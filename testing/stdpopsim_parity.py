#!/usr/bin/env python3
"""Compare discoal and msprime on demes models.

Models come from three sources:

  - every demographic model in the stdpopsim catalog;
  - random instances of stdpopsim's generic models PiecewiseConstantSize
    and IsolationWithMigration (--generic N of each, drawn with --seed);
  - every *.yaml demes file under each --demes-dir, for example the
    examples and test cases of the demes specification:
        git clone https://github.com/popsim-consortium/demes-spec
        python testing/stdpopsim_parity.py \\
            --demes-dir demes-spec/examples --demes-dir demes-spec/test-cases/valid

For each model:

  1. Write it as a demes file.
  2. Simulate `reps` genealogies with discoal (through a -Y config that
     loads the demes file, with tree sequence output) and with msprime
     (Demography.from_demes), with recombination if --rho > 0. --rho is
     4 N_max r L, where N_max is the largest deme size in the model; discoal
     gets the same r as its own rho = 4 N_ref r L. Where that would give
     fewer than about MIN_TREES recombination breakpoints per replicate (a
     very large size somewhere), r is raised to reach it.
  3. Compute branch-mode statistics on each tree sequence: diversity,
     segregating sites (total branch length) and Tajima's D over all
     samples, diversity within each set of samples ("diversity[A]" for deme
     A at time 0, "diversity[A@320]" for deme A 320 generations ago), and,
     with recombination, the number of distinct trees along the sequence.
  4. Compare the two distributions of each statistic with a two-sample KS
     test, Bonferroni-corrected over all tests in the run.

The comparison is on absolute scale. discoal writes node times in units of
2 N_ref generations, where N_ref is the present-day size of the first
present-day deme; they are converted to generations before the statistics
are computed. A time-unit or size-scaling error therefore shows up as a
shift in diversity and segregating sites, not only in Tajima's D.

Samples: sets of `--samples` haploid genomes, each from one deme at one
time; both simulators get the same sets.

  - Every deme that exists at time 0 gives a set at time 0.
  - In a stdpopsim catalog model, a population with a default sampling time
    later than 0 is sampled at that time instead (ancient DNA samples, for
    example the Neanderthal in HomSap/AncientEurasia_9K19). This is the
    sampling its authors give, including for extinct demes.
  - Each generic model and demes file also runs as a second model,
    "<label>+ancient", which adds one set from every deme, extinct demes
    included, at a random time inside the deme's lifetime (away from its
    ends; for a deme with no start time, up to twice its oldest epoch's
    size in generations after the end of that epoch). The times are drawn from a seed made
    from the label, so they are the same in every run.

discoal gets the sets at times after 0 as YAML ancient_samples. Its
deme_sample_size is indexed by the order of demes in the graph, extinct
demes included, and counts the ancient samples too. The statistics find
each set's samples by population and node time, and stop with an error if
a set does not hold exactly `--samples` nodes. A model is skipped,
with the reason, if it uses features that discoal does not support
(selfing, cloning, size functions other than constant and exponential), has
no deme at time 0, or if msprime cannot simulate it (stranded lineages, or
node times beyond floating-point precision). In a model that can
leave lineages in demes that never exchange migrants again, a rare
replicate strands them: msprime then reports an infinite waiting time and
discoal stops with an error. Either one skips the model.

A model with a deme size below 100 is first scaled up: sizes and times are
multiplied by c and migration rates divided by c, which leaves the
coalescent unchanged. discoal's output does not depend on this scaling, but
msprime's does at sizes of about one individual (several demes
specification test cases use sizes of 0.1 to 1); see rescale().

Usage:
    python testing/stdpopsim_parity.py [--reps 1000] [--samples 10]
        [--workers 24] [--generic 20] [--seed 1] [--demes-dir DIR ...]
        [--rho 0] [--sites 1000] [--models LABEL,...]
        [--sample K --sample-seed S] [--timeout 600]
        [--cache-dir DIR | --no-cache] [--tsv out.tsv] [--jobs-tsv jobs.tsv]

Each simulation's statistics are cached (default testing/.parity_cache), so
a rerun repeats only what changed: after a discoal rebuild, only discoal
runs again. Each model's simulation times are printed as it finishes,
summarized at the end, and written with --jobs-tsv.

Needs build/discoal, stdpopsim, msprime, demes, tskit and scipy.
Exit code 0 if no test rejects and no model errors, 1 otherwise.
"""

import argparse
import concurrent.futures
import contextlib
import glob
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
import warnings
import zlib

import demes
import msprime
import numpy as np
import stdpopsim
import tskit
from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
ALPHA = 0.01


def present_demes(graph):
    """Names of the demes that exist at time 0, in graph order."""
    return [d.name for d in graph.demes if d.end_time == 0]


def sample_sets(kind, payload, graph, scale, rng):
    """(deme name, time in generations) for each set of samples; see the
    module docstring. `scale` is rescale()'s factor, which also applies to
    stdpopsim's sampling times. `rng` draws a random ancient set from every
    deme, or is None."""
    times = {}
    if kind == "catalog":
        species, model_id = payload
        model = stdpopsim.get_species(species).get_demographic_model(model_id)
        times = {pop.name: pop.default_sampling_time * scale
                 for pop in model.populations
                 if pop.default_sampling_time}
    sets = [(d.name, times.get(d.name, 0.0)) for d in graph.demes
            if d.name in times or d.end_time == 0]
    if rng is not None:
        # Graph times are in its time_units (generations or years); sample
        # times are in generations.
        gt = graph.generation_time
        for d in graph.demes:
            lo = d.end_time
            # With no start time, go past the end of the oldest epoch by
            # twice its size: the oldest size alone can be tiny next to
            # N_ref (valid/growth_rates_03.yaml).
            hi = (d.start_time if np.isfinite(d.start_time)
                  else max(lo, d.epochs[0].end_time)
                  + 2 * d.epochs[0].start_size * gt)
            span = hi - lo
            t = rng.uniform(lo + 0.01 * span, hi - 0.01 * span)
            sets.append((d.name, float(t / gt)))
    return sets


def set_label(name, t):
    """Statistic label of a set: the deme name, with "@time" if ancient."""
    return name if t == 0 else f"{name}@{t:.6g}"


def unsupported(graph):
    """Why discoal cannot run this graph, or None."""
    epochs = [e for d in graph.demes for e in d.epochs]
    if any(e.size_function not in ("constant", "exponential") for e in epochs):
        return "size function other than constant or exponential"
    if any(e.selfing_rate > 0 or e.cloning_rate > 0 for e in epochs):
        return "selfing or cloning"
    if not present_demes(graph):
        return "no deme exists at time 0"
    return None


MIN_SIZE = 100


def rescale(graph):
    """Return (graph, c): if some deme size is below MIN_SIZE, the graph with
    sizes and times multiplied by c and migration rates divided by c, so the
    smallest size is MIN_SIZE. The coalescent does not change under this
    scaling, and discoal's output does not depend on it, but msprime's does
    for sizes of about one individual or less."""
    sizes = [x for d in graph.demes for e in d.epochs
             for x in (e.start_size, e.end_size)]
    c = MIN_SIZE / min(sizes)
    if c <= 1:
        return graph, 1.0
    g = graph.asdict()
    for d in g["demes"]:
        d["start_time"] = d["start_time"] * c
        for e in d["epochs"]:
            e["end_time"] *= c
            e["start_size"] *= c
            e["end_size"] *= c
    for m in g.get("migrations", []):
        m["start_time"] = m["start_time"] * c
        m["end_time"] *= c
        m["rate"] /= c
    for pulse in g.get("pulses", []):
        pulse["time"] *= c
    return demes.Builder.fromdict(g).resolve(), c


def build_graph(kind, payload):
    """The demes graph for one job."""
    if kind == "catalog":
        species, model_id = payload
        model = stdpopsim.get_species(species).get_demographic_model(model_id)
        return model.model.to_demes()
    if kind == "pcs":
        return stdpopsim.PiecewiseConstantSize(
            payload["N0"], *payload["changes"]).model.to_demes()
    if kind == "im":
        return stdpopsim.IsolationWithMigration(**payload).model.to_demes()
    return demes.load(payload)


def random_generic(rng, n):
    """n random PiecewiseConstantSize and n IsolationWithMigration models:
    sizes log-uniform in [1e2, 1e5], times log-uniform in [50, 2e5]
    generations, migration rates 0 (one time in three) or log-uniform in
    [1e-6, 1e-3] per generation."""
    def size():
        return float(10 ** rng.uniform(2, 5))

    def time():
        return float(10 ** rng.uniform(np.log10(50), np.log10(2e5)))

    def mig():
        return 0.0 if rng.random() < 1 / 3 else float(10 ** rng.uniform(-6, -3))

    jobs = []
    for i in range(n):
        times = sorted(time() for _ in range(rng.integers(1, 5)))
        jobs.append(("pcs", f"generic/PiecewiseConstantSize_{i}",
                     dict(N0=size(), changes=[(t, size()) for t in times])))
    for i in range(n):
        jobs.append(("im", f"generic/IsolationWithMigration_{i}",
                     dict(NA=size(), N1=size(), N2=size(), T=time(),
                          M12=mig(), M21=mig())))
    return jobs


def set_nodes(ts, set_ids, n):
    """Sample node arrays, one per (population id, node time) in set_ids;
    each must hold exactly n nodes."""
    out = []
    node_time = ts.tables.nodes.time
    for p, t in set_ids:
        s = ts.samples(population=p)
        # No absolute tolerance: an ancient time near 0 is not time 0.
        s = s[np.isclose(node_time[s], t, rtol=1e-6, atol=0)]
        if len(s) != n:
            raise RuntimeError(f"{len(s)} sample nodes in population {p} at "
                               f"time {t:g}, expected {n}")
        out.append(s)
    return out


def collect(tree_seqs, set_ids, labels, n, time_scale, count_trees):
    """Statistics for each tree sequence; returns {stat name: array}.
    set_ids holds (population id, node time in the tree sequence's units)
    for each set of samples. With count_trees (recombination on), every
    replicate also records its number of distinct trees, which tracks the
    recombination rate."""
    out = {}
    for ts in tree_seqs:
        sets = set_nodes(ts, set_ids, n)
        allset = np.concatenate(sets)
        vals = {
            "diversity": ts.diversity(allset, mode="branch") * time_scale,
            "segregating_sites":
                ts.segregating_sites(allset, mode="branch") * time_scale,
            "tajimas_d": ts.Tajimas_D(allset, mode="branch"),
        }
        if count_trees:
            vals["num_trees"] = ts.num_trees
        if len(sets) > 1:
            for name, s in zip(labels, sets):
                vals[f"diversity[{name}]"] = (
                    ts.diversity(s, mode="branch") * time_scale)
        for k, v in vals.items():
            out.setdefault(k, []).append(float(v))
    return {k: np.array(v) for k, v in out.items()}


def run_discoal(demes_path, graph, sets, n, reps, seed, tmp, timeout,
                rho, sites):
    """Run discoal; return ({stat: array}, N_ref as discoal reports it)."""
    index = {d.name: i for i, d in enumerate(graph.demes)}
    sizes = [n * sum(name == d.name for name, _ in sets) for d in graph.demes]
    # YAML times are in units of 4 N_ref generations.
    n_ref = reference_size(graph)
    ancient = ", ".join(
        f"{{time: {float(t / (4 * n_ref))!r}, population: {index[name]}, "
        f"sample_size: {n}}}" for name, t in sets if t > 0)
    events = (f", demographic_events: {{ancient_samples: [{ancient}]}}"
              if ancient else "")
    cfg = os.path.join(tmp, "config.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {sum(sizes)}, "
                f"num_replicates: {reps}, num_sites: {sites}, "
                f"seed: [{seed}, {seed + 1}]}}\n"
                # Branch-mode statistics need no mutations; with tall trees a
                # positive rate can exceed discoal's mutation limit.
                f"genetics: {{mutation_rate: 0, recombination_rate: {rho}}}\n"
                f"demography: {{deme_sample_size: {sizes}, "
                f"demes_filename: \"{demes_path}\"{events}}}\n")
    prefix = os.path.join(tmp, "d.trees")
    r = subprocess.run([DISCOAL, "-v", "-Y", cfg, "-ts", prefix],
                       capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"discoal exit {r.returncode}: "
                           f"{r.stderr.strip()[-300:]}")
    m = re.search(r"as reference with N=([0-9.eE+-]+)", r.stderr)
    if not m:
        raise RuntimeError("could not find N_ref in discoal -v output")
    # The log prints N_ref to 6 significant digits; audit_model checks it
    # against reference_size(), and the exact value scales the times here.
    reported_n_ref = float(m.group(1))
    files = [os.path.join(tmp, f"d_rep{i}.trees") for i in range(1, reps + 1)]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        # discoal can exit 0 after a failed tree sequence write.
        errs = [l for l in r.stderr.splitlines() if "rror" in l]
        raise RuntimeError(f"{len(missing)} of {reps} tree files missing; "
                           f"{errs[-1] if errs else 'no error message'}")
    # Node times are in units of 2 N_ref generations.
    set_ids = [(index[name], t / (2 * n_ref)) for name, t in sets]
    res = collect((tskit.load(f) for f in files), set_ids,
                  [set_label(*x) for x in sets], n, 2 * n_ref, rho > 0)
    for f in files:
        os.remove(f)
    return res, reported_n_ref


def reference_size(graph):
    """discoal's N_ref: the present-day size of the first deme at time 0."""
    return next(d.epochs[-1].end_size for d in graph.demes if d.end_time == 0)


def largest_size(graph):
    return max(x for d in graph.demes for e in d.epochs
               for x in (e.start_size, e.end_size))


def recombination(graph, rho, sites):
    """(rate per generation for msprime, rho for discoal). --rho is
    4 N_max r L, with N_max the largest deme size in the model, so the
    number of trees stays moderate even where some size is far above
    discoal's N_ref. discoal's rho is 4 N_ref r L over its L = `sites`
    sites (units_convention_check.py pins this)."""
    if rho <= 0:
        return 0.0, 0.0
    r = rho / (4 * largest_size(graph) * sites)
    return r, 4 * reference_size(graph) * r * sites


def run_msprime(graph, sets, n, reps, seed, r, sites):
    """msprime with the same samples, recombination rate r per generation
    on a discrete genome of length `sites`."""
    dem = msprime.Demography.from_demes(graph)
    samples = [msprime.SampleSet(n, population=name, ploidy=1, time=t)
               for name, t in sets]
    set_ids = [(dem[name].id, t) for name, t in sets]
    tss = msprime.sim_ancestry(samples=samples, demography=dem, ploidy=2,
                               sequence_length=sites, recombination_rate=r,
                               num_replicates=reps, random_seed=seed)
    return collect(tss, set_ids, [set_label(*x) for x in sets], n, 1.0, r > 0)


# demes specification test cases with near-zero migration between demes:
# the trees are extremely tall, and with recombination neither simulator
# finishes 2 replicates in 2 minutes (valid/migration_05.yaml, --rho 100
# --sites 10000). Without recombination they take under a second each and
# pass, so they are skipped only when --rho > 0.
SLOW_WITH_RECOMBINATION = {
    "valid/basic_resolution_04.yaml", "valid/basic_resolution_05.yaml",
    "valid/migration_05.yaml", "valid/migration_06.yaml",
    "valid/migration_07.yaml", "valid/migration_08.yaml",
    "valid/structure_02.yaml", "valid/structure_06.yaml",
}


HARNESS_ID = (hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
              msprime.__version__, tskit.__version__, demes.__version__,
              stdpopsim.__version__)

# With recombination, a model whose rho (scaled to its largest deme size)
# would give fewer than about this many recombination breakpoints per
# replicate gets its recombination rate raised, so that the tree-count test
# means something (a very large size somewhere makes rho tiny elsewhere).
MIN_TREES = 20


class TimeLimit(Exception):
    pass


@contextlib.contextmanager
def time_limit(seconds):
    """Raise TimeLimit in this process after `seconds` (msprime checks for
    signals while it runs, so this interrupts a long simulation)."""
    def on_alarm(signum, frame):
        raise TimeLimit()
    old = signal.signal(signal.SIGALRM, on_alarm)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def cache_key(*parts):
    """Key for a cached simulation. It includes this script's code and the
    versions of the libraries that build the models and compute the
    statistics, so a change to any of them is not served old results."""
    return hashlib.sha256(json.dumps((HARNESS_ID, parts), sort_keys=True,
                                     default=str).encode()).hexdigest()


def run_cached(cache_dir, key, limit, compute, deterministic):
    """Run compute() (which returns {stat: array} or ({stat: array}, extra))
    under the time limit, or reuse a cached result. Returns (result, cached),
    where result has "stats" and "extra", or "timed_out", or "failed", and
    always "seconds" (the time of the original run). A cached time-out is
    reused only if the limit has not been raised since. A failure is cached
    only if deterministic(exception) says it would happen again (a model
    msprime cannot simulate, a discoal error); others, such as running out of
    memory, are retried on the next run."""
    path = os.path.join(cache_dir, key[:2], key + ".json") if cache_dir else None
    if path and os.path.exists(path):
        with open(path) as f:
            hit = json.load(f)
        if not (hit.get("timed_out") and limit > hit["limit"]):
            return hit, True
    t0 = time.perf_counter()
    try:
        with time_limit(limit):
            value = compute()
        stats_, extra = value if isinstance(value, tuple) else (value, None)
        out = {"stats": {k: v.tolist() for k, v in stats_.items()},
               "extra": extra}
    except (TimeLimit, subprocess.TimeoutExpired):
        out = {"timed_out": True, "limit": limit}
    except Exception as e:
        out = {"failed": f"{type(e).__name__}: {e}"[:300]}
        if not deterministic(e):
            out["seconds"] = time.perf_counter() - t0
            return out, False
    out["seconds"] = time.perf_counter() - t0
    if path:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + f".{os.getpid()}"
        with open(tmp, "w") as f:
            json.dump(out, f)
        os.replace(tmp, path)
    return out, False


def audit_model(job):
    """Run one model. Returns a dict with rows of test results, an error, or
    a reason to skip it, and the time of each simulation."""
    (kind, label, payload, ancient, reps, n, seed, timeout, rho, sites,
     cache_dir, discoal_id) = job
    timings = {}
    try:
        if rho > 0 and label.removesuffix("+ancient") in SLOW_WITH_RECOMBINATION:
            return dict(model=label, timings=timings, skipped="too slow with "
                        "recombination (near-zero migration, very tall trees)")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graph = build_graph(kind, payload)
        reason = unsupported(graph)
        if reason:
            return dict(model=label, timings=timings, skipped=reason)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graph, scale = rescale(graph)
            model_id = demes.dumps(graph)
        rng = (np.random.default_rng(zlib.crc32(label.encode()))
               if ancient else None)
        sets = sample_sets(kind, payload, graph, scale, rng)
        r, discoal_rho = recombination(graph, rho, sites)
        ms_id = ("msprime", msprime.__version__)

        def simulate(name, key_parts, compute, deterministic):
            res, cached = run_cached(cache_dir, cache_key(name, *key_parts),
                                     timeout, compute, deterministic)
            timings[name] = (res["seconds"], cached)
            return res

        def quiet_msprime(rate, n_reps):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return run_msprime(graph, sets, n, n_reps, seed, rate, sites)

        def msprime_cannot(e):
            # Stranded lineages, or node times that floating point cannot
            # tell apart (a parent at the same time as its child), which
            # happens with ancient samples at times near 1e25 generations
            # in models with sizes near 1e-20 (valid/growth_rates_03.yaml).
            return ("Infinite waiting time" in str(e)
                    or "time value <= to its child" in str(e))

        def msprime_result(name, rate, n_reps):
            """Run (or reuse) an msprime simulation; returns (stats, outcome)
            where outcome is None, or a skip or an error for the model."""
            res = simulate(name, (model_id, sets, n, n_reps, seed, sites,
                                  rate, ms_id),
                           lambda: quiet_msprime(rate, n_reps), msprime_cannot)
            if "stats" in res:
                return {k: np.array(v) for k, v in res["stats"].items()}, None
            if res.get("timed_out") or msprime_cannot(res["failed"]):
                return None, dict(model=label, timings=timings,
                                  skipped=skip_reason("msprime", res, timeout))
            return None, dict(model=label, timings=timings,
                              error="msprime failed: " + res["failed"])

        rho_factor = 1.0
        if r > 0:
            # A model that can strand lineages makes msprime report "infinite
            # waiting time" without recombination; with recombination the
            # stranded lineages keep recombining and msprime never stops.
            # Find these models first.
            check, outcome = msprime_result("msprime_check", 0.0, reps)
            if outcome:
                return outcome
            # The expected number of recombination breakpoints on the
            # genealogy is about r * L * T, with T the total branch length in
            # generations (branch-mode segregating sites without
            # recombination). Raise r to give at least MIN_TREES.
            r_min = MIN_TREES / (sites * check["segregating_sites"].mean())
            if r < r_min:
                rho_factor = r_min / r
                r *= rho_factor
                discoal_rho *= rho_factor
        m, outcome = msprime_result("msprime", r, reps)
        if outcome:
            return outcome

        def run_d():
            with tempfile.TemporaryDirectory() as tmp:
                path = os.path.join(tmp, "model.demes.yaml")
                demes.dump(graph, path)
                return run_discoal(path, graph, sets, n, reps, seed, tmp,
                                   timeout, discoal_rho, sites)
        # discoal errors (non-zero exit, missing tree files) repeat for the
        # same binary and inputs; other exceptions are retried.
        res = simulate("discoal", (model_id, sets, n, reps, seed, sites,
                                   discoal_rho, discoal_id), run_d,
                       lambda e: isinstance(e, RuntimeError))
        if res.get("timed_out"):
            return dict(model=label, timings=timings,
                        error=f"discoal did not finish in {timeout} s")
        if "failed" in res:
            # In a model where lineages can be left in demes that never
            # exchange migrants again, a rare replicate strands them. discoal
            # stops there as msprime does ("infinite waiting time"), which
            # skips the model above; a run of msprime's replicates can miss
            # a replicate that discoal's run hits, so skip it here too.
            if "can never coalesce" in res["failed"]:
                return dict(model=label, timings=timings,
                            skipped="discoal: lineages stranded in one "
                                    "replicate (msprime: infinite waiting "
                                    "time)")
            return dict(model=label, timings=timings, error=res["failed"])
        d = {k: np.array(v) for k, v in res["stats"].items()}
        n_ref = res["extra"]
        if not np.isclose(n_ref, reference_size(graph)):
            raise RuntimeError(f"discoal N_ref {n_ref} differs from "
                               f"{reference_size(graph)}")
        missing = sorted(set(d) ^ set(m))
        if missing:
            return dict(model=label, timings=timings,
                        error="statistics missing from one simulator: "
                              + ", ".join(missing))
        rows = []
        for stat in sorted(d):
            ks = stats.ks_2samp(d[stat], m[stat])
            rows.append(dict(model=label, stat=stat,
                             discoal_mean=d[stat].mean(),
                             msprime_mean=m[stat].mean(),
                             ks_d=ks.statistic, p=ks.pvalue))
        return dict(model=label, n_ref=n_ref,
                    demes=[set_label(*x) for x in sets], rows=rows,
                    scale=scale, rho_factor=rho_factor, timings=timings)
    except Exception as e:  # report and keep going with other models
        return dict(model=label, timings=timings,
                    error=f"{type(e).__name__}: {e}")


def skip_reason(simulator, res, timeout):
    if res.get("timed_out"):
        return f"{simulator} did not finish in {timeout} s"
    return f"{simulator} cannot simulate it: {res['failed']}"[:200]


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def timing_text(timings):
    return ", ".join(f"{name} {sec:.1f} s{' (cached)' if cached else ''}"
                     for name, (sec, cached) in timings.items())


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--samples", type=int, default=10,
                    help="haploid samples per sampled deme")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--generic", type=int, default=20,
                    help="random instances of each generic stdpopsim model")
    ap.add_argument("--seed", type=int, default=1,
                    help="seed for the random generic models")
    ap.add_argument("--demes-dir", action="append", default=[],
                    help="audit every *.yaml demes file under this directory "
                         "(may be given more than once)")
    ap.add_argument("--models", default="",
                    help="comma-separated model labels (as printed) to run; "
                         "default all")
    ap.add_argument("--sample", type=int, default=0,
                    help="also run K models drawn at random from the "
                         "others (with --models), or from all")
    ap.add_argument("--sample-seed", type=int, default=0,
                    help="seed for --sample; the same seed draws the same "
                         "models")
    ap.add_argument("--rho", type=float, default=0.0,
                    help="recombination: 4 N_max r L over the whole locus, "
                         "N_max the largest deme size in each model")
    ap.add_argument("--sites", type=int, default=1000,
                    help="number of sites (discoal) and sequence length "
                         "(msprime)")
    ap.add_argument("--timeout", type=int, default=600,
                    help="seconds allowed for each simulation (msprime or "
                         "discoal); a model over the limit is skipped "
                         "(msprime) or reported as an error (discoal)")
    ap.add_argument("--cache-dir",
                    default=os.path.join(REPO, "testing", ".parity_cache"),
                    help="reuse simulation results stored here; the key "
                         "includes the model, settings, the msprime version "
                         "and a hash of the discoal binary")
    ap.add_argument("--no-cache", action="store_true",
                    help="run every simulation and store nothing")
    ap.add_argument("--tsv", help="write per-test results to this file")
    ap.add_argument("--jobs-tsv",
                    help="write per-model status and simulation times to "
                         "this file")
    args = ap.parse_args()
    if args.reps < 2:
        sys.exit("--reps must be at least 2")
    cache_dir = None if args.no_cache else args.cache_dir
    discoal_id = ("discoal", file_sha256(DISCOAL))

    # (kind, label, payload, random ancient samples)
    specs = [("catalog", f"{sp.id}/{mdl.id}", (sp.id, mdl.id), False)
             for sp in stdpopsim.all_species() for mdl in sp.demographic_models]
    others = random_generic(np.random.default_rng(args.seed), args.generic)
    for d in args.demes_dir:
        root = os.path.dirname(os.path.abspath(d))
        for path in sorted(glob.glob(os.path.join(d, "**", "*.yaml"),
                                     recursive=True)):
            others.append(("file", os.path.relpath(path, root), path))
    specs += [(kind, label, payload, False) for kind, label, payload in others]
    specs += [(kind, label + "+ancient", payload, True)
              for kind, label, payload in others]
    wanted = set(filter(None, args.models.split(",")))
    unknown = wanted - {label for _, label, _, _ in specs}
    if unknown:
        sys.exit(f"unknown model labels: {', '.join(sorted(unknown))}")
    if args.sample > 0:
        rest = [label for _, label, _, _ in specs if label not in wanted]
        drawn = np.random.default_rng(args.sample_seed).choice(
            len(rest), size=min(args.sample, len(rest)), replace=False)
        sample = sorted(rest[i] for i in drawn)
        print(f"Random sample of {len(sample)} models (--sample "
              f"{args.sample} --sample-seed {args.sample_seed}): "
              f"{', '.join(sample)}\n", flush=True)
        wanted |= set(sample)
    jobs = [(kind, label, payload, ancient, args.reps, args.samples,
             1000 + 10 * i, args.timeout, args.rho, args.sites, cache_dir,
             discoal_id)
            for i, (kind, label, payload, ancient) in enumerate(specs)
            if not wanted or label in wanted]
    # Report each model as it finishes; long runs then show progress.
    results = []
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(audit_model, job) for job in jobs]
        for k, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            res = fut.result()
            results.append(res)
            state = ("skipped" if "skipped" in res else
                     "error" if "error" in res else "done")
            print(f"[{k}/{len(jobs)}] {res['model']}: {state}"
                  f"{' (' + timing_text(res['timings']) + ')' if res['timings'] else ''}",
                  file=sys.stderr, flush=True)
    order = {job[1]: i for i, job in enumerate(jobs)}
    results.sort(key=lambda res: order[res["model"]])

    rows = [r for res in results if "rows" in res for r in res["rows"]]
    errors = [res for res in results if "error" in res]
    skipped = [res for res in results if "skipped" in res]
    n_tests = len(rows)
    threshold = ALPHA / max(n_tests, 1)
    print(f"{len(jobs)} models, {n_tests} tests, {args.reps} reps, "
          f"Bonferroni threshold p < {threshold:.2e}\n")
    print(f"{'model':<55} {'stat':<28} {'discoal':>11} {'msprime':>11} "
          f"{'ratio':>6} {'KS D':>6} {'p':>9}")
    rejected = []
    for r in rows:
        ratio = (r["discoal_mean"] / r["msprime_mean"]
                 if r["msprime_mean"] != 0 else float("nan"))
        flag = r["p"] < threshold
        if flag:
            rejected.append(r)
        print(f"{r['model']:<55} {r['stat']:<28} {r['discoal_mean']:>11.4g} "
              f"{r['msprime_mean']:>11.4g} {ratio:>6.3f} {r['ks_d']:>6.3f} "
              f"{r['p']:>9.2e}{'  REJECT' if flag else ''}")
    for e in errors:
        print(f"{e['model']:<55} ERROR {e['error']}")
    for s in skipped:
        print(f"{s['model']:<55} SKIPPED {s['skipped']}")
    if args.tsv:
        with open(args.tsv, "w") as f:
            f.write("model\tstat\tdiscoal_mean\tmsprime_mean\tks_d\tp\t"
                    "reject\n")
            for r in rows:
                f.write(f"{r['model']}\t{r['stat']}\t{r['discoal_mean']:.6g}"
                        f"\t{r['msprime_mean']:.6g}\t{r['ks_d']:.4f}\t"
                        f"{r['p']:.3e}\t{int(r['p'] < threshold)}\n")
    if args.jobs_tsv:
        with open(args.jobs_tsv, "w") as f:
            f.write("model\tstatus\tmsprime_check_s\tmsprime_s\tdiscoal_s\t"
                    "cached\trho_factor\treason\n")
            for res in results:
                t = res["timings"]
                status = ("skipped" if "skipped" in res else
                          "error" if "error" in res else "done")
                secs = [f"{t[k][0]:.2f}" if k in t else ""
                        for k in ("msprime_check", "msprime", "discoal")]
                cached = ",".join(k for k, (sec, c) in t.items() if c)
                reason = res.get("skipped", res.get("error", ""))
                factor = f"{res['rho_factor']:.3g}" if "rho_factor" in res else ""
                f.write(f"{res['model']}\t{status}\t" + "\t".join(secs)
                        + f"\t{cached}\t{factor}\t{reason}\n")
    computed = {}
    for res in results:
        for name, (sec, cached) in res["timings"].items():
            tot = computed.setdefault(name, [0.0, 0, 0])
            tot[0] += sec
            tot[1 if cached else 2] += 1
    print("\nSimulation time (original run times, cached or not):")
    for name, (sec, n_cached, n_run) in computed.items():
        print(f"  {name:<14} {sec / 60:8.1f} min total; {n_run} run now, "
              f"{n_cached} from cache")
    slowest = sorted(results, key=lambda res: -sum(
        sec for sec, _ in res["timings"].values()))[:5]
    print("  slowest models: " + "; ".join(
        f"{res['model']} {sum(sec for sec, _ in res['timings'].values()):.0f} s"
        for res in slowest if res["timings"]))
    trees = [r["msprime_mean"] for r in rows if r["stat"] == "num_trees"]
    if trees:
        raised = sum(1 for res in results if res.get("rho_factor", 1) > 1)
        print(f"Recombination: {len(trees)} models compared on the number of "
              f"trees; msprime mean trees per replicate min {min(trees):.1f}, "
              f"median {np.median(trees):.0f}; rate raised (to about "
              f"{MIN_TREES} breakpoints per replicate) in {raised} models.")
    bad = {r["model"] for r in rejected} | {e["model"] for e in errors}
    tested = len(results) - len(skipped)
    print(f"\n{tested - len(bad)} of {tested} tested models pass; "
          f"{len(rejected)} tests reject; {len(errors)} models error; "
          f"{len(skipped)} skipped.")
    sys.exit(0 if not bad else 1)


if __name__ == "__main__":
    main()
