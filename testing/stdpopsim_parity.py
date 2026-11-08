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
     (Demography.from_demes), with recombination if --rho > 0.
  3. Compute branch-mode statistics on each tree sequence: diversity,
     segregating sites (total branch length) and Tajima's D over all
     samples, diversity within each sampled deme, and, with recombination,
     the number of distinct trees along the sequence.
  4. Compare the two distributions of each statistic with a two-sample KS
     test, Bonferroni-corrected over all tests in the run.

The comparison is on absolute scale. discoal writes node times in units of
2 N_ref generations, where N_ref is the present-day size of the first
present-day deme; they are converted to generations before the statistics
are computed. A time-unit or size-scaling error therefore shows up as a
shift in diversity and segregating sites, not only in Tajima's D.

Samples: `--samples` haploid genomes from every deme that exists at time 0,
whatever stdpopsim's default sampling time; both simulators get the same
samples. discoal's deme_sample_size is indexed by the order of demes in the
graph, extinct demes included, so other demes get 0. A model is skipped,
with the reason, if it uses features that discoal does not support
(selfing, cloning, size functions other than constant and exponential), has
no deme at time 0, or if msprime cannot simulate it.

A model with a deme size below 100 is first scaled up: sizes and times are
multiplied by c and migration rates divided by c, which leaves the
coalescent unchanged. discoal's output does not depend on this scaling, but
msprime's does at sizes of about one individual (several demes
specification test cases use sizes of 0.1 to 1); see rescale().

Usage:
    python testing/stdpopsim_parity.py [--reps 1000] [--samples 10]
        [--workers 24] [--generic 20] [--seed 1] [--demes-dir DIR ...]
        [--rho 0] [--sites 1000] [--models LABEL,...] [--timeout 600]
        [--tsv out.tsv]

Needs build/discoal, stdpopsim, msprime, demes, tskit and scipy.
Exit code 0 if no test rejects and no model errors, 1 otherwise.
"""

import argparse
import concurrent.futures
import glob
import os
import re
import subprocess
import sys
import tempfile
import warnings

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


def per_deme_sets(ts, pop_ids):
    """Sample node arrays, one per population id, in the given order."""
    return [ts.samples(population=p) for p in pop_ids]


def collect(tree_seqs, pop_ids, names, time_scale):
    """Statistics for each tree sequence; returns {stat name: array}."""
    out = {}
    for ts in tree_seqs:
        sets = per_deme_sets(ts, pop_ids)
        allset = np.concatenate(sets)
        vals = {
            "diversity": ts.diversity(allset, mode="branch") * time_scale,
            "segregating_sites":
                ts.segregating_sites(allset, mode="branch") * time_scale,
            "tajimas_d": ts.Tajimas_D(allset, mode="branch"),
        }
        if ts.num_trees > 1:
            # Distinct trees along the sequence: tracks the recombination rate.
            vals["num_trees"] = ts.num_trees
        if len(sets) > 1:
            for name, s in zip(names, sets):
                vals[f"diversity[{name}]"] = (
                    ts.diversity(s, mode="branch") * time_scale)
        for k, v in vals.items():
            out.setdefault(k, []).append(float(v))
    return {k: np.array(v) for k, v in out.items()}


def run_discoal(demes_path, graph, names, n, reps, seed, tmp, timeout,
                rho, sites):
    """Run discoal; return ({stat: array}, N_ref)."""
    index = {d.name: i for i, d in enumerate(graph.demes)}
    sizes = [n if d.name in names else 0 for d in graph.demes]
    cfg = os.path.join(tmp, "config.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {sum(sizes)}, "
                f"num_replicates: {reps}, num_sites: {sites}, "
                f"seed: [{seed}, {seed + 1}]}}\n"
                # Branch-mode statistics need no mutations; with tall trees a
                # positive rate can exceed discoal's mutation limit.
                f"genetics: {{mutation_rate: 0, recombination_rate: {rho}}}\n"
                f"demography: {{deme_sample_size: {sizes}, "
                f"demes_filename: \"{demes_path}\"}}\n")
    prefix = os.path.join(tmp, "d.trees")
    r = subprocess.run([DISCOAL, "-v", "-Y", cfg, "-ts", prefix],
                       capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"discoal exit {r.returncode}: "
                           f"{r.stderr.strip()[-300:]}")
    m = re.search(r"as reference with N=([0-9.eE+-]+)", r.stderr)
    if not m:
        raise RuntimeError("could not find N_ref in discoal -v output")
    n_ref = float(m.group(1))
    files = [os.path.join(tmp, f"d_rep{i}.trees") for i in range(1, reps + 1)]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        # discoal can exit 0 after a failed tree sequence write.
        errs = [l for l in r.stderr.splitlines() if "rror" in l]
        raise RuntimeError(f"{len(missing)} of {reps} tree files missing; "
                           f"{errs[-1] if errs else 'no error message'}")
    pop_ids = [index[name] for name in names]
    # Node times are in units of 2 N_ref generations.
    res = collect((tskit.load(f) for f in files), pop_ids, names, 2 * n_ref)
    for f in files:
        os.remove(f)
    return res, n_ref


def reference_size(graph):
    """discoal's N_ref: the present-day size of the first deme at time 0."""
    return next(d.epochs[-1].end_size for d in graph.demes if d.end_time == 0)


def run_msprime(graph, names, n, reps, seed, rho, sites):
    """msprime with the same samples. discoal's rho is 4 N_ref r L over its
    L = `sites` sites (units_convention_check.py pins this), so msprime gets
    r = rho / (4 N_ref L) on a discrete genome of length L."""
    dem = msprime.Demography.from_demes(graph)
    samples = [msprime.SampleSet(n, population=name, ploidy=1, time=0)
               for name in names]
    pop_ids = [dem[name].id for name in names]
    r = rho / (4 * reference_size(graph) * sites) if rho > 0 else 0.0
    tss = msprime.sim_ancestry(samples=samples, demography=dem, ploidy=2,
                               sequence_length=sites, recombination_rate=r,
                               num_replicates=reps, random_seed=seed)
    return collect(tss, pop_ids, names, 1.0)


def audit_model(job):
    """Run one model. Returns a dict with rows of test results, an error, or
    a reason to skip it."""
    kind, label, payload, reps, n, seed, timeout, rho, sites = job
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graph = build_graph(kind, payload)
        reason = unsupported(graph)
        if reason:
            return dict(model=label, skipped=reason)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graph, scale = rescale(graph)
        names = present_demes(graph)
        # msprime first: a graph it cannot simulate is skipped, not an error.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = run_msprime(graph, names, n, reps, seed, rho, sites)
        except Exception as e:
            return dict(model=label,
                        skipped=f"msprime cannot simulate it: {e}"[:200])
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "model.demes.yaml")
            demes.dump(graph, path)
            d, n_ref = run_discoal(path, graph, names, n, reps, seed, tmp,
                                   timeout, rho, sites)
        if not np.isclose(n_ref, reference_size(graph)):
            raise RuntimeError(f"discoal N_ref {n_ref} differs from "
                               f"{reference_size(graph)}")
        rows = []
        for stat in d:
            ks = stats.ks_2samp(d[stat], m[stat])
            rows.append(dict(model=label, stat=stat,
                             discoal_mean=d[stat].mean(),
                             msprime_mean=m[stat].mean(),
                             ks_d=ks.statistic, p=ks.pvalue))
        return dict(model=label, n_ref=n_ref, demes=names, rows=rows,
                    scale=scale)
    except subprocess.TimeoutExpired:
        return dict(model=label, error=f"discoal did not finish in {timeout} s")
    except Exception as e:  # report and keep going with other models
        return dict(model=label, error=f"{type(e).__name__}: {e}")


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
    ap.add_argument("--rho", type=float, default=0.0,
                    help="recombination: 4 N_ref r over the whole locus")
    ap.add_argument("--sites", type=int, default=1000,
                    help="number of sites (discoal) and sequence length "
                         "(msprime)")
    ap.add_argument("--timeout", type=int, default=600,
                    help="seconds allowed for each discoal run")
    ap.add_argument("--tsv", help="write per-test results to this file")
    args = ap.parse_args()
    if args.reps < 2:
        sys.exit("--reps must be at least 2")

    specs = [("catalog", f"{sp.id}/{mdl.id}", (sp.id, mdl.id))
             for sp in stdpopsim.all_species() for mdl in sp.demographic_models]
    specs += random_generic(np.random.default_rng(args.seed), args.generic)
    for d in args.demes_dir:
        root = os.path.dirname(os.path.abspath(d))
        for path in sorted(glob.glob(os.path.join(d, "**", "*.yaml"),
                                     recursive=True)):
            specs.append(("file", os.path.relpath(path, root), path))
    wanted = set(filter(None, args.models.split(",")))
    jobs = [(kind, label, payload, args.reps, args.samples, 1000 + 10 * i,
             args.timeout, args.rho, args.sites)
            for i, (kind, label, payload) in enumerate(specs)
            if not wanted or label in wanted]
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(audit_model, jobs))

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
    bad = {r["model"] for r in rejected} | {e["model"] for e in errors}
    tested = len(results) - len(skipped)
    print(f"\n{tested - len(bad)} of {tested} tested models pass; "
          f"{len(rejected)} tests reject; {len(errors)} models error; "
          f"{len(skipped)} skipped.")
    sys.exit(0 if not bad else 1)


if __name__ == "__main__":
    main()
