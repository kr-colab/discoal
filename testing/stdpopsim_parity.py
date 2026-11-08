#!/usr/bin/env python3
"""Compare discoal and msprime on every demographic model in stdpopsim.

For each (species, model) in the stdpopsim catalog:

  1. Export the model to a demes file.
  2. Simulate `reps` genealogies with discoal (through a -Y config that
     loads the demes file, with tree sequence output) and with msprime
     (Demography.from_demes), with no recombination.
  3. Compute branch-mode statistics on each tree: diversity, segregating
     sites (total branch length) and Tajima's D over all samples, and
     diversity within each sampled deme.
  4. Compare the two distributions of each statistic with a two-sample KS
     test, Bonferroni-corrected over all tests in the run.

The comparison is on absolute scale. discoal writes node times in units of
2 N_ref generations, where N_ref is the present-day size of the first
present-day deme; they are converted to generations before the statistics
are computed. A time-unit or size-scaling error therefore shows up as a
shift in diversity and segregating sites, not only in Tajima's D.

Samples: `--samples` haploid genomes from each deme that exists at time 0
and that stdpopsim samples at time 0. discoal's deme_sample_size is indexed
by the order of demes in the graph, extinct demes included, so demes that
are not sampled get 0.

Usage:
    python testing/stdpopsim_parity.py [--reps 1000] [--samples 10]
        [--workers 24] [--models HomSap/OutOfAfrica_3G09,...] [--tsv out.tsv]

Needs build/discoal, stdpopsim, msprime, demes, tskit and scipy.
Exit code 0 if no test rejects, 1 otherwise.
"""

import argparse
import concurrent.futures
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


def sampled_demes(model, graph):
    """Names of demes to sample: present at time 0 and sampled at time 0 by
    stdpopsim."""
    sampling_time = {p.name: p.default_sampling_time
                     for p in model.populations}
    return [d.name for d in graph.demes
            if d.end_time == 0 and sampling_time.get(d.name) in (None, 0)]


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
        if len(sets) > 1:
            for name, s in zip(names, sets):
                vals[f"diversity[{name}]"] = (
                    ts.diversity(s, mode="branch") * time_scale)
        for k, v in vals.items():
            out.setdefault(k, []).append(float(v))
    return {k: np.array(v) for k, v in out.items()}


def run_discoal(demes_path, graph, names, n, reps, seed, tmp, timeout):
    """Run discoal; return ({stat: array}, N_ref)."""
    index = {d.name: i for i, d in enumerate(graph.demes)}
    sizes = [n if d.name in names else 0 for d in graph.demes]
    cfg = os.path.join(tmp, "config.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {sum(sizes)}, "
                f"num_replicates: {reps}, num_sites: 1000, "
                f"seed: [{seed}, {seed + 1}]}}\n"
                "genetics: {mutation_rate: 1.0, recombination_rate: 0}\n"
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


def run_msprime(graph, names, n, reps, seed):
    dem = msprime.Demography.from_demes(graph)
    samples = [msprime.SampleSet(n, population=name, ploidy=1, time=0)
               for name in names]
    pop_ids = [dem[name].id for name in names]
    tss = msprime.sim_ancestry(samples=samples, demography=dem, ploidy=2,
                               num_replicates=reps, random_seed=seed)
    return collect(tss, pop_ids, names, 1.0)


def audit_model(job):
    """Run one model. Returns a dict with rows of test results or an error."""
    species, model_id, reps, n, seed, timeout = job
    label = f"{species}/{model_id}"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = stdpopsim.get_species(species).get_demographic_model(
                model_id)
            graph = model.model.to_demes()
        names = sampled_demes(model, graph)
        if not names:
            return dict(model=label, skipped="no deme is sampled at time 0 "
                        "(discoal YAML has no ancient samples)")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "model.demes.yaml")
            demes.dump(graph, path)
            d, n_ref = run_discoal(path, graph, names, n, reps, seed, tmp,
                                   timeout)
        m = run_msprime(graph, names, n, reps, seed)
        rows = []
        for stat in d:
            ks = stats.ks_2samp(d[stat], m[stat])
            rows.append(dict(model=label, stat=stat,
                             discoal_mean=d[stat].mean(),
                             msprime_mean=m[stat].mean(),
                             ks_d=ks.statistic, p=ks.pvalue))
        return dict(model=label, n_ref=n_ref, demes=names, rows=rows)
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
    ap.add_argument("--models", default="",
                    help="comma-separated Species/Model ids; default all")
    ap.add_argument("--timeout", type=int, default=600,
                    help="seconds allowed for each discoal run")
    ap.add_argument("--tsv", help="write per-test results to this file")
    args = ap.parse_args()
    if args.reps < 2:
        sys.exit("--reps must be at least 2")

    jobs = []
    wanted = set(filter(None, args.models.split(",")))
    for sp in stdpopsim.all_species():
        for mdl in sp.demographic_models:
            if not wanted or f"{sp.id}/{mdl.id}" in wanted:
                jobs.append((sp.id, mdl.id, args.reps, args.samples,
                             1000 + 10 * len(jobs), args.timeout))
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
