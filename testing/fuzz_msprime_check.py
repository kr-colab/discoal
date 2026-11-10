#!/usr/bin/env python3
"""Compare discoal with msprime on random command-line demographic models.

Each model is drawn at random from a seed: 1 to 3 populations, their
sizes, size changes (-en), exponential growth (-eg) that a size change
ends, constant and changing migration (-m, -em), an admixture (-ea, with
3 populations), ancient samples (-A), and recombination in some models.
All populations finally merge (-ed), and a merged or admixed population
stops all its migration (-em ... 0), so no model can strand lineages.
The same model is built in msprime, with the unit conventions that
testing/units_convention_check.py pins: times in 4N0 generations, -eg
alpha per 4N0 generations, migration rates 4N0 m with -m i j moving
lineages from i to j backward in time.

For each model, both simulators give the branch-mode diversity of every
sample set (each population at time 0, and each ancient sample), the
divergence of every pair of sets, and the total branch length, in
generations. Each statistic is compared with a two-sample KS test; a
model fails if any p-value is below 0.001 / (number of tests of the run)
(Bonferroni). The output lists each model, its command line, and its
smallest p-value.

Usage:
    python testing/fuzz_msprime_check.py [--models 30] [--reps 300]
        [--seed 1] [--workers 8] [--discoal PATH]
Needs build/discoal, msprime, tskit, numpy and scipy.
"""

import argparse
import concurrent.futures
import glob
import itertools
import os
import subprocess
import sys
import tempfile

import msprime
import numpy as np
import tskit
from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N0 = 10000
L = 1000
ALPHA = 0.001


def draw_model(seed):
    """A random model: discoal arguments (after the sample size, number of
    replicates and sites) and a function that builds the msprime
    Demography and sample sets. Times are in 4N0 generations."""
    rng = np.random.default_rng(seed)
    npops = int(rng.choice([1, 2, 3], p=[0.25, 0.4, 0.35]))
    n_present = [int(rng.integers(2, 7)) for _ in range(npops)]
    sizes = [float(np.exp(rng.uniform(np.log(0.2), np.log(5.0))))
             for _ in range(npops)]
    rho = float(rng.choice([0.0, 0.0, 10.0]))

    # Each population's end: an admixture (population 2, with 3
    # populations) or a merge; the last population (0) never ends.
    end = {0: np.inf}
    admix = None
    if npops == 3 and rng.random() < 0.5:
        t_admix = float(rng.uniform(0.3, 0.9))
        admix = (t_admix, 2, 0, 1, float(rng.uniform(0.1, 0.9)))
        end[2] = t_admix
    merge_times = sorted(rng.uniform(1.0, 2.0, size=npops - 1))
    merges = []
    alive = [p for p in range(npops) if p == 0 or p not in end]
    for t in merge_times:
        if len(alive) < 2:
            break
        src = alive[-1]
        dst = alive[int(rng.integers(0, len(alive) - 1))]
        merges.append((float(t), src, dst))
        end[src] = float(t)
        alive.remove(src)

    events = []          # (time, discoal args, msprime adder)

    def add(t, args, fn):
        events.append((t, args, fn))

    # Size changes and growth that ends with a size change.
    for p in range(npops):
        for _ in range(int(rng.integers(0, 3))):
            t = float(rng.uniform(0.01, min(end[p], 1.0)))
            s = float(np.exp(rng.uniform(np.log(0.2), np.log(5.0))))
            add(t, ["-en", f"{t!r}", str(p), f"{s!r}"],
                lambda d, t=t, p=p, s=s: d.add_population_parameters_change(
                    time=t * 4 * N0, initial_size=s * N0, growth_rate=0,
                    population=p))
        if rng.random() < 0.4:
            t1 = float(rng.uniform(0.0, min(end[p], 1.0) * 0.6))
            t2 = float(rng.uniform(t1 + 0.05, min(end[p], 1.0)))
            a = float(rng.uniform(-5.0, 5.0))
            s = float(np.exp(rng.uniform(np.log(0.2), np.log(5.0))))
            add(t1, ["-eg", f"{t1!r}", str(p), f"{a!r}"],
                lambda d, t=t1, p=p, a=a: d.add_population_parameters_change(
                    time=t * 4 * N0, growth_rate=a / (4 * N0), population=p))
            add(t2, ["-en", f"{t2!r}", str(p), f"{s!r}"],
                lambda d, t=t2, p=p, s=s: d.add_population_parameters_change(
                    time=t * 4 * N0, initial_size=s * N0, growth_rate=0,
                    population=p))

    # Migration: constant rates and changes while both populations exist.
    mig0 = {}
    for i, j in itertools.permutations(range(npops), 2):
        if rng.random() < 0.5:
            mig0[(i, j)] = float(rng.uniform(0.1, 5.0))
        if rng.random() < 0.3:
            t = float(rng.uniform(0.01, min(end[i], end[j], 1.0)))
            r = float(rng.choice([0.0, rng.uniform(0.1, 5.0)]))
            add(t, ["-em", f"{t!r}", str(i), str(j), f"{r!r}"],
                lambda d, t=t, i=i, j=j, r=r: d.add_migration_rate_change(
                    time=t * 4 * N0, rate=r / (4 * N0), source=i, dest=j))

    def stop_migration(t, p):
        for q in range(npops):
            if q == p:
                continue
            for i, j in ((p, q), (q, p)):
                add(t, ["-em", f"{t!r}", str(i), str(j), "0"],
                    lambda d, t=t, i=i, j=j: d.add_migration_rate_change(
                        time=t * 4 * N0, rate=0, source=i, dest=j))

    if admix:
        t, dau, f1, f2, prop = admix
        def admix_fn(d, t=t, dau=dau, f1=f1, f2=f2, prop=prop):
            d.add_mass_migration(time=t * 4 * N0, source=dau, dest=f1,
                                 proportion=prop)
            d.add_mass_migration(time=t * 4 * N0, source=dau, dest=f2,
                                 proportion=1.0)
        add(t, ["-ea", f"{t!r}", str(dau), str(f1), str(f2),
                f"{prop!r}"], admix_fn)
        stop_migration(t, dau)
    for t, src, dst in merges:
        add(t, ["-ed", f"{t!r}", str(src), str(dst)],
            lambda d, t=t, src=src, dst=dst: d.add_mass_migration(
                time=t * 4 * N0, source=src, dest=dst, proportion=1.0))
        stop_migration(t, src)

    # Ancient samples, taken while their population exists.
    ancient = []
    for p in range(npops):
        if rng.random() < 0.3:
            t = float(rng.uniform(0.02, min(end[p], 1.5) * 0.9))
            ancient.append((p, int(rng.integers(2, 4)), t))

    events.sort(key=lambda e: e[0])   # stable: same-time order is kept
    n_total = sum(n_present) + sum(k for _, k, _ in ancient)
    args = ["-t", "0", "-r", f"{rho:g}"]
    if npops > 1:
        args += ["-p", str(npops)] + [
            str(n_present[p] + sum(k for q, k, _ in ancient if q == p))
            for p in range(npops)]
    for (i, j), r in mig0.items():
        args += ["-m", str(i), str(j), f"{r!r}"]
    for p in range(npops):
        args += ["-en", "0", str(p), f"{sizes[p]!r}"]
    for p, k, t in ancient:
        args += ["-A", str(k), str(p), f"{t!r}"]
    for _, a, _ in events:
        args += a

    def build():
        d = msprime.Demography()
        for p in range(npops):
            d.add_population(name=f"p{p}", initial_size=sizes[p] * N0)
        for (i, j), r in mig0.items():
            d.set_migration_rate(source=i, dest=j, rate=r / (4 * N0))
        for _, _, fn in events:
            fn(d)
        d.sort_events()
        sets = [msprime.SampleSet(n_present[p], population=p, ploidy=1)
                for p in range(npops)]
        sets += [msprime.SampleSet(k, population=p, time=t * 4 * N0,
                                   ploidy=1) for p, k, t in ancient]
        return d, sets

    set_keys = [(p, 0.0) for p in range(npops)] + [
        (p, t * 4 * N0) for p, _, t in ancient]
    return dict(n=n_total, args=args, build=build, rho=rho,
                sets=set_keys, sizes=[n_present[p] for p in range(npops)]
                + [k for _, k, _ in ancient])


def sample_sets(ts, keys, sizes):
    """The sample nodes of each (population, time in generations) set."""
    times = ts.nodes_time
    out = []
    for (p, t), k in zip(keys, sizes):
        s = ts.samples(population=p)
        s = s[np.isclose(times[s], t, rtol=1e-6, atol=1e-6)]
        if len(s) != k:
            raise RuntimeError(f"set {p}@{t}: {len(s)} samples, not {k}")
        out.append(s)
    return out


def statistics(ts, keys, sizes):
    sets = sample_sets(ts, keys, sizes)
    v = [ts.diversity(s, mode="branch") for s in sets]
    if len(sets) > 1:
        v += list(ts.divergence(
            sets, indexes=list(itertools.combinations(range(len(sets)), 2)),
            mode="branch"))
    v.append(ts.segregating_sites(np.concatenate(sets), mode="branch"))
    return v


def run_model(job):
    """Runs one model; returns (seed, discoal arguments, (KS p-values,
    discoal/msprime mean ratios) or None, error text)."""
    seed, reps, discoal = job
    m = draw_model(seed)
    try:
        return compare(seed, reps, discoal, m)
    except (RuntimeError, subprocess.TimeoutExpired) as e:
        return seed, m["args"], None, str(e)[-300:]


def compare(seed, reps, discoal, m):
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        cmd = [discoal, str(m["n"]), str(reps), str(L)] + m["args"] + [
            "-d", str(seed), str(seed + 1), "-ts", prefix]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
        if r.returncode != 0:
            return seed, m["args"], None, r.stderr.strip()[-300:]
        files = sorted(glob.glob(os.path.join(tmp, "d_rep*.trees")))
        # discoal node times are in units of 2N0 generations.
        d = np.array([statistics(tskit.load(f), [
            (p, t / (2 * N0)) for p, t in m["sets"]], m["sizes"])
            for f in files]) * 2 * N0
    dem, sets = m["build"]()
    ms = np.array([statistics(ts, m["sets"], m["sizes"])
                   for ts in msprime.sim_ancestry(
                       samples=sets, demography=dem, ploidy=2,
                       sequence_length=L, discrete_genome=True,
                       recombination_rate=m["rho"] / (4 * N0 * L),
                       num_replicates=reps, random_seed=seed)])
    pvals = [stats.ks_2samp(d[:, k], ms[:, k]).pvalue
             for k in range(d.shape[1])]
    ratios = d.mean(0) / ms.mean(0)
    return seed, m["args"], (pvals, ratios), ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--models", type=int, default=30)
    ap.add_argument("--reps", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--discoal", default=DISCOAL,
                    help="the discoal binary (default: build/discoal)")
    args = ap.parse_args()
    seeds = [args.seed * 100000 + k for k in range(args.models)]
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(run_model, [(s, args.reps, args.discoal)
                                            for s in seeds]))
    n_tests = sum(len(r[2][0]) for r in results if r[2])
    threshold = ALPHA / max(n_tests, 1)
    failed = 0
    for seed, margs, res, err in results:
        line = " ".join(margs)
        if res is None:
            failed += 1
            print(f"ERROR model {seed}: discoal {line}\n      {err}")
            continue
        pvals, ratios = res
        worst = int(np.argmin(pvals))
        ok = min(pvals) >= threshold
        failed += not ok
        print(f"{'PASS' if ok else 'FAIL'}  model {seed}: min KS p "
              f"{min(pvals):.2g} (statistic {worst}, discoal/msprime mean "
              f"{ratios[worst]:.3f})  discoal ... {line}", flush=True)
    print(f"\n{len(results) - failed} of {len(results)} models pass; "
          f"{n_tests} tests, Bonferroni threshold {threshold:.2g}.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
