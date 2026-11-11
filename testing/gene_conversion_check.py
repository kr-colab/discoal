#!/usr/bin/env python3
"""Check gene conversion against theory and against msprime.

discoal's gene conversion is the model of ms and msprime (Wiuf and Hein
2000): tracts start at rate g per site and run to the right, with a length
that is geometric on 1, 2, 3, ... with mean lam. The chromosome goes on to
the left of the locus, so a tract can also start left of a lineage's first
site and run into its material. A tract that changes no lineage's material
(it falls in a gap, or covers all of the material) is not an event.

1. Two-locus theory. With two samples, an event that separates two sites
   acts as a recombination between them, so the correlation of the sites'
   TMRCAs is (R + 18) / (R^2 + 13 R + 18), where R = 2 g lam (1 - pc^d)
   (4N-scaled, pc = 1 - 1/lam) is the rate of tracts that cover one of two
   sites d apart and not the other. R does not depend on where the pair is,
   so pairs at the left end, in the middle and at the right end must all
   agree with it (|z| < 4). Without tracts that come in from the left, the
   pairs at the left end are too correlated (z > 15 here).

2. msprime. Ten samples, gene conversion and crossover: the number of
   trees and the TMRCA at three sites must agree in mean (|z| < 4) and in
   distribution (two-sample KS p > 1e-4).

Units: discoal -g gamma is 4 N g nSites (as -r rho is 4 N r nSites), and
node times are in units of 2N generations. msprime: discrete genome of
nSites positions, haploid samples from a diploid population of size N.

Needs build/discoal, msprime, tskit, numpy and scipy. Fixed seeds, so the
result does not change from run to run. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
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
N = 10_000
CHUNK = 500

# Two-locus part: 1001 sites, gene conversion only.
TL_SITES = 1001
TL_CASES = {  # name: (lam, g per site, 4N-scaled)
    "mean tract 100": (100, 0.025),
    "mean tract 10": (10, 0.25),
}
TL_PAIRS = []
for _d in (5, 20, 100):
    TL_PAIRS += [("left", 0, _d), ("middle", 300, 300 + _d),
                 ("right", TL_SITES - 1 - _d, TL_SITES - 1)]
TL_SITES_USED = sorted({p[1] for p in TL_PAIRS} | {p[2] for p in TL_PAIRS})

# msprime part: ten samples, gene conversion and crossover.
MS_SAMPLES = 10
MS_SITES = 1000
MS_RHO = 10.0
MS_GAMMA = 20.0
MS_LAM = 50
MS_AT = [0, 500, MS_SITES - 1]
MS_STATS = ["trees"] + [f"TMRCA site {s}" for s in MS_AT]


def run_discoal(binary, n, reps, sites, args, seed, tmp):
    prefix = os.path.join(tmp, "d.trees")
    subprocess.run([binary, str(n), str(reps), str(sites), "-t", "0"] + args
                   + ["-d", str(seed), str(seed + 1), "-ts", prefix],
                   check=True, capture_output=True)
    return [tskit.load(os.path.join(tmp, f"d_rep{k}.trees"))
            for k in range(1, reps + 1)]


def root_time(ts, site):
    tree = ts.at(site)
    return tree.time(tree.root)


def two_locus_chunk(job):
    binary, lam, g, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        return [[root_time(ts, s) for s in TL_SITES_USED]
                for ts in run_discoal(binary, 2, CHUNK, TL_SITES,
                                      ["-r", "0", "-g", repr(g * TL_SITES),
                                       str(lam)], seed, tmp)]


def ms_stats(ts, scale):
    return [ts.num_trees] + [scale * root_time(ts, s) for s in MS_AT]


def discoal_chunk(job):
    binary, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        # discoal node times are in units of 2N generations.
        return [ms_stats(ts, 2 * N) for ts in run_discoal(
            binary, MS_SAMPLES, CHUNK, MS_SITES,
            ["-r", str(MS_RHO), "-g", str(MS_GAMMA), str(MS_LAM)], seed, tmp)]


def msprime_chunk(seed):
    return [ms_stats(ts, 1) for ts in msprime.sim_ancestry(
        [msprime.SampleSet(MS_SAMPLES, ploidy=1)], ploidy=2,
        population_size=N, sequence_length=MS_SITES, discrete_genome=True,
        recombination_rate=MS_RHO / (4 * N * MS_SITES),
        gene_conversion_rate=MS_GAMMA / (4 * N * MS_SITES),
        gene_conversion_tract_length=MS_LAM,
        num_replicates=CHUNK, random_seed=seed)]


def jackknife_corr(x, y, groups):
    """Correlation of x and y, and its standard error by the jackknife over
    groups of rows. TMRCAs are far from normal, so the normal-theory
    standard error of a correlation does not apply."""
    r = np.corrcoef(x, y)[0, 1]
    parts = np.array_split(np.arange(len(x)), groups)
    loo = np.array([np.corrcoef(np.delete(x, idx), np.delete(y, idx))[0, 1]
                    for idx in parts])
    se = np.sqrt((groups - 1) / groups * ((loo - loo.mean()) ** 2).sum())
    return r, se


def theory(lam, g, d):
    pc = 1 - 1 / lam
    r = 2 * g * lam * (1 - pc ** d)
    return (r + 18) / (r * r + 13 * r + 18)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=20000,
                        help="replicates of each two-locus case")
    parser.add_argument("--ms-reps", type=int, default=5000,
                        help="replicates of the msprime comparison")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for c, (name, (lam, g)) in enumerate(TL_CASES.items()):
            t = np.array([row for rows in pool.map(two_locus_chunk, [
                (opts.discoal, lam, g, 1000 * (c + 1) + 2 * k)
                for k in range(max(1, opts.reps // CHUNK))]) for row in rows])
            for where, i, j in TL_PAIRS:
                expected = theory(lam, g, j - i)
                r, se = jackknife_corr(t[:, TL_SITES_USED.index(i)],
                                       t[:, TL_SITES_USED.index(j)], 40)
                z = (r - expected) / se
                good = abs(z) < 4
                ok &= good
                print(f"{'ok  ' if good else 'FAIL'} two-locus, {name:14s} "
                      f"{where:6s} sites {i:4d},{j:4d}: corr(TMRCA) {r:.4f} "
                      f"theory {expected:.4f}  z={z:+.2f}", flush=True)

        chunks = max(1, opts.ms_reps // CHUNK)
        d = np.array([row for rows in pool.map(discoal_chunk, [
            (opts.discoal, 9000 + 2 * k) for k in range(chunks)])
            for row in rows])
        m = np.array([row for rows in pool.map(msprime_chunk, [
            9500 + k for k in range(chunks)]) for row in rows])
        for k, stat in enumerate(MS_STATS):
            a, b = d[:, k], m[:, k]
            se = np.sqrt(a.var() / len(a) + b.var() / len(b))
            z = (a.mean() - b.mean()) / se if se > 0 else 0.0
            p = stats.ks_2samp(a, b).pvalue
            good = abs(z) < 4 and p > 1e-4
            ok &= good
            print(f"{'ok  ' if good else 'FAIL'} msprime, {stat:15s} "
                  f"discoal {a.mean():10.2f}  msprime {b.mean():10.2f}  "
                  f"z={z:+.2f}  KS p={p:.3g}", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
