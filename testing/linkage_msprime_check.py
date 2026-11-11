#!/usr/bin/env python3
"""Check two-site linkage against msprime.

In Hudson's model (ms, msprime) a breakpoint anywhere in a lineage's span
splits it, also in a gap whose sites have found their most recent common
ancestor (MRCA). If such breakpoints are rejected, distant sites stay
linked for too long, while single-site statistics and the number of trees
do not change. This check measures linkage directly, for pairs of sites
100, 400 and 800 apart on a locus of 1000 sites:

  - P(same MRCA): the two sites have their MRCA in the same ancestor;
  - corr(TMRCA): the correlation of the two sites' TMRCAs.

Each must agree with msprime (|z| < 4). With the rejection, P(same MRCA)
for distant sites is 2 to 4 times msprime's (z > 10 here).

msprime: discrete genome, haploid samples from a diploid population of
size N, recombination rate rho / (4 N L) per bp. TMRCAs are compared only
within each simulator, so the time units do not matter.

Needs build/discoal, msprime, tskit and numpy. Fixed seeds, so the result
does not change from run to run. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
import glob
import os
import subprocess
import sys
import tempfile

import msprime
import numpy as np
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N = 10000
L = 1000
N_SAMPLE = 10
RHO = 20
PAIRS = [(450, 550), (300, 700), (100, 900)]
CHUNK = 500


def constant():
    d = msprime.Demography()
    d.add_population(name="a", initial_size=N)
    return d


def growth():
    # discoal -eg 0 0 2.0 -en 0.5 0 0.37: growth at rate 2 per 4N
    # generations, then (when the size is exp(-1)) a constant size.
    d = msprime.Demography()
    d.add_population(name="a", initial_size=N, growth_rate=2.0 / (4 * N))
    d.add_population_parameters_change(time=0.5 * 4 * N, population="a",
                                       initial_size=0.37 * N, growth_rate=0)
    return d


MODELS = {
    "constant size": ([], constant),
    "exponential growth": (["-eg", "0", "0", "2.0", "-en", "0.5", "0", "0.37"],
                           growth),
}


def root_times(ts):
    out = []
    for a, b in PAIRS:
        ta, tb = ts.at(a), ts.at(b)
        out += [ta.time(ta.root), tb.time(tb.root)]
    return out


def discoal_chunk(job):
    binary, args, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, str(N_SAMPLE), str(CHUNK), str(L), "-t", "0",
                        "-r", str(RHO)] + args
                       + ["-d", str(seed), str(seed + 1), "-ts", prefix],
                       check=True, capture_output=True)
        return [root_times(tskit.load(f))
                for f in glob.glob(os.path.join(tmp, "d_rep*.trees"))]


def msprime_chunk(job):
    model, seed = job
    return [root_times(ts) for ts in msprime.sim_ancestry(
        [msprime.SampleSet(N_SAMPLE, ploidy=1)], ploidy=2,
        demography=MODELS[model][1](), sequence_length=L, discrete_genome=True,
        recombination_rate=RHO / (4 * N * L), num_replicates=CHUNK,
        random_seed=seed)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=20000)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()
    chunks = max(1, opts.reps // CHUNK)

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for model, (args, _) in MODELS.items():
            d = np.array([row for rows in pool.map(
                discoal_chunk, [(opts.discoal, args, 1000 + 2 * k)
                                for k in range(chunks)]) for row in rows])
            m = np.array([row for rows in pool.map(
                msprime_chunk, [(model, 5000 + k) for k in range(chunks)])
                for row in rows])
            for k, (a, b) in enumerate(PAIRS):
                same_d = np.isclose(d[:, 2 * k], d[:, 2 * k + 1]).mean()
                same_m = np.isclose(m[:, 2 * k], m[:, 2 * k + 1]).mean()
                se = np.sqrt(same_d * (1 - same_d) / len(d)
                             + same_m * (1 - same_m) / len(m))
                z_same = (same_d - same_m) / se if se > 0 else 0.0
                r_d = np.corrcoef(d[:, 2 * k], d[:, 2 * k + 1])[0, 1]
                r_m = np.corrcoef(m[:, 2 * k], m[:, 2 * k + 1])[0, 1]
                z_corr = ((np.arctanh(r_d) - np.arctanh(r_m))
                          / np.sqrt(1 / (len(d) - 3) + 1 / (len(m) - 3)))
                good = abs(z_same) < 4 and abs(z_corr) < 4
                ok &= good
                print(f"{'ok  ' if good else 'FAIL'} {model:18s} sites {a},{b}: "
                      f"P(same MRCA) discoal {same_d:.4f} msprime {same_m:.4f} "
                      f"(z={z_same:+.1f}); corr(TMRCA) {r_d:.3f} vs {r_m:.3f} "
                      f"(z={z_corr:+.1f})", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
