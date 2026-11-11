#!/usr/bin/env python3
"""Check hard sweeps against msprime.

discoal's hard sweep (-ws) and msprime's SweepGenicSelection model the same
process: a new allele with additive selection goes from frequency 1/(2N) to
1 - 1/(2N) along a diffusion trajectory, and lineages coalesce and recombine
between the two backgrounds. With discoal's alpha = 2Ns (s per allele copy)
msprime's s is alpha / N (fitness 1 + s for the homozygote).

For each case both simulators give 20 genomes on a locus of 10^4 sites with
the sweep site in the middle, and the statistics per replicate are, per
base pair and in generations: total branch length and diversity (branch
mode, so with no mutation noise). They must agree in mean (|z| < 4) and in
distribution (two-sample KS p > 1e-4).

discoal: -ws tau is in units of 4N generations and node times in units of
2N generations. msprime: diploid population of size N, haploid samples,
recombination rate rho / (4 N L) per bp per generation.

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
L = 10_000
GENOMES = 20
CHUNK = 50
STATS = ["branch length", "diversity"]

# name: (alpha = 2Ns, rho = 4NrL, tau = end of the sweep, in 4N generations)
CASES = {
    "alpha=1000 rho=20": (1000, 20, 0.0),
    "alpha=200 rho=20": (200, 20, 0.0),
    "alpha=2000 rho=100": (2000, 100, 0.0),
    "alpha=1000 rho=20, ended 0.01 x 4N ago": (1000, 20, 0.01),
}


def branch_stats(ts):
    return (ts.segregating_sites(mode="branch"), ts.diversity(mode="branch"))


def run_discoal(job):
    binary, alpha, rho, tau, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, str(GENOMES), str(CHUNK), str(L), "-t", "0",
                        "-r", str(rho), "-ws", str(tau), "-a", str(alpha),
                        "-x", "0.5", "-N", str(N), "-d", str(seed),
                        str(seed + 1), "-ts", prefix],
                       check=True, capture_output=True)
        # discoal node times are in units of 2N generations.
        return [tuple(2 * N * v for v in branch_stats(
                    tskit.load(os.path.join(tmp, f"d_rep{k}.trees"))))
                for k in range(1, CHUNK + 1)]


def run_msprime(job):
    alpha, rho, tau, seed = job
    sweep = msprime.SweepGenicSelection(
        position=L // 2, start_frequency=1 / (2 * N),
        end_frequency=1 - 1 / (2 * N), s=alpha / N, dt=1e-6)
    model = [sweep, msprime.StandardCoalescent()]
    if tau > 0:
        model = [msprime.StandardCoalescent(duration=4 * N * tau)] + model
    return [branch_stats(ts) for ts in msprime.sim_ancestry(
        GENOMES // 2, population_size=N, sequence_length=L,
        recombination_rate=rho / (4 * N * L), model=model,
        num_replicates=CHUNK, random_seed=seed)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()
    chunks = max(1, opts.reps // CHUNK)

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for c, (name, (alpha, rho, tau)) in enumerate(CASES.items()):
            d = np.array([row for rows in pool.map(run_discoal, [
                (opts.discoal, alpha, rho, tau, 1000 * (c + 1) + 2 * k)
                for k in range(chunks)]) for row in rows])
            m = np.array([row for rows in pool.map(run_msprime, [
                (alpha, rho, tau, 5000 * (c + 1) + k)
                for k in range(chunks)]) for row in rows])
            for j, stat in enumerate(STATS):
                a, b = d[:, j], m[:, j]
                se = np.sqrt(a.var() / len(a) + b.var() / len(b))
                z = (a.mean() - b.mean()) / se if se > 0 else 0.0
                p = stats.ks_2samp(a, b).pvalue
                good = abs(z) < 4 and p > 1e-4
                ok &= good
                print(f"{'ok  ' if good else 'FAIL'} {name:40s} {stat:13s} "
                      f"discoal {a.mean():9.4f}  msprime {b.mean():9.4f}  "
                      f"z={z:+.2f}  KS p={p:.3g}", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
