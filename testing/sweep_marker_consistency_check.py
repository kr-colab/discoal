#!/usr/bin/env python3
"""Check that diversity at a marker after a hard sweep does not depend on the
sample size.

The genealogy at a site is a property of the process, not of how many
chromosomes are sampled: the mean coalescence time of a pair at a marker is
the same whether the sample is 2 or 20 chromosomes (the mean over all pairs
of a larger sample estimates the same quantity). In a sweep that holds only if
a lineage's class (carrying the selected allele or not) changes at the rate
set by the distance between the selected site and the marker, whatever else
the lineage carries.

For a lineage whose material is only on one side of the selected site, the
class changes when a crossover (or a gene conversion tract) falls between the
site and the material, and that does not cut the material. Without those
events a lineage that has coalesced with others, and so carries more material
and still spans the site, changes class at a different rate from one that has
not, and the mean diversity at the marker depends on the sample size. At the
parameters below the difference is about 10% for crossovers and long tracts,
and about 6% for short tracts, between 2 and 20 samples.

Each case simulates the sweep with 2 and with 20 samples and compares the mean
branch diversity (branch mode, so with no mutation noise) in a one-site window
at the marker: they must agree (|z| < 4).

Locus: 10000 sites, the selected site 5% in (site 500, -x 0.05), the marker at
site 3500, so a = rho / alpha = 0.3 between them, where the difference is
largest. Hard sweep ending now with alpha = 2Ns = 100 (-a 50), N = 10000.
Gene conversion only (-r 0) with rho_total = alpha for crossovers, and with
tract lengths chosen so that the rate of tracts that separate the site from the
marker is the same: R = 2 g lambda (1 - pc^d) = 30 (4N-scaled; see
gene_conversion_check.py), with g = gamma / nSites and pc = 1 - 1/lambda.

Needs build/discoal, tskit and numpy. Fixed seeds, so the result does not
change from run to run. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
import os
import subprocess
import sys
import tempfile

import numpy as np
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N = 10_000
L = 10_000
SITE = 0.05
MARKER = 3500
CHUNK = 100
SAMPLES = (2, 20)
SWEEP = ["-t", "0", "-N", str(N), "-ws", "0", "-a", "50", "-x", str(SITE)]

# name: (recombination arguments, replicates of 2 samples, replicates of 20)
CASES = {
    "crossover (rho = 100)": (["-r", "100"], 16000, 3200),
    "gene conversion, tracts of mean 10^7": (["-r", "0", "-g", "50", "10000000"], 16000, 3200),
    "gene conversion, tracts of mean 1000": (["-r", "0", "-g", "158", "1000"], 50000, 10000),
}


def run_chunk(job):
    binary, n, args, reps, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, str(n), str(reps), str(L)] + SWEEP + args
                       + ["-d", str(seed), str(seed + 1), "-ts", prefix],
                       check=True, capture_output=True)
        return [tskit.load(os.path.join(tmp, f"d_rep{k}.trees")).diversity(
                    mode="branch", windows=[0, MARKER, MARKER + 1, L],
                    span_normalise=True)[1]
                for k in range(1, reps + 1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    parser.add_argument("--scale", type=float, default=1.0,
                        help="multiply the numbers of replicates")
    opts = parser.parse_args()

    jobs, keys = [], []
    for c, (name, (args, reps2, reps20)) in enumerate(CASES.items()):
        for i, (n, reps) in enumerate(zip(SAMPLES, (reps2, reps20))):
            reps = max(CHUNK, int(reps * opts.scale) // CHUNK * CHUNK)
            for k in range(reps // CHUNK):
                jobs.append((opts.discoal, n, args, CHUNK, 1000 + 100000 * c + 10000 * i + 13 * k))
                keys.append((name, n))

    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        results = list(pool.map(run_chunk, jobs))

    values = {}
    for key, rows in zip(keys, results):
        values.setdefault(key, []).extend(rows)

    ok = True
    for name in CASES:
        stats = []
        for n in SAMPLES:
            v = np.array(values[(name, n)])
            stats.append((v.mean(), v.std() / np.sqrt(len(v)), len(v)))
        (m2, s2, r2), (m20, s20, r20) = stats
        z = (m2 - m20) / np.sqrt(s2 ** 2 + s20 ** 2)
        good = abs(z) < 4
        ok &= good
        print(f"{'ok  ' if good else 'FAIL'} {name:38s} "
              f"n=2 {m2:.3f} (se {s2:.3f}, {r2} reps)  "
              f"n=20 {m20:.3f} (se {s20:.3f}, {r20} reps)  "
              f"difference {100 * (m2 / m20 - 1):+5.1f}%  z={z:+.2f}", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
