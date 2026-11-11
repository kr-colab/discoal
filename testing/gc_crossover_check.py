#!/usr/bin/env python3
"""Check gene conversion against crossover.

With very long tracts (-g gamma 10^7), gene conversion acts as crossover at
twice its rate. A tract that starts in a lineage's span runs past the end of
the locus, and splits the lineage at its start. A tract that starts left of
the lineage's first site and ends in its span splits it at its end: these
arrive at gamma / nSites x mean length per lineage, and end at a position of
the span with probability about 1 / mean each, so at gamma / nSites per
position too (less 0.1% here, the share of tracts shorter than the locus).
So a simulation with only such gene conversion (-r 0 -g gamma 10^7) must give
the same distribution of output as one with only crossovers (-r rho) with
gamma = rho / 2.

The check covers the choice of the lineage and of the position (both by
the lineage's span), in the neutral phase and in the sweep phase. In the
sweep phase it also covers the class of the parent that does not get the
sweep site: a random draw from the population, so it joins the lineage's
class with that class's frequency (x for the beneficial class, 1 - x for
the wild-type class).

Statistics per replicate: segregating sites, pi, distinct haplotypes, and
pi in the left, middle and right thirds of the locus. Each must agree in
mean (|z| < 4) and in distribution (two-sample KS p > 1e-4).

Needs build/discoal, numpy and scipy. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
import os
import subprocess
import sys

import numpy as np
from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N_SAMPLE = 20
SITES = 10000
RHO = 200.0
GAMMA = RHO / 2
TRACT_MEAN = "10000000"  # P(tract < 10^4 sites) = 0.1%
STATS = ["segsites", "pi", "haplotypes", "pi_left", "pi_mid", "pi_right"]
MODELS = {
    "neutral": ["-t", "40"],
    "sweep": ["-t", "40", "-ws", "0", "-a", "1000", "-x", "0.5"],
}


def summarize(output):
    rows = []
    for block in output.split("//")[1:]:
        lines = block.strip().split("\n")
        segsites = int(lines[0].split()[1])
        if segsites == 0:
            rows.append([0, 0, 1, 0, 0, 0])
            continue
        pos = np.array([float(v) for v in lines[1].split()[1:]])
        hap = np.array([[c == "1" for c in line.strip()]
                        for line in lines[2:2 + N_SAMPLE]])
        p = hap.mean(0)
        pi = 2 * p * (1 - p) * N_SAMPLE / (N_SAMPLE - 1)
        rows.append([segsites, pi.sum(),
                     len({row.tobytes() for row in hap}),
                     pi[pos < 1 / 3].sum(),
                     pi[(pos >= 1 / 3) & (pos < 2 / 3)].sum(),
                     pi[pos >= 2 / 3].sum()])
    return rows


def run_chunk(job):
    binary, args, reps, seed = job
    out = subprocess.run([binary, str(N_SAMPLE), str(reps), str(SITES)] + args
                         + ["-d", str(seed), str(seed + 1)],
                         capture_output=True, text=True, check=True).stdout
    return summarize(out)


def simulate(pool, binary, args, reps, chunks, seed):
    per = reps // chunks
    jobs = [(binary, args, per, seed + 2 * k) for k in range(chunks)]
    return np.array([row for rows in pool.map(run_chunk, jobs) for row in rows],
                    float)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=6000)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()
    chunks = max(1, min(opts.workers * 2, opts.reps // 100))

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for name, model in MODELS.items():
            crossover = simulate(pool, opts.discoal, model + ["-r", str(RHO)],
                                 opts.reps, chunks, 1000)
            conversion = simulate(pool, opts.discoal, model + ["-r", "0", "-g", str(GAMMA),
                                                 TRACT_MEAN],
                                  opts.reps, chunks, 5000)
            for j, stat in enumerate(STATS):
                a, b = conversion[:, j], crossover[:, j]
                se = np.sqrt(a.var() / len(a) + b.var() / len(b))
                z = (a.mean() - b.mean()) / se if se > 0 else 0.0
                p = stats.ks_2samp(a, b).pvalue
                good = abs(z) < 4 and p > 1e-4
                ok &= good
                print(f"{'ok  ' if good else 'FAIL'} {name:8s} {stat:10s} "
                      f"conversion {a.mean():8.3f}  crossover {b.mean():8.3f}  "
                      f"z={z:+.2f}  KS p={p:.3g}", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
