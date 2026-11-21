#!/usr/bin/env python3
"""Check diversity near a soft sweep against a forward diffusion.

A sweep from standing variation (-f f0) became beneficial when its
frequency first reached f0, forward in time, and stayed beneficial after
that, also when it fell below f0 again. Back in time, the beneficial phase
therefore ends at the last passage of the frequency through f0, not the
first. If it ends at the first, the trajectory spends too little time in
the beneficial phase near f0, and diversity near the sweep is about 1.5%
too low at the markers below.

The model: N = 1000, alpha = 100 (-a 100), a sweep that has just ended
(-ws 0), with f0 = 0.05 and 0.1. Markers to the right of the selected site
at crossover rates r = 5 and 10 from it (discoal's units, 2N generations).
For 20 samples, diversity at each marker (branch mode) is compared with the
reference of testing/soft_sweep_reference.py: the model simulated forward
in time as a diffusion (neutral until the allele first reaches f0, then
selected; lost paths start again), with the exact expected coalescence
time of two lineages on each path. The reference has 2 million paths. Each
value must agree with it (|z| < 4).

Needs build/discoal, tskit and numpy. Fixed seeds, so the result does not
change from run to run. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
CHUNK = 500
SITES = 10000
RHO = 40  # rho / (2 SITES) = 0.002 per site
MARKERS = {5.0: 2500, 10.0: 5000}
# testing/soft_sweep_reference.py: mean and standard error
REFERENCE = {
    0.05: {5.0: (0.82454, 0.00005), 10.0: (1.20985, 0.00005)},
    0.1: {5.0: (0.95100, 0.00007), 10.0: (1.33889, 0.00006)},
}


def run_chunk(job):
    binary, f0, seed = job
    positions = list(MARKERS.values())
    windows = [0]
    for p in positions:
        windows += [p, p + 1]
    windows.append(SITES)
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, "20", str(CHUNK), str(SITES), "-t", "0",
                        "-r", str(RHO), "-x", "0", "-N", "1000", "-ws", "0",
                        "-a", "100", "-f", str(f0),
                        "-d", str(seed), str(seed + 1), "-ts", prefix],
                       check=True, capture_output=True)
        rows = []
        for k in range(1, CHUNK + 1):
            d = tskit.load(os.path.join(tmp, f"d_rep{k}.trees")).diversity(
                mode="branch", windows=windows)
            rows.append([d[windows.index(p)] for p in positions])
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=80000)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for i, (f0, ref) in enumerate(REFERENCE.items()):
            values = np.array([row for rows in pool.map(run_chunk, [
                (opts.discoal, f0, 1000 + 1000000 * i + 2 * k)
                for k in range(max(1, opts.reps // CHUNK))]) for row in rows])
            for j, r in enumerate(MARKERS):
                mean = values[:, j].mean()
                se = values[:, j].std() / math.sqrt(len(values))
                target, target_se = ref[r]
                z = (mean - target) / math.hypot(se, target_se)
                good = abs(z) < 4
                ok &= good
                print(f"{'ok  ' if good else 'FAIL'} f0={f0:<5g} r={r:<3g} diversity "
                      f"{mean:.4f} +- {se:.4f}  reference {target:.4f}  "
                      f"{100 * (mean / target - 1):+.2f}%  z={z:+.2f}", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
