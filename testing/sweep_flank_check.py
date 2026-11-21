#!/usr/bin/env python3
"""Check diversity at the flank of a hard sweep against exact values.

Three parts of the sweep model meet here (issue #164):

  - a lineage changes background (beneficial or wild type) by crossing
    over between the selected site and its material, also when the site
    lies outside the material it carries;
  - the sweep's wait for its next event ends at the last step of the
    trajectory, and no event may follow unless one was drawn;
  - at the end of its trajectory the allele is at its origin, one copy, so
    the lineages that carry it have one ancestor.

Each error makes lineages hitchhike too strongly or too weakly by a few
percent. Two checks tell them apart from the correct model:

1. Exact, with no sampling error in the reference: a deterministic sweep
   (-wd 0, alpha = 50, N = 1000) on a locus of two sites, the selected site
   at site 0 and the marker at site 1, rho = 40 (a background redraw rate of
   10 for the marker). The expected pairwise coalescence time at the marker
   comes from the three-state chain of two lineages (both beneficial, one
   of each, both wild type) on discoal's trajectory and time grid. Diversity
   at the marker must agree with it (|z| < 4).

2. The model of issue #164: a hard sweep that has just ended (-ws 0),
   alpha = 2 N_chromosomes s = 100 (discoal -a 50, -N 1000), and two markers
   to the right of a selected site at the left end of the locus, at map
   distances a1 = 0.1 and a3 = 0.2 in units of s. For 20 samples the
   statistics are ratios, sweep over neutral, of branch-mode expectations:
   diversity at each marker (h1, h3) and pi2 between them. The targets come
   from an exact solver of the two-locus dynamics under the conditioned
   trajectory (issue #164, N. Pope); forward simulations in SLiM agree with
   it. The neutral expectation of diversity is exactly 2; that of pi2 is
   the mean of 1.44 million neutral replicates (discoal, and a separate
   two-locus simulator). Each ratio must agree with its target (|z| < 4).

Needs build/discoal, tskit, numpy and scipy. Fixed seeds, so the result
does not change from run to run. Exit code 0 if every check passes.
"""

import argparse
import concurrent.futures
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.linalg
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
CHUNK = 500

# Check 1
ALPHA_DET = 50.0
N = 1000
DELTA_T_MOD = 40  # discoal's default -i
TWO_SITE_RHO = 40
REDRAW_RATE = TWO_SITE_RHO * 0.5 / 2  # rho / 2 per locus, one of two sites

# Check 2: rho over the locus is alpha x a per map distance a, so 20 for
# a3 = 0.2 at the far end; a1 = 0.1 is at half the locus.
SITES = 10000
RHO = 20
MARKERS = (SITES // 2, SITES - 2)
TARGETS = {"h1": 0.5587, "h3": 0.7492, "pi2": 0.4436}
NEUTRAL = {"h1": (2.0, 0.0), "h3": (2.0, 0.0), "pi2": (1.0149, 0.0007)}


def det_freq(t):
    """discoal's deterministic trajectory (detSweepFreq)."""
    eps = 0.05 / ALPHA_DET
    ts = -2.0 * math.log(eps) / ALPHA_DET
    return eps / (eps + (1.0 - eps) * math.exp(ALPHA_DET * (t - ts)))


def exact_pair_time():
    """E[T] for two marker lineages: a three-state chain on discoal's grid
    (x from step k = 1 until x <= 1/(2N), all samples beneficial at the
    start; at the origin two beneficial lineages are one, and the others
    go on in the neutral coalescent, one unit of time on average)."""
    dt = 1.0 / (DELTA_T_MOD * N)
    r = REDRAW_RATE
    p = np.array([1.0, 0.0, 0.0])  # both B, one of each, both b
    expected = 0.0
    k, x = 0, 1.0
    while x > 1.0 / (2 * N):
        k += 1
        x = det_freq(k * dt)
        q = np.array([[-(1 / x + 2 * r * (1 - x)), 2 * r * (1 - x), 0.0],
                      [r * x, -r, r * (1 - x)],
                      [0.0, 2 * r * x, -(1 / (1 - x) + 2 * r * x)]])
        # The step's survival integral with the state: one exponential of
        # the generator extended by an accumulator.
        m = np.zeros((4, 4))
        m[:3, :3] = q.T
        m[3, :3] = 1.0
        v = scipy.linalg.expm(m * dt) @ np.append(p, 0.0)
        expected += v[3]
        p = v[:3]
    return expected + p[1] + p[2]


def run_two_site(job):
    binary, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, "20", str(CHUNK), "2", "-t", "0",
                        "-r", str(TWO_SITE_RHO), "-wd", "0", "-a", str(ALPHA_DET),
                        "-x", "0", "-N", str(N), "-d", str(seed), str(seed + 1),
                        "-ts", prefix], check=True, capture_output=True)
        return [tskit.load(os.path.join(tmp, f"d_rep{k}.trees")).diversity(
                    mode="branch", windows=[0, 1, 2])[1]
                for k in range(1, CHUNK + 1)]


def flank_stats(ts):
    x1, x3 = MARKERS
    # Only the trees at the two markers matter; the rest of the locus is
    # cut away first, which makes the statistics much faster.
    ts = ts.keep_intervals([[x1, x1 + 1], [x3, x3 + 1]], simplify=False)
    h = [ts.diversity(mode="branch", windows=[0, x, x + 1, SITES])[1]
         for x in (x1, x3)]
    pi2 = ts.ld_matrix(mode="branch", stat="pi2_unbiased",
                       positions=[[x1], [x3]])[0, 0]
    return [h[0], h[1], pi2]


def run_flank(job):
    binary, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, "20", str(CHUNK), str(SITES), "-t", "0",
                        "-r", str(RHO), "-N", str(N), "-ws", "0", "-a", "50",
                        "-x", "0", "-d", str(seed), str(seed + 1),
                        "-ts", prefix], check=True, capture_output=True)
        return [flank_stats(tskit.load(os.path.join(tmp, f"d_rep{k}.trees")))
                for k in range(1, CHUNK + 1)]


def report(name, value, se, target, extra=""):
    z = (value - target) / se
    good = abs(z) < 4
    print(f"{'ok  ' if good else 'FAIL'} {name:20s} {value:.4f} +- {se:.4f}  "
          f"exact {target:.4f}  z={z:+.2f}{extra}", flush=True)
    return good


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reps", type=int, default=120000,
                        help="replicates of the issue #164 model")
    parser.add_argument("--two-site-reps", type=int, default=40000)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    opts = parser.parse_args()

    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        two = np.array([v for vs in pool.map(run_two_site, [
            (opts.discoal, 1000 + 2 * k)
            for k in range(max(1, opts.two_site_reps // CHUNK))]) for v in vs])
        flank = np.array([row for rows in pool.map(run_flank, [
            (opts.discoal, 500000 + 2 * k)
            for k in range(max(1, opts.reps // CHUNK))]) for row in rows])

    ok = report("two sites, -wd", two.mean(), two.std() / math.sqrt(len(two)),
                2.0 * exact_pair_time())
    for j, (name, target) in enumerate(TARGETS.items()):
        mean = flank[:, j].mean()
        se_mean = flank[:, j].std() / math.sqrt(len(flank))
        ref, ref_se = NEUTRAL[name]
        ratio = mean / ref
        se = ratio * math.hypot(se_mean / mean, ref_se / ref)
        ok &= report(f"{name} sweep/neutral", ratio, se, target, "  (solver)")
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
