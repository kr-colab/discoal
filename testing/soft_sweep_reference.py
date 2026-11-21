#!/usr/bin/env python3
"""Reference values for testing/soft_sweep_check.py (not run in CI).

The model is the one of testing/slim_sweep_check.py, as a diffusion,
simulated forward in time: the allele starts as one copy (x = 1/(2N)), is
neutral until its frequency first reaches f0, and then has the selection
drift alpha x (1 - x) (discoal's alpha, its -a); a path that is lost
(x <= 0) starts again, and a path ends when x reaches 1 - 1/(2N), the time
of sampling (-ws 0). Steps of discoal's size, 1/(40 N), with Gaussian
increments. No conditioning is used: the restarts give it.

On each path, back in time, the expected coalescence time of two lineages
at a marker comes from the chain of their backgrounds (both beneficial,
one of each, both wild type): a lineage changes background at rate
r (1 - x) or r x, where r is the crossover rate between the selected site
and the marker (discoal's units), and two lineages on one background
coalesce at rate 1/x or 1/(1 - x). At the origin the beneficial lineages
are one; the others go on in the neutral coalescent (one unit of time on
average). Diversity at the marker is 2 E[T].

Needs numba. Prints the mean and its standard error for each case.
"""

import argparse
import math

import numba as nb
import numpy as np

N = 1000
DT = 1.0 / (40 * N)
MAXS = 400000


@nb.njit(cache=True)
def forward_path(alpha, f0, buf):
    """A path that fixes, stored back in time (buf[0] at the end)."""
    tmp = np.empty(MAXS)
    while True:
        x = 1.0 / (2 * N)
        selected = f0 <= 0.0
        k = 0
        ok = True
        while x < 1.0 - 1.0 / (2 * N):
            if not selected and x >= f0:
                selected = True
            tmp[k] = x
            k += 1
            if k >= MAXS:
                ok = False
                break
            pq = x * (1 - x)
            x += ((alpha * pq) if selected else 0.0) * DT \
                + math.sqrt(pq * DT) * np.random.normal()
            if x <= 0.0:
                ok = False
                break
        if ok:
            for i in range(k):
                buf[i] = tmp[k - 1 - i]
            return k


@nb.njit(cache=True)
def pair_time(buf, length, r):
    pBB, pBb, pbb = 1.0, 0.0, 0.0
    t = 0.0
    for i in range(length):
        x = buf[i]
        t += (pBB + pBb + pbb) * DT
        toB_b = 1.0 - math.exp(-r * (1 - x) * DT)
        tob_B = 1.0 - math.exp(-r * x * DT)
        nBB = pBB * (1 - toB_b) ** 2 + pBb * tob_B * (1 - toB_b) + pbb * tob_B ** 2
        nbb = pbb * (1 - tob_B) ** 2 + pBb * toB_b * (1 - tob_B) + pBB * toB_b ** 2
        nBb = (pBB + pBb + pbb) - nBB - nbb
        pBB = nBB * math.exp(-DT / x)
        pbb = nbb * math.exp(-DT / (1 - x))
        pBb = nBb
    return t + pBb + pbb


@nb.njit(parallel=True, cache=True)
def run(reps, alpha, f0, rs, seed):
    out = np.zeros((reps, len(rs)))
    for q in nb.prange(reps):
        np.random.seed(seed + q)
        buf = np.zeros(MAXS)
        length = forward_path(alpha, f0, buf)
        for j in range(len(rs)):
            out[q, j] = 2.0 * pair_time(buf, length, rs[j])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--paths", type=int, default=2000000)
    opts = parser.parse_args()
    rs = np.array([5.0, 10.0])
    for f0 in (0.0, 0.05, 0.1):
        a = run(opts.paths, 100.0, f0, rs, 20261121)
        mean = a.mean(0)
        se = a.std(0) / math.sqrt(len(a))
        print(f"f0={f0}: " + "  ".join(
            f"r={r:g}: {m:.5f} +- {s:.5f}" for r, m, s in zip(rs, mean, se)),
            flush=True)


if __name__ == "__main__":
    main()
