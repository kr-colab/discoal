#!/usr/bin/env python3
"""Check the selected SNP of partial sweeps (-c) in discoal's output.

At a partial sweep that ends at the present (-ws 0 ... -c c), each sampled
lineage is in the beneficial class with probability c, so the number of
samples that carry the selected allele is Binomial(n, c), and the selected
site is at the sweep position given by -x. The checks:

  1. the selected site is at -x, and its derived count has the binomial
     mean and variance (two final frequencies);
  2. -h hides the selected site;
  3. recurrent partial sweeps (-R ... -c) leave one site per sweep, so a
     replicate can hold more than one (with a tiny theta, every site is a
     selected site);
  4. in tree sequence output (-ts) the site is at -x times the number of
     sites;
  5. ancient samples taken during the sweep carry the allele with the
     frequency it has then (here a few generations back, so about c).

Needs build/discoal, numpy and tskit. Exit code 0 if every check passes.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N_SAMPLES = 20
SITES = 10000
Z_MAX = 4.0
# Population size during sweeps (-N). The sweep phase takes time in
# proportion to it; the checks do not depend on it, and 1e4 keeps the
# script at seconds rather than minutes (the default is 1e6).
N_SWEEP = 10000


def run(args):
    return subprocess.run([DISCOAL] + [str(a) for a in args],
                          capture_output=True, text=True, check=True).stdout


def replicates(out, n_samples=N_SAMPLES):
    """(positions, haplotype rows) for each replicate of ms output."""
    reps = []
    for block in out.split("//")[1:]:
        lines = block.strip().splitlines()
        nseg = int(lines[0].split()[1])
        if nseg == 0:
            reps.append(([], []))
            continue
        positions = [float(x) for x in lines[1].split()[1:]]
        reps.append((positions, lines[2:2 + n_samples]))
    return reps


def derived_count_at(rep, x):
    """Samples carrying the derived allele at position x, or None if the
    replicate has no site there."""
    positions, rows = rep
    for k, p in enumerate(positions):
        if abs(p - x) < 1e-6:
            return sum(row[k] == "1" for row in rows)
    return None


def main():
    ok = True

    def report(passed, text):
        nonlocal ok
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'}  {text}", flush=True)

    reps = 4000
    x = 0.3
    for c in (0.5, 0.2):
        out = run([N_SAMPLES, reps, SITES, "-t", 0.0001, "-r", 20, "-ws", 0,
                   "-a", 500, "-N", N_SWEEP, "-x", x, "-c", c, "-d", 11, 12])
        found = [derived_count_at(r, x) for r in replicates(out)]
        # No site means that no sampled lineage carries the allele.
        counts = np.array([0 if k is None else k for k in found])
        mean, var = N_SAMPLES * c, N_SAMPLES * c * (1 - c)
        z = (counts.mean() - mean) / np.sqrt(var / len(counts))
        ratio = counts.var() / var
        report(abs(z) < Z_MAX and 0.85 < ratio < 1.15
               and all(k is None or k > 0 for k in found),
               f"-c {c}: derived count mean {counts.mean():.2f} (expected "
               f"{mean:.1f}, z={z:+.1f}), variance ratio {ratio:.3f}; every "
               f"carried site at x={x}")

    # A partial sweep that ended in the past, with much recombination:
    # when the sweep starts, a lineage can be a recombination parent or a
    # common ancestor whose material at the selected site comes from one
    # child. The selected allele must still reach every sample below it.
    # After the sweep the allele drifts neutrally, so its mean count stays
    # N_SAMPLES * c.
    out = run([N_SAMPLES, reps, SITES, "-t", 0.0001, "-r", 200, "-ws", 0.05,
               "-a", 500, "-N", N_SWEEP, "-x", x, "-c", 0.5, "-d", 13, 14])
    counts = np.array([derived_count_at(r, x) or 0 for r in replicates(out)])
    z = (counts.mean() - N_SAMPLES * 0.5) / (counts.std() / np.sqrt(len(counts)))
    report(abs(z) < Z_MAX, f"partial sweep in the past with recombination: "
                           f"derived count mean {counts.mean():.2f} (expected "
                           f"{N_SAMPLES * 0.5:.1f}, z={z:+.1f})")

    # The same with the full ARG (-F), where each lineage has its own node.
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "f.trees")
        run([N_SAMPLES, 500, SITES, "-t", 0.0001, "-r", 200, "-ws", 0.05,
             "-a", 500, "-N", N_SWEEP, "-x", x, "-c", 0.5, "-d", 13, 14,
             "-ts", prefix, "-F"])
        counts = []
        for k in range(1, 501):
            ts = tskit.load(os.path.join(tmp, f"f_rep{k}.trees"))
            counts.append(next((int((v.genotypes > 0).sum())
                                for v in ts.variants()
                                if abs(v.site.position - x * SITES) < 1e-6), 0))
    counts = np.array(counts)
    z = (counts.mean() - N_SAMPLES * 0.5) / (counts.std() / np.sqrt(len(counts)))
    report(abs(z) < Z_MAX, f"the same with the full ARG (-F): derived count "
                           f"mean {counts.mean():.2f} (expected "
                           f"{N_SAMPLES * 0.5:.1f}, z={z:+.1f})")

    out = run([N_SAMPLES, 500, SITES, "-t", 0.0001, "-r", 20, "-ws", 0,
               "-a", 500, "-N", N_SWEEP, "-x", x, "-c", 0.5, "-h", "-d", 11, 12])
    hidden = sum(derived_count_at(r, x) is not None for r in replicates(out))
    report(hidden == 0, f"-h hides the selected site ({hidden} replicates "
                        "show it)")

    out = run([N_SAMPLES, 2000, SITES, "-t", 0.0001, "-r", 20, "-R", 4,
               "-a", 500, "-N", N_SWEEP, "-c", 0.5, "-d", 11, 12])
    nsites = np.array([len(r[0]) for r in replicates(out)])
    report(nsites.mean() > 2.0 and nsites.max() > 3,
           f"recurrent partial sweeps: {nsites.mean():.2f} selected sites "
           f"per replicate on average, maximum {nsites.max()}")

    # Ancient samples come first in their population's rows. Taken 0.4
    # generations before the end of a partial sweep to c = 0.5, each
    # carries the allele with probability about 0.5.
    n_anc, n_all = 10, N_SAMPLES + 10
    out = run([n_all, reps, SITES, "-t", 0.0001, "-r", 20, "-ws", 0,
               "-a", 500, "-N", N_SWEEP, "-x", x, "-c", 0.5,
               "-A", n_anc, 0, 0.00001, "-d", 11, 12])
    anc = []
    for positions, rows in replicates(out, n_all):
        k = next((i for i, p in enumerate(positions) if abs(p - x) < 1e-6),
                 None)
        anc.append(0 if k is None else sum(r[k] == "1" for r in rows[:n_anc]))
    anc = np.array(anc)
    z = (anc.mean() - n_anc * 0.5) / np.sqrt(n_anc * 0.25 / len(anc))
    report(abs(z) < Z_MAX, f"ancient samples during the sweep: derived count "
                           f"mean {anc.mean():.2f} of {n_anc} (about 5, z={z:+.1f})")

    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "p.trees")
        run([N_SAMPLES, 1, SITES, "-t", 0.0001, "-r", 20, "-ws", 0,
             "-a", 500, "-N", N_SWEEP, "-x", x, "-c", 0.5, "-d", 5, 6, "-ts", prefix])
        ts = tskit.load(prefix)
        positions = [s.position for s in ts.sites()]
        report(x * SITES in positions,
               f"-ts: selected site at {x * SITES:g} (sites at "
               f"{positions[:5]})")

    print("\nAll checks pass." if ok else "\nOne or more checks FAILED.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
