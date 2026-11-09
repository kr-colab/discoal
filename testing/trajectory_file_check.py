#!/usr/bin/env python3
"""Check sweeps with a trajectory read from a file (-S) in discoal.

The checks:

  1. a file holding discoal's deterministic trajectory gives the same
     genealogies as -wd: the mean total branch length and diversity
     (branch mode, from -ts output) agree within Z_MAX standard errors;
  2. a file whose first frequency is 0.5 gives a partial sweep: the number
     of samples that carry the selected SNP at -x is Binomial(n, 0.5);
  3. frequencies of 1 are used as 1 - 1/(2N), where a simulated sweep
     starts: a file that holds the allele at 1 before the sweep (going
     back in time) gives the same genealogies as the file without that
     plateau (at exactly 1 the wild-type coalescence rate is 0/0);
  4. -S runs with ancient samples taken during the sweep.

The comparison with -wd cannot be byte for byte: a simulated trajectory
uses a random number for its acceptance step, and a trajectory from a
file does not.

Needs build/discoal, numpy and tskit. Exit code 0 if every check passes.
"""

import glob
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import tskit

from sweep_snp_check import derived_count_at, replicates

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
Z_MAX = 4.0
N, ALPHA, TAU = 10000, 1000, 0.01   # -N, -a, -wd (4N0 generations)
N_SAMPLES, SITES, RHO = 10, 10000, 10
REPS = 1000


def write_deterministic(path):
    """discoal's deterministic trajectory (detSweepFreq), at its own time
    step 1/(40 N) in 2N0 generations, back from fixation at TAU until the
    frequency is 1/(2N), then 0. File times are in 4N0 generations."""
    eps = 0.05 / ALPHA
    ts = -2 * math.log(eps) / ALPHA
    dt = 1 / (40 * N)
    with open(path, "w") as f:
        f.write("# time (4N0 generations)  frequency\n")
        k = 0
        while True:
            t = k * dt
            x = 1.0 if k == 0 else eps / (eps + (1 - eps) * math.exp(ALPHA * (t - ts)))
            if x <= 1 / (2 * N):
                break
            f.write(f"{TAU + t / 2:.10g} {x:.10g}\n")
            k += 1
        f.write(f"{TAU + t / 2:.10g} 0\n")


def write_partial(path):
    """A logistic trajectory that is at 0.5 at the present and falls to 0
    going back in time."""
    with open(path, "w") as f:
        for i in range(400):
            t = i * 2.5e-5
            x = 1 / (1 + math.exp(ALPHA * 2 * t))
            f.write(f"{t:.10g} {x:.10g}\n")
        f.write(f"{400 * 2.5e-5:.10g} 0\n")


def write_rows(path, rows):
    with open(path, "w") as f:
        f.writelines(f"{t:g} {x:g}\n" for t, x in rows)


def compare(report, what, a, b):
    """Report whether two branch_stats results agree."""
    for name, (ma, sa, na), (mb, sb, nb) in zip(
            ("total branch length", "diversity"), a, b):
        z = (ma - mb) / math.sqrt(sa ** 2 + sb ** 2)
        report(abs(z) < Z_MAX and na == nb == REPS,
               f"{what}: {name} {ma:.4f} vs {mb:.4f} (z={z:+.1f}, "
               f"{na} replicates)")


def branch_stats(tmp, tag, sweep_args):
    """Mean total branch length and diversity per replicate, with their
    standard errors."""
    prefix = os.path.join(tmp, tag + ".trees")
    subprocess.run([DISCOAL, str(N_SAMPLES), str(REPS), str(SITES),
                    "-t", "1", "-r", str(RHO), "-N", str(N), "-ts", prefix,
                    "-d", "3", "4"] + sweep_args,
                   capture_output=True, text=True, check=True)
    length, pi = [], []
    for f in glob.glob(os.path.join(tmp, tag + "_rep*.trees")):
        ts = tskit.load(f)
        length.append(ts.segregating_sites(mode="branch"))
        pi.append(ts.diversity(mode="branch"))
    return [(np.mean(v), np.std(v, ddof=1) / math.sqrt(len(v)), len(v))
            for v in (length, pi)]


def main():
    ok = True

    def report(passed, text):
        nonlocal ok
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'}  {text}", flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        det = os.path.join(tmp, "det.traj")
        write_deterministic(det)
        compare(report, "-S with -wd's trajectory vs -wd",
                branch_stats(tmp, "file", ["-S", det]),
                branch_stats(tmp, "wd", ["-wd", str(TAU), "-a", str(ALPHA)]))

        plateau = os.path.join(tmp, "plateau.traj")
        no_plateau = os.path.join(tmp, "no_plateau.traj")
        sweep = [(0.05, 1), (0.051, 0.5), (0.052, 0)]
        write_rows(plateau, [(0, 1)] + sweep)
        write_rows(no_plateau, sweep)
        compare(report, "-S with a plateau at 1 vs without",
                branch_stats(tmp, "plateau", ["-S", plateau]),
                branch_stats(tmp, "no_plateau", ["-S", no_plateau]))

        part = os.path.join(tmp, "partial.traj")
        write_partial(part)
        x = 0.3
        out = subprocess.run(
            [DISCOAL, str(N_SAMPLES), str(REPS), str(SITES), "-t", "0.0001",
             "-r", str(RHO), "-N", str(N), "-S", part, "-x", str(x),
             "-d", "5", "6"],
            capture_output=True, text=True, check=True).stdout
        counts = np.array([derived_count_at(rep, x) or 0
                           for rep in replicates(out, N_SAMPLES)])
        z = (counts.mean() - N_SAMPLES * 0.5) / math.sqrt(
            N_SAMPLES * 0.25 / len(counts))
        report(abs(z) < Z_MAX,
               f"-S starting at 0.5: selected SNP carried by "
               f"{counts.mean():.2f} of {N_SAMPLES} samples on average "
               f"(about 5, z={z:+.1f})")

        for extra, what in ((["-A", "5", "0", str(TAU + 0.002)],
                             "ancient samples during the sweep"),):
            r = subprocess.run(
                [DISCOAL, str(N_SAMPLES), "20", str(SITES), "-t", "5",
                 "-r", str(RHO), "-N", str(N), "-S", det, "-d", "7", "8"] + extra,
                capture_output=True, text=True)
            report(r.returncode == 0 and r.stdout.count("//") == 20,
                   f"-S with {what}: exit {r.returncode}, "
                   f"{r.stdout.count('//')} replicates")

    print("\nAll checks pass." if ok else "\nOne or more checks FAILED.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
