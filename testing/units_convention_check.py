#!/usr/bin/env python3
"""Pin the units of discoal's user-facing time and rate parameters.

Each case runs discoal and msprime on the same model and compares the mean
number of segregating sites. discoal reports site counts; msprime reports
the branch-mode expectation (total branch length times the per-generation
mutation rate), which has no mutational noise.

For each case the documented convention must match (|z| < Z_MATCH), and
every listed alternative convention must be rejected (|z| > Z_REJECT).
The second check shows that the test can tell the conventions apart, so a
change of scale cannot pass unseen.

Conventions under test (N0 is the reference size, theta = 4 N0 mu):
  - times (-en, -ed, -em, -eM): units of 4N0 generations
  - -eg / -eG alpha: growth rate per 4N0 generations (as in ms)
  - -el / -eL gamma: change in relative size per 4N0 generations
  - -M / -m / -em rates: 4N0 m; -m i j moves lineages from i to j
    backward in time (msprime's convention)
  - demes import: demes generations and per-generation rates

Usage: python testing/units_convention_check.py [--reps N]
Needs build/discoal, msprime and demes.
"""

import argparse
import math
import os
import subprocess
import sys
import tempfile

import demes
import msprime
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
N0 = 10000
THETA = 20.0
MU = THETA / (4 * N0)  # per-locus, per-generation mutation rate
Z_MATCH = 4.0
Z_REJECT = 6.0


def mean_se(x):
    x = np.asarray(x, dtype=float)
    return x.mean(), x.std(ddof=1) / math.sqrt(len(x))


def z(a, b):
    """z-score of the difference of two (mean, standard error) pairs."""
    return (a[0] - b[0]) / math.hypot(a[1], b[1])


def run_discoal(argv, reps):
    """Run discoal and return (mean segsites, standard error)."""
    out = subprocess.run([DISCOAL] + argv, capture_output=True, text=True,
                         check=True).stdout
    s = [int(l.split()[1]) for l in out.splitlines()
         if l.startswith("segsites:")]
    if len(s) != reps:
        sys.exit(f"discoal returned {len(s)} of {reps} replicates: {argv}")
    return mean_se(s)


def discoal_argv(case, reps, seed, tmp):
    """Command-line arguments for a case. A demes case runs through a -Y
    config that loads the demes file."""
    sizes = list(case["samples"].values())
    if "demes" not in case:
        return ([str(sum(sizes)), str(reps), "1000", "-t", str(THETA)]
                + case["flags"] + ["-d", str(seed), str(seed + 1)])
    cfg = os.path.join(tmp, "config.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {sum(sizes)}, "
                f"num_replicates: {reps}, num_sites: 1000, "
                f"seed: [{seed}, {seed + 1}]}}\n"
                f"genetics: {{mutation_rate: {THETA}, recombination_rate: 0}}\n"
                f"demography: {{deme_sample_size: {sizes}, "
                f"demes_filename: \"{case['demes']}\"}}\n")
    return ["-Y", cfg]


def msprime_mean(demography, samples, reps, seed):
    """Expected segsites under msprime: mean branch length times MU."""
    sets = [msprime.SampleSet(k, population=p, ploidy=1)
            for p, k in samples.items()]
    return mean_se([ts.segregating_sites(mode="branch") * MU
                    for ts in msprime.sim_ancestry(
                        samples=sets, demography=demography, ploidy=2,
                        num_replicates=reps, random_seed=seed)])


def one_pop(changes):
    """Single population of size N0 with (time_gens, kwargs) changes."""
    d = msprime.Demography()
    d.add_population(name="A", initial_size=N0)
    for t, kw in changes:
        d.add_population_parameters_change(time=t, population="A", **kw)
    return d


def linear_pop(t1_gens, gamma_rel_per_gen, final_rel, steps=400):
    """Linear size change with forward-time slope gamma_rel_per_gen
    (relative size per generation) from the present back to t1_gens, then
    constant at final_rel * N0. Approximated by `steps` constant pieces."""
    changes = [(t1_gens * i / steps,
                {"initial_size": N0 * (1.0 - gamma_rel_per_gen
                                       * t1_gens * (i + 0.5) / steps)})
               for i in range(steps)]
    changes.append((t1_gens, {"initial_size": final_rel * N0}))
    return one_pop(changes)


def two_pop(rate, merge_gens, source, dest):
    """Pops P0 and P1; lineages move from source to dest (backward in time)
    only; both merge into ANC at merge_gens. All samples come from P0, so
    the direction of migration changes the model."""
    d = msprime.Demography()
    for name in ("P0", "P1", "ANC"):
        d.add_population(name=name, initial_size=N0)
    d.set_migration_rate(source=source, dest=dest, rate=rate)
    d.add_population_split(time=merge_gens, derived=["P0", "P1"],
                           ancestral="ANC")
    return d


DEMES_GROWTH = """time_units: generations
demes:
  - name: A
    epochs:
      - end_time: 20000
        start_size: 1000
      - end_time: 0
        start_size: 1000
        end_size: 10000
        size_function: exponential
"""


# A growth epoch between two constant epochs: at 1000 generations the importer
# emits a size change and a growth change at the same time.
DEMES_GROWTH_MIDDLE = """time_units: generations
demes:
  - name: A
    epochs:
      - end_time: 20000
        start_size: 1000
      - end_time: 1000
        start_size: 1000
        end_size: 10000
        size_function: exponential
      - end_time: 0
        start_size: 10000
"""


def cases(tmp):
    """Each case has: name, samples ({population: count}), either discoal
    `flags` or a `demes` file, the msprime `model` for the documented
    convention, and `alts`, msprime models for conventions that must be
    rejected. k is the number of N0 in one time unit (4 documented, 2 not)."""
    alpha, gamma, m = 10.0, 0.8, 2.0
    s20 = {"A": 20}

    def en(k):
        return one_pop([(0.5 * k * N0, {"initial_size": 0.1 * N0})])

    def growth(k):
        return one_pop([(0, {"growth_rate": alpha / (k * N0)}),
                        (0.5 * 4 * N0, {"growth_rate": 0})])

    def linear(k):
        return linear_pop(0.5 * 4 * N0, gamma / (k * N0), 1 - gamma * 0.5)

    growth_file = os.path.join(tmp, "growth.demes.yaml")
    with open(growth_file, "w") as f:
        f.write(DEMES_GROWTH)
    middle_file = os.path.join(tmp, "growth_middle.demes.yaml")
    with open(middle_file, "w") as f:
        f.write(DEMES_GROWTH_MIDDLE)
    example = os.path.join(REPO, "config_examples", "demes_example.demes.yaml")

    return [
        dict(name="-en time in 4N0 generations", samples=s20,
             flags=["-en", "0.5", "0", "0.1"],
             model=en(4), alts={"2N0 generations": en(2)}),
        dict(name="-eg alpha per 4N0 generations", samples=s20,
             flags=["-eg", "0", "0", str(alpha), "-eg", "0.5", "0", "0"],
             model=growth(4), alts={"per 2N0 generations": growth(2)}),
        dict(name="-eG alpha per 4N0 generations", samples=s20,
             flags=["-eG", "0", str(alpha), "-eG", "0.5", "0"],
             model=growth(4), alts={"per 2N0 generations": growth(2)}),
        dict(name="-el gamma per 4N0 generations", samples=s20,
             flags=["-el", "0", "0", str(gamma),
                    "-en", "0.5", "0", str(1 - gamma * 0.5)],
             model=linear(4), alts={"per 2N0 generations": linear(2)}),
        dict(name="-m i j rate is 4N0 m, lineages move i -> j backward",
             samples={"P0": 20},
             flags=["-p", "2", "20", "0", "-m", "0", "1", str(m),
                    "-ed", "1.0", "1", "0"],
             model=two_pop(m / (4 * N0), 4 * N0, "P0", "P1"),
             alts={"rate 2N0 m": two_pop(m / (2 * N0), 4 * N0, "P0", "P1"),
                   "lineages move j -> i":
                       two_pop(m / (4 * N0), 4 * N0, "P1", "P0")}),
        dict(name="demes exponential epoch", samples=s20, demes=growth_file,
             model=msprime.Demography.from_demes(demes.load(growth_file)),
             alts={}),
        dict(name="demes exponential epoch between constant epochs",
             samples=s20, demes=middle_file,
             model=msprime.Demography.from_demes(demes.load(middle_file)),
             alts={}),
        dict(name="demes splits, sizes, migration windows",
             samples={"A": 10, "B": 5, "C": 5}, demes=example,
             model=msprime.Demography.from_demes(demes.load(example)),
             alts={}),
    ]


def report(name, d, ref, alts):
    zr = z(d, ref)
    ok = abs(zr) < Z_MATCH
    line = f"{'PASS' if ok else 'FAIL'}  {name}: discoal {d[0]:.2f}, msprime {ref[0]:.2f} (z={zr:+.1f})"
    for alt_name, a in alts.items():
        za = z(d, a)
        rejected = abs(za) > Z_REJECT
        ok = ok and rejected
        line += f"\n      alternative '{alt_name}': msprime {a[0]:.2f} (z={za:+.1f}) {'rejected' if rejected else 'NOT REJECTED'}"
    print(line, flush=True)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=4000)
    reps = ap.parse_args().reps
    all_ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for i, case in enumerate(cases(tmp)):
            seed = 1010 + 10 * i
            d = run_discoal(discoal_argv(case, reps, seed, tmp), reps)
            ref = msprime_mean(case["model"], case["samples"], reps, seed)
            alts = {k: msprime_mean(mm, case["samples"], reps, seed + 1)
                    for k, mm in case["alts"].items()}
            all_ok &= report(case["name"], d, ref, alts)

    print("\nAll unit conventions hold." if all_ok
          else "\nOne or more unit conventions FAILED.")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
