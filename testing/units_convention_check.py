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
  - -r rho: 4 N_ref r L over all L sites, r per generation whatever the
    population size; compared on the number of distinct trees in tree
    sequence output (N_ref is N0 on the command line, the present size of
    the first deme for a demes file)

Usage: python testing/units_convention_check.py [--reps N]
Needs build/discoal, msprime, demes and tskit.
"""

import argparse
import glob
import math
import os
import re
import subprocess
import sys
import tempfile
import warnings

import demes
import msprime
import numpy as np
import tskit

# The same-time pulse case is deliberate; demes warns about it.
warnings.filterwarnings("ignore", message="Multiple pulses are defined")

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


def run_discoal(argv, reps, timeout=300):
    """Run discoal and return (mean segsites, standard error)."""
    try:
        out = subprocess.run([DISCOAL] + argv, capture_output=True,
                             text=True, check=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        sys.exit(f"discoal did not finish in {timeout} s: {argv}")
    s = [int(l.split()[1]) for l in out.splitlines()
         if l.startswith("segsites:")]
    if len(s) != reps:
        sys.exit(f"discoal returned {len(s)} of {reps} replicates: {argv}")
    return mean_se(s)


def discoal_argv(case, reps, seed, tmp):
    """Command-line arguments for a case. A demes case runs through a -Y
    config that loads the demes file; deme_sample_size follows the order of
    demes in the file, extinct demes included. Ancient samples of a demes
    case go in the config's ancient_samples, with times converted from
    generations to units of 4 N_ref generations (N_ref: the present size of
    the first deme that exists at time 0)."""
    n = sum(case["samples"].values())
    if "demes" not in case:
        return ([str(n), str(reps), "1000", "-t", str(THETA)]
                + case["flags"] + ["-d", str(seed), str(seed + 1)])
    graph = case["graph"]
    names = [d.name for d in graph.demes]
    sizes = [case["samples"].get(name, 0) for name in names]
    n_ref = next(d.epochs[-1].end_size for d in graph.demes
                 if d.end_time == 0)
    events = ""
    if case.get("ancient"):
        entries = ", ".join(
            f"{{time: {t / (4 * n_ref)!r}, population: {names.index(p)}, "
            f"sample_size: {k}}}" for p, k, t in case["ancient"])
        events = f", demographic_events: {{ancient_samples: [{entries}]}}"
    cfg = os.path.join(tmp, "config.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {n}, "
                f"num_replicates: {reps}, num_sites: 1000, "
                f"seed: [{seed}, {seed + 1}]}}\n"
                f"genetics: {{mutation_rate: {THETA}, recombination_rate: 0}}\n"
                f"demography: {{deme_sample_size: {sizes}, "
                f"demes_filename: \"{case['demes']}\"{events}}}\n")
    return ["-Y", cfg]


def msprime_mean(demography, samples, reps, seed, ancient=()):
    """Expected segsites under msprime: mean branch length times MU.
    `ancient` lists (population, count, time in generations); as with
    discoal's -A, these samples are part of the counts in `samples`.
    `demography` may also be a function of a numpy Generator that returns a
    Demography; it is then called once per replicate (a model with priors)."""
    present = dict(samples)
    for p, k, _ in ancient:
        present[p] -= k
    sets = [msprime.SampleSet(k, population=p, ploidy=1)
            for p, k in present.items() if k > 0]
    sets += [msprime.SampleSet(k, population=p, time=t, ploidy=1)
             for p, k, t in ancient]
    if callable(demography):
        rng = np.random.default_rng(seed)
        return mean_se([msprime.sim_ancestry(
            samples=sets, demography=demography(rng), ploidy=2,
            random_seed=int(rng.integers(1, 2**31))).segregating_sites(
                mode="branch") * MU for _ in range(reps)])
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


def two_pop(rate, merge_gens, source, dest, growth_p1=0.0):
    """Pops P0 and P1; lineages move from source to dest (backward in time)
    only; both merge into ANC at merge_gens. All samples come from P0, so
    the direction of migration changes the model. growth_p1 is P1's growth
    rate per generation."""
    d = msprime.Demography()
    for name in ("P0", "P1", "ANC"):
        d.add_population(name=name, initial_size=N0,
                         growth_rate=growth_p1 if name == "P1" else 0.0)
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


# Two present-day demes that merge into an extinct ancestor of a different
# size. The ancestor's size must be set where its last epoch ends.
DEMES_EXTINCT_ANCESTOR = """time_units: generations
demes:
  - name: ANC
    epochs:
      - end_time: 20000
        start_size: 2000
  - name: A
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 10000
  - name: B
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 10000
migrations:
  - demes: [A, B]
    rate: 1.0e-5
"""


# A deme formed from three ancestors of different sizes. The first
# present-day deme (A) sets discoal's N_ref, so it has size N0.
DEMES_THREE_WAY_ADMIXTURE = """time_units: generations
demes:
  - name: ANC
    epochs:
      - end_time: 20000
        start_size: 10000
  - name: A
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 10000
  - name: B
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 1000
  - name: C
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 50000
  - name: D
    ancestors: [A, B, C]
    proportions: [0.2, 0.3, 0.5]
    start_time: 100
    epochs:
      - end_time: 0
        start_size: 10000
"""


# A pulse from two sources into A. Lineages not moved must stay in A.
DEMES_TWO_SOURCE_PULSE = """time_units: generations
demes:
  - name: ANC
    epochs:
      - end_time: 20000
        start_size: 10000
  - name: A
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 10000
  - name: B
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 1000
  - name: C
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 50000
pulses:
  - sources: [B, C]
    dest: A
    time: 500
    proportions: [0.2, 0.3]
"""


# One-way migration from A into B, forward in time. Backward in time,
# lineages sampled in B move to A.
DEMES_ONE_WAY_MIGRATION = """time_units: generations
demes:
  - name: A
    epochs:
      - end_time: 0
        start_size: 10000
  - name: B
    ancestors: [A]
    start_time: 20000
    epochs:
      - end_time: 0
        start_size: 2000
migrations:
  - {source: A, dest: B, rate: 1.0e-4}
"""


# An extinct deme (MID) that is not population 0 and has no lineages until
# B merges into it at 10000 generations. Before then its size is unset (0).
DEMES_EMPTY_EXTINCT_DEME = """time_units: generations
demes:
  - name: X
    epochs:
      - end_time: 0
        start_size: 10000
  - name: MID
    ancestors: [X]
    start_time: 30000
    epochs:
      - end_time: 10000
        start_size: 3000
  - name: B
    ancestors: [MID]
    start_time: 10000
    epochs:
      - end_time: 0
        start_size: 5000
migrations:
  - {demes: [X, B], rate: 1.0e-5}
"""


# Two pulses at the same time. demes applies them in the order given,
# forward in time, so backward in time C->B runs before B->A.
DEMES_SAME_TIME_PULSES = """time_units: generations
demes:
  - name: ANC
    epochs: [{start_size: 10000, end_time: 20000}]
  - name: A
    ancestors: [ANC]
    epochs: [{start_size: 10000}]
  - name: B
    ancestors: [ANC]
    epochs: [{start_size: 10000}]
  - name: C
    ancestors: [ANC]
    epochs: [{start_size: 200000}]
pulses:
  - {sources: [B], dest: A, time: 1000, proportions: [0.8]}
  - {sources: [C], dest: B, time: 1000, proportions: [0.8]}
"""

# A pulse into A at the time D forms from A and B. Backward in time the pulse
# runs first, so lineages that just left D take no part in it.
DEMES_PULSE_AND_ADMIXTURE = """time_units: generations
demes:
  - name: ANC
    epochs: [{start_size: 10000, end_time: 20000}]
  - name: A
    ancestors: [ANC]
    epochs: [{start_size: 10000}]
  - name: B
    ancestors: [ANC]
    epochs: [{start_size: 10000}]
  - name: X
    ancestors: [ANC]
    epochs: [{start_size: 100000}]
  - name: D
    ancestors: [A, B]
    proportions: [0.5, 0.5]
    start_time: 1000
    epochs: [{start_size: 10000}]
pulses:
  - {sources: [X], dest: A, time: 1000, proportions: [0.9]}
"""


# Two root demes joined only by a pulse that moves all lineages of B into A
# (backward in time): no split and no migration connects them.
DEMES_PULSE_ONLY = """time_units: generations
demes:
  - name: A
    epochs: [{start_size: 10000}]
  - name: B
    epochs: [{start_size: 2000}]
pulses:
  - {sources: [A], dest: B, time: 1000, proportions: [1.0]}
"""


# A present-day deme A and a ghost deme X that went extinct 2000
# generations ago, both from ANC; they exchange migrants while both exist.
# All of X's samples are ancient. N_ref is A's size.
DEMES_EXTINCT_GHOST = """time_units: generations
demes:
  - name: ANC
    epochs:
      - end_time: 8000
        start_size: 10000
  - name: A
    ancestors: [ANC]
    epochs:
      - end_time: 0
        start_size: 10000
  - name: X
    ancestors: [ANC]
    epochs:
      - end_time: 2000
        start_size: 5000
migrations:
  - demes: [A, X]
    rate: 1e-4
"""


# An ancient sample of C taken at the time of a pulse that moves all of
# C's lineages into A (backward in time). msprime takes the sample after
# the pulse, so the sample stays in C, and so must discoal.
DEMES_SAMPLE_AT_PULSE = """time_units: generations
demes:
  - name: ANC
    epochs: [{start_size: 10000, end_time: 20000}]
  - name: C
    ancestors: [ANC]
    epochs: [{start_size: 10000}]
  - name: A
    ancestors: [ANC]
    epochs: [{start_size: 1000}]
pulses:
  - {sources: [A], dest: C, time: 500, proportions: [1.0]}
"""


def demes_case(tmp, name, samples, text=None, path=None, ancient=()):
    """A case that loads a demes model, given as YAML text or a file path,
    and compares it with msprime's Demography.from_demes. `ancient` lists
    (deme, count, time in generations), as in msprime_mean."""
    if path is None:
        path = os.path.join(tmp, re.sub(r"\W+", "_", name) + ".demes.yaml")
        with open(path, "w") as f:
            f.write(text)
    graph = demes.load(path)
    return dict(name=name, samples=samples, demes=path, graph=graph, alts={},
                ancient=list(ancient),
                model=msprime.Demography.from_demes(graph))


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

    def prior_en(k):
        """-Pe1 0.2 0.6 0.1 1.0: time uniform in [0.2, 0.6] (k N0
        generations per unit) and size uniform in [0.1, 1.0] N0, drawn
        for each replicate."""
        return lambda rng: one_pop([(rng.uniform(0.2, 0.6) * k * N0,
                                     {"initial_size":
                                      rng.uniform(0.1, 1.0) * N0})])

    def admixed(k, p):
        """P0 is admixed at 0.1 (k N0 generations per unit): its lineages
        go to P1 (size 0.1 N0) with probability p, else to P2; P2 merges
        into P1 at 0.5."""
        d = msprime.Demography()
        d.add_population(name="P0", initial_size=N0)
        d.add_population(name="P1", initial_size=0.1 * N0,
                         initially_active=True)
        d.add_population(name="P2", initial_size=N0, initially_active=True)
        d.add_admixture(time=0.1 * k * N0, derived="P0",
                        ancestral=["P1", "P2"], proportions=[p, 1 - p])
        d.add_population_split(time=0.5 * k * N0, derived=["P2"],
                               ancestral="P1")
        return d

    def ancient(k, sets):
        """-A sample sets (population, count, time) with times in k N0."""
        return [(pop, n, t * k * N0) for pop, n, t in sets]

    # Lineages in P1 move to P0 backward in time; P1 merges into P0.
    one_way = msprime.Demography()
    # msprime keeps the ancestral population of a split inactive until the
    # split unless told otherwise; P0 has samples and receives migrants.
    one_way.add_population(name="P0", initial_size=N0, initially_active=True)
    one_way.add_population(name="P1", initial_size=N0)
    one_way.set_migration_rate(source="P1", dest="P0", rate=m / (4 * N0))
    one_way.add_population_split(time=4 * N0, derived=["P1"], ancestral="P0")
    anc1 = [("A", 3, 0.25), ("A", 2, 1.0)]
    anc2 = [("P1", 4, 0.3), ("P0", 3, 0.6), ("P1", 2, 0.6)]

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
        # Growth makes the sampler use the time-varying (NHPP) path, which
        # draws migration times separately from the constant-rate path. With
        # alpha < 0 the past was larger, so migration out of P1 matters.
        dict(name="-m with -eg: migration rate while sizes change",
             samples={"P1": 20},
             flags=["-p", "2", "0", "20", "-m", "1", "0", str(m),
                    "-eg", "0", "1", "-5", "-ed", "1.0", "1", "0"],
             model=two_pop(m / (4 * N0), 4 * N0, "P1", "P0",
                           growth_p1=-5 / (4 * N0)),
             alts={"rate 2 x 4N0 m": two_pop(2 * m / (4 * N0), 4 * N0,
                                            "P1", "P0",
                                            growth_p1=-5 / (4 * N0))}),
        # Recurrent-sweep mode at a rate so low that no sweep happens must
        # give the neutral coalescent: the recurrent sampler used to leave
        # the edges of the last coalescence unwritten, which removed the
        # root branch of every genealogy.
        dict(name="-R at a near-zero rate gives the neutral coalescent",
             samples=s20, flags=["-R", "0.000001", "-a", "500"],
             model=one_pop([]), alts={}),
        dict(name="-ea proportion is the share that goes to the first source",
             samples={"P0": 20},
             flags=["-p", "3", "20", "0", "0", "-en", "0", "1", "0.1",
                    "-ea", "0.1", "0", "1", "2", "0.3", "-ed", "0.5", "2", "1"],
             model=admixed(4, 0.3),
             alts={"share to the second source": admixed(4, 0.7),
                   "2N0 generations": admixed(2, 0.3)}),
        dict(name="-Pe1: size change with priors on time and size",
             samples=s20, flags=["-Pe1", "0.2", "0.6", "0.1", "1.0"],
             model=prior_en(4), alts={"2N0 generations": prior_en(2)}),
        # The size change with a prior must not replace another event (it
        # used to write to events[1], which here is the -en event).
        dict(name="-Pe1 with a fixed -en at an earlier time",
             samples=s20,
             flags=["-en", "0.1", "0", "3.0", "-Pe1", "0.5", "0.5", "0.2", "0.2"],
             model=one_pop([(0.1 * 4 * N0, {"initial_size": 3 * N0}),
                            (0.5 * 4 * N0, {"initial_size": 0.2 * N0})]),
             alts={"-en replaced":
                   one_pop([(0.5 * 4 * N0, {"initial_size": 0.2 * N0})])}),
        dict(name="-A time in 4N0 generations, two ancient sets",
             samples=s20,
             flags=["-A", "3", "0", "0.25", "-A", "2", "0", "1.0"],
             model=one_pop([]), ancient=ancient(4, anc1),
             alts={"2N0 generations": (one_pop([]), ancient(2, anc1))}),
        dict(name="-A in two populations with migration and a merge",
             samples={"P0": 10, "P1": 10},
             flags=["-p", "2", "10", "10", "-m", "1", "0", str(m),
                    "-ed", "1.0", "1", "0",
                    "-A", "4", "1", "0.3", "-A", "3", "0", "0.6",
                    "-A", "2", "1", "0.6"],
             model=one_way, ancient=ancient(4, anc2),
             alts={"2N0 generations": (one_way, ancient(2, anc2))}),
        demes_case(tmp, "demes exponential epoch", s20, DEMES_GROWTH),
        demes_case(tmp, "demes exponential epoch between constant epochs",
                   s20, DEMES_GROWTH_MIDDLE),
        demes_case(tmp, "demes extinct ancestor with its own size",
                   {"A": 10, "B": 10}, DEMES_EXTINCT_ANCESTOR),
        demes_case(tmp, "demes deme with three ancestors", {"D": 20},
                   DEMES_THREE_WAY_ADMIXTURE),
        demes_case(tmp, "demes pulse with two sources", {"A": 20},
                   DEMES_TWO_SOURCE_PULSE),
        demes_case(tmp, "demes one-way migration direction", {"B": 20},
                   DEMES_ONE_WAY_MIGRATION),
        demes_case(tmp, "demes extinct deme that is empty at first",
                   {"X": 10, "B": 10}, DEMES_EMPTY_EXTINCT_DEME),
        demes_case(tmp, "demes pulses at the same time", {"A": 20},
                   DEMES_SAME_TIME_PULSES),
        demes_case(tmp, "demes pulse at the time of an admixture", {"D": 20},
                   DEMES_PULSE_AND_ADMIXTURE),
        demes_case(tmp, "demes joined only by a pulse", {"A": 10, "B": 10},
                   DEMES_PULSE_ONLY),
        demes_case(tmp, "demes ancient samples, one from an extinct deme",
                   {"A": 12, "X": 8}, DEMES_EXTINCT_GHOST,
                   ancient=[("X", 8, 3000), ("A", 4, 1000)]),
        demes_case(tmp, "demes ancient sample at the time of a pulse",
                   {"C": 8, "A": 4}, DEMES_SAMPLE_AT_PULSE,
                   ancient=[("C", 4, 500)]),
        demes_case(tmp, "demes splits, sizes, migration windows",
                   {"A": 10, "B": 5, "C": 5},
                   path=os.path.join(REPO, "config_examples",
                                     "demes_example.demes.yaml")),
    ]


def discoal_tree_mean(argv, tmp):
    """Run discoal with tree sequence output; return (mean, standard error)
    of the number of distinct trees along the sequence."""
    for f in glob.glob(os.path.join(tmp, "t_rep*.trees")):
        os.remove(f)
    subprocess.run([DISCOAL] + argv + ["-ts", os.path.join(tmp, "t.trees")],
                   capture_output=True, check=True, timeout=600)
    return mean_se([tskit.load(f).num_trees
                    for f in glob.glob(os.path.join(tmp, "t_rep*.trees"))])


def msprime_tree_mean(demography, samples, sites, r, reps, seed):
    """Mean number of distinct trees under msprime on a discrete genome;
    samples is {population: haploid sample count}."""
    return mean_se([ts.num_trees for ts in msprime.sim_ancestry(
        samples=[msprime.SampleSet(k, population=p, ploidy=1)
                 for p, k in samples.items()],
        demography=demography, ploidy=2, sequence_length=sites,
        recombination_rate=r, num_replicates=reps, random_seed=seed)])


def tree_cases(tmp, reps):
    """Recombination cases, compared on the number of distinct trees. Each
    case has discoal arguments, the msprime model and population, the
    number of sites L, the per-generation rate r for the documented
    convention (rho = 4 N_ref r L), rates that must be rejected, and its
    number of replicates. Only the L versus L - 1 case needs `reps`; the
    others reject their alternatives by a wide margin with fewer, and their
    many trees are slow to write."""
    rho, n = 10.0, 10
    few = max(2, reps // 10)

    def r(N, L):
        return rho / (4 * N * L)

    grown = msprime.Demography.isolated_model([N0])
    grown.add_population_parameters_change(time=0.05 * 4 * N0,
                                           initial_size=100 * N0,
                                           population="pop_0")
    # Two present-day demes of different sizes, and a much larger past size:
    # rho must be relative to the present size of the first deme (A, 300),
    # not to B (3000) or to A's past size (90000). B merges into A at 1000
    # generations, so lineages coalesce in A long before the large epoch.
    text = ("time_units: generations\ndemes:\n"
            "  - name: A\n    epochs:\n"
            "      - {end_time: 50000, start_size: 90000}\n"
            "      - {end_time: 0, start_size: 300}\n"
            "  - name: B\n    ancestors: [A]\n    start_time: 1000\n"
            "    epochs: [{start_size: 3000}]\n")
    path = os.path.join(tmp, "two_present_demes.demes.yaml")
    with open(path, "w") as f:
        f.write(text)
    half = max(2, reps // 2)
    cfg = os.path.join(tmp, "two_present_demes.yaml")
    with open(cfg, "w") as f:
        f.write(f"simulation: {{sample_size: {n}, num_replicates: {half}, "
                f"num_sites: 1000, seed: [7, 8]}}\n"
                f"genetics: {{mutation_rate: 0, recombination_rate: {rho}}}\n"
                f"demography: {{deme_sample_size: [{n // 2}, {n // 2}], "
                f"demes_filename: \"{path}\"}}\n")
    return [
        dict(name="-r rho is 4N0 r over all L sites (L = 10)", reps=reps,
             argv=[str(n), str(reps), "10", "-t", "0", "-r", str(rho),
                   "-d", "3", "4"],
             model=msprime.Demography.isolated_model([N0]),
             samples={"pop_0": n}, sites=10, r=r(N0, 10),
             alts={"over L - 1 sites": rho / (4 * N0 * 9),
                   "2 rho": 2 * r(N0, 10)}),
        dict(name="-r with a past size 100 x N0: r stays per generation",
             reps=few,
             argv=[str(n), str(few), "1000", "-t", "0", "-r", str(rho),
                   "-en", "0.05", "0", "100", "-d", "5", "6"],
             model=grown, samples={"pop_0": n}, sites=1000, r=r(N0, 1000),
             alts={"rho relative to the past size": r(100 * N0, 1000)}),
        dict(name="demes: rho is relative to the present size of the first deme",
             reps=half, argv=["-Y", cfg], model=msprime.Demography.from_demes(
                 demes.loads(text)), samples={"A": n // 2, "B": n // 2},
             sites=1000, r=r(300, 1000),
             alts={"rho relative to deme B (3000)": r(3000, 1000),
                   "rho relative to A's past size (90000)": r(90000, 1000)}),
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
            anc = case.get("ancient", ())
            ref = msprime_mean(case["model"], case["samples"], reps, seed,
                               anc)
            # An alternative is a model, or a (model, ancient samples) pair.
            alts = {}
            for k, mm in case["alts"].items():
                model, alt_anc = mm if isinstance(mm, tuple) else (mm, anc)
                alts[k] = msprime_mean(model, case["samples"], reps,
                                       seed + 1, alt_anc)
            all_ok &= report(case["name"], d, ref, alts)

        for i, case in enumerate(tree_cases(tmp, reps)):
            seed = 2010 + 10 * i
            d = discoal_tree_mean(case["argv"], tmp)
            args = (case["model"], case["samples"], case["sites"])
            ref = msprime_tree_mean(*args, case["r"], case["reps"], seed)
            alts = {k: msprime_tree_mean(*args, rk, case["reps"], seed + 1)
                    for k, rk in case["alts"].items()}
            all_ok &= report(case["name"] + " (number of trees)", d, ref,
                             alts)

    print("\nAll unit conventions hold." if all_ok
          else "\nOne or more unit conventions FAILED.")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
