#!/usr/bin/env python3
"""Check diversity at a marker after a hard sweep against forward simulations
in SLiM, with crossovers and with gene conversion.

For a lineage whose material is only on one side of the selected site, the
class of the lineage (carrying the selected allele or not) changes when a
crossover or a gene conversion tract falls between the site and the material,
even if it does not cut the material. These events are what this check looks
at: mean pair coalescence time at a marker, after the sweep, relative to
neutral. It is compared with SLiM, a different simulator, so that an error
that discoal shares with a model it was written from would show.

Both simulate a locus of 100000 sites with the selected site at site 5000 and
one marker at site 35000. Only the recombination between the site and the
marker, and the material around them, matter.

discoal: 2 samples (the sample size where the error is largest), hard sweep
ending now with alpha = 2Ns = 100 (-a 50), N = 10000, -i 160.

SLiM (5.x): haploid Wright-Fisher, 2000 chromosomes, selection coefficient
s = alpha / 2000 for the allele, a new mutation at the selected site,
conditioned on fixation (restarts on loss). 20 chromosomes are drawn at
fixation, and their ancestry before the sweep is completed backward in time
with msprime, at the same recombination rate. Recombination: uniform over the
locus, with a total map length of rho / alpha times s (here 1, so the
recombination parameter of the locus equals alpha).

Cases (both at the same scaled distance a = rho * dx / alpha):

  crossover:  -r 100 in discoal; the SLiM rate per site is s / 100000.
  gene conversion, tracts of mean lambda = 1000 sites, no crossover: tracts
  that cover the site and not the marker, or the marker and not the site,
  arrive at R = 2 g lambda (1 - pc^d) (4N-scaled, g = gamma / nSites,
  pc = 1 - 1/lambda, d = 30000, see gene_conversion_check.py), the same as a
  crossover at an effective scaled distance R / alpha = 0.1. gamma is
  computed from that in discoal, and in SLiM the rate of conversion tracts per
  site (initializeGeneConversion, all events conversion) is computed to
  give the same tract rate per generation: 0.1 s = 2 r lambda (1 - pc^d).

The statistic is branch diversity (branch mode, so with no mutation noise) in a
one-site window at the marker; discoal's and SLiM's means over replicates are
each divided by the neutral mean (discoal with no sweep, and msprime for the
SLiM side). They must agree within 3% plus the standard error of the
difference (about 1.2% with the defaults). 3% allows for SLiM being a finite
population of 2000 chromosomes with discrete generations (its sweep ratios
scatter by 1 to 2% from the exact diffusion value), and the standard error
for the noise of the two simulations. Without the draws of crossovers and
tracts between the site and a lineage's material, discoal's ratio is 6 to 12%
below SLiM's.

Needs build/discoal, slim (5.x) on the PATH, msprime, tskit and numpy. With the
defaults, about 50 minutes on 20 workers for the SLiM sweeps (done once with
--slim-cache) and about 3 minutes for discoal. Fixed seeds, so the result does not change from run to run.
Exit code 0 if every check passes.

Usage:
    python testing/sweep_marker_slim_check.py [--workers 10] [--sweeps 8000]
        [--reps 60000] [--discoal PATH] [--slim-cache DIR]

--slim-cache DIR keeps the SLiM replicates in DIR, so that several discoal
builds can be compared without running SLiM again; more replicates are added
to what is there if --sweeps is raised.
"""

import argparse
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings

import msprime
import numpy as np
import tskit

# A SLiM tick is a generation in this model, as msprime assumes.
warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
L = 100_000
SITE = 5_000
MARKER = 35_000
N = 10_000
ALPHA = 100.0
SIZE = 2_000              # SLiM chromosomes
S = ALPHA / (2 * SIZE)    # SLiM selection coefficient
SAMPLES = 20              # SLiM samples
TRACT = 1000              # gene conversion mean tract length
A_EFF = 0.1               # gene conversion: effective scaled distance
PC = 1 - 1 / TRACT
D = MARKER - SITE
SEPARATION = 1 - PC ** D  # share of tracts long enough to reach from one site to the other
GAMMA = round(A_EFF * ALPHA * L / (2 * TRACT * SEPARATION))
SLIM_GC_RATE = A_EFF * S / (2 * TRACT * SEPARATION)
TOLERANCE = 0.03
CHUNK = 100

# name: (discoal arguments, SLiM recombination rate per site, SLiM gene conversion mean tract (0 for none))
CASES = {
    "crossover (rho = alpha)": (["-r", "100"], S / L, 0),
    f"gene conversion (mean tract {TRACT})": (["-r", "0", "-g", str(GAMMA), str(TRACT)], SLIM_GC_RATE, TRACT),
}

SLIM_MODEL = """
initialize() {
    setSeed(SEED);
    initializeTreeSeq();
    initializeMutationType("m2", 1.0, "f", SELECTION);
    initializeGenomicElementType("g1", m2, 1.0);
    initializeChromosome(1, LENGTH, type="H");
    initializeGenomicElement(g1);
    initializeMutationRate(0.0);
    initializeRecombinationRate(RATE);
    if (MEANLEN > 0) initializeGeneConversion(1.0, MEANLEN, 1.0);
}

1 early() {
    sim.addSubpop("p1", POPULATION);
}

1 late() {
    target = sample(p1.haplosomes, 1);
    target.addNewDrawnMutation(m2, SITE);
    defineConstant("SAVED", SAVE_PATH);
    sim.outputFull(SAVED);
}

1:100000000 late() {
    living = sim.mutationsOfType(m2);
    if (size(living) == 0) {
        if (sum(sim.substitutions.mutationType == m2) == 1) {
            sim.treeSeqOutput(READ_PATH);
            sim.simulationFinished();
        } else {
            sim.readFromPopulationFile(SAVED);
        }
    }
}
"""


def marker_diversity(ts):
    return ts.diversity(mode="branch", windows=[0, MARKER, MARKER + 1, L],
                        span_normalise=True)[1]


def msprime_rates(rate, mean_tract):
    if mean_tract:
        return dict(gene_conversion_rate=rate, gene_conversion_tract_length=mean_tract)
    return dict(recombination_rate=rate)


def slim_neutral(job):
    rate, mean_tract, seed = job
    rng = np.random.default_rng(seed)
    ts = msprime.sim_ancestry(SAMPLES, ploidy=1, population_size=SIZE, sequence_length=L,
                              discrete_genome=True, random_seed=int(rng.integers(1, 2**31)),
                              **msprime_rates(rate, mean_tract))
    return marker_diversity(ts)


def slim_sweep(job):
    rate, mean_tract, seed = job
    rng = np.random.default_rng(seed)
    with tempfile.TemporaryDirectory() as tmp:
        saved, read = os.path.join(tmp, "saved.trees"), os.path.join(tmp, "read.trees")
        constants = dict(SEED=seed, SAVE_PATH=f"'{saved}'", READ_PATH=f"'{read}'", SELECTION=S,
                         POPULATION=SIZE, LENGTH=L, RATE=rate, MEANLEN=mean_tract, SITE=SITE)
        args = ["slim"]
        for name, value in constants.items():
            args += ["-d", f"{name}={value}"]
        done = subprocess.run(args, input=SLIM_MODEL.encode("ascii"), capture_output=True)
        if done.returncode or not os.path.exists(read):
            raise RuntimeError(done.stderr.decode()[-2000:] + done.stdout.decode()[-2000:])
        traced = tskit.load(read)
    live = [u for u in traced.samples() if traced.node(u).metadata["is_vacant"][0] == 0]
    drawn = sorted(rng.choice(live, size=SAMPLES, replace=False).tolist())
    reduced = traced.simplify(drawn, keep_input_roots=True)
    demography = msprime.Demography.from_tree_sequence(reduced)
    for population in demography.populations:
        population.initial_size = SIZE
    full = msprime.sim_ancestry(initial_state=reduced, demography=demography, ploidy=1,
                                discrete_genome=True, random_seed=int(rng.integers(1, 2**31)),
                                **msprime_rates(rate, mean_tract))
    return marker_diversity(full)


def discoal_chunk(job):
    binary, args, sweep, reps, seed = job
    flags = ["-ws", "0", "-a", str(ALPHA / 2), "-x", str(SITE / L), "-i", "160"] if sweep else []
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([binary, "2", str(reps), str(L), "-t", "0", "-N", str(N)] + flags + args
                       + ["-d", str(seed), str(seed + 1), "-ts", prefix],
                       check=True, capture_output=True)
        return [marker_diversity(tskit.load(os.path.join(tmp, f"d_rep{k}.trees")))
                for k in range(1, reps + 1)]


def ratio(sweep, neutral):
    """Mean diversity after the sweep over neutral, and its standard error."""
    r = sweep.mean() / neutral.mean()
    se = r * np.sqrt(sweep.var() / len(sweep) / sweep.mean() ** 2
                     + neutral.var() / len(neutral) / neutral.mean() ** 2)
    return r, se


def slim_values(pool, name, rate, mean_tract, sweeps, cache):
    """SLiM sweeps and neutral replicates for a case. With a cache directory, what
    is there is reused and only the missing replicates are run."""
    path = os.path.join(cache, f"{name.split()[0]}.npz") if cache else None
    sweep, neutral = np.array([]), np.array([])
    if path and os.path.exists(path):
        z = np.load(path)
        sweep, neutral = z["sweep"], z["neutral"]
    seed = 100_000 * (1 + list(CASES).index(name))
    n_neutral = max(sweeps, 2000)
    if len(sweep) < sweeps:
        more = list(pool.map(slim_sweep, [(rate, mean_tract, seed + i) for i in range(len(sweep), sweeps)],
                             chunksize=1))
        sweep = np.concatenate([sweep, np.array(more)])
    if len(neutral) < n_neutral:
        more = list(pool.map(slim_neutral, [(rate, mean_tract, seed + 50_000 + i)
                                            for i in range(len(neutral), n_neutral)], chunksize=20))
        neutral = np.concatenate([neutral, np.array(more)])
    if path:
        os.makedirs(cache, exist_ok=True)
        np.savez(path, sweep=sweep, neutral=neutral)
    return sweep[:sweeps], neutral[:n_neutral]


def discoal_values(pool, binary, name, args, reps):
    seed = 7_000_000 + 100_000 * list(CASES).index(name)
    out = []
    for sweep in (True, False):
        jobs = [(binary, args, sweep, CHUNK, seed + (0 if sweep else 50_000) + 13 * k)
                for k in range(max(1, reps // CHUNK))]
        out.append(np.array([v for chunk in pool.map(discoal_chunk, jobs) for v in chunk]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--discoal", default=DISCOAL, help="binary to test")
    parser.add_argument("--sweeps", type=int, default=8000, help="SLiM sweeps per case")
    parser.add_argument("--reps", type=int, default=60000, help="discoal replicates per case")
    parser.add_argument("--slim-cache", default=None, help="directory for the SLiM replicates")
    opts = parser.parse_args()

    if shutil.which("slim") is None:
        print("slim was not found on the PATH; this check needs SLiM 5")
        return 1
    version = subprocess.run(["slim", "-v"], capture_output=True, text=True).stdout
    m = re.search(r"SLiM version (\d+)", version)
    if m is None or int(m.group(1)) < 5:
        print(f"this check needs SLiM 5 (haplosomes); found: {version.strip()}")
        return 1

    ok = True
    with concurrent.futures.ProcessPoolExecutor(opts.workers) as pool:
        for name, (args, rate, mean_tract) in CASES.items():
            ss, sn = slim_values(pool, name, rate, mean_tract, opts.sweeps, opts.slim_cache)
            ds, dn = discoal_values(pool, opts.discoal, name, args, opts.reps)
            (rs, ses), (rd, sed) = ratio(ss, sn), ratio(ds, dn)
            diff = rd / rs - 1
            se = (rd / rs) * np.sqrt((ses / rs) ** 2 + (sed / rd) ** 2)
            limit = TOLERANCE + se
            good = abs(diff) < limit
            ok &= good
            print(f"{'ok  ' if good else 'FAIL'} {name:34s} discoal {rd:.4f} (se {sed:.4f})  "
                  f"SLiM {rs:.4f} (se {ses:.4f}, {len(ss)} sweeps)  "
                  f"difference {100 * diff:+5.1f}% (limit +-{100 * limit:.1f}%)", flush=True)
    print("all checks passed" if ok else "some checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
