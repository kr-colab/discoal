#!/usr/bin/env python3
"""Compare discoal's sweeps with forward simulations in SLiM.

The SLiM script is the one stdpopsim writes for its constant-size model
(PiecewiseConstantSize) on a generic contig with one drawn mutation. Only
the sweep logic is replaced, so that it follows discoal's -ws / -f / -c /
-wn model, which stdpopsim's fixed-time sweep events cannot express:

  - the drawn mutation is neutral until its frequency first reaches f0
    (for a hard sweep, f0 = 0: selected from its origin; for -wn it stays
    neutral);
  - it then has selection coefficient alpha/N with heterozygote fitness
    1 + alpha/(2N). stdpopsim and SLiM use msprime's convention (W_BB =
    1 + s), and discoal's s is half of that: alpha = 2Ns with per-copy
    advantage s;
  - if it is lost, the run restarts from the state saved just before the
    mutation was drawn (stdpopsim's save and restore);
  - the run ends when the allele fixes (or reaches c, for -c), as -ws 0.

The statistics are branch-mode total branch length and diversity in
generations, for 20 genomes over the whole contig (no mutation noise).
msprime's SweepGenicSelection gives a diffusion reference for the hard
sweep.

What to expect: discoal agrees with msprime's hard sweep. SLiM is a
Wright-Fisher population with discrete generations, and its values are a
few percent above both diffusion models (about 3% at N = 1000, 4 to 5% at
N = 5000); discoal's soft and partial sweeps should differ from SLiM by
about as much as its hard sweep does. Before the soft-sweep fix, discoal's
soft sweep with f0 = 0.05 was 23% above SLiM.

Needs build/discoal, stdpopsim, SLiM (slim on the PATH), pyslim, msprime,
tskit and numpy. stdpopsim 0.3 writes scripts for SLiM 4; for SLiM 5 the
one call that uses genomes is changed to haplosomes. A run of all cases
with 1000 replicates takes about 7 minutes on 120 cores at N = 1000.

Usage:
    python testing/slim_sweep_check.py [--reps 1000] [--workers 100]
        [--n 1000] [--cases "hard,soft f0=0.05"]
"""

import argparse
import concurrent.futures
import contextlib
import io
import os
import subprocess
import tempfile

import msprime
import numpy as np
import pyslim
import stdpopsim
import tskit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
L = 100_000
N_GENOMES = 20
SWEEP_SITE = L // 2

SWEEP_LOGIC = r"""
// discoal's sweep model in place of stdpopsim's fixed-time sweep events.
// Phase 0: the drawn mutation is neutral; phase 1: selected, once its
// frequency first reaches F0. The run ends when it reaches FINAL (as -ws 0).
function (void)sweep_step(void) {
    f = af(m1, p0);
    phase = sim.getValue("phase");
    if (f == 0.0) {
        sim.setValue("phase", 0);
        restore();
        return;
    }
    if (phase == 0 & f >= F0) {
        sim.mutationsOfType(m1).setSelectionCoeff(S_HOM);
        sim.setValue("phase", 1);
        phase = 1;
    }
    if (phase == 1 & f >= FINAL) {
        metadata.setValue("sweep_done", 1);
        end();
    }
}

1 early() {
    sim.setValue("phase", 0);
    g_mut = time_to_tick(drawn_mutations[0, 0]);
    // G0 is the last tick (the present) of this constant-size model.
    community.registerLateEvent(NULL, "{sweep_step();}", g_mut + 1, G0);
}
"""


def stdpopsim_script(n, rho, t_mut):
    """stdpopsim's SLiM script with one mutation drawn t_mut generations
    ago, and the state saved before the draw."""
    sp = stdpopsim.get_species("HomSap")
    model = stdpopsim.PiecewiseConstantSize(n)
    contig = sp.get_contig(length=L, mutation_rate=0,
                           recombination_rate=rho / (4 * n * L))
    contig.add_single_site(id="sweep", coordinate=SWEEP_SITE)
    # A frequency condition makes stdpopsim save the state before the draw,
    # so restore() can restart from it. "Not lost one generation after the
    # draw" is also applied by sweep_step() (with the phase still 0), so it
    # adds nothing.
    events = [stdpopsim.DrawMutation(time=t_mut, single_site_id="sweep",
                                     population="pop_0"),
              stdpopsim.ConditionOnAlleleFrequency(
                  start_time=t_mut - 1, end_time=t_mut - 1,
                  single_site_id="sweep", population="pop_0", op=">",
                  allele_frequency=0.0)]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        stdpopsim.get_engine("slim").simulate(
            model, contig, {"pop_0": N_GENOMES // 2}, extended_events=events,
            slim_script=True, slim_burn_in=t_mut / n + 0.5, seed=1)
    return buf.getvalue()


def sweep_script(base, f0, s_hom, final, trees_file):
    # SLiM 5 calls genomes haplosomes; this is the script's only use.
    s = base.replace("sample(pop.genomes, 1)", "sample(pop.haplosomes, 1)")
    s = s.replace('defineConstant("trees_file", ',
                  'defineConstant("trees_file_unused", ', 1)
    head = (f'    defineConstant("trees_file", "{trees_file}");\n'
            f'    defineConstant("F0", {float(f0)!r});\n'
            f'    defineConstant("S_HOM", {float(s_hom)!r});\n'
            f'    defineConstant("FINAL", {float(final)!r});\n')
    i = s.index('    defineConstant("pop_names"')
    return s[:i] + head + s[i:] + SWEEP_LOGIC


def stats(ts):
    """(total branch length, diversity) per base pair, in the tree
    sequence's time units."""
    return (ts.segregating_sites(mode="branch"), ts.diversity(mode="branch"))


def run_slim(job):
    base, n, rho, f0, s_hom, final, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        trees = os.path.join(tmp, "out.trees")
        script = os.path.join(tmp, "sweep.slim")
        with open(script, "w") as fh:
            fh.write(sweep_script(base, f0, s_hom, final, trees))
        subprocess.run(["slim", "-s", str(seed), script], check=True,
                       capture_output=True)
        ts = tskit.load(trees)
    # A run that reaches the present before the sweep ends is not used.
    if ts.metadata["SLiM"]["user_metadata"].get("sweep_done") != [1]:
        return None
    ts = pyslim.recapitate(ts, ancestral_Ne=n,
                           recombination_rate=rho / (4 * n * L),
                           random_seed=seed)
    rng = np.random.default_rng(seed)
    inds = rng.choice(pyslim.individuals_alive_at(ts, 0), N_GENOMES // 2,
                      replace=False)
    return stats(ts.simplify(np.concatenate([ts.individual(i).nodes
                                             for i in inds])))


def run_discoal(job):
    args, n, reps, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "d.trees")
        subprocess.run([DISCOAL, str(N_GENOMES), str(reps), str(L)] + args +
                       ["-N", str(n), "-d", str(seed), str(seed + 1),
                        "-ts", prefix], check=True, capture_output=True)
        # discoal node times are in units of 2N generations.
        return [tuple(2 * n * v for v in
                      stats(tskit.load(os.path.join(tmp, f"d_rep{k}.trees"))))
                for k in range(1, reps + 1)]


def run_msprime(job):
    n, rho, alpha, reps, seed = job
    sweep = msprime.SweepGenicSelection(
        position=SWEEP_SITE, start_frequency=1 / (2 * n),
        end_frequency=1 - 1 / (2 * n), s=alpha / n, dt=1e-6)
    return [stats(ts) for ts in msprime.sim_ancestry(
        samples=N_GENOMES // 2, population_size=n, sequence_length=L,
        recombination_rate=rho / (4 * n * L),
        model=[sweep, msprime.StandardCoalescent()],
        num_replicates=reps, random_seed=seed)]


def mean_se(vals):
    a = np.array(vals)
    return a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a)), len(a)


# name: (f0, final frequency, alpha, rho, generations before the present at
# which SLiM draws the mutation, in units of N). f0 = 1 with alpha = 0 is a
# neutral fixation (-wn).
CASES = {
    "hard": (0.0, 1.0, 100, 20, 6),
    "soft f0=0.02": (0.02, 1.0, 100, 20, 6),
    "soft f0=0.05": (0.05, 1.0, 100, 20, 6),
    "soft f0=0.1": (0.1, 1.0, 100, 20, 6),
    "soft f0=0.2": (0.2, 1.0, 100, 20, 6),
    "soft f0=0.05 rho=100": (0.05, 1.0, 100, 100, 6),
    "partial soft f0=0.05 c=0.7": (0.05, 0.7, 100, 20, 6),
    "neutral fixation (-wn)": (1.0, 1.0, 0, 20, 20),
}


def discoal_args(f0, final, alpha, rho):
    """discoal's arguments for a case (sweep ending at the present)."""
    args = ["-t", "0", "-r", str(rho), "-x", "0.5", "-a", str(alpha)]
    if alpha == 0:
        return args + ["-wn", "0"]
    args += ["-ws", "0"]
    if 0 < f0 < 1:
        args += ["-f", str(f0)]
    if final < 1:
        args += ["-c", str(final), "-h"]
    return args


def row(name, ref, ref_se, val, val_se, extra=""):
    """One line of the table: a reference (SLiM or msprime) and discoal."""
    z = (val - ref) / np.hypot(ref_se, val_se)
    return (f"{name:28s} {ref[0]:8.0f}({ref_se[0]:4.0f}) {val[0]:8.0f}({val_se[0]:4.0f})"
            f" {100 * (val[0] / ref[0] - 1):+6.1f}% {z[0]:+6.1f}"
            f"   {ref[1]:7.0f}({ref_se[1]:3.0f}) {val[1]:7.0f}({val_se[1]:3.0f})"
            f" {z[1]:+6.1f}{extra}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=100)
    ap.add_argument("--n", type=int, default=1000,
                    help="population size (diploid individuals)")
    ap.add_argument("--cases", default="",
                    help="comma-separated case names; default all")
    args = ap.parse_args()
    n, reps = args.n, args.reps
    wanted = set(filter(None, args.cases.split(",")))
    print(f"N = {n}; mean (standard error) over {reps} replicates; "
          "branch length and diversity in generations")
    print(f"{'case':28s} {'reference':>15s} {'discoal':>15s} {'diff':>7s} {'z':>6s}"
          f"   {'ref. div':>13s} {'discoal div':>13s} {'z':>6s}"
          "   (reference: SLiM, or msprime where named)", flush=True)
    with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
        for name, (f0, final, alpha, rho, t_mut) in CASES.items():
            if wanted and name not in wanted:
                continue
            base = stdpopsim_script(n, rho, t_mut * n)
            jobs = [(base, n, rho, f0, alpha / n, final, 1000 + k)
                    for k in range(int(reps * 1.05))]
            sl = [v for v in pool.map(run_slim, jobs) if v is not None][:reps]
            chunks = max(1, reps // 25)
            dc = [v for c in pool.map(run_discoal,
                                      [(discoal_args(f0, final, alpha, rho), n,
                                        reps // chunks, 7000 + 10 * k)
                                       for k in range(chunks)]) for v in c]
            (ms, ss, ns), (md, sd, _) = mean_se(sl), mean_se(dc)
            print(row(name, ms, ss, md, sd, f"   [SLiM n={ns}]"), flush=True)
            if name == "hard":
                mm, sm, _ = mean_se([v for c in pool.map(
                    run_msprime, [(n, rho, alpha, 20, 100 + k)
                                  for k in range(max(1, reps // 20))])
                    for v in c])
                print(row("  discoal against msprime", mm, sm, md, sd),
                      flush=True)

if __name__ == "__main__":
    main()
