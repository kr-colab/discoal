"""test_issue82_fixture.py — multi-window migration parity (closes issue #82).

Loads config_examples/demes_example.demes.yaml in both discoal (via the
-Y flag with the wrapper YAML config) and msprime (via Demography.from_demes),
runs N replicates each at matched parameters, and compares summary statistics
via Bonferroni-corrected Kolmogorov-Smirnov tests at alpha = 0.01.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import demes
import msprime
import numpy as np
from scipy import stats

# Reuse Phase 4b parity utilities
sys.path.insert(0, str(Path(__file__).parent.parent / "phase4b_msprime"))
from parity_utils import (
    DEFAULT_NE,
    Convention,
    collect_stats,
    msprime_to_ms,
    parse_ms,
)


# Use Convention B (ploidy=2, popsize=Ne) — the validated convention from
# Phase 4b. We work around the ploidy=2 sample-count divisibility issue by
# passing per-population SampleSets with ploidy=1, which makes msprime treat
# each sample as a single haploid genome regardless of the global ploidy.
# Coalescent rate and migration scaling continue to use the global ploidy=2.
PHASE7_CONVENTION = Convention(
    name="B: ploidy=2, popsize=Ne (Phase 7 issue #82 test)",
    ploidy=2, population_size_factor=1.0,
)


REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
DISCOAL = REPO_ROOT / "build" / "discoal"
WRAPPER_YAML = REPO_ROOT / "config_examples" / "demes_example.yaml"
DEMES_FILE = REPO_ROOT / "config_examples" / "demes_example.demes.yaml"


def make_temp_yaml(reps: int) -> Path:
    """Create a temp YAML config that overrides num_replicates."""
    base = WRAPPER_YAML.read_text()
    # Replace the line `  num_replicates: 10` with the requested reps.
    overridden = base.replace("num_replicates: 10", f"num_replicates: {reps}")
    if "num_replicates: " + str(reps) not in overridden:
        raise RuntimeError("failed to override num_replicates in YAML")
    tmp = Path(tempfile.mkstemp(suffix=".yaml", prefix="phase7_")[1])
    tmp.write_text(overridden)
    return tmp


def run_discoal_via_yaml(yaml_path: Path) -> str:
    """Run discoal with -Y <yaml>, return stdout."""
    proc = subprocess.run(
        [str(DISCOAL), "-Y", str(yaml_path)],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"discoal exited {proc.returncode}: {proc.stderr[-500:]}")
    return proc.stdout


def run_msprime_from_demes(reps: int, n_per_pop: dict, theta: float, rho: float, L: int):
    """Run msprime with Demography.from_demes against the issue 82 fixture.

    n_per_pop is {pop_name: n_haploid_samples}.
    Uses Convention B: ploidy=2, popsize unchanged. Per-population SampleSets
    declared ploidy=1 to extract haploid samples without requiring even sample
    counts (works around the [10, 5, 5] divisibility issue while keeping
    Convention B's validated coalescent semantics).
    """
    graph = demes.load(str(DEMES_FILE))
    demography = msprime.Demography.from_demes(graph)

    # Apply the population_size_factor (1.0 under Convention B; identity).
    if PHASE7_CONVENTION.population_size_factor != 1.0:
        for pop in demography.populations:
            pop.initial_size = pop.initial_size * PHASE7_CONVENTION.population_size_factor

    # Discoal effective N for the fixture: reference N is the present-day pop A size.
    # From the fixture: A's present epoch has start_size 10000. So Ne_discoal = 10000.
    # Convert theta and rho to per-gen per-site rates using this Ne.
    ne_discoal = 10000  # present-day pop A
    mu = theta / (4.0 * ne_discoal * L)
    r = rho / (4.0 * ne_discoal * L)

    samples = [
        msprime.SampleSet(num_samples=n, population=pop_name, ploidy=1)
        for pop_name, n in n_per_pop.items()
    ]

    out = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples=samples,
            sequence_length=L,
            recombination_rate=r,
            demography=demography,
            ploidy=PHASE7_CONVENTION.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(
            ts, rate=mu,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            random_seed=seed,
        )
        out.append(msprime_to_ms(ts))
    return out


def main(reps: int = 500):
    print(f"=== Phase 7: issue #82 fixture parity (REPS={reps}) ===")
    print(f"Convention: {PHASE7_CONVENTION.name}")

    # discoal run
    yaml_path = make_temp_yaml(reps)
    print(f"discoal: -Y {yaml_path} (overridden num_replicates={reps})")
    t0 = time.time()
    discoal_out = run_discoal_via_yaml(yaml_path)
    print(f"  discoal: {time.time()-t0:.1f}s")
    discoal_reps = parse_ms(discoal_out)
    yaml_path.unlink()  # cleanup

    # msprime run with matched parameters
    n_per_pop = {"A": 10, "B": 5, "C": 5}  # haploid samples per fixture YAML
    L = 50000
    theta = 20.0
    rho = 15.0

    print("msprime: Demography.from_demes against the fixture")
    t0 = time.time()
    msp_reps = run_msprime_from_demes(reps, n_per_pop, theta, rho, L)
    print(f"  msprime: {time.time()-t0:.1f}s")

    n_total = 20  # 10 + 5 + 5
    s_discoal = collect_stats(discoal_reps, n_total)
    s_msp = collect_stats(msp_reps, n_total)

    print()
    print(f"discoal: mean ss={s_discoal['ss'].mean():.2f}, mean pi={s_discoal['pi'].mean():.4f}, "
          f"mean tajD={s_discoal['td'].mean():.4f}")
    print(f"msprime: mean ss={s_msp['ss'].mean():.2f}, mean pi={s_msp['pi'].mean():.4f}, "
          f"mean tajD={s_msp['td'].mean():.4f}")

    # Bonferroni-corrected KS comparisons
    comparisons = []
    for stat_name in ("ss", "pi", "td", "wtheta", "hapdiv", "nhap"):
        D, p = stats.ks_2samp(s_discoal[stat_name], s_msp[stat_name])
        comparisons.append((f"issue82_{stat_name}", D, p))

    n_comp = len(comparisons)
    bonf_alpha = 0.01 / n_comp
    print()
    print(f"{n_comp} comparisons, Bonferroni alpha = {bonf_alpha:.2e}")
    print(f"{'comparison':<24} {'D':>10} {'p':>10} {'reject?'}")
    rejected = 0
    for name, D, p in comparisons:
        flag = "  **" if p < bonf_alpha else ""
        if p < bonf_alpha:
            rejected += 1
        print(f"{name:<24} {D:>10.4f} {p:>10.2e}{flag}")

    print()
    if rejected == 0:
        print(f"PASS: all {n_comp} comparisons within Bonferroni-corrected p > {bonf_alpha:.2e}")
        print()
        print("Issue #82 closed: discoal under the rewritten importer matches")
        print("msprime running the same demes graph at statistical parity.")
        sys.exit(0)
    else:
        print(f"FAIL: {rejected}/{n_comp} comparisons reject equality at Bonferroni p < {bonf_alpha:.2e}")
        sys.exit(1)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    main(reps=reps)
