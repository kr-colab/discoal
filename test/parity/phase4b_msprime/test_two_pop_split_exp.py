"""test_two_pop_split_exp.py — two-population split + EXP growth parity.

Configuration: pop0 + pop1, split (msprime POPULATION_SPLIT) at time
t_split = 1.0 (CLI). Pop0 experiences exponential growth from t=0 backward to
t_growth_end_cli = 0.3, then constant before that. No migration.

discoal CLI: -p 2 n n -ed t_split 0 1 -eg t_growth_end_cli 0 alpha
msprime: Demography with pop_split + add_population_parameters_change.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import msprime
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from parity_utils import (
    ACTIVE_CONVENTION,
    DEFAULT_NE,
    collect_stats,
    discoal_alpha_to_msp_growth,
    discoal_rho_to_msp_r,
    discoal_theta_to_msp_mu,
    discoal_time_to_msp_gen,
    msprime_to_ms,
    parse_ms,
    run_discoal,
    validate_active_convention,
)


def run_msprime_two_pop_split_exp(
    n_per_pop: int, L: int, theta: float, rho: float,
    eg_time_cli: float, eg_alpha: float, t_split_cli: float,
    Ne: int, reps: int,
):
    if ACTIVE_CONVENTION is None:
        raise RuntimeError("ACTIVE_CONVENTION is None. Run parity_utils.py first.")
    if n_per_pop % ACTIVE_CONVENTION.ploidy != 0:
        raise ValueError(
            f"n_per_pop={n_per_pop} must be divisible by ploidy={ACTIVE_CONVENTION.ploidy}"
        )
    n_individuals_per_pop = n_per_pop // ACTIVE_CONVENTION.ploidy
    mu = discoal_theta_to_msp_mu(theta, Ne, L)
    r = discoal_rho_to_msp_r(rho, Ne, L)
    t_growth_gen = discoal_time_to_msp_gen(eg_time_cli, Ne)
    g = discoal_alpha_to_msp_growth(eg_alpha, Ne)
    t_split_gen = discoal_time_to_msp_gen(t_split_cli, Ne)
    pop_size = Ne * ACTIVE_CONVENTION.population_size_factor

    demography = msprime.Demography()
    demography.add_population(name="pop0", initial_size=pop_size, growth_rate=0)
    demography.add_population(name="pop1", initial_size=pop_size, growth_rate=0)
    demography.add_population(name="ancestral", initial_size=pop_size, growth_rate=0)
    demography.add_population_parameters_change(
        time=t_growth_gen, population="pop0", growth_rate=g,
    )
    demography.add_population_split(
        time=t_split_gen, derived=["pop0", "pop1"], ancestral="ancestral",
    )

    out = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples={"pop0": n_individuals_per_pop, "pop1": n_individuals_per_pop},
            sequence_length=L,
            recombination_rate=r,
            demography=demography,
            ploidy=ACTIVE_CONVENTION.ploidy,
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


def main(reps: int = 1000):
    print("=== Phase 4b: two-pop split + EXP parity ===")
    print(f"reps per simulator = {reps}")
    validate_active_convention(reps=200)

    n_per_pop = 6  # must be divisible by ploidy
    L = 10000
    theta = 10.0
    rho = 10.0
    Ne = DEFAULT_NE
    eg_time_cli = 0.3
    eg_alpha = 50.0
    t_split_cli = 1.0  # both pops share an ancestor at this time

    discoal_args = [
        str(2 * n_per_pop), str(reps), str(L),
        "-t", str(theta), "-r", str(rho),
        "-p", "2", str(n_per_pop), str(n_per_pop),
        "-eg", str(eg_time_cli), "0", str(eg_alpha),
        "-ed", str(t_split_cli), "0", "1",
        "-d", "12345", "67890",
    ]
    print(f"discoal: {' '.join(['./build/discoal'] + discoal_args)}")
    t0 = time.time()
    discoal_out = run_discoal(discoal_args)
    print(f"  discoal: {time.time()-t0:.1f}s")
    discoal_reps = parse_ms(discoal_out)

    print("msprime: matched two-pop split + exp growth")
    t0 = time.time()
    msp_reps = run_msprime_two_pop_split_exp(
        n_per_pop, L, theta, rho, eg_time_cli, eg_alpha, t_split_cli, Ne, reps,
    )
    print(f"  msprime: {time.time()-t0:.1f}s")

    s_discoal = collect_stats(discoal_reps, 2 * n_per_pop)
    s_msp = collect_stats(msp_reps, 2 * n_per_pop)

    print()
    print(f"discoal: mean ss={s_discoal['ss'].mean():.2f}, mean pi={s_discoal['pi'].mean():.4f}, "
          f"mean tajD={s_discoal['td'].mean():.4f}")
    print(f"msprime: mean ss={s_msp['ss'].mean():.2f}, mean pi={s_msp['pi'].mean():.4f}, "
          f"mean tajD={s_msp['td'].mean():.4f}")

    comparisons = []
    for stat_name in ("ss", "pi", "td"):
        D, p = stats.ks_2samp(s_discoal[stat_name], s_msp[stat_name])
        comparisons.append((f"two_pop_split_exp_{stat_name}", D, p))

    n_comp = len(comparisons)
    alpha = 0.01 / n_comp
    print()
    print(f"{n_comp} comparisons, Bonferroni alpha = {alpha:.2e}")
    rejected = 0
    for name, D, p in comparisons:
        flag = "  **" if p < alpha else ""
        if p < alpha:
            rejected += 1
        print(f"{name:<32} {D:>10.4f} {p:>10.2e}{flag}")

    print()
    if rejected == 0:
        print(f"PASS: all {n_comp} comparisons within Bonferroni-corrected p > {alpha:.2e}")
        sys.exit(0)
    else:
        print(f"FAIL: {rejected}/{n_comp} comparisons reject equality at Bonferroni p < {alpha:.2e}")
        sys.exit(1)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    main(reps=reps)
