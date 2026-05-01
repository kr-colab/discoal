"""test_single_pop_exp.py — single-population exponential growth parity.

Run discoal with `-eg time 0 alpha` and msprime with a matched
add_population_parameters_change at the converted generation time
and per-generation growth rate. Compare summary statistics via
Bonferroni-corrected Kolmogorov-Smirnov tests.
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


def run_msprime_single_pop_exp(
    n: int, L: int, theta: float, rho: float,
    eg_time_cli: float, eg_alpha: float, Ne: int, reps: int,
):
    """Run msprime with a single population that experiences exponential growth
    starting at the converted generation time, with per-generation growth rate
    derived from discoal's 4N-scaled alpha. Uses ACTIVE_CONVENTION (set in
    parity_utils.py after the constant-parity trial) for ploidy/popsize."""
    if ACTIVE_CONVENTION is None:
        raise RuntimeError("ACTIVE_CONVENTION is None. Run parity_utils.py first.")
    mu = discoal_theta_to_msp_mu(theta, Ne, L)
    r = discoal_rho_to_msp_r(rho, Ne, L)
    t_gen = discoal_time_to_msp_gen(eg_time_cli, Ne)
    g = discoal_alpha_to_msp_growth(eg_alpha, Ne)
    pop_size = Ne * ACTIVE_CONVENTION.population_size_factor

    demography = msprime.Demography()
    demography.add_population(name="pop0", initial_size=pop_size, growth_rate=0)
    demography.add_population_parameters_change(
        time=t_gen, population="pop0", growth_rate=g,
    )

    # discoal's `n` is the haploid sample count. msprime's `samples={pop: k}`
    # with ploidy=p produces k*p haploid samples. Pass n // ploidy individuals
    # so the total haploid sample count matches discoal.
    n_individuals = n // ACTIVE_CONVENTION.ploidy
    if n_individuals * ACTIVE_CONVENTION.ploidy != n:
        raise RuntimeError(
            f"n={n} not divisible by ploidy={ACTIVE_CONVENTION.ploidy}; "
            "cannot match discoal's haploid sample count exactly."
        )

    out = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples={"pop0": n_individuals},
            sequence_length=L,
            recombination_rate=r,
            demography=demography,
            ploidy=ACTIVE_CONVENTION.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(
            ts, rate=mu, random_seed=seed,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
        )
        out.append(msprime_to_ms(ts))
    return out


def main(reps: int = 1000):
    print("=== Phase 4b: single-pop EXP growth parity ===")
    print(f"reps per simulator = {reps}")
    validate_active_convention(reps=200)

    n = 10
    L = 10000
    theta = 10.0
    rho = 10.0
    Ne = DEFAULT_NE
    eg_time_cli = 0.5  # 2N units; -> 4N=1.0 internal -> 4Ne gen
    eg_alpha = 50.0    # 4N-scaled

    discoal_args = [
        str(n), str(reps), str(L),
        "-t", str(theta), "-r", str(rho),
        "-eg", str(eg_time_cli), "0", str(eg_alpha),
        "-d", "12345", "67890",
    ]
    print(f"discoal: {' '.join(['./build/discoal'] + discoal_args)}")
    t0 = time.time()
    discoal_out = run_discoal(discoal_args)
    print(f"  discoal: {time.time()-t0:.1f}s")
    discoal_reps = parse_ms(discoal_out)

    print("msprime: matched single-pop exp growth")
    t0 = time.time()
    msp_reps = run_msprime_single_pop_exp(n, L, theta, rho, eg_time_cli, eg_alpha, Ne, reps)
    print(f"  msprime: {time.time()-t0:.1f}s")

    s_discoal = collect_stats(discoal_reps, n)
    s_msp = collect_stats(msp_reps, n)

    print()
    print(f"discoal: mean ss={s_discoal['ss'].mean():.2f}, mean pi={s_discoal['pi'].mean():.4f}, "
          f"mean tajD={s_discoal['td'].mean():.4f}")
    print(f"msprime: mean ss={s_msp['ss'].mean():.2f}, mean pi={s_msp['pi'].mean():.4f}, "
          f"mean tajD={s_msp['td'].mean():.4f}")

    # KS tests across statistics
    comparisons = []
    for stat_name in ("ss", "pi", "td"):
        D, p = stats.ks_2samp(s_discoal[stat_name], s_msp[stat_name])
        comparisons.append((f"single_pop_exp_{stat_name}", D, p))

    # Per-bin chi-squared on SFS (folded)
    sfs_d = s_discoal["sfs"].sum(axis=0)
    sfs_m = s_msp["sfs"].sum(axis=0)
    if sfs_d.sum() > 0 and sfs_m.sum() > 0 and len(sfs_d) == len(sfs_m):
        scale = sfs_m.sum() / sfs_d.sum()
        expected = sfs_d * scale
        mask = expected > 5
        if mask.sum() >= 2:
            chi2 = float(np.sum((sfs_m[mask] - expected[mask]) ** 2 / expected[mask]))
            df = mask.sum() - 1
            p_chi = float(1 - stats.chi2.cdf(chi2, df))
            comparisons.append(("single_pop_exp_sfs_chi2", chi2, p_chi))

    n_comp = len(comparisons)
    alpha = 0.01 / n_comp
    print()
    print(f"{n_comp} comparisons, Bonferroni alpha = {alpha:.2e}")
    print(f"{'comparison':<32} {'D/chi2':>10} {'p':>10} {'reject?'}")
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
