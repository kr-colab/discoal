"""test_single_pop_exp.py — single-population exponential growth parity sweep.

Sweeps a grid of (alpha, eg_time_cli) configurations, runs REPS replicates of
discoal -eg and msprime add_population_parameters_change for each cell, and
runs Bonferroni-corrected KS tests across multiple summary statistics:
ss, pi, Tajima's D, Watterson's theta, haplotype diversity, number of
haplotypes, and per-bin folded SFS counts.

Reports overall PASS/FAIL plus the full per-cell-per-stat table.
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
    if ACTIVE_CONVENTION is None:
        raise RuntimeError("ACTIVE_CONVENTION is None. Run parity_utils.py first.")
    if n % ACTIVE_CONVENTION.ploidy != 0:
        raise ValueError(f"n={n} must be divisible by ploidy={ACTIVE_CONVENTION.ploidy}")
    n_individuals = n // ACTIVE_CONVENTION.ploidy
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
            ts, rate=mu,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            random_seed=seed,
        )
        out.append(msprime_to_ms(ts))
    return out


def run_one_cell(n: int, L: int, theta: float, rho: float,
                 eg_time_cli: float, eg_alpha: float, Ne: int, reps: int):
    """Run discoal and msprime for one (alpha, eg_time) cell. Return collected stats."""
    discoal_args = [
        str(n), str(reps), str(L),
        "-t", str(theta), "-r", str(rho),
        "-eg", str(eg_time_cli), "0", str(eg_alpha),
        "-d", "12345", "67890",
    ]
    discoal_out = run_discoal(discoal_args)
    discoal_reps = parse_ms(discoal_out)
    msp_reps = run_msprime_single_pop_exp(n, L, theta, rho, eg_time_cli, eg_alpha, Ne, reps)
    return collect_stats(discoal_reps, n), collect_stats(msp_reps, n)


def main(reps: int = 5000):
    print(f"=== Phase 4b: single-pop EXP growth parity sweep (REPS={reps}) ===")
    validate_active_convention(reps=200)

    n = 10
    L = 10000
    theta = 10.0
    rho = 10.0
    Ne = DEFAULT_NE

    alphas = [50.0, 200.0, 1000.0]
    eg_times = [0.1, 0.5, 1.0]

    # Run the grid, collect all comparison results
    comparisons = []  # (cell_label, stat_name, D_or_chi2, p)
    cell_means = []   # (cell_label, sim, stat, mean) for the diagnostic table

    for alpha in alphas:
        for t_cli in eg_times:
            cell = f"a{int(alpha)}_t{t_cli}"
            t0 = time.time()
            s_d, s_m = run_one_cell(n, L, theta, rho, t_cli, alpha, Ne, reps)
            elapsed = time.time() - t0
            print(f"  cell {cell}: {elapsed:.1f}s")

            # Scalar stats
            for stat_name in ("ss", "pi", "td", "wtheta", "hapdiv", "nhap"):
                D, p = stats.ks_2samp(s_d[stat_name], s_m[stat_name])
                comparisons.append((cell, stat_name, D, p))

            # Per-bin SFS
            n_bins = s_d["sfs"].shape[1]
            for b in range(n_bins):
                D, p = stats.ks_2samp(s_d["sfs"][:, b], s_m["sfs"][:, b])
                comparisons.append((cell, f"sfs_bin{b+1}", D, p))

            for stat_name in ("ss", "pi", "td", "wtheta", "hapdiv", "nhap"):
                cell_means.append((cell, "discoal", stat_name, s_d[stat_name].mean()))
                cell_means.append((cell, "msprime", stat_name, s_m[stat_name].mean()))

    n_comp = len(comparisons)
    bonf_alpha = 0.01 / n_comp
    print()
    print(f"{n_comp} comparisons, Bonferroni alpha = {bonf_alpha:.2e}")
    print()

    # Print mean comparisons (compact)
    print(f"{'cell':<14} {'stat':<10} {'discoal':>12} {'msprime':>12} {'rel_err':>10}")
    for cell in [f"a{int(a)}_t{t}" for a in alphas for t in eg_times]:
        for stat in ("ss", "pi", "td", "wtheta", "hapdiv", "nhap"):
            d_mean = next((m for c, sim, s, m in cell_means if c == cell and sim == "discoal" and s == stat), 0.0)
            m_mean = next((m for c, sim, s, m in cell_means if c == cell and sim == "msprime" and s == stat), 0.0)
            denom = max(abs(d_mean), abs(m_mean), 1e-12)
            rel_err = abs(d_mean - m_mean) / denom
            print(f"{cell:<14} {stat:<10} {d_mean:>12.4f} {m_mean:>12.4f} {rel_err:>10.3f}")
        print()

    # Print per-comparison table, flagging rejections
    print(f"{'cell':<14} {'stat':<14} {'D/chi2':>10} {'p':>10} {'reject?'}")
    rejected = 0
    for cell, stat_name, D, p in comparisons:
        flag = "  **" if p < bonf_alpha else ""
        if p < bonf_alpha:
            rejected += 1
        print(f"{cell:<14} {stat_name:<14} {D:>10.4f} {p:>10.2e}{flag}")

    print()
    if rejected == 0:
        print(f"PASS: all {n_comp} comparisons within Bonferroni-corrected p > {bonf_alpha:.2e}")
        sys.exit(0)
    else:
        print(f"FAIL: {rejected}/{n_comp} comparisons reject equality at Bonferroni p < {bonf_alpha:.2e}")
        sys.exit(1)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    main(reps=reps)
