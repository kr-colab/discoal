"""parity_utils.py — discoal vs msprime parameter conversion + statistic extraction.

discoal and msprime use different time and rate conventions. The mapping is
encoded in a Convention struct that captures msprime's ploidy and a multiplier
for the effective population size. check_constant_parity() trials three
plausible conventions against a no-demography baseline and reports which one
agrees with discoal within 10% relative error on mean pi.

The validated convention is locked into ACTIVE_CONVENTION at module level.
The EXP-parity tests in Tasks 2-3 use ACTIVE_CONVENTION exclusively.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import msprime


REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
DISCOAL = REPO_ROOT / "build" / "discoal"

# Default discoal effective population size (matches discoal_multipop.c default
# EFFECTIVE_POPN_SIZE = 1_000_000).
DEFAULT_NE = 1_000_000


@dataclass(frozen=True)
class Convention:
    """A ploidy / population_size combination for matching discoal in msprime.

    For each convention, time conversion is fixed: discoal CLI time (in 2N units)
    becomes msprime generations as `t_cli * 4 * N_e`, and discoal alpha
    (4N-scaled) becomes per-generation growth as `alpha / (4 * N_e)`. The
    variable is what to pass msprime as ploidy and population_size — these
    determine msprime's per-pair coalescent rate per generation, which must
    match discoal's coalescent timescale to produce the same pi distribution.
    """
    name: str
    ploidy: int
    population_size_factor: float  # passed to msprime as ploidy * Ne * factor


CONVENTION_A = Convention(name="A: ploidy=1, popsize=Ne", ploidy=1, population_size_factor=1.0)
CONVENTION_B = Convention(name="B: ploidy=2, popsize=Ne", ploidy=2, population_size_factor=1.0)
CONVENTION_C = Convention(name="C: ploidy=1, popsize=2*Ne", ploidy=1, population_size_factor=2.0)
CONVENTIONS = [CONVENTION_A, CONVENTION_B, CONVENTION_C]

# Set this to the validated Convention after running check_constant_parity once.
# Default value is None; tests in Tasks 2-3 call validate_active_convention() at
# startup and raise if it's None or doesn't match within tolerance.
ACTIVE_CONVENTION = CONVENTION_B


@dataclass
class MsHaplotypes:
    """One simulation replicate parsed from ms-format output."""
    positions: np.ndarray  # in [0, 1]
    haplotypes: np.ndarray  # shape (n, S), 0/1


def parse_ms(text: str) -> List[MsHaplotypes]:
    """Parse discoal/ms-format text into a list of replicate haplotype matrices.

    Skips DEBUG lines on stderr; works on stdout. Tolerates the discoal
    command-line echo on line 1 and seeds line on line 2.
    """
    reps: List[MsHaplotypes] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("//"):
            i += 1
            assert i < len(lines), "truncated ms output"
            m = re.match(r"segsites:\s*(\d+)", lines[i].strip())
            assert m, f"expected segsites line, got: {lines[i]}"
            S = int(m.group(1))
            if S == 0:
                reps.append(MsHaplotypes(np.zeros(0), np.zeros((0, 0))))
                i += 1
                continue
            i += 1
            assert i < len(lines), "truncated ms output"
            m = re.match(r"positions:\s*(.*)", lines[i].strip())
            assert m, f"expected positions line, got: {lines[i]}"
            pos = np.array([float(x) for x in m.group(1).split()])
            assert pos.shape == (S,), f"expected {S} positions, got {pos.shape}"
            i += 1
            haps = []
            while i < len(lines):
                ln = lines[i].strip()
                if not ln or ln.startswith("//"):
                    break
                if all(c in "01" for c in ln):
                    haps.append([int(c) for c in ln])
                    i += 1
                else:
                    break
            haps_arr = np.array(haps, dtype=np.int8)
            reps.append(MsHaplotypes(pos, haps_arr))
        else:
            i += 1
    return reps


def msprime_to_ms(ts) -> MsHaplotypes:
    """Convert an msprime tree-sequence to ms-format MsHaplotypes."""
    n = ts.num_samples
    S = ts.num_sites
    if S == 0:
        return MsHaplotypes(np.zeros(0), np.zeros((n, 0), dtype=np.int8))
    pos = np.array([s.position for s in ts.sites()]) / ts.sequence_length
    haps = np.zeros((n, S), dtype=np.int8)
    for var in ts.variants():
        haps[:, var.site.id] = var.genotypes
    return MsHaplotypes(pos, haps)


def segsites(rep: MsHaplotypes) -> int:
    return rep.haplotypes.shape[1] if rep.haplotypes.ndim == 2 else 0


def pi(rep: MsHaplotypes) -> float:
    """Nucleotide diversity (mean pairwise differences per site)."""
    if rep.haplotypes.size == 0:
        return 0.0
    n, S = rep.haplotypes.shape
    if n < 2:
        return 0.0
    p = rep.haplotypes.mean(axis=0)
    h = 2 * p * (1 - p) * n / (n - 1)
    return float(h.sum())


def tajima_D(rep: MsHaplotypes) -> float:
    """Tajima's D. Returns 0.0 if undefined (no sites or n<2)."""
    if rep.haplotypes.size == 0:
        return 0.0
    n, S = rep.haplotypes.shape
    if S == 0 or n < 2:
        return 0.0
    a1 = sum(1.0 / k for k in range(1, n))
    a2 = sum(1.0 / (k * k) for k in range(1, n))
    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    var = e1 * S + e2 * S * (S - 1)
    if var <= 0:
        return 0.0
    pi_hat = pi(rep)
    theta_w = S / a1
    return float((pi_hat - theta_w) / np.sqrt(var))


def sfs(rep: MsHaplotypes) -> np.ndarray:
    """Folded SFS (no outgroup distinction). Bins are (1..n//2)."""
    if rep.haplotypes.size == 0:
        return np.zeros(0)
    n, S = rep.haplotypes.shape
    if S == 0:
        return np.zeros(n // 2)
    counts = rep.haplotypes.sum(axis=0)
    folded = np.minimum(counts, n - counts)
    return np.bincount(folded, minlength=n // 2 + 1)[1:n // 2 + 1].astype(float)


def watterson_theta(rep: MsHaplotypes) -> float:
    """Watterson's theta = S / sum_{k=1}^{n-1}(1/k)."""
    if rep.haplotypes.size == 0:
        return 0.0
    n, S = rep.haplotypes.shape
    if S == 0 or n < 2:
        return 0.0
    a1 = sum(1.0 / k for k in range(1, n))
    return float(S / a1)


def haplotype_diversity(rep: MsHaplotypes) -> float:
    """Haplotype diversity h = (n / (n-1)) * (1 - sum p_i^2) where p_i is the
    frequency of haplotype i. Returns 0.0 if n < 2."""
    if rep.haplotypes.size == 0:
        return 0.0
    n, S = rep.haplotypes.shape
    if n < 2:
        return 0.0
    if S == 0:
        return 0.0  # all identical
    # Convert each row to a tuple, count occurrences
    rows = [tuple(row) for row in rep.haplotypes]
    counts = {}
    for r in rows:
        counts[r] = counts.get(r, 0) + 1
    p = np.array(list(counts.values()), dtype=float) / n
    return float(n / (n - 1) * (1.0 - np.sum(p * p)))


def num_haplotypes(rep: MsHaplotypes) -> int:
    """Number of distinct haplotype patterns."""
    if rep.haplotypes.size == 0:
        return 1
    rows = {tuple(row) for row in rep.haplotypes}
    return len(rows)


def collect_stats(reps: Iterable[MsHaplotypes], n: int) -> dict:
    """Compute per-replicate stats; return arrays."""
    ss_list, pi_list, td_list, wt_list, hd_list, nh_list, sfs_list = [], [], [], [], [], [], []
    n_bins = n // 2
    for r in reps:
        ss_list.append(segsites(r))
        pi_list.append(pi(r))
        td_list.append(tajima_D(r))
        wt_list.append(watterson_theta(r))
        hd_list.append(haplotype_diversity(r))
        nh_list.append(num_haplotypes(r))
        sfs_row = sfs(r)
        # Pad/truncate empty-replicate SFS (shape 0,) to length n_bins for stacking.
        if sfs_row.shape[0] != n_bins:
            sfs_row = np.zeros(n_bins)
        sfs_list.append(sfs_row)
    sfs_arr = np.array(sfs_list) if sfs_list else np.zeros((0, n_bins))
    return dict(
        ss=np.array(ss_list, dtype=float),
        pi=np.array(pi_list),
        td=np.array(td_list),
        wtheta=np.array(wt_list),
        hapdiv=np.array(hd_list),
        nhap=np.array(nh_list, dtype=float),
        sfs=sfs_arr,
    )


def run_discoal(args: List[str], capture_stderr: bool = False) -> str:
    """Run discoal with given args, return stdout."""
    if not DISCOAL.exists():
        raise RuntimeError(f"discoal not built: {DISCOAL}")
    proc = subprocess.run(
        [str(DISCOAL)] + args,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"discoal exited {proc.returncode}: {proc.stderr[-500:]}")
    return proc.stdout


def discoal_time_to_msp_gen(t_cli: float, Ne: int) -> float:
    """Convert discoal CLI time (2N units) to msprime generations.

    discoal CLI time is doubled internally (-> 4N units), and 4N units = 4*Ne
    generations.
    """
    return t_cli * 4 * Ne


def discoal_alpha_to_msp_growth(alpha: float, Ne: int) -> float:
    """Convert discoal alpha to msprime per-generation growth rate.

    discoal stores `alpha` as the rate_param of an exponential shape applied in
    discoal *internal* time units. From shapes.c:
        size(t_internal) = anchor_value * exp(-alpha * (t_internal - t0_internal))
    Discoal internal time runs at the standard pair-coalescent rate of 1.0 per
    unit (see neutralPhase: cRate = n*(n-1)/2 / sizeRatio), so 1 internal unit
    equals 2*Ne_diploid generations under the WF coalescent. Therefore:

        alpha * t_internal = g * t_gen
        g = alpha / (2 * Ne)

    Note: the discoal CLI multiplies the user-supplied event time by 2.0 to
    convert it to internal units, so the time conversion (CLI -> generations)
    has an extra factor of 2 (4*Ne). The growth rate does not; it is anchored
    in internal units directly.
    """
    return alpha / (2 * Ne)


def discoal_theta_to_msp_mu(theta: float, Ne: int, L: int) -> float:
    """Convert discoal theta (4*Ne*mu*L) to msprime per-generation per-site mutation rate."""
    return theta / (4 * Ne * L)


def discoal_rho_to_msp_r(rho: float, Ne: int, L: int) -> float:
    """Convert discoal rho (4*Ne*r*L) to msprime per-generation per-site recombination rate."""
    return rho / (4 * Ne * L)


def _msp_constant_parity_pi(n: int, L: int, theta: float, rho: float,
                              Ne: int, reps: int, conv: Convention) -> float:
    """Mean pi from msprime no-demography simulation under the given convention.

    `n` is the desired haploid sample count (matches discoal's first positional
    arg). msprime takes `samples` as the number of *individuals* of the given
    ploidy, so we pass `n // ploidy` to get exactly `n` haploid output samples.
    Requires n % ploidy == 0.
    """
    if n % conv.ploidy != 0:
        raise ValueError(
            f"n={n} must be divisible by ploidy={conv.ploidy} for sample-count parity"
        )
    n_individuals = n // conv.ploidy
    mu = discoal_theta_to_msp_mu(theta, Ne, L)
    r = discoal_rho_to_msp_r(rho, Ne, L)
    pop_size = Ne * conv.population_size_factor
    pis = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples=n_individuals,
            sequence_length=L,
            recombination_rate=r,
            population_size=pop_size,
            ploidy=conv.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(
            ts, rate=mu,
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            random_seed=seed,
        )
        pis.append(pi(msprime_to_ms(ts)))
    return float(np.mean(pis))


def check_constant_parity(n: int = 6, theta: float = 5.0, rho: float = 5.0,
                          L: int = 1000, Ne: int = DEFAULT_NE,
                          reps: int = 200, tolerance: float = 0.10,
                          ) -> Convention:
    """Trial all three conventions against a no-demography baseline.

    Returns the Convention whose mean pi agrees with discoal within tolerance.
    Raises RuntimeError if no convention matches.
    """
    discoal_args = [
        str(n), str(reps), str(L),
        "-t", str(theta), "-r", str(rho),
        "-d", "12345", "67890",
    ]
    discoal_out = run_discoal(discoal_args)
    discoal_reps = parse_ms(discoal_out)
    discoal_pi = float(np.mean([pi(r) for r in discoal_reps]))
    print(f"discoal mean pi = {discoal_pi:.4f} (n={n}, theta={theta}, rho={rho}, reps={reps})")

    results = []
    for conv in CONVENTIONS:
        msp_pi = _msp_constant_parity_pi(n, L, theta, rho, Ne, reps, conv)
        rel_err = abs(discoal_pi - msp_pi) / max(abs(discoal_pi), abs(msp_pi), 1e-12)
        marker = "  <-- best so far" if not results or rel_err < min(r[2] for r in results) else ""
        print(f"  {conv.name:<28} msprime mean pi = {msp_pi:.4f}, rel_err = {rel_err:.3f}{marker}")
        results.append((conv, msp_pi, rel_err))

    best = min(results, key=lambda x: x[2])
    if best[2] > tolerance:
        raise RuntimeError(
            f"No convention matches within tolerance={tolerance}. "
            f"Best was {best[0].name} with rel_err={best[2]:.3f}. "
            f"Either the parameter conversions in parity_utils.py are wrong, "
            f"or msprime's underlying coalescent semantics differ from discoal "
            f"in a way the conventions don't capture. Halt and investigate."
        )
    print(f"\nBest convention: {best[0].name} (rel_err = {best[2]:.3f})")
    return best[0]


def validate_active_convention(reps: int = 100, tolerance: float = 0.10) -> None:
    """Sanity-check that ACTIVE_CONVENTION is set and still produces matching pi.

    Called at the start of each EXP-parity test in Tasks 2-3. Halts if the
    constant-N check fails — guards against subtle conversion regressions.
    """
    if ACTIVE_CONVENTION is None:
        raise RuntimeError(
            "ACTIVE_CONVENTION is None. Run check_constant_parity() once to "
            "discover the right convention, then set ACTIVE_CONVENTION at the "
            "top of parity_utils.py."
        )
    discoal_args = [
        "6", str(reps), "1000",
        "-t", "5.0", "-r", "5.0",
        "-d", "12345", "67890",
    ]
    discoal_out = run_discoal(discoal_args)
    discoal_pi = float(np.mean([pi(r) for r in parse_ms(discoal_out)]))
    msp_pi = _msp_constant_parity_pi(6, 1000, 5.0, 5.0, DEFAULT_NE, reps, ACTIVE_CONVENTION)
    rel_err = abs(discoal_pi - msp_pi) / max(abs(discoal_pi), abs(msp_pi), 1e-12)
    if rel_err > tolerance:
        raise RuntimeError(
            f"validate_active_convention FAILED: rel_err = {rel_err:.3f}. "
            f"ACTIVE_CONVENTION = {ACTIVE_CONVENTION.name}. The conversion is "
            f"no longer valid; re-run check_constant_parity() to find the right one."
        )
    print(f"validate_active_convention PASS: rel_err = {rel_err:.3f} "
          f"(convention = {ACTIVE_CONVENTION.name})")


if __name__ == "__main__":
    best = check_constant_parity()
    print(f"\nTo lock in the best convention, edit parity_utils.py:")
    print(f"  ACTIVE_CONVENTION = {best!r}")
