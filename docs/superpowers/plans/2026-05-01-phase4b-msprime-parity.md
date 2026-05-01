# Phase 4b: msprime Parity Validation for SHAPE_EXPONENTIAL

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Validate that the `SHAPE_EXPONENTIAL` support added in Phase 4 produces statistically correct results by comparing discoal (`-eg`) output against msprime running the same demographic model. The test passes if Bonferroni-corrected Kolmogorov-Smirnov tests on summary statistics show no significant distributional differences at p > 0.01 across all comparisons.

**Architecture:** A Python harness in `test/parity/phase4b_msprime/` builds matched single-pop and two-pop EXP-growth configurations, runs both simulators with $10^4$ replicates each, parses ms-format output to extract per-replicate summary statistics (segregating sites, π, Tajima's D, SFS bins), and runs statistical comparisons. The discoal-internal-time vs msprime-generations conversion is handled in a small `parity_utils.py` module (CLI 2N units → internal 4N units → generations × $N_e$, alpha 4N-scaled → per-generation × $1/(4N_e)$).

**Tech Stack:** Python 3, `msprime 1.4.1`, `numpy 2.4.4`, `scipy 1.17.1`, bash.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §6.2 (msprime parity test design) + Phase 4 in §5.

**Deliverables:**
- `test/parity/phase4b_msprime/parity_utils.py` — parameter-conversion helpers and statistic extraction.
- `test/parity/phase4b_msprime/test_single_pop_exp.py` — single-pop EXP growth parity (Phase 4 deliverable per spec §5).
- `test/parity/phase4b_msprime/test_two_pop_split_exp.py` — two-pop split + EXP growth in one branch parity.
- `test/parity/phase4b_msprime/run_parity.sh` — orchestration script that runs both tests and emits PASS/FAIL.
- `test/parity/phase4b_msprime/results/.gitignore` — keeps raw outputs out of git.
- A committed analysis summary at `test/parity/phase4b_msprime/results/analysis.txt`.
- A spec addendum recording the Phase 4b outcome.
- Tag `phase4b-msprime-parity-complete` (or `phase4b-msprime-parity-fail` if the result is FAIL — the result determines the tag).

**Out of scope:**
- `SHAPE_LINEAR` parity (Phase 6).
- Sweep parity (Phase 5b after sweep wiring lands).
- Full demes-graph parity (Phase 7 — closes issue #82).

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push (push is a separate manual step at end).

---

## Conversion math (used in every task)

**discoal vs msprime are not the same simulator and have non-obvious unit conventions. Get this wrong and the parity tests will reject for spurious reasons.**

### Two independent conversion knobs

There are TWO independent choices that determine the conversion:

1. **discoal's internal time unit.** The docs say discoal uses 4N internal time units (with CLI time in 2N), but this needs to be empirically confirmed — it could be 2N if the codebase treats N as haploid. The critical relation is: `1 discoal-internal-time-unit = K * N_e` generations for some constant K ∈ {1, 2, 4}.

2. **msprime ploidy and the matched `population_size`.** msprime's per-pair coalescent rate is `1 / (ploidy * population_size)` per generation. For matching, we choose `(ploidy, population_size)` so the per-pair rate yields a coalescent timescale of `K * N_e` generations per discoal-internal-time-unit.

The combinations that *could* match (depending on K):

| If discoal `K` (gen per discoal-time-unit) is... | Use msprime ploidy | And `population_size` |
|---|---|---|
| 4 (the documented "4N convention") | 2 (diploid) | $N_e$ — gives expected pairwise coal time = $2 N_e$ generations per discoal-time-unit; then `t_gen = t_internal * 2 * N_e = t_cli * 4 * N_e` ✗ wait — the math here is inconsistent |
| 4 | 1 (haploid) | $2 N_e$ — gives expected pairwise coal time of $2 N_e$ generations | 

Actually the cleanest way to state it: msprime's coalescent rate for a pair of haploid samples is `1 / population_size` per generation, while for diploid it's `1 / (2 * population_size)` per generation. Either way, the *effective population size in haploid terms* is `ploidy * population_size`. We want msprime's effective haploid size to equal whatever discoal's `EFFECTIVE_POPN_SIZE` is (in haploid terms).

**Pragmatic approach:** treat conversion as an empirical question, not a derivation. Try the three plausible combinations of (ploidy, time_factor) until one passes the constant-parity sanity check at < 10% relative error on mean π. The combinations to try, in order:

| Try # | msprime `ploidy` | msprime `population_size` | `t_gen = t_cli * ?` | `alpha_per_gen = alpha_discoal / ?` |
|---|---|---|---|---|
| A | 1 | $N_e$ | $4 \cdot N_e$ | $4 \cdot N_e$ |
| B | 2 | $N_e$ | $4 \cdot N_e$ | $4 \cdot N_e$ |
| C | 1 | $2 \cdot N_e$ | $4 \cdot N_e$ | $4 \cdot N_e$ |

(All three keep `t_gen = t_cli * 4 * N_e` and `alpha_per_gen = alpha_discoal / (4 * N_e)` because the docs explicitly say "command line uses 2N, internal uses 4N". The variable is what to feed msprime as `population_size` × `ploidy`.)

### Mutation and recombination rates

discoal `-t theta` is $\theta = 4 N_e \mu L$ regardless of ploidy convention (it's the diversity-scaling rate). msprime `mutation_rate` is per-generation per-site. Convert:

```
mu_per_gen = theta / (4 * N_e * L)
```

This is independent of the ploidy choice because `theta = 4*N_e*mu*L` is the standard population-scaled mutation rate definition, and msprime expects raw per-generation `mu`.

Recombination: `r_per_gen = rho / (4 * N_e * L)` for the same reason.

### Empirical validation strategy

`parity_utils.py` includes `check_constant_parity()` which:

1. Runs discoal with `-t theta -r rho` (no demography) at $\sim 200$ replicates.
2. Runs msprime with each of the three trial conversions in the table above.
3. Reports `discoal_pi`, `msp_pi_A`, `msp_pi_B`, `msp_pi_C` and the pairwise relative errors.
4. **The implementer picks the conversion that minimizes the relative error** and updates `discoal_*_to_msp_*` to use that combination. The final implementation has only one conversion path; the trial code is for diagnostic purposes during the initial conversion-discovery phase, then can be simplified.

If NONE of the three combinations match within ~10%, something deeper is wrong (e.g., msprime's `recombination_rate` semantics differ, or there's a mutation-model difference). Halt and investigate before proceeding to EXP-parity tests.

### After the conversion is locked in

The parity tests in Tasks 2 and 3 use the *single* validated conversion (no trial logic). They re-run `check_constant_parity` at startup as a regression guard — if a future change breaks the conversion, the tests fail loudly at the sanity check rather than misreport an EXP-parity failure.

---

## Task 1: Parity infrastructure (`parity_utils.py`)

**Files:**
- Create: `test/parity/phase4b_msprime/parity_utils.py`
- Create: `test/parity/phase4b_msprime/__init__.py` (empty, makes the dir a package for test imports)
- Create: `test/parity/phase4b_msprime/results/.gitignore`

- [ ] **Step 1: Create the directory**

```bash
mkdir -p test/parity/phase4b_msprime/results
touch test/parity/phase4b_msprime/__init__.py
echo '*.ms' > test/parity/phase4b_msprime/results/.gitignore
echo '*.err' >> test/parity/phase4b_msprime/results/.gitignore
echo '*.npz' >> test/parity/phase4b_msprime/results/.gitignore
echo '*.csv' >> test/parity/phase4b_msprime/results/.gitignore
```

- [ ] **Step 2: Write `parity_utils.py`**

The library exposes a `Convention` named-tuple plus three pre-defined conventions A/B/C as described in the conversion math above. `check_constant_parity()` tries all three and prints the relative errors so the implementer can pick the right one. After picking, set `ACTIVE_CONVENTION` in the module to lock it in for the actual EXP tests.

```python
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
ACTIVE_CONVENTION: Optional[Convention] = None


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
            # Read segsites
            i += 1
            assert i < len(lines), "truncated ms output"
            m = re.match(r"segsites:\s*(\d+)", lines[i].strip())
            assert m, f"expected segsites line, got: {lines[i]}"
            S = int(m.group(1))
            if S == 0:
                # No variation; emit empty rep
                reps.append(MsHaplotypes(np.zeros(0), np.zeros((0, 0))))
                i += 1
                continue
            i += 1
            # Read positions
            assert i < len(lines), "truncated ms output"
            m = re.match(r"positions:\s*(.*)", lines[i].strip())
            assert m, f"expected positions line, got: {lines[i]}"
            pos = np.array([float(x) for x in m.group(1).split()])
            assert pos.shape == (S,), f"expected {S} positions, got {pos.shape}"
            i += 1
            # Read haplotypes
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
    p = rep.haplotypes.mean(axis=0)  # per-site allele freq
    h = 2 * p * (1 - p) * n / (n - 1)  # corrected per-site heterozygosity
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
    counts = rep.haplotypes.sum(axis=0)  # derived count per site
    folded = np.minimum(counts, n - counts)
    return np.bincount(folded, minlength=n // 2 + 1)[1:n // 2 + 1].astype(float)


def collect_stats(reps: Iterable[MsHaplotypes], n: int) -> dict:
    """Compute per-replicate stats; return arrays."""
    ss_list, pi_list, td_list, sfs_list = [], [], [], []
    for r in reps:
        ss_list.append(segsites(r))
        pi_list.append(pi(r))
        td_list.append(tajima_D(r))
        sfs_list.append(sfs(r))
    sfs_arr = np.array(sfs_list) if sfs_list else np.zeros((0, n // 2))
    return dict(
        ss=np.array(ss_list, dtype=float),
        pi=np.array(pi_list),
        td=np.array(td_list),
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
    generations for diploid Ne.
    """
    return t_cli * 4 * Ne


def discoal_alpha_to_msp_growth(alpha: float, Ne: int) -> float:
    """Convert discoal alpha (4N-scaled) to msprime per-generation growth rate."""
    return alpha / (4 * Ne)


def discoal_theta_to_msp_mu(theta: float, Ne: int, L: int) -> float:
    """Convert discoal theta (4*Ne*mu*L) to msprime per-generation per-site mutation rate."""
    return theta / (4 * Ne * L)


def discoal_rho_to_msp_r(rho: float, Ne: int, L: int) -> float:
    """Convert discoal rho (4*Ne*r*L) to msprime per-generation per-site recombination rate."""
    return rho / (4 * Ne * L)


def _msp_constant_parity_pi(n: int, L: int, theta: float, rho: float,
                              Ne: int, reps: int, conv: Convention) -> float:
    """Mean pi from msprime no-demography simulation under the given convention."""
    mu = discoal_theta_to_msp_mu(theta, Ne, L)
    r = discoal_rho_to_msp_r(rho, Ne, L)
    pop_size = Ne * conv.population_size_factor
    pis = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples=n,
            sequence_length=L,
            recombination_rate=r,
            population_size=pop_size,
            ploidy=conv.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
        pis.append(pi(msprime_to_ms(ts)))
    return float(np.mean(pis))


def check_constant_parity(n: int = 6, theta: float = 5.0, rho: float = 5.0,
                          L: int = 1000, Ne: int = DEFAULT_NE,
                          reps: int = 200, tolerance: float = 0.10,
                          ) -> Convention:
    """Trial all three conventions against a no-demography baseline.

    Returns the Convention whose mean pi agrees with discoal within tolerance.
    Raises RuntimeError if no convention matches.

    Use the returned Convention to seed ACTIVE_CONVENTION:
        ACTIVE_CONVENTION = check_constant_parity()
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
```

- [ ] **Step 3: Discover and lock in the right convention**

```bash
make discoal
python3 test/parity/phase4b_msprime/parity_utils.py
```

Expected output structure:

```
discoal mean pi = X.XXXX (n=6, theta=5.0, rho=5.0, reps=200)
  A: ploidy=1, popsize=Ne       msprime mean pi = X.XXXX, rel_err = 0.XXX  <-- best so far
  B: ploidy=2, popsize=Ne       msprime mean pi = X.XXXX, rel_err = 0.XXX
  C: ploidy=1, popsize=2*Ne     msprime mean pi = X.XXXX, rel_err = 0.XXX

Best convention: <name> (rel_err = 0.XXX)

To lock in the best convention, edit parity_utils.py:
  ACTIVE_CONVENTION = Convention(...)
```

Pick the BEST convention (lowest rel_err) **if it is < 10%**. Edit `parity_utils.py`:

Replace:
```python
ACTIVE_CONVENTION: Optional[Convention] = None
```

with one of:
```python
ACTIVE_CONVENTION = CONVENTION_A
ACTIVE_CONVENTION = CONVENTION_B
ACTIVE_CONVENTION = CONVENTION_C
```

depending on which won.

Then run validation:
```bash
cd test/parity/phase4b_msprime
python3 -c "import parity_utils; parity_utils.validate_active_convention()"
cd -
```

Expected: `validate_active_convention PASS: rel_err = 0.XXX (convention = ...)`.

**Halt conditions** (report BLOCKED):
- If NONE of the three conventions match within 10% (best rel_err > 10%), don't commit. Possible deeper issues:
  - Recombination model semantics differ between simulators.
  - Mutation model differences (infinite-sites vs JC — discoal uses infinite-sites; verify msprime's `sim_mutations` default matches).
  - Sample-count or population-namespacing mismatch.
  Print all three rel_err values and your analysis.
- If two conventions match closely (e.g., A=0.04 and C=0.03), pick the lower but flag it as suspicious — that's an unexpected coincidence and may indicate a near-identity in your specific test parameters that wouldn't generalize.

- [ ] **Step 4: Commit**

```bash
git add test/parity/phase4b_msprime/parity_utils.py \
        test/parity/phase4b_msprime/__init__.py \
        test/parity/phase4b_msprime/results/.gitignore
git commit -m "$(cat <<'EOF'
Add Phase 4b msprime parity utility library

parity_utils.py provides:
- ms-format parser for discoal stdout (skips DEBUG lines / cmdline echo)
- msprime tree-sequence to ms-format converter
- Per-replicate stats: segregating sites, pi, Tajima's D, folded SFS
- Parameter conversions (discoal CLI 2N time -> msprime generations,
  alpha 4N-scaled -> per-generation, theta/rho -> per-gen-per-site)
- check_constant_parity() sanity test that runs no-demography discoal
  and matched msprime, comparing mean pi within 10% tolerance

The constant-parity check runs as the script's __main__ entry to
validate the conversion functions before any actual EXP parity test
is run. Halt on conversion mismatch rather than misreport later.
EOF
)"
```

---

## Task 2: Single-pop EXP growth parity test

**Files:**
- Create: `test/parity/phase4b_msprime/test_single_pop_exp.py`

- [ ] **Step 1: Write the test**

```python
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

    out = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples={"pop0": n},
            sequence_length=L,
            recombination_rate=r,
            demography=demography,
            ploidy=ACTIVE_CONVENTION.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
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
```

- [ ] **Step 2: Smoke-run with small reps**

```bash
python3 test/parity/phase4b_msprime/test_single_pop_exp.py 100
```

Expected output: harness runs both simulators (each takes a few seconds for 100 reps), prints comparison table, and prints PASS or FAIL.

If FAIL with reps=100: this is a tight test and a small N may have spurious significance. Re-run with `reps=1000` to confirm.

If PASS with reps=100: proceed to a real run with 1000+ reps.

If the smoke run errors (e.g., msprime API mismatch, missing fields): debug the harness before attempting larger runs.

- [ ] **Step 3: Commit**

```bash
git add test/parity/phase4b_msprime/test_single_pop_exp.py
git commit -m "Add single-pop EXP growth msprime parity test

Runs discoal -eg 0.5 0 50 and msprime with matched
add_population_parameters_change at the converted generation time
and per-generation growth rate. Computes segregating sites, pi,
Tajima's D, and folded SFS across replicates, then runs Bonferroni-
corrected Kolmogorov-Smirnov / chi-squared tests at alpha = 0.01.

Usage: python3 test_single_pop_exp.py [reps]
Default reps=1000; use 100 for a quick smoke run."
```

---

## Task 3: Two-pop split + EXP growth parity test

**Files:**
- Create: `test/parity/phase4b_msprime/test_two_pop_split_exp.py`

- [ ] **Step 1: Write the test**

```python
"""test_two_pop_split_exp.py — two-population split + EXP growth parity.

Configuration: pop0 + pop1, split (msprime POPULATION_SPLIT) at time
t_split = 1.0 (CLI). Pop0 experiences exponential growth from t=0 backward to
t_growth_end_cli = 0.5, then constant before that. No migration.

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
            samples={"pop0": n_per_pop, "pop1": n_per_pop},
            sequence_length=L,
            recombination_rate=r,
            demography=demography,
            ploidy=ACTIVE_CONVENTION.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
        out.append(msprime_to_ms(ts))
    return out


def main(reps: int = 1000):
    print("=== Phase 4b: two-pop split + EXP parity ===")
    print(f"reps per simulator = {reps}")
    validate_active_convention(reps=200)

    n_per_pop = 5
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
```

- [ ] **Step 2: Smoke-run with small reps**

```bash
python3 test/parity/phase4b_msprime/test_two_pop_split_exp.py 100
```

Expected output: ms-format outputs match between simulators, comparison table, PASS/FAIL.

- [ ] **Step 3: Commit**

```bash
git add test/parity/phase4b_msprime/test_two_pop_split_exp.py
git commit -m "Add two-pop split + EXP growth msprime parity test

Two populations with EXP growth in pop0 and a split at t_split.
Compares discoal (-p 2 n n -eg t pop alpha -ed t_split 0 1) and
msprime (Demography with pop_split + add_population_parameters_change)
on segregating sites, pi, Tajima's D via Bonferroni-corrected KS."
```

---

## Task 4: Run + analyze full parity suite

**Files:**
- Create: `test/parity/phase4b_msprime/run_parity.sh`

- [ ] **Step 1: Write the orchestrator**

```bash
#!/usr/bin/env bash
# Phase 4b: run msprime parity tests at full replicate count, log results.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/results"
mkdir -p "$OUT"

REPS="${REPS:-1000}"

cd "$ROOT"
make discoal >/dev/null 2>&1 || { echo "discoal build failed"; exit 1; }

echo "=== Phase 4b msprime parity tests (REPS=$REPS) ===" | tee "$OUT/analysis.txt"

set +e
python3 "$HERE/test_single_pop_exp.py" "$REPS" 2>&1 | tee -a "$OUT/analysis.txt"
single_rc=$?
echo | tee -a "$OUT/analysis.txt"
python3 "$HERE/test_two_pop_split_exp.py" "$REPS" 2>&1 | tee -a "$OUT/analysis.txt"
two_pop_rc=$?
set -e

echo | tee -a "$OUT/analysis.txt"
if [[ $single_rc -eq 0 && $two_pop_rc -eq 0 ]]; then
  echo "OVERALL: PASS" | tee -a "$OUT/analysis.txt"
  exit 0
else
  echo "OVERALL: FAIL (single=$single_rc, two_pop=$two_pop_rc)" | tee -a "$OUT/analysis.txt"
  exit 1
fi
```

Make executable: `chmod +x test/parity/phase4b_msprime/run_parity.sh`.

- [ ] **Step 2: Run with reps=1000**

```bash
./test/parity/phase4b_msprime/run_parity.sh
```

Expected: takes ~5-10 minutes. Outputs per-test PASS/FAIL and OVERALL PASS/FAIL.

- [ ] **Step 3: Save analysis.txt to commit**

`results/analysis.txt` is generated by the orchestrator. The `.gitignore` from Task 1 excludes it (it has `*.txt` exclusion? Actually no — only `.ms`, `.err`, `.npz`, `.csv`. `.txt` is committed.)

Verify:
```bash
cat test/parity/phase4b_msprime/results/.gitignore
```

If `*.txt` is not in the gitignore, the analysis.txt commits cleanly. If it IS, remove that line — we want analysis.txt committed.

- [ ] **Step 4: Commit results**

```bash
git add test/parity/phase4b_msprime/run_parity.sh \
        test/parity/phase4b_msprime/results/analysis.txt
git commit -m "$(cat <<'EOF'
Add Phase 4b orchestrator and run results

run_parity.sh runs both single-pop and two-pop EXP parity tests
at the configured REPS (default 1000), logs to results/analysis.txt,
and exits with overall PASS/FAIL.

The committed analysis.txt records the outcome at the SHA where
Phase 4b was run.
EOF
)"
```

---

## Task 5: Spec addendum recording the outcome

**Files:**
- Modify: `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md`

- [ ] **Step 1: Append addendum**

At the end of the spec file, append a new section. Use the actual outcome from Task 4's analysis.txt.

If PASS:

```markdown

## Addendum (2026-MM-DD): Phase 4b msprime Parity Result

The Phase 4 SHAPE_EXPONENTIAL implementation was validated against msprime
1.4.1 via the parity harness at `test/parity/phase4b_msprime/`. Both
single-pop (`-eg 0.5 0 50`) and two-pop split + EXP (`-eg 0.3 0 50 -ed 1.0 0 1`)
configurations were run at $10^3$ replicates each, with summary statistics
compared via Bonferroni-corrected Kolmogorov-Smirnov tests at $\alpha = 0.01$.

**Result: PASS.** All comparisons across single-pop {ss, pi, tajD, sfs} and
two-pop {ss, pi, tajD} statistics yielded $p > \alpha_{\text{Bonf}}$, indicating
the discoal SHAPE_EXPONENTIAL implementation is statistically indistinguishable
from msprime's `add_population_parameters_change` under the tested conditions.

Detailed per-comparison output: `test/parity/phase4b_msprime/results/analysis.txt`.
```

If FAIL: substitute the appropriate FAIL paragraph identifying which statistics rejected and at what $D$ / $p$.

- [ ] **Step 2: Commit**

```bash
git add docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md
git commit -m "Spec addendum: Phase 4b msprime parity result"
```

---

## Task 6: Final regression sweep + tag

- [ ] **Step 1: Confirm Phase 3 / Phase 4 regressions still pass**

```bash
make discoal discoal_pre_phase3
./test/parity/phase3_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
make run_tests
```

Expected: all PASS.

- [ ] **Step 2: Tag**

If Phase 4b PASSed:

```bash
git tag -a phase4b-msprime-parity-complete -m "$(cat <<'TAGEOF'
Phase 4b of issue-82 design complete

msprime parity for SHAPE_EXPONENTIAL validated. Single-pop and
two-pop split + EXP configurations both pass Bonferroni-corrected
KS tests at alpha = 0.01 across summary statistics (ss, pi, tajD, sfs).

The discoal -eg implementation is statistically indistinguishable
from msprime's add_population_parameters_change under the tested
conditions.

Detailed result: docs/superpowers/specs/...-design.md addendum and
test/parity/phase4b_msprime/results/analysis.txt.
TAGEOF
)"
```

If FAIL: tag as `phase4b-msprime-parity-fail` with diagnostic notes in the message.

- [ ] **Step 3: List tags**

```bash
git tag --list 'phase*'
```

Expected: 5 tags including `phase4b-msprime-parity-complete` (or `-fail`).

---

## Self-Review Checklist

- [ ] Every task has explicit file paths.
- [ ] Every code step contains the actual code an engineer needs.
- [ ] All Python and bash code is included verbatim.
- [ ] Conversion math is documented in `parity_utils.py` with a runtime sanity check.
- [ ] Bit-equality regression (Phase 3) and `-eg` smoke test (Phase 4) still PASS at the end.
- [ ] No emojis, no Claude/AI references.
- [ ] All work lands on `feature/issue-82-time-varying-demography`.

## What's Next After This Plan

- **Plan: Phase 5** — Sweep accessor wiring with closed-form `detSweepFreqGeneral` per the spec revision. Removes `--det-sweep-mode` flag and `detSweepFreqEuler` function as side effect.
- **Plan: Phase 6** — Add `SHAPE_LINEAR` support, msprime parity (uses the harness from this plan).
- **Plan: Phase 7** — `'em'` migration shape events, importer rewrite, removal of back-derivation hack. **Closes issue #82.**
- **Plan: Phase 8** — Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep.
