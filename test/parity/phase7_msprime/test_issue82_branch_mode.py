"""test_issue82_branch_mode.py - branch-mode parity diagnostic for issue #82 fixture.

Compares tskit branch-mode statistics computed on tree sequences from discoal
and msprime on the issue #82 demes fixture. Branch-mode stats bypass the
mutation process: they're computed directly from tree topology and branch lengths.

If branch-mode stats MATCH, the residual site-mode divergence (ss, Watterson
theta, Tajima's D) is a mutation-placement problem (e.g., a theta-to-mu
convention mismatch).

If branch-mode stats DIFFER, the trees themselves differ (a coalescent
simulator divergence).

Tajima's D in branch mode is scale-invariant and used as the primary diagnostic.
For diversity and segregating-sites-branch (which differ by an overall time-scale
factor), each replicate's value is divided by the across-replicate mean to
produce a unitless distribution before KS comparison.
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
import tskit
from scipy import stats

# Reuse Phase 4b parity utilities for msprime setup.
sys.path.insert(0, str(Path(__file__).parent.parent / "phase4b_msprime"))
from parity_utils import Convention

PHASE7_CONVENTION = Convention(
    name="C: ploidy=1, popsize=2*Ne (Phase 7 branch-mode)",
    ploidy=1, population_size_factor=2.0,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
DISCOAL = REPO_ROOT / "build" / "discoal"
WRAPPER_YAML = REPO_ROOT / "config_examples" / "demes_example.yaml"
DEMES_FILE = REPO_ROOT / "config_examples" / "demes_example.demes.yaml"


def make_temp_yaml(reps: int) -> Path:
    base = WRAPPER_YAML.read_text()
    overridden = base.replace("num_replicates: 10", f"num_replicates: {reps}")
    if f"num_replicates: {reps}" not in overridden:
        raise RuntimeError("failed to override num_replicates in YAML")
    tmp = Path(tempfile.mkstemp(suffix=".yaml", prefix="phase7_branch_")[1])
    tmp.write_text(overridden)
    return tmp


def run_discoal_trees(reps: int, prefix: Path) -> list[Path]:
    """Run discoal and return list of tree-sequence file paths.

    Discoal writes <prefix_basename_rep{N}.{ext}> for each replicate (1-indexed).
    """
    yaml_path = make_temp_yaml(reps)
    cmd = [str(DISCOAL), "-Y", str(yaml_path), "-ts", str(prefix)]
    # Run from REPO_ROOT so the relative demes_filename in the wrapper YAML
    # ("config_examples/demes_example.demes.yaml") resolves correctly.
    proc = subprocess.run(
        cmd, capture_output=True, text=True, check=False, cwd=str(REPO_ROOT),
    )
    yaml_path.unlink()
    if proc.returncode != 0:
        raise RuntimeError(f"discoal exited {proc.returncode}: {proc.stderr[-500:]}")
    # discoal writes <stem>_rep{N}<ext>; e.g., for prefix='/tmp/d.trees',
    # outputs are /tmp/d_rep1.trees, /tmp/d_rep2.trees, ...
    stem = prefix.stem
    ext = prefix.suffix
    out_dir = prefix.parent
    files = sorted(
        out_dir.glob(f"{stem}_rep*{ext}"),
        key=lambda p: int(p.stem.rsplit("rep", 1)[1])
    )
    if len(files) != reps:
        raise RuntimeError(f"expected {reps} tree files, found {len(files)}")
    return files


def run_msprime_trees(reps: int, n_per_pop: dict, L: int) -> list[tskit.TreeSequence]:
    graph = demes.load(str(DEMES_FILE))
    demography = msprime.Demography.from_demes(graph)
    for pop in demography.populations:
        pop.initial_size = pop.initial_size * PHASE7_CONVENTION.population_size_factor
    ne_discoal = 10000
    L_int = int(L)
    rho = 15.0
    r = rho / (4.0 * ne_discoal * L_int)
    samples = [
        msprime.SampleSet(num_samples=n, population=pop_name, ploidy=1)
        for pop_name, n in n_per_pop.items()
    ]
    out = []
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples=samples, sequence_length=L_int,
            recombination_rate=r, demography=demography,
            ploidy=PHASE7_CONVENTION.ploidy, random_seed=seed,
        )
        out.append(ts)
    return out


def compute_branch_stats(ts: tskit.TreeSequence) -> dict:
    """Compute branch-mode stats for a single tree sequence.

    Returns dict with keys: diversity, segregating_sites, tajimas_d.
    All computed in mode="branch" (operates on tree topology + branch lengths,
    independent of mutations).
    """
    samples = ts.samples()
    return {
        "diversity": ts.diversity(sample_sets=[samples], mode="branch")[0],
        "segregating_sites": ts.segregating_sites(sample_sets=[samples], mode="branch")[0],
        "tajimas_d": ts.Tajimas_D(sample_sets=[samples], mode="branch")[0],
    }


def collect_branch_stats(tree_sources, source_label: str) -> dict:
    """tree_sources is a list of either Path (load with tskit.load) or TreeSequence."""
    rows = []
    for src in tree_sources:
        ts = tskit.load(str(src)) if isinstance(src, Path) else src
        rows.append(compute_branch_stats(ts))
    arr = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    return arr


def main(reps: int = 500):
    print(f"=== Phase 7 branch-mode parity diagnostic (REPS={reps}) ===")
    print(f"Convention: {PHASE7_CONVENTION.name}")
    print()

    # Discoal: simulate to tree-sequence files.
    tmpdir = Path(tempfile.mkdtemp(prefix="phase7_branch_"))
    try:
        prefix = tmpdir / "discoal.trees"
        print(f"discoal -> {prefix} ({reps} reps)")
        t0 = time.time()
        discoal_files = run_discoal_trees(reps, prefix)
        print(f"  discoal: {time.time()-t0:.1f}s, wrote {len(discoal_files)} tree files")

        # msprime: simulate tree sequences in memory.
        n_per_pop = {"A": 10, "B": 5, "C": 5}
        L = 50000
        print(f"msprime Demography.from_demes ({reps} reps)")
        t0 = time.time()
        msp_ts_list = run_msprime_trees(reps, n_per_pop, L)
        print(f"  msprime: {time.time()-t0:.1f}s")

        # Compute branch-mode stats.
        s_discoal = collect_branch_stats(discoal_files, "discoal")
        s_msp = collect_branch_stats(msp_ts_list, "msprime")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    print()
    print("Raw branch-mode stats (note units differ by scale factor):")
    print(f"discoal mean: diversity={s_discoal['diversity'].mean():.4e}, "
          f"seg_sites={s_discoal['segregating_sites'].mean():.4e}, "
          f"tajD={s_discoal['tajimas_d'].mean():.4f}")
    print(f"msprime mean: diversity={s_msp['diversity'].mean():.4e}, "
          f"seg_sites={s_msp['segregating_sites'].mean():.4e}, "
          f"tajD={s_msp['tajimas_d'].mean():.4f}")

    # Scale ratio (informational only).
    if s_msp['diversity'].mean() > 0 and s_discoal['diversity'].mean() > 0:
        scale_ratio = s_msp['diversity'].mean() / s_discoal['diversity'].mean()
        print(f"  scale ratio msprime/discoal (diversity): {scale_ratio:.4f} "
              f"(expected ~20000 if discoal is in coal units, msprime in gens)")

    # Mean-normalize before comparison so the time-scale mismatch doesn't
    # cause spurious rejection.
    print()
    print("KS comparisons on mean-normalized distributions:")
    print(f"  Tajima's D is scale-invariant; compared directly.")
    print()

    comparisons = []
    for stat_name in ("diversity", "segregating_sites", "tajimas_d"):
        if stat_name == "tajimas_d":
            d_arr = s_discoal[stat_name]
            m_arr = s_msp[stat_name]
        else:
            d_arr = s_discoal[stat_name] / s_discoal[stat_name].mean()
            m_arr = s_msp[stat_name] / s_msp[stat_name].mean()
        D, p = stats.ks_2samp(d_arr, m_arr)
        comparisons.append((f"branch_{stat_name}", D, p))

    n_comp = len(comparisons)
    bonf_alpha = 0.01 / n_comp
    print(f"{n_comp} comparisons, Bonferroni alpha = {bonf_alpha:.2e}")
    print(f"{'comparison':<28} {'D':>10} {'p':>12}  reject?")
    rejected = 0
    for name, D, p in comparisons:
        flag = "  **" if p < bonf_alpha else ""
        if p < bonf_alpha:
            rejected += 1
        print(f"{name:<28} {D:>10.4f} {p:>12.2e}{flag}")
    print()
    if rejected == 0:
        print(f"PASS: branch-mode parity holds across all {n_comp} statistics.")
        print()
        print("Diagnostic conclusion: the trees themselves match in distribution.")
        print("The residual site-mode divergence (ss, wtheta, td in test_issue82_fixture.py)")
        print("must therefore arise in MUTATION PLACEMENT, not in the coalescent simulation.")
        print("Investigate: theta-to-mu conversion, mutation-placement model differences.")
        sys.exit(0)
    else:
        print(f"FAIL: {rejected}/{n_comp} branch-mode comparisons reject equality.")
        print()
        print("Diagnostic conclusion: the TREES differ in distribution.")
        print("This points to a coalescent-simulation divergence between discoal")
        print("and msprime running the same demes graph.")
        sys.exit(1)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    main(reps=reps)
