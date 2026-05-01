#!/usr/bin/env python3
"""YAML parser fuzzer for discoal.

Generates many edge-case variants of a small set of valid baseline configs
and runs each through `./build/discoal -Y`. Flags any outcome that isn't a
clean exit 0 (parsed and ran) or exit 1 (parser rejected with a message).

Usage:
    python3 testing/yaml_parser_fuzz.py [--seed N] [--cases N] [--timeout S]

This is exploratory, not deterministic in scope. Not run by CI. Findings
should be triaged by hand and, when reproducible, promoted to a
`testing/debug_yaml/bug_*.yaml` reproducer.
"""
import argparse
import copy
import json
import math
import random
import shutil
import string
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DISCOAL = REPO_ROOT / "build" / "discoal"
FINDINGS_DIR = REPO_ROOT / "testing" / "fuzz_findings"

# Mirrors `#define MAXPOPS 121` in src/core/discoal.h
MAXPOPS = 121

# Exit codes we consider "clean": parser handled the input gracefully one
# way or the other. Anything else (signals, asserts, allocator aborts) is
# a real defect.
CLEAN_EXITS = {0, 1, 2}

# Numeric edge values to drop into double-typed fields.
DOUBLE_EDGES = [
    float("nan"),
    float("inf"),
    float("-inf"),
    0.0,
    -0.0,
    1.0,
    -1.0,
    1e-300,
    -1e-300,
    1e300,
    -1e300,
    2**31 - 1,
    -(2**31),
    2**53,  # largest exact integer in double
    -(2**53),
]

# Edge values for int-typed fields (cyaml will reject most non-int yaml
# scalars at parse time, but these probe the parser's int range checks).
INT_EDGES = [
    0,
    -1,
    1,
    2**31 - 1,
    -(2**31),
    2**31,
    -(2**31) - 1,
]

# Some fields directly drive simulator runtime (loops scaled by their
# value); a huge but well-formed value parses instantly and then runs
# forever. We still want to probe those fields for sign/zero/NaN edges
# but not for "huge but valid". Keep the small/nonsense edges only.
RUNTIME_FIELDS = {
    ("simulation", "num_replicates"),
    ("simulation", "num_sites"),
    ("genetics", "recombination_rate"),
    ("genetics", "gene_conversion_rate"),
}
SAFE_DOUBLE_EDGES = [
    float("nan"), float("inf"), float("-inf"),
    0.0, -0.0, -1.0, 1e-300, -1e-300, -1e300, -(2**31), -(2**53),
]
SAFE_INT_EDGES = [0, -1, -(2**31), -(2**31) - 1]

# Pop-index edge values, parameterized by num_demes the baseline declares.
def pop_index_edges(num_demes):
    return [-1, 0, num_demes - 1, num_demes, num_demes + 1, MAXPOPS, MAXPOPS + 1, -(2**31), 2**31 - 1]


# Fields that are sized populations and where we have a strong "must be
# rejected" prior. These map a (block, key) to a predicate the value must
# satisfy. Anything failing the predicate that nonetheless gets exit 0
# is flagged ACCEPT_SUSPECT.
EXPECT_REJECT = {
    ("genetics", "mutation_rate"):       lambda v: not (math.isnan(v) or math.isinf(v) or v < 0),
    ("genetics", "recombination_rate"):  lambda v: not (math.isnan(v) or math.isinf(v) or v < 0),
    ("genetics", "gene_conversion_rate"):lambda v: not (math.isnan(v) or math.isinf(v) or v < 0),
    ("simulation", "sample_size"):       lambda v: isinstance(v, int) and 0 < v < 2**31,
    ("simulation", "num_replicates"):    lambda v: isinstance(v, int) and 0 < v < 2**31,
    ("simulation", "num_sites"):         lambda v: isinstance(v, int) and 0 < v < 2**31,
    # Optional field; cyaml zero-fills, so 0 is the "not set" sentinel and
    # legitimately accepted. Anything else must be a positive int in range.
    ("demography", "effective_population_size"):
        lambda v: v == 0 or (0 < v <= 2**31 - 1 and float(v).is_integer()),
}


BASELINES = {
    "single_pop": {
        "simulation": {"sample_size": 10, "num_replicates": 1, "num_sites": 1000, "seed": [1, 2]},
        "genetics":   {"mutation_rate": 0.1, "recombination_rate": 0.1},
        "demography": {"deme_sample_size": [10], "effective_population_size": 10000},
    },
    "two_pop_migration": {
        "simulation": {"sample_size": 20, "num_replicates": 1, "num_sites": 1000, "seed": [1, 2]},
        "genetics":   {"mutation_rate": 0.1, "recombination_rate": 0.1},
        "demography": {
            "deme_sample_size": [10, 10],
            "effective_population_size": 10000,
            "migration_matrix": [
                {"row": [0.0, 0.1]},
                {"row": [0.1, 0.0]},
            ],
        },
    },
    "single_pop_sweep": {
        "simulation": {"sample_size": 10, "num_replicates": 1, "num_sites": 1000, "seed": [1, 2]},
        "genetics":   {"mutation_rate": 0.1, "recombination_rate": 0.1},
        "demography": {"deme_sample_size": [10], "effective_population_size": 10000},
        "selection":  {
            "sweep_mode": "stochastic",
            "selection_coefficient": 100,
            "sweep_position": 0.5,
            "fixation_time_ago": 0.01,
        },
    },
}


# --- mutators ----------------------------------------------------------------

def random_string(rng, n):
    return "".join(rng.choices(string.printable, k=n))


def numeric_field_mutations(rng, baseline_name, baseline):
    """Yield (label, mutated, expectation) for one-field numeric replacements."""
    for block_name, block in baseline.items():
        if not isinstance(block, dict):
            continue
        for key, val in list(block.items()):
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                runtime = (block_name, key) in RUNTIME_FIELDS
                if isinstance(val, float):
                    edges = SAFE_DOUBLE_EDGES if runtime else DOUBLE_EDGES
                else:
                    edges = (SAFE_INT_EDGES + SAFE_DOUBLE_EDGES) if runtime \
                        else (INT_EDGES + DOUBLE_EDGES)
                for edge in edges:
                    mutated = copy.deepcopy(baseline)
                    mutated[block_name][key] = edge
                    pred = EXPECT_REJECT.get((block_name, key))
                    expect_reject = pred is not None and not pred(edge)
                    yield (
                        f"{baseline_name}/{block_name}.{key}={repr(edge)}",
                        mutated,
                        expect_reject,
                    )


def pop_index_mutations(rng, baseline_name, baseline):
    """Inject demographic_events with edge pop indices."""
    if "demography" not in baseline:
        return
    n_demes = len(baseline["demography"]["deme_sample_size"])
    for pop_edge in pop_index_edges(n_demes):
        mutated = copy.deepcopy(baseline)
        mutated["demography"]["demographic_events"] = {
            "population_size_changes": [
                {"time": 0.1, "size": 0.5, "population": pop_edge}
            ],
        }
        # Anything outside [0, n_demes) ought to be rejected.
        expect_reject = pop_edge < 0 or pop_edge >= n_demes
        yield (f"{baseline_name}/size_change.population={pop_edge}", mutated, expect_reject)

    if n_demes >= 2:
        for derived in pop_index_edges(n_demes):
            for ancestral in pop_index_edges(n_demes):
                mutated = copy.deepcopy(baseline)
                mutated["demography"]["demographic_events"] = {
                    "population_splits": [
                        {"time": 0.1, "derived": derived, "ancestral": ancestral}
                    ],
                }
                expect_reject = (
                    derived < 0 or derived >= n_demes
                    or ancestral < 0 or ancestral >= n_demes
                    or derived == ancestral
                )
                yield (
                    f"{baseline_name}/split.d={derived},a={ancestral}",
                    mutated,
                    expect_reject,
                )


def deme_sample_size_mutations(rng, baseline_name, baseline):
    """Perturb the deme_sample_size array to wrong length, negatives, etc."""
    if "demography" not in baseline:
        return
    sample_size = baseline["simulation"]["sample_size"]
    sizes_to_try = [
        [],                            # empty
        [sample_size],                 # right total but wrong length if multi-pop
        [-1, sample_size + 1],         # negative entry, sum still matches
        [2**31 - 1, -(2**31 - 1) + sample_size],  # overflow but sum matches
        [sample_size // 2] * (MAXPOPS + 1),       # too many
        [sample_size, sample_size],    # sum > sample_size
        [0, 0],                        # sum < sample_size
    ]
    for sizes in sizes_to_try:
        mutated = copy.deepcopy(baseline)
        mutated["demography"]["deme_sample_size"] = sizes
        # No ground-truth for whether each gets rejected, only a no-crash
        # contract. Set expect_reject False so we don't false-flag accepts.
        yield (f"{baseline_name}/deme_sample_size={sizes}", mutated, False)


def migration_matrix_mutations(rng, baseline_name, baseline):
    """Stress the migration_matrix shape and values."""
    if "demography" not in baseline or "migration_matrix" not in baseline["demography"]:
        return
    n = len(baseline["demography"]["migration_matrix"])
    options = [
        [{"row": []}] * n,                        # empty rows
        [{"row": [0.0] * (n + 1)}] * n,           # too many cols
        [{"row": [0.0] * (n - 1)}] * n if n > 1 else None,  # too few cols
        [{"row": [float("nan")] * n}] * n,        # NaN
        [{"row": [-1.0] * n}] * n,                # all-negative
        [{"row": [0.5] * n}] * n,                 # non-zero diagonal
        [{"row": [0.0] * n}] * (n + 1),           # extra row
    ]
    for opt in options:
        if opt is None:
            continue
        mutated = copy.deepcopy(baseline)
        mutated["demography"]["migration_matrix"] = opt
        yield (f"{baseline_name}/migmat-shape", mutated, False)


def string_field_mutations(rng, baseline_name, baseline):
    """Stress tree_sequence_filename: empty, very long, control chars, traversal."""
    long_str = "A" * 4096
    weird_str = "\x00\x01\x02embedded-control"
    traversal = "../" * 200 + "etc/passwd"
    for s in ["", long_str, weird_str, traversal, "\n", "\t", " "]:
        mutated = copy.deepcopy(baseline)
        mutated.setdefault("output", {})
        mutated["output"]["output_type"] = "tree_sequence"
        mutated["output"]["tree_sequence_filename"] = s
        # Empty string we know should be rejected (#67); others should
        # at least not crash.
        yield (f"{baseline_name}/tree_seq_filename={s!r}", mutated, s == "")


def structural_mutations(rng, baseline_name, baseline):
    """Drop required fields and inject unknown fields."""
    # Drop one required key at a time
    required = [
        ("simulation", "sample_size"),
        ("simulation", "num_replicates"),
        ("simulation", "num_sites"),
        ("genetics", "mutation_rate"),
        ("genetics", "recombination_rate"),
        ("demography", "deme_sample_size"),
    ]
    for block, key in required:
        if block in baseline and key in baseline[block]:
            mutated = copy.deepcopy(baseline)
            del mutated[block][key]
            yield (f"{baseline_name}/drop.{block}.{key}", mutated, False)

    # Inject one unknown field at a time
    for block in ["simulation", "genetics", "demography"]:
        if block in baseline:
            mutated = copy.deepcopy(baseline)
            mutated[block]["totally_made_up_key"] = "value"
            yield (f"{baseline_name}/unknown.{block}", mutated, False)

    # Inject unknown top-level block
    mutated = copy.deepcopy(baseline)
    mutated["this_block_does_not_exist"] = {"foo": 1}
    yield (f"{baseline_name}/unknown_top_block", mutated, False)


def random_combo_mutations(rng, baseline_name, baseline, n=200):
    """Apply N independent random numeric mutations per call."""
    mutators = list(numeric_field_mutations(rng, baseline_name, baseline))
    if not mutators:
        return
    for _ in range(n):
        # Pick a random mutation as the basis, then chain a second random
        # field perturbation on top.
        a = rng.choice(mutators)
        _, base, _ = a
        b = rng.choice(mutators)
        _, mut2, _ = b
        # Trivial merge: take base, copy b's deltas. b is also a fully
        # mutated config so just swap one block.
        block = rng.choice(list(base.keys()))
        if block in mut2:
            merged = copy.deepcopy(base)
            merged[block] = mut2[block]
            yield (f"combo/{a[0]}+{b[0]}", merged, False)


# --- runner ------------------------------------------------------------------

import yaml as pyyaml


def yaml_dump(cfg):
    # default_flow_style=False keeps it readable for triage. allow_unicode
    # and explicit start mark intentional.
    return pyyaml.safe_dump(cfg, default_flow_style=False, allow_unicode=True)


def run_one(content, timeout_s):
    """Write the YAML to stdin via a temp file is awkward; use process pipe."""
    p = subprocess.Popen(
        [str(DISCOAL), "-Y", "/dev/stdin"],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    try:
        _, err = p.communicate(content.encode("utf-8"), timeout=timeout_s)
        return p.returncode, err.decode("utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        p.kill()
        p.communicate()
        return None, "<timeout>"


def classify(rc, expect_reject):
    if rc is None:
        return "HANG"
    if rc not in CLEAN_EXITS:
        return "CRASH"
    if expect_reject and rc == 0:
        return "ACCEPT_SUSPECT"
    return "OK_REJECT" if rc != 0 else "OK_RUN"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cases", type=int, default=0,
                    help="cap total cases; 0 = exhaustive enumeration of all generators")
    ap.add_argument("--timeout", type=float, default=5.0)
    ap.add_argument("--keep-findings", action="store_true",
                    help="don't wipe testing/fuzz_findings before running")
    args = ap.parse_args()

    if not DISCOAL.exists():
        print(f"discoal binary not found at {DISCOAL}; run `make discoal` first.",
              file=sys.stderr)
        sys.exit(2)

    rng = random.Random(args.seed)

    if FINDINGS_DIR.exists() and not args.keep_findings:
        shutil.rmtree(FINDINGS_DIR)
    FINDINGS_DIR.mkdir(parents=True, exist_ok=True)

    cases = []
    for name, baseline in BASELINES.items():
        cases.extend(numeric_field_mutations(rng, name, baseline))
        cases.extend(pop_index_mutations(rng, name, baseline))
        cases.extend(deme_sample_size_mutations(rng, name, baseline))
        cases.extend(migration_matrix_mutations(rng, name, baseline))
        cases.extend(string_field_mutations(rng, name, baseline))
        cases.extend(structural_mutations(rng, name, baseline))
        cases.extend(random_combo_mutations(rng, name, baseline, n=300))

    if args.cases:
        rng.shuffle(cases)
        cases = cases[: args.cases]

    counts = {"OK_RUN": 0, "OK_REJECT": 0, "ACCEPT_SUSPECT": 0, "CRASH": 0, "HANG": 0}
    findings = []

    print(f"Fuzzing {len(cases)} cases against {DISCOAL.name} (timeout={args.timeout}s)")
    for i, (label, cfg, expect_reject) in enumerate(cases):
        try:
            content = yaml_dump(cfg)
        except Exception as e:
            # Some mutated configs are invalid Python (e.g. NaN can't dump
            # to YAML safely depending on dialect). Skip but count.
            counts.setdefault("DUMP_ERROR", 0)
            counts["DUMP_ERROR"] += 1
            continue
        rc, err = run_one(content, args.timeout)
        kind = classify(rc, expect_reject)
        counts[kind] += 1
        if kind in ("CRASH", "HANG", "ACCEPT_SUSPECT"):
            findings.append({
                "label": label,
                "kind": kind,
                "exit_code": rc,
                "stderr_tail": err[-500:],
                "expect_reject": expect_reject,
            })
            # Save reproducer
            safe_label = label.replace("/", "__").replace(" ", "_")[:120]
            (FINDINGS_DIR / f"{kind}_{i:05d}_{safe_label}.yaml").write_text(content)
        if (i + 1) % 200 == 0:
            print(f"  ... {i+1}/{len(cases)} done; counts so far: {counts}")

    print()
    print("=== Fuzz summary ===")
    for k, v in counts.items():
        print(f"  {k:<15} {v}")
    print()
    if findings:
        print(f"Findings written to {FINDINGS_DIR.relative_to(REPO_ROOT)}/")
        manifest = FINDINGS_DIR / "manifest.json"
        manifest.write_text(json.dumps(findings, indent=2, default=str))
        # Brief preview
        for f in findings[:10]:
            print(f"  [{f['kind']}] {f['label']} exit={f['exit_code']}")
        if len(findings) > 10:
            print(f"  ... and {len(findings) - 10} more in {manifest.name}")
        sys.exit(1)
    print("No crashes, hangs, or unexpected accepts.")
    sys.exit(0)


if __name__ == "__main__":
    main()
