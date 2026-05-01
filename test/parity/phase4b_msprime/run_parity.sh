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
