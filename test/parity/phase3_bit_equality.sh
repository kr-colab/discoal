#!/usr/bin/env bash
# Phase 3 regression: verify build/discoal output is byte-identical to
# build/discoal_pre_phase3 across a configuration grid. Excludes the
# line-1 command-line echo (which differs because the binary paths differ).

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

PRE="$ROOT/build/discoal_pre_phase3"
CUR="$ROOT/build/discoal"
[[ -x "$PRE" ]] || { echo "FAIL: build/discoal_pre_phase3 missing (make discoal_pre_phase3)"; exit 1; }
[[ -x "$CUR" ]] || { echo "FAIL: build/discoal missing (make discoal)"; exit 1; }

# Configuration grid. Each line is: tag <discoal args>
CONFIGS=(
  "single_pop_neutral|6 5 1000 -t 5 -r 5 -d 12345 67890"
  "single_pop_size_change|6 5 1000 -t 5 -r 5 -en 0.5 0 0.1 -en 1.0 0 1.0 -d 12345 67890"
  "two_pop_constant_mig|8 5 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.3 -d 12345 67890"
  "two_pop_split|8 5 1000 -t 5 -r 5 -p 2 4 4 -ed 0.5 0 1 -d 12345 67890"
  "two_pop_split_mig|8 5 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.3 -ed 1.0 0 1 -d 12345 67890"
  "three_pop_combo|9 5 1000 -t 5 -r 5 -p 3 3 3 3 -m 0 1 0.5 -ed 0.5 0 1 -ed 1.5 1 2 -d 12345 67890"
  "single_pop_with_admixture|8 5 1000 -t 5 -r 5 -p 2 4 4 -ea 0.2 0 0 1 0.7 -ed 1.0 0 1 -d 12345 67890"
  "ancient_sample|6 5 1000 -t 5 -r 5 -p 2 4 0 -A 2 1 0.3 -ed 0.5 0 1 -d 12345 67890"
)

OUT="$HERE/phase3_bit_equality_results"
mkdir -p "$OUT"

fails=0
for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # Run both; record exit codes. Don't abort on non-zero so identical
  # crashes (e.g. a pre-existing bug exercised by the same input) still
  # get diffed against each other.
  # shellcheck disable=SC2086
  "$PRE" $args > "$OUT/${tag}_pre.ms" 2>"$OUT/${tag}_pre.err" && pre_rc=$? || pre_rc=$?
  # shellcheck disable=SC2086
  "$CUR" $args > "$OUT/${tag}_cur.ms" 2>"$OUT/${tag}_cur.err" && cur_rc=$? || cur_rc=$?
  # diff excluding line 1 (the command-line echo, which differs by binary path)
  if [[ "$pre_rc" != "$cur_rc" ]]; then
    echo "  FAIL  (exit codes differ: pre=$pre_rc cur=$cur_rc)"
    fails=$((fails+1))
    continue
  fi
  if diff <(sed '1d' "$OUT/${tag}_pre.ms") <(sed '1d' "$OUT/${tag}_cur.ms") > "$OUT/${tag}.diff"; then
    echo "  PASS"
  else
    echo "  FAIL  (see $OUT/${tag}.diff)"
    head -20 "$OUT/${tag}.diff" || true
    fails=$((fails+1))
  fi
done

if [[ $fails -eq 0 ]]; then
  echo
  echo "PASS: all ${#CONFIGS[@]} configurations are byte-equal."
  exit 0
else
  echo
  echo "FAIL: $fails of ${#CONFIGS[@]} configurations differ."
  exit 1
fi
