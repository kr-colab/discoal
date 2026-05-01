#!/usr/bin/env bash
# Phase 5 regression: verify build/discoal sweep output is byte-identical to
# build/discoal_pre_phase5 for SHAPE_CONSTANT-only configurations.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

PRE="$ROOT/build/discoal_pre_phase5"
CUR="$ROOT/build/discoal"
[[ -x "$PRE" ]] || { echo "FAIL: build/discoal_pre_phase5 missing (make discoal_pre_phase5)"; exit 1; }
[[ -x "$CUR" ]] || { echo "FAIL: build/discoal missing (make discoal)"; exit 1; }

# Sweep configurations covering deterministic / stochastic-forward / neutral-stochastic
# modes, varying alpha and tau, single-pop and two-pop.
CONFIGS=(
  "single_pop_neutral_sweep|6 5 1000 -t 5 -r 5 -wn 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_stoch_sweep|6 5 1000 -t 5 -r 5 -ws 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep_alpha50|6 5 1000 -t 5 -r 5 -wd 0.5 -a 50 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep_size_change|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -en 0.1 0 0.5 -d 12345 67890"
  "two_pop_det_sweep|8 5 1000 -t 5 -r 5 -p 2 4 4 -wd 0.05 -a 200 -x 0.5 -ed 0.5 0 1 -d 12345 67890"
)

OUT="$HERE/phase5_sweep_bit_equality_results"
mkdir -p "$OUT"

fails=0
for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # shellcheck disable=SC2086
  "$PRE" $args > "$OUT/${tag}_pre.ms" 2>"$OUT/${tag}_pre.err"; pre_rc=$?
  # shellcheck disable=SC2086
  "$CUR" $args > "$OUT/${tag}_cur.ms" 2>"$OUT/${tag}_cur.err"; cur_rc=$?
  if [[ $pre_rc -ne $cur_rc ]]; then
    echo "  FAIL exit code mismatch: pre=$pre_rc cur=$cur_rc"
    fails=$((fails+1))
    continue
  fi
  if diff <(sed '1d' "$OUT/${tag}_pre.ms") <(sed '1d' "$OUT/${tag}_cur.ms") > "$OUT/${tag}.diff"; then
    echo "  PASS"
  else
    echo "  FAIL  (see $OUT/${tag}.diff)"
    head -10 "$OUT/${tag}.diff" || true
    fails=$((fails+1))
  fi
done

if [[ $fails -eq 0 ]]; then
  echo
  echo "PASS: all ${#CONFIGS[@]} sweep configurations are byte-equal."
  exit 0
else
  echo
  echo "FAIL: $fails of ${#CONFIGS[@]} sweep configurations differ."
  exit 1
fi
