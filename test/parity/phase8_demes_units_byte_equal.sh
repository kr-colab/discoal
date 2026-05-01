#!/usr/bin/env bash
# Phase 8 regression: a demes file with time_units=years and one with
# time_units=generations describing the same demographic model must
# produce byte-equal output. Demes is the same simulation; only the
# time-unit convention differs.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

OUT="$HERE/phase8_demes_units/results"
mkdir -p "$OUT"

cd "$ROOT"
"$DISCOAL" -Y "$HERE/phase8_demes_units/years.yaml" > "$OUT/years.ms" 2>"$OUT/years.err" || {
    echo "FAIL: discoal exited non-zero on years fixture"
    cat "$OUT/years.err"
    exit 1
}
"$DISCOAL" -Y "$HERE/phase8_demes_units/gens.yaml" > "$OUT/gens.ms" 2>"$OUT/gens.err" || {
    echo "FAIL: discoal exited non-zero on gens fixture"
    cat "$OUT/gens.err"
    exit 1
}

# Diff excluding line 1 (the command-line echo, which differs).
if diff <(sed '1d' "$OUT/years.ms") <(sed '1d' "$OUT/gens.ms") > "$OUT/diff.out"; then
    echo "PASS: years-based and generations-based demes fixtures produce byte-equal output."
    exit 0
else
    echo "FAIL: output differs between time_units=years and time_units=generations"
    head -20 "$OUT/diff.out"
    exit 1
fi
