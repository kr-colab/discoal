#!/usr/bin/env bash
# Q1 verification: detSweepFreq closed-form vs Euler step under constant N.
#
# Runs discoal in two modes across a grid of (alpha, tau) configs,
# generates niceStats summary statistics, and compares distributions
# with KS tests via q1_detsweep_analyze.py.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="$HERE/q1_results"
mkdir -p "$OUT"

DISCOAL="$ROOT/build/discoal"
NICESTATS="$ROOT/build/niceStats"
[[ -x "$DISCOAL" ]] || { echo "build discoal first"; exit 1; }
[[ -x "$NICESTATS" ]] || { echo "build niceStats first"; exit 1; }

REPS=10000
N=10
NSITES=10000
THETA=10
RHO=10

ALPHAS=(50 200 1000)
TAUS=(0.01 0.1 0.5)

for alpha in "${ALPHAS[@]}"; do
  for tau in "${TAUS[@]}"; do
    cfg="alpha${alpha}_tau${tau}"
    for mode in closed euler; do
      out="$OUT/${cfg}_${mode}.ms"
      stats="$OUT/${cfg}_${mode}.stats"
      echo "=== $cfg $mode ==="
      "$DISCOAL" "$N" "$REPS" "$NSITES" \
        -t "$THETA" -r "$RHO" \
        -wd "$tau" -a "$alpha" -x 0.5 \
        -d 12345 67890 \
        --det-sweep-mode "$mode" > "$out"
      "$NICESTATS" "$N" "$REPS" < "$out" > "$stats"
    done
  done
done

echo
echo "All configurations run. Pass result directory to q1_detsweep_analyze.py:"
echo "  python3 $HERE/q1_detsweep_analyze.py $OUT"
