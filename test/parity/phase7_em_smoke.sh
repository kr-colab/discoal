#!/usr/bin/env bash
# Phase 7 smoke: confirm -em produces distributionally different output
# from comparable constant-migration configs.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

NREPS=20
configs=(
  "no_migration|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -ed 1.0 0 1 -d 12345 67890"
  "constant_migration|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.5 -ed 1.0 0 1 -d 12345 67890"
  "with_em_off|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.5 -em 0.3 0 1 0.0 -em 0.3 1 0 0.0 -ed 1.0 0 1 -d 12345 67890"
)

OUT="$HERE/phase7_em_smoke_results"
mkdir -p "$OUT"

for cfg in "${configs[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # shellcheck disable=SC2086
  "$DISCOAL" $args > "$OUT/${tag}.ms" 2>"$OUT/${tag}.err" || {
    echo "  FAIL: discoal exited non-zero (see $OUT/${tag}.err)"
    exit 1
  }
done

if diff <(sed '1d' "$OUT/no_migration.ms") <(sed '1d' "$OUT/constant_migration.ms") > /dev/null; then
  echo "FAIL: no-migration and constant-migration configs are byte-identical (migration appears not to be exercised)"
  exit 1
fi

if diff <(sed '1d' "$OUT/constant_migration.ms") <(sed '1d' "$OUT/with_em_off.ms") > /dev/null; then
  echo "FAIL: constant-migration and -em-off configs are byte-identical (-em appears to be a no-op)"
  exit 1
fi

echo
echo "PASS: -em produces distributionally different output from constant-migration baseline."
