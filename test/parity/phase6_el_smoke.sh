#!/usr/bin/env bash
# Phase 6 smoke: confirm -el produces distributionally different output
# from comparable -en and -eg configurations. Sanity check that the
# SHAPE_LINEAR path is exercised.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

NREPS=20
configs=(
  "no_demography|6 $NREPS 1000 -t 5 -r 5 -d 12345 67890"
  "with_el_gamma_pos|6 $NREPS 1000 -t 5 -r 5 -el 0.5 0 0.5 -d 12345 67890"
  "with_en_size_change|6 $NREPS 1000 -t 5 -r 5 -en 0.5 0 0.5 -d 12345 67890"
  "with_eg_alpha|6 $NREPS 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890"
)

OUT="$HERE/phase6_el_smoke_results"
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

if diff <(sed '1d' "$OUT/no_demography.ms") <(sed '1d' "$OUT/with_el_gamma_pos.ms") > /dev/null; then
  echo "FAIL: -el config produced byte-identical output to no-demography (el appears to be a no-op)"
  exit 1
fi

if diff <(sed '1d' "$OUT/with_el_gamma_pos.ms") <(sed '1d' "$OUT/with_en_size_change.ms") > /dev/null; then
  echo "FAIL: -el and -en configs produced byte-identical output (el may be silently treated as en)"
  exit 1
fi

if diff <(sed '1d' "$OUT/with_el_gamma_pos.ms") <(sed '1d' "$OUT/with_eg_alpha.ms") > /dev/null; then
  echo "FAIL: -el and -eg configs produced byte-identical output (el may be silently treated as eg)"
  exit 1
fi

echo
echo "PASS: -el produces distributionally different output from no-demography, -en, and -eg."
