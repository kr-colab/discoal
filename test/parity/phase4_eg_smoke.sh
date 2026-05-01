#!/usr/bin/env bash
# Phase 4 smoke test: confirm -eg produces distributionally different output
# from comparable -en configurations. This is a sanity check that the EXP
# shape is actually being used by the simulation engine.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

# Three pairs of configs; each pair (constant vs eg) should yield different output
# at the same RNG seeds.

NREPS=20
configs=(
  "no_demography|6 $NREPS 1000 -t 5 -r 5 -d 12345 67890"
  "with_eg_alpha50|6 $NREPS 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890"
  "with_en_size_change|6 $NREPS 1000 -t 5 -r 5 -en 0.5 0 0.5 -d 12345 67890"
)

OUT="$HERE/phase4_eg_smoke_results"
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

# Confirm with_eg_alpha50 differs from no_demography (sed 1d strips the cmdline echo)
if diff <(sed '1d' "$OUT/no_demography.ms") <(sed '1d' "$OUT/with_eg_alpha50.ms") > /dev/null; then
  echo "FAIL: -eg config produced byte-identical output to no-demography config (eg appears to be a no-op)"
  exit 1
fi

# Confirm with_eg_alpha50 differs from with_en_size_change
if diff <(sed '1d' "$OUT/with_eg_alpha50.ms") <(sed '1d' "$OUT/with_en_size_change.ms") > /dev/null; then
  echo "FAIL: -eg and -en configs produced byte-identical output (eg may be silently treated as en)"
  exit 1
fi

echo
echo "PASS: -eg produces distributionally different output from constant-N and from -en configurations."
