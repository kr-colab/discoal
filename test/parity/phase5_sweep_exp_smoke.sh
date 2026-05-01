#!/usr/bin/env bash
# Phase 5 smoke: confirm sweep + -eg runs to completion.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"
[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

CONFIGS=(
  "det_sweep_with_eg|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -eg 0.5 0 50 -d 12345 67890"
  "stoch_sweep_with_eg|6 5 1000 -t 5 -r 5 -ws 0.05 -a 200 -x 0.5 -eg 0.5 0 50 -d 12345 67890"
)

for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  if "$DISCOAL" $args > /tmp/$tag.ms 2>/tmp/$tag.err; then
    n_segsites=$(grep -c '^segsites:' /tmp/$tag.ms || echo 0)
    if [[ $n_segsites -ge 1 ]]; then
      echo "  PASS ($n_segsites segsites entries)"
    else
      echo "  FAIL: no segsites lines"
      exit 1
    fi
  else
    echo "  FAIL: discoal exited nonzero"
    cat /tmp/$tag.err >&2
    exit 1
  fi
done

echo
echo "PASS: all sweep+eg smoke configurations ran to completion."
