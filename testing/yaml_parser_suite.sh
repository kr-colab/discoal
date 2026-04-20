#!/bin/bash
#
# YAML parser regression suite.
#
# Walks every fixture in config_examples/debug_yaml/ exercising it
# against `./build/discoal -Y`:
#
#   error_*.yaml  must exit non-zero AND emit the error message that the
#                 fixture header documents (the line of the form
#                 `#   "Error parsing config: ..."`)
#   other *.yaml  must exit 0
#
# This is the end-to-end counterpart to test_config_interface.c, which
# exercises the parser via its API rather than the discoal binary.

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
DISCOAL=${DISCOAL:-$REPO_ROOT/build/discoal}
FIXTURES="$REPO_ROOT/config_examples/debug_yaml"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

PASS=0
FAIL=0
FAIL_NAMES=()

# Pull the documented expected error message out of a fixture header.
# Every error_*.yaml uses the convention `#   "<message>"` on its own
# line. Trailing " ..." (continuation) and "<placeholder>" suffixes are
# stripped so the result is a literal substring suitable for `grep -F`.
extract_expected() {
    grep -m1 '^#   "' "$1" \
        | sed -e 's/^#   "//' \
              -e 's/"$//' \
              -e 's/ \.\.\.$//' \
              -e 's/<.*//'
}

run_error_case() {
    local f=$1
    local name
    name=$(basename "$f")
    local expected
    expected=$(extract_expected "$f")
    local out
    out=$("$DISCOAL" -Y "$f" 2>&1 >/dev/null)
    local rc=$?
    if [ "$rc" -eq 0 ]; then
        echo -e "${RED}FAIL${NC} $name (expected non-zero exit, got 0)"
        FAIL_NAMES+=("$name")
        FAIL=$((FAIL + 1))
        return
    fi
    if [ -n "$expected" ] && ! printf '%s\n' "$out" | grep -qF -- "$expected"; then
        echo -e "${RED}FAIL${NC} $name"
        echo "  expected substring: $expected"
        echo "  got: $(printf '%s\n' "$out" | head -1)"
        FAIL_NAMES+=("$name")
        FAIL=$((FAIL + 1))
        return
    fi
    echo -e "${GREEN}PASS${NC} $name"
    PASS=$((PASS + 1))
}

run_success_case() {
    local f=$1
    local name
    name=$(basename "$f")
    local out
    out=$("$DISCOAL" -Y "$f" 2>&1)
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo -e "${RED}FAIL${NC} $name (expected exit 0, got $rc)"
        echo "  last stderr line: $(printf '%s\n' "$out" | tail -1)"
        FAIL_NAMES+=("$name")
        FAIL=$((FAIL + 1))
        return
    fi
    echo -e "${GREEN}PASS${NC} $name"
    PASS=$((PASS + 1))
}

echo -e "${BLUE}=== YAML parser regression suite ===${NC}"
echo "discoal:  $DISCOAL"
echo "fixtures: $FIXTURES"
echo

if [ ! -x "$DISCOAL" ]; then
    echo -e "${RED}discoal binary not found at $DISCOAL. Run 'make discoal' first.${NC}" >&2
    exit 2
fi
if [ ! -d "$FIXTURES" ]; then
    echo -e "${RED}fixtures directory $FIXTURES not found.${NC}" >&2
    exit 2
fi

echo -e "${YELLOW}-- Error cases (must exit != 0 with expected message) --${NC}"
for f in "$FIXTURES"/error_*.yaml; do
    [ -e "$f" ] || continue
    run_error_case "$f"
done
echo
echo -e "${YELLOW}-- Success cases (must exit 0) --${NC}"
for f in "$FIXTURES"/*.yaml; do
    [ -e "$f" ] || continue
    case $(basename "$f") in
        error_*) ;;
        *) run_success_case "$f" ;;
    esac
done
echo
echo -e "${BLUE}=== Summary ===${NC}"
echo "Pass: $PASS"
echo "Fail: $FAIL"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed fixtures: ${FAIL_NAMES[*]}"
    exit 1
fi
