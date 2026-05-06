#!/bin/bash
#
# Preflight validation regression suite.
#
# Drives the discoal binary with a battery of CLI invocations that
# exercise validateShapeTrajectories from src/core/shapes.c. Each case
# is one line in run_cases below: a short name, expected exit-status
# class (PASS or FAIL), an expected-stderr substring (FAIL only) and the
# CLI args. PASS cases must exit 0 and produce non-empty output. FAIL
# cases must exit non-zero AND emit the documented error substring.
#
# This is the end-to-end counterpart to the unit tests in
# test/unit/test_shapes.c -- the unit tests cover the validator in
# isolation; this suite confirms the CLI rejects bad inputs at the
# right time and accepts well-formed ones.

set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
DISCOAL=${DISCOAL:-$REPO_ROOT/build/discoal}

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

PASS=0
FAIL=0
FAIL_NAMES=()

run_pass_case() {
    local name=$1; shift
    local out rc
    out=$("$DISCOAL" "$@" 2>&1)
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo -e "${RED}FAIL${NC} $name (expected exit 0, got $rc)"
        echo "  last stderr line: $(printf '%s\n' "$out" | tail -1)"
        FAIL_NAMES+=("$name"); FAIL=$((FAIL + 1)); return
    fi
    if ! printf '%s\n' "$out" | grep -q '^//$'; then
        echo -e "${RED}FAIL${NC} $name (exit 0 but no '//' segment header in output)"
        FAIL_NAMES+=("$name"); FAIL=$((FAIL + 1)); return
    fi
    echo -e "${GREEN}PASS${NC} $name"
    PASS=$((PASS + 1))
}

run_fail_case() {
    local name=$1 expect=$2; shift 2
    local out rc
    out=$("$DISCOAL" "$@" 2>&1)
    rc=$?
    if [ "$rc" -eq 0 ]; then
        echo -e "${RED}FAIL${NC} $name (expected non-zero exit, got 0)"
        FAIL_NAMES+=("$name"); FAIL=$((FAIL + 1)); return
    fi
    if ! printf '%s\n' "$out" | grep -qF -- "$expect"; then
        echo -e "${RED}FAIL${NC} $name (exit $rc but expected substring missing)"
        echo "  expected: $expect"
        echo "  last stderr line: $(printf '%s\n' "$out" | tail -1)"
        FAIL_NAMES+=("$name"); FAIL=$((FAIL + 1)); return
    fi
    echo -e "${GREEN}PASS${NC} $name"
    PASS=$((PASS + 1))
}

echo -e "${BLUE}=== Preflight validation regression suite ===${NC}"
echo "discoal: $DISCOAL"
echo

if [ ! -x "$DISCOAL" ]; then
    echo -e "${RED}discoal binary not found at $DISCOAL. Run 'make discoal' first.${NC}" >&2
    exit 2
fi

# ----- Cases that should pass -----
echo -e "${YELLOW}-- Success cases (must exit 0 and produce '//' output) --${NC}"

run_pass_case "constant_single_pop" \
    5 1 100 -t 5

run_pass_case "constant_two_pop_with_migration" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -m 1 0 0.5

run_pass_case "en_positive_size" \
    5 1 100 -t 5 -en 0.1 0 0.5

run_pass_case "eg_grows_backward" \
    5 1 100 -t 5 -eg 0.1 0 0.5

run_pass_case "el_short_interval_no_zero_crossing" \
    5 1 100 -t 5 -el 0.1 0 0.5 -en 0.2 0 1.0

run_pass_case "em_zero_turns_migration_off" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -em 0.1 0 1 0.0

run_pass_case "em_after_pop_merge_skipped" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -ed 0.05 1 0 -em 0.1 0 1 0.5

# ----- Cases that should fail with a documented error -----
echo
echo -e "${YELLOW}-- Failure cases (must exit non-zero with expected message) --${NC}"

run_fail_case "en_zero_size" \
    "size must be strictly positive" \
    5 1 100 -t 5 -en 0.1 0 0.0

run_fail_case "en_negative_size" \
    "size must be strictly positive" \
    5 1 100 -t 5 -en 0.1 0 -0.5

run_fail_case "el_positive_gamma_zero_crossing_before_next_event" \
    "drives population 0 size to zero" \
    5 1 100 -t 5 -el 0.05 0 5.0 -en 1.0 0 1.0

run_fail_case "m_initial_negative" \
    "initial migration rate" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 -0.5

run_fail_case "M_initial_negative" \
    "initial migration rate" \
    5 1 100 -t 5 -p 2 5 0 -M -0.5

run_fail_case "em_negative_rate" \
    "rate must be non-negative" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -em 0.1 0 1 -0.2

echo
echo -e "${BLUE}=== Summary ===${NC}"
echo "Pass: $PASS"
echo "Fail: $FAIL"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed cases: ${FAIL_NAMES[*]}"
    exit 1
fi
