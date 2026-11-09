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

# The -ed brings back lineages that moved to population 1 before
# migration stops; without it they could never coalesce, and discoal
# would stop with an error.
run_pass_case "em_zero_turns_migration_off" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -em 0.1 0 1 0.0 -ed 0.2 1 0

run_pass_case "em_after_pop_merge_skipped" \
    5 1 100 -t 5 -p 2 5 0 -m 0 1 0.5 -ed 0.05 1 0 -em 0.1 0 1 0.5

# Coalescence that takes longer than the old time cap of 100000 units (2N):
# discoal used to stop there and write the replicate unfinished.
run_pass_case "coalescence_after_old_time_cap" \
    2 3 100 -t 0.00001 -en 0.01 0 1e6

# An event later than the old time cap; the stranded-lineage check must not
# treat the interval before it as the last one.
run_pass_case "event_after_old_time_cap" \
    4 2 100 -t 0.00001 -p 2 2 2 -ed 60000 1 0

# After migration stops, one population holds only a lineage without
# ancestral material; that must not count as stranded. (Other seeds of this
# model do strand lineages and stop with an error.)
run_pass_case "lineage_without_ancestral_material_is_not_stranded" \
    2 1 2 -t 1 -r 20 -p 2 2 0 -m 0 1 2 -m 1 0 2 -em 1.0 0 1 0 -em 1.0 1 0 0 -d 66 73

# Migration that starts later, with -eM and no -m, is still migration for
# the check that a model with several populations can coalesce.
run_pass_case "eM_alone_counts_as_migration" \
    4 1 100 -t 2 -p 2 2 2 -eM 0.1 1 -d 1 2

# Size changes with priors (#96): the events hold the lower bounds until
# the draw, so the check on sizes passes.
run_pass_case "Pe1_and_Pe2_run" \
    20 5 10000 -t 10 -Pe1 0.1 0.2 0.1 1.0 -Pe2 0.3 0.4 0.5 2.0

# Recurrent sweeps that happen: the sampler used to read a sweep
# trajectory that was never mapped, and crashed.
run_pass_case "R_recurrent_sweeps_happen" \
    20 3 10000 -t 10 -r 10 -R 5 -a 500 -d 12345 10102

run_pass_case "L_recurrent_sweeps_left_of_locus_happen" \
    20 3 10000 -t 10 -r 10 -L 5 -a 500 -d 12345 10102

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

run_fail_case "Pe1_zero_size_bound" \
    "0 < size low <= size high" \
    20 5 10000 -t 10 -Pe1 0.1 0.2 0 1.0

run_fail_case "Px_with_sweep_outside_locus" \
    "cannot be combined with a sweep outside the locus" \
    20 2 1000 -t 5 -ls 0.01 10 -a 100 -Px 0.4 0.6

run_fail_case "Pu_with_recurrent_sweeps" \
    "needs a single sweep, not recurrent sweeps" \
    20 2 1000 -t 5 -R 0.1 -a 100 -Pu 0 0.1

run_fail_case "A_negative_time" \
    "the time must be >= 0" \
    4 2 100 -t 2 -A 1 0 -0.1

echo
echo -e "${BLUE}=== Summary ===${NC}"
echo "Pass: $PASS"
echo "Fail: $FAIL"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed cases: ${FAIL_NAMES[*]}"
    exit 1
fi
