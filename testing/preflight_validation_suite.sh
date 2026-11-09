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
# A case that runs longer than this fails, so a hang cannot stall CI.
CASE_TIMEOUT=${CASE_TIMEOUT:-60}

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
    out=$(timeout "$CASE_TIMEOUT" "$DISCOAL" "$@" 2>&1)
    rc=$?
    if [ "$rc" -ne 0 ]; then
        [ "$rc" -eq 124 ] && rc="124, timed out after ${CASE_TIMEOUT} s"
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
    out=$(timeout "$CASE_TIMEOUT" "$DISCOAL" "$@" 2>&1)
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

# A partial sweep must end above one copy of the allele, 1/(2N) with
# N = 10^6 by default.
run_fail_case "c_below_one_copy" \
    "must be above 1/(2N)" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -c 0.0000001

run_fail_case "Pc_high_below_one_copy" \
    "must be above 1/(2N)" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pc 0.00000001 0.0000001

# A prior with a low bound below one copy is allowed: a draw below one
# copy gives a sweep that ends at once.
run_pass_case "Pc_low_bound_zero" \
    20 20 1000 -t 5 -r 5 -ws 0 -a 500 -Pc 0 0.000001 -d 1 2

# Above one copy at N0, but below it at the sweep, after a size change to
# 0.01: the trajectory has no steps, and the sweep ends at once.
run_pass_case "c_below_one_copy_after_size_change" \
    20 5 10000 -t 1 -r 10 -ws 0.01 -a 500 -c 0.00001 -en 0 0 0.01 -d 1 2

run_pass_case "R_c_below_one_copy_after_size_change" \
    20 20 10000 -t 1 -r 10 -R 2 -a 500 -c 0.00001 -en 0 0 0.01 -d 1 2

# Sweeps outside the locus with two populations: left-recombination
# events must be chosen among the sweep events only, or a draw can take a
# lineage from an empty beneficial class.
run_pass_case "ls_two_populations" \
    20 10 10000 -t 1 -r 10 -p 2 10 10 -ed 5 1 0 -ls 0.005 100 -a 500 -N 10000 -d 1 2

run_pass_case "L_two_populations" \
    20 10 10000 -t 10 -p 2 10 10 -m 0 1 2 -ed 0.3 1 0 -L 1 -a 500 -d 20 1020

# A size change in population 1 must not set the sweep trajectory's N:
# here it would end the trajectory at once while population 0 is still
# sweeping.
run_pass_case "en_other_population_during_sweep" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -ws 0.01 -a 500 -c 0.3 -en 0.01 1 0.00015 -d 1 2

run_pass_case "en_other_population_recurrent_sweeps" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -R 50 -a 500 -c 0.3 -en 0.001 1 0.00015 -d 1 2

# With growth (-eg) the trajectory's end follows N, so it can end before
# the sampler's fixed end frequency; the sweep then ends with the
# trajectory.
run_pass_case "eg_during_sweep" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 200 -d 1 2

run_pass_case "eg_during_partial_sweep" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -c 0.001 -eg 0 0 200 -d 1 2

run_pass_case "eg_during_sweep_default_N" \
    20 20 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 2 -d 3 4

# Exponential growth over a long time: after the sweep one lineage is left
# in each population, and while they wait to migrate the size underflows to
# 0. Coalescence must then still happen (at once).
run_pass_case "eG_size_underflow_two_populations" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -ws 0.01 -a 500 -eG 0 200 -d 3 6

echo
echo -e "${BLUE}=== Summary ===${NC}"
echo "Pass: $PASS"
echo "Fail: $FAIL"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed cases: ${FAIL_NAMES[*]}"
    exit 1
fi
