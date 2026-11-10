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

# The sweep is outside the locus only with -ls, -ld, -ln or -L, not
# because of the value of -x: with -Px, -x is a placeholder.
run_pass_case "x_placeholder_with_Px" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -x -1 -Px 0 1 -d 1 2

run_fail_case "x_minus_one_without_sweep_outside_locus" \
    "-x (-1) must be in [0, 1]" \
    20 2 1000 -t 5 -r 5 -R 2 -a 500 -x -1

run_fail_case "Px_with_L" \
    "cannot be combined with a sweep outside the locus" \
    20 2 1000 -t 5 -r 5 -L 1 -a 500 -Px 0 1

run_fail_case "x_inside_locus_with_L" \
    "cannot be combined with a sweep outside the locus" \
    20 2 1000 -t 5 -r 5 -L 1 -a 500 -x 0.3

# A YAML sweep_position with -L from the command line.
run_fail_case "yaml_position_with_cli_L" \
    "YAML sweep_position) cannot be combined with a sweep outside the locus" \
    -Y "$REPO_ROOT/testing/debug_yaml/recurrent_sweeps_with_position.yaml" -L 1

run_fail_case "x_negative_with_ls" \
    "cannot be combined with a sweep outside the locus" \
    20 2 1000 -t 5 -r 5 -ls 0.01 10 -a 500 -x -0.5

# Recurrent sweeps check, as the neutral sampler does, that the lineages
# can still meet after the last event (here migration from population 2
# into population 1 strands lineages after population 1 merges).
run_fail_case "R_lineages_can_never_meet" \
    "they can never coalesce" \
    20 2 10000 -t 10 -r 10 -R 0.5 -a 500 -p 3 7 7 6 -M 1 -ed 2 1 0 -ed 3 2 0

# One single sweep at most: the sweep sampler follows one sweep at a time.
run_fail_case "two_ws" \
    "only one single sweep can be given" \
    20 5 1000 -t 5 -r 5 -ws 0.01 -ws 0.012 -a 500

run_fail_case "ws_and_ls" \
    "only one single sweep can be given" \
    20 5 1000 -t 5 -r 5 -ws 0.01 -ls 0.05 5 -a 500

run_fail_case "yaml_sweep_and_cli_ws" \
    "only one single sweep can be given" \
    -Y "$REPO_ROOT/testing/debug_yaml/sweep_with_tau.yaml" -ws 0.05

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

# Sweeps are not modeled under growth: population 0 may not grow (-eg,
# -el, -eG, -eL) at or after the earliest start of a sweep, and no
# population may grow with recurrent sweeps.
run_fail_case "eg_during_sweep" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 200

run_fail_case "eg_during_partial_sweep" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -c 0.001 -eg 0 0 200

run_fail_case "eg_starting_after_sweep_start" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0.05 0 2

run_fail_case "el_during_sweep" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -el 0 0 2

run_fail_case "eG_during_sweep" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -ws 0.01 -a 500 -eG 0 2

run_fail_case "Pu_reaching_into_growth" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 5 -en 0.005 0 1 -Pu 0.001 0.02

run_fail_case "Pe1_may_end_growth_after_sweep_start" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 5 -Pe1 0.004 0.03 0.5 2

run_fail_case "R_with_growth_in_population_1" \
    "not modeled under growth" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -R 2 -a 500 -eg 0 1 2

# A size change with a prior on its time ends growth only if it surely
# comes after the growth starts and surely before the sweep.
run_fail_case "Pe1_may_come_before_growth" \
    "not modeled under growth" \
    10 3 1000 -N 10000 -t 5 -r 5 -a 200 -ws 0.05 -eg 0.01 0 30 -Pe1 0.001 0.0101 0.5 1

run_pass_case "Pe1_ends_growth_before_sweep" \
    10 3 1000 -N 10000 -t 5 -r 5 -a 200 -ws 0.05 -eg 0.01 0 30 -Pe1 0.02 0.04 0.5 1 -d 1 2

# An -en at the sweep's time ends growth before the sweep only if it is
# given before the sweep (events at one time apply in command-line order).
run_fail_case "en_at_sweep_time_after_sweep" \
    "not modeled under growth" \
    10 2 1000 -N 10000 -t 5 -r 5 -a 200 -ws 0.05 -eg 0.01 0 30 -en 0.05 0 1

run_pass_case "en_at_sweep_time_before_sweep" \
    10 2 1000 -N 10000 -t 5 -r 5 -a 200 -eg 0.01 0 30 -en 0.05 0 1 -ws 0.05 -d 1 2

# With -Pu, only its low bound counts as the sweep's start.
run_pass_case "Pu_after_growth" \
    10 2 1000 -N 10000 -t 5 -r 5 -a 200 -ws 0.01 -Pu 0.05 0.08 -eg 0.02 0 5 -en 0.03 0 1 -d 1 2

# A trajectory from a file (-S) does not use the acceptance step, and its
# sampler follows the size at each step, so growth is allowed with it.
run_pass_case "S_with_growth" \
    20 2 10000 -t 1 -r 10 -N 10000 -S "$REPO_ROOT/config_examples/trajectory_example.traj" -eg 0 0 1 -d 1 2

run_pass_case "growth_ending_before_sweep" \
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -eg 0 0 5 -en 0.005 0 1 -d 1 2

run_pass_case "growth_in_population_1_with_sweep" \
    20 5 10000 -N 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 1 -m 1 0 1 -ws 0.01 -a 500 -eg 0 1 2 -d 1 2

# Priors on the sweep's frequencies: 0 <= low <= high < 1, and a soft
# sweep that is partial needs every f0 below every final frequency.
run_fail_case "Pc_negative_low" \
    "the bounds must satisfy 0 <= low <= high < 1" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pc -0.5 0.5

run_fail_case "Pc_high_above_one" \
    "the bounds must satisfy 0 <= low <= high < 1" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pc 0.5 1.5

run_fail_case "Pf_low_above_high" \
    "the bounds must satisfy 0 <= low <= high < 1" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pf 0.5 0.4 -c 0.6

run_pass_case "f_with_Pc" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -f 0.05 -Pc 0.3 0.6 -d 1 2

run_fail_case "f_above_Pc_low" \
    "must be below the final frequency" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -f 0.5 -Pc 0.3 0.6

run_pass_case "Pf_with_c" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pf 0.01 0.1 -c 0.5 -d 1 2

run_fail_case "Pf_high_above_c" \
    "must be below the final frequency" \
    20 2 1000 -t 5 -r 5 -ws 0 -a 500 -Pf 0.01 0.6 -c 0.5

# A size change of population 0 long after the sweep: the trajectory has
# ended by then, so the change must not lower the acceptance of every
# proposed trajectory (the run took longer than a minute).
run_pass_case "en_long_after_sweep" \
    20 2 10000 -t 1 -r 10 -ws 0.01 -a 500 -en 1.0 0 0.00001 -d 1 2

# Recurrent soft sweeps are not implemented; -Pf makes a sweep soft too.
run_fail_case "R_with_Pf" \
    "recurrent soft sweeps" \
    20 2 1000 -t 5 -r 5 -R 1 -a 500 -Pf 0.01 0.1

# A sweep needs at least one individual in population 0; every size of
# population 0 (an -en, or the low bound of a size prior) is checked.
run_fail_case "sweep_in_population_below_one" \
    "at least one individual" \
    20 2 1000 -t 5 -r 5 -en 0 0 0.0000001 -ws 0 -a 500

run_fail_case "sweep_with_size_prior_below_one" \
    "at least one individual" \
    10 5 1000 -N 10000 -t 5 -r 5 -a 200 -Pe1 0.01 0.02 0.000001 0.001 -ws 0.05

# Priors (-P flags) and fixed values follow the rules of their YAML
# fields; a negative rate could hang the run (-Pr, -g, -R).
run_fail_case "Pt_negative_low" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -Pt -5 5

run_fail_case "Pt_low_above_high" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -Pt 10 5

run_fail_case "Pr_negative_low" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -Pr -5 5

run_fail_case "Pre_negative_mean" \
    "the mean and the upper bound must be > 0" \
    20 2 1000 -t 5 -r 5 -Pre -1 5

run_fail_case "Pa_negative_low" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -Pa -100 100

run_fail_case "Pu_negative_low" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -Pu -0.1 0.1

run_fail_case "PuA_negative_low" \
    "the bounds must satisfy 0 <= low <= high" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -PuA -1 1

run_fail_case "Px_outside_locus" \
    "the bounds must satisfy 0 <= low <= high <= 1" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -Px -0.5 1.5

run_fail_case "t_negative" \
    "-t (-5) must be >= 0" \
    20 2 1000 -t -5 -r 5

run_fail_case "r_negative" \
    "-r (-5) must be >= 0" \
    20 2 1000 -t 5 -r -5

run_fail_case "g_negative_rate" \
    "-g (rate) (-1) must be >= 0" \
    20 2 1000 -t 5 -r 5 -g -1 100

run_fail_case "g_zero_tract_length" \
    "tract length (YAML gene_conversion_tract_length) (0) must be >= 1" \
    20 2 1000 -t 5 -r 5 -g 10 0

run_fail_case "gr_zero_tract_length" \
    "tract length (YAML gene_conversion_tract_length) (0) must be >= 1" \
    20 2 1000 -t 5 -r 5 -gr 0.5 0

run_pass_case "g_zero_rate_zero_tract_length" \
    20 2 1000 -t 5 -r 5 -g 0 0

run_fail_case "sample_size_above_limit" \
    "the sample size (70000) must be <= 65535" \
    70000 1 1000 -t 5

run_fail_case "num_sites_above_maxsites" \
    "the number of sites (2e+08) must be <= 100000000" \
    20 1 200000000 -t 5

run_fail_case "a_negative" \
    "-a (-100) must be >= 0" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a -100

run_fail_case "x_outside_locus" \
    "-x (1.5) must be in [0, 1]" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -x 1.5

run_fail_case "f_above_one" \
    "-f (1.5) must be in [0, 1)" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -f 1.5

run_fail_case "i_zero" \
    "-i (0) must be > 0" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -i 0

run_fail_case "uA_negative" \
    "-uA (-1) must be >= 0" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -uA -1

run_fail_case "ws_negative_time" \
    "a sweep time" \
    20 2 1000 -t 5 -r 5 -ws -0.1 -a 500

run_fail_case "ls_negative_leftRho" \
    "(leftRho) (-5) must be >= 0" \
    20 2 1000 -t 5 -r 5 -ls 0.01 -5 -a 500

run_fail_case "R_negative_rate" \
    "-R or -L (rate) (-1) must be > 0" \
    20 2 1000 -t 5 -r 5 -R -1 -a 500

run_fail_case "ea_proportion_above_one" \
    "-ea (proportion) (1.5) must be in [0, 1]" \
    20 2 1000 -t 5 -r 5 -p 3 10 5 5 -ea 0.1 0 1 2 1.5 -ed 0.5 1 0 -ed 0.6 2 0

run_fail_case "p_negative_sample_size" \
    "a -p sample size (-5) must be >= 0" \
    20 2 1000 -t 5 -r 5 -p 2 25 -5 -ed 0.5 1 0

run_fail_case "p_sample_sizes_do_not_add_up" \
    "the sum of the -p sample sizes (15) must be equal to the sample size" \
    20 2 1000 -t 5 -r 5 -p 2 10 5 -ed 0.5 1 0

run_fail_case "p_one_population_wrong_size" \
    "the sum of the -p sample sizes (10) must be equal to the sample size" \
    20 1 1000 -t 5 -r 5 -p 1 10

run_fail_case "ws_nan_with_Pu" \
    "a sweep time (-ws, -wd, -wn, -ls, -ld, -ln) (nan) must be >= 0" \
    20 2 1000 -t 5 -r 5 -ws nan -a 500 -Pu 0.05 0.1

run_fail_case "x_nan" \
    "-x (nan) must be in [0, 1]" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -x nan

# A fixed value that a prior replaces is not used, so it is not checked.
run_pass_case "placeholders_replaced_by_priors" \
    20 2 1000 -t -1 -r -1 -ws -1 -a -1 -Pt 1 10 -Pr 0 5 -Pu 0 0.1 -Pa 100 1000 -d 1 2

run_pass_case "priors_inside_their_ranges" \
    20 2 1000 -t 5 -r 5 -ws 0.01 -a 500 -Pt 0 10 -Pr 0 20 -Pa 100 1000 -Pu 0 0.05 -PuA 0 0.01 -Px 0 1 -d 1 2

# Exponential growth over a long time: one lineage is left in each
# population, and while they wait to migrate the size underflows to 0.
# Coalescence must then still happen (at once).
run_pass_case "eG_size_underflow_two_populations" \
    20 5 10000 -t 1 -r 10 -p 2 10 10 -m 0 1 0.1 -m 1 0 0.1 -eG 0 200 -d 3 6

echo
echo -e "${BLUE}=== Summary ===${NC}"
echo "Pass: $PASS"
echo "Fail: $FAIL"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed cases: ${FAIL_NAMES[*]}"
    exit 1
fi
