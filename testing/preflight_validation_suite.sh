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

# The last interval has no end, so a linear decline back in time with no
# later event reaches size zero; it must be rejected, not run with a
# negative size.
run_fail_case "el_zero_crossing_in_last_interval" \
    "drives population 1 size to zero" \
    4 2 100 -t 2 -p 2 2 2 -M 1 -el 0.0 1 10

# An admixture moves every lineage out of the admixed population, so its
# size after the admixture does not matter.
run_pass_case "el_in_population_ended_by_admixture" \
    6 5 100 -t 2 -p 3 2 2 2 -el 0 2 1 -ea 0.1 2 0 1 0.5 -ed 1.0 1 0

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

# Demographic events with a population outside 0 .. npops-1 were ignored; a
# negative time stopped inside tskit; a NaN rate made the run hang.
run_fail_case "eg_population_out_of_range" \
    "population 5 must be between 0 and 1" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -eg 0 5 1.0

run_fail_case "en_negative_time" \
    "-en (YAML population_size_changes): the time (-1) must be a finite number >= 0" \
    4 1 100 -t 1 -en -1 0 0.5

run_fail_case "eG_nan_rate" \
    "Invalid numeric argument for -eG: 'nan'" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -eG 0 nan

run_fail_case "em_same_population" \
    "the two populations must be different (population 1 is given twice)" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -em 0 1 1 1

run_fail_case "ed_into_itself" \
    "the two populations must be different (population 0 is given twice)" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -ed 0.5 0 0

# An -ea source can be the admixed population: a pulse.
run_pass_case "ea_pulse_source_is_admixed_population" \
    4 2 100 -t 2 -p 2 2 2 -M 1 -ea 0.5 0 1 0 0.3 -d 1 2

run_fail_case "ed_population_out_of_range" \
    "population 4 must be between 0 and 1" \
    4 1 100 -t 1 -p 2 2 2 -ed 1 0 4

# The command line: a missing or non-numeric value, or an option name that
# does not exist, stops with a message. These crashed, ran with the value
# read as 0, or were ignored.
run_fail_case "en_missing_value" \
    "option -en needs more values" \
    10 1 100 -t 5 -en 0.1 0

run_fail_case "ws_missing_value" \
    "option -ws needs more values" \
    10 1 100 -t 5 -ws

run_fail_case "A_missing_value" \
    "option -A needs more values" \
    10 1 100 -t 5 -A 2 0

run_fail_case "en_non_numeric_time" \
    "Invalid numeric argument for -en: 'abc'" \
    10 1 100 -t 5 -en abc 0 0.5

run_fail_case "en_non_numeric_population" \
    "Invalid integer argument for -en: 'zero'" \
    10 1 100 -t 5 -en 0.1 zero 0.5

run_fail_case "ws_non_numeric" \
    "Invalid numeric argument for -ws: 'abc'" \
    10 1 100 -t 5 -ws abc -a 100

run_fail_case "positional_non_numeric" \
    "Invalid integer argument for numReplicates: 'x'" \
    10 x 100 -t 5

run_fail_case "positional_not_an_integer" \
    "Invalid integer argument for nSites: '1e3'" \
    10 1 1e3 -t 5

run_fail_case "unknown_e_option" \
    "Unknown option '-ez'" \
    10 1 100 -t 5 -ez 1 2 3

run_fail_case "option_prefix_alone" \
    "Unknown option '-P'" \
    10 1 100 -t 5 -P

run_fail_case "unknown_u_option" \
    "Unknown option '-uX'" \
    10 1 100 -t 5 -uX 3

# Only the exact name -Y loads a file: -Yzzz is not a YAML run (it gets
# the usage text, as any call with too few arguments), and after a full
# command line it is an unknown option.
# Values that are not finite, population indices of -m outside the model,
# and sizes below the smallest model stop with a message. (-R inf ran
# without end; -m 10 1000 1 wrote outside the migration matrix; 0
# replicates gave no output and no message.)
run_fail_case "R_infinite" \
    "Invalid numeric argument for -R: 'inf'" \
    4 1 1000 -R inf -a 100

run_fail_case "t_overflow" \
    "Invalid numeric argument for -t: '1e400'" \
    4 1 100 -t 1e400

run_fail_case "m_population_out_of_range" \
    "-m 10 1000: the two populations must be different and between 0 and 1" \
    4 1 100 -t 1 -p 2 2 2 -m 10 1000 1

run_fail_case "m_same_population" \
    "-m 1 1: the two populations must be different" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -m 1 1 2

run_fail_case "zero_replicates" \
    "the number of replicates (0) must be >= 1" \
    10 0 100 -t 5

run_fail_case "zero_sites" \
    "the number of sites (0) must be >= 1" \
    10 1 0 -t 5

run_fail_case "one_sample" \
    "the sample size (1) must be >= 2" \
    1 1 100 -t 5

run_fail_case "p_zero_populations" \
    "-p needs at least 1 population" \
    4 1 100 -t 1 -p 0

run_fail_case "c_out_of_range" \
    "-c (1.5) must be above 0 and below 1" \
    10 1 100 -t 5 -ws 0 -a 100 -c 1.5

run_fail_case "F_before_ts" \
    "-F needs -ts before it" \
    10 1 100 -t 5 -F

run_fail_case "Y_name_with_suffix" \
    "usage: discoal" \
    -Yzzz "$REPO_ROOT/config_examples/basic_example.yaml"

run_fail_case "Y_name_with_suffix_after_arguments" \
    "Unknown option '-Yzzz'" \
    10 1 100 -t 5 -Yzzz "$REPO_ROOT/config_examples/basic_example.yaml"

run_fail_case "Y_twice" \
    "-Y can be given only once" \
    -Y "$REPO_ROOT/config_examples/basic_example.yaml" -Y "$REPO_ROOT/config_examples/basic_example.yaml"

run_fail_case "ts_directory_missing" \
    "cannot write the tree sequence" \
    10 1 100 -t 5 -ts /nonexistent-directory/x.trees

# Long option names. Each command below runs with the short names and with
# the long names, which the help text gives ("  -t, --theta ..."), and the
# two outputs must be identical. So the help text, the table of long names
# and the parser agree. Together the commands use every option.
declare -A LONG_NAME
while read -r short long; do
    LONG_NAME[$short]=$long
done < <("$DISCOAL" --help | sed -n 's/^  \(-[A-Za-z0-9]*\), \(--[a-z-]*\).*/\1 \2/p')
LONG_NAME[-ej]=--join
USED_LONG=""
check_long_names() {
    local name=$1 tok; shift
    local long=()
    for tok in "$@"; do
        if [ -n "${LONG_NAME[$tok]+set}" ]; then
            long+=("${LONG_NAME[$tok]}"); USED_LONG="$USED_LONG ${LONG_NAME[$tok]}"
        else
            long+=("$tok")
        fi
    done
    local outa outb
    # stdout without its first line, the command line; an error gives none.
    outa=$(timeout "$CASE_TIMEOUT" "$DISCOAL" "$@" 2>/dev/null | tail -n +2)
    outb=$(timeout "$CASE_TIMEOUT" "$DISCOAL" "${long[@]}" 2>/dev/null | tail -n +2)
    if [ -n "$outa" ] && [ "$outa" = "$outb" ]; then
        echo -e "${GREEN}PASS${NC} $name"; PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL${NC} $name (short and long names differ, or no output)"
        FAIL=$((FAIL + 1)); FAIL_NAMES+=("$name")
    fi
}
check_long_names "long_names_neutral" \
    20 2 10000 -t 10 -r 10 -g 5 50 -p 2 10 10 -M 1 -m 0 1 2 -en 0.1 0 0.5 \
    -eg 0.2 1 2 -en 0.5 1 1 -em 0.6 1 0 0.5 -A 2 0 0.05 -ed 1.0 1 0 -d 1 2
check_long_names "long_names_all_populations" \
    20 2 10000 -t 10 -r 10 -gr 0.5 50 -p 2 10 10 -M 1 -eG 0 2 -eL 0.3 -0.5 \
    -eM 0.5 0.2 -el 0.6 0 0 -ej 1.0 1 0 -v -d 1 2
check_long_names "long_names_admixture" \
    20 2 10000 -t 10 -p 3 10 5 5 -ea 0.1 0 1 2 0.3 -ed 1.0 2 1 -d 1 2
check_long_names "long_names_sweep" \
    20 2 10000 -t 10 -r 10 -ws 0.01 -a 500 -x 0.5 -f 0.05 -c 0.8 -h \
    -N 10000 -i 40 -d 1 2
check_long_names "long_names_recurrent_mutation" \
    20 2 10000 -t 10 -r 10 -wd 0.01 -a 500 -x 0.5 -uA 0.0001 -d 1 2
check_long_names "long_names_neutral_fixation" \
    20 2 10000 -t 10 -r 10 -wn 0.01 -x 0.5 -d 1 2
check_long_names "long_names_linked_sweep" \
    20 2 10000 -t 10 -r 10 -ls 0.01 20 -a 500 -d 1 2
check_long_names "long_names_linked_deterministic" \
    20 2 10000 -t 10 -r 10 -ld 0.01 20 -a 500 -d 1 2
check_long_names "long_names_linked_neutral" \
    20 2 10000 -t 10 -r 10 -ln 0.01 20 -d 1 2
check_long_names "long_names_recurrent" \
    20 2 10000 -t 10 -r 10 -R 2 -a 500 -d 1 2
check_long_names "long_names_recurrent_linked" \
    20 2 10000 -t 10 -r 10 -L 2 -a 500 -d 1 2
check_long_names "long_names_trajectory_file" \
    20 2 10000 -t 10 -r 10 -N 10000 -S "$REPO_ROOT/config_examples/trajectory_example.traj" -x 0.3 -d 1 2
check_long_names "long_names_priors" \
    20 3 10000 -Pt 5 15 -Pr 5 15 -ws 0.01 -Pa 100 500 -Pu 0 0.05 -Px 0.2 0.8 \
    -Pf 0.01 0.05 -Pc 0.6 0.9 -PuA 0 0.0001 -Pe1 0.1 0.2 0.5 1 -Pe2 0.3 0.4 0.5 1 -d 1 2
check_long_names "long_names_prior_rho_exponential" \
    20 3 10000 -t 10 -Pre 10 100 -d 1 2
check_long_names "long_names_config" \
    -Y "$REPO_ROOT/config_examples/basic_example.yaml" -d 1 2
LONG_TS_DIR=$(mktemp -d)
check_long_names "long_names_tree_sequence" \
    10 1 1000 -t 5 -r 5 -ts "$LONG_TS_DIR/a.trees" -F -d 1 2
rm -r "$LONG_TS_DIR"
# Every long name of the help text was used above.
missing=""
for long in "${LONG_NAME[@]}"; do
    case " $USED_LONG " in *" $long "*) ;; *) missing="$missing $long" ;; esac
done
if [ "${#LONG_NAME[@]}" -ge 50 ] && [ -z "$missing" ]; then
    echo -e "${GREEN}PASS${NC} long_names_all_used"; PASS=$((PASS + 1))
else
    echo -e "${RED}FAIL${NC} long_names_all_used (${#LONG_NAME[@]} names; not used:$missing)"
    FAIL=$((FAIL + 1)); FAIL_NAMES+=("long_names_all_used")
fi

run_fail_case "unknown_long_option" \
    "Unknown option '--thetas'" \
    10 1 100 --thetas 5

# --help and --version write to stdout and end with status 0; a write error
# on stdout ends with status 1.
check_status_stdout() {
    local name=$1 want_rc=$2 want_out=$3; shift 3
    local out rc
    out=$(timeout "$CASE_TIMEOUT" "$@" 2>/dev/null); rc=$?
    if [ "$rc" -eq "$want_rc" ] && { [ -z "$want_out" ] || printf '%s' "$out" | grep -qF -- "$want_out"; }; then
        echo -e "${GREEN}PASS${NC} $name"; PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL${NC} $name (exit $rc)"; FAIL=$((FAIL + 1)); FAIL_NAMES+=("$name")
    fi
}
check_status_stdout "help_to_stdout" 0 "usage: discoal" "$DISCOAL" --help
check_status_stdout "version_to_stdout" 0 "discoal " "$DISCOAL" --version
if [ -w /dev/full ]; then
    check_status_stdout "stdout_write_error" 1 "" sh -c "\"$DISCOAL\" 10 5 1000 -t 50 > /dev/full"
fi

# A -ts filename longer than the buffer stops with a message.
run_fail_case "ts_filename_too_long" \
    "the -ts filename is too long" \
    4 1 100 -t 1 -ts "$(printf 'a%.0s' $(seq 1 1100)).trees"

# The most negative int is a population index like any other, not a mark
# for "all populations".
run_fail_case "A_population_int_min" \
    "ancient sample of 1 lineages from population -2147483648" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -A 1 -2147483648 0.1

run_fail_case "en_population_int_min" \
    "population -2147483648 must be between 0 and 1" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -en 0.1 -2147483648 0.5

run_fail_case "em_population_int_min" \
    "population -2147483648 must be between 0 and 1" \
    4 1 100 -t 1 -p 2 2 2 -M 1 -em 0.5 0 -2147483648 1

# Times in messages are in units of 4N0 generations, as given.
run_fail_case "el_message_in_user_time" \
    "drives population 0 size to zero at time 0.25" \
    5 1 100 -t 5 -el 0.05 0 5.0 -en 1.0 0 1.0

# -eG, -eL and -eM apply to every population, also when -p comes after
# them: the output must equal that of the same options after -p.
check_same_output() {
    local name=$1 a=$2 b=$3
    local outa outb
    outa=$("$DISCOAL" $a 2>&1 | tail -n +2)
    outb=$("$DISCOAL" $b 2>&1 | tail -n +2)
    if [ -n "$outa" ] && [ "$outa" = "$outb" ]; then
        echo -e "${GREEN}PASS${NC} $name"
        PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL${NC} $name (outputs differ)"
        FAIL=$((FAIL + 1))
        FAIL_NAMES+=("$name")
    fi
}
check_same_output "eG_eL_eM_before_p" \
    "20 2 10000 -t 10 -r 10 -eG 0 2 -eL 0.3 -0.5 -eM 0.5 0.2 -p 2 10 10 -M 1 -d 1 2" \
    "20 2 10000 -t 10 -r 10 -p 2 10 10 -M 1 -eG 0 2 -eL 0.3 -0.5 -eM 0.5 0.2 -d 1 2"

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
# into population 1 strands lineages after population 1 merges). A fixed
# seed: with some seeds all lineages coalesce before the last event, and
# the run is then correct.
run_fail_case "R_lineages_can_never_meet" \
    "they can never coalesce" \
    20 2 10000 -t 10 -r 10 -R 0.5 -a 500 -p 3 7 7 6 -M 1 -ed 2 1 0 -ed 3 2 0 -d 1 7

# One single sweep at most: the sweep sampler follows one sweep at a time.
run_fail_case "two_ws" \
    "only one single sweep can be given" \
    20 5 1000 -t 5 -r 5 -ws 0.01 -ws 0.012 -a 500

run_fail_case "ws_and_ls" \
    "only one single sweep can be given" \
    20 5 1000 -t 5 -r 5 -ws 0.01 -ls 0.05 5 -a 500

# A demes model gives the demography and the YAML its samples, checked
# against the demes' lifetimes; options that change either would bypass
# these checks.
run_fail_case "demes_model_with_cli_p" \
    "-p cannot be used with a demes model" \
    -Y "$REPO_ROOT/config_examples/demes_example.yaml" -p 2 10 10

run_fail_case "demes_model_with_cli_A" \
    "-A cannot be used with a demes model" \
    -Y "$REPO_ROOT/config_examples/demes_example.yaml" -A 2 1 0.5

run_fail_case "demes_model_with_cli_ed" \
    "-ed cannot be used with a demes model" \
    -Y "$REPO_ROOT/config_examples/demes_example.yaml" -ed 0.1 1 0

run_fail_case "yaml_sweep_and_cli_ws" \
    "only one single sweep can be given" \
    -Y "$REPO_ROOT/testing/debug_yaml/sweep_with_tau.yaml" -ws 0.05

run_fail_case "Pu_with_recurrent_sweeps" \
    "needs a single sweep, not recurrent sweeps" \
    20 2 1000 -t 5 -R 0.1 -a 100 -Pu 0 0.1

run_fail_case "A_negative_time" \
    "the time must be >= 0" \
    4 2 100 -t 2 -A 1 0 -0.1

# Seeds outside the range of the generator's setall() stop with a message
# that names the option, not inside the generator.
run_fail_case "seed1_zero" \
    "the first seed (-d or YAML seed) (0) must be in [1, 2147483562]" \
    4 2 100 -t 2 -d 0 5

run_fail_case "seed2_above_bound" \
    "the second seed (-d or YAML seed) (2147483399) must be in [1, 2147483398]" \
    4 2 100 -t 2 -d 5 2147483399

run_pass_case "seeds_at_bounds" \
    4 2 100 -t 2 -d 2147483562 2147483398

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

# A soft sweep (-f) under many size changes of population 0: the model of
# issue #30, cut to its size changes before time 1.0. Its trajectory below
# f0 could fail to end, and the run stopped at the trajectory-length limit.
# It must run. With 8 replicates, the bug shows for most seeds.
run_pass_case "soft_sweep_many_size_changes" \
    10 8 1000 -t 1 -r 1 -ws 0 -a 800 -f 0.1 -x 0.5 \
    -en 0.0005124 0 1 -en 0.000605 0 1 \
    -en 0.0007144 0 1 -en 0.0008435 0 1.003 \
    -en 0.0009961 0 1.013 -en 0.001176 0 1.039 \
    -en 0.001389 0 1.086 -en 0.00164 0 1.163 \
    -en 0.001936 0 1.288 -en 0.002286 0 1.488 \
    -en 0.0027 0 1.77 -en 0.003188 0 2.122 \
    -en 0.003764 0 2.525 -en 0.004445 0 2.971 \
    -en 0.005249 0 3.457 -en 0.006198 0 3.981 \
    -en 0.007318 0 4.538 -en 0.008641 0 5.12 \
    -en 0.0102 0 5.721 -en 0.01205 0 6.332 \
    -en 0.01423 0 6.944 -en 0.0168 0 7.547 \
    -en 0.01984 0 8.13 -en 0.02342 0 8.685 \
    -en 0.02766 0 9.202 -en 0.03266 0 9.673 \
    -en 0.03856 0 10.09 -en 0.04554 0 10.45 \
    -en 0.05377 0 10.74 -en 0.06349 0 10.97 \
    -en 0.07497 0 11.13 -en 0.08852 0 11.21 \
    -en 0.1045 0 11.23 -en 0.1234 0 11.19 \
    -en 0.1457 0 11.08 -en 0.1721 0 10.92 \
    -en 0.2032 0 10.71 -en 0.2399 0 10.45 \
    -en 0.2833 0 10.15 -en 0.3346 0 9.817 \
    -en 0.395 0 9.459 -en 0.4665 0 9.08 \
    -en 0.5508 0 8.686 -en 0.6504 0 8.284 \
    -en 0.768 0 7.877 -en 0.9068 0 7.471 \
    -d 2 99

# More samples than the 1000 lineages that nodes[] starts with: the samples
# overflowed the array (glibc then stopped with "realloc(): invalid old
# size", or the output was wrong).
run_pass_case "more_than_1000_samples" \
    1500 2 1000 -t 5 -r 5 -d 3 4

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
    20 5 10000 -N 10000 -t 1 -r 10 -ws 0.01 -a 500 -el 0 0 -2

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
    "the number of sites (200000000) must be <= 100000000" \
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
    "Invalid numeric argument for -ws: 'nan'" \
    20 2 1000 -t 5 -r 5 -ws nan -a 500 -Pu 0.05 0.1

run_fail_case "x_nan" \
    "Invalid numeric argument for -x: 'nan'" \
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
