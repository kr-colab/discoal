#!/usr/bin/env bash
#
# YAML example regression suite.
#
# Runs every example YAML in config_examples/ (except all_options.yaml,
# which is a non-runnable parser-options index) and asserts:
#
#   1. discoal exits 0.
#   2. For ms-style examples: stdout contains exactly num_replicates
#      `//` blocks, each followed by a `segsites:` line. For replicates
#      with segsites > 0: a `positions:` line with exactly that many
#      floats in [0, 1], followed by sample_size genotype rows of length
#      segsites over {0, 1}. Upper bound is inclusive: at `%6.6f` print
#      precision, true values just under 1 can round up to "1.000000".
#   3. For tree-sequence examples: the declared .trees output file is
#      created and non-empty.
#
# Exit code: 0 if every fixture passes, 1 otherwise. Never use set -e;
# we want all fixtures to run so the final report is a complete picture.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DISCOAL="${REPO_ROOT}/build/discoal"
EXAMPLES_DIR="${REPO_ROOT}/config_examples"
WORK_DIR="$(mktemp -d -t yaml_regression.XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT

PER_FIXTURE_TIMEOUT=60

RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
YELLOW=$'\033[1;33m'
NC=$'\033[0m'

pass_count=0
fail_count=0
failed_fixtures=()

# Declare each fixture's expected behavior. The harness is the source of
# truth for what "correct" means per fixture; if you change a YAML's
# sample_size / num_replicates, update the matching entry here.
#
# Format (space-separated, per line, matching $fixture):
#   FIXTURE_MODE:       ms | trees
#   FIXTURE_SAMPLE:     per-replicate sample count
#   FIXTURE_REPS:       expected number of `//` blocks (ms) or 1 (trees)
#   FIXTURE_TREES_FILE: .trees output path (trees mode only, otherwise "")
declare -A FIXTURE_MODE=(
    [basic_example.yaml]="ms"
    [sweep_example.yaml]="ms"
    [demographic_example.yaml]="ms"
    [demes_example.yaml]="ms"
    [tree_sequence_example.yaml]="trees"
)
declare -A FIXTURE_SAMPLE=(
    [basic_example.yaml]=10
    [sweep_example.yaml]=15
    [demographic_example.yaml]=20
    [demes_example.yaml]=20
    [tree_sequence_example.yaml]=20
)
declare -A FIXTURE_REPS=(
    [basic_example.yaml]=2
    [sweep_example.yaml]=10
    [demographic_example.yaml]=2
    [demes_example.yaml]=10
    [tree_sequence_example.yaml]=5
)
declare -A FIXTURE_TREES_FILE=(
    [tree_sequence_example.yaml]="tree_sequence_example.trees"
)
# Auxiliary files referenced by a fixture YAML via a `config_examples/<name>`
# relative path. The harness copies each into ${fixture_workdir}/config_examples/
# so the YAML's relative path resolves inside the per-fixture tempdir.
# Values are space-separated basenames under EXAMPLES_DIR.
declare -A FIXTURE_AUX_FILES=(
    [demes_example.yaml]="demes_example.demes.yaml"
)

fail() {
    local fixture="$1" reason="$2"
    echo "${RED}FAIL${NC} ${fixture}: ${reason}"
    fail_count=$((fail_count + 1))
    failed_fixtures+=("${fixture}")
}

pass() {
    local fixture="$1"
    echo "${GREEN}PASS${NC} ${fixture}"
    pass_count=$((pass_count + 1))
}

# Verify one ms-style replicate starting at the line index after `//`.
# Arguments: fixture_name, stdout_file, line_index_of_slashslash, sample_size
# Sets REPLICATE_NEXT_LINE to the line after this replicate's last genotype row.
# Returns 0 on success, non-zero with an error message printed on failure.
check_ms_replicate() {
    local fixture="$1" stdout_file="$2" start_idx="$3" sample_size="$4"
    local segsites_line positions_line segsites positions_count
    local -a genotype_rows

    # Line $start_idx is `//`. Line $start_idx + 1 must be `segsites: N`.
    segsites_line=$(sed -n "$((start_idx + 1))p" "${stdout_file}")
    if [[ ! "${segsites_line}" =~ ^segsites:\ ([0-9]+)$ ]]; then
        fail "${fixture}" "expected 'segsites: N' after '//' at line $((start_idx + 1)), got: '${segsites_line}'"
        return 1
    fi
    segsites="${BASH_REMATCH[1]}"

    if [[ "${segsites}" -eq 0 ]]; then
        # segsites=0 replicates have no positions/genotype lines.
        REPLICATE_NEXT_LINE=$((start_idx + 2))
        return 0
    fi

    positions_line=$(sed -n "$((start_idx + 2))p" "${stdout_file}")
    if [[ ! "${positions_line}" =~ ^positions:\  ]]; then
        fail "${fixture}" "expected 'positions: ...' at line $((start_idx + 2)), got: '${positions_line:0:80}'"
        return 1
    fi

    # Validate position count and that each is in [0, 1].
    # strip leading "positions: ", then count whitespace-separated tokens.
    # Upper bound is inclusive because discoal prints positions at %6.6f,
    # so a true value slightly under 1 can display as "1.000000".
    local positions_values
    positions_values="${positions_line#positions: }"
    positions_count=$(echo "${positions_values}" | awk '{print NF}')
    if [[ "${positions_count}" -ne "${segsites}" ]]; then
        fail "${fixture}" "positions count ${positions_count} != segsites ${segsites}"
        return 1
    fi
    if ! echo "${positions_values}" | awk '
        { for (i = 1; i <= NF; i++) {
            if ($i + 0 < 0 || $i + 0 > 1) { print "out_of_range:" $i; exit 1 }
        }}' >/dev/null; then
        fail "${fixture}" "positions contain value(s) outside [0, 1]"
        return 1
    fi

    # sample_size genotype rows, each length segsites, over {0, 1}.
    local row_idx row_start row
    row_start=$((start_idx + 3))
    for (( row_idx = 0; row_idx < sample_size; row_idx++ )); do
        row=$(sed -n "$((row_start + row_idx))p" "${stdout_file}")
        if [[ "${#row}" -ne "${segsites}" ]]; then
            fail "${fixture}" "genotype row $((row_idx + 1)) length ${#row} != segsites ${segsites}"
            return 1
        fi
        if [[ "${row}" =~ [^01] ]]; then
            fail "${fixture}" "genotype row $((row_idx + 1)) contains non-{0,1} chars: '${row}'"
            return 1
        fi
    done

    REPLICATE_NEXT_LINE=$((row_start + sample_size))
    return 0
}

check_ms_fixture() {
    local fixture="$1" stdout_file="$2" sample_size="$3" expected_reps="$4"

    # Count `//` lines; must equal expected_reps. grep -c exits 1 when the
    # match count is zero but still prints "0", so || : absorbs that under set -u.
    local actual_reps
    actual_reps=$(grep -c '^//$' "${stdout_file}" || :)
    if [[ "${actual_reps}" -ne "${expected_reps}" ]]; then
        fail "${fixture}" "found ${actual_reps} '//' blocks, expected ${expected_reps}"
        return 1
    fi

    # Walk each replicate in order.
    local -a slash_line_numbers
    mapfile -t slash_line_numbers < <(grep -n '^//$' "${stdout_file}" | cut -d: -f1)

    local idx
    for idx in "${slash_line_numbers[@]}"; do
        if ! check_ms_replicate "${fixture}" "${stdout_file}" "${idx}" "${sample_size}"; then
            return 1
        fi
    done

    return 0
}

check_trees_fixture() {
    local fixture="$1" trees_base_path="$2" expected_reps="$3"

    # discoal writes one tree-sequence file per replicate, derived from
    # tree_sequence_filename by stripping the ".trees" suffix and appending
    # "_rep<N>.trees" for N = 1..num_replicates. Require that every expected
    # per-replicate file exists and is non-empty.
    local base="${trees_base_path%.trees}"
    local n
    for (( n = 1; n <= expected_reps; n++ )); do
        local per_rep="${base}_rep${n}.trees"
        if [[ ! -f "${per_rep}" ]]; then
            fail "${fixture}" "expected tree-sequence file '${per_rep}' was not created"
            return 1
        fi
        if [[ ! -s "${per_rep}" ]]; then
            fail "${fixture}" "tree-sequence file '${per_rep}' is empty"
            return 1
        fi
    done

    return 0
}

# Pull the CLI invocation declared by a fixture YAML's load-bearing
# header. The YAML must contain exactly one line of the form
#   # Equivalent to: discoal <args...>
# We strip the prefix (including the literal "discoal " token, since
# the harness substitutes the absolute binary path) and print what is
# left. An empty result means the marker was not found.
parse_cli_args_from_yaml() {
    local yaml_path="$1"
    grep -m1 '^# Equivalent to: discoal ' "$yaml_path" \
        | sed 's|^# Equivalent to: discoal ||'
}

# Diff the YAML and CLI runs' stdouts. discoal echoes its own argv on
# the first stdout line, which is naturally different between -Y and
# the positional invocation, so we skip line 1 on both sides. Any
# remaining difference is a parity failure.
# Per-replicate trees parity. Walks each expected replicate file
# under the YAML and CLI workdirs in parallel and shells out to
# cmp_trees.py, which loads each .trees with tskit and compares
# table collections with provenance ignored (discoal stamps a
# timestamp into provenance, so byte-equal would never hold).
# Returns 0 if every replicate matches, 1 on first mismatch.
# Caller is responsible for deciding whether to invoke this — when
# tskit is not installed, the suite skips it altogether.
check_parity_trees() {
    local fixture="$1" yaml_dir="$2" cli_dir="$3"
    local expected_reps="$4" trees_basename="$5"
    local base="${trees_basename%.trees}"
    local n yaml_rep cli_rep cmp_out cmp_rc
    for (( n = 1; n <= expected_reps; n++ )); do
        yaml_rep="${yaml_dir}/${base}_rep${n}.trees"
        cli_rep="${cli_dir}/${base}_rep${n}.trees"
        cmp_out=$(python3 "${SCRIPT_DIR}/cmp_trees.py" \
                    "${yaml_rep}" "${cli_rep}" 2>&1)
        cmp_rc=$?
        if [[ "${cmp_rc}" -ne 0 ]]; then
            local detail="${cmp_out:-no comparator output}"
            fail "${fixture}" "trees parity mismatch on replicate ${n} (rc=${cmp_rc}: ${detail})"
            return 1
        fi
    done
    return 0
}

check_parity_ms_stdout() {
    local fixture="$1" yaml_stdout="$2" cli_stdout="$3"
    local diff_output
    diff_output=$(diff <(tail -n +2 "${yaml_stdout}") \
                      <(tail -n +2 "${cli_stdout}"))
    if [[ -z "${diff_output}" ]]; then
        return 0
    fi
    # Compact summary: total diff line count plus the first three diff
    # lines truncated to 60 chars each. Genotype/positions rows are
    # hundreds of chars wide, so we never let a single line dominate.
    local diff_lines diff_head
    diff_lines=$(echo "${diff_output}" | wc -l)
    diff_head=$(echo "${diff_output}" | head -n 3 | cut -c1-60 | tr '\n' '|')
    fail "${fixture}" "YAML/CLI stdout parity mismatch (${diff_lines} diff lines; first: ${diff_head})"
    return 1
}

run_fixture() {
    local fixture="$1"
    local yaml_path="${EXAMPLES_DIR}/${fixture}"
    local mode="${FIXTURE_MODE[${fixture}]}"
    local sample_size="${FIXTURE_SAMPLE[${fixture}]}"
    local expected_reps="${FIXTURE_REPS[${fixture}]}"

    if [[ ! -f "${yaml_path}" ]]; then
        fail "${fixture}" "fixture YAML missing at ${yaml_path}"
        return
    fi

    # The CLI form is the one declared in the YAML header — single
    # source of truth, no parallel array in this script.
    local cli_args
    cli_args=$(parse_cli_args_from_yaml "${yaml_path}")
    if [[ -z "${cli_args}" ]]; then
        fail "${fixture}" "YAML missing '# Equivalent to: discoal ...' header"
        return
    fi

    # Run -Y and the equivalent CLI invocation in separate workdirs so
    # tree-sequence fixtures' .trees files don't collide.
    local fixture_workdir="${WORK_DIR}/${fixture%.yaml}"
    local yaml_dir="${fixture_workdir}/yaml"
    local cli_dir="${fixture_workdir}/cli"
    mkdir -p "${yaml_dir}" "${cli_dir}"

    # Aux files that a fixture YAML cites by relative path (e.g. a
    # demes file) must be reachable from each run's CWD. The CLI form
    # generally doesn't need them, but copying into both keeps the
    # setup symmetric.
    local aux_files="${FIXTURE_AUX_FILES[${fixture}]:-}"
    if [[ -n "${aux_files}" ]]; then
        local d
        for d in "${yaml_dir}" "${cli_dir}"; do
            mkdir -p "${d}/config_examples"
            for aux in ${aux_files}; do
                cp "${EXAMPLES_DIR}/${aux}" "${d}/config_examples/${aux}"
            done
        done
    fi

    local yaml_stdout="${yaml_dir}/stdout.txt"
    local yaml_stderr="${yaml_dir}/stderr.txt"
    local cli_stdout="${cli_dir}/stdout.txt"
    local cli_stderr="${cli_dir}/stderr.txt"

    (
        cd "${yaml_dir}" && \
        timeout "${PER_FIXTURE_TIMEOUT}" "${DISCOAL}" -Y "${yaml_path}" \
            > "${yaml_stdout}" 2> "${yaml_stderr}"
    )
    local yaml_rc=$?
    if [[ "${yaml_rc}" -ne 0 ]]; then
        local stderr_head
        stderr_head=$(head -n 3 "${yaml_stderr}" | tr '\n' ' ')
        fail "${fixture}" "discoal -Y exited ${yaml_rc} (stderr: ${stderr_head})"
        return
    fi

    # Word-split ${cli_args} on whitespace into argv. None of the
    # current fixtures need quoted arguments; if a future fixture
    # does, this is where to revisit.
    (
        cd "${cli_dir}" && \
        timeout "${PER_FIXTURE_TIMEOUT}" "${DISCOAL}" ${cli_args} \
            > "${cli_stdout}" 2> "${cli_stderr}"
    )
    local cli_rc=$?
    if [[ "${cli_rc}" -ne 0 ]]; then
        local stderr_head
        stderr_head=$(head -n 3 "${cli_stderr}" | tr '\n' ' ')
        fail "${fixture}" "discoal (CLI form) exited ${cli_rc} (stderr: ${stderr_head})"
        return
    fi

    # Run all relevant checks; only mark the fixture as passed if every
    # one of them passes.
    local fixture_failed=0
    case "${mode}" in
        ms)
            if ! check_parity_ms_stdout "${fixture}" \
                "${yaml_stdout}" "${cli_stdout}"; then
                fixture_failed=1
            fi
            if ! check_ms_fixture "${fixture}" "${yaml_stdout}" \
                "${sample_size}" "${expected_reps}"; then
                fixture_failed=1
            fi
            ;;
        trees)
            local trees_basename="${FIXTURE_TREES_FILE[${fixture}]}"
            local yaml_trees="${yaml_dir}/${trees_basename}"
            if ! check_trees_fixture "${fixture}" "${yaml_trees}" \
                "${expected_reps}"; then
                fixture_failed=1
            fi
            # Parity is only meaningful if both the YAML run and the
            # CLI run produced the expected files; the existing
            # non-empty check above already covers the YAML side.
            if [[ "${TSKIT_AVAILABLE}" -eq 1 ]]; then
                if ! check_parity_trees "${fixture}" \
                    "${yaml_dir}" "${cli_dir}" \
                    "${expected_reps}" "${trees_basename}"; then
                    fixture_failed=1
                fi
            else
                echo "${YELLOW}SKIP${NC} ${fixture}: trees-parity (tskit unavailable)"
            fi
            ;;
        *)
            fail "${fixture}" "unknown mode '${mode}' in harness config"
            return
            ;;
    esac

    if [[ "${fixture_failed}" -eq 0 ]]; then
        pass "${fixture}"
    fi
}

# Build discoal if it doesn't exist yet so the harness is self-bootstrapping.
if [[ ! -x "${DISCOAL}" ]]; then
    echo "${YELLOW}build/discoal not found; running make...${NC}"
    ( cd "${REPO_ROOT}" && make ) || {
        echo "${RED}make failed; aborting${NC}"
        exit 1
    }
fi

# Probe for the trees-parity comparator. tskit is the heavy dep —
# environments without it just skip the trees parity check rather
# than failing the whole suite.
TSKIT_AVAILABLE=0
if command -v python3 >/dev/null 2>&1 \
   && python3 -c "import tskit" >/dev/null 2>&1 \
   && [[ -f "${SCRIPT_DIR}/cmp_trees.py" ]]; then
    TSKIT_AVAILABLE=1
fi

echo "Running YAML example regression suite against ${DISCOAL}"
echo

# Sort fixture names for stable report ordering. We must not pipe the
# loop through `sort`, because that would run the loop in a subshell
# and discard pass_count / fail_count / failed_fixtures.
readarray -t sorted_fixtures < <(printf '%s\n' "${!FIXTURE_MODE[@]}" | sort)
for fixture in "${sorted_fixtures[@]}"; do
    run_fixture "${fixture}"
done

echo
echo "---"
echo "${pass_count} passed, ${fail_count} failed"
if [[ "${fail_count}" -gt 0 ]]; then
    echo "Failed fixtures:"
    for f in "${failed_fixtures[@]}"; do
        echo "  - ${f}"
    done
    exit 1
fi
exit 0
