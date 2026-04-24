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

    pass "${fixture}"
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

    pass "${fixture}"
    return 0
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

    local fixture_workdir="${WORK_DIR}/${fixture%.yaml}"
    mkdir -p "${fixture_workdir}"
    local stdout_file="${fixture_workdir}/stdout.txt"
    local stderr_file="${fixture_workdir}/stderr.txt"

    local aux_files="${FIXTURE_AUX_FILES[${fixture}]:-}"
    if [[ -n "${aux_files}" ]]; then
        mkdir -p "${fixture_workdir}/config_examples"
        for aux in ${aux_files}; do
            cp "${EXAMPLES_DIR}/${aux}" "${fixture_workdir}/config_examples/${aux}"
        done
    fi

    # Tree-sequence outputs land next to the working directory so we can
    # pick them up by the filename declared in the YAML.
    (
        cd "${fixture_workdir}" && \
        timeout "${PER_FIXTURE_TIMEOUT}" "${DISCOAL}" -Y "${yaml_path}" \
            > "${stdout_file}" 2> "${stderr_file}"
    )
    local rc=$?
    if [[ "${rc}" -ne 0 ]]; then
        local stderr_head
        stderr_head=$(head -n 3 "${stderr_file}" | tr '\n' ' ')
        fail "${fixture}" "discoal exited ${rc} (stderr: ${stderr_head})"
        return
    fi

    case "${mode}" in
        ms)
            check_ms_fixture "${fixture}" "${stdout_file}" \
                "${sample_size}" "${expected_reps}"
            ;;
        trees)
            local trees_file="${fixture_workdir}/${FIXTURE_TREES_FILE[${fixture}]}"
            check_trees_fixture "${fixture}" "${trees_file}" "${expected_reps}"
            ;;
        *)
            fail "${fixture}" "unknown mode '${mode}' in harness config"
            ;;
    esac
}

# Build discoal if it doesn't exist yet so the harness is self-bootstrapping.
if [[ ! -x "${DISCOAL}" ]]; then
    echo "${YELLOW}build/discoal not found; running make...${NC}"
    ( cd "${REPO_ROOT}" && make ) || {
        echo "${RED}make failed; aborting${NC}"
        exit 1
    }
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
