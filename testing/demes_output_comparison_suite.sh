#!/bin/bash

# Demes output comparison test suite
# Tests that demographic models loaded via demes files produce identical results
# to the same models specified via command line arguments

# Time units: command-line times are in units of 4N0 generations, and
# demes times in generations, so a demes time of g generations is the
# command-line time g / (4 N0). N0 is the present-day size of the first
# deme. The demes models go through a YAML config (-Y with
# demography.demes_filename), written by demes_run.

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
DISCOAL=${DISCOAL:-$SCRIPT_DIR/../build/discoal}
source "$SCRIPT_DIR/demes_run.sh"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test files go in a temporary directory, removed on exit.
TEST_DIR=$(mktemp -d -t demes_comparison.XXXXXX)
trap 'rm -rf "$TEST_DIR"' EXIT
mkdir -p "$TEST_DIR/demes_files"
mkdir -p "$TEST_DIR/outputs"

echo "=== DEMES OUTPUT COMPARISON SUITE ==="
echo "Testing equivalence between demes files and command-line demographic models"
echo "Test directory: $TEST_DIR"
echo ""

# Counter for tests
PASSED=0
FAILED=0

# Function to compare outputs
compare_outputs() {
    local test_name=$1
    local demes_cmd=$2
    local cmdline_cmd=$3
    
    echo "=== Test: $test_name ==="
    echo "  Demes command: $demes_cmd"
    echo "  Cmdline command: $cmdline_cmd"
    
    # Run with demes file
    if ! $demes_cmd > "$TEST_DIR/outputs/${test_name}_demes.out" 2> "$TEST_DIR/outputs/${test_name}_demes.err"; then
        echo -e "  ${RED}FAILED${NC}: Demes command failed"
        cat "$TEST_DIR/outputs/${test_name}_demes.err"
        ((FAILED++))
        return
    fi
    
    # Run with command line args
    if ! $cmdline_cmd > "$TEST_DIR/outputs/${test_name}_cmdline.out" 2> "$TEST_DIR/outputs/${test_name}_cmdline.err"; then
        echo -e "  ${RED}FAILED${NC}: Command line version failed"
        cat "$TEST_DIR/outputs/${test_name}_cmdline.err"
        ((FAILED++))
        return
    fi
    
    # The outputs must be identical apart from line 1, the echoed command.
    local demes_out="$TEST_DIR/outputs/${test_name}_demes.out"
    local cmdline_out="$TEST_DIR/outputs/${test_name}_cmdline.out"
    if diff -q <(sed 1d "$demes_out") <(sed 1d "$cmdline_out") > /dev/null; then
        echo -e "  ${GREEN}PASSED${NC}: Outputs are identical"
        ((PASSED++))
    else
        echo -e "  ${RED}FAILED${NC}: Outputs differ"
        echo "  First 10 lines of diff:"
        diff <(sed 1d "$demes_out") <(sed 1d "$cmdline_out") | head -10
        ((FAILED++))
    fi
    echo ""
}

# Fixed seeds for reproducibility
SEED1=12345
SEED2=67890

# Common simulation parameters
NSAMP=20
NREPS=5
NSITES=10000
THETA=10
RHO=5

# Test 1: Single population, constant size
# The trick here is that discoal's default EFFECTIVE_POPN_SIZE = 1000000
# So for demes to match, we need to use that same size
cat > "$TEST_DIR/demes_files/single_pop_constant.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000000}
EOF

compare_outputs "single_pop_constant" \
    "demes_run $TEST_DIR/demes_files/single_pop_constant.yaml $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 1 $NSAMP -d $SEED1 $SEED2" \
    "$DISCOAL $NSAMP $NREPS $NSITES -t $THETA -r $RHO -d $SEED1 $SEED2"

# Test 2: Population size change
# For this to work, we need the same effective population size trajectory
# Without demes: starts at size 1.0, changes to 2.0 at time 0.01
# With demes: we need present size that when divided by itself gives 1.0
cat > "$TEST_DIR/demes_files/size_change.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 40000, start_size: 2000000}
  - {end_time: 0, start_size: 1000000}
EOF

# With demes: N0 = 1000000 (present), ancient = 2000000
# Time: 40000/(4*1000000) = 0.01
# Relative ancient size: 2000000/1000000 = 2.0
compare_outputs "size_change" \
    "demes_run $TEST_DIR/demes_files/size_change.yaml $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 1 $NSAMP -d $SEED1 $SEED2" \
    "$DISCOAL $NSAMP $NREPS $NSITES -t $THETA -r $RHO -en 0.01 0 2.0 -d $SEED1 $SEED2"

# Test 3: Two populations with migration
cat > "$TEST_DIR/demes_files/two_pop_migration.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 1000000}
migrations:
- {demes: [pop0, pop1], rate: 1e-6}
EOF

# Migration rate: 1e-6 * 4 * 1000000 = 4.0
compare_outputs "two_pop_migration" \
    "demes_run $TEST_DIR/demes_files/two_pop_migration.yaml $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -d $SEED1 $SEED2" \
    "$DISCOAL $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -M 4.0 -d $SEED1 $SEED2"

# Test 4: Population split: deme1 splits from deme0, whose size changes then
cat > "$TEST_DIR/demes_files/pop_split_corrected.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: deme0
  epochs:
  - {end_time: 200000, start_size: 2000000}
  - {end_time: 0, start_size: 1000000}
- name: deme1
  ancestors: [deme0]
  start_time: 200000
  epochs:
  - {end_time: 0, start_size: 1000000}
EOF

# The demes model: deme0 has size 2M until 200k gen ago, then 1M; deme1 splits at 200k gen ago
# The command line: at time 0.05, pop1 joins pop0, and pop0 size changes to 2.0
# Time conversion: 200,000 generations / (4 N0) = 0.05 command-line units
compare_outputs "pop_split_corrected" \
    "demes_run $TEST_DIR/demes_files/pop_split_corrected.yaml $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -d $SEED1 $SEED2" \
    "$DISCOAL $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -ed 0.05 1 0 -en 0.05 0 2.0 -d $SEED1 $SEED2"

# Test 5: Asymmetric migration
cat > "$TEST_DIR/demes_files/asymmetric_migration.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 1000000}
migrations:
- {source: pop0, dest: pop1, rate: 2e-6}
- {source: pop1, dest: pop0, rate: 1e-6}
EOF

# A demes migration from source to dest is forward in time, so lineages
# move from dest to source going back, which is -m dest source:
# pop0 -> pop1 at 2e-6 is -m 1 0 8.0 (2e-6 * 4 * 1000000), and
# pop1 -> pop0 at 1e-6 is -m 0 1 4.0.
compare_outputs "asymmetric_migration" \
    "demes_run $TEST_DIR/demes_files/asymmetric_migration.yaml $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -d $SEED1 $SEED2" \
    "$DISCOAL $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -m 1 0 8.0 -m 0 1 4.0 -d $SEED1 $SEED2"

# Summary
echo "=== TEST SUMMARY ==="
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}All output comparison tests passed!${NC}"
    echo "Demes file loading produces identical results to command line specification"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi