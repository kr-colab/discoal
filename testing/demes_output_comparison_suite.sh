#!/bin/bash

# Demes output comparison test suite
# Tests that demographic models loaded via demes files produce identical results
# to the same models specified via command line arguments

# IMPORTANT NOTE ON TIME CONVERSIONS:
# Discoal has a confusing time scaling convention:
# - Command line times are specified in units of 2N generations
# - Internally, times are stored in units of 4N generations
# - Command line parser multiplies times by 2.0 (e.g., -en 0.05 becomes 0.1 internally)
# - Demes times are in generations, so we convert: generations / (2N) to get internal units
# 
# Example: 200,000 generations with N=1,000,000
# - Demes: 200,000 generations
# - Internal: 200,000 / (2 * 1,000,000) = 0.1
# - Command line equivalent: 0.1 / 2 = 0.05 (because parser will multiply by 2)

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

TEST_DIR="demes_comparison_$(date +%Y%m%d_%H%M%S)"
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
    local demes_file=$2
    local demes_cmd=$3
    local cmdline_cmd=$4
    
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
    
    # Remove the command line and debug output for comparison
    grep -v "^Loaded .* populations" "$TEST_DIR/outputs/${test_name}_demes.out" | \
        grep -v "^Note: You must use" | \
        grep -v "^Example:" | \
        grep -v "^===" | \
        grep -v "DEMES" | \
        grep -v "DISCOAL" | \
        grep -v "Reference N" | \
        grep -v "Number of populations" | \
        grep -v "Generation time" | \
        grep -v "Population [0-9]" | \
        grep -v "Start time" | \
        grep -v "epochs" | \
        grep -v "Ancestors" | \
        grep -v "Migrations:" | \
        grep -v "Pulses:" | \
        grep -v "Present-day" | \
        grep -v "currentSize" | \
        grep -v "Total events" | \
        grep -v "Has migration" | \
        grep -v "Validating" | \
        sed '1,2d' > "$TEST_DIR/outputs/${test_name}_demes_filtered.out"
    
    sed '1,2d' "$TEST_DIR/outputs/${test_name}_cmdline.out" > "$TEST_DIR/outputs/${test_name}_cmdline_filtered.out"
    
    # Compare outputs
    if diff -q "$TEST_DIR/outputs/${test_name}_demes_filtered.out" "$TEST_DIR/outputs/${test_name}_cmdline_filtered.out" > /dev/null; then
        echo -e "  ${GREEN}PASSED${NC}: Outputs are identical"
        ((PASSED++))
    else
        echo -e "  ${RED}FAILED${NC}: Outputs differ"
        echo "  First 10 lines of diff:"
        diff "$TEST_DIR/outputs/${test_name}_demes_filtered.out" "$TEST_DIR/outputs/${test_name}_cmdline_filtered.out" | head -10
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
    "$TEST_DIR/demes_files/single_pop_constant.yaml" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -D $TEST_DIR/demes_files/single_pop_constant.yaml -p 1 $NSAMP -d $SEED1 $SEED2" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -d $SEED1 $SEED2"

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
    "$TEST_DIR/demes_files/size_change.yaml" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -D $TEST_DIR/demes_files/size_change.yaml -p 1 $NSAMP -d $SEED1 $SEED2" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -en 0.01 0 2.0 -d $SEED1 $SEED2"

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
    "$TEST_DIR/demes_files/two_pop_migration.yaml" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -D $TEST_DIR/demes_files/two_pop_migration.yaml -p 2 10 10 -d $SEED1 $SEED2" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -M 4.0 -d $SEED1 $SEED2"

# Test 4: Population join (modeled as pulse in demes)
# The command line -ed 0.05 1 0 means at time 0.05, all lineages from pop1 move to pop0
# In demes, we model this as a pulse migration of proportion 1.0
cat > "$TEST_DIR/demes_files/pop_join.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 1000000}
pulses:
- {time: 200000, sources: [pop1], dest: pop0, proportions: [1.0]}
EOF

# After the pulse, we need to change pop0's size to reflect the ancestral size
# But demes doesn't allow size changes at pulse times, so we use a different approach
# Actually, let's model this correctly with splits instead
cat > "$TEST_DIR/demes_files/pop_split_correct.yaml" << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 200000, start_size: 1000000}
  - {end_time: 0, start_size: 2000000}
- name: pop1
  ancestors: [pop0]
  start_time: 200000
  epochs:
  - {end_time: 0, start_size: 1000000}
EOF

# N0 = 1000000 (present-day size of pop0 from epoch 1, but we use size from epoch 0 = 2000000)
# Actually wait, that's not right either. Let me think...

# The command line model is: pop0 and pop1 exist, at time 0.05 all lineages from pop1 go to pop0,
# and pop0's size changes to 2.0. This is effectively modeling a split backwards in time.
# In forward time: one population of size 2N splits into two populations of size N each.

# Test 4: Population split
# In demes: ancestor splits into pop0 and pop1
# In discoal command line: pop1 joins pop0 going backwards
# The key is ensuring the population indices match what discoal expects

# Test 4: Population split (corrected version)
# Model where deme0 exists throughout, deme1 splits from deme0
# This ensures population indices match what discoal expects
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
# Time conversion: 200,000 generations → 0.1 internal units → 0.05 command line units
# (200,000 / (2*N) = 0.1 internal, divided by 2 for command line = 0.05)
compare_outputs "pop_split_corrected" \
    "$TEST_DIR/demes_files/pop_split_corrected.yaml" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -D $TEST_DIR/demes_files/pop_split_corrected.yaml -p 2 10 10 -d $SEED1 $SEED2" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -ed 0.05 1 0 -en 0.05 0 2.0 -d $SEED1 $SEED2"

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

# Rate from 0->1: 2e-6 * 4 * 1000000 = 8.0
# Rate from 1->0: 1e-6 * 4 * 1000000 = 4.0
compare_outputs "asymmetric_migration" \
    "$TEST_DIR/demes_files/asymmetric_migration.yaml" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -D $TEST_DIR/demes_files/asymmetric_migration.yaml -p 2 10 10 -d $SEED1 $SEED2" \
    "../build/discoal $NSAMP $NREPS $NSITES -t $THETA -r $RHO -p 2 10 10 -m 0 1 8.0 -m 1 0 4.0 -d $SEED1 $SEED2"

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
    echo "Check the outputs in: $TEST_DIR/outputs/"
    exit 1
fi