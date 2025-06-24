#!/bin/bash

# Demes validation test suite for discoal
# Tests various demographic models specified in demes format

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create temporary directory for test files
TEST_DIR="demes_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p $TEST_DIR/demes_files

echo "=== DISCOAL DEMES VALIDATION SUITE ==="
echo "Creating test files in $TEST_DIR/"
echo
echo "IMPORTANT NOTES:"
echo "1. Current implementation requires -p flag even for single population models"
echo "   This should be improved in future versions"
echo "2. The -D flag must come BEFORE the -p flag in the command line"
echo "3. When time_units='generations', generation_time must be 1 (demes spec requirement)"
echo

# Counter for tests
PASSED=0
FAILED=0

# Function to run a test
run_test() {
    local test_name=$1
    local demes_file=$2
    local discoal_args=$3
    local expect_fail=$4
    
    echo -n "Testing $test_name... "
    
    # Parse the arguments to separate basic args from -p flag
    local basic_args=$(echo "$discoal_args" | sed 's/-p.*//')
    local p_flag=$(echo "$discoal_args" | grep -o -- '-p.*' || echo "")
    
    # Construct command with -D before -p
    local cmd="../build/discoal $basic_args -D $demes_file $p_flag"
    
    if OUTPUT=$($cmd 2>&1); then
        if [ "$expect_fail" = "true" ]; then
            echo -e "${RED}FAILED${NC} (expected to fail but succeeded)"
            ((FAILED++))
        else
            echo -e "${GREEN}PASSED${NC}"
            ((PASSED++))
        fi
    else
        if [ "$expect_fail" = "true" ]; then
            echo -e "${GREEN}PASSED${NC} (correctly rejected)"
            ((PASSED++))
        else
            echo -e "${RED}FAILED${NC}"
            echo "  Command: $cmd"
            echo "  Error output:"
            echo "$OUTPUT" | grep -E "(Error:|line|unknown)" | head -5
            ((FAILED++))
        fi
    fi
}

# Test 1: Single population, constant size
cat > $TEST_DIR/demes_files/single_pop_constant.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
EOF

run_test "Single population constant size" \
    "$TEST_DIR/demes_files/single_pop_constant.yaml" \
    "20 1 1000 -t 10 -p 1 20" \
    false

# Test 2: Two populations with migration
cat > $TEST_DIR/demes_files/two_pop_migration.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 15000}
migrations:
- {demes: [pop0, pop1], rate: 1e-4}
EOF

run_test "Two populations with migration" \
    "$TEST_DIR/demes_files/two_pop_migration.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    false

# Test 3: Two populations without migration (should fail)
cat > $TEST_DIR/demes_files/two_pop_isolated.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 15000}
EOF

run_test "Two isolated populations (should fail)" \
    "$TEST_DIR/demes_files/two_pop_isolated.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    true

# Test 4: Population split
cat > $TEST_DIR/demes_files/pop_split.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: ancestor
  epochs:
  - {end_time: 1000, start_size: 20000}
- name: pop0
  ancestors: [ancestor]
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  ancestors: [ancestor]
  epochs:
  - {end_time: 0, start_size: 15000}
EOF

run_test "Population split" \
    "$TEST_DIR/demes_files/pop_split.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    false

# Test 5: Size change event
cat > $TEST_DIR/demes_files/size_change.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 1000, start_size: 20000}
  - {end_time: 0, start_size: 5000}
EOF

run_test "Population size change" \
    "$TEST_DIR/demes_files/size_change.yaml" \
    "20 1 1000 -t 10 -p 1 20" \
    false

# Test 6: Complex three population model (like Out-of-Africa)
cat > $TEST_DIR/demes_files/three_pop_ooa.yaml << EOF
time_units: years
generation_time: 25
demes:
- name: YRI
  epochs:
  - {end_time: 0, start_size: 14474}
- name: OOA
  ancestors: [YRI]
  start_time: 52500
  epochs:
  - {end_time: 27500, start_size: 1861}
  - {end_time: 0, start_size: 12300}
- name: CEU
  ancestors: [OOA]
  start_time: 27500
  epochs:
  - {end_time: 0, start_size: 1032}
- name: CHB
  ancestors: [OOA]
  start_time: 27500
  epochs:
  - {end_time: 0, start_size: 554}
migrations:
- {demes: [YRI, OOA], rate: 25e-5, start_time: 52500, end_time: 27500}
- {demes: [YRI, CEU], rate: 3e-5, start_time: 27500}
- {demes: [YRI, CHB], rate: 1.9e-5, start_time: 27500}
- {demes: [CEU, CHB], rate: 9.6e-5, start_time: 27500}
EOF

run_test "Three population Out-of-Africa model" \
    "$TEST_DIR/demes_files/three_pop_ooa.yaml" \
    "30 1 1000 -t 10 -p 4 8 8 7 7" \
    false

# Test 7: Admixture event
cat > $TEST_DIR/demes_files/admixture.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: src1
  epochs:
  - {end_time: 500, start_size: 10000}
- name: src2
  epochs:
  - {end_time: 500, start_size: 15000}
- name: admixed
  ancestors: [src1, src2]
  proportions: [0.3, 0.7]
  start_time: 500
  epochs:
  - {end_time: 0, start_size: 12000}
EOF

run_test "Admixture model" \
    "$TEST_DIR/demes_files/admixture.yaml" \
    "20 1 1000 -t 10 -p 1 20" \
    false

# Test 8: Pulse migration (admixture)
cat > $TEST_DIR/demes_files/pulse_migration.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 15000}
migrations:
- {demes: [pop0, pop1], rate: 1e-4}
pulses:
- {sources: [pop1], dest: pop0, proportions: [0.1], time: 500}
EOF

run_test "Pulse migration" \
    "$TEST_DIR/demes_files/pulse_migration.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    false

# Test 9: Exponential growth (should fail)
cat > $TEST_DIR/demes_files/exp_growth.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000, end_size: 10000, size_function: exponential}
EOF

run_test "Exponential growth (should fail)" \
    "$TEST_DIR/demes_files/exp_growth.yaml" \
    "20 1 1000 -t 10" \
    true

# Test 10: Selfing rate (should fail)
cat > $TEST_DIR/demes_files/selfing.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000, selfing_rate: 0.5}
EOF

run_test "Selfing rate (should fail)" \
    "$TEST_DIR/demes_files/selfing.yaml" \
    "20 1 1000 -t 10" \
    true

# Test 11: Complex migration changes
cat > $TEST_DIR/demes_files/migration_changes.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 15000}
migrations:
- {demes: [pop0, pop1], rate: 1e-3, start_time: .inf, end_time: 1000}
- {demes: [pop0, pop1], rate: 5e-4, start_time: 1000, end_time: 500}
- {demes: [pop0, pop1], rate: 1e-4, start_time: 500, end_time: 0}
EOF

run_test "Time-varying migration rates" \
    "$TEST_DIR/demes_files/migration_changes.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    false

# Test 12: Multiple sequential splits
cat > $TEST_DIR/demes_files/sequential_splits.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: root
  epochs:
  - {end_time: 3000, start_size: 30000}
- name: branch1
  ancestors: [root]
  start_time: 3000
  epochs:
  - {end_time: 2000, start_size: 20000}
- name: pop0
  ancestors: [branch1]
  start_time: 2000
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  ancestors: [branch1]
  start_time: 2000
  epochs:
  - {end_time: 0, start_size: 15000}
- name: pop2
  ancestors: [root]
  start_time: 3000
  epochs:
  - {end_time: 0, start_size: 12000}
EOF

run_test "Sequential population splits" \
    "$TEST_DIR/demes_files/sequential_splits.yaml" \
    "30 1 1000 -t 10 -p 3 10 10 10" \
    false

# Test 13: Island model (all populations connected)
cat > $TEST_DIR/demes_files/island_model.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop2
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop3
  epochs:
  - {end_time: 0, start_size: 10000}
migrations:
- {demes: [pop0, pop1], rate: 1e-4}
- {demes: [pop0, pop2], rate: 1e-4}
- {demes: [pop0, pop3], rate: 1e-4}
- {demes: [pop1, pop2], rate: 1e-4}
- {demes: [pop1, pop3], rate: 1e-4}
- {demes: [pop2, pop3], rate: 1e-4}
EOF

run_test "Four population island model" \
    "$TEST_DIR/demes_files/island_model.yaml" \
    "40 1 1000 -t 10 -p 4 10 10 10 10" \
    false

# Test 14: Asymmetric migration
cat > $TEST_DIR/demes_files/asymmetric_migration.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 5000}
migrations:
- {source: pop0, dest: pop1, rate: 2e-4}
- {source: pop1, dest: pop0, rate: 1e-4}
EOF

run_test "Asymmetric migration" \
    "$TEST_DIR/demes_files/asymmetric_migration.yaml" \
    "20 1 1000 -t 10 -p 2 10 10" \
    false

# Test 15: Complex pulse admixture
cat > $TEST_DIR/demes_files/complex_pulses.yaml << EOF
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 10000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 15000}
- name: pop2
  epochs:
  - {end_time: 0, start_size: 12000}
migrations:
- {demes: [pop0, pop1], rate: 1e-4}
- {demes: [pop1, pop2], rate: 1e-4}
pulses:
- {sources: [pop0, pop1], dest: pop2, proportions: [0.2, 0.3], time: 800}
- {sources: [pop2], dest: pop0, proportions: [0.1], time: 400}
EOF

run_test "Complex pulse admixture" \
    "$TEST_DIR/demes_files/complex_pulses.yaml" \
    "30 1 1000 -t 10 -p 3 10 10 10" \
    false

# Summary
echo
echo "=== TEST SUMMARY ==="
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"
echo

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi