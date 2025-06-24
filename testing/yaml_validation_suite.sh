#!/bin/bash

# YAML Configuration Validation Suite
# Tests YAML configurations against equivalent command line arguments
# Ensures YAML config produces identical results to command line usage

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
DISCOAL_YAML="../../build/discoal"
DISCOAL_CMD="../../build/discoal"
TEST_DIR="yaml_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

echo -e "${BLUE}=== YAML Configuration Validation Suite ===${NC}"
echo -e "${BLUE}Testing YAML configs vs equivalent command line arguments${NC}"
echo -e "${BLUE}Test directory: $TEST_DIR${NC}"
echo

# Function to run a test case
run_test() {
    local test_name="$1"
    local yaml_config="$2"
    local cmd_args="$3"
    local description="$4"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo -e "${YELLOW}Test $TOTAL_TESTS: $test_name${NC}"
    echo "Description: $description"
    echo "Command line: $cmd_args"
    
    # Create YAML config file
    echo "$yaml_config" > "test_${TOTAL_TESTS}.yaml"
    
    # Run YAML version
    $DISCOAL_YAML -Y "test_${TOTAL_TESTS}.yaml" > "yaml_output_${TOTAL_TESTS}.txt" 2>&1
    yaml_exit_code=$?
    
    # Run command line version
    eval "$DISCOAL_CMD $cmd_args" > "cmd_output_${TOTAL_TESTS}.txt" 2>&1
    cmd_exit_code=$?
    
    # Check exit codes
    if [ $yaml_exit_code -ne 0 ] || [ $cmd_exit_code -ne 0 ]; then
        echo -e "${RED}FAIL: Exit code mismatch (YAML: $yaml_exit_code, CMD: $cmd_exit_code)${NC}"
        echo "YAML stderr:"
        cat "yaml_output_${TOTAL_TESTS}.txt"
        echo "CMD stderr:"
        cat "cmd_output_${TOTAL_TESTS}.txt"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        echo
        return 1
    fi
    
    # Extract simulation results (skip debug output)
    grep -A 1000 "^//" "yaml_output_${TOTAL_TESTS}.txt" | head -n 50 > "yaml_clean_${TOTAL_TESTS}.txt" 2>/dev/null || echo "No simulation output found" > "yaml_clean_${TOTAL_TESTS}.txt"
    grep -A 1000 "^//" "cmd_output_${TOTAL_TESTS}.txt" | head -n 50 > "cmd_clean_${TOTAL_TESTS}.txt" 2>/dev/null || echo "No simulation output found" > "cmd_clean_${TOTAL_TESTS}.txt"
    
    # Compare outputs
    if diff -q "yaml_clean_${TOTAL_TESTS}.txt" "cmd_clean_${TOTAL_TESTS}.txt" > /dev/null 2>&1; then
        echo -e "${GREEN}PASS: Outputs match${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}FAIL: Outputs differ${NC}"
        echo "YAML output preview:"
        head -5 "yaml_clean_${TOTAL_TESTS}.txt"
        echo "CMD output preview:"
        head -5 "cmd_clean_${TOTAL_TESTS}.txt"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    echo
}

# Test 1: Basic simulation parameters
run_test "basic_simulation" \
"simulation:
  sample_size: 20
  num_replicates: 3
  num_sites: 1000
  seed: [12345, 67890]

genetics:
  mutation_rate: 10.0

populations:
  count: 1
  sample_sizes: [20]" \
"20 3 1000 -t 10.0 -d 12345 67890" \
"Basic simulation with mutation rate"

# Test 2: Recombination rate
run_test "with_recombination" \
"simulation:
  sample_size: 15
  num_replicates: 2
  num_sites: 2000
  seed: [11111, 22222]

genetics:
  mutation_rate: 8.0
  recombination_rate: 4.0

populations:
  count: 1
  sample_sizes: [15]" \
"15 2 2000 -t 8.0 -r 4.0 -d 11111 22222" \
"Simulation with mutation and recombination"

# Test 3: Gene conversion
run_test "gene_conversion" \
"simulation:
  sample_size: 25
  num_replicates: 2
  num_sites: 1500
  seed: [98765, 54321]

genetics:
  mutation_rate: 12.0
  recombination_rate: 6.0
  gene_conversion:
    rate: 2.0
    tract_length: 500

populations:
  count: 1
  sample_sizes: [25]" \
"25 2 1500 -t 12.0 -r 6.0 -g 2.0 500 -d 98765 54321" \
"Simulation with gene conversion"

# Test 4: Gene conversion with crossover ratio
run_test "gene_conversion_co_ratio" \
"simulation:
  sample_size: 30
  num_replicates: 2
  num_sites: 3000
  seed: [13579, 24680]

genetics:
  mutation_rate: 15.0
  recombination_rate: 8.0
  gene_conversion:
    tract_length: 300
    crossover_ratio: 0.5" \
"30 2 3000 -t 15.0 -r 8.0 -gr 0.5 300 -d 13579 24680" \
"Gene conversion with crossover ratio"

# Test 5: Population size change
run_test "population_size_change" \
"simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 2000
  seed: [55555, 66666]

genetics:
  mutation_rate: 10.0

populations:
  count: 1
  sample_sizes: [20]

events:
  demographic:
    - type: size_change
      time: 0.1
      population: 0
      size: 2.0" \
"20 2 2000 -t 10.0 -en 0.1 0 2.0 -d 55555 66666" \
"Population size change event"

# Test 6: Selection sweep - stochastic
run_test "stochastic_sweep" \
"simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 5000
  seed: [77777, 88888]

genetics:
  mutation_rate: 15.0
  recombination_rate: 8.0

populations:
  count: 1
  sample_sizes: [20]

events:
  selection:
    - type: sweep
      mode: stochastic
      time: 0.01
      selection_coeff: 100.0
      position: 0.5" \
"20 2 5000 -t 15.0 -r 8.0 -ws 0.01 -a 100.0 -x 0.5 -d 77777 88888" \
"Stochastic selection sweep"

# Test 7: Selection sweep - deterministic
run_test "deterministic_sweep" \
"simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 4000
  seed: [11223, 33445]

genetics:
  mutation_rate: 12.0
  recombination_rate: 6.0

populations:
  count: 1
  sample_sizes: [20]

events:
  selection:
    - type: sweep
      mode: deterministic
      time: 0.02
      selection_coeff: 150.0
      position: 0.3" \
"20 2 4000 -t 12.0 -r 6.0 -wd 0.02 -a 150.0 -x 0.3 -d 11223 33445" \
"Deterministic selection sweep"

# Test 8: Soft sweep (initial frequency)
run_test "soft_sweep" \
"simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 3000
  seed: [99887, 77665]

genetics:
  mutation_rate: 10.0
  recombination_rate: 5.0

populations:
  count: 1
  sample_sizes: [20]

events:
  selection:
    - type: sweep
      mode: stochastic
      time: 0.015
      selection_coeff: 120.0
      position: 0.7
      initial_freq: 0.1" \
"20 2 3000 -t 10.0 -r 5.0 -ws 0.015 -a 120.0 -x 0.7 -f 0.1 -d 99887 77665" \
"Soft sweep with initial frequency"

# Test 9: Effective population size
run_test "effective_popn_size" \
"simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 2000
  seed: [44444, 55555]
  effective_popn_size: 500000

genetics:
  mutation_rate: 8.0

populations:
  count: 1
  sample_sizes: [20]

events:
  selection:
    - type: sweep
      mode: stochastic
      time: 0.01
      selection_coeff: 80.0" \
"20 2 2000 -t 8.0 -ws 0.01 -a 80.0 -N 500000 -d 44444 55555" \
"Custom effective population size with sweep"

# Test 10: Complex scenario with multiple features
run_test "complex_scenario" \
"simulation:
  sample_size: 40
  num_replicates: 2
  num_sites: 5000
  seed: [12121, 21212]

genetics:
  mutation_rate: 20.0
  recombination_rate: 10.0
  gene_conversion:
    rate: 4.0
    tract_length: 400

populations:
  count: 1
  sample_sizes: [40]

events:
  demographic:
    - type: size_change
      time: 0.05
      population: 0
      size: 1.5
  selection:
    - type: sweep
      mode: stochastic
      time: 0.02
      selection_coeff: 50.0
      position: 0.4" \
"40 2 5000 -t 20.0 -r 10.0 -g 4.0 400 -en 0.05 0 1.5 -ws 0.02 -a 50.0 -x 0.4 -d 12121 21212" \
"Complex scenario with gene conversion, size change, and sweep"

# Summary
echo -e "${BLUE}=== Test Summary ===${NC}"
echo -e "Total tests: $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
echo -e "${RED}Failed: $FAILED_TESTS${NC}"

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}All tests passed! YAML configuration system working correctly.${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed. YAML configuration needs debugging.${NC}"
    echo -e "${YELLOW}Check individual test outputs in $TEST_DIR for details.${NC}"
    exit 1
fi