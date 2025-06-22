#!/bin/bash

# demographic_events_validation.sh
# Statistical validation of demographic event handling (admixture, ancient samples, etc.)
# Compares tskit-integrated version against legacy baseline using nicestats

set -e

# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"
WORK_DIR="$ROOT_DIR/work/demographic_validation_$(date +%Y%m%d_%H%M%S)"
NREPS=${1:-100}  # Number of replicates (default 100)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Ensure we're in the root directory
cd "$ROOT_DIR"

echo "=========================================="
echo "Demographic Events Statistical Validation"
echo "=========================================="
echo "Replicates: $NREPS"
echo "Work directory: $WORK_DIR"
echo

# Create work directory
mkdir -p "$WORK_DIR"

# Build executables
echo "Building executables..."
make discoal_edited > /dev/null 2>&1
make discoal_legacy_backup > /dev/null 2>&1
make niceStats > /dev/null 2>&1

# Copy executables to work directory for consistency
cp build/discoal_edited "$WORK_DIR/"
cp build/discoal_legacy_backup "$WORK_DIR/"
cp build/niceStats "$WORK_DIR/"

cd "$WORK_DIR"

# Function to run simulation and analyze with nicestats
run_and_analyze() {
    local exec_name=$1
    local output_prefix=$2
    local params="${@:3}"
    
    echo -n "Running $exec_name with: $params ... "
    
    # Run simulation
    ./$exec_name $params > ${output_prefix}.ms 2>${output_prefix}.err
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}FAILED${NC}"
        echo "Error output:"
        cat ${output_prefix}.err
        return 1
    fi
    
    # Check for -1 values (lost ancestral material)
    if grep -q "^-1" ${output_prefix}.ms; then
        echo -e "${RED}FAILED${NC} - Found -1 values (lost ancestral material)"
        return 1
    fi
    
    # Run nicestats
    ./niceStats < ${output_prefix}.ms > ${output_prefix}.stats 2>${output_prefix}.stats.err
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}FAILED${NC} - nicestats error"
        cat ${output_prefix}.stats.err
        return 1
    fi
    
    echo -e "${GREEN}OK${NC}"
    return 0
}

# Function to compare statistics
compare_stats() {
    local test_name=$1
    local legacy_stats=$2
    local edited_stats=$3
    local output_file=$4
    
    echo "=== $test_name ===" >> $output_file
    echo >> $output_file
    
    # Extract all statistics from nicestats output
    # nicestats outputs: pi, ss, D, thetaH, H, tajd, fusf, fusfs, rm, rmmg, nhaps, hdiv, wallb, wallq, rosasR, sweepK
    
    # Create temporary files for comparison
    grep -E "^pi:|^ss:|^D:|^thetaH:|^H:|^tajd:|^fusf:|^fusfs:|^rm:|^rmmg:|^nhaps:|^hdiv:|^wallb:|^wallq:|^rosasR:|^sweepK:" $legacy_stats | sort > legacy_sorted.tmp
    grep -E "^pi:|^ss:|^D:|^thetaH:|^H:|^tajd:|^fusf:|^fusfs:|^rm:|^rmmg:|^nhaps:|^hdiv:|^wallb:|^wallq:|^rosasR:|^sweepK:" $edited_stats | sort > edited_sorted.tmp
    
    # Compare each statistic
    local all_match=true
    
    while IFS=: read -r stat_name legacy_values; do
        # Get corresponding values from edited version
        edited_line=$(grep "^$stat_name:" edited_sorted.tmp || echo "")
        
        if [ -z "$edited_line" ]; then
            echo "ERROR: Statistic '$stat_name' missing in edited version" >> $output_file
            all_match=false
            continue
        fi
        
        edited_values=$(echo "$edited_line" | cut -d: -f2-)
        
        # Clean up whitespace
        legacy_values=$(echo "$legacy_values" | xargs)
        edited_values=$(echo "$edited_values" | xargs)
        
        # Compare values
        if [ "$legacy_values" = "$edited_values" ]; then
            echo "✓ $stat_name: MATCH" >> $output_file
        else
            echo "✗ $stat_name: DIFFER" >> $output_file
            echo "  Legacy: $legacy_values" >> $output_file
            echo "  Edited: $edited_values" >> $output_file
            all_match=false
            
            # Try to compute relative difference for numeric values
            # This is a simple comparison - could be enhanced with proper statistical tests
            if [[ "$legacy_values" =~ ^[0-9.]+$ ]] && [[ "$edited_values" =~ ^[0-9.]+$ ]]; then
                rel_diff=$(awk "BEGIN {if ($legacy_values != 0) print abs(($edited_values - $legacy_values) / $legacy_values * 100); else print 0}")
                echo "  Relative difference: ${rel_diff}%" >> $output_file
            fi
        fi
    done < legacy_sorted.tmp
    
    # Clean up temp files
    rm -f legacy_sorted.tmp edited_sorted.tmp
    
    echo >> $output_file
    
    if [ "$all_match" = true ]; then
        echo "Result: ALL STATISTICS MATCH ✓" >> $output_file
        return 0
    else
        echo "Result: STATISTICS DIFFER ✗" >> $output_file
        return 1
    fi
}

# Test suite
echo
echo "Running demographic event tests..."
echo

# Summary file
SUMMARY_FILE="summary.txt"
echo "Demographic Events Validation Summary" > $SUMMARY_FILE
echo "====================================" >> $SUMMARY_FILE
echo "Date: $(date)" >> $SUMMARY_FILE
echo "Replicates: $NREPS" >> $SUMMARY_FILE
echo >> $SUMMARY_FILE

# Track overall success
all_tests_passed=true

# Test 1: Basic admixture event
echo -e "\n${YELLOW}Test 1: Basic admixture event${NC}"
echo "Parameters: 6 $NREPS 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2 -ea 0.02 0 1 0 0.15 -d 1300526888 1325944824"

if run_and_analyze "discoal_legacy_backup" "test1_legacy" 6 $NREPS 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2 -ea 0.02 0 1 0 0.15 -d 1300526888 1325944824 &&
   run_and_analyze "discoal_edited" "test1_edited" 6 $NREPS 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2 -ea 0.02 0 1 0 0.15 -d 1300526888 1325944824; then
    compare_stats "Test 1: Admixture Event" "test1_legacy.stats" "test1_edited.stats" "$SUMMARY_FILE"
    if [ $? -ne 0 ]; then all_tests_passed=false; fi
else
    echo "Test 1: FAILED - Simulation error" >> $SUMMARY_FILE
    all_tests_passed=false
fi

# Test 2: Ancient samples
echo -e "\n${YELLOW}Test 2: Ancient samples${NC}"
echo "Parameters: 10 $NREPS 100 -t 2 -r 2 -p 2 5 5 -M 1.0 -A 2 0 0.1 -A 2 1 0.2 -d 123 456"

if run_and_analyze "discoal_legacy_backup" "test2_legacy" 10 $NREPS 100 -t 2 -r 2 -p 2 5 5 -M 1.0 -A 2 0 0.1 -A 2 1 0.2 -d 123 456 &&
   run_and_analyze "discoal_edited" "test2_edited" 10 $NREPS 100 -t 2 -r 2 -p 2 5 5 -M 1.0 -A 2 0 0.1 -A 2 1 0.2 -d 123 456; then
    compare_stats "Test 2: Ancient Samples" "test2_legacy.stats" "test2_edited.stats" "$SUMMARY_FILE"
    if [ $? -ne 0 ]; then all_tests_passed=false; fi
else
    echo "Test 2: FAILED - Simulation error" >> $SUMMARY_FILE
    all_tests_passed=false
fi

# Test 3: Population merge
echo -e "\n${YELLOW}Test 3: Population merge${NC}"
echo "Parameters: 8 $NREPS 100 -t 3 -r 1.5 -p 2 4 4 -ed 0.5 0 1 -d 999 888"

if run_and_analyze "discoal_legacy_backup" "test3_legacy" 8 $NREPS 100 -t 3 -r 1.5 -p 2 4 4 -ed 0.5 0 1 -d 999 888 &&
   run_and_analyze "discoal_edited" "test3_edited" 8 $NREPS 100 -t 3 -r 1.5 -p 2 4 4 -ed 0.5 0 1 -d 999 888; then
    compare_stats "Test 3: Population Merge" "test3_legacy.stats" "test3_edited.stats" "$SUMMARY_FILE"
    if [ $? -ne 0 ]; then all_tests_passed=false; fi
else
    echo "Test 3: FAILED - Simulation error" >> $SUMMARY_FILE
    all_tests_passed=false
fi

# Test 4: Complex demographic scenario
echo -e "\n${YELLOW}Test 4: Complex demographic scenario${NC}"
echo "Parameters: 12 $NREPS 100 -t 4 -r 3 -p 4 3 3 3 3 -M 0.5 -en 0.1 0 2.0 -en 0.2 1 0.5 -ed 0.3 2 0 -ed 0.4 3 1 -ea 0.15 1 0 2 0.3 -d 111 222"

if run_and_analyze "discoal_legacy_backup" "test4_legacy" 12 $NREPS 100 -t 4 -r 3 -p 4 3 3 3 3 -M 0.5 -en 0.1 0 2.0 -en 0.2 1 0.5 -ed 0.3 2 0 -ed 0.4 3 1 -ea 0.15 1 0 2 0.3 -d 111 222 &&
   run_and_analyze "discoal_edited" "test4_edited" 12 $NREPS 100 -t 4 -r 3 -p 4 3 3 3 3 -M 0.5 -en 0.1 0 2.0 -en 0.2 1 0.5 -ed 0.3 2 0 -ed 0.4 3 1 -ea 0.15 1 0 2 0.3 -d 111 222; then
    compare_stats "Test 4: Complex Demographics" "test4_legacy.stats" "test4_edited.stats" "$SUMMARY_FILE"
    if [ $? -ne 0 ]; then all_tests_passed=false; fi
else
    echo "Test 4: FAILED - Simulation error" >> $SUMMARY_FILE
    all_tests_passed=false
fi

# Test 5: Migration with admixture
echo -e "\n${YELLOW}Test 5: Migration with admixture${NC}"
echo "Parameters: 8 $NREPS 100 -t 2.5 -r 2 -p 3 3 3 2 -m 0 1 2.0 -m 1 0 2.0 -m 1 2 1.0 -m 2 1 1.0 -ea 0.1 0 1 2 0.4 -d 777 666"

if run_and_analyze "discoal_legacy_backup" "test5_legacy" 8 $NREPS 100 -t 2.5 -r 2 -p 3 3 3 2 -m 0 1 2.0 -m 1 0 2.0 -m 1 2 1.0 -m 2 1 1.0 -ea 0.1 0 1 2 0.4 -d 777 666 &&
   run_and_analyze "discoal_edited" "test5_edited" 8 $NREPS 100 -t 2.5 -r 2 -p 3 3 3 2 -m 0 1 2.0 -m 1 0 2.0 -m 1 2 1.0 -m 2 1 1.0 -ea 0.1 0 1 2 0.4 -d 777 666; then
    compare_stats "Test 5: Migration + Admixture" "test5_legacy.stats" "test5_edited.stats" "$SUMMARY_FILE"
    if [ $? -ne 0 ]; then all_tests_passed=false; fi
else
    echo "Test 5: FAILED - Simulation error" >> $SUMMARY_FILE
    all_tests_passed=false
fi

# Final summary
echo >> $SUMMARY_FILE
echo "====================================" >> $SUMMARY_FILE
if [ "$all_tests_passed" = true ]; then
    echo "OVERALL RESULT: ALL TESTS PASSED ✓" >> $SUMMARY_FILE
    echo -e "\n${GREEN}All demographic event tests passed!${NC}"
else
    echo "OVERALL RESULT: SOME TESTS FAILED ✗" >> $SUMMARY_FILE
    echo -e "\n${RED}Some demographic event tests failed!${NC}"
fi

# Display summary
echo
echo "====================================="
cat $SUMMARY_FILE
echo "====================================="

# Save detailed results
echo -e "\nDetailed results saved in: $WORK_DIR"
echo "Summary file: $WORK_DIR/summary.txt"

# Return appropriate exit code
if [ "$all_tests_passed" = true ]; then
    exit 0
else
    exit 1
fi