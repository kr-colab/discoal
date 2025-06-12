#!/bin/bash

# msprime comparison test suite for discoal
# This suite compares discoal output with msprime simulations to ensure
# statistical equivalence for neutral and selection models

# Ensure we're in the testing directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check for required tools
if ! command -v /home/adkern/miniforge3/envs/discoal_dev/bin/python &> /dev/null; then
    echo "ERROR: /home/adkern/miniforge3/envs/discoal_dev/bin/python is required but not found"
    exit 1
fi

# Check if msprime is available
if ! /home/adkern/miniforge3/envs/discoal_dev/bin/python -c "import msprime" 2>/dev/null; then
    echo "ERROR: msprime Python package is required but not found"
    echo "Please activate the discoal_dev conda environment: conda activate discoal_dev"
    exit 1
fi

# Check if scipy is available
if ! /home/adkern/miniforge3/envs/discoal_dev/bin/python -c "import scipy" 2>/dev/null; then
    echo "ERROR: scipy Python package is required but not found"
    echo "Please install scipy: pip install scipy"
    exit 1
fi

# Create test directory
TEST_DIR="msprime_comparison_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "=== MSPRIME COMPARISON TEST SUITE FOR DISCOAL ==="
echo "Test directory: $TEST_DIR"
echo ""

# Build discoal if needed
cd ..
if [ ! -f "discoal" ]; then
    echo "Building discoal..."
    make discoal
fi
DISCOAL_BIN="$PWD/discoal"
cd "$SCRIPT_DIR"

# Test configuration
N_REPLICATES=100
SIGNIFICANCE_LEVEL=0.05
NE=1  # Use Ne=1 for msprime with haploid samples
REFSIZE=10000  # Reference population size for sweep scaling

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to extract segregating sites from discoal output
extract_seg_sites() {
    grep "segsites:" | awk '{print $2}'
}

# Function to run comparison test
run_comparison() {
    local test_name="$1"
    local discoal_cmd="$2"
    local msprime_args="$3"
    local description="$4"
    
    echo ""
    echo "Running test: $test_name"
    echo "Description: $description"
    echo "Discoal command: $discoal_cmd"
    
    # Run discoal simulations with timing
    echo -n "Running discoal simulations..."
    local discoal_start=$(date +%s.%N)
    eval "$discoal_cmd" 2>/dev/null | extract_seg_sites > "$TEST_DIR/${test_name}_discoal.txt"
    local discoal_end=$(date +%s.%N)
    local discoal_time=$(echo "$discoal_end - $discoal_start" | bc)
    local discoal_count=$(wc -l < "$TEST_DIR/${test_name}_discoal.txt")
    echo " done ($discoal_count replicates, ${discoal_time}s)"
    
    # Run msprime simulations with timing
    echo -n "Running msprime simulations..."
    local msprime_start=$(date +%s.%N)
    /home/adkern/miniforge3/envs/discoal_dev/bin/python msprime_comparison.py $msprime_args > "$TEST_DIR/${test_name}_msprime.txt"
    local msprime_end=$(date +%s.%N)
    local msprime_time=$(echo "$msprime_end - $msprime_start" | bc)
    local msprime_count=$(wc -l < "$TEST_DIR/${test_name}_msprime.txt")
    echo " done ($msprime_count replicates, ${msprime_time}s)"
    
    # Statistical comparison
    echo "Performing statistical comparison..."
    /home/adkern/miniforge3/envs/discoal_dev/bin/python - <<EOF > "$TEST_DIR/${test_name}_comparison.txt"
import numpy as np
from scipy import stats

# Read data
discoal_data = np.loadtxt("$TEST_DIR/${test_name}_discoal.txt")
msprime_data = np.loadtxt("$TEST_DIR/${test_name}_msprime.txt")

# Calculate statistics
d_mean = np.mean(discoal_data)
d_std = np.std(discoal_data)
m_mean = np.mean(msprime_data)
m_std = np.std(msprime_data)

# Statistical tests
ks_stat, ks_pval = stats.ks_2samp(discoal_data, msprime_data)
mw_stat, mw_pval = stats.mannwhitneyu(discoal_data, msprime_data, alternative='two-sided')

# Output results
print(f"Discoal: mean={d_mean:.2f}, std={d_std:.2f}")
print(f"msprime: mean={m_mean:.2f}, std={m_std:.2f}")
print(f"Relative difference: {abs(d_mean - m_mean) / max(d_mean, m_mean) * 100:.1f}%")
print(f"KS test: statistic={ks_stat:.4f}, p-value={ks_pval:.4f}")
print(f"MW test: statistic={mw_stat:.1f}, p-value={mw_pval:.4f}")

# Save detailed results
with open("$TEST_DIR/${test_name}_results.txt", "w") as f:
    f.write(f"TEST: ${test_name}\\n")
    f.write(f"DISCOAL_MEAN: {d_mean}\\n")
    f.write(f"DISCOAL_STD: {d_std}\\n")
    f.write(f"MSPRIME_MEAN: {m_mean}\\n")
    f.write(f"MSPRIME_STD: {m_std}\\n")
    f.write(f"KS_PVALUE: {ks_pval}\\n")
    f.write(f"MW_PVALUE: {mw_pval}\\n")
    f.write(f"SIMILAR: {'YES' if ks_pval > $SIGNIFICANCE_LEVEL else 'NO'}\\n")
    f.write(f"DISCOAL_TIME: ${discoal_time}\\n")
    f.write(f"MSPRIME_TIME: ${msprime_time}\\n")
EOF
    
    # Display results
    cat "$TEST_DIR/${test_name}_comparison.txt"
    
    # Performance comparison
    local speedup=$(echo "scale=2; $msprime_time / $discoal_time" | bc)
    echo "Performance: discoal=${discoal_time}s, msprime=${msprime_time}s (discoal is ${speedup}x faster)"
    
    # Save timing data
    echo "${test_name},${discoal_time},${msprime_time},${speedup}" >> "$TEST_DIR/timing_summary.csv"
    
    # Check if distributions are similar
    local ks_pval=$(grep "KS_PVALUE:" "$TEST_DIR/${test_name}_results.txt" | awk '{print $2}')
    if (( $(echo "$ks_pval > $SIGNIFICANCE_LEVEL" | bc -l) )); then
        echo -e "${GREEN}✓ PASS${NC}: Distributions are statistically similar (p > $SIGNIFICANCE_LEVEL)"
        echo "PASS" >> "$TEST_DIR/test_results.txt"
    else
        echo -e "${RED}✗ FAIL${NC}: Distributions differ significantly (p = $ks_pval)"
        echo "FAIL: $test_name (p = $ks_pval)" >> "$TEST_DIR/test_results.txt"
    fi
}

# Initialize timing summary file
echo "test_name,discoal_time,msprime_time,speedup" > "$TEST_DIR/timing_summary.csv"

# Test cases
echo "=== TEST CASES ==="
echo "Note: Using Ne=$NE for parameter scaling (discoal's scaled parameters assume this)"
echo ""

# 1. Basic neutral model
run_comparison \
    "neutral_basic" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10" \
    "--mode neutral --samples 20 --sites 10000 --theta 10 --Ne $NE --replicates $N_REPLICATES" \
    "Basic neutral coalescent without recombination"

# 2. Neutral with recombination
run_comparison \
    "neutral_recomb" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20" \
    "--mode recombination --samples 20 --sites 10000 --theta 10 --rho 20 --Ne $NE --replicates $N_REPLICATES" \
    "Neutral coalescent with recombination"

# 3. Higher recombination rate
run_comparison \
    "neutral_high_recomb" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 20 -r 100" \
    "--mode recombination --samples 20 --sites 10000 --theta 20 --rho 100 --Ne $NE --replicates $N_REPLICATES" \
    "Neutral coalescent with high recombination"

# 4. Large sample size
run_comparison \
    "neutral_large_sample" \
    "$DISCOAL_BIN 50 $N_REPLICATES 10000 -t 20" \
    "--mode neutral --samples 50 --sites 10000 --theta 20 --Ne $NE --replicates $N_REPLICATES" \
    "Neutral coalescent with larger sample size"

# 5. Small sample size
run_comparison \
    "neutral_small_sample" \
    "$DISCOAL_BIN 10 $N_REPLICATES 10000 -t 5" \
    "--mode neutral --samples 10 --sites 10000 --theta 5 --Ne $NE --replicates $N_REPLICATES" \
    "Neutral coalescent with small sample size"

# 6. Hard sweep (stochastic)
run_comparison \
    "sweep_stochastic" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20 -ws 0.0 -a 1000 -x 0.5 -N $REFSIZE" \
    "--mode sweep --samples 20 --sites 10000 --theta 10 --rho 20 --alpha 1000 --tau 0.0 --Ne $NE --refsize $REFSIZE --replicates $N_REPLICATES" \
    "Stochastic hard sweep (alpha=1000, tau=0.0)"

# 7. Strong sweep
run_comparison \
    "sweep_strong" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20 -ws 0.00 -a 2000 -x 0.5 -N $REFSIZE" \
    "--mode sweep --samples 20 --sites 10000 --theta 10 --rho 20 --alpha 2000 --tau 0.0 --Ne $NE --refsize $REFSIZE --replicates $N_REPLICATES" \
    "Strong hard sweep (alpha=2000, tau=0.0)"

# 8. Weak sweep  
run_comparison \
    "sweep_weak" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20 -ws 0.0 -a 500 -x 0.5 -N $REFSIZE" \
    "--mode sweep --samples 20 --sites 10000 --theta 10 --rho 20 --alpha 500 --tau 0.0 --Ne $NE --refsize $REFSIZE --replicates $N_REPLICATES" \
    "Weak hard sweep (alpha=500, tau=0.0)"

# 9. Recent sweep
run_comparison \
    "sweep_recent" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20 -ws 0.00 -a 1000 -x 0.5 -N $REFSIZE" \
    "--mode sweep --samples 20 --sites 10000 --theta 10 --rho 20 --alpha 1000 --tau 0.0 --Ne $NE --refsize $REFSIZE --replicates $N_REPLICATES" \
    "Recent hard sweep (alpha=1000, tau=0.0)"

# 10. Old sweep
run_comparison \
    "sweep_old" \
    "$DISCOAL_BIN 20 $N_REPLICATES 10000 -t 10 -r 20 -ws 0.05 -a 1000 -x 0.5 -N $REFSIZE" \
    "--mode sweep --samples 20 --sites 10000 --theta 10 --rho 20 --alpha 1000 --tau 0.05 --Ne $NE --refsize $REFSIZE --replicates $N_REPLICATES" \
    "Old hard sweep (alpha=1000, tau=0.05)"

echo ""
echo "=== SUMMARY ==="

# Count passes and fails
TOTAL_TESTS=$(grep -c "." "$TEST_DIR/test_results.txt" 2>/dev/null || echo 0)
PASSED_TESTS=$(grep -c "PASS" "$TEST_DIR/test_results.txt" 2>/dev/null || echo 0)
FAILED_TESTS=$(grep -c "FAIL" "$TEST_DIR/test_results.txt" 2>/dev/null || echo 0)

echo "Total tests: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed: ${RED}$FAILED_TESTS${NC}"

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "\n${GREEN}All tests passed!${NC} Discoal and msprime produce statistically equivalent results."
else
    echo -e "\n${YELLOW}Some tests failed.${NC} Check $TEST_DIR for detailed results."
    echo "Failed tests:"
    grep "FAIL:" "$TEST_DIR/test_results.txt"
fi

echo ""
echo "=== PERFORMANCE SUMMARY ==="
echo ""
echo "Runtime comparison (100 replicates each):"
echo ""
printf "%-20s %15s %15s %10s\n" "Test" "Discoal (s)" "msprime (s)" "Speedup"
printf "%-20s %15s %15s %10s\n" "----" "-----------" "-----------" "-------"
tail -n +2 "$TEST_DIR/timing_summary.csv" | while IFS=',' read -r test discoal msprime speedup; do
    printf "%-20s %15s %15s %10sx\n" "$test" "$discoal" "$msprime" "$speedup"
done

echo ""
# Calculate average speedup
avg_speedup=$(tail -n +2 "$TEST_DIR/timing_summary.csv" | awk -F',' '{sum+=$4; count++} END {printf "%.2f", sum/count}')
echo "Average speedup: discoal is ${avg_speedup}x faster than msprime"

echo ""
echo "Detailed results saved in: $TEST_DIR/"
echo ""
echo "To examine specific comparisons:"
echo "  - Segregating sites distributions: *_discoal.txt and *_msprime.txt"
echo "  - Statistical test results: *_results.txt"
echo "  - Comparison summaries: *_comparison.txt"
echo "  - Performance data: timing_summary.csv"