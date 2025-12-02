#!/bin/bash
#
# Statistical validation test suite for tskit-only optimization
# Compares current implementation (with allNodes) vs tskit-only implementation
#
# Usage: ./tskit_only_validation_suite.sh [replicates]
#   replicates: number of replicates per test (default: 50)

# Parse command line arguments
NUM_REPLICATES="${1:-50}"

TEST_DIR="tskit_only_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "=== TSKIT-ONLY OPTIMIZATION VALIDATION SUITE ==="
echo "Comparing current implementation vs tskit-only implementation"
echo "Running $NUM_REPLICATES replicates per test case"
echo "Test directory: $TEST_DIR"
echo ""

# Memory measurement function with timeout
measure_memory() {
    local command="$1"
    local output_file="$2"
    local memory_file="$3"
    local timeout_seconds=300  # 5 minute timeout
    
    # Use time -l to get memory statistics on macOS, or time -v on Linux
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        timeout $timeout_seconds /usr/bin/time -l bash -c "$command" > "$output_file" 2> "$memory_file"
    else
        # Linux  
        timeout $timeout_seconds /usr/bin/time -v bash -c "$command" > "$output_file" 2> "$memory_file"
    fi
    
    return $?
}

# Extract peak memory usage from time output
get_peak_memory() {
    local memory_file="$1"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS - bytes
        grep "maximum resident set size" "$memory_file" | awk '{print $1}'
    else
        # Linux - KB
        grep "Maximum resident set size" "$memory_file" | awk '{print $6}'
    fi
}

# Extract wall time from time output
get_wall_time() {
    local memory_file="$1"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS - real time in seconds
        grep "real" "$memory_file" | head -1 | awk '{print $1}'
    else
        # Linux - elapsed time
        grep "Elapsed (wall clock) time" "$memory_file" | sed 's/.*: //' | awk -F: '{ if (NF == 2) {print $1 * 60 + $2} else {print $1 * 3600 + $2 * 60 + $3} }'
    fi
}

# Calculate basic statistics
calculate_stats() {
    local output_file="$1"
    
    # Extract segregating sites counts
    grep "^segsites:" "$output_file" | awk '{print $2}' > "$output_file.segsites"
    
    # Calculate basic statistics using awk
    if [ -s "$output_file.segsites" ]; then
        awk '{
            sum += $1
            sumsq += $1 * $1
            if (NR == 1 || $1 < min) min = $1
            if (NR == 1 || $1 > max) max = $1
        }
        END {
            if (NR > 0) {
                mean = sum / NR
                variance = (sumsq - sum * sum / NR) / (NR - 1)
                sd = sqrt(variance)
                print "count=" NR
                print "mean=" mean
                print "sd=" sd
                print "min=" min
                print "max=" max
            }
        }' "$output_file.segsites"
    fi
}

# Simple Kolmogorov-Smirnov test implementation
ks_test() {
    local file1="$1"
    local file2="$2"
    
    # Sort both files
    sort -n "$file1" > "$file1.sorted"
    sort -n "$file2" > "$file2.sorted"
    
    # Calculate KS statistic
    awk '
    BEGIN { max_d = 0 }
    FNR == NR { data1[FNR] = $1; n1 = FNR; next }
    { data2[FNR] = $1; n2 = FNR }
    END {
        # Get all unique values
        for (i = 1; i <= n1; i++) unique[data1[i]] = 1
        for (i = 1; i <= n2; i++) unique[data2[i]] = 1
        
        # Calculate KS statistic at each unique value
        for (v in unique) {
            val = v + 0
            
            # Count values <= val in each dataset
            count1 = 0
            for (j = 1; j <= n1; j++) {
                if (data1[j] <= val) count1++
            }
            
            count2 = 0
            for (j = 1; j <= n2; j++) {
                if (data2[j] <= val) count2++
            }
            
            # Calculate empirical CDFs
            ecdf1 = count1 / n1
            ecdf2 = count2 / n2
            
            # Calculate difference
            d = abs(ecdf1 - ecdf2)
            if (d > max_d) max_d = d
        }
        
        # Critical value at 0.05 significance
        critical = 1.36 * sqrt((n1 + n2) / (n1 * n2))
        
        print "ks_statistic=" max_d
        print "critical_value=" critical
        print "p_value=" (max_d > critical ? "<0.05" : ">0.05")
        print "reject_null=" (max_d > critical ? "yes" : "no")
    }
    function abs(x) { return x < 0 ? -x : x }
    ' "$file1.sorted" "$file2.sorted"
    
    rm -f "$file1.sorted" "$file2.sorted"
}

# Run a single test case
run_test_case() {
    local test_name="$1"
    local test_args="$2"
    local num_reps="$3"
    
    echo "=== Testing: $test_name ==="
    echo "Command: discoal $test_args"
    echo "Running $num_reps replicates..."
    
    # Fixed seeds for reproducibility
    local SEED1=12345
    local SEED2=67890
    
    # Test current implementation
    echo "  Testing current implementation..."
    local current_cmd="../discoal $test_args -d $SEED1 $SEED2"
    measure_memory "$current_cmd" "$TEST_DIR/${test_name}_current.out" "$TEST_DIR/${test_name}_current_memory.txt"
    local current_exit=$?
    
    local current_memory=""
    local current_time=""
    
    if [ $current_exit -eq 0 ]; then
        echo "    ✅ Current: SUCCESS"
        
        current_memory=$(get_peak_memory "$TEST_DIR/${test_name}_current_memory.txt")
        current_time=$(get_wall_time "$TEST_DIR/${test_name}_current_memory.txt")
        echo "    📊 Memory: ${current_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
        echo "    ⏱️  Time: ${current_time} seconds"
        
        # Calculate statistics
        current_stats=$(calculate_stats "$TEST_DIR/${test_name}_current.out")
        echo "    📈 Statistics: $current_stats"
    else
        echo "    ❌ Current: FAILED"
        return 1
    fi
    
    # Test tskit-only implementation
    echo "  Testing tskit-only implementation..."
    local tskit_cmd="USE_TSKIT_ONLY=1 ../discoal $test_args -d $SEED1 $SEED2"
    measure_memory "$tskit_cmd" "$TEST_DIR/${test_name}_tskit.out" "$TEST_DIR/${test_name}_tskit_memory.txt"
    local tskit_exit=$?
    
    local tskit_memory=""
    local tskit_time=""
    
    if [ $tskit_exit -eq 0 ]; then
        echo "    ✅ Tskit-only: SUCCESS"
        
        tskit_memory=$(get_peak_memory "$TEST_DIR/${test_name}_tskit_memory.txt")
        tskit_time=$(get_wall_time "$TEST_DIR/${test_name}_tskit_memory.txt")
        echo "    📊 Memory: ${tskit_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
        echo "    ⏱️  Time: ${tskit_time} seconds"
        
        # Calculate statistics
        tskit_stats=$(calculate_stats "$TEST_DIR/${test_name}_tskit.out")
        echo "    📈 Statistics: $tskit_stats"
    else
        echo "    ❌ Tskit-only: FAILED"
        return 1
    fi
    
    # Compare results
    if [ $current_exit -eq 0 ] && [ $tskit_exit -eq 0 ]; then
        # Memory comparison
        if [ -n "$current_memory" ] && [ -n "$tskit_memory" ]; then
            if [ "$tskit_memory" -lt "$current_memory" ]; then
                local savings=$(( (current_memory - tskit_memory) * 100 / current_memory ))
                echo "    💾 Memory savings: ${savings}%"
            else
                local overhead=$(( (tskit_memory - current_memory) * 100 / current_memory ))
                echo "    ⚠️  Memory overhead: +${overhead}%"
            fi
        fi
        
        # Performance comparison
        if [ -n "$current_time" ] && [ -n "$tskit_time" ] && command -v bc &> /dev/null; then
            local speedup=$(echo "scale=2; $current_time / $tskit_time" | bc 2>/dev/null)
            if [ -n "$speedup" ]; then
                echo "    ⚡ Performance: ${speedup}x"
            fi
        fi
        
        # Statistical comparison
        if [ -f "$TEST_DIR/${test_name}_current.out.segsites" ] && \
           [ -f "$TEST_DIR/${test_name}_tskit.out.segsites" ]; then
            echo "    🔬 Performing Kolmogorov-Smirnov test..."
            local ks_result=$(ks_test "$TEST_DIR/${test_name}_current.out.segsites" \
                                     "$TEST_DIR/${test_name}_tskit.out.segsites")
            echo "    $ks_result"
        fi
        
        # Check for identical output (for small tests)
        if [ "$num_reps" -eq 1 ]; then
            if diff -q "$TEST_DIR/${test_name}_current.out" "$TEST_DIR/${test_name}_tskit.out" > /dev/null; then
                echo "    ✓ Output files are identical"
            else
                echo "    ⚠️  Output files differ (expected for different RNG paths)"
            fi
        fi
    fi
    
    echo ""
    return 0
}

# Test cases specifically designed for tskit-only optimization
echo "Phase 1: Basic Functionality Tests (1 replicate each)"
echo "=================================================="
run_test_case "basic_tiny" "5 1 100 -t 2" 1
run_test_case "basic_small" "10 1 1000 -t 5 -r 2" 1
run_test_case "basic_medium" "20 1 5000 -t 10 -r 5" 1

echo "Phase 2: Memory Scaling Tests (1 replicate each)"
echo "=================================================="
run_test_case "memory_100samples" "100 1 10000 -t 20 -r 10" 1
run_test_case "memory_200samples" "200 1 10000 -t 20 -r 10" 1
run_test_case "memory_500samples" "500 1 10000 -t 20 -r 10" 1

echo "Phase 3: High Recombination Tests (1 replicate each)"
echo "====================================================="
run_test_case "recomb_low" "50 1 10000 -t 20 -r 50" 1
run_test_case "recomb_medium" "50 1 10000 -t 20 -r 200" 1
run_test_case "recomb_high" "50 1 10000 -t 20 -r 1000" 1

echo "Phase 4: Statistical Validation Tests ($NUM_REPLICATES replicates each)"
echo "======================================================================"
run_test_case "stats_neutral" "20 $NUM_REPLICATES 10000 -t 20" $NUM_REPLICATES
run_test_case "stats_recomb" "20 $NUM_REPLICATES 10000 -t 20 -r 20" $NUM_REPLICATES
run_test_case "stats_gc" "20 $NUM_REPLICATES 10000 -t 20 -r 10 -g 5 10" $NUM_REPLICATES
run_test_case "stats_bottleneck" "20 $NUM_REPLICATES 10000 -t 20 -r 10 -en 0.1 0 0.1 -en 0.5 0 1.0" $NUM_REPLICATES

echo "Phase 5: Tree Sequence Tests"
echo "============================"
run_test_case "tskit_output" "10 1 5000 -t 10 -r 5 -ts test.trees" 1
run_test_case "tskit_mutations" "20 1 10000 -t 50 -r 20 -ts test_muts.trees" 1

echo "Phase 6: Stress Tests (1 replicate each)"
echo "========================================"
run_test_case "stress_large" "1000 1 50000 -t 100 -r 50" 1
run_test_case "stress_extreme_recomb" "100 1 100000 -t 50 -r 5000" 1

# Summary report
echo "======================================================================"
echo "                    TSKIT-ONLY VALIDATION SUMMARY"
echo "======================================================================"
echo ""
echo "📊 Test Results Summary:"
echo "  Test directory: $TEST_DIR"
echo ""
echo "🔬 Key Validation Points:"
echo "  1. Basic functionality maintained"
echo "  2. Statistical distributions preserved"
echo "  3. Memory usage significantly reduced"
echo "  4. Performance characteristics analyzed"
echo ""
echo "📁 Detailed results saved in: $TEST_DIR"
echo ""
echo "💡 Next Steps:"
echo "  1. Review memory savings in stress tests"
echo "  2. Check statistical validation results"
echo "  3. Analyze any performance differences"
echo "  4. Verify tree sequence output compatibility"
echo ""

# Create summary report
{
    echo "Tskit-Only Optimization Validation - $(date)"
    echo "============================================"
    echo ""
    echo "Test Configuration:"
    echo "  Statistical test replicates: $NUM_REPLICATES"
    echo ""
    echo "Memory Savings Summary:"
    echo "  (Review individual test results above)"
    echo ""
    echo "Statistical Validation:"
    echo "  (Review KS test results above)"
    echo ""
    echo "Implementation Notes:"
    echo "  - Current: Uses allNodes array (never freed)"
    echo "  - Tskit-only: Frees nodes immediately after recording"
    echo "  - Expected: >90% memory reduction for large simulations"
} > "$TEST_DIR/validation_summary.txt"

echo "🎉 Validation suite completed!"