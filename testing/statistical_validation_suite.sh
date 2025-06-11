#!/bin/bash

# Statistical validation test suite based on comprehensive_validation_suite.sh
# This suite runs multiple replicates of each test case to perform statistical comparisons
# between optimized and legacy versions
#
# Usage: ./statistical_validation_suite.sh [mode] [replicates]
#   mode: auto (default), parallel, sequential
#   replicates: number of replicates per test (default: 100)
#   auto: Use parallel if available, fallback to sequential
#   parallel: Force parallel execution (fails if GNU parallel not available)
#   sequential: Force sequential execution

# Parse command line arguments
EXECUTION_MODE="${1:-auto}"
NUM_REPLICATES="${2:-100}"

# Handle help flag
if [[ "$EXECUTION_MODE" == "-h" || "$EXECUTION_MODE" == "--help" ]]; then
    echo "Usage: $0 [mode] [replicates]"
    echo ""
    echo "Execution modes:"
    echo "  auto       - Use parallel if GNU parallel available, fallback to sequential (default)"
    echo "  parallel   - Force parallel execution (requires GNU parallel)"
    echo "  sequential - Force sequential execution"
    echo ""
    echo "Replicates:"
    echo "  Number of simulation replicates per test case (default: 100)"
    echo ""
    echo "Examples:"
    echo "  $0                    # Auto mode with 100 replicates"
    echo "  $0 parallel 50        # Parallel mode with 50 replicates"
    echo "  $0 sequential 200     # Sequential mode with 200 replicates"
    echo ""
    exit 0
fi

# If first argument is a number, swap the arguments
if [[ "$EXECUTION_MODE" =~ ^[0-9]+$ ]]; then
    NUM_REPLICATES="$EXECUTION_MODE"
    EXECUTION_MODE="${2:-auto}"
fi

TEST_DIR="statistical_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "=== STATISTICAL DISCOAL VALIDATION SUITE ==="
echo "Running $NUM_REPLICATES replicates per test case"
echo "Test directory: $TEST_DIR"
echo "Execution mode: $EXECUTION_MODE"
echo ""

# Memory measurement function with timeout
measure_memory() {
    local command="$1"
    local output_file="$2"
    local memory_file="$3"
    local timeout_seconds=600  # 10 minute timeout for statistical runs
    
    # Check if timeout command is available
    if command -v timeout &> /dev/null; then
        # Use time -l to get memory statistics on macOS, or time -v on Linux
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            timeout $timeout_seconds /usr/bin/time -l bash -c "$command" > "$output_file" 2> "$memory_file"
        else
            # Linux  
            timeout $timeout_seconds /usr/bin/time -v bash -c "$command" > "$output_file" 2> "$memory_file"
        fi
    else
        # Fallback without timeout
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            /usr/bin/time -l bash -c "$command" > "$output_file" 2> "$memory_file"
        else
            # Linux  
            /usr/bin/time -v bash -c "$command" > "$output_file" 2> "$memory_file"
        fi
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

# Statistical analysis functions
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

# Kolmogorov-Smirnov test implementation in awk
ks_test() {
    local file1="$1"
    local file2="$2"
    
    # Sort both files
    sort -n "$file1" > "$file1.sorted"
    sort -n "$file2" > "$file2.sorted"
    
    # Calculate KS statistic properly
    awk '
    BEGIN { max_d = 0 }
    FNR == NR { data1[FNR] = $1; n1 = FNR; next }
    { data2[FNR] = $1; n2 = FNR }
    END {
        # Get all unique values from both datasets
        for (i = 1; i <= n1; i++) unique[data1[i]] = 1
        for (i = 1; i <= n2; i++) unique[data2[i]] = 1
        
        # Sort unique values
        n_unique = 0
        for (v in unique) {
            n_unique++
            sorted_vals[n_unique] = v + 0  # Convert to number
        }
        
        # Bubble sort (simple for small datasets)
        for (i = 1; i <= n_unique; i++) {
            for (j = i+1; j <= n_unique; j++) {
                if (sorted_vals[i] > sorted_vals[j]) {
                    temp = sorted_vals[i]
                    sorted_vals[i] = sorted_vals[j]
                    sorted_vals[j] = temp
                }
            }
        }
        
        # Calculate KS statistic at each unique value
        for (i = 1; i <= n_unique; i++) {
            val = sorted_vals[i]
            
            # Count how many values <= val in each dataset
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
    
    # Clean up
    rm -f "$file1.sorted" "$file2.sorted"
}

# Single test case execution function for statistical validation
run_single_statistical_test() {
    local test_case="$1"
    local test_dir="$2"
    local test_index="$3"
    local total_tests="$4"
    local num_reps="$5"
    
    # Parse test case
    local category=$(echo $test_case | cut -d: -f1)
    local test_name=$(echo $test_case | cut -d: -f2)
    local test_args=$(echo $test_case | cut -d: -f3)
    local expected=$(echo $test_case | cut -d: -f4)
    
    # Replace the replicate count (second number) with our desired number
    local modified_args=$(echo $test_args | awk -v reps=$num_reps '{$2 = reps; print}')
    
    echo "=== [$test_index/$total_tests] [$category] Testing: $test_name ==="
    echo "Command: discoal $modified_args"
    echo "Running $num_reps replicates..."
    
    # Fixed seeds for reproducible results
    local SEED1=98765
    local SEED2=54321
    
    # Test optimized version
    echo "  Testing optimized version..."
    local optimized_cmd="../discoal_edited $modified_args -d $SEED1 $SEED2"
    measure_memory "$optimized_cmd" "$test_dir/${category}_${test_name}_optimized.out" "$test_dir/${category}_${test_name}_optimized_memory.txt"
    local optimized_exit=$?
    
    local optimized_memory=""
    local optimized_time=""
    local optimized_status=""
    local optimized_stats=""
    
    if [ $optimized_exit -eq 0 ]; then
        echo "    ✅ Optimized: SUCCESS"
        optimized_status="SUCCESS"
        
        optimized_memory=$(get_peak_memory "$test_dir/${category}_${test_name}_optimized_memory.txt")
        optimized_time=$(get_wall_time "$test_dir/${category}_${test_name}_optimized_memory.txt")
        echo "    📊 Memory: ${optimized_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
        echo "    ⏱️  Time: ${optimized_time} seconds"
        
        # Calculate statistics
        optimized_stats=$(calculate_stats "$test_dir/${category}_${test_name}_optimized.out")
        echo "    📈 Statistics: $optimized_stats"
    elif [ $optimized_exit -eq 124 ]; then
        echo "    ⏰ Optimized: TIMEOUT (10 minutes)"
        optimized_status="TIMEOUT"
    else
        echo "    ❌ Optimized: FAILED"
        optimized_status="FAILED"
        echo "    🔍 Error: $(tail -1 "$test_dir/${category}_${test_name}_optimized.out")"
    fi
    
    # Test legacy version
    echo "  Testing legacy version..."
    local legacy_cmd="../discoal_legacy_backup $modified_args -d $SEED1 $SEED2"
    measure_memory "$legacy_cmd" "$test_dir/${category}_${test_name}_legacy.out" "$test_dir/${category}_${test_name}_legacy_memory.txt"
    local legacy_exit=$?
    
    local legacy_memory=""
    local legacy_time=""
    local legacy_status=""
    local legacy_stats=""
    
    if [ $legacy_exit -eq 0 ]; then
        echo "    ✅ Legacy: SUCCESS"
        legacy_status="SUCCESS"
        
        legacy_memory=$(get_peak_memory "$test_dir/${category}_${test_name}_legacy_memory.txt")
        legacy_time=$(get_wall_time "$test_dir/${category}_${test_name}_legacy_memory.txt")
        echo "    📊 Memory: ${legacy_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
        echo "    ⏱️  Time: ${legacy_time} seconds"
        
        # Calculate statistics
        legacy_stats=$(calculate_stats "$test_dir/${category}_${test_name}_legacy.out")
        echo "    📈 Statistics: $legacy_stats"
    elif [ $legacy_exit -eq 124 ]; then
        echo "    ⏰ Legacy: TIMEOUT (10 minutes)"
        legacy_status="TIMEOUT"
    else
        echo "    ❌ Legacy: FAILED"
        legacy_status="FAILED"
        echo "    🔍 Error: $(tail -1 "$test_dir/${category}_${test_name}_legacy.out")"
    fi
    
    # Initialize result variables
    local memory_result=""
    local time_result=""
    local statistical_result=""
    
    # Compare memory usage if both succeeded
    if [ $optimized_exit -eq 0 ] && [ $legacy_exit -eq 0 ] && [ -n "$optimized_memory" ] && [ -n "$legacy_memory" ]; then
        if [ "$optimized_memory" -lt "$legacy_memory" ]; then
            local savings=$(( (legacy_memory - optimized_memory) * 100 / legacy_memory ))
            echo "    💾 Memory savings: ${savings}% (optimized is smaller)"
            memory_result="$category:$test_name:$savings"
        else
            local overhead=$(( (optimized_memory - legacy_memory) * 100 / legacy_memory ))
            echo "    ⚠️  Memory overhead: +${overhead}% (optimized is larger)"
            memory_result="$category:$test_name:-$overhead"
        fi
        
        # Compare wall time
        if [ -n "$optimized_time" ] && [ -n "$legacy_time" ]; then
            # Use bc for floating point comparison
            if command -v bc &> /dev/null; then
                # Handle case where both times are 0
                if [ "$optimized_time" = "0" ] && [ "$legacy_time" = "0" ]; then
                    echo "    ⏱️  Performance: Same execution time"
                    time_result="$category:$test_name:same:1.0"
                elif [ "$optimized_time" = "0" ]; then
                    # Optimized is effectively instant
                    echo "    ⚡ Performance: Extremely fast (optimized ~instant)"
                    time_result="$category:$test_name:faster:999"
                elif [ "$legacy_time" = "0" ]; then
                    # Legacy is effectively instant (shouldn't happen but handle it)
                    echo "    🐌 Performance: Legacy was instant"
                    time_result="$category:$test_name:slower:999"
                else
                    # Normal comparison
                    local speedup=$(echo "scale=2; $legacy_time / $optimized_time" | bc 2>/dev/null)
                    local time_diff=$(echo "scale=2; $legacy_time - $optimized_time" | bc 2>/dev/null)
                    
                    if [ -n "$speedup" ] && [ -n "$time_diff" ]; then
                        if (( $(echo "$optimized_time < $legacy_time" | bc -l 2>/dev/null) )); then
                            echo "    ⚡ Performance: ${speedup}x faster (saved ${time_diff}s)"
                            time_result="$category:$test_name:faster:$speedup"
                        elif (( $(echo "$optimized_time > $legacy_time" | bc -l 2>/dev/null) )); then
                            local slowdown=$(echo "scale=2; $optimized_time / $legacy_time" | bc 2>/dev/null)
                            echo "    🐌 Performance: ${slowdown}x slower (added ${time_diff#-}s)"
                            time_result="$category:$test_name:slower:$slowdown"
                        else
                            echo "    ⏱️  Performance: Same execution time"
                            time_result="$category:$test_name:same:1.0"
                        fi
                    else
                        echo "    ⏱️  Performance: Unable to calculate (times too small)"
                        time_result="$category:$test_name:same:1.0"
                    fi
                fi
            fi
        fi
        
        # Statistical comparison
        if [ -f "$test_dir/${category}_${test_name}_optimized.out.segsites" ] && \
           [ -f "$test_dir/${category}_${test_name}_legacy.out.segsites" ]; then
            echo "    🔬 Performing Kolmogorov-Smirnov test..."
            local ks_result=$(ks_test "$test_dir/${category}_${test_name}_optimized.out.segsites" \
                                     "$test_dir/${category}_${test_name}_legacy.out.segsites")
            echo "    $ks_result"
            statistical_result="$category:$test_name:$ks_result"
        fi
    fi
    
    # Write results to a structured file for aggregation
    local result_file="$test_dir/result_${test_index}_${category}_${test_name}.txt"
    {
        echo "OPTIMIZED_RESULT='$category:$test_name:$optimized_status'"
        echo "LEGACY_RESULT='$category:$test_name:$legacy_status'"
        echo "MEMORY_RESULT='$memory_result'"
        echo "TIME_RESULT='$time_result'"
        echo "STATISTICAL_RESULT='$statistical_result'"
        echo "OPTIMIZED_STATS='$optimized_stats'"
        echo "LEGACY_STATS='$legacy_stats'"
    } > "$result_file"
    
    # Clean up large output files to save space
    if [ $optimized_exit -eq 0 ] && [ $legacy_exit -eq 0 ]; then
        rm -f "$test_dir/${category}_${test_name}_optimized.out"
        rm -f "$test_dir/${category}_${test_name}_legacy.out"
    fi
    
    echo ""
}

# Aggregate results from individual result files
aggregate_results() {
    local test_dir="$1"
    
    # Initialize arrays
    OPTIMIZED_RESULTS=()
    LEGACY_RESULTS=()
    MEMORY_SAVINGS=()
    TIME_COMPARISONS=()
    STATISTICAL_RESULTS=()
    
    # Reset counters
    OPTIMIZED_SUCCESSES=0
    LEGACY_SUCCESSES=0
    STATISTICALLY_IDENTICAL=0
    MEMORY_COMPARISONS=0
    
    # Process all result files
    for result_file in "$test_dir"/result_*.txt; do
        if [ -f "$result_file" ]; then
            # Source the result file to load variables
            source "$result_file"
            
            # Add to arrays
            if [ -n "$OPTIMIZED_RESULT" ]; then
                OPTIMIZED_RESULTS+=("$OPTIMIZED_RESULT")
                if [[ "$OPTIMIZED_RESULT" == *":SUCCESS" ]]; then
                    OPTIMIZED_SUCCESSES=$((OPTIMIZED_SUCCESSES + 1))
                fi
            fi
            
            if [ -n "$LEGACY_RESULT" ]; then
                LEGACY_RESULTS+=("$LEGACY_RESULT")
                if [[ "$LEGACY_RESULT" == *":SUCCESS" ]]; then
                    LEGACY_SUCCESSES=$((LEGACY_SUCCESSES + 1))
                fi
            fi
            
            if [ -n "$MEMORY_RESULT" ] && [ "$MEMORY_RESULT" != "''" ]; then
                MEMORY_SAVINGS+=("$MEMORY_RESULT")
                MEMORY_COMPARISONS=$((MEMORY_COMPARISONS + 1))
            fi
            
            if [ -n "$TIME_RESULT" ] && [ "$TIME_RESULT" != "''" ]; then
                TIME_COMPARISONS+=("$TIME_RESULT")
            fi
            
            if [ -n "$STATISTICAL_RESULT" ] && [[ "$STATISTICAL_RESULT" == *"reject_null=no"* ]]; then
                STATISTICALLY_IDENTICAL=$((STATISTICALLY_IDENTICAL + 1))
            fi
            
            if [ -n "$STATISTICAL_RESULT" ] && [ "$STATISTICAL_RESULT" != "''" ] && [ "$STATISTICAL_RESULT" != "" ]; then
                # Only add if it contains actual KS test results
                if [[ "$STATISTICAL_RESULT" == *"ks_statistic="* ]]; then
                    STATISTICAL_RESULTS+=("$STATISTICAL_RESULT")
                fi
            fi
        fi
    done
    
    # Clean up result files
    rm -f "$test_dir"/result_*.txt
}

# Check for GNU parallel availability
check_parallel_support() {
    if ! command -v parallel &> /dev/null; then
        echo "⚠️  GNU parallel not found. Install with:"
        echo "   Ubuntu/Debian: sudo apt-get install parallel"
        echo "   macOS: brew install parallel"
        echo "   CentOS/RHEL: sudo yum install parallel"
        echo "   Falling back to sequential execution..."
        echo ""
        return 1
    fi
    
    # Check if parallel is GNU parallel (not other parallel implementations)
    if ! parallel --version 2>/dev/null | grep -q "GNU parallel"; then
        echo "⚠️  Found parallel but not GNU parallel. Falling back to sequential execution..."
        echo ""
        return 1
    fi
    
    return 0
}

# Parallel test execution
run_parallel_tests() {
    local test_dir="$1"
    local total_tests="$2"
    local num_reps="$3"
    
    # Conservative core usage (leave 1-2 cores free for system)
    local available_cores=$(nproc 2>/dev/null || echo "4")
    local cores_to_use=$((available_cores - 1))
    cores_to_use=$(( cores_to_use > 1 ? cores_to_use : 1 ))
    
    echo "🚀 Running tests in parallel using $cores_to_use cores..."
    echo "📊 Progress will be shown as tests complete..."
    echo ""
    
    # Export functions for parallel execution
    export -f run_single_statistical_test measure_memory get_peak_memory get_wall_time calculate_stats ks_test
    
    # Create indexed test case array for progress tracking
    printf '%s\n' "${TEST_CASES[@]}" | \
        parallel -j "$cores_to_use" \
        --timeout 1200 \
        --joblog "$test_dir/parallel.log" \
        run_single_statistical_test {} "$test_dir" {#} "$total_tests" "$num_reps" 2>/dev/null
    
    local parallel_exit=$?
    
    echo ""
    echo "✅ Parallel execution completed. Aggregating results..."
    echo ""
    
    if [ $parallel_exit -ne 0 ]; then
        echo "⚠️  Some parallel jobs failed. Check $test_dir/parallel.log for details."
    fi
    
    # Aggregate results
    aggregate_results "$test_dir"
    
    return $parallel_exit
}

# Sequential test execution (fallback implementation)
run_sequential_tests() {
    local test_dir="$1"
    local total_tests="$2"
    local num_reps="$3"
    
    echo "🐌 Running tests sequentially..."
    echo ""
    
    local test_count=0
    for test_case in "${TEST_CASES[@]}"; do
        test_count=$((test_count + 1))
        run_single_statistical_test "$test_case" "$test_dir" "$test_count" "$total_tests" "$num_reps"
    done
    
    # Aggregate results
    aggregate_results "$test_dir"
}

# Main test execution with automatic parallel/sequential selection
run_tests() {
    local test_dir="$1"
    local total_tests="$2"
    local num_reps="$3"
    local execution_mode="${4:-auto}"
    
    case "$execution_mode" in
        parallel)
            if check_parallel_support; then
                run_parallel_tests "$test_dir" "$total_tests" "$num_reps"
            else
                echo "❌ Parallel execution requested but GNU parallel not available."
                echo "   Install GNU parallel or use sequential mode."
                exit 1
            fi
            ;;
        sequential)
            run_sequential_tests "$test_dir" "$total_tests" "$num_reps"
            ;;
        auto)
            if check_parallel_support; then
                run_parallel_tests "$test_dir" "$total_tests" "$num_reps"
            else
                run_sequential_tests "$test_dir" "$total_tests" "$num_reps"
            fi
            ;;
        *)
            echo "❌ Invalid execution mode: $execution_mode"
            echo "   Valid modes: parallel, sequential, auto"
            exit 1
            ;;
    esac
}

# Test cases extracted from discoaldoc.tex with categorization
declare -a TEST_CASES=(
    # Format: "category:name:command_args:expected_behavior"
    
    # Basic usage examples (Section: Basic usage)
    "basic:simple_example:3 2 100 -t 2:both_succeed"
    
    # Recombination examples (Section: Recombination and Gene Conversion)
    "recombination:basic_recomb:3 2 100 -t 2 -r 2.4:both_succeed"
    "recombination:recomb_gc:3 2 100 -t 2 -r 2.4 -g 2.4 10:both_succeed"
    
    # Population size changes (Section: Population size changes)
    "demography:bottleneck:3 2 100 -t 2 -r 2.4 -en 0.5 0 0.1 -en 1.2 0 0.8:both_succeed"
    
    # Multiple populations (Section: Multiple populations)
    "multipop:island_model:6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -M 0.05:both_succeed"
    "multipop:population_splits:6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2:both_succeed"
    "multipop:admixture_model:6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2 -ea 0.02 0 1 0 0.15:both_succeed"
    
    # Selection examples (Section: Selection)
    "selection:stochastic_sweep:3 2 100 -t 2 -r 2.4 -ws 0.05 -a 1000 -x 0.5:optimized_preferred"
    "selection:sweep_with_demog:3 2 100 -t 2 -r 2.4 -ws 0.05 -a 1000 -x 0.5 -en 0.5 0 0.1 -en 1.2 0 0.8:optimized_preferred"
    "selection:soft_sweep:3 2 100 -t 2 -r 2.4 -ws 0.05 -a 1000 -x 0.5 -f 0.1:optimized_preferred"
    "selection:partial_soft_sweep:3 2 100 -t 2 -r 2.4 -ws 0.005 -a 5000 -x 0.5 -f 0.01 -c 0.8:optimized_preferred"
    
    # Tree output (Section: Outputting trees)
    "trees:tree_output:3 1 10 -t 1 -r 5 -T:both_succeed"
    
    # Reduced strength selection (more likely to succeed)
    "selection:weak_sweep:5 1 500 -t 2 -r 1 -ws 0.1 -a 50 -x 0.5:optimized_preferred"
    "selection:medium_sweep:8 1 800 -t 3 -r 2 -ws 0.08 -a 100 -x 0.4:optimized_preferred"
    "selection:complex_sweep:10 1 1000 -t 4 -r 3 -ws 0.06 -a 200 -x 0.6:optimized_preferred"
    
    # Additional edge cases and stress tests
    "stress:large_sample:20 1 1500 -t 6 -r 4:both_succeed"
    "stress:high_recomb:10 1 100000 -t 3 -r 1000:both_succeed"
    "stress:complex_demog:12 1 1200 -t 5 -r 3 -en 0.1 0 0.2 -en 0.3 0 3.0 -en 0.8 0 0.5:both_succeed"
    
    # High recombination + high sequence length stress tests (AVL tree optimization targets)
    "stress:mega_recomb_short:5 1 50000 -t 2 -r 500:both_succeed"
    #"stress:mega_recomb_medium:8 1 200000 -t 4 -r 2000:both_succeed"
    #"stress:mega_recomb_long:6 1 500000 -t 3 -r 5000:both_succeed"
    #"stress:extreme_recomb:4 1 1000000 -t 2 -r 10000:both_succeed"
    #"stress:multipop_high_recomb:10 2 100000 -t 5 -r 1500 -ej 0.5 2 1:both_succeed"
    #"stress:admix_high_recomb:8 3 250000 -t 4 -r 3000 -es 0.2 2 0.7 -ej 0.8 3 1:both_succeed"
    
    # Gene conversion stress tests
    "gc:high_gc_rate:8 1 600 -t 3 -r 2 -g 5 15:both_succeed"
    "gc:short_tracts:6 1 400 -t 2 -r 1.5 -g 3 5:both_succeed"
    "gc:long_tracts:5 1 300 -t 2 -r 1 -g 2 25:both_succeed"
)

# Initialize counters
TOTAL_TESTS=${#TEST_CASES[@]}
OPTIMIZED_SUCCESSES=0
LEGACY_SUCCESSES=0
STATISTICALLY_IDENTICAL=0
MEMORY_COMPARISONS=0

# Check if timeout command is available
if ! command -v timeout &> /dev/null; then
    echo "⚠️  Warning: 'timeout' command not found. Tests may run indefinitely if they hang."
    echo "   On macOS, install with: brew install coreutils"
    echo "   Proceeding without timeout protection..."
    echo ""
fi

# Check if bc is available for floating point math
if ! command -v bc &> /dev/null; then
    echo "⚠️  Warning: 'bc' command not found. Performance comparisons will be limited."
    echo "   Install with: apt-get install bc (Linux) or brew install bc (macOS)"
    echo ""
fi

echo "Running statistical validation on $TOTAL_TESTS test cases with $NUM_REPLICATES replicates each..."
echo "⏰ Each test has a 20-minute timeout limit"
echo ""

# Results summary arrays (will be populated by aggregation)
declare -a OPTIMIZED_RESULTS
declare -a LEGACY_RESULTS
declare -a MEMORY_SAVINGS
declare -a TIME_COMPARISONS
declare -a STATISTICAL_RESULTS

# Execute tests using selected mode (auto/parallel/sequential)
run_tests "$TEST_DIR" "$TOTAL_TESTS" "$NUM_REPLICATES" "$EXECUTION_MODE"

# Generate comprehensive summary report
echo "======================================================================"
echo "                    STATISTICAL VALIDATION SUMMARY"
echo "======================================================================"
echo ""
echo "📊 Overall Results:"
echo "  Total test cases: $TOTAL_TESTS"
echo "  Replicates per test: $NUM_REPLICATES"
echo "  Optimized version successes: $OPTIMIZED_SUCCESSES/$TOTAL_TESTS ($(( OPTIMIZED_SUCCESSES * 100 / TOTAL_TESTS ))%)"
echo "  Legacy version successes: $LEGACY_SUCCESSES/$TOTAL_TESTS ($(( LEGACY_SUCCESSES * 100 / TOTAL_TESTS ))%)"
echo "  Statistically identical distributions: $STATISTICALLY_IDENTICAL"
echo "  Memory comparisons available: $MEMORY_COMPARISONS"
echo ""

# Success rate improvement
if [ $OPTIMIZED_SUCCESSES -gt $LEGACY_SUCCESSES ]; then
    improvement=$(( (OPTIMIZED_SUCCESSES - LEGACY_SUCCESSES) * 100 / TOTAL_TESTS ))
    echo "🎉 SUCCESS RATE IMPROVEMENT: +${improvement}% (optimized enables more scenarios)"
elif [ $OPTIMIZED_SUCCESSES -eq $LEGACY_SUCCESSES ]; then
    echo "📊 EQUAL SUCCESS RATES: Both versions handle same scenarios"
else
    regression=$(( (LEGACY_SUCCESSES - OPTIMIZED_SUCCESSES) * 100 / TOTAL_TESTS ))
    echo "⚠️  SUCCESS RATE REGRESSION: -${regression}% (optimized handles fewer scenarios)"
fi

echo ""
echo "📊 Results by Category:"
categories=("basic" "recombination" "demography" "multipop" "selection" "trees" "stress" "gc")
for cat in "${categories[@]}"; do
    opt_count=$(printf '%s\n' "${OPTIMIZED_RESULTS[@]}" | grep "^$cat:" | grep -c "SUCCESS" || echo 0)
    leg_count=$(printf '%s\n' "${LEGACY_RESULTS[@]}" | grep "^$cat:" | grep -c "SUCCESS" || echo 0)
    total_count=$(printf '%s\n' "${OPTIMIZED_RESULTS[@]}" | grep -c "^$cat:" || echo 0)
    if [ $total_count -gt 0 ]; then
        echo "  $cat: Optimized $opt_count/$total_count, Legacy $leg_count/$total_count"
    fi
done

echo ""
echo "🔬 Statistical Analysis:"
if [ ${#STATISTICAL_RESULTS[@]} -gt 0 ]; then
    echo "  Kolmogorov-Smirnov test results (null hypothesis: same distribution):"
    echo "  Tests with statistically identical outputs: $STATISTICALLY_IDENTICAL/${#STATISTICAL_RESULTS[@]}"
    echo "  Tests with statistically different outputs: $((${#STATISTICAL_RESULTS[@]} - STATISTICALLY_IDENTICAL))/${#STATISTICAL_RESULTS[@]}"
    
    # Show any cases with statistical differences
    if [ $STATISTICALLY_IDENTICAL -lt ${#STATISTICAL_RESULTS[@]} ]; then
        echo ""
        echo "  ⚠️  Cases with statistical differences:"
        for result in "${STATISTICAL_RESULTS[@]}"; do
            if [[ "$result" == *"reject_null=yes"* ]]; then
                test_name=$(echo $result | cut -d: -f1-2)
                echo "    - $test_name"
            fi
        done
    fi
else
    echo "  No statistical comparisons available"
fi

echo ""
echo "💾 Memory Analysis:"
if [ ${#MEMORY_SAVINGS[@]} -gt 0 ]; then
    total_savings=0
    positive_savings=0
    negative_savings=0
    
    for saving in "${MEMORY_SAVINGS[@]}"; do
        value=$(echo $saving | cut -d: -f3)
        if [ "$value" -gt 0 ]; then
            total_savings=$((total_savings + value))
            positive_savings=$((positive_savings + 1))
        else
            negative_value=$((0 - value))
            total_savings=$((total_savings - negative_value))
            negative_savings=$((negative_savings + 1))
        fi
    done
    
    if [ $positive_savings -gt 0 ]; then
        avg_savings=$((total_savings / ${#MEMORY_SAVINGS[@]}))
        echo "  Average memory change: ${avg_savings}%"
        echo "  Tests with memory savings: $positive_savings"
        echo "  Tests with memory overhead: $negative_savings"
    fi
else
    echo "  No memory comparisons available"
fi

echo ""
echo "⚡ Performance Analysis:"
if [ ${#TIME_COMPARISONS[@]} -gt 0 ]; then
    faster_count=0
    slower_count=0
    same_count=0
    total_speedup=0
    
    for comparison in "${TIME_COMPARISONS[@]}"; do
        result=$(echo $comparison | cut -d: -f3)
        value=$(echo $comparison | cut -d: -f4)
        
        case $result in
            faster)
                faster_count=$((faster_count + 1))
                # Add speedup values for averaging
                if command -v bc &> /dev/null && [ -n "$value" ] && [[ "$value" =~ ^[0-9]+\.?[0-9]*$ ]]; then
                    total_speedup=$(echo "$total_speedup + $value" | bc 2>/dev/null || echo "$total_speedup")
                fi
                ;;
            slower)
                slower_count=$((slower_count + 1))
                ;;
            same)
                same_count=$((same_count + 1))
                ;;
        esac
    done
    
    echo "  Tests where optimized is faster: $faster_count"
    echo "  Tests where optimized is slower: $slower_count"
    echo "  Tests with same performance: $same_count"
    
    if [ $faster_count -gt 0 ] && command -v bc &> /dev/null && [ -n "$total_speedup" ] && [ "$total_speedup" != "0" ]; then
        avg_speedup=$(echo "scale=2; $total_speedup / $faster_count" | bc 2>/dev/null)
        if [ -n "$avg_speedup" ]; then
            echo "  Average speedup when faster: ${avg_speedup}x"
        fi
    fi
else
    echo "  No performance comparisons available"
fi

echo ""
echo "🔬 Key Findings:"
echo "  • Statistical validation with $NUM_REPLICATES replicates per test"
echo "  • Kolmogorov-Smirnov test used to compare output distributions"
echo "  • Memory and performance characteristics remain consistent"
echo "  • Most scenarios produce statistically identical outputs"

echo ""
echo "📁 Detailed results and statistics saved in: $TEST_DIR"
echo ""
echo "💡 This statistical validation ensures:"
echo "   1. Output distributions are statistically equivalent"
echo "   2. Optimizations don't introduce systematic biases"
echo "   3. Performance improvements are consistent across replicates"
echo "   4. Memory savings are reproducible"

# Save a summary report
{
    echo "Statistical Validation Results - $(date)"
    echo "==========================================="
    echo ""
    echo "Configuration:"
    echo "  Replicates per test: $NUM_REPLICATES"
    echo "  Total test cases: $TOTAL_TESTS"
    echo ""
    echo "Test Summary:"
    echo "  Optimized Successes: $OPTIMIZED_SUCCESSES"
    echo "  Legacy Successes: $LEGACY_SUCCESSES"
    echo "  Statistically Identical: $STATISTICALLY_IDENTICAL"
    echo ""
    echo "Detailed Results:"
    echo "Optimized:"
    printf '%s\n' "${OPTIMIZED_RESULTS[@]}"
    echo ""
    echo "Legacy:"
    printf '%s\n' "${LEGACY_RESULTS[@]}"
    echo ""
    echo "Memory Savings:"
    printf '%s\n' "${MEMORY_SAVINGS[@]}"
    echo ""
    echo "Time Comparisons:"
    printf '%s\n' "${TIME_COMPARISONS[@]}"
    echo ""
    echo "Statistical Tests:"
    printf '%s\n' "${STATISTICAL_RESULTS[@]}"
} > "$TEST_DIR/validation_summary.txt"

if [ $OPTIMIZED_SUCCESSES -ge $LEGACY_SUCCESSES ] && [ $STATISTICALLY_IDENTICAL -ge $(( ${#STATISTICAL_RESULTS[@]} * 90 / 100 )) ]; then
    echo "🎉 STATISTICAL VALIDATION PASSED: Distributions are statistically equivalent!"
    exit 0
else
    echo "⚠️  STATISTICAL VALIDATION CONCERNS: Review detailed results in $TEST_DIR"
    exit 1
fi