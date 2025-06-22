#!/bin/bash

# Benchmark script to compare segment pool vs standard allocation
# Run after implementing segment pool to measure improvement

echo "=== Segment Pool Performance Comparison ==="
echo "Date: $(date)"
echo ""

# Build both versions
echo "Building standard version..."
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1
mv discoal discoal_standard

echo "Building pooled version..."
make clean > /dev/null 2>&1
USE_SEGMENT_POOL=1 make discoal > /dev/null 2>&1
mv discoal discoal_pooled

# Test configurations focusing on allocation-heavy scenarios
declare -a configs=(
    "20 1000 10000 -t 20 -r 500"          # Baseline
    "20 1000 10000 -t 20 -r 5000"         # High recombination
    "20 5000 10000 -t 20 -r 1000"         # Many replicates  
    "50 500 10000 -t 50 -r 2000"          # Complex scenario
    "20 10000 5000 -t 10 -r 100"          # Many replicates, low recomb
)

declare -a descriptions=(
    "Baseline (moderate recombination)"
    "High recombination (10x)"
    "Many replicates (5x)"
    "Complex (more samples, high recomb)"
    "Many replicates, low recomb"
)

# Results file
RESULTS="pool_comparison_$(date +%Y%m%d_%H%M%S).csv"

echo "Config,Standard_Time,Pooled_Time,Speedup,Standard_Memory,Pooled_Memory,Memory_Reduction" > "$RESULTS"

# Function to get median of three runs
get_median_time() {
    local cmd="$1"
    local times=()
    
    for i in 1 2 3; do
        if command -v gtime &> /dev/null; then
            TIME_CMD="gtime"
        else
            TIME_CMD="/usr/bin/time"
        fi
        
        local output=$($TIME_CMD -f "%e,%M" $cmd 2>&1 > /dev/null | tail -1)
        local time=$(echo "$output" | cut -d',' -f1)
        local mem=$(echo "$output" | cut -d',' -f2)
        times+=("$time,$mem")
    done
    
    # Return the middle run
    echo "${times[1]}"
}

# Run comparisons
for i in "${!configs[@]}"; do
    config="${configs[$i]}"
    desc="${descriptions[$i]}"
    
    echo ""
    echo "Testing: $desc"
    echo "Command: ./discoal $config"
    
    # Test standard version
    echo -n "  Standard version: "
    standard_result=$(get_median_time "./discoal_standard $config")
    standard_time=$(echo "$standard_result" | cut -d',' -f1)
    standard_mem=$(echo "$standard_result" | cut -d',' -f2)
    echo "${standard_time}s, ${standard_mem}KB"
    
    # Test pooled version
    echo -n "  Pooled version:   "
    pooled_result=$(get_median_time "./discoal_pooled $config")
    pooled_time=$(echo "$pooled_result" | cut -d',' -f1)
    pooled_mem=$(echo "$pooled_result" | cut -d',' -f2)
    echo "${pooled_time}s, ${pooled_mem}KB"
    
    # Calculate improvements
    speedup=$(echo "scale=2; $standard_time / $pooled_time" | bc)
    mem_reduction=$(echo "scale=1; 100 - ($pooled_mem * 100 / $standard_mem)" | bc)
    
    echo "  Speedup: ${speedup}x"
    echo "  Memory reduction: ${mem_reduction}%"
    
    # Save results
    echo "\"$desc\",$standard_time,$pooled_time,$speedup,$standard_mem,$pooled_mem,$mem_reduction" >> "$RESULTS"
done

# Summary
echo ""
echo "=== Summary ==="
echo ""
column -t -s ',' "$RESULTS" | head -20

# Calculate average improvement
avg_speedup=$(tail -n +2 "$RESULTS" | awk -F',' '{sum+=$4; count++} END {print sum/count}')
avg_mem=$(tail -n +2 "$RESULTS" | awk -F',' '{sum+=$7; count++} END {print sum/count}')

echo ""
echo "Average speedup: ${avg_speedup}x"
echo "Average memory reduction: ${avg_mem}%"

# Profile the high recombination case
echo ""
echo "Profiling high recombination scenario..."
perf record -g ./discoal_standard 20 100 10000 -t 20 -r 5000 > /dev/null 2>&1
perf report --stdio --no-children | grep -E "(malloc|free|calloc)" > standard_malloc_profile.txt

perf record -g ./discoal_pooled 20 100 10000 -t 20 -r 5000 > /dev/null 2>&1
perf report --stdio --no-children | grep -E "(allocSegment|freeSegment|resetSegment)" > pooled_alloc_profile.txt

echo ""
echo "Results saved to: $RESULTS"
echo "Profile data saved to: standard_malloc_profile.txt, pooled_alloc_profile.txt"

# Cleanup
rm -f discoal_standard discoal_pooled perf.data perf.data.old