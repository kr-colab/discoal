#!/bin/bash

# Baseline performance benchmark for discoal
# Run before implementing optimizations to establish baseline

echo "=== Discoal Performance Baseline Benchmark ==="
echo "Date: $(date)"
echo "Commit: $(git rev-parse HEAD)"
echo ""

# Ensure we have a clean build
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1

# Test configurations
declare -a configs=(
    "20 1000 10000 -t 20 -r 500"          # Original profiling params
    "20 1000 10000 -t 20 -r 5000"         # 10x recombination
    "20 1000 10000 -t 100 -r 500"         # 5x mutation  
    "100 100 10000 -t 20 -r 500"          # More samples
    "20 1000 100000 -t 20 -r 5000"        # 10x sites
    "20 10000 10000 -t 20 -r 500"         # 10x replicates
)

declare -a descriptions=(
    "Standard parameters"
    "High recombination"
    "High mutation"
    "Many samples"
    "Long sequences"
    "Many replicates"
)

# Results file
RESULTS="baseline_performance_$(date +%Y%m%d_%H%M%S).txt"

echo "Configuration|Description|Real Time (s)|User Time (s)|Max Memory (KB)" > "$RESULTS"
echo "-------------|-----------|--------------|--------------|---------------" >> "$RESULTS"

# Run benchmarks
for i in "${!configs[@]}"; do
    config="${configs[$i]}"
    desc="${descriptions[$i]}"
    
    echo "Running: $desc ($config)"
    
    # Use GNU time for detailed stats
    if command -v gtime &> /dev/null; then
        TIME_CMD="gtime"
    else
        TIME_CMD="/usr/bin/time"
    fi
    
    # Run 3 times and take median
    times=()
    for run in 1 2 3; do
        # Capture time output
        time_output=$($TIME_CMD -f "real:%e user:%U maxrss:%M" ./discoal $config 2>&1 > /dev/null | tail -1)
        real_time=$(echo "$time_output" | grep -o "real:[0-9.]*" | cut -d: -f2)
        user_time=$(echo "$time_output" | grep -o "user:[0-9.]*" | cut -d: -f2)
        max_rss=$(echo "$time_output" | grep -o "maxrss:[0-9]*" | cut -d: -f2)
        
        times+=("$real_time")
        
        if [ $run -eq 2 ]; then  # Use middle run
            echo "$config|$desc|$real_time|$user_time|$max_rss" >> "$RESULTS"
        fi
    done
    
    # Calculate average
    avg=$(echo "${times[@]}" | awk '{sum=0; for(i=1;i<=NF;i++)sum+=$i; print sum/NF}')
    echo "  Average time: ${avg}s"
    echo ""
done

echo "Results saved to: $RESULTS"
echo ""
echo "Summary:"
column -t -s "|" "$RESULTS"

# Create performance profile for hottest functions
echo ""
echo "Creating detailed performance profile..."
perf record -g ./discoal 20 1000 10000 -t 20 -r 500 > /dev/null 2>&1
perf report --stdio --no-children | head -30 > "baseline_profile_$(date +%Y%m%d_%H%M%S).txt"

echo "Baseline benchmarking complete!"