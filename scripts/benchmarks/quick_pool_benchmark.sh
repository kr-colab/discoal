#!/bin/bash

echo "=== Quick Pool Performance Test ==="
echo "Date: $(date)"
echo ""

# Function to time a command
time_command() {
    local cmd="$1"
    local start=$(date +%s.%N)
    $cmd > /dev/null 2>&1
    local end=$(date +%s.%N)
    echo "$end - $start" | bc
}

# Test configurations
declare -a configs=(
    "20 10 10000 -t 20 -r 500"          # Baseline
    "20 10 10000 -t 20 -r 5000"         # High recombination  
    "20 100 10000 -t 20 -r 500"         # Many replicates
)

declare -a descriptions=(
    "Baseline"
    "High recombination (10x)"
    "Many replicates (10x)"
)

echo "Config,Standard Time,Pooled Time,Speedup"

for i in "${!configs[@]}"; do
    config="${configs[$i]}"
    desc="${descriptions[$i]}"
    
    # Time standard version (3 runs, take middle)
    times=()
    for run in 1 2 3; do
        t=$(time_command "./discoal_standard $config")
        times+=($t)
    done
    std_time=${times[1]}
    
    # Time pooled version (3 runs, take middle)  
    times=()
    for run in 1 2 3; do
        t=$(time_command "./discoal_pooled $config")
        times+=($t)
    done
    pool_time=${times[1]}
    
    # Calculate speedup
    speedup=$(echo "scale=2; $std_time / $pool_time" | bc)
    
    echo "$desc,$std_time,$pool_time,$speedup"
done