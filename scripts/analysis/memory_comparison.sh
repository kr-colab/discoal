#!/bin/bash

echo "=== Memory Usage Comparison: Standard vs Pooled ==="
echo "Date: $(date)"
echo ""
echo "Recombination,Standard Memory (MB),Pooled Memory (MB),Memory Overhead"

# Test different recombination rates
for rho in 500 1000 2500 5000 10000; do
    # Standard version
    std_output=$(/usr/bin/time -v ./discoal_standard 20 10 10000 -t 20 -r $rho 2>&1 > /dev/null)
    std_mem=$(echo "$std_output" | grep "Maximum resident" | awk '{print $6}')
    std_mem_mb=$(echo "scale=1; $std_mem / 1024" | bc)
    
    # Pooled version
    pool_output=$(/usr/bin/time -v ./discoal_pooled 20 10 10000 -t 20 -r $rho 2>&1 > /dev/null)
    pool_mem=$(echo "$pool_output" | grep "Maximum resident" | awk '{print $6}')
    pool_mem_mb=$(echo "scale=1; $pool_mem / 1024" | bc)
    
    # Calculate overhead
    overhead=$(echo "scale=1; ($pool_mem - $std_mem) * 100 / $std_mem" | bc)
    
    echo "$rho,$std_mem_mb,$pool_mem_mb,${overhead}%"
done