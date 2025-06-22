#!/bin/bash

# Memory profiling script for discoal with recombination
# Tests various recombination rates to understand memory usage

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="memory_profile_recomb_${TIMESTAMP}"
mkdir -p "$OUTDIR"

echo "🔬 Memory Profiling discoal with Recombination"
echo "=============================================="
echo "Output directory: $OUTDIR"
echo ""

# Test parameters
SAMPLE_SIZE=50
NUM_REPS=1
NUM_SITES=50000

# Function to run valgrind memory profiling
run_memory_profile() {
    local test_name=$1
    local rho=$2
    local extra_args=$3
    
    echo "📊 Running $test_name (ρ=$rho)..."
    
    # Run with massif for detailed heap profiling
    valgrind --tool=massif \
             --massif-out-file="$OUTDIR/massif_${test_name}.out" \
             --detailed-freq=1 \
             --max-snapshots=100 \
             --depth=40 \
             ../discoal $SAMPLE_SIZE $NUM_REPS $NUM_SITES -t 20 -r $rho $extra_args \
             > "$OUTDIR/output_${test_name}.ms" 2> "$OUTDIR/valgrind_${test_name}.log"
    
    # Extract peak memory usage
    local peak_mem=$(grep "mem_heap_B" "$OUTDIR/massif_${test_name}.out" | awk -F= '{print $2}' | sort -n | tail -1)
    local peak_mb=$(echo "scale=2; $peak_mem / 1048576" | bc)
    echo "  Peak memory: ${peak_mb} MB"
    
    # Run ms_print for detailed analysis
    ms_print "$OUTDIR/massif_${test_name}.out" > "$OUTDIR/massif_${test_name}.txt"
    
    # Also run with memcheck for leak detection
    valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             ../discoal $SAMPLE_SIZE $NUM_REPS $NUM_SITES -t 20 -r $rho $extra_args \
             > /dev/null 2> "$OUTDIR/memcheck_${test_name}.log"
    
    # Extract leak summary
    local definitely_lost=$(grep "definitely lost:" "$OUTDIR/memcheck_${test_name}.log" | awk '{print $4}')
    local indirectly_lost=$(grep "indirectly lost:" "$OUTDIR/memcheck_${test_name}.log" | awk '{print $4}')
    echo "  Memory leaks: definitely lost: $definitely_lost, indirectly lost: $indirectly_lost"
    echo ""
}

# Test different recombination rates
echo "🧬 Testing different recombination rates..."
echo "Sample size: $SAMPLE_SIZE, Sites: $NUM_SITES"
echo ""

run_memory_profile "no_recomb" 0
run_memory_profile "low_recomb" 10
run_memory_profile "medium_recomb" 50
run_memory_profile "high_recomb" 200
run_memory_profile "very_high_recomb" 1000
run_memory_profile "extreme_recomb" 5000

# Test with tree sequence output
echo "🌳 Testing with tree sequence output..."
run_memory_profile "tskit_low_recomb" 10 "-ts $OUTDIR/test_low.trees"
run_memory_profile "tskit_high_recomb" 200 "-ts $OUTDIR/test_high.trees"

# Test with larger sample size
echo "👥 Testing with larger sample size..."
SAMPLE_SIZE=200
run_memory_profile "large_sample_low" 10
run_memory_profile "large_sample_high" 200

# Generate summary report
echo "📝 Generating summary report..."
cat > "$OUTDIR/summary.txt" << EOF
Memory Profile Summary
=====================
Date: $(date)
discoal version: $(../discoal 2>&1 | head -1)

Test Parameters:
- Sample sizes: 50, 200
- Number of sites: $NUM_SITES
- Theta: 20

Peak Memory Usage:
EOF

# Extract peak memory for all tests
for massif_file in "$OUTDIR"/massif_*.out; do
    if [ -f "$massif_file" ]; then
        test_name=$(basename "$massif_file" .out | sed 's/massif_//')
        peak_mem=$(grep "mem_heap_B" "$massif_file" | awk -F= '{print $2}' | sort -n | tail -1)
        peak_mb=$(echo "scale=2; $peak_mem / 1048576" | bc)
        printf "%-25s %8.2f MB\n" "$test_name:" "$peak_mb" >> "$OUTDIR/summary.txt"
    fi
done

echo "" >> "$OUTDIR/summary.txt"
echo "Memory Leak Summary:" >> "$OUTDIR/summary.txt"
for memcheck_file in "$OUTDIR"/memcheck_*.log; do
    if [ -f "$memcheck_file" ]; then
        test_name=$(basename "$memcheck_file" .log | sed 's/memcheck_//')
        definitely_lost=$(grep "definitely lost:" "$memcheck_file" | awk '{print $4}')
        printf "%-25s %s bytes\n" "$test_name:" "$definitely_lost" >> "$OUTDIR/summary.txt"
    fi
done

echo ""
echo "✅ Memory profiling complete!"
echo "📁 Results saved in: $OUTDIR"
echo "📊 Summary available in: $OUTDIR/summary.txt"
echo ""
echo "To view detailed heap profiles:"
echo "  ms_print $OUTDIR/massif_<test_name>.out | less"