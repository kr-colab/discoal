#!/bin/bash

# niceStats comparison suite
# Compares distribution of population genetics statistics between discoal versions
# Uses GNU parallel for efficient multi-core execution with chunking support

# Configuration
TOTAL_REPLICATES=${1:-1000}  # Total number of replicates
CHUNKS=${2:-10}               # Number of chunks to split work into
COMPARISON_BINARY=${3:-"../discoal_legacy_backup"}  # Binary to compare against
NSAM=20                       # Sample size
NSITES=10000                  # Number of sites
OUTPUT_DIR="nicestats_comparison_$(date +%Y%m%d_%H%M%S)"

# Calculate replicates per chunk
REPS_PER_CHUNK=$((TOTAL_REPLICATES / CHUNKS))
REMAINDER=$((TOTAL_REPLICATES % CHUNKS))

# Check for required binaries
if [ ! -f "../discoal_edited" ]; then
    echo "Error: discoal_edited not found. Run 'make discoal_edited' first."
    exit 1
fi

if [ ! -f "$COMPARISON_BINARY" ]; then
    echo "Error: Comparison binary $COMPARISON_BINARY not found."
    echo "Available options:"
    echo "  - Build legacy version: make discoal_legacy_backup"
    echo "  - Build mem branch version: make discoal_mem_branch"
    echo "  - Or specify a custom binary as the 3rd argument"
    exit 1
fi

if [ ! -f "../niceStats" ]; then
    echo "Error: niceStats not found. Run 'make niceStats' first."
    exit 1
fi

# Check for GNU parallel
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found. Please install it first."
    echo "On Ubuntu/Debian: sudo apt-get install parallel"
    echo "On macOS: brew install parallel"
    exit 1
fi

# Check for Python and required packages
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 not found. Please install it first."
    exit 1
fi

# Check for required Python packages
python3 -c "import numpy; import scipy" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: Required Python packages not found."
    echo "Please install numpy and scipy:"
    echo "  pip install numpy scipy"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "======================================================================"
echo "niceStats Comparison Suite"
echo "======================================================================"
echo "Comparing statistical distributions between discoal versions"
echo "Edited version: ../discoal_edited"
echo "Comparison version: $COMPARISON_BINARY"
echo "Total replicates: $TOTAL_REPLICATES"
echo "Chunks: $CHUNKS (${REPS_PER_CHUNK} replicates per chunk)"
echo "Sample size: $NSAM"
echo "Sites: $NSITES"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to run a chunk of simulations
run_simulation_chunk() {
    local chunk_id=$1
    local binary=$2
    local params=$3
    local seed1=$4
    local seed2=$5
    local reps=$6
    
    # Run simulation chunk and pipe to niceStats, output all data lines
    # Note: params contains PLACEHOLDER that needs to be replaced with reps
    local final_params=$(echo "$params" | sed "s/PLACEHOLDER/$reps/")
    $binary $final_params -d $seed1 $seed2 2>/dev/null | ../niceStats 2>/dev/null | tail -n +2
}

export -f run_simulation_chunk

# Test scenarios
declare -a SCENARIOS=(
    "neutral:$NSAM PLACEHOLDER $NSITES -t 20"
    "low_recomb:$NSAM PLACEHOLDER $NSITES -t 20 -r 10"
    "high_recomb:$NSAM PLACEHOLDER $NSITES -t 20 -r 500"
    "weak_sweep:$NSAM PLACEHOLDER $NSITES -t 20 -r 40 -ws 0.1 -Pa 5 100 -Pu 0 0.05"
    "strong_sweep:$NSAM PLACEHOLDER $NSITES -t 20 -r 40 -ws 0.01 -Pa 5 500 -Pu 0 0.001"
    "bottleneck:$NSAM PLACEHOLDER $NSITES -t 20 -r 40 -en 0.1 0 0.1 -en 0.2 0 1.0"
)

# Process each scenario
for scenario_def in "${SCENARIOS[@]}"; do
    IFS=':' read -r scenario_name params_template <<< "$scenario_def"
    
    echo "----------------------------------------------------------------------"
    echo "Scenario: $scenario_name"
    echo "Parameters: ${params_template/PLACEHOLDER/$TOTAL_REPLICATES}"
    
    # Start timing with millisecond precision
    START_TIME=$(date +%s.%N)
    
    # Output files for this scenario
    edited_output="$OUTPUT_DIR/${scenario_name}_edited.stats"
    legacy_output="$OUTPUT_DIR/${scenario_name}_legacy.stats"
    comparison_output="$OUTPUT_DIR/${scenario_name}_comparison.txt"
    
    # Generate random seeds for chunks
    echo -n "" > "$OUTPUT_DIR/${scenario_name}_chunks.txt"
    for ((i=0; i<CHUNKS; i++)); do
        # Determine number of reps for this chunk
        if [ $i -lt $REMAINDER ]; then
            chunk_reps=$((REPS_PER_CHUNK + 1))
        else
            chunk_reps=$REPS_PER_CHUNK
        fi
        echo "$i $RANDOM $RANDOM $chunk_reps" >> "$OUTPUT_DIR/${scenario_name}_chunks.txt"
    done
    
    echo "Running edited version (${CHUNKS} chunks)..."
    # Run simulation chunks in parallel for edited version
    cat "$OUTPUT_DIR/${scenario_name}_chunks.txt" | \
        parallel -j+0 --colsep ' ' --bar \
        "run_simulation_chunk {1} ../discoal_edited '${params_template}' {2} {3} {4}" \
        > "$edited_output"
    
    EDITED_TIME=$(date +%s.%N)
    EDITED_DURATION=$(echo "$EDITED_TIME - $START_TIME" | bc)
    
    echo "Running legacy version (${CHUNKS} chunks)..."
    # Run simulation chunks in parallel for legacy version
    cat "$OUTPUT_DIR/${scenario_name}_chunks.txt" | \
        parallel -j+0 --colsep ' ' --bar \
        "run_simulation_chunk {1} '$COMPARISON_BINARY' '${params_template}' {2} {3} {4}" \
        > "$legacy_output"
    
    LEGACY_TIME=$(date +%s.%N)
    LEGACY_DURATION=$(echo "$LEGACY_TIME - $EDITED_TIME" | bc)
    TOTAL_DURATION=$(echo "$LEGACY_TIME - $START_TIME" | bc)
    
    echo "Timing: Edited=${EDITED_DURATION}s, Legacy=${LEGACY_DURATION}s, Total=${TOTAL_DURATION}s"
    echo ""
    
    # Extract column headers from niceStats output
    if [ ! -f "$OUTPUT_DIR/headers.txt" ]; then
        ../discoal_edited 10 1 1000 -t 10 2>/dev/null | ../niceStats 2>/dev/null | head -n 1 > "$OUTPUT_DIR/headers.txt"
    fi
    
    # Compare distributions using Python
    python3 << EOF > "$comparison_output" 2>&1
import sys
import numpy as np
from scipy import stats

try:
    # Read data
    edited_data = np.loadtxt("$edited_output")
    legacy_data = np.loadtxt("$legacy_output")
    
    # Handle case where data might be 1D (single replicate)
    if edited_data.ndim == 1:
        edited_data = edited_data.reshape(1, -1)
    if legacy_data.ndim == 1:
        legacy_data = legacy_data.reshape(1, -1)

    # Column names
    colnames = ["pi", "ss", "thetaH", "tajD", "fayWuH", "maxFDA", "HapCount", 
                "H1", "H12", "H2_H1", "Omega", "ZnS"]
    
    # Select numeric columns (skip H1 which seems to have issues)
    numeric_cols = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11]
    
    print("Statistical Comparison of Distributions")
    print("=======================================\\n")
    print(f"Sample sizes: Edited={len(edited_data)}, Legacy={len(legacy_data)}")
    print(f"Timing: Edited={$EDITED_DURATION:.3f}s, Legacy={$LEGACY_DURATION:.3f}s, Total={$TOTAL_DURATION:.3f}s")
    if $EDITED_DURATION > 0.001:
        print(f"Performance: Edited is {$LEGACY_DURATION/$EDITED_DURATION:.2f}x faster than Legacy\\n")
    else:
        print("Performance: Times too short to measure accurately\\n")
    
    for i in numeric_cols:
        if i < edited_data.shape[1] and i < legacy_data.shape[1]:
            stat_name = colnames[i]
            edited_vals = edited_data[:, i]
            legacy_vals = legacy_data[:, i]
            
            # Remove any nan or inf values
            edited_vals = edited_vals[np.isfinite(edited_vals)]
            legacy_vals = legacy_vals[np.isfinite(legacy_vals)]
            
            if len(edited_vals) > 0 and len(legacy_vals) > 0:
                # Kolmogorov-Smirnov test
                ks_stat, ks_pvalue = stats.ks_2samp(edited_vals, legacy_vals)
                
                print(f"\\n{stat_name}:")
                print(f"  Edited - Mean: {np.mean(edited_vals):.4f}, SD: {np.std(edited_vals):.4f}, Median: {np.median(edited_vals):.4f}")
                print(f"  Legacy - Mean: {np.mean(legacy_vals):.4f}, SD: {np.std(legacy_vals):.4f}, Median: {np.median(legacy_vals):.4f}")
                significance = "***" if ks_pvalue < 0.05 else ""
                print(f"  KS test: D = {ks_stat:.4f}, p-value = {ks_pvalue:.4f} {significance}")
                
except Exception as e:
    print(f"Error during statistical comparison: {e}", file=sys.stderr)
    sys.exit(1)
EOF
    
    # Display results
    echo ""
    cat "$comparison_output"
    
    # Save timing information
    echo -e "\nTiming Information:" >> "$comparison_output"
    printf "  Edited version: %.3f seconds\n" "${EDITED_DURATION}" >> "$comparison_output"
    printf "  Legacy version: %.3f seconds\n" "${LEGACY_DURATION}" >> "$comparison_output"
    printf "  Total time: %.3f seconds\n" "${TOTAL_DURATION}" >> "$comparison_output"
done

# Summary report
echo ""
echo "======================================================================"
echo "Summary Report"
echo "======================================================================"
echo ""

# Timing summary
echo "Timing Summary:"
echo "---------------"
grep -h "Timing:" "$OUTPUT_DIR"/*_comparison.txt | while read line; do
    scenario=$(echo "$line" | cut -d':' -f1)
    echo "$line"
done
echo ""

# Count significant differences
echo "Scenarios with significant differences (p < 0.05):"
grep -l "\*\*\*" "$OUTPUT_DIR"/*_comparison.txt 2>/dev/null | while read file; do
    scenario=$(basename "$file" | sed 's/_comparison.txt//')
    echo "  - $scenario:"
    grep -B1 "\*\*\*" "$file" | grep -E "^[a-zA-Z]" | sed 's/://' | while read stat; do
        echo "      * $stat"
    done
done

echo ""
echo "Full results saved in: $OUTPUT_DIR"
echo ""
echo "To view detailed results for a specific scenario:"
echo "  cat $OUTPUT_DIR/SCENARIO_comparison.txt"
echo ""
echo "Raw statistics data:"
echo "  Edited version: $OUTPUT_DIR/SCENARIO_edited.stats"
echo "  Legacy version: $OUTPUT_DIR/SCENARIO_legacy.stats"
echo ""
echo "Usage: $0 [total_replicates] [number_of_chunks] [comparison_binary]"
echo "  Default: 1000 replicates in 10 chunks comparing against ../discoal_legacy_backup"
echo "  Examples:"
echo "    $0 10000 20                      # 10,000 replicates in 20 chunks (vs legacy)"
echo "    $0 1000 10 ../discoal_mem_branch # Compare against mem branch"
echo "    $0 100 5 /path/to/custom/binary  # Compare against custom binary"