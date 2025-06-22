#!/bin/bash

# Growth comparison suite - Compare discoal -eg against ms -eg baseline
# Uses niceStats to compare summary statistics distributions

# Configuration
TOTAL_REPLICATES=${1:-1000}
CHUNKS=${2:-10}
NSAM=20
NSITES=10000
OUTPUT_DIR="growth_comparison_$(date +%Y%m%d_%H%M%S)"

# Calculate replicates per chunk
REPS_PER_CHUNK=$((TOTAL_REPLICATES / CHUNKS))
REMAINDER=$((TOTAL_REPLICATES % CHUNKS))

# Check for required binaries
echo "Checking for required binaries..."
if [ ! -f "../discoal" ]; then
    echo "Error: discoal not found. Run 'make discoal' first."
    exit 1
fi

if [ ! -f "../extern/ms" ]; then
    echo "Error: ms not found. Run 'make extern/ms' first."
    exit 1
fi

if [ ! -f "../niceStats" ]; then
    echo "Error: niceStats not found. Run 'make niceStats' first."
    exit 1
fi

# Check for GNU parallel
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found. Install with: sudo apt-get install parallel"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Growth comparison suite: $TOTAL_REPLICATES replicates in $CHUNKS chunks"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run simulation chunk and pipe to niceStats
run_simulation_chunk() {
    local chunk_id=$1
    local binary=$2
    local params=$3
    local seed1=$4
    local seed2=$5
    local reps=$6
    
    # Replace PLACEHOLDER with actual reps
    local final_params=$(echo "$params" | sed "s/PLACEHOLDER/$reps/")
    
    # Execute command with explicit path and error handling
    # Note: discoal uses -d for seeds, ms uses a seedms file
    if [[ "$binary" == *"ms"* ]]; then
        # For ms, set up random seed file for each chunk
        local seedfile="seedms_${chunk_id}_$$"
        echo "$seed1 $seed2" > "$seedfile"
        timeout 300 $binary $final_params 2>/dev/null | timeout 30 ../niceStats 2>/dev/null | tail -n +2
        rm -f "$seedfile"
    else
        # For discoal, use -d flag
        timeout 300 $binary $final_params -d $seed1 $seed2 2>/dev/null | timeout 30 ../niceStats 2>/dev/null | tail -n +2
    fi
    
    # Check for timeout or other failures
    if [ ${PIPESTATUS[0]} -eq 124 ] || [ ${PIPESTATUS[1]} -eq 124 ]; then
        echo "Error: Command timed out" >&2
        return 1
    fi
}

export -f run_simulation_chunk

# Growth test scenarios (ms format: nsam nreps -t theta -r rho nsites [growth_options])
# Note: ms uses 1-indexed populations, discoal uses 1-indexed populations internally
declare -a SCENARIOS=(
    "no_growth:$NSAM PLACEHOLDER -t 20 -r 40 $NSITES"
    "pop_growth_single:$NSAM PLACEHOLDER -t 20 -r 40 $NSITES -eg 0.0 1 1.0"
    "pop_growth_structured:$NSAM PLACEHOLDER -t 20 -r 40 $NSITES -I 2 10 10 1.0 -eg 0.0 1 1.0"
)

# Convert discoal params to ms params
convert_to_ms_params() {
    local params="$1"
    # Convert discoal syntax to ms syntax
    # discoal: -p 2 10 10 -> ms: -I 2 10 10 (same syntax)
    # Note: ms uses 1-indexed populations, discoal uses 0-indexed populations
    echo "$params" | sed -E 's/-p ([0-9]+) /-I \1 /g'
}

# Convert ms params to discoal params
convert_to_discoal_params() {
    local params="$1"
    # ms: nsam nreps -t theta -r rho nsites [options]
    # discoal: nsam nreps nsites -t theta -r rho [options]
    # Also convert: -I -> -p, -eG -> -G, -eg popID -> -eg popID+1 (1-indexed)
    
    # Extract nsites from -r rho nsites pattern and move it
    # ms: nsam nreps -t theta -r rho nsites [options]
    # discoal: nsam nreps nsites -t theta -r rho [options]
    local converted="$params"
    
    # Use awk to properly parse and rearrange parameters
    if echo "$params" | grep -q '\-r'; then
        converted=$(echo "$params" | awk '{
            nsam = $1; nreps = $2
            # Find -r flag position
            for (i = 3; i <= NF; i++) {
                if ($i == "-r") {
                    rho = $(i+1)
                    nsites = $(i+2)
                    # Build the new command
                    result = nsam " " nreps " " nsites
                    # Add everything before -r
                    for (j = 3; j < i; j++) result = result " " $j
                    # Add -r rho
                    result = result " -r " rho
                    # Add everything after nsites
                    for (j = i+3; j <= NF; j++) result = result " " $j
                    print result
                    exit
                }
            }
            # If no -r found, print as is
            print $0
        }')
    fi
    
    # Convert population structure: ms -I npop n1 n2 ... migrate -> discoal -p npop n1 n2 ... -M migrate
    # Use sed to handle the conversion more reliably
    converted=$(echo "$converted" | sed -E 's/-I ([0-9]+) ([0-9 ]+) ([0-9.]+)/-p \1 \2 -M \3/g')
    # Also handle simpler case without migration rate
    converted=$(echo "$converted" | sed -E 's/-I ([0-9]+) /-p \1 /g')
    
    # Convert population indices: ms uses 1-indexed, discoal uses 0-indexed
    # -eg time popID alpha -> -eg time (popID-1) alpha
    converted=$(echo "$converted" | awk '{
        for(i=1; i<=NF; i++) {
            if($i == "-eg" && i+2 <= NF) {
                printf "%s %s %d %s ", $i, $(i+1), $(i+2)-1, $(i+3)
                i += 3
            } else if($i == "-en" && i+3 <= NF) {
                printf "%s %s %d %s ", $i, $(i+1), $(i+2)-1, $(i+3)
                i += 3
            } else if($i == "-ej" && i+3 <= NF) {
                printf "%s %s %d %d ", $i, $(i+1), $(i+2)-1, $(i+3)-1
                i += 3
            } else if($i == "-em" && i+3 <= NF) {
                printf "%s %s %d %d %s ", $i, $(i+1), $(i+2)-1, $(i+3)-1, $(i+4)
                i += 4
            } else if($i == "-eM" && i+2 <= NF) {
                printf "%s %s %s ", $i, $(i+1), $(i+2)
                i += 2
            } else {
                printf "%s ", $i
            }
        }
        printf "\n"
    }')
    
    echo "$converted"
}

# Process each scenario
for scenario_def in "${SCENARIOS[@]}"; do
    IFS=':' read -r scenario_name ms_params_template <<< "$scenario_def"
    
    echo "----------------------------------------------------------------------"
    echo "Scenario: $scenario_name"
    echo "MS Parameters: ${ms_params_template/PLACEHOLDER/$TOTAL_REPLICATES}"
    
    START_TIME=$(date +%s.%N)
    
    # Output files
    discoal_output="$OUTPUT_DIR/${scenario_name}_discoal.stats"
    ms_output="$OUTPUT_DIR/${scenario_name}_ms.stats"
    comparison_output="$OUTPUT_DIR/${scenario_name}_comparison.txt"
    
    # Generate random seeds for chunks
    echo -n "" > "$OUTPUT_DIR/${scenario_name}_chunks.txt"
    for ((i=0; i<CHUNKS; i++)); do
        if [ $i -lt $REMAINDER ]; then
            chunk_reps=$((REPS_PER_CHUNK + 1))
        else
            chunk_reps=$REPS_PER_CHUNK
        fi
        echo "$i $RANDOM $RANDOM $chunk_reps" >> "$OUTPUT_DIR/${scenario_name}_chunks.txt"
    done
    
    # Convert ms parameters to discoal parameters
    discoal_params=$(convert_to_discoal_params "$ms_params_template")
    
    echo "Running discoal (${CHUNKS} chunks)..."
    echo "Discoal params: ${discoal_params/PLACEHOLDER/$TOTAL_REPLICATES}"
    cat "$OUTPUT_DIR/${scenario_name}_chunks.txt" | \
        parallel --no-notice -j+0 --colsep ' ' \
        "run_simulation_chunk {1} ../discoal '${discoal_params}' {2} {3} {4}" \
        > "$discoal_output" 2>/dev/null
    
    DISCOAL_TIME=$(date +%s.%N)
    DISCOAL_DURATION=$(echo "$DISCOAL_TIME - $START_TIME" | bc)
    
    echo "Running ms baseline (${CHUNKS} chunks)..."
    cat "$OUTPUT_DIR/${scenario_name}_chunks.txt" | \
        parallel --no-notice -j+0 --colsep ' ' \
        "run_simulation_chunk {1} ../extern/ms '${ms_params_template}' {2} {3} {4}" \
        > "$ms_output" 2>/dev/null
    
    MS_TIME=$(date +%s.%N)
    MS_DURATION=$(echo "$MS_TIME - $DISCOAL_TIME" | bc)
    TOTAL_DURATION=$(echo "$MS_TIME - $START_TIME" | bc)
    
    echo "Timing: discoal=${DISCOAL_DURATION}s, ms=${MS_DURATION}s, Total=${TOTAL_DURATION}s"
    echo ""
    
    # Compare distributions using detailed statistical tests
    python3 << EOF 2>/dev/null || echo "Python comparison failed"
import sys
import numpy as np
from scipy import stats

# Read data
discoal_data = np.loadtxt('$discoal_output')
ms_data = np.loadtxt('$ms_output')

if discoal_data.size == 0 or ms_data.size == 0:
    print('Error: Empty data files')
    sys.exit(1)

# Ensure 2D arrays
if discoal_data.ndim == 1:
    discoal_data = discoal_data.reshape(1, -1)
if ms_data.ndim == 1:
    ms_data = ms_data.reshape(1, -1)

print('Scenario: $scenario_name')
print('Discoal shape:', discoal_data.shape, 'MS shape:', ms_data.shape)
print('')

# Get column names from niceStats header
headers = ['SS', 'H', 'Pi', 'ThetaW', 'ThetaH', 'ThetaL', 'TajimasD', 'H1', 'H12', 'H2H1', 'VarPi', 'VarH', 'VarThetaW', 'VarThetaH', 'VarThetaL']

# Print detailed comparison to console
print('Detailed Statistical Comparison:')
print('=' * 80)
print(f"{'Statistic':<12} {'Discoal':<25} {'MS':<25} {'KS_pval':<10} {'Status':<6}")
print(f"{'':^12} {'Mean+/-SD [Median]':<25} {'Mean+/-SD [Median]':<25} {'':^10} {'':^6}")
print('-' * 80)

with open('$comparison_output', 'w') as f:
    f.write('Statistic\\tDiscoal_Mean\\tDiscoal_Median\\tDiscoal_StdDev\\tMS_Mean\\tMS_Median\\tMS_StdDev\\tDifference\\tPercent_Diff\\tKS_pvalue\\tStatus\\n')
    
    n_stats = min(discoal_data.shape[1], ms_data.shape[1], len(headers))
    all_pass = True
    
    for i in range(n_stats):
        stat_name = headers[i] if i < len(headers) else f'Stat_{i}'
        
        # Calculate detailed statistics
        discoal_vals = discoal_data[:, i]
        ms_vals = ms_data[:, i]
        
        # Remove any nan or inf values
        discoal_vals = discoal_vals[np.isfinite(discoal_vals)]
        ms_vals = ms_vals[np.isfinite(ms_vals)]
        
        if len(discoal_vals) == 0 or len(ms_vals) == 0:
            continue
        
        discoal_mean = np.mean(discoal_vals)
        discoal_median = np.median(discoal_vals)
        discoal_std = np.std(discoal_vals)
        
        ms_mean = np.mean(ms_vals)
        ms_median = np.median(ms_vals)
        ms_std = np.std(ms_vals)
        
        diff = discoal_mean - ms_mean
        
        if ms_mean != 0:
            pct_diff = abs(diff / ms_mean) * 100
        else:
            pct_diff = abs(diff) * 100 if diff != 0 else 0
        
        # Handle cases where all values are the same
        if discoal_std == 0 and ms_std == 0:
            if discoal_mean == ms_mean:
                ks_pval = 1.0
            else:
                ks_pval = 0.0
        else:
            _, ks_pval = stats.ks_2samp(discoal_vals, ms_vals)
        
        # Status: PASS if KS p-value > 0.01 (distributions not significantly different)
        status = 'PASS' if ks_pval > 0.01 else 'FAIL'
        if status == 'FAIL':
            all_pass = False
        
        # Print to console
        discoal_summary = f'{discoal_mean:.3f}+/-{discoal_std:.3f} [{discoal_median:.3f}]'
        ms_summary = f'{ms_mean:.3f}+/-{ms_std:.3f} [{ms_median:.3f}]'
        print(f'{stat_name:<12} {discoal_summary:<25} {ms_summary:<25} {ks_pval:<10.4f} {status:<6}')
        
        # Write to file
        f.write(f'{stat_name}\\t{discoal_mean:.6f}\\t{discoal_median:.6f}\\t{discoal_std:.6f}\\t{ms_mean:.6f}\\t{ms_median:.6f}\\t{ms_std:.6f}\\t{diff:.6f}\\t{pct_diff:.2f}%\\t{ks_pval:.6f}\\t{status}\\n')
    
    # Overall result
    overall = 'PASS' if all_pass else 'FAIL'
    f.write(f'OVERALL\\t-\\t-\\t-\\t-\\t-\\t-\\t-\\t-\\t-\\t{overall}\\n')
    
    print('-' * 80)
    print(f'Overall result: {overall}')
    print('')

EOF

done

echo ""
echo "========================================================================"
echo "Growth comparison suite completed"
echo "Results saved in: $OUTPUT_DIR"
echo ""
echo "Summary:"
find "$OUTPUT_DIR" -name "*_comparison.txt" -exec echo "{}:" \; -exec tail -n 1 {} \;