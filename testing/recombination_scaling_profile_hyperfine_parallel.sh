#!/bin/bash

# Recombination Scaling Profile with Hyperfine and GNU Parallel
# Tests memory and time scaling across recombination rates
# Uses hyperfine for accurate timing and GNU parallel for efficiency

# Set up paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
EDITED_BIN="$PROJECT_ROOT/discoal_edited"
LEGACY_BIN="$PROJECT_ROOT/discoal_legacy_backup"

# Create output directory
OUTPUT_DIR="$SCRIPT_DIR/recomb_scaling_hyperfine_parallel_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

echo "=== RECOMBINATION SCALING PROFILER (HYPERFINE + PARALLEL) ==="
echo "Testing memory and time scaling across recombination rates"
echo "Output directory: $OUTPUT_DIR"
echo

# Check dependencies
if ! command -v hyperfine &> /dev/null; then
    echo "Error: hyperfine not found. Install with:"
    echo "  Ubuntu/Debian: wget https://github.com/sharkdp/hyperfine/releases/download/v1.16.1/hyperfine_1.16.1_amd64.deb && sudo dpkg -i hyperfine_1.16.1_amd64.deb"
    echo "  macOS: brew install hyperfine"
    echo "  Or: cargo install hyperfine"
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel not found. Install with:"
    echo "  Ubuntu/Debian: sudo apt-get install parallel"
    echo "  macOS: brew install parallel"
    exit 1
fi

# Check if binaries exist
if [ ! -f "$EDITED_BIN" ]; then
    echo "Error: $EDITED_BIN not found. Please build with 'make discoal_edited'"
    exit 1
fi

if [ ! -f "$LEGACY_BIN" ]; then
    echo "Error: $LEGACY_BIN not found. Please build with 'make discoal_legacy_backup'"
    exit 1
fi

# Base command parameters
if [ $# -eq 0 ]; then
    # Default parameters for measurable timing
    SAMPLES=50
    REPS=1
    SITES=200000
    THETA=20
    BASE_CMD="$SAMPLES $REPS $SITES -t $THETA"
else
    # Use provided command line
    BASE_CMD="$@"
    SAMPLES=$(echo "$@" | awk '{print $1}')
    SITES=$(echo "$@" | awk '{print $3}')
fi

# Define recombination rates to test
REC_RATES=(10 20 50 100 200 500 1000)

# Output files
RESULTS_FILE="$OUTPUT_DIR/scaling_results.csv"
SUMMARY_FILE="$OUTPUT_DIR/scaling_summary.txt"
TEMP_DIR="$OUTPUT_DIR/temp"
mkdir -p "$TEMP_DIR"

# Write CSV header
echo "rec_rate,version,mean_time,stddev_time,min_time,max_time,memory_kb" > "$RESULTS_FILE"

# Function to run hyperfine for a single recombination rate
run_single_benchmark() {
    local rec_rate=$1
    local output_dir=$2
    local edited_bin=$3
    local legacy_bin=$4
    local base_cmd="$5"
    local temp_dir=$6
    
    echo "Testing r=$rec_rate..."
    
    # Commands to test
    local edited_cmd="$edited_bin $base_cmd -r $rec_rate -d 12345 67890"
    local legacy_cmd="$legacy_bin $base_cmd -r $rec_rate -d 12345 67890"
    
    # Run hyperfine benchmark
    local hyperfine_output="$output_dir/hyperfine_r${rec_rate}.json"
    
    hyperfine \
        --warmup 3 \
        --min-runs 10 \
        --export-json "$hyperfine_output" \
        --command-name "edited" "$edited_cmd" \
        --command-name "legacy" "$legacy_cmd" \
        > "$output_dir/hyperfine_r${rec_rate}.log" 2>&1
    
    # Get memory usage
    echo "  Measuring memory usage..."
    local mem_output_edited="$temp_dir/mem_edited_r${rec_rate}.txt"
    local mem_output_legacy="$temp_dir/mem_legacy_r${rec_rate}.txt"
    
    /usr/bin/time -v $edited_cmd > /dev/null 2> "$mem_output_edited"
    local edited_mem=$(grep "Maximum resident" "$mem_output_edited" | awk '{print $6}')
    
    /usr/bin/time -v $legacy_cmd > /dev/null 2> "$mem_output_legacy"
    local legacy_mem=$(grep "Maximum resident" "$mem_output_legacy" | awk '{print $6}')
    
    # Debug: show memory values
    echo "    Edited memory: ${edited_mem:-0} KB"
    echo "    Legacy memory: ${legacy_mem:-0} KB"
    
    # Set default values if empty
    edited_mem=${edited_mem:-0}
    legacy_mem=${legacy_mem:-0}
    
    # Parse hyperfine JSON and create temporary result file
    if [ -f "$hyperfine_output" ]; then
        python3 << EOF > "$temp_dir/result_${rec_rate}.csv"
import json
with open('$hyperfine_output') as f:
    data = json.load(f)
    for result in data['results']:
        if result['command'] == 'edited':
            print(f"$rec_rate,edited,{result['mean']},{result['stddev']},{result['min']},{result['max']},$edited_mem")
        elif result['command'] == 'legacy':
            print(f"$rec_rate,legacy,{result['mean']},{result['stddev']},{result['min']},{result['max']},$legacy_mem")
EOF
    fi
}

export -f run_single_benchmark
export OUTPUT_DIR
export EDITED_BIN
export LEGACY_BIN
export BASE_CMD
export TEMP_DIR

echo "Running benchmarks in parallel..."
echo "Testing ${#REC_RATES[@]} recombination rates"
echo

# Run benchmarks in parallel
printf '%s\n' "${REC_RATES[@]}" | parallel -j 4 --progress \
    "run_single_benchmark {} '$OUTPUT_DIR' '$EDITED_BIN' '$LEGACY_BIN' '$BASE_CMD' '$TEMP_DIR'"

# Combine results
echo "Combining results..."
for rec_rate in "${REC_RATES[@]}"; do
    if [ -f "$TEMP_DIR/result_${rec_rate}.csv" ]; then
        cat "$TEMP_DIR/result_${rec_rate}.csv" >> "$RESULTS_FILE"
    fi
done

# Sort results by recombination rate and version
echo "Sorting results..."
head -1 "$RESULTS_FILE" > "$OUTPUT_DIR/scaling_results_sorted.csv"
tail -n +2 "$RESULTS_FILE" | sort -t, -k1,1n -k2,2 >> "$OUTPUT_DIR/scaling_results_sorted.csv"
mv "$OUTPUT_DIR/scaling_results_sorted.csv" "$RESULTS_FILE"

# Clean up temp directory
rm -rf "$TEMP_DIR"

# Create Python script for plotting
cat > "$OUTPUT_DIR/plot_scaling.py" << 'EOF'
#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Read the CSV file
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, 'scaling_results.csv')
df = pd.read_csv(csv_file)

# Convert memory from KB to MB
df['memory_mb'] = df['memory_kb'] / 1024

# Filter out rows with zero memory (measurement failures)
df_valid_memory = df[df['memory_kb'] > 0].copy()

# Create figure with subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# Color scheme
edited_color = '#2E86AB'
legacy_color = '#E63946'
reduction_color = '#06D6A0'

# Plot 1: Memory usage
edited_data = df[df['version'] == 'edited']
legacy_data = df[df['version'] == 'legacy']
edited_mem_data = df_valid_memory[df_valid_memory['version'] == 'edited']
legacy_mem_data = df_valid_memory[df_valid_memory['version'] == 'legacy']

if not edited_mem_data.empty:
    ax1.plot(edited_mem_data['rec_rate'], edited_mem_data['memory_mb'], 
             marker='o', linewidth=2.5, markersize=8,
             label='Edited version', color=edited_color)
if not legacy_mem_data.empty:
    ax1.plot(legacy_mem_data['rec_rate'], legacy_mem_data['memory_mb'], 
             marker='s', linewidth=2.5, markersize=8,
             label='Legacy version', color=legacy_color)

ax1.set_xscale('log')
if df['memory_mb'].max() / df['memory_mb'].min() > 10:
    ax1.set_yscale('log')
ax1.set_xlabel('Recombination Rate', fontsize=12)
ax1.set_ylabel('Memory Usage (MB)', fontsize=12)
ax1.set_title('Memory Scaling with Recombination Rate', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3, which='both')
ax1.legend(fontsize=10)

# Add memory reduction annotations
if not edited_mem_data.empty and not legacy_mem_data.empty:
    for i in range(len(edited_mem_data)):
        rec_rate = edited_mem_data.iloc[i]['rec_rate']
        if rec_rate in legacy_mem_data['rec_rate'].values:
            edited_mem = edited_mem_data.iloc[i]['memory_mb']
            legacy_mem = legacy_mem_data[legacy_mem_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
            if legacy_mem > 0:
                reduction = (1 - edited_mem / legacy_mem) * 100
                if abs(reduction) > 5:  # Only show significant differences
                    ax1.annotate(f'{reduction:.0f}%', 
                               xy=(rec_rate, edited_mem),
                               xytext=(0, -15), textcoords='offset points',
                               ha='center', fontsize=8, color=reduction_color if reduction > 0 else 'red')

# Plot 2: Execution time with error bars
ax2.errorbar(edited_data['rec_rate'], edited_data['mean_time'], 
             yerr=edited_data['stddev_time'],
             marker='o', linewidth=2.5, markersize=8, capsize=5,
             label='Edited version', color=edited_color)
ax2.errorbar(legacy_data['rec_rate'], legacy_data['mean_time'], 
             yerr=legacy_data['stddev_time'],
             marker='s', linewidth=2.5, markersize=8, capsize=5,
             label='Legacy version', color=legacy_color)

ax2.set_xscale('log')
if df['mean_time'].max() > 0 and df['mean_time'].max() / df['mean_time'].min() > 10:
    ax2.set_yscale('log')
ax2.set_xlabel('Recombination Rate', fontsize=12)
ax2.set_ylabel('Execution Time (seconds)', fontsize=12)
ax2.set_title('Execution Time Scaling', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(fontsize=10)

# Add speedup annotations
for i in range(len(edited_data)):
    rec_rate = edited_data.iloc[i]['rec_rate']
    if rec_rate in legacy_data['rec_rate'].values:
        edited_time = edited_data.iloc[i]['mean_time']
        legacy_time = legacy_data[legacy_data['rec_rate'] == rec_rate]['mean_time'].values[0]
        if edited_time > 0:
            speedup = legacy_time / edited_time
            if speedup > 1.2:
                ax2.annotate(f'{speedup:.1f}x', 
                           xy=(rec_rate, edited_time),
                           xytext=(5, -15), textcoords='offset points',
                           fontsize=8, color='blue')

# Plot 3: Speedup bar chart
merged = pd.merge(edited_data[['rec_rate', 'mean_time']], 
                  legacy_data[['rec_rate', 'mean_time']], 
                  on='rec_rate', suffixes=('_edited', '_legacy'))
merged['speedup'] = merged.apply(lambda r: r['mean_time_legacy'] / r['mean_time_edited'] 
                                if r['mean_time_edited'] > 0 else 0, axis=1)

bars = ax3.bar(range(len(merged)), merged['speedup'], color='#4CAF50', alpha=0.8)
ax3.set_xlabel('Recombination Rate', fontsize=12)
ax3.set_ylabel('Speedup Factor', fontsize=12)
ax3.set_title('Performance Speedup (Legacy / Edited)', fontsize=14, fontweight='bold')
ax3.set_xticks(range(len(merged)))
ax3.set_xticklabels([str(int(r)) for r in merged['rec_rate']], rotation=45)
ax3.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
ax3.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, val in zip(bars, merged['speedup']):
    if val > 0:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f'{val:.1f}x', ha='center', va='bottom', fontsize=9)

# Plot 4: Summary table
ax4.axis('tight')
ax4.axis('off')

# Create summary table
table_data = []
headers = ['Rec Rate', 'Edited Time', 'Legacy Time', 'Speedup', 'Edited Mem', 'Legacy Mem', 'Mem Diff']

for i, row in merged.iterrows():
    rec_rate = int(row['rec_rate'])
    edited_time = edited_data[edited_data['rec_rate'] == rec_rate]['mean_time'].values[0]
    legacy_time = legacy_data[legacy_data['rec_rate'] == rec_rate]['mean_time'].values[0]
    edited_mem = edited_data[edited_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
    legacy_mem = legacy_data[legacy_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
    
    mem_diff = (1 - edited_mem/legacy_mem) * 100 if legacy_mem > 0 else 0
    
    table_data.append([
        f"{rec_rate}",
        f"{edited_time:.3f}s",
        f"{legacy_time:.3f}s",
        f"{row['speedup']:.1f}x" if row['speedup'] > 0 else "N/A",
        f"{edited_mem:.1f} MB",
        f"{legacy_mem:.1f} MB",
        f"{mem_diff:+.0f}%" if abs(mem_diff) > 1 else "~0%"
    ])

table = ax4.table(cellText=table_data, colLabels=headers, 
                  cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)

# Style the header
for i in range(len(headers)):
    table[(0, i)].set_facecolor('#4CAF50')
    table[(0, i)].set_text_props(weight='bold', color='white')

plt.tight_layout()

# Save plots
plot_file = os.path.join(script_dir, 'recombination_scaling_hyperfine.png')
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {plot_file}")

pdf_file = os.path.join(script_dir, 'recombination_scaling_hyperfine.pdf')
plt.savefig(pdf_file, bbox_inches='tight')
print(f"PDF saved to: {pdf_file}")

# Print detailed statistics
print("\n=== PERFORMANCE COMPARISON ===")
print(f"{'Rec Rate':<10} {'Edited (s)':<15} {'Legacy (s)':<15} {'Speedup':<10} {'Edited Mem':<15} {'Legacy Mem':<15} {'Mem Diff':<10}")
print("-" * 95)

for i, row in merged.iterrows():
    rec_rate = int(row['rec_rate'])
    edited_time = edited_data[edited_data['rec_rate'] == rec_rate]['mean_time'].values[0]
    legacy_time = legacy_data[legacy_data['rec_rate'] == rec_rate]['mean_time'].values[0]
    edited_mem = edited_data[edited_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
    legacy_mem = legacy_data[legacy_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
    
    mem_diff = (1 - edited_mem/legacy_mem) * 100 if legacy_mem > 0 else 0
    speedup = row['speedup'] if row['speedup'] > 0 else 0
    
    print(f"{rec_rate:<10} {edited_time:<15.4f} {legacy_time:<15.4f} "
          f"{speedup:<9.1f}x {edited_mem:<14.1f}MB {legacy_mem:<14.1f}MB "
          f"{mem_diff:+9.0f}%")

# Calculate averages
valid_speedups = merged[merged['speedup'] > 0]['speedup']
if not valid_speedups.empty:
    avg_speedup = valid_speedups.mean()
    print(f"\nAverage speedup: {avg_speedup:.1f}x")

# Check for significant memory differences
mem_diffs = []
for i in range(len(edited_data)):
    rec_rate = edited_data.iloc[i]['rec_rate']
    if rec_rate in legacy_data['rec_rate'].values:
        edited_mem = edited_data.iloc[i]['memory_mb']
        legacy_mem = legacy_data[legacy_data['rec_rate'] == rec_rate]['memory_mb'].values[0]
        if legacy_mem > 0:
            diff = (1 - edited_mem/legacy_mem) * 100
            if abs(diff) > 1:
                mem_diffs.append(diff)

if mem_diffs:
    print(f"Average memory difference: {np.mean(mem_diffs):+.1f}%")
else:
    print("No significant memory differences detected with current parameters")

plt.show()
EOF

# Make Python script executable
chmod +x "$OUTPUT_DIR/plot_scaling.py"

# Generate plots
echo
echo "Generating plots..."
cd "$OUTPUT_DIR"
python3 plot_scaling.py | tee -a "$SUMMARY_FILE"

# Final summary
echo
echo "📊 Results saved to: $OUTPUT_DIR"
echo "📈 Plots saved as: recombination_scaling_hyperfine.png and .pdf"
echo "📋 Raw data in: scaling_results.csv"
echo "🔍 Detailed hyperfine results in: hyperfine_r*.json"
echo
echo "To regenerate plots: cd $OUTPUT_DIR && python3 plot_scaling.py"
echo
echo "Note: Running benchmarks in parallel with -j 4 to balance accuracy and speed"
echo "      Adjust -j parameter in the script if needed (current: 4 concurrent jobs)"