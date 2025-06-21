#!/bin/bash
#
# Profile mutation placement algorithms comparing:
# 1. Legacy dropMutations() 
# 2. tskit edge-based algorithm
# 3. tskit node-based algorithm
#

set -e

echo "========================================"
echo "Mutation Algorithm Profiling"
echo "========================================"

# Create output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="mutation_profiling_${TIMESTAMP}"
mkdir -p "$OUTDIR"

# Build discoal
echo "Building discoal..."
cd ..
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

# Test parameters - larger scale for profiling
PARAMS_SMALL="50 10 10000 -t 100 -r 50"
PARAMS_MEDIUM="100 10 50000 -t 200 -r 100"
PARAMS_LARGE="200 5 100000 -t 400 -r 200"

profile_run() {
    local label=$1
    local params=$2
    local extra_args=$3
    local output_prefix=$4
    
    echo -n "  $label: "
    
    # Use perf stat for basic metrics
    if command -v perf &> /dev/null; then
        perf stat -o "testing/${OUTDIR}/${output_prefix}_perf.txt" \
            ./discoal $params $extra_args > "testing/${OUTDIR}/${output_prefix}.ms" 2>&1
    else
        # Fallback to time command
        /usr/bin/time -v -o "testing/${OUTDIR}/${output_prefix}_time.txt" \
            ./discoal $params $extra_args > "testing/${OUTDIR}/${output_prefix}.ms" 2>&1
    fi
    
    # Extract mutation count
    local muts=$(grep -m1 "segsites:" "testing/${OUTDIR}/${output_prefix}.ms" | awk '{print $2}')
    echo "$muts mutations"
}

# Small scale test
echo -e "\nSmall scale test ($PARAMS_SMALL)"
echo "----------------------------------------"
profile_run "Legacy" "$PARAMS_SMALL" "-d 12345 67890" "small_legacy"
export DISCOAL_TSKIT_MUTATION_ALGO=node
profile_run "Node-based" "$PARAMS_SMALL" "-d 12345 67890 -ts testing/${OUTDIR}/small_node.trees" "small_node"
export DISCOAL_TSKIT_MUTATION_ALGO=edge
profile_run "Edge-based" "$PARAMS_SMALL" "-d 12345 67890 -ts testing/${OUTDIR}/small_edge.trees" "small_edge"

# Medium scale test
echo -e "\nMedium scale test ($PARAMS_MEDIUM)"
echo "----------------------------------------"
profile_run "Legacy" "$PARAMS_MEDIUM" "-d 12345 67890" "medium_legacy"
export DISCOAL_TSKIT_MUTATION_ALGO=node
profile_run "Node-based" "$PARAMS_MEDIUM" "-d 12345 67890 -ts testing/${OUTDIR}/medium_node.trees" "medium_node"
export DISCOAL_TSKIT_MUTATION_ALGO=edge
profile_run "Edge-based" "$PARAMS_MEDIUM" "-d 12345 67890 -ts testing/${OUTDIR}/medium_edge.trees" "medium_edge"

# Large scale test
echo -e "\nLarge scale test ($PARAMS_LARGE)"
echo "----------------------------------------"
profile_run "Legacy" "$PARAMS_LARGE" "-d 12345 67890" "large_legacy"
export DISCOAL_TSKIT_MUTATION_ALGO=node
profile_run "Node-based" "$PARAMS_LARGE" "-d 12345 67890 -ts testing/${OUTDIR}/large_node.trees" "large_node"
export DISCOAL_TSKIT_MUTATION_ALGO=edge
profile_run "Edge-based" "$PARAMS_LARGE" "-d 12345 67890 -ts testing/${OUTDIR}/large_edge.trees" "large_edge"

# Analyze results
cd testing/
echo -e "\nGenerating performance report..."

# Create analysis script
cat > "${OUTDIR}/analyze_profile.py" << 'EOF'
#!/usr/bin/env python3
import os
import re
import sys

def parse_perf_stat(filename):
    """Parse perf stat output"""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract key metrics
    metrics = {}
    
    # Time
    time_match = re.search(r'(\d+\.\d+)\s+seconds time elapsed', content)
    if time_match:
        metrics['time'] = float(time_match.group(1))
    
    # Instructions
    inst_match = re.search(r'([\d,]+)\s+instructions', content)
    if inst_match:
        metrics['instructions'] = int(inst_match.group(1).replace(',', ''))
    
    # Cycles
    cycles_match = re.search(r'([\d,]+)\s+cycles', content)
    if cycles_match:
        metrics['cycles'] = int(cycles_match.group(1).replace(',', ''))
    
    return metrics

def parse_time_output(filename):
    """Parse GNU time output"""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        content = f.read()
    
    metrics = {}
    
    # User time
    time_match = re.search(r'User time \(seconds\):\s+(\d+\.\d+)', content)
    if time_match:
        metrics['time'] = float(time_match.group(1))
    
    # Maximum resident set size
    mem_match = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', content)
    if mem_match:
        metrics['max_memory_kb'] = int(mem_match.group(1))
    
    return metrics

def get_mutations(ms_file):
    """Extract mutation count from MS file"""
    if not os.path.exists(ms_file):
        return 0
    
    with open(ms_file, 'r') as f:
        for line in f:
            if line.startswith('segsites:'):
                return int(line.split()[1])
    return 0

# Analyze each scale
scales = ['small', 'medium', 'large']
algorithms = ['legacy', 'node', 'edge']

print("\nPERFORMANCE COMPARISON")
print("=" * 70)

for scale in scales:
    print(f"\n{scale.upper()} SCALE:")
    print("-" * 40)
    
    results = {}
    for algo in algorithms:
        prefix = f"{scale}_{algo}"
        
        # Try perf first, then time
        metrics = parse_perf_stat(f"{prefix}_perf.txt")
        if not metrics:
            metrics = parse_time_output(f"{prefix}_time.txt")
        
        if metrics:
            metrics['mutations'] = get_mutations(f"{prefix}.ms")
            results[algo] = metrics
    
    # Print comparison
    if 'legacy' in results:
        legacy_time = results['legacy'].get('time', 1)
        
        for algo in algorithms:
            if algo in results:
                time = results[algo].get('time', 0)
                muts = results[algo].get('mutations', 0)
                speedup = legacy_time / time if time > 0 else 0
                
                print(f"  {algo:10}: {time:6.2f}s ({speedup:4.2f}x) - {muts} mutations", end="")
                
                if 'max_memory_kb' in results[algo]:
                    mem_mb = results[algo]['max_memory_kb'] / 1024
                    print(f" - {mem_mb:.1f} MB", end="")
                
                print()

print("\nNOTES:")
print("- Times include full simulation, not just mutation placement")
print("- Edge-based includes tskit tree sequence construction overhead")
print("- Node-based attempts to maintain RNG compatibility")
EOF

chmod +x "${OUTDIR}/analyze_profile.py"
python3 "${OUTDIR}/analyze_profile.py"

echo -e "\nDetailed results saved in: testing/${OUTDIR}/"
echo "To examine specific metrics, check the *_perf.txt or *_time.txt files"