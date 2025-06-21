#!/bin/bash
#
# Profile mutation algorithms to identify performance hotspots
# Uses perf record/report for detailed analysis
#

set -e

echo "========================================"
echo "Mutation Algorithm Hotspot Profiling"
echo "========================================"

# Check for perf
if ! command -v perf &> /dev/null; then
    echo "Error: perf not found. Install linux-tools-generic or equivalent."
    exit 1
fi

# Create output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="mutation_hotspots_${TIMESTAMP}"
mkdir -p "$OUTDIR"

# Build with debug symbols
echo "Building discoal with debug symbols..."
cd ..
make clean > /dev/null 2>&1
make CFLAGS="-g -O2 -march=native" discoal > /dev/null 2>&1

# Test parameters for hotspot analysis
PARAMS="100 20 50000 -t 200 -r 100 -d 12345 67890"

echo -e "\nProfiling parameters: $PARAMS"
echo "This will generate ~200 mutations per replicate"

# Profile legacy dropMutations
echo -e "\n1. Profiling legacy dropMutations()..."
perf record -g -o "testing/${OUTDIR}/legacy.perf.data" \
    ./discoal $PARAMS > "testing/${OUTDIR}/legacy.ms" 2>&1
echo "   Generated $(grep -c "segsites:" "testing/${OUTDIR}/legacy.ms") replicates"

# Profile edge-based
echo -e "\n2. Profiling edge-based algorithm..."
export DISCOAL_TSKIT_MUTATION_ALGO=edge
perf record -g -o "testing/${OUTDIR}/edge.perf.data" \
    ./discoal $PARAMS -ts "testing/${OUTDIR}/edge.trees" > "testing/${OUTDIR}/edge.ms" 2>&1
echo "   Generated $(grep -c "segsites:" "testing/${OUTDIR}/edge.ms") replicates"

# Profile node-based
echo -e "\n3. Profiling node-based algorithm..."
export DISCOAL_TSKIT_MUTATION_ALGO=node
perf record -g -o "testing/${OUTDIR}/node.perf.data" \
    ./discoal $PARAMS -ts "testing/${OUTDIR}/node.trees" > "testing/${OUTDIR}/node.ms" 2>&1
echo "   Generated $(grep -c "segsites:" "testing/${OUTDIR}/node.ms") replicates"

cd testing/

# Generate reports
echo -e "\nGenerating performance reports..."

# Top functions for each algorithm
echo -e "\n=== TOP FUNCTIONS BY ALGORITHM ===" > "${OUTDIR}/hotspots_summary.txt"

echo -e "\nLEGACY dropMutations():" >> "${OUTDIR}/hotspots_summary.txt"
perf report --stdio -i "${OUTDIR}/legacy.perf.data" 2>/dev/null | \
    grep -A 20 "Overhead" | grep -E "dropMutations|addMutation|isAncestral|ignpoi|genunf" \
    >> "${OUTDIR}/hotspots_summary.txt" || true

echo -e "\nEDGE-BASED:" >> "${OUTDIR}/hotspots_summary.txt"
perf report --stdio -i "${OUTDIR}/edge.perf.data" 2>/dev/null | \
    grep -A 20 "Overhead" | grep -E "tskit_place_mutations_edge|tskit_add_site|tskit_add_mutation|ignpoi|genunf" \
    >> "${OUTDIR}/hotspots_summary.txt" || true

echo -e "\nNODE-BASED:" >> "${OUTDIR}/hotspots_summary.txt"
perf report --stdio -i "${OUTDIR}/node.perf.data" 2>/dev/null | \
    grep -A 20 "Overhead" | grep -E "tskit_place_mutations_node|tskit_add_site|tskit_add_mutation|ignpoi|genunf" \
    >> "${OUTDIR}/hotspots_summary.txt" || true

# Generate full reports
perf report --stdio -i "${OUTDIR}/legacy.perf.data" > "${OUTDIR}/legacy_full_report.txt" 2>&1
perf report --stdio -i "${OUTDIR}/edge.perf.data" > "${OUTDIR}/edge_full_report.txt" 2>&1
perf report --stdio -i "${OUTDIR}/node.perf.data" > "${OUTDIR}/node_full_report.txt" 2>&1

# Create comparison script
cat > "${OUTDIR}/compare_hotspots.py" << 'EOF'
#!/usr/bin/env python3
import re
import sys

def parse_perf_report(filename):
    """Extract function percentages from perf report"""
    functions = {}
    
    try:
        with open(filename, 'r') as f:
            in_overhead = False
            for line in f:
                if "Overhead" in line and "Symbol" in line:
                    in_overhead = True
                    continue
                
                if in_overhead and line.strip() == "":
                    break
                
                if in_overhead:
                    match = re.match(r'\s*(\d+\.\d+)%.*\[.\]\s+(\S+)', line)
                    if match:
                        pct = float(match.group(1))
                        func = match.group(2)
                        # Group by base function name
                        base_func = func.split('.')[0]
                        if base_func in functions:
                            functions[base_func] += pct
                        else:
                            functions[base_func] = pct
    except:
        pass
    
    return functions

# Parse reports
print("\nMUTATION ALGORITHM PERFORMANCE HOTSPOTS")
print("=" * 60)

algorithms = [
    ("Legacy", "legacy_full_report.txt"),
    ("Edge-based", "edge_full_report.txt"),
    ("Node-based", "node_full_report.txt")
]

all_functions = {}
for name, file in algorithms:
    funcs = parse_perf_report(file)
    if funcs:
        all_functions[name] = funcs
        print(f"\n{name} Top Functions:")
        print("-" * 40)
        sorted_funcs = sorted(funcs.items(), key=lambda x: x[1], reverse=True)[:10]
        for func, pct in sorted_funcs:
            print(f"  {func:30} {pct:5.1f}%")

# Find mutation-specific functions
print("\n\nMUTATION-SPECIFIC FUNCTIONS:")
print("-" * 60)
mutation_funcs = [
    'dropMutations', 'addMutation', 'isAncestralHere',
    'tskit_place_mutations_edge', 'tskit_place_mutations_node',
    'tskit_add_site', 'tskit_add_mutation', 'ignpoi', 'genunf'
]

for func in mutation_funcs:
    print(f"\n{func}:")
    for name, funcs in all_functions.items():
        if func in funcs:
            print(f"  {name:15} {funcs[func]:5.1f}%")
EOF

chmod +x "${OUTDIR}/compare_hotspots.py"
python3 "${OUTDIR}/compare_hotspots.py"

cat "${OUTDIR}/hotspots_summary.txt"

echo -e "\nDetailed profiling results saved in: testing/${OUTDIR}/"
echo "View interactive reports with:"
echo "  perf report -i testing/${OUTDIR}/legacy.perf.data"
echo "  perf report -i testing/${OUTDIR}/edge.perf.data"
echo "  perf report -i testing/${OUTDIR}/node.perf.data"