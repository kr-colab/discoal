#!/bin/bash
#
# Comprehensive validation suite for tskit mutation placement algorithms
# Tests both node-based (RNG-compatible) and edge-based algorithms
#

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test parameters
REPLICATES=100
SAMPLE_SIZE=20
SEQUENCE_LENGTH=10000
THETA=20
RECOMB=20

echo "========================================"
echo "TSKIT Mutation Algorithm Validation Suite"
echo "========================================"
echo ""

# Create output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="tskit_mutation_validation_${TIMESTAMP}"
mkdir -p "$OUTDIR"

# Build discoal
echo "Building discoal..."
cd ..
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Build successful${NC}"

# Function to run test
run_test() {
    local test_name=$1
    local algo=$2
    local extra_args=$3
    
    echo -n "  Testing $test_name... "
    
    if [ "$algo" == "node" ]; then
        export DISCOAL_TSKIT_MUTATION_ALGO=node
    elif [ "$algo" == "edge" ]; then
        export DISCOAL_TSKIT_MUTATION_ALGO=edge
    else
        unset DISCOAL_TSKIT_MUTATION_ALGO
    fi
    
    # Run with fixed seed for reproducibility
    ./discoal $SAMPLE_SIZE $REPLICATES $SEQUENCE_LENGTH -t $THETA -r $RECOMB \
        -d 12345 67890 $extra_args -ts "testing/$OUTDIR/${test_name}_${algo}.trees" \
        > "testing/$OUTDIR/${test_name}_${algo}.ms" 2> "testing/$OUTDIR/${test_name}_${algo}.err"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓${NC}"
        return 0
    else
        echo -e "${RED}✗${NC}"
        return 1
    fi
}

# Test 1: Basic comparison
echo ""
echo "Test 1: Basic mutation placement"
echo "--------------------------------"
run_test "basic" "node" ""
run_test "basic" "edge" ""
run_test "basic" "default" ""

# Test 2: No recombination
echo ""
echo "Test 2: No recombination (single tree)"
echo "--------------------------------------"
run_test "no_recomb" "node" "-r 0"
run_test "no_recomb" "edge" "-r 0"

# Test 3: High recombination
echo ""
echo "Test 3: High recombination rate"
echo "-------------------------------"
run_test "high_recomb" "node" "-r 100"
run_test "high_recomb" "edge" "-r 100"

# # Test 4: Multiple populations
echo ""
echo "Test 4: Multiple populations"
echo "----------------------------"
run_test "multipop" "node" "-I 2 10 10 0.1"
run_test "multipop" "edge" "-I 2 10 10 0.1"

# # Test 5: Selection
echo ""
echo "Test 5: Selection (sweep)"
echo "------------------------"
run_test "selection" "node" "-ws 0.1 -a 100 -x 0.5"
run_test "selection" "edge" "-ws 0.1 -a 100 -x 0.5"

# Compare outputs
echo ""
echo "Analyzing results..."
echo "===================="

cd testing/

# Python script to analyze results (without tskit dependency)
cat > "$OUTDIR/analyze_results.py" << 'EOF'
#!/usr/bin/env python3
import sys
import os
import re

def analyze_ms_file(filename):
    """Extract segregating sites from MS file"""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract all segsites values
    segsites = re.findall(r'segsites:\s*(\d+)', content)
    if segsites:
        values = [int(s) for s in segsites]
        total = sum(values)
        avg = total / len(values) if values else 0
        return {
            "total_sites": total,
            "average_sites": avg,
            "num_replicates": len(values),
            "min_sites": min(values) if values else 0,
            "max_sites": max(values) if values else 0
        }
    return None

def analyze_test(test_name, outdir):
    """Analyze a specific test comparing node vs edge algorithms"""
    results = {}
    
    for algo in ["node", "edge", "default"]:
        ms_file = f"{outdir}/{test_name}_{algo}.ms"
        ms_results = analyze_ms_file(ms_file)
        if ms_results:
            results[algo] = ms_results
    
    return results

def print_comparison(test_name, results):
    """Print comparison results"""
    print(f"\n{test_name.upper()} TEST RESULTS:")
    print("-" * 50)
    
    if not results:
        print("  No results found")
        return
    
    # Check if default matches node-based
    if "default" in results and "node" in results:
        if results["default"]["total_sites"] == results["node"]["total_sites"]:
            print("  ✓ Default correctly uses node-based algorithm")
        else:
            print("  ✗ WARNING: Default doesn't match node-based!")
    
    # Compare node vs edge
    if "node" in results and "edge" in results:
        node_avg = results["node"]["average_sites"]
        edge_avg = results["edge"]["average_sites"]
        node_total = results["node"]["total_sites"]
        edge_total = results["edge"]["total_sites"]
        
        diff_pct = (edge_avg - node_avg) / max(node_avg, 1) * 100
        
        print(f"  Node-based: {results['node']['num_replicates']} replicates")
        print(f"    Average mutations: {node_avg:.1f}")
        print(f"    Range: [{results['node']['min_sites']}, {results['node']['max_sites']}]")
        print(f"  Edge-based: {results['edge']['num_replicates']} replicates")  
        print(f"    Average mutations: {edge_avg:.1f}")
        print(f"    Range: [{results['edge']['min_sites']}, {results['edge']['max_sites']}]")
        print(f"  Difference: {diff_pct:+.1f}%")
        
        # Statistical check - they should be similar on average
        if abs(diff_pct) > 50:
            print("  ⚠ Large difference detected - may indicate an issue")

# Main analysis
if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    
    tests = ["basic", "no_recomb", "high_recomb", "multipop", "selection"]
    
    print("\nTSKIT MUTATION ALGORITHM COMPARISON")
    print("=" * 50)
    
    all_results = {}
    for test in tests:
        results = analyze_test(test, outdir)
        if results:
            all_results[test] = results
            print_comparison(test, results)
    
    # Overall summary
    print("\n\nOVERALL SUMMARY:")
    print("=" * 50)
    
    total_tests = len(all_results)
    defaults_match = sum(1 for r in all_results.values() 
                        if "default" in r and "node" in r 
                        and r["default"]["total_sites"] == r["node"]["total_sites"])
    
    print(f"Total tests run: {total_tests}")
    print(f"Default algorithm matches node-based: {defaults_match}/{total_tests}")
    
    # Check for systematic differences
    diffs = []
    for test, results in all_results.items():
        if "node" in results and "edge" in results:
            node_avg = results["node"]["average_sites"]
            edge_avg = results["edge"]["average_sites"]
            if node_avg > 0:
                diffs.append((edge_avg - node_avg) / node_avg)
    
    if diffs:
        mean_diff = sum(diffs) / len(diffs) * 100
        print(f"\nMean difference (edge vs node): {mean_diff:.1f}%")
        if abs(mean_diff) < 5:
            print("✓ Algorithms produce statistically similar results")
        else:
            print("⚠ Algorithms show systematic differences")
EOF

chmod +x "$OUTDIR/analyze_results.py"

# Run analysis
python3 "$OUTDIR/analyze_results.py" "$OUTDIR"

# Create summary report
echo ""
echo "Summary report saved to: testing/$OUTDIR/summary.txt"
python3 "$OUTDIR/analyze_results.py" "$OUTDIR" > "$OUTDIR/summary.txt" 2>&1

echo ""
echo -e "${GREEN}Validation complete!${NC}"
echo "Results saved in: testing/$OUTDIR/"
echo ""
echo "Analysis based on MS output files (segregating sites counts)."
echo "Tree sequence files saved for future analysis with tskit."