#!/bin/bash
#
# Benchmark how mutation algorithms scale with different parameters
#

set -e

echo "========================================"
echo "Mutation Algorithm Scaling Benchmark"
echo "========================================"

# Create output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="mutation_scaling_${TIMESTAMP}"
mkdir -p "$OUTDIR"

# Build discoal
echo "Building discoal..."
cd ..
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1

# Function to run benchmark
benchmark() {
    local label=$1
    local algo=$2
    local params=$3
    local output=$4
    
    if [ "$algo" == "legacy" ]; then
        unset DISCOAL_TSKIT_MUTATION_ALGO
        extra_args=""
    else
        export DISCOAL_TSKIT_MUTATION_ALGO=$algo
        extra_args="-ts testing/${OUTDIR}/temp.trees"
    fi
    
    # Run 3 times and take average
    total_time=0
    for i in 1 2 3; do
        start=$(date +%s.%N)
        ./discoal $params -d 12345 67890 $extra_args > /dev/null 2>&1
        end=$(date +%s.%N)
        elapsed=$(echo "$end - $start" | bc)
        total_time=$(echo "$total_time + $elapsed" | bc)
    done
    
    avg_time=$(echo "scale=3; $total_time / 3" | bc)
    echo "$label,$avg_time" >> "testing/${OUTDIR}/$output"
    echo "  $label: ${avg_time}s"
}

# Test 1: Scaling with mutation rate (theta)
echo -e "\n1. Scaling with mutation rate (constant n=50, L=10000)"
echo "-----------------------------------------------------"

for algo in legacy node edge; do
    echo -e "\n$algo algorithm:"
    echo "theta,time" > "testing/${OUTDIR}/theta_scaling_${algo}.csv"
    
    for theta in 10 50 100 200 400 800; do
        benchmark "theta=$theta" "$algo" "50 5 10000 -t $theta" "theta_scaling_${algo}.csv"
    done
done

# Test 2: Scaling with sample size
echo -e "\n2. Scaling with sample size (constant theta=100, L=10000)"
echo "--------------------------------------------------------"

for algo in legacy node edge; do
    echo -e "\n$algo algorithm:"
    echo "samples,time" > "testing/${OUTDIR}/sample_scaling_${algo}.csv"
    
    for n in 10 20 50 100 200; do
        benchmark "n=$n" "$algo" "$n 5 10000 -t 100" "sample_scaling_${algo}.csv"
    done
done

# Test 3: Scaling with sequence length
echo -e "\n3. Scaling with sequence length (constant n=50, theta=100)"
echo "---------------------------------------------------------"

for algo in legacy node edge; do
    echo -e "\n$algo algorithm:"
    echo "length,time" > "testing/${OUTDIR}/length_scaling_${algo}.csv"
    
    for L in 1000 5000 10000 50000 100000; do
        benchmark "L=$L" "$algo" "50 5 $L -t 100" "length_scaling_${algo}.csv"
    done
done

# Test 4: Scaling with recombination rate
echo -e "\n4. Scaling with recombination rate (n=50, theta=100, L=10000)"
echo "------------------------------------------------------------"

for algo in legacy node edge; do
    echo -e "\n$algo algorithm:"
    echo "rho,time" > "testing/${OUTDIR}/recomb_scaling_${algo}.csv"
    
    for rho in 0 10 50 100 200 500; do
        benchmark "rho=$rho" "$algo" "50 5 10000 -t 100 -r $rho" "recomb_scaling_${algo}.csv"
    done
done

cd testing/

# Create analysis script
cat > "${OUTDIR}/plot_scaling.py" << 'EOF'
#!/usr/bin/env python3
import os
import sys
import csv

def read_csv(filename):
    """Read benchmark CSV file"""
    if not os.path.exists(filename):
        return None, None
    
    x_vals = []
    y_vals = []
    
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            x_vals.append(float(row[0]))
            y_vals.append(float(row[1]))
    
    return x_vals, y_vals

def print_comparison(title, param_name, legacy_file, node_file, edge_file):
    """Print scaling comparison"""
    print(f"\n{title}")
    print("=" * 60)
    
    legacy_x, legacy_y = read_csv(legacy_file)
    node_x, node_y = read_csv(node_file)
    edge_x, edge_y = read_csv(edge_file)
    
    if not legacy_x:
        return
    
    print(f"\n{param_name:>10} | {'Legacy':>10} | {'Node':>10} | {'Edge':>10} | {'Node/Leg':>10} | {'Edge/Leg':>10}")
    print("-" * 75)
    
    for i in range(len(legacy_x)):
        param = legacy_x[i]
        leg_time = legacy_y[i]
        node_time = node_y[i] if i < len(node_y) else 0
        edge_time = edge_y[i] if i < len(edge_y) else 0
        
        node_ratio = node_time / leg_time if leg_time > 0 else 0
        edge_ratio = edge_time / leg_time if leg_time > 0 else 0
        
        print(f"{param:10.0f} | {leg_time:10.3f} | {node_time:10.3f} | {edge_time:10.3f} | {node_ratio:10.2f}x | {edge_ratio:10.2f}x")

# Analyze results
print("\nMUTATION ALGORITHM SCALING ANALYSIS")

print_comparison(
    "Scaling with Mutation Rate (theta)",
    "theta",
    "theta_scaling_legacy.csv",
    "theta_scaling_node.csv", 
    "theta_scaling_edge.csv"
)

print_comparison(
    "Scaling with Sample Size",
    "n",
    "sample_scaling_legacy.csv",
    "sample_scaling_node.csv",
    "sample_scaling_edge.csv"
)

print_comparison(
    "Scaling with Sequence Length",
    "L",
    "length_scaling_legacy.csv",
    "length_scaling_node.csv",
    "length_scaling_edge.csv"
)

print_comparison(
    "Scaling with Recombination Rate",
    "rho",
    "recomb_scaling_legacy.csv",
    "recomb_scaling_node.csv",
    "recomb_scaling_edge.csv"
)

print("\n\nKEY OBSERVATIONS:")
print("-" * 60)
print("1. Both tskit algorithms include tree sequence construction overhead")
print("2. Legacy algorithm benefits from simpler data structures")
print("3. Edge-based may be more efficient for high recombination")
print("4. Node-based maintains conceptual compatibility with legacy")
EOF

chmod +x "${OUTDIR}/plot_scaling.py"
python3 "${OUTDIR}/plot_scaling.py"

echo -e "\nScaling benchmark results saved in: testing/${OUTDIR}/"
echo "CSV files contain raw timing data for further analysis"