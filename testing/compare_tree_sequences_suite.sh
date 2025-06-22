#!/bin/bash
#
# Tree sequence comparison suite between discoal and msprime
# Uses proper parameter conversion from discoal's coalescent units to msprime's physical units
#

set -e  # Exit on error

# Configuration
REPLICATES=10
SAMPLE_SIZE=20
SEQUENCE_LENGTH=100000
OUTPUT_DIR="tree_sequence_comparison_$(date +%Y%m%d_%H%M%S)"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Tree Sequence Comparison Suite"
echo "========================================"
echo "Comparing discoal and msprime tree sequences"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Build discoal if needed
cd ..
if [ ! -f "discoal" ]; then
    echo "Building discoal..."
    make discoal
fi
cd testing

# Python script to generate msprime tree sequences with proper parameter conversion
cat > "$OUTPUT_DIR/generate_msprime_trees.py" << 'EOF'
#!/usr/bin/env python3
"""
Generate msprime tree sequences with parameters matching discoal conventions.

Key parameter conversions:
- discoal uses haploid samples, msprime uses diploid individuals
- theta = 4*Ne*mu*L (scaled mutation rate)
- rho = 4*Ne*r*L (scaled recombination rate)
- Ne = 0.5 for msprime with haploid samples (following mspms convention)
"""

import msprime
import sys
import numpy as np

def generate_msprime_trees(n_samples, sequence_length, theta, rho, n_replicates, output_prefix, seed=None):
    """Generate tree sequences with msprime using discoal parameter conventions."""
    
    if seed is not None:
        np.random.seed(seed)
    
    # Key insight from mspms: Use Ne=0.5 for haploid samples
    Ne = 0.5
    
    # Convert scaled parameters to per-bp rates
    mutation_rate = theta / (4 * Ne * sequence_length)
    recombination_rate = rho / (4 * Ne * sequence_length) if rho > 0 else 0
    
    print(f"msprime parameters:")
    print(f"  Samples: {n_samples} haploid")
    print(f"  Ne: {Ne}")
    print(f"  Sequence length: {sequence_length}")
    print(f"  Mutation rate: {mutation_rate:.2e} per bp")
    print(f"  Recombination rate: {recombination_rate:.2e} per bp")
    
    for i in range(n_replicates):
        # Simulate ancestry - use haploid samples directly
        ts = msprime.sim_ancestry(
            samples=n_samples,
            ploidy=1,  # Haploid
            sequence_length=sequence_length,
            population_size=Ne,
            recombination_rate=recombination_rate,
            random_seed=seed + i if seed else None
        )
        
        # Add mutations
        ts = msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            random_seed=(seed + 1000 + i) if seed else None
        )
        
        # Save tree sequence
        filename = f"{output_prefix}_rep{i+1}.ts"
        ts.dump(filename)
        print(f"  Saved replicate {i+1}: {ts.num_trees} trees, {ts.num_sites} sites")

if __name__ == "__main__":
    # Parse arguments
    if len(sys.argv) != 8:
        print("Usage: generate_msprime_trees.py n_samples seq_length theta rho n_reps output_prefix seed")
        sys.exit(1)
    
    n_samples = int(sys.argv[1])
    sequence_length = int(sys.argv[2])
    theta = float(sys.argv[3])
    rho = float(sys.argv[4])
    n_replicates = int(sys.argv[5])
    output_prefix = sys.argv[6]
    seed = int(sys.argv[7])
    
    generate_msprime_trees(n_samples, sequence_length, theta, rho, n_replicates, output_prefix, seed)
EOF

chmod +x "$OUTPUT_DIR/generate_msprime_trees.py"

# Function to run comparison for a specific scenario
run_comparison() {
    local scenario_name=$1
    local theta=$2
    local rho=$3
    local extra_args=$4
    local description=$5
    
    echo "----------------------------------------------------------------------"
    echo "Scenario: $scenario_name"
    echo "Description: $description"
    echo "Parameters: theta=$theta, rho=$rho"
    
    # Create subdirectory for this scenario
    mkdir -p "$OUTPUT_DIR/$scenario_name"
    
    # Generate discoal tree sequences
    echo "Generating discoal tree sequences..."
    ../discoal $SAMPLE_SIZE $REPLICATES $SEQUENCE_LENGTH \
        -t $theta -r $rho $extra_args \
        -d 12345 67890 \
        -ts "$OUTPUT_DIR/$scenario_name/discoal.trees" \
        > "$OUTPUT_DIR/$scenario_name/discoal.ms" 2> "$OUTPUT_DIR/$scenario_name/discoal.err"
    
    # Generate msprime tree sequences
    echo "Generating msprime tree sequences..."
    /home/adkern/miniforge3/envs/discoal_dev/bin/python "$OUTPUT_DIR/generate_msprime_trees.py" \
        $SAMPLE_SIZE $SEQUENCE_LENGTH $theta $rho $REPLICATES \
        "$OUTPUT_DIR/$scenario_name/msprime" 12345
    
    # Compare using our tree sequence comparison tool
    echo "Comparing tree sequences..."
    /home/adkern/miniforge3/envs/discoal_dev/bin/python ../compare_tree_sequences.py compare \
        "$OUTPUT_DIR/$scenario_name/msprime" \
        "$OUTPUT_DIR/$scenario_name/discoal" \
        $REPLICATES > "$OUTPUT_DIR/$scenario_name/comparison.txt"
    
    # Show summary
    echo ""
    tail -n 20 "$OUTPUT_DIR/$scenario_name/comparison.txt" | head -n 15
    echo ""
}

# Test scenarios (starting with just a few models)

# 1. Neutral model without recombination
run_comparison "neutral_no_recomb" 20 0 "" \
    "Neutral evolution without recombination"

# 2. Neutral model with moderate recombination
run_comparison "neutral_with_recomb" 20 50 "" \
    "Neutral evolution with moderate recombination"

# 3. Neutral model with high recombination
run_comparison "neutral_high_recomb" 20 200 "" \
    "Neutral evolution with high recombination"

# Summary
echo "======================================================================"
echo "Summary of comparisons"
echo "======================================================================"
echo ""

for scenario in neutral_no_recomb neutral_with_recomb neutral_high_recomb; do
    if [ -f "$OUTPUT_DIR/$scenario/comparison.txt" ]; then
        echo "Scenario: $scenario"
        grep "p-value" "$OUTPUT_DIR/$scenario/comparison.txt" | grep -E "trees|height|diversity" | head -3
        echo ""
    fi
done

echo "Full results saved in: $OUTPUT_DIR"
echo ""
echo "To examine tree sequences in detail:"
echo "  import tskit"
echo "  ts_discoal = tskit.load('$OUTPUT_DIR/SCENARIO/discoal_rep1.trees')"
echo "  ts_msprime = tskit.load('$OUTPUT_DIR/SCENARIO/msprime_rep1.ts')"