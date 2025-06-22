#!/bin/bash
# Test script to compare minimal vs full ARG mode

# Parameters
SAMPLE_SIZE=20
NUM_REPS=1
NUM_SITES=10000
RHO=40  # 4Nr = 40, so r = 10/N per site
SEED1=12345
SEED2=67890

echo "Testing full ARG mode (-F) vs minimal mode"
echo "Parameters: n=$SAMPLE_SIZE, sites=$NUM_SITES, rho=$RHO"
echo ""

# Build discoal if needed
if [ ! -f ./discoal ]; then
    echo "Building discoal..."
    make discoal
fi

# Run minimal mode (default)
echo "Running minimal mode..."
./discoal $SAMPLE_SIZE $NUM_REPS $NUM_SITES -r $RHO -d $SEED1 $SEED2 -ts minimal_mode.trees > minimal_mode.out 2>&1

# Run full ARG mode
echo "Running full ARG mode (-F)..."
./discoal $SAMPLE_SIZE $NUM_REPS $NUM_SITES -r $RHO -d $SEED1 $SEED2 -ts full_mode.trees -F > full_mode.out 2>&1

# Analyze the tree sequences
echo ""
echo "=== MINIMAL MODE ==="
python3 measure_tree_roots.py minimal_mode.trees minimal_mode

echo ""
echo "=== FULL ARG MODE ==="
python3 measure_tree_roots.py full_mode.trees full_mode

# Also check tree sequence properties
echo ""
echo "=== Tree sequence comparison ==="
echo "Minimal mode:"
python3 -c "import tskit; ts=tskit.load('minimal_mode.trees'); print(f'  Trees: {ts.num_trees}, Nodes: {ts.num_nodes}, Edges: {ts.num_edges}')"
echo "Full mode:"
python3 -c "import tskit; ts=tskit.load('full_mode.trees'); print(f'  Trees: {ts.num_trees}, Nodes: {ts.num_nodes}, Edges: {ts.num_edges}')"