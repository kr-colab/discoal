#!/bin/bash
# Test script to verify full ARG produces multi-root trees

# Parameters
SAMPLE_SIZE=10
NUM_REPS=1
NUM_SITES=1000
RHO=100  # Very high recombination
SEED1=12345
SEED2=67890

echo "Testing full ARG with high recombination to see multi-root trees"
echo "Parameters: n=$SAMPLE_SIZE, sites=$NUM_SITES, rho=$RHO"
echo ""

# Run full ARG mode with high recombination
./discoal $SAMPLE_SIZE $NUM_REPS $NUM_SITES -r $RHO -d $SEED1 $SEED2 -ts high_recomb_full.trees -F

# Analyze
python3 -c "
import tskit
ts = tskit.load('high_recomb_full.trees')
print(f'Trees: {ts.num_trees}')
print(f'Nodes: {ts.num_nodes}')
print(f'Edges: {ts.num_edges}')

# Check for multi-root trees
multi_root_count = 0
max_roots = 0
for tree in ts.trees():
    num_roots = len(tree.roots)
    if num_roots > 1:
        multi_root_count += 1
    max_roots = max(max_roots, num_roots)

print(f'\\nTrees with multiple roots: {multi_root_count} / {ts.num_trees}')
print(f'Maximum roots in any tree: {max_roots}')

# Show some recombination nodes
print('\\nChecking for recombination nodes:')
for i in range(min(20, ts.num_nodes)):
    node = ts.node(i)
    if node.flags == 0:  # Non-sample
        # Count how many trees this node appears in
        tree_count = 0
        child_count = 0
        for tree in ts.trees():
            if tree.parent(i) != tskit.NULL:
                tree_count += 1
            for c in tree.children(i):
                child_count += 1
        if child_count > 0:
            print(f'  Node {i}: appears as parent in {tree_count} trees, time={node.time:.4f}')
"