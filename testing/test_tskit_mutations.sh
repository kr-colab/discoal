#!/bin/bash
# Test script to compare node-based and edge-based mutation algorithms

echo "Building discoal..."
cd ..
make discoal

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo -e "\n=== Testing Node-Based Algorithm (RNG-compatible) ==="
echo "This should match the traditional dropMutations() output"
export DISCOAL_TSKIT_MUTATION_ALGO=node
./discoal 10 1 1000 -t 10 -r 10 -d 123 456 -ts testing/test_node_based.trees
echo "Tree sequence saved to testing/test_node_based.trees"

echo -e "\n=== Testing Edge-Based Algorithm (msprime-style) ==="
echo "This uses independent mutation placement per edge"
export DISCOAL_TSKIT_MUTATION_ALGO=edge
./discoal 10 1 1000 -t 10 -r 10 -d 123 456 -ts testing/test_edge_based.trees
echo "Tree sequence saved to testing/test_edge_based.trees"

echo -e "\n=== Testing Default (should use node-based) ==="
unset DISCOAL_TSKIT_MUTATION_ALGO
./discoal 10 1 1000 -t 10 -r 10 -d 123 456 -ts testing/test_default.trees
echo "Tree sequence saved to testing/test_default.trees"

echo -e "\n=== Comparing outputs ==="
echo "You can use tskit Python API to examine these files:"
echo "  import tskit"
echo "  ts_node = tskit.load('testing/test_node_based.trees')"
echo "  ts_edge = tskit.load('testing/test_edge_based.trees')"
echo "  print(f'Node-based: {ts_node.num_mutations} mutations')"
echo "  print(f'Edge-based: {ts_edge.num_mutations} mutations')"

# Return to testing directory
cd testing/