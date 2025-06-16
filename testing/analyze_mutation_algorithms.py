#!/usr/bin/env python3
"""
Analyze and compare mutation placement algorithms in discoal tskit output
"""
import tskit
import numpy as np
import sys

def analyze_mutations(ts, name):
    """Analyze mutation patterns in a tree sequence"""
    print(f"\n=== {name} ===")
    print(f"Number of mutations: {ts.num_mutations}")
    print(f"Number of sites: {ts.num_sites}")
    print(f"Sequence length: {ts.sequence_length}")
    
    # Analyze mutation distribution across the genome
    if ts.num_sites > 0:
        positions = [site.position for site in ts.sites()]
        print(f"Position range: [{min(positions):.3f}, {max(positions):.3f}]")
        
        # Mutations per edge
        edge_mutations = {}
        for mut in ts.mutations():
            # Find which edge this mutation is on
            tree = ts.at(ts.site(mut.site).position)
            node = mut.node
            parent = tree.parent(node)
            if parent != tskit.NULL:
                # Create edge key
                edge_key = (parent, node)
                edge_mutations[edge_key] = edge_mutations.get(edge_key, 0) + 1
        
        if edge_mutations:
            mut_counts = list(edge_mutations.values())
            print(f"Edges with mutations: {len(edge_mutations)}")
            print(f"Mutations per edge: mean={np.mean(mut_counts):.2f}, "
                  f"std={np.std(mut_counts):.2f}, max={max(mut_counts)}")
        
        # Check for nodes with multiple edges
        node_edge_count = {}
        for edge in ts.edges():
            child = edge.child
            if child not in node_edge_count:
                node_edge_count[child] = []
            node_edge_count[child].append((edge.left, edge.right))
        
        multi_edge_nodes = {n: e for n, e in node_edge_count.items() if len(e) > 1}
        print(f"Nodes spanning multiple edges: {len(multi_edge_nodes)}")
        
        if multi_edge_nodes and edge_mutations:
            # Check mutation distribution on multi-edge nodes
            multi_edge_muts = 0
            for node_id, intervals in multi_edge_nodes.items():
                for tree in ts.trees():
                    parent = tree.parent(node_id)
                    if parent != tskit.NULL:
                        edge_key = (parent, node_id)
                        if edge_key in edge_mutations:
                            multi_edge_muts += edge_mutations[edge_key]
            print(f"Mutations on multi-edge nodes: {multi_edge_muts}")

def compare_algorithms(ts1, ts2, name1="Algorithm 1", name2="Algorithm 2"):
    """Compare mutation patterns between two algorithms"""
    print(f"\n=== Comparison: {name1} vs {name2} ===")
    print(f"Mutation count difference: {ts1.num_mutations} vs {ts2.num_mutations}")
    print(f"Site count difference: {ts1.num_sites} vs {ts2.num_sites}")
    
    # Compare site positions
    if ts1.num_sites > 0 and ts2.num_sites > 0:
        pos1 = sorted([site.position for site in ts1.sites()])
        pos2 = sorted([site.position for site in ts2.sites()])
        
        # Check for exact position matches (unlikely due to different RNG)
        common_positions = set(pos1) & set(pos2)
        print(f"Common mutation positions: {len(common_positions)}")
        
        # Compare position distributions
        print(f"Position statistics:")
        print(f"  {name1}: mean={np.mean(pos1):.3f}, std={np.std(pos1):.3f}")
        print(f"  {name2}: mean={np.mean(pos2):.3f}, std={np.std(pos2):.3f}")

if __name__ == "__main__":
    # Default files to analyze
    files = {
        "Node-based (RNG-compatible)": "test_node_based.trees",
        "Edge-based (msprime-style)": "test_edge_based.trees",
        "Default": "test_default.trees"
    }
    
    # Load and analyze each file
    tree_sequences = {}
    for name, filename in files.items():
        try:
            ts = tskit.load(filename)
            tree_sequences[name] = ts
            analyze_mutations(ts, name)
        except FileNotFoundError:
            print(f"\nWarning: {filename} not found. Run test_tskit_mutations.sh first.")
    
    # Compare algorithms if we have multiple files
    if len(tree_sequences) >= 2:
        ts_list = list(tree_sequences.items())
        if "Node-based (RNG-compatible)" in tree_sequences and "Edge-based (msprime-style)" in tree_sequences:
            compare_algorithms(
                tree_sequences["Node-based (RNG-compatible)"],
                tree_sequences["Edge-based (msprime-style)"],
                "Node-based", "Edge-based"
            )
        
        # Check if default matches node-based
        if "Default" in tree_sequences and "Node-based (RNG-compatible)" in tree_sequences:
            ts_default = tree_sequences["Default"]
            ts_node = tree_sequences["Node-based (RNG-compatible)"]
            if ts_default.num_mutations == ts_node.num_mutations:
                print("\n✓ Default algorithm correctly uses node-based approach")
            else:
                print("\n✗ Warning: Default algorithm doesn't match node-based!")