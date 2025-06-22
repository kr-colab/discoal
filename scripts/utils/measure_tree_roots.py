#!/usr/bin/env python3
"""
Measure the number of roots in each tree across the chromosome
for tree sequences produced by discoal in minimal vs full ARG mode.
"""

import tskit
import sys
import numpy as np
import matplotlib.pyplot as plt

def analyze_tree_roots(ts_file, output_prefix):
    """Analyze and plot the number of roots per tree."""
    ts = tskit.load(ts_file)
    
    # Collect data about each tree
    tree_data = []
    for tree in ts.trees():
        num_roots = len(tree.roots)
        tree_data.append({
            'interval': tree.interval,
            'num_roots': num_roots,
            'span': tree.span
        })
    
    # Print summary statistics
    print(f"\nTree sequence file: {ts_file}")
    print(f"Number of trees: {ts.num_trees}")
    print(f"Number of nodes: {ts.num_nodes}")
    print(f"Number of edges: {ts.num_edges}")
    print(f"Number of sites: {ts.num_sites}")
    print(f"Number of mutations: {ts.num_mutations}")
    
    # Root count statistics
    root_counts = [t['num_roots'] for t in tree_data]
    print(f"\nRoot count statistics:")
    print(f"  Min roots: {min(root_counts)}")
    print(f"  Max roots: {max(root_counts)}")
    print(f"  Mean roots: {np.mean(root_counts):.2f}")
    print(f"  Median roots: {np.median(root_counts):.0f}")
    print(f"  Trees with >1 root: {sum(1 for r in root_counts if r > 1)} ({100*sum(1 for r in root_counts if r > 1)/len(root_counts):.1f}%)")
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Root count across chromosome
    positions = []
    root_values = []
    for t in tree_data:
        # Add points at start and end of each interval
        positions.extend([t['interval'][0], t['interval'][1]])
        root_values.extend([t['num_roots'], t['num_roots']])
    
    ax1.plot(positions[:-1], root_values[:-1], 'b-', linewidth=0.5)
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Number of Roots')
    ax1.set_title(f'Number of Tree Roots Across Chromosome - {ts_file}')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Histogram of root counts
    ax2.hist(root_counts, bins=range(1, max(root_counts)+2), alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Number of Roots')
    ax2.set_ylabel('Number of Trees')
    ax2.set_title('Distribution of Root Counts')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_roots.png", dpi=150)
    print(f"\nPlot saved to: {output_prefix}_roots.png")
    
    # Also output detailed data for first few multi-root trees
    print("\nFirst 5 trees with multiple roots:")
    count = 0
    for i, t in enumerate(tree_data):
        if t['num_roots'] > 1:
            print(f"  Tree {i}: interval {t['interval']}, {t['num_roots']} roots")
            count += 1
            if count >= 5:
                break

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python measure_tree_roots.py <tree_sequence.trees> [output_prefix]")
        sys.exit(1)
    
    ts_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else ts_file.rsplit('.', 1)[0]
    
    analyze_tree_roots(ts_file, output_prefix)