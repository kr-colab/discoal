#!/usr/bin/env python3
"""
Tree sequence comparison utility extracted from notebooks/tree_peak.ipynb
Provides comprehensive analysis and comparison of tree sequences from discoal and msprime.
"""

import tskit
import numpy as np
import pandas as pd
from scipy import stats
import sys
import os

def analyze_tree_sequences(output_prefix, reps, source="discoal"):
    """
    Analyze multiple replicate tree sequences and compute summary statistics.
    
    Parameters:
    - output_prefix: Base name for output files (e.g., "discoal_output" or "msprime_output")
    - reps: Number of replicates to analyze
    - source: Either "discoal" or "msprime" to handle file naming conventions
    """
    
    all_stats = []
    
    for i in range(1, reps + 1):
        # Load tree sequence
        if source == "discoal":
            tree_file_path = f"{output_prefix}_rep{i}.trees"
        else:  # msprime
            tree_file_path = f"{output_prefix}_rep{i}.ts"
        
        try:
            ts = tskit.load(tree_file_path)
        except FileNotFoundError:
            print(f"Warning: Could not find {tree_file_path}, skipping...")
            continue
        
        # Calculate tree heights
        tree_heights = []
        tree_spans = []
        num_roots_per_tree = []
        
        for tree in ts.trees():
            # Get all roots (there may be multiple in case of non-coalescence)
            roots = tree.roots
            num_roots_per_tree.append(len(roots))
            
            # For each root, calculate the tree height
            if len(roots) > 0:
                max_height = max(tree.time(root) for root in roots)
                tree_heights.append(max_height)
                tree_spans.append(tree.span)
        
        # Calculate diversity metrics
        diversity = ts.diversity(mode="branch") if ts.num_samples > 1 else 0
        
        # Calculate segregating sites
        segregating_sites = ts.segregating_sites(mode="branch")
        
        # Store statistics for this replicate
        stats = {
            "replicate": i,
            "num_samples": ts.num_samples,
            "num_trees": ts.num_trees,
            "num_sites": ts.num_sites,
            "num_mutations": ts.num_mutations,
            "num_nodes": ts.num_nodes,
            "num_edges": ts.num_edges,
            "sequence_length": ts.sequence_length,
            "mean_tree_height": np.average(tree_heights, weights=tree_spans) if tree_heights else 0,
            "max_tree_height": max(tree_heights) if tree_heights else 0,
            "min_tree_height": min(tree_heights) if tree_heights else 0,
            "diversity": diversity,
            "segregating_sites": segregating_sites,
            "mean_roots_per_tree": np.mean(num_roots_per_tree),
            "max_roots_per_tree": max(num_roots_per_tree) if num_roots_per_tree else 0,
            "fully_coalesced": all(n == 1 for n in num_roots_per_tree)
        }
        
        all_stats.append(stats)
    
    # Convert to DataFrame for easy analysis
    df = pd.DataFrame(all_stats)
    
    if len(df) == 0:
        print(f"No valid tree sequences found for {source}")
        return df
    
    # Print summary across all replicates
    print(f"\n{'='*60}")
    print(f"Summary statistics for {source} ({len(df)} replicates analyzed)")
    print(f"{'='*60}\n")
    
    # Basic counts
    print("Basic metrics (mean ± std):")
    print(f"  Samples: {df['num_samples'].iloc[0] if len(df) > 0 else 'N/A'}")  # Should be constant
    print(f"  Trees: {df['num_trees'].mean():.1f} ± {df['num_trees'].std():.1f}")
    print(f"  Sites: {df['num_sites'].mean():.1f} ± {df['num_sites'].std():.1f}")
    print(f"  Mutations: {df['num_mutations'].mean():.1f} ± {df['num_mutations'].std():.1f}")
    print(f"  Nodes: {df['num_nodes'].mean():.1f} ± {df['num_nodes'].std():.1f}")
    print(f"  Edges: {df['num_edges'].mean():.1f} ± {df['num_edges'].std():.1f}")
    
    # Tree structure
    print("\nTree structure:")
    print(f"  Mean tree height: {df['mean_tree_height'].mean():.4f} ± {df['mean_tree_height'].std():.4f}")
    print(f"  Max tree height: {df['max_tree_height'].mean():.4f} ± {df['max_tree_height'].std():.4f}")
    print(f"  Mean roots per tree: {df['mean_roots_per_tree'].mean():.3f} ± {df['mean_roots_per_tree'].std():.3f}")
    print(f"  Fully coalesced replicates: {df['fully_coalesced'].sum()}/{len(df)} ({100*df['fully_coalesced'].mean():.1f}%)")
    
    # Genetic diversity
    print("\nGenetic diversity:")
    print(f"  Pairwise diversity (π): {df['diversity'].mean():.6f} ± {df['diversity'].std():.6f}")
    print(f"  Segregating sites (S): {df['segregating_sites'].mean():.6f} ± {df['segregating_sites'].std():.6f}")
    
    return df

def compare_simulators(df_msprime, df_discoal):
    """Compare summary statistics between msprime and discoal outputs."""
    
    print(f"\n{'='*60}")
    print("Comparison between msprime and discoal")
    print(f"{'='*60}\n")
    
    # Metrics to compare
    metrics = [
        ("num_trees", "Number of trees"),
        ("num_sites", "Number of sites"), 
        ("num_mutations", "Number of mutations"),
        ("mean_tree_height", "Mean tree height"),
        ("diversity", "Pairwise diversity"),
        ("segregating_sites", "Segregating sites"),
        ("mean_roots_per_tree", "Mean roots per tree")
    ]
    
    for metric, label in metrics:
        msp_mean = df_msprime[metric].mean()
        dis_mean = df_discoal[metric].mean()
        
        # Perform t-test
        t_stat, p_value = stats.ttest_ind(df_msprime[metric], df_discoal[metric])
        
        print(f"{label}:")
        print(f"  msprime: {msp_mean:.6f}")
        print(f"  discoal: {dis_mean:.6f}")
        print(f"  Difference: {abs(msp_mean - dis_mean):.6f} ({100*abs(msp_mean - dis_mean)/msp_mean if msp_mean != 0 else 0:.1f}%)")
        print(f"  t-test p-value: {p_value:.4f} {'*' if p_value < 0.05 else ''}")
        print()

def summarize_single_ts(ts, label):
    """Summarize a single tree sequence."""
    print(f"\nSummary statistics for {label}:")
    print(f"  Number of samples: {ts.num_samples}")
    print(f"  Number of trees: {ts.num_trees}")
    print(f"  Number of sites: {ts.num_sites}")
    print(f"  Number of mutations: {ts.num_mutations}")
    print(f"  Number of nodes: {ts.num_nodes}")
    print(f"  Number of edges: {ts.num_edges}")
    print(f"  Sequence length: {ts.sequence_length}")
    
    # Calculate tree heights
    tree_heights = []
    tree_spans = []
    for tree in ts.trees():
        if tree.num_roots == 1:
            tree_heights.append(tree.time(tree.root))
            tree_spans.append(tree.span)
    
    if tree_heights:
        mean_height = np.average(tree_heights, weights=tree_spans)
        print(f"  Mean tree height: {mean_height:.4f}")
    
    # Diversity
    if ts.num_samples > 1:
        print(f"  Mean pairwise diversity (π): {ts.diversity(mode='branch'):.6f}")
        print(f"  Segregating sites (S): {ts.segregating_sites(mode='branch'):.6f}")

def main():
    """Main function for command-line usage."""
    if len(sys.argv) < 2:
        print("Usage: python compare_tree_sequences.py <command> [args]")
        print("\nCommands:")
        print("  single <file.trees>        - Summarize a single tree sequence")
        print("  analyze <prefix> <reps>    - Analyze multiple replicates")
        print("  compare <prefix1> <prefix2> <reps> - Compare two sets of replicates")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "single":
        if len(sys.argv) < 3:
            print("Error: Please provide a tree sequence file")
            sys.exit(1)
        ts = tskit.load(sys.argv[2])
        summarize_single_ts(ts, sys.argv[2])
    
    elif command == "analyze":
        if len(sys.argv) < 4:
            print("Error: Please provide prefix and number of replicates")
            sys.exit(1)
        prefix = sys.argv[2]
        reps = int(sys.argv[3])
        source = sys.argv[4] if len(sys.argv) > 4 else "discoal"
        analyze_tree_sequences(prefix, reps, source)
    
    elif command == "compare":
        if len(sys.argv) < 5:
            print("Error: Please provide two prefixes and number of replicates")
            sys.exit(1)
        prefix1 = sys.argv[2]
        prefix2 = sys.argv[3]
        reps = int(sys.argv[4])
        
        # Analyze both sets
        df1 = analyze_tree_sequences(prefix1, reps, "msprime")
        df2 = analyze_tree_sequences(prefix2, reps, "discoal")
        
        # Compare if both have data
        if len(df1) > 0 and len(df2) > 0:
            compare_simulators(df1, df2)
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

if __name__ == "__main__":
    main()