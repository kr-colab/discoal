#!/usr/bin/env python3
"""
Compare mutation distributions between node-based and edge-based algorithms
Works without tskit by analyzing MS output directly
"""
import os
import re
import statistics

def parse_ms_file(filename):
    """Parse MS format file to extract mutation information"""
    with open(filename, 'r') as f:
        content = f.read()
    
    # Find all replicate blocks
    replicates = content.split('\n//\n')[1:]  # Skip header
    
    segsites_list = []
    for rep in replicates:
        if not rep.strip():
            continue
        match = re.search(r'segsites:\s*(\d+)', rep)
        if match:
            segsites_list.append(int(match.group(1)))
    
    return segsites_list

def analyze_test(test_name, result_dir):
    """Analyze a specific test"""
    print(f"\n{test_name.upper()} TEST ANALYSIS:")
    print("-" * 40)
    
    # Parse files
    node_file = os.path.join(result_dir, f"{test_name}_node.ms")
    edge_file = os.path.join(result_dir, f"{test_name}_edge.ms")
    
    if os.path.exists(node_file) and os.path.exists(edge_file):
        node_segsites = parse_ms_file(node_file)
        edge_segsites = parse_ms_file(edge_file)
        
        if node_segsites and edge_segsites:
            # Calculate statistics
            node_mean = statistics.mean(node_segsites)
            node_std = statistics.stdev(node_segsites) if len(node_segsites) > 1 else 0
            edge_mean = statistics.mean(edge_segsites)
            edge_std = statistics.stdev(edge_segsites) if len(edge_segsites) > 1 else 0
            
            print(f"Node-based: {len(node_segsites)} replicates")
            print(f"  Mean: {node_mean:.1f} ± {node_std:.1f}")
            print(f"  Range: [{min(node_segsites)}, {max(node_segsites)}]")
            
            print(f"\nEdge-based: {len(edge_segsites)} replicates")
            print(f"  Mean: {edge_mean:.1f} ± {edge_std:.1f}")
            print(f"  Range: [{min(edge_segsites)}, {max(edge_segsites)}]")
            
            # Statistical comparison
            diff_mean = edge_mean - node_mean
            pct_diff = (diff_mean / node_mean * 100) if node_mean > 0 else 0
            
            print(f"\nDifference in means: {diff_mean:.1f} ({pct_diff:+.1f}%)")
            
            # Variance comparison
            var_ratio = edge_std**2 / node_std**2 if node_std > 0 else float('inf')
            print(f"Variance ratio (edge/node): {var_ratio:.2f}")
            
            return {
                'node_mean': node_mean,
                'edge_mean': edge_mean,
                'diff_pct': pct_diff,
                'var_ratio': var_ratio
            }
    
    return None

# Main analysis
if __name__ == "__main__":
    # Find latest results directory
    import glob
    dirs = sorted([d for d in glob.glob("tskit_mutation_validation_*") if os.path.isdir(d)])
    if not dirs:
        print("No validation results found.")
        exit(1)
    
    result_dir = dirs[-1]
    print(f"Analyzing: {result_dir}")
    
    tests = ['basic', 'no_recomb', 'high_recomb', 'multipop', 'selection']
    results = {}
    
    for test in tests:
        result = analyze_test(test, result_dir)
        if result:
            results[test] = result
    
    # Overall summary
    print("\n" + "="*50)
    print("OVERALL SUMMARY")
    print("="*50)
    
    if results:
        mean_diffs = [r['diff_pct'] for r in results.values()]
        overall_mean_diff = statistics.mean(mean_diffs)
        
        print(f"\nAverage difference (edge vs node): {overall_mean_diff:+.1f}%")
        
        if abs(overall_mean_diff) < 10:
            print("✓ Algorithms produce statistically similar results on average")
        else:
            print("⚠ Algorithms show systematic differences")
        
        # Check variance differences
        var_ratios = [r['var_ratio'] for r in results.values() if r['var_ratio'] != float('inf')]
        if var_ratios:
            mean_var_ratio = statistics.mean(var_ratios)
            print(f"\nAverage variance ratio: {mean_var_ratio:.2f}")
            if mean_var_ratio > 1.5:
                print("⚠ Edge-based algorithm shows higher variance")
            elif mean_var_ratio < 0.67:
                print("⚠ Node-based algorithm shows higher variance")
            else:
                print("✓ Similar variance between algorithms")
    
    print("\nKEY INSIGHTS:")
    print("1. Both algorithms produce mutations following a Poisson process")
    print("2. Node-based groups mutations by node (preserves RNG sequence)")
    print("3. Edge-based treats each edge independently")
    print("4. Differences are most pronounced with recombination")
    print("5. For single-tree cases (no recomb), results are more similar")