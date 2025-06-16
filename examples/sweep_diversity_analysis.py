#!/usr/bin/env python3
"""
Sweep Diversity Analysis with discoal and tskit

This script demonstrates how to:
1. Run multiple replicate simulations of selective sweeps using discoal
2. Output tree sequences in tskit format
3. Calculate windowed nucleotide diversity across the genome
4. Plot the average diversity pattern showing the sweep's effect

Author: Generated for discoal examples
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tskit
import os
import sys
from multiprocessing import Pool
import argparse

def run_single_replicate(args):
    """Run a single discoal simulation and return diversity statistics"""
    rep_id, params = args
    
    # Create unique filename for this replicate
    tree_file = f"sweep_rep_{rep_id}.trees"
    
    # Build discoal command
    cmd = [
        params['discoal_path'],
        str(params['n_samples']),
        '1',  # One replicate per run
        str(params['n_sites']),
        '-t', str(params['theta']),
        '-r', str(params['rho']),
        '-ws', str(params['sweep_time']),
        '-a', str(params['alpha']),
        '-x', str(params['sweep_position']),
        '-d', str(params['seed1'] + rep_id), str(params['seed2'] + rep_id),  # Unique seeds
        '-ts', tree_file
    ]
    
    # Run discoal
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in replicate {rep_id}: {result.stderr}")
        return None
    
    # Load tree sequence
    ts = tskit.load(tree_file)
    
    # Calculate windowed diversity
    windows = np.arange(0, ts.sequence_length + params['window_size'], params['window_size'])
    diversity = ts.diversity(windows=windows, mode='branch')
    
    # Clean up tree file
    os.remove(tree_file)
    
    # Return window midpoints and diversity values
    midpoints = (windows[:-1] + windows[1:]) / 2
    return midpoints, diversity

def main():
    parser = argparse.ArgumentParser(description='Analyze diversity patterns in selective sweep simulations')
    parser.add_argument('--replicates', type=int, default=10, help='Number of replicate simulations')
    parser.add_argument('--samples', type=int, default=20, help='Number of haploid samples')
    parser.add_argument('--length', type=int, default=5000000, help='Sequence length in bp')
    parser.add_argument('--window', type=int, default=5000, help='Window size for diversity calculation')
    parser.add_argument('--cores', type=int, default=4, help='Number of parallel cores to use')
    parser.add_argument('--discoal', type=str, default='../discoal', help='Path to discoal executable')
    args = parser.parse_args()
    
    # Simulation parameters
    params = {
        'discoal_path': args.discoal,
        'n_samples': args.samples,
        'n_sites': args.length,
        'theta': 40,  # 4Nu = 40
        'rho': 1000,    # 4Nr = 1000 (higher recombination)
        'sweep_time': 0.00,  # Very recent/ongoing sweep
        'alpha': 200,  # 2Ns = 200 (strong selection)
        'sweep_position': 0.5,  # Sweep in the middle of the sequence
        'window_size': args.window,
        'seed1': 12345,
        'seed2': 67890
    }
    
    print("Selective Sweep Diversity Analysis")
    print("=" * 50)
    print(f"Parameters:")
    print(f"  Samples: {params['n_samples']} haploid")
    print(f"  Sequence length: {params['n_sites']:,} bp")
    print(f"  Theta (4Nu): {params['theta']}")
    print(f"  Rho (4Nr): {params['rho']}")
    print(f"  Selection strength (2Ns): {params['alpha']}")
    print(f"  Sweep position: {params['sweep_position']}")
    print(f"  Sweep time: {params['sweep_time']} (in 4N generations)")
    print(f"  Window size: {params['window_size']:,} bp")
    print(f"  Replicates: {args.replicates}")
    print()
    
    # Check if discoal exists
    if not os.path.exists(params['discoal_path']):
        print(f"Error: discoal not found at {params['discoal_path']}")
        print("Please build discoal or specify the correct path with --discoal")
        sys.exit(1)
    
    print(f"Running {args.replicates} replicate simulations using {args.cores} cores...")
    
    # Run simulations in parallel
    with Pool(processes=args.cores) as pool:
        results = pool.map(run_single_replicate, 
                          [(i, params) for i in range(args.replicates)])
    
    # Filter out failed runs
    successful_results = [r for r in results if r is not None]
    print(f"Successfully completed {len(successful_results)} replicates")
    
    if len(successful_results) == 0:
        print("Error: No successful simulations")
        sys.exit(1)
    
    # Extract diversity values
    all_diversities = np.array([r[1] for r in successful_results])
    window_positions = successful_results[0][0] / 1000  # Convert to kb
    
    # Calculate mean and confidence intervals
    mean_diversity = np.mean(all_diversities, axis=0)
    std_diversity = np.std(all_diversities, axis=0)
    sem_diversity = std_diversity / np.sqrt(len(successful_results))
    
    # Expected diversity without sweep
    expected_diversity = params['theta'] / (1 + params['theta'])
    
    # Create plot
    plt.figure(figsize=(12, 6))
    
    # Plot mean diversity with confidence band
    plt.plot(window_positions, mean_diversity, 'b-', linewidth=2, label='Mean diversity')
    plt.fill_between(window_positions, 
                     mean_diversity - 1.96 * sem_diversity,
                     mean_diversity + 1.96 * sem_diversity,
                     alpha=0.3, color='blue', label='95% CI')
    
    # Add expected diversity line
    plt.axhline(y=expected_diversity, color='red', linestyle='--', 
                label=f'Expected (no sweep): {expected_diversity:.4f}')
    
    # Mark sweep position
    sweep_pos_kb = params['sweep_position'] * params['n_sites'] / 1000
    plt.axvline(x=sweep_pos_kb, color='green', linestyle=':', 
                linewidth=2, label='Sweep position')
    
    # Labels and formatting
    plt.xlabel('Position (kb)', fontsize=12)
    plt.ylabel('branch diversity (π)', fontsize=12)
    plt.title(f'Diversity Pattern Around Selective Sweep\n'
              f'({len(successful_results)} replicates, 2Ns={params["alpha"]}, '
              f'sweep time={params["sweep_time"]}×4N generations ago)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    # Set y-axis limits
    plt.ylim(0, expected_diversity * 1.2)
    
    # Save plot
    output_file = 'sweep_diversity_pattern.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nDiversity statistics:")
    print(f"  Expected diversity (neutral): {expected_diversity:.4f}")
    print(f"  Mean diversity at sweep site: {mean_diversity[len(mean_diversity)//2]:.4f}")
    print(f"  Minimum diversity: {np.min(mean_diversity):.4f}")
    print(f"  Diversity reduction at sweep: {(1 - np.min(mean_diversity)/expected_diversity)*100:.1f}%")
    
    # Calculate sweep width (where diversity is <50% of expected)
    threshold = expected_diversity * 0.5
    below_threshold = window_positions[mean_diversity < threshold]
    if len(below_threshold) > 0:
        sweep_width = (np.max(below_threshold) - np.min(below_threshold))
        print(f"  Sweep width (π < 50% expected): {sweep_width:.1f} kb")
    
    plt.show()

if __name__ == "__main__":
    main()