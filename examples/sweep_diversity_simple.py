#!/usr/bin/env python3
"""
Simple Sweep Diversity Analysis

A simplified version that runs fewer replicates for quick demonstration.
Shows the characteristic valley of diversity around a selective sweep.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tskit
import os

# Parameters
n_replicates = 10
n_samples = 20
sequence_length = 5000000  # 200 kb
window_size = 5000  # 5 kb windows
theta = 40  # 4Nu
rho = 1000    # 4Nr  
alpha = 200  # 2Ns (strong selection)
sweep_position = 0.5  # Middle of sequence
sweep_time = 0.00  # Recent sweep

print("Running selective sweep simulations...")
print(f"Parameters: {n_samples} samples, {sequence_length:,} bp, {n_replicates} replicates")

# Store diversity results
all_diversities = []

for rep in range(n_replicates):
    # Run discoal
    tree_file = f"temp_sweep_{rep}.trees"
    cmd = [
        '../discoal', str(n_samples), '1', str(sequence_length),
        '-t', str(theta), '-r', str(rho),
        '-ws', str(sweep_time), '-a', str(alpha), '-x', str(sweep_position),
        '-d', str(12345 + rep), str(67890 + rep),
        '-ts', tree_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in replicate {rep}: {result.stderr}")
        continue
    
    # Load and analyze
    ts = tskit.load(tree_file)
    
    # Calculate windowed diversity
    windows = np.arange(0, ts.sequence_length + window_size, window_size)
    diversity = ts.diversity(windows=windows, mode='branch')
    all_diversities.append(diversity)
    
    # Clean up
    os.remove(tree_file)
    
    print(f"  Completed replicate {rep + 1}/{n_replicates}")

# Calculate mean diversity
all_diversities = np.array(all_diversities)
mean_diversity = np.mean(all_diversities, axis=0)
window_positions = np.arange(window_size/2, sequence_length, window_size) / 1000  # kb

# Expected diversity
expected_diversity = theta / (theta + 1) 

# Plot
plt.figure(figsize=(10, 6))
plt.plot(window_positions, mean_diversity, 'b-', linewidth=2, label='Mean diversity')
plt.axhline(y=expected_diversity, color='red', linestyle='--', 
            label=f'Expected (neutral): {expected_diversity:.3f}')
plt.axvline(x=sweep_position * sequence_length / 1000, color='green', 
            linestyle=':', linewidth=2, label='Sweep position')

plt.xlabel('Position (kb)')
plt.ylabel('Nucleotide diversity (π)')
plt.title(f'Diversity Pattern Around Selective Sweep (2Ns={alpha})')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig('simple_sweep_diversity.png', dpi=150)
print("\nPlot saved to: simple_sweep_diversity.png")

# Print statistics
min_div = np.min(mean_diversity)
print(f"\nDiversity at sweep: {min_div:.4f} (reduction: {(1-min_div/expected_diversity)*100:.1f}%)")

plt.show()