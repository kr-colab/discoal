#!/usr/bin/env python3
"""
Profile discoal vs mspms performance across different recombination rates.
Generates publication-quality figures comparing execution time and memory usage.

USAGE:
    # First activate the discoal_dev conda environment:
    conda activate discoal_dev
    
    # Then run the script:
    python profile_discoal_vs_mspms.py
"""

import subprocess
import time
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import pandas as pd
import argparse

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.minor.size'] = 2

def run_simulation(program, n_samples, seq_length, recomb_rate, mut_rate, replicate, output_dir):
    """Run a single simulation and measure time/memory."""
    
    # Construct command based on program
    if program == 'discoal':
        cmd = ['./discoal', str(n_samples), '1', str(seq_length), 
               '-r', str(recomb_rate), '-t', str(mut_rate)]
        cmd_str = ' '.join(cmd)
    else:  # mspms
        cmd = ['mspms', str(n_samples), '1', 
               '-r', str(recomb_rate), str(seq_length), '-t', str(mut_rate)]
        cmd_str = ' '.join(cmd)
    
    # Output file for this run
    output_file = os.path.join(output_dir, f"{program}_r{recomb_rate}_rep{replicate}.out")
    
    # Use GNU time for memory tracking if available
    if subprocess.run(['which', 'gtime'], capture_output=True).returncode == 0:
        time_cmd = 'gtime'
    elif subprocess.run(['which', '/usr/bin/time'], capture_output=True).returncode == 0:
        time_cmd = '/usr/bin/time'
    else:
        # Fallback to Python timing (no memory info)
        start_time = time.time()
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            elapsed_time = time.time() - start_time
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            return {
                'program': program,
                'recomb_rate': recomb_rate,
                'replicate': replicate,
                'time': elapsed_time,
                'memory': None,
                'status': 'success',
                'cmd': cmd_str
            }
        except subprocess.CalledProcessError as e:
            return {
                'program': program,
                'recomb_rate': recomb_rate,
                'replicate': replicate,
                'time': None,
                'memory': None,
                'status': 'failed',
                'error': str(e),
                'cmd': cmd_str
            }
    
    # Run with time command for memory tracking
    full_cmd = [time_cmd, '-v'] + cmd
    
    try:
        start_time = time.time()
        result = subprocess.run(full_cmd, capture_output=True, text=True, check=True)
        elapsed_time = time.time() - start_time
        
        # Parse GNU time output
        stderr_lines = result.stderr.split('\n')
        max_memory = None
        user_time = None
        sys_time = None
        
        for line in stderr_lines:
            if 'Maximum resident set size' in line:
                # Extract memory in KB
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        max_memory = float(parts[1].strip()) / 1024  # Convert to MB
                    except:
                        pass
            elif 'User time (seconds):' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        user_time = float(parts[1].strip())
                    except:
                        pass
            elif 'System time (seconds):' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        sys_time = float(parts[1].strip())
                    except:
                        pass
        
        # Use CPU time if available, otherwise wall time
        if user_time is not None and sys_time is not None:
            cpu_time = user_time + sys_time
        else:
            cpu_time = elapsed_time
        
        # Save output
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        
        return {
            'program': program,
            'recomb_rate': recomb_rate,
            'replicate': replicate,
            'time': cpu_time,
            'wall_time': elapsed_time,
            'memory': max_memory,
            'status': 'success',
            'cmd': cmd_str
        }
    
    except subprocess.CalledProcessError as e:
        return {
            'program': program,
            'recomb_rate': recomb_rate,
            'replicate': replicate,
            'time': None,
            'memory': None,
            'status': 'failed',
            'error': str(e),
            'cmd': cmd_str
        }

def run_profiling(n_samples=100, seq_length=100000, recomb_rates=None, mut_rate=10, 
                  n_replicates=5, n_workers=10, output_dir='profiling_results'):
    """Run profiling comparison between discoal and mspms."""
    
    if recomb_rates is None:
        # Default range of recombination rates
        recomb_rates = [0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000]
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = os.path.join(output_dir, f'run_{timestamp}')
    os.makedirs(run_dir, exist_ok=True)
    
    # Save parameters
    params = {
        'n_samples': n_samples,
        'seq_length': seq_length,
        'recomb_rates': recomb_rates,
        'mut_rate': mut_rate,
        'n_replicates': n_replicates,
        'timestamp': timestamp
    }
    with open(os.path.join(run_dir, 'parameters.json'), 'w') as f:
        json.dump(params, f, indent=2)
    
    print(f"Running profiling comparison")
    print(f"Sample size: {n_samples}")
    print(f"Sequence length: {seq_length:,} bp")
    print(f"Mutation rate (theta): {mut_rate}")
    print(f"Recombination rates: {recomb_rates}")
    print(f"Replicates per condition: {n_replicates}")
    print(f"Output directory: {run_dir}")
    print(f"Using {n_workers} parallel workers")
    print()
    
    # Check if programs exist
    if not os.path.exists('./discoal'):
        print("ERROR: ./discoal not found!")
        return None, None
    
    if subprocess.run(['which', 'mspms'], capture_output=True).returncode != 0:
        print("ERROR: mspms not found in PATH!")
        return None, None
    
    # Prepare all tasks
    tasks = []
    for program in ['discoal', 'mspms']:
        for recomb_rate in recomb_rates:
            for replicate in range(n_replicates):
                tasks.append((program, n_samples, seq_length, recomb_rate, 
                            mut_rate, replicate, run_dir))
    
    # Run tasks in parallel
    results = []
    total_tasks = len(tasks)
    completed = 0
    
    print(f"Running {total_tasks} simulations...")
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(run_simulation, *task): task 
            for task in tasks
        }
        
        # Process results as they complete
        for future in as_completed(future_to_task):
            result = future.result()
            results.append(result)
            completed += 1
            
            # Progress update
            if completed % 10 == 0 or completed == total_tasks:
                print(f"Progress: {completed}/{total_tasks} ({100*completed/total_tasks:.1f}%)")
    
    # Convert to DataFrame for analysis
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(run_dir, 'raw_results.csv'), index=False)
    
    # Check for failures
    failed = df[df['status'] == 'failed']
    if len(failed) > 0:
        print(f"\nWARNING: {len(failed)} simulations failed!")
        print(failed[['program', 'recomb_rate', 'replicate', 'error']])
    
    # Filter successful runs
    df_success = df[df['status'] == 'success'].copy()
    
    return df_success, run_dir

def generate_figures(df, output_dir):
    """Generate publication-quality figures."""
    
    # Set up the plotting style
    sns.set_style("whitegrid")
    sns.set_context("paper")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle('Performance Comparison: discoal vs mspms', fontsize=14, y=0.98)
    
    # Prepare summary statistics
    summary = df.groupby(['program', 'recomb_rate']).agg({
        'time': ['mean', 'std', 'count'],
        'memory': ['mean', 'std', 'count']
    }).reset_index()
    
    # Flatten column names
    summary.columns = ['program', 'recomb_rate', 'time_mean', 'time_std', 'time_count',
                      'memory_mean', 'memory_std', 'memory_count']
    
    # Colors for programs
    colors = {'discoal': '#1f77b4', 'mspms': '#ff7f0e'}
    
    # 1. Time vs Recombination Rate (linear scale)
    ax = axes[0, 0]
    for program in ['discoal', 'mspms']:
        data = summary[summary['program'] == program]
        ax.errorbar(data['recomb_rate'], data['time_mean'], 
                   yerr=data['time_std'], 
                   label=program, color=colors[program],
                   marker='o', capsize=5, capthick=1, linewidth=1.5)
    ax.set_xlabel('Recombination rate (ρ)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Execution Time vs Recombination Rate')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Time vs Recombination Rate (log scale)
    ax = axes[0, 1]
    for program in ['discoal', 'mspms']:
        data = summary[summary['program'] == program]
        # Add small offset to handle r=0
        r_plot = data['recomb_rate'].values
        r_plot[r_plot == 0] = 0.1
        ax.errorbar(r_plot, data['time_mean'], 
                   yerr=data['time_std'], 
                   label=program, color=colors[program],
                   marker='o', capsize=5, capthick=1, linewidth=1.5)
    ax.set_xlabel('Recombination rate (ρ)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Execution Time (log scale)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Memory vs Recombination Rate
    ax = axes[1, 0]
    # Check if memory data is available
    if df['memory'].notna().any():
        for program in ['discoal', 'mspms']:
            data = summary[summary['program'] == program]
            ax.errorbar(data['recomb_rate'], data['memory_mean'], 
                       yerr=data['memory_std'], 
                       label=program, color=colors[program],
                       marker='o', capsize=5, capthick=1, linewidth=1.5)
        ax.set_xlabel('Recombination rate (ρ)')
        ax.set_ylabel('Memory (MB)')
        ax.set_title('Peak Memory Usage vs Recombination Rate')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'Memory data not available\n(GNU time not found)', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel('Recombination rate (ρ)')
        ax.set_ylabel('Memory (MB)')
        ax.set_title('Peak Memory Usage vs Recombination Rate')
    
    # 4. Speedup ratio
    ax = axes[1, 1]
    # Calculate speedup (mspms time / discoal time)
    discoal_data = summary[summary['program'] == 'discoal'].set_index('recomb_rate')
    mspms_data = summary[summary['program'] == 'mspms'].set_index('recomb_rate')
    
    speedup = mspms_data['time_mean'] / discoal_data['time_mean']
    ax.plot(speedup.index, speedup.values, 'ko-', linewidth=1.5, markersize=6)
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Recombination rate (ρ)')
    ax.set_ylabel('Speedup (mspms time / discoal time)')
    ax.set_title('discoal Speedup Relative to mspms')
    ax.grid(True, alpha=0.3)
    
    # Add text annotation for average speedup
    avg_speedup = speedup.mean()
    ax.text(0.02, 0.98, f'Average speedup: {avg_speedup:.1f}x', 
            transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save figure
    fig_path = os.path.join(output_dir, 'performance_comparison.pdf')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    fig_path_png = os.path.join(output_dir, 'performance_comparison.png')
    plt.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    print(f"\nFigures saved to:")
    print(f"  - {fig_path}")
    print(f"  - {fig_path_png}")
    
    # Generate detailed statistics table
    stats_file = os.path.join(output_dir, 'performance_statistics.txt')
    with open(stats_file, 'w') as f:
        f.write("Performance Comparison: discoal vs mspms\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Parameters:\n")
        params_file = os.path.join(output_dir, 'parameters.json')
        if os.path.exists(params_file):
            with open(params_file, 'r') as pf:
                params = json.load(pf)
                for key, value in params.items():
                    f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        f.write("Summary Statistics:\n")
        f.write("-" * 50 + "\n")
        for _, row in summary.iterrows():
            f.write(f"\n{row['program']} - ρ={row['recomb_rate']}:\n")
            f.write(f"  Time: {row['time_mean']:.3f} ± {row['time_std']:.3f} seconds\n")
            if pd.notna(row['memory_mean']):
                f.write(f"  Memory: {row['memory_mean']:.1f} ± {row['memory_std']:.1f} MB\n")
        
        f.write("\n\nSpeedup Analysis:\n")
        f.write("-" * 50 + "\n")
        for r in speedup.index:
            f.write(f"ρ={r}: {speedup[r]:.2f}x faster\n")
        f.write(f"\nAverage speedup: {avg_speedup:.2f}x\n")
    
    print(f"  - {stats_file}")
    
    return fig

def main():
    parser = argparse.ArgumentParser(
        description='Profile discoal vs mspms performance across recombination rates')
    parser.add_argument('-n', '--samples', type=int, default=100,
                        help='Number of samples (default: 100)')
    parser.add_argument('-l', '--length', type=int, default=100000,
                        help='Sequence length in bp (default: 100000)')
    parser.add_argument('-t', '--theta', type=float, default=10,
                        help='Mutation rate theta (default: 10)')
    parser.add_argument('--reps', type=int, default=5,
                        help='Number of replicates per condition (default: 5)')
    parser.add_argument('--workers', type=int, default=12,
                        help='Number of parallel workers (default: 4)')
    parser.add_argument('-r', '--recomb-rates', nargs='+', type=float,
                        help='Recombination rates to test (default: 0 1 5 10 20 50 100 200 500 1000)')
    parser.add_argument('-o', '--output', type=str, default='profiling_results',
                        help='Output directory (default: profiling_results)')
    
    args = parser.parse_args()
    
    # Run profiling
    df, output_dir = run_profiling(
        n_samples=args.samples,
        seq_length=args.length,
        recomb_rates=args.recomb_rates,
        mut_rate=args.theta,
        n_replicates=args.reps,
        n_workers=args.workers,
        output_dir=args.output
    )
    
    if df is not None and len(df) > 0:
        # Generate figures
        generate_figures(df, output_dir)
        print("\nProfiling complete!")
    else:
        print("\nProfiling failed - no results to plot")
        sys.exit(1)

if __name__ == '__main__':
    main()