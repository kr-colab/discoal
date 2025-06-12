#!/usr/bin/env python3
"""
msprime comparison utilities for discoal testing

This script provides functions to run equivalent simulations in msprime
and compare the results with discoal output.
"""

import msprime
import numpy as np
import sys
import argparse
from scipy import stats
import re


def parse_ms_output(ms_output):
    """Parse ms-format output to extract segregating sites and haplotypes."""
    lines = ms_output.strip().split('\n')
    
    results = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('//'):
            i += 1
            # Parse segsites line
            if i < len(lines) and lines[i].startswith('segsites:'):
                n_sites = int(lines[i].split()[1])
                i += 1
                
                # Parse positions line if there are segregating sites
                if n_sites > 0 and i < len(lines) and lines[i].startswith('positions:'):
                    positions = list(map(float, lines[i].split()[1:]))
                    i += 1
                    
                    # Parse haplotypes
                    haplotypes = []
                    while i < len(lines) and not lines[i].startswith('//') and lines[i].strip():
                        haplotypes.append(lines[i].strip())
                        i += 1
                    
                    results.append({
                        'n_sites': n_sites,
                        'positions': positions,
                        'haplotypes': haplotypes
                    })
                else:
                    # No segregating sites
                    results.append({
                        'n_sites': 0,
                        'positions': [],
                        'haplotypes': []
                    })
        else:
            i += 1
    
    return results


def run_neutral_comparison(n_samples, n_sites, theta, rho=0, n_replicates=100, seed=None, Ne=0.5):
    """
    Run neutral coalescent simulations in msprime and return summary statistics.
    
    Parameters matching discoal conventions:
    - n_samples: number of haploid samples (chromosomes)
    - n_sites: number of sites (L)
    - theta: population mutation rate (4*Ne*mu*L)
    - rho: population recombination rate (4*Ne*r*(L-1))
    - n_replicates: number of replicates
    - seed: random seed
    - Ne: effective population size (default 1000)
    
    Note: discoal uses haploid samples, while msprime typically uses diploid individuals
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Convert from discoal's scaled parameters to per-generation rates
    # In discoal, theta = 4*Ne*mu*L where mu is per-basepair per-generation mutation rate
    # and L is the sequence length
    # So for msprime, we need mu = theta / (4*Ne*L)
    mutation_rate = theta / (4 * Ne * n_sites)
    
    # Similarly, rho = 4*Ne*r*L where r is per-basepair per-generation recombination rate
    if rho > 0:
        recombination_rate = rho / (4 * Ne * n_sites)
    else:
        recombination_rate = 0
    
    seg_sites = []
    pi_values = []
    tajimas_d = []
    
    for rep in range(n_replicates):
        # Simulate ancestry
        # discoal uses haploid samples, so we need to convert to diploid individuals for msprime
        # n_samples haploid chromosomes = n_samples/2 diploid individuals
        ts = msprime.sim_ancestry(
            samples=n_samples // 2,  # Convert haploid samples to diploid individuals
            sequence_length=n_sites,
            recombination_rate=recombination_rate,
            population_size=Ne,
            ploidy=2,  # Use diploid (default)
            random_seed=np.random.randint(1, 2**31)
        )
        
        # Add mutations
        mts = msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            random_seed=np.random.randint(1, 2**31)
        )
        
        # Calculate summary statistics
        n_seg_sites = mts.num_sites
        seg_sites.append(n_seg_sites)
        
        # Calculate nucleotide diversity (pi)
        if n_seg_sites > 0:
            div = mts.diversity()
            pi_values.append(div)
            
            # Calculate Tajima's D
            # This is a simplified calculation
            n = n_samples
            a1 = sum(1/i for i in range(1, n))
            S = n_seg_sites
            theta_w = S / a1
            theta_pi = div * n_sites
            
            if theta_w > 0:
                # Tajima's D calculation
                a2 = sum(1/i**2 for i in range(1, n))
                b1 = (n + 1) / (3 * (n - 1))
                b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
                c1 = b1 - 1/a1
                c2 = b2 - (n + 2)/(a1 * n) + a2/a1**2
                e1 = c1/a1
                e2 = c2/(a1**2 + a2)
                
                var_d = e1 * S + e2 * S * (S - 1)
                if var_d > 0:
                    d = (theta_pi - theta_w) / np.sqrt(var_d)
                    tajimas_d.append(d)
        else:
            pi_values.append(0)
    
    return {
        'seg_sites': seg_sites,
        'pi': pi_values,
        'tajimas_d': tajimas_d,
        'mean_seg_sites': np.mean(seg_sites),
        'std_seg_sites': np.std(seg_sites),
        'mean_pi': np.mean(pi_values),
        'std_pi': np.std(pi_values)
    }


def run_sweep_comparison(n_samples, n_sites, theta, rho, alpha, sweep_position=0.5, 
                        tau=0.01, n_replicates=100, seed=None, Ne=0.5, refsize=10000):
    """
    Run selective sweep simulations in msprime.
    
    Parameters matching discoal conventions:
    - alpha: selection strength (2*Ne*s)
    - tau: time since sweep in units of 4*Ne generations
    - sweep_position: position of selected site (0-1)
    
    Note: msprime's sweep models work differently than discoal's, so this provides
    an approximation using structured coalescent with selection.
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Convert parameters following discoal conventions exactly as in reference script
    # Initial calculation
    s = alpha / (2.0 * refsize)
    s = s * 2  # msprime fitness model adjustment
    
    # Handle time scaling exactly as in reference script
    original_refsize = refsize  # Keep original for frequency calculations
    if tau > 0:
        # When sweep_mod_time > 0, rescale to Ne=0.25 for consistent time scales
        refsize = 0.25
        s = alpha / refsize  # Recalculate s with new refsize
        Ne_sim = refsize
        sweep_time = tau  # tau is already in coalescent units
    else:
        # When tau=0, use the original refsize
        Ne_sim = refsize
        sweep_time = 0
    
    # Convert rates using the simulation Ne (refsize)
    recomb_rate = rho / (4 * refsize * (n_sites - 1))
    mutation_rate = theta / (4 * refsize * n_sites)
    
    seg_sites = []
    
    for rep in range(n_replicates):
        # Create a sweep model following reference script exactly
        # From reference script: start = 1/(2*original_refsize), end = 1 - 1/(2*original_refsize)
        # Use original refsize for frequency bounds to ensure valid [0,1] range
        start_freq = 1.0 / (2 * original_refsize)
        end_freq = 1.0 - (1.0 / (2 * original_refsize))
        
        sweep_model = msprime.SweepGenicSelection(
            position=np.floor(sweep_position * n_sites),
            start_frequency=start_freq,
            end_frequency=end_freq,
            s=s,
            dt=1e-6  # Small time step
        )
        
        # Simulate with sweep
        # If tau=0, the sweep just completed so we only run the sweep model
        # If tau>0, we run standard coalescent for sweep_time, then the sweep
        if sweep_time > 0:
            models = [
                msprime.StandardCoalescent(duration=sweep_time),
                sweep_model,
                msprime.StandardCoalescent()
            ]
        else:
            models = [
                sweep_model,
                msprime.StandardCoalescent()
                ]
            
        # Use diploid samples for consistency with neutral model
        ts = msprime.sim_ancestry(
            samples=n_samples // 2,  # Convert haploid samples to diploid individuals
            sequence_length=n_sites,
            recombination_rate=recomb_rate,
            population_size=refsize,  # Use refsize as population size
            ploidy=2,  # Use diploid
            model=models,
            discrete_genome=False,  # As in reference script
            random_seed=np.random.randint(1, 2**31)
        )
        
        # Add mutations
        mts = msprime.sim_mutations(
            ts,
            rate=mutation_rate,
            random_seed=np.random.randint(1, 2**31)
        )
        
        seg_sites.append(mts.num_sites)
    
    return {
        'seg_sites': seg_sites,
        'mean_seg_sites': np.mean(seg_sites),
        'std_seg_sites': np.std(seg_sites)
    }


def compare_distributions(discoal_stats, msprime_stats, test_name):
    """
    Compare distributions from discoal and msprime using statistical tests.
    """
    results = {}
    
    # Kolmogorov-Smirnov test
    ks_stat, ks_pval = stats.ks_2samp(discoal_stats, msprime_stats)
    results['ks_statistic'] = ks_stat
    results['ks_pvalue'] = ks_pval
    
    # Mann-Whitney U test
    mw_stat, mw_pval = stats.mannwhitneyu(discoal_stats, msprime_stats, alternative='two-sided')
    results['mw_statistic'] = mw_stat
    results['mw_pvalue'] = mw_pval
    
    # Summary statistics
    results['discoal_mean'] = np.mean(discoal_stats)
    results['discoal_std'] = np.std(discoal_stats)
    results['msprime_mean'] = np.mean(msprime_stats)
    results['msprime_std'] = np.std(msprime_stats)
    
    # Check if distributions are statistically similar
    results['similar'] = ks_pval > 0.05 and mw_pval > 0.05
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Run msprime simulations for comparison with discoal')
    parser.add_argument('--mode', choices=['neutral', 'recombination', 'sweep'], 
                       default='neutral', help='Simulation mode')
    parser.add_argument('--samples', type=int, default=20, help='Number of haploid samples (chromosomes)')
    parser.add_argument('--sites', type=int, default=10000, help='Number of sites')
    parser.add_argument('--theta', type=float, default=10.0, help='Population mutation rate (4*Ne*mu*L)')
    parser.add_argument('--rho', type=float, default=0.0, help='Population recombination rate (4*Ne*r*L)')
    parser.add_argument('--alpha', type=float, default=1000.0, help='Selection coefficient (2*Ne*s)')
    parser.add_argument('--tau', type=float, default=0.01, help='Time since sweep (in 4Ne generations)')
    parser.add_argument('--Ne', type=float, default=1, help='Effective population size')
    parser.add_argument('--refsize', type=float, default=10000, help='Reference population size for sweeps')
    parser.add_argument('--replicates', type=int, default=100, help='Number of replicates')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    parser.add_argument('--output-stats', action='store_true', 
                       help='Output summary statistics instead of raw simulations')
    
    args = parser.parse_args()
    
    if args.mode == 'neutral' or args.mode == 'recombination':
        stats = run_neutral_comparison(
            n_samples=args.samples,
            n_sites=args.sites,
            theta=args.theta,
            rho=args.rho,
            n_replicates=args.replicates,
            seed=args.seed,
            Ne=args.Ne
        )
        
        if args.output_stats:
            print(f"Mean segregating sites: {stats['mean_seg_sites']:.2f} (SD: {stats['std_seg_sites']:.2f})")
            print(f"Mean pi: {stats['mean_pi']:.6f} (SD: {stats['std_pi']:.6f})")
        else:
            # Output segregating sites counts for comparison
            for s in stats['seg_sites']:
                print(s)
    
    elif args.mode == 'sweep':
        stats = run_sweep_comparison(
            n_samples=args.samples,
            n_sites=args.sites,
            theta=args.theta,
            rho=args.rho,
            alpha=args.alpha,
            tau=args.tau,
            n_replicates=args.replicates,
            seed=args.seed,
            Ne=args.Ne,
            refsize=args.refsize
        )
        
        if args.output_stats:
            print(f"Mean segregating sites: {stats['mean_seg_sites']:.2f} (SD: {stats['std_seg_sites']:.2f})")
        else:
            for s in stats['seg_sites']:
                print(s)


if __name__ == '__main__':
    main()