#!/usr/bin/env python3
"""
Statistical comparison of summary statistics between discoal versions.
Compares means, medians, standard deviations, and performs statistical tests.
"""

import sys
import numpy as np
from scipy import stats
import argparse
from collections import defaultdict
import pandas as pd

def parse_nicestats_output(filename):
    """Parse nicestats output file and extract statistics for each replicate."""
    stats_dict = defaultdict(list)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
        if not lines:
            return stats_dict
            
        # First line should be headers
        headers = lines[0].strip().split('\t')
        
        # Parse data lines
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
                
            values = line.split('\t')
            if len(values) != len(headers):
                continue
                
            # Store each statistic
            for i, header in enumerate(headers):
                try:
                    value = float(values[i])
                    stats_dict[header].append(value)
                except ValueError:
                    # Non-numeric value, skip
                    pass
    
    return stats_dict

def calculate_summary_stats(values):
    """Calculate summary statistics for a list of values."""
    arr = np.array(values)
    return {
        'mean': np.mean(arr),
        'median': np.median(arr),
        'std': np.std(arr),
        'min': np.min(arr),
        'max': np.max(arr),
        'q25': np.percentile(arr, 25),
        'q75': np.percentile(arr, 75),
        'n': len(arr)
    }

def compare_distributions(legacy_values, edited_values, stat_name, alpha=0.05):
    """Compare distributions between legacy and edited versions."""
    result = {
        'stat_name': stat_name,
        'legacy': calculate_summary_stats(legacy_values),
        'edited': calculate_summary_stats(edited_values)
    }
    
    # Calculate relative differences
    if result['legacy']['mean'] != 0:
        result['rel_diff_mean'] = ((result['edited']['mean'] - result['legacy']['mean']) / 
                                   abs(result['legacy']['mean']) * 100)
    else:
        result['rel_diff_mean'] = 0 if result['edited']['mean'] == 0 else float('inf')
    
    if result['legacy']['median'] != 0:
        result['rel_diff_median'] = ((result['edited']['median'] - result['legacy']['median']) / 
                                     abs(result['legacy']['median']) * 100)
    else:
        result['rel_diff_median'] = 0 if result['edited']['median'] == 0 else float('inf')
    
    # Statistical tests
    # Kolmogorov-Smirnov test for distribution equality
    ks_stat, ks_pvalue = stats.ks_2samp(legacy_values, edited_values)
    result['ks_statistic'] = ks_stat
    result['ks_pvalue'] = ks_pvalue
    
    # Mann-Whitney U test (non-parametric test for medians)
    mw_stat, mw_pvalue = stats.mannwhitneyu(legacy_values, edited_values, alternative='two-sided')
    result['mw_statistic'] = mw_stat
    result['mw_pvalue'] = mw_pvalue
    
    # Levene's test for equality of variances
    levene_stat, levene_pvalue = stats.levene(legacy_values, edited_values)
    result['levene_statistic'] = levene_stat
    result['levene_pvalue'] = levene_pvalue
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt((result['legacy']['std']**2 + result['edited']['std']**2) / 2)
    if pooled_std > 0:
        result['cohens_d'] = (result['edited']['mean'] - result['legacy']['mean']) / pooled_std
    else:
        result['cohens_d'] = 0
    
    return result

def print_comparison_table(results):
    """Print a formatted comparison table."""
    print("\n" + "="*120)
    print(f"{'Statistic':<15} {'Legacy Mean':<12} {'Edited Mean':<12} {'Diff %':<8} "
          f"{'Legacy Std':<11} {'Edited Std':<11} {'KS p-val':<10} {'MW p-val':<10}")
    print("="*120)
    
    for r in results:
        if 'legacy' in r and 'edited' in r:
            print(f"{r['stat_name']:<15} "
                  f"{r['legacy']['mean']:>11.4f} "
                  f"{r['edited']['mean']:>11.4f} "
                  f"{r['rel_diff_mean']:>7.1f}% "
                  f"{r['legacy']['std']:>10.4f} "
                  f"{r['edited']['std']:>10.4f} "
                  f"{r['ks_pvalue']:>9.4f} "
                  f"{r['mw_pvalue']:>9.4f}")

def main():
    parser = argparse.ArgumentParser(description='Statistical comparison of nicestats outputs')
    parser.add_argument('legacy_file', help='Legacy nicestats output file')
    parser.add_argument('edited_file', help='Edited nicestats output file')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance level (default: 0.05)')
    parser.add_argument('--output', help='Output CSV file for detailed results')
    parser.add_argument('--verbose', action='store_true', help='Show detailed statistics')
    
    args = parser.parse_args()
    
    # Parse both files
    print(f"Parsing {args.legacy_file}...")
    legacy_stats = parse_nicestats_output(args.legacy_file)
    
    print(f"Parsing {args.edited_file}...")
    edited_stats = parse_nicestats_output(args.edited_file)
    
    # Get all unique statistic names
    all_stats = set(legacy_stats.keys()) & set(edited_stats.keys())  # Only stats in both
    
    # Compare each statistic
    results = []
    for stat_name in sorted(all_stats):
        if len(legacy_stats[stat_name]) > 0 and len(edited_stats[stat_name]) > 0:
            result = compare_distributions(legacy_stats[stat_name], edited_stats[stat_name], 
                                         stat_name, args.alpha)
            results.append(result)
    
    # Print summary table
    print_comparison_table(results)
    
    # Detailed output if requested
    if args.verbose:
        print("\n" + "="*80)
        print("DETAILED COMPARISON")
        print("="*80)
        
        for r in results:
            print(f"\n{r['stat_name']}:")
            print(f"  Legacy: mean={r['legacy']['mean']:.4f}, median={r['legacy']['median']:.4f}, "
                  f"std={r['legacy']['std']:.4f}, n={r['legacy']['n']}")
            print(f"  Edited: mean={r['edited']['mean']:.4f}, median={r['edited']['median']:.4f}, "
                  f"std={r['edited']['std']:.4f}, n={r['edited']['n']}")
            print(f"  Relative difference: mean={r['rel_diff_mean']:.2f}%, median={r['rel_diff_median']:.2f}%")
            print(f"  Cohen's d (effect size): {r['cohens_d']:.4f}")
            print(f"  KS test: p={r['ks_pvalue']:.4f} {'*' if r['ks_pvalue'] < args.alpha else ''}")
            print(f"  MW test: p={r['mw_pvalue']:.4f} {'*' if r['mw_pvalue'] < args.alpha else ''}")
            print(f"  Levene test: p={r['levene_pvalue']:.4f} {'*' if r['levene_pvalue'] < args.alpha else ''}")
    
    # Summary statistics
    significant_ks = sum(1 for r in results if r['ks_pvalue'] < args.alpha)
    significant_mw = sum(1 for r in results if r['mw_pvalue'] < args.alpha)
    large_mean_diff = sum(1 for r in results if abs(r['rel_diff_mean']) > 5.0)
    large_effect = sum(1 for r in results if abs(r['cohens_d']) > 0.8)
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total statistics compared: {len(results)}")
    print(f"Significant KS tests (p < {args.alpha}): {significant_ks}")
    print(f"Significant MW tests (p < {args.alpha}): {significant_mw}")
    print(f"Large mean differences (>5%): {large_mean_diff}")
    print(f"Large effect sizes (|d| > 0.8): {large_effect}")
    
    # Save detailed results if requested
    if args.output:
        df_data = []
        for r in results:
            row = {
                'statistic': r['stat_name'],
                'legacy_mean': r['legacy']['mean'],
                'legacy_median': r['legacy']['median'],
                'legacy_std': r['legacy']['std'],
                'legacy_n': r['legacy']['n'],
                'edited_mean': r['edited']['mean'],
                'edited_median': r['edited']['median'],
                'edited_std': r['edited']['std'],
                'edited_n': r['edited']['n'],
                'rel_diff_mean_%': r['rel_diff_mean'],
                'rel_diff_median_%': r['rel_diff_median'],
                'cohens_d': r['cohens_d'],
                'ks_pvalue': r['ks_pvalue'],
                'mw_pvalue': r['mw_pvalue'],
                'levene_pvalue': r['levene_pvalue']
            }
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        df.to_csv(args.output, index=False)
        print(f"\nDetailed results saved to: {args.output}")
    
    # Overall assessment
    print("\nOVERALL ASSESSMENT:")
    if significant_ks == 0 and significant_mw == 0 and large_mean_diff == 0:
        print("✓ EXCELLENT: No significant differences detected between versions")
    elif large_mean_diff == 0 and large_effect == 0:
        print("✓ GOOD: Some statistical differences detected, but effect sizes are small")
    elif large_mean_diff <= 2:
        print("⚠ WARNING: Some notable differences detected between versions")
    else:
        print("✗ CONCERN: Substantial differences detected between versions")

if __name__ == '__main__':
    main()