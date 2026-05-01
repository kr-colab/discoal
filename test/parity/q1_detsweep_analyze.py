#!/usr/bin/env python3
"""Compare niceStats output distributions between closed-form and Euler
deterministic-sweep modes. Bonferroni-corrected KS at p > 0.01.

Usage: q1_detsweep_analyze.py <results_dir>
"""
import sys
import os
import re
from pathlib import Path

import numpy as np
from scipy import stats

def parse_stats(path):
    cols = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for h in header:
            cols[h] = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != len(header):
                continue
            for h, p in zip(header, parts):
                try:
                    cols[h].append(float(p))
                except ValueError:
                    pass
    return cols

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    results = Path(sys.argv[1])
    configs = set()
    for p in results.glob('*_closed.stats'):
        configs.add(p.name.replace('_closed.stats', ''))
    configs = sorted(configs)

    # collect all (config, statistic) pairs and run KS
    comparisons = []
    for cfg in configs:
        a = parse_stats(results / f'{cfg}_closed.stats')
        b = parse_stats(results / f'{cfg}_euler.stats')
        common = sorted(set(a.keys()) & set(b.keys()))
        for stat in common:
            arr_a = np.asarray(a[stat], dtype=float)
            arr_b = np.asarray(b[stat], dtype=float)
            if len(arr_a) < 100 or len(arr_b) < 100:
                continue
            if np.all(arr_a == arr_a[0]) and np.all(arr_b == arr_b[0]):
                continue
            D, p = stats.ks_2samp(arr_a, arr_b)
            comparisons.append((cfg, stat, D, p, len(arr_a), len(arr_b)))

    if not comparisons:
        print('no comparisons made')
        sys.exit(1)

    n = len(comparisons)
    alpha = 0.01 / n  # Bonferroni
    print(f'{n} comparisons, Bonferroni alpha = {alpha:.2e}')
    print(f"{'config':<24} {'stat':<16} {'D':>10} {'p':>10} {'reject?'}")
    rejected = 0
    for cfg, stat, D, p, na, nb in comparisons:
        flag = '  **' if p < alpha else ''
        if p < alpha:
            rejected += 1
        print(f'{cfg:<24} {stat:<16} {D:>10.4f} {p:>10.2e}{flag}')

    print()
    if rejected == 0:
        print(f'PASS: all {n} comparisons within Bonferroni-corrected p > {alpha:.2e}')
        print('CONCLUSION: detSweepFreq and Euler are statistically indistinguishable')
        print('            under constant N for the tested grid.')
        sys.exit(0)
    else:
        print(f'FAIL: {rejected}/{n} comparisons reject equality at Bonferroni p < {alpha:.2e}')
        sys.exit(1)

if __name__ == '__main__':
    main()
