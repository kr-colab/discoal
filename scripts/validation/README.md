# Validation Scripts

This directory contains validation and testing scripts for discoal.

## Building niceStats

Many validation scripts require the `niceStats` utility:

```bash
# From project root
make niceStats
```

## Statistical Comparison Tools

### compare_nicestats_distributions.py

Compares summary statistics distributions between different versions of discoal using Kolmogorov-Smirnov tests.

**Usage:**
```bash
python3 compare_nicestats_distributions.py legacy_file.stats edited_file.stats
```

**Features:**
- Parses nicestats output files
- Performs KS tests for distribution equality
- Reports mean, SD, median for each statistic
- Flags significant differences

## Validation Scripts

### demographic_events_validation.sh

Tests demographic event handling, ensuring population events (admixture, migration, merges) work correctly.

**Usage:**
```bash
./demographic_events_validation.sh [num_replicates]
```

### test_full_arg_mode.sh

Tests the `-F` flag for full ARG mode vs minimal tree sequence mode.

### test_full_arg_roots.sh

Validates multi-root tree generation with high recombination rates.

## Notes

- Different RNGs between versions will produce statistical variation, but distributional properties should be preserved
- Use higher replicate counts (100+) for more robust statistical comparisons