# Validation Scripts

This directory contains validation and testing scripts for discoal.

## Statistical Comparison Tools

### statistical_comparison.py

A comprehensive tool for comparing summary statistics between different versions of discoal. This script parses nicestats output files and performs statistical tests to assess whether two versions produce similar results.

**Usage:**
```bash
python3 statistical_comparison.py legacy_file.stats edited_file.stats [options]
```

**Features:**
- Compares means, medians, and standard deviations for all nicestats output
- Performs statistical tests:
  - Kolmogorov-Smirnov test (distribution equality)
  - Mann-Whitney U test (median comparison)
  - Levene's test (variance equality)
  - Cohen's d (effect size calculation)
- Provides detailed summary with relative differences
- Can output results to CSV for further analysis

**Options:**
- `--alpha`: Significance level for tests (default: 0.05)
- `--output`: Save detailed results to CSV file
- `--verbose`: Show detailed per-statistic comparisons

**Example:**
```bash
# Basic comparison
python3 statistical_comparison.py test1_legacy.stats test1_edited.stats

# Detailed comparison with CSV output
python3 statistical_comparison.py test1_legacy.stats test1_edited.stats --verbose --output results.csv
```

## Validation Test Suites

### demographic_events_validation.sh

Tests demographic event handling across different discoal versions, ensuring population events (admixture, migration, merges) work correctly.

**Usage:**
```bash
./demographic_events_validation.sh [num_replicates]
```

Default: 100 replicates

### Other Validation Scripts

- `comprehensive_validation_suite.sh` - Full test suite comparing outputs
- `statistical_validation_suite.sh` - Statistical properties validation
- `msprime_comparison_suite.sh` - Compare discoal with msprime
- `focused_validation_suite.sh` - Quick validation for common use cases

## Notes

- The statistical comparison tool expects nicestats output format (12 columns)
- Different RNGs between versions will produce variation, but statistical properties should be preserved
- Use higher replicate counts (100+) for more robust statistical comparisons