# Discoal Testing Framework

This directory contains the comprehensive testing framework for validating discoal optimizations and ensuring functional compatibility.

## Quick Start

```bash
# Navigate to testing directory
cd testing/

# Run comprehensive validation (all documented examples)
./comprehensive_validation_suite.sh

# Run quick validation (core functionality)
./focused_validation_suite.sh

# Compare with msprime simulator
./msprime_comparison_suite.sh

# Run nicestats comparison (requires niceStats - see below)
./nicestats_comparison_suite.sh
```

## Building niceStats

Several test suites require the `niceStats` utility for computing population genetics summary statistics:

```bash
# From project root
make niceStats
```

This builds `niceStats` from the msUtils library and creates a symlink in the project root.

## Test Suites

| Script | Purpose | Runtime |
|--------|---------|---------|
| `comprehensive_validation_suite.sh` | Complete validation of all documented functionality | 15-30 min |
| `focused_validation_suite.sh` | Quick regression testing of core features | 3-5 min |
| `nicestats_comparison_suite.sh` | Statistical comparison of summary statistics | 5-15 min |
| `msprime_comparison_suite.sh` | Statistical comparison with msprime simulator | 5-10 min |
| `demes_validation_suite.sh` | Demes models through a YAML config (`-Y`): each is accepted or, where marked, rejected (in CI) | under 1 s |
| `demes_output_comparison_suite.sh` | Demes models through `-Y` give output identical to the equivalent command line (in CI) | under 1 s |
| `yaml_validation_suite.sh` | YAML configuration file validation | 2 min |
| `yaml_parser_suite.sh` | YAML parser regression suite (walks `testing/debug_yaml/`) | < 5 sec |
| `units_convention_check.py` | Checks the units of times, growth rates, and migration rates against msprime (needs msprime and demes) | under 1 min |
| `slim_sweep_check.py` | Compares discoal's hard, soft (`-f`), partial (`-c`) and neutral (`-wn`) sweeps with forward simulations in SLiM (stdpopsim's constant-size script with discoal's sweep logic), and the hard sweep with msprime; run by hand (needs stdpopsim, SLiM, pyslim, msprime) | about 7 min with 120 workers |
| `sweep_snp_check.py` | Checks the selected SNP of partial sweeps (`-c`): its position, the binomial count of its carriers, `-h`, one site per recurrent sweep, and `-ts` output (needs numpy, tskit) | about 15 s |
| `trajectory_file_check.py` | Checks `-S`: a file holding discoal's deterministic trajectory gives the same branch length and diversity as `-wd`; a file that starts at frequency 0.5 gives a partial sweep whose selected SNP is carried by Binomial(n, 0.5) samples; `-S` runs with ancient samples (needs numpy, tskit) | about 3 s |
| `recurrent_sweep_check.py` | Checks that recurrent sweeps (`-R`) follow migration changes (`-eM`, `-em`) between sweeps: at a negligible sweep rate the output must equal the neutral output | about 1 s |
| `fuzz_msprime_check.py` | Random command-line demographic models (sizes, growth, migration, admixture, merges, ancient samples, recombination) against msprime, compared with KS tests on branch statistics; `--models`, `--reps`, `--seed`, `--discoal` | a few seconds |
| `stdpopsim_parity.py` | Compares discoal and msprime branch-mode statistics on every stdpopsim demographic model, random generic stdpopsim models (`--generic`), and demes files (`--demes-dir`, e.g. the examples and test cases of the demes specification), with ancient samples at stdpopsim's default sampling times and, in "+ancient" variants of the generic models and demes files, from every deme at random times (needs stdpopsim, msprime, demes, tskit, scipy). `--sample K --sample-seed S` adds K models drawn at random. CI runs 10 core models and 80 random ones, seeded by the commit, on every push and PR; `.github/workflows/parity-weekly.yml` runs everything weekly | about 15 min with 24 workers for the catalog and generic models |

## Documentation

- **`TEST_SUITE_DOCUMENTATION.md`**: Complete documentation of the testing framework
- **`README.md`**: This quick reference guide

## Test Results

Test results are stored in timestamped directories:
```
comprehensive_validation_YYYYMMDD_HHMMSS/
├── validation_summary.txt           # Overall results
├── *_optimized.out                  # Optimized outputs  
├── *_legacy.out                     # Legacy outputs
├── *_optimized_memory.txt           # Memory profiles
└── *_legacy_memory.txt              # Memory profiles
```

## Usage Guidelines

### Before Optimization
1. Run comprehensive suite to establish baseline
2. Document current success rates and memory usage

### During Development
1. Use focused suite for rapid iteration testing
2. Run relevant subset tests for specific optimizations

### After Optimization
1. Run comprehensive suite to validate no regressions
2. Document improvements in success rates and memory usage
3. Update test cases if new capabilities are added

## Requirements

- **Git setup**: Access to both current feature branch and master branch
- **Build system**: Automated binary generation via `make test_binaries`
  - `discoal_edited`: Built from current working branch
  - `discoal_legacy_backup`: Built from master branch for stable baseline comparison
- **System**: BSD `time` command for memory profiling
- **Optional**: `timeout` command (install with `brew install coreutils` on macOS)

## msprime Comparison Tests

The `msprime_comparison_suite.sh` compares discoal with the msprime coalescent simulator to ensure statistical equivalence:

```bash
# Run comparison tests (requires conda environment)
conda activate discoal_dev
./msprime_comparison_suite.sh

# Analyze results from a previous run
./analyze_msprime_comparison.py msprime_comparison_20240611_120000/
```

The suite tests:
- Neutral coalescent with various sample sizes
- Models with recombination at different rates
- Statistical equivalence using Kolmogorov-Smirnov and Mann-Whitney U tests

## Adding New Tests

Extend the `TEST_CASES` array in any test script:

```bash
declare -a TEST_CASES=(
    # Format: "category:name:command_args:expected_behavior"
    "new_category:test_name:discoal_args:both_succeed"
)
```

Expected behaviors:
- `both_succeed`: Both versions should succeed with identical output
- `optimized_preferred`: Optimized version should succeed, legacy may fail
- `legacy_only`: Legacy version should succeed (regression detection)

## Maintenance

- **Monthly**: Run comprehensive suite on latest codebase
- **Pre-release**: Full validation with all test suites  
- **Post-optimization**: Update documentation with new capabilities

For detailed information, see `TEST_SUITE_DOCUMENTATION.md`.