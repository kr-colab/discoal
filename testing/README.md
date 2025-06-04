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

# Test trajectory optimization specifically
./test_trajectory_memory_comparison.sh
```

## Test Suites

| Script | Purpose | Runtime | Test Cases |
|--------|---------|---------|------------|
| `comprehensive_validation_suite.sh` | Complete validation of all documented functionality | 15-30 min | 21 tests |
| `focused_validation_suite.sh` | Quick regression testing of core features | 3-5 min | 10 tests |
| `test_trajectory_memory_comparison.sh` | Trajectory optimization validation | 2-3 min | 5 tests |
| `test_trajectory_lazy_optimization.sh` | Legacy trajectory validation | 5-10 min | 10 tests |

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
  - `discoal_trajectory_optimized`: Built from current working branch
  - `discoal_legacy_backup`: Built from master branch for stable baseline comparison
- **System**: BSD `time` command for memory profiling
- **Optional**: `timeout` command (install with `brew install coreutils` on macOS)

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