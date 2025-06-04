# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

- `make discoal` - Build the main discoal executable
- `make clean` - Remove all build artifacts
- `make run_tests` - Build and run all unit tests
- Individual unit tests: `make test_node`, `make test_event`, `make test_mutations`, `make test_node_operations`

## Testing

**IMPORTANT: Always run the comprehensive testing suite when making changes to ensure no regressions.**

### Comprehensive Testing Framework (PRIMARY)
The comprehensive testing suite in `testing/` directory is based on all documented examples and must be run for any changes:

**Quick validation (use during development):**
```bash
cd testing/
./focused_validation_suite.sh
```

**Full validation (required before commits):**
```bash
cd testing/
./comprehensive_validation_suite.sh
```

The testing framework:
- Tests all documented functionality from discoaldoc.tex
- Compares optimized vs legacy versions for identical output
- Profiles memory usage and detects regressions
- Provides 21 comprehensive test cases covering all features
- Ensures trajectory optimization maintains compatibility

### Unit Tests
Unit tests use a custom Unity testing framework in `test/unit/unity.h`. Tests are located in `test/unit/`:
- `test_node.c` - Tests rootedNode structure operations
- `test_mutations.c` - Tests mutation handling
- `test_event.c` - Tests event structure
- `test_node_operations.c` - Tests node operations

### Legacy Integration Tests
Bash scripts in root directory test different aspects:
- Performance: `test_performance_*.sh`
- Memory: `test_memory_*.sh`, `test_memory_leaks.sh`
- Features: `test_recomb_*.sh`, `test_muts_*.sh`, `test_gene_conversion.sh`
- Stress: `test_stress_muts.sh`, `test_extreme_stress_muts.sh`

## Architecture

### Core Components

1. **discoal_multipop.c** - Main program entry point with command-line parsing and simulation loop
2. **discoalFunctions.c** - Core simulation functions (coalescence, recombination, mutations)
3. **alleleTraj.c** - Allele trajectory calculations for sweep simulations
4. **ranlibComplete.c** - Random number generation library

### Key Data Structures

- **rootedNode** (`discoal.h:34-52`) - Core coalescent tree node with:
  - Recombination support (left/right parents and children)
  - Dynamic arrays for mutations (`muts`) and ancestral sites (`ancSites`)
  - Memory management with capacity tracking
  - Population and sweep population tracking

- **event** (`discoal.h:60-68`) - Demographic events (population size changes, splits, admixture)

### Memory Management

The codebase has been optimized for memory efficiency:
- Dynamic allocation for breakpoints, mutations, and ancestral sites
- Capacity tracking with `*Capacity` fields in structures
- Memory allocation functions: `initialize*()`, `ensure*Capacity()`, `cleanup*()`
- Multiple executable versions exist for comparison (legacy vs optimized)

### Simulation Architecture

1. **Parameter parsing** in main()
2. **Initialization** of global arrays and structures
3. **Event-driven simulation** processing demographic events chronologically
4. **Coalescent process** with recombination and mutation
5. **Output generation** in ms format or tree format (with `-T` flag)

### Key Global Variables

- `nodes[]` - Active nodes in coalescent tree
- `breakPoints[]` - Recombination breakpoints (dynamically allocated)
- `activeMaterial[]` - Tracks which sites have ancestral material
- `events[]` - Array of demographic events sorted by time

## Code Conventions

- C style with global variables for simulation state
- Memory-efficient dynamic allocation patterns
- Extensive use of function pointers for different simulation modes
- Consistent naming: `*AtTime*()` for time-based operations, `*Popn*()` for population-specific operations

## Development Workflow

### Required Testing Workflow
1. **During development**: Run `cd testing/ && ./focused_validation_suite.sh` frequently
2. **Before commits**: Run `cd testing/ && ./comprehensive_validation_suite.sh` to ensure no regressions
3. **After major changes**: Run full comprehensive suite and document any performance improvements
4. **Never commit without validation**: The testing suite must pass before any code changes are committed

### Optimization Guidelines
- Always establish baseline performance with testing suite before optimization
- Use focused testing suite for rapid iteration during optimization work
- Document memory savings and success rate improvements in commit messages
- Update test cases when adding new functionality
- Maintain backward compatibility verified through identical output comparison

## Git Workflow and Branch Management

### Branch Strategy
- **master branch**: Protected baseline containing stable legacy code for comparison
- **Feature branches** (e.g., `mem`, `optimization`, etc.): Active development branches
- **NEVER commit to master**: All development work must be done on feature branches

### Legacy Testing Convention
- Legacy binary (`discoal_legacy_backup`) is automatically built from master branch
- This ensures consistent baseline for functional compatibility testing
- Master branch serves as the "source of truth" for original functionality
- All optimizations must maintain identical output to master branch version

### Commit Guidelines
- Only commit to feature branches, never to master
- When making commits, don't add Claude as a co-author
- Run testing suite before all commits to ensure no regressions
- Document optimization benefits and performance improvements in commit messages

### Testing Requirements
- Legacy version always built from master branch via `git show master:filename`
- Optimized version built from current working branch
- All tests compare feature branch (optimized) vs master branch (legacy)
- This guarantees that optimizations maintain compatibility with baseline functionality