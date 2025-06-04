# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

- `make discoal` - Build the main discoal executable
- `make clean` - Remove all build artifacts
- `make run_tests` - Build and run all unit tests
- Individual unit tests: `make test_node`, `make test_event`, `make test_mutations`, `make test_node_operations`

## Testing

### Unit Tests
Unit tests use a custom Unity testing framework in `test/unit/unity.h`. Tests are located in `test/unit/`:
- `test_node.c` - Tests rootedNode structure operations
- `test_mutations.c` - Tests mutation handling
- `test_event.c` - Tests event structure
- `test_node_operations.c` - Tests node operations

### Integration Tests
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