# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

discoal is a coalescent simulation program for generating genetic data under models with recombination, population structure, and selection. It extends Richard Hudson's ms program with the ability to simulate selective sweeps and other forms of natural selection.

## Build Commands

```bash
# Build the main executable
make discoal

# Build with specific compiler
make CC=gcc discoal
make CC=clang discoal

# Clean build artifacts
make clean

# Build and run all tests
make test

# Run specific test suites
make run_tests              # Unit tests only
make comprehensive_validation_test
make focused_validation_test
make memory_comparison_test
make msvalidate             # Statistical validation against msprime
```

## Testing

### Unit Tests
```bash
# Run all unit tests
make run_tests

# Run a single unit test (from project root)
./build/test_ancestrySegment
./build/test_activeSegment
```

### Integration Tests
```bash
# Run comprehensive validation (tests all examples from documentation)
make comprehensive_validation_test

# Quick validation of core functionality
make focused_validation_test

# Memory usage comparison
make memory_comparison_test
```

## Code Architecture

### Core Components

- **src/core/discoal_multipop.c**: Main program entry point and command-line parsing
- **src/core/discoalFunctions.c/h**: Core simulation functions, coalescent process, recombination, selection
- **src/core/ancestrySegment.c/h**: Tracks ancestry segments along chromosomes
- **src/core/activeSegment.c/h**: Manages active lineages during simulation
- **src/core/alleleTraj.c/h**: Generates allele frequency trajectories for selection models

### Key Data Structures

- **struct ancestrySegment**: Represents a contiguous genomic segment with shared ancestry
- **struct activeSegment**: Tracks actively coalescing lineages
- **struct event**: Represents demographic/selection events
- **struct chunk**: Memory management for ancestry segments

### Memory Management

The codebase uses custom memory management for ancestry segments:
- Pre-allocated chunks of segments to reduce malloc overhead
- Recycling of freed segments
- Global pools for efficient allocation/deallocation

### Random Number Generation

Two RNG implementations available:
- Legacy L'Ecuyer RNG (default for backward compatibility)
- Modern xoshiro256++ RNG (use with --xoshiro flag)

### Tree Sequence Integration

- Optional tskit output format via --tree-sequence flag
- Interface in src/tskit/tskitInterface.c
- Converts internal ancestry representation to tskit tables

## Development Guidelines

### Memory Optimization

Recent optimizations have reduced memory usage by 15-99% across scenarios. Key strategies:
- Use chunk-based allocation for ancestry segments
- Minimize per-segment memory overhead
- Efficient data structure packing

### Testing Requirements

Before committing changes:
1. Run unit tests: `make run_tests`
2. Run focused validation: `make focused_validation_test`
3. Check for memory leaks: `make valgrind_test` (if available)

### Adding New Features

1. Update parameter parsing in discoal_multipop.c
2. Implement logic in discoalFunctions.c
3. Add unit tests in test/unit/
4. Add integration tests in testing/
5. Update documentation in docs/discoaldoc.tex

### Common Development Tasks

```bash
# Debug build
make DEBUG=1 discoal

# Profile memory usage
./scripts/run_memory_profiling.sh

# Compare output with msprime
make msvalidate

# Generate comparison plots
python scripts/validation/compare_nicestats_distributions.py
```

## Important Notes

- Maximum sample size: 65,535 (increased from 254)
- Support for up to 2,000 populations
- Backward compatibility maintained with ms output format
- Tree sequence output requires --tree-sequence flag
- Gene conversion uses -gc flag (fraction of recombination events)

## Parameter System Refactoring (dsc_params branch)

### Overview
We are transitioning from scattered global variables to a centralized parameter management system that supports both command-line and YAML configuration.

### New Components

#### Parameter Structure (src/core/params.h/c)
- Hierarchical organization of all simulation parameters
- Main struct: `SimulationParams` containing:
  - `SimulationCore`: samples, replicates, sites, populations
  - `EvolutionaryForces`: mutation, recombination, migration
  - `Demographics`: population sizes, events, ancient samples
  - `SelectionParams`: sweeps, selection coefficients
  - `OutputConfig`: output format settings
  - `PriorConfig`: prior distributions
  - `DebugConfig`: debugging options
  - `RandomConfig`: RNG seeds

#### YAML Support (src/core/yaml_loader.h/c)
- Uses libyaml for robust YAML parsing
- Functions: `yaml_load_params()`, `yaml_save_params()`
- Example configuration:
```yaml
simulation:
  samples: 100
  replicates: 1000
  sites: 10000

evolution:
  mutation:
    theta: 100.0
  recombination:
    rho: 100.0
    gene_conversion:
      ratio: 0.1
      tract_length: 500

selection:
  coefficient: 100.0
  position: 0.5
  mode: stochastic
  time: 0.01
```

### Migration Status
- ✓ Parameter structures defined
- ✓ Command-line parser implemented (backward compatible)
- ✓ YAML loader implemented (using libyaml)
- ✓ Parameter validation framework
- ✓ Comprehensive unit tests
- ⏳ Refactor main code to use SimulationParams
- ⏳ Integrate demes-c for demographic models
- ⏳ Update documentation

### Dependencies
- **libyaml**: Required for YAML configuration support
  - Install: `sudo apt-get install libyaml-dev` (Ubuntu/Debian)
  - The build system links with `-lyaml`

### Future Integration: demes-c
- Located in `reference_code/demes-c/`
- Will provide standardized demographic model support
- Also requires libyaml (dependency alignment achieved)

## Demes Integration (2024-12-22)

### Overview
Successfully integrated demes-c library to support loading demographic models from Demes specification files. This allows discoal to use standardized demographic models that are compatible with other population genetics software.

### Key Implementation Details

1. **Time Units Requirement**: 
   - discoal requires Demes files to use `time_units: 4N` (coalescent time scaling)
   - Also requires `generation_time: 1` when using 4N units

2. **Population Size Scaling**:
   - All population sizes are relative to population 0's present-day size
   - This matches discoal's internal representation where pop 0 is the reference

3. **Epoch Transitions**:
   - Successfully implemented conversion of Demes epoch transitions to discoal's `-en` events
   - Each epoch boundary with a size change creates an EVENT_SIZE_CHANGE

4. **Supported Demes Features**:
   - Single and multiple populations
   - Population splits (ancestors → descendants)
   - Discrete size changes (step functions)
   - Symmetric migration between populations
   - Basic demographic events

5. **Unsupported Demes Features** (documented in `docs/demes_unsupported_features.md`):
   - Exponential growth (discoal only supports step functions)
   - Selfing and cloning rates
   - Pulse admixture events
   - Multi-way admixture (more than 2 ancestors)
   - Time-bounded migration
   - Complex migration matrices

### Test Files
Example Demes files in `test/demes_examples/`:
- `single_pop.yaml`: Simple single population
- `bottleneck.yaml`: Population with bottleneck
- `two_pop_migration.yaml`: Two populations with migration
- `size_changes.yaml`: Multiple epoch transitions

### Usage (pending CLI integration)
```bash
# Currently can be tested with:
./build/test_demes_run test/demes_examples/bottleneck.yaml

# Future usage will be:
./build/discoal -d demographic_model.yaml [other options]
```

## Sample Specification Support (2025-06-23)

### Overview
Implemented support for specifying sample distributions across populations when using Demes files. This allows users to define which populations to sample from and how many samples to take from each, using population names from the Demes file.

### Implementation Details

1. **SampleSpec Structure** (params.h):
   ```c
   typedef struct {
       char    *population;    /* Population name (for Demes) */
       int     size;          /* Number of samples */
       double  time;          /* Sampling time (default 0 = present) */
   } SampleSpec;
   ```

2. **YAML Configuration**:
   - Sample specifications are defined in the `demographics->samples` section
   - Each spec includes: population name, sample size, and optional time
   - Times are in coalescent units (4N) and default to 0 (present day)

3. **Name Resolution**:
   - `demes_apply_sample_specs()` maps population names from Demes to indices
   - Validates that all specified populations exist in the Demes model
   - Updates `params->core.sample_sizes[]` based on specifications

4. **Integration Flow**:
   - `yaml_load_params_with_demes()` loads both Demes file and sample specs
   - Maintains name mapping during Demes conversion
   - Applies sample specifications after demographic model is loaded

### Example Usage

```yaml
# discoal_config.yaml
simulation:
  replicates: 100
  sites: 100000

evolution:
  mutation:
    theta: 10.0
  recombination:
    rho: 10.0

demographics:
  demes_file: model.yaml
  samples:
    - population: deme1
      size: 20
      time: 0      # Present-day samples (default)
    - population: deme2  
      size: 30
    - population: ancient_pop
      size: 10
      time: 0.5    # Ancient samples (parsed but not fully implemented)

output:
  format: ms
```

### Current Limitations
- Ancient samples (time > 0) are parsed but not yet implemented in the simulation engine
- Sample specifications require using the YAML configuration approach
- The main discoal program hasn't been refactored to use the new parameter system yet

### Testing
- Created `test_sample_specification.c` to validate the implementation
- Example configurations in `test/configs/`
- Test Demes files in `test/examples/demes/`

## Minimal Phase 1 Parameter Refactoring (2025-06-23)

### Overview
Completed minimal Phase 1 of refactoring global variables to use SimulationParams structure. Added essential missing parameters while maintaining backward compatibility.

### Parameters Added

1. **SelectionParams**:
   - `recurrent_sweep_rate` - Rate of recurrent sweeps (for -R option)
   - Maps to global `recurSweepRate` variable

2. **OutputConfig**:
   - `mask` - Masking parameter (reserved for future use)
   - Maps to global `mask` variable

3. **Verified existing support**:
   - `Ne` (effective population size) - Already in SimulationCore
   - Maps to global `EFFECTIVE_POPN_SIZE`

### Implementation Details

- Updated `params_create()` to initialize new fields with appropriate defaults
- Modified `params_load_from_args()` to handle -R flag correctly
- Extended YAML parser to support:
  - `selection->recurrent_sweep->rate` 
  - `output->mask`
  - `simulation->effective_size` (already supported)
- Created comprehensive tests for both command-line and YAML loading

### Next Steps
With these parameters in place, the main refactoring task can proceed to replace global variable usage with SimulationParams throughout the codebase.