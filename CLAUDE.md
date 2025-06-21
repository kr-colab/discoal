# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## When you open this file
- **Read the entire file** to understand the current state of the codebase, recent changes, and testing requirements.
- **Do not skip any sections** as they contain important information about the code structure, memory optimizations, and testing procedures.
- **Follow the instructions carefully** to ensure you maintain the integrity of the codebase and adhere to the development workflow.
- **Pay attention to the directory structure** and build commands, as they are crucial for running tests and building the project correctly.
- **Note the current work status** to understand recent changes and ongoing optimizations.
- **Be aware of the testing requirements** to ensure all changes are validated against the comprehensive testing suite.
- **Read the makefile** to understand how to build the project and run tests.
- **read the last few git commits** to see the latest changes and optimizations made to the codebase.
- **read any changes on the current git branch that are uncommitted** to understand ongoing work and ensure you do not duplicate efforts.
- **review the directory structure** to understand where files are located and how to navigate the project.
- **read the readme** to familiarize yourself with the project, its components, and how to run it.
- **acknowledge that you have read this**

## CRITICAL: File Deletion Warning
**NEVER DELETE TEST SUITE FILES**: When cleaning up artifacts, be extremely careful:
- ✅ OK to delete: Test result directories (e.g., `comprehensive_validation_TIMESTAMP/`)
- ✅ OK to delete: Temporary files (`.tmp`, `.out`, `perf.data`)
- ❌ NEVER delete: Test scripts (`*_suite.sh`, `*.sh` in testing/)
- ❌ NEVER delete: Source code or configuration files
- Always use `git status` before cleanup to ensure important files aren't deleted
- If unsure, ask before deleting ANY file

## Current Work Status (as of latest fixes on `mem` branch)

### Recent Memory Optimizations
1. **Dynamic Memory Allocation** (previous commits)
   - Converted fixed arrays to dynamic allocation for nodes, breakpoints, mutations
   - Reduced baseline memory from ~400MB to ~10MB
   - Achieved 15-93% memory savings across test cases

2. **Memory-Mapped Trajectory Storage** (commit df6368b)
   - Implemented file-backed mmap for sweep trajectories
   - Trajectories write directly to temp files during generation
   - Accepted trajectories are memory-mapped for reading
   - Reduced memory usage from ~2GB to ~1.4MB for medium sweeps (99.9% reduction)
   - Enables previously impossible large sweep simulations

3. **Fixed RNG Corruption Bug** (previous fix)
   - Changed `malloc` to `calloc` in `initializeAncSites()` function
   - Resolved output differences between optimized and legacy versions
   - Fixed memory explosion in admixture_model test (4.6GB → 3.4MB)
   - All 21 comprehensive tests now produce identical output

4. **Ancestry Segment Tree Implementation** (completed)
   - Implemented succinct tree structure for tracking ancestral segments
   - Proper interval merging with sweep-line algorithm and automatic coalescing
   - Reference counting infrastructure for future segment sharing
   - Wrapper functions enable seamless transition from arrays
   - Tree-only mode tested with 80% memory reduction for large simulations

5. **Complete Removal of ancSites Arrays** (previous commit)
   - Removed all ancSites array references from codebase
   - Eliminated conditional compilation (USE_ANCESTRY_TREE_ONLY)
   - All ancestry tracking now uses segment trees exclusively
   - Maintained full backward compatibility - all 21 tests pass
   - Memory savings of 80% for large simulations (50 samples, 50k sites)

6. **Segment Sharing with Reference Counting** (previous commit)
   - Implemented reference counting for all ancestry segments
   - Smart copying shares immutable segments instead of deep copying
   - Optimized merge operations retain child segments
   - Efficient split operations share segments when possible
   - Additional memory savings: 10-16% on top of previous optimizations

7. **AVL Tree Indexing for Ancestry Lookups** (completed)
   - Implemented self-balancing AVL trees for O(log n) ancestry lookups
   - Optimized threshold: builds AVL tree when ≥3 segments exist
   - Performance profiling showed 21% speedup for high recombination (r=200)
   - Enables extreme recombination scenarios with 97% memory reduction
   - Trade-off: slower for very high recombination (r>1000) but prevents memory exhaustion

8. **Active Material Segment-Based Tracking** (completed)
   - Replaced legacy `activeMaterial[]` array with segment-based structure
   - Tracks active genomic regions as intervals instead of per-site arrays
   - Integrated with AVL tree indexing for O(log n) lookups when needed
   - Memory savings: 70-99% for high recombination scenarios
   - Performance improvements: Up to 47x faster for memory-intensive simulations
   - Maintains 100% output compatibility with legacy implementation

9. **msprime Comparison Suite** (completed)
   - Created comprehensive comparison framework with msprime coalescent simulator
   - Implemented proper parameter scaling between discoal and msprime conventions
   - All 10 test cases pass statistical validation (KS test p>0.05)
   - Validates that memory optimizations maintain statistical correctness
   - Performance tracking integrated (discoal averages 56x faster than msprime)
   - Documentation includes detailed parameter conversion guide

10. **Mutation Handling Optimization - Phase 1** (completed)
   - Identified O(n²) bottleneck in duplicate mutation detection
   - Implemented hash table for O(n) duplicate detection
   - Hash table size increased to 40,009 (prime > MAXMUTS) to handle extreme cases
   - Performance improvements:
     - θ=1,000: 1.50x speedup
     - θ=5,000: 1.97x speedup  
     - θ=10,000: 2.16x speedup
   - Added high mutation rate tests to comprehensive suite

11. **Mutation Handling Optimization - Phase 2** (completed)
   - Implemented sorting of mutations after placement for binary search
   - Modified hasMutation() to use binary search for arrays >10 elements
   - sortAllMutations() called within makeGametesMS after all mutations placed
   - Combined Phase 1+2 performance improvements:
     - θ=1,000: 1.95x speedup (0.04s vs 0.14s)
     - θ=5,000: 6.46x speedup (0.21s vs 2.82s)
     - θ=10,000: >50x speedup (0.37s vs >30s timeout)
   - Maintains 100% output compatibility with legacy version

12. **Mutation Handling Optimization - Phase 3** (completed)
   - Pre-compute mutation presence matrix for output generation
   - Eliminates O(n×m) ancestry lookups during output phase
   - Single pre-computation pass stores results in char matrix
   - Performance improvements (on top of Phase 1+2):
     - 100 samples, θ=1000: 6.4x additional speedup
     - 200 samples, θ=500: 6.2x additional speedup
     - Most effective with large sample sizes and many segregating sites
   - Minimal memory overhead (temporary sampleSize × mutNumber bytes)
   - Memory overhead minimal (<1-3%)

13. **Comprehensive Unit Test Framework** (completed)
   - Migrated from custom Unity implementation to official Unity framework
   - Created unified test runner for all unit tests
   - Added 5 new test suites with 46 additional tests:
     - Ancestry segment tests (13 tests): Tree operations, reference counting, splitting/merging
     - Active segment tests (12 tests): Active material tracking, fixed region removal
     - Trajectory tests (12 tests): Memory-mapped files, cleanup, persistence
     - Coalescence/recombination tests (11 tests): Core simulation operations
     - Memory management tests (17 tests): Dynamic arrays, capacity growth, stress testing
   - Total test coverage: 77 unit tests across 9 test files
   - All tests integrated into Makefile with individual and unified execution options
   - Updated development documentation with comprehensive testing guide

### Test Results Summary
- **Success Rate**: 31/31 tests pass (100%) for both versions (includes 4 new high mutation tests)
- **Output Compatibility**: 31/31 tests produce identical output
- **Memory Improvements**: Average 23% reduction across all tests
  - Multipop models: up to 58% reduction
  - Selection sweeps: up to 52% reduction
  - High recombination: 97-99% reduction
  - Large simulations: 80% reduction (tree-only implementation)
  - No recombination: Additional 10.8% with segment sharing
  - Large sample sizes: Additional 16.1% with segment sharing
  - Extreme recombination (r=10000): 98% reduction with 33-47x speedup
- **Performance Improvements** (mutation optimization Phase 1+2+3):
  - High mutation (θ=1000): Up to 12.4x speedup with large samples
  - Very high mutation (θ=5000): Up to 40x speedup with large samples
  - Extreme mutation (θ=10000): >50x speedup (legacy times out)
  - Phase 3 adds 1.14x-6.4x speedup depending on sample size

### Active TODOs
- [x] Complete removal of ancSites array from main codebase ✓
- [x] Implement segment sharing using reference counting ✓
- [x] Add AVL tree indexing for high-recombination scenarios ✓
- [x] Replace activeMaterial array with segment-based tracking ✓
- [x] Add msprime comparison suite for statistical validation ✓
- [x] Optimize mutation duplicate detection with hash table (Phase 1) ✓
- [x] Optimize hasMutation() with binary search (Phase 2) ✓
- [x] Pre-compute mutation presence matrix for output (Phase 3) ✓
- [ ] Document memory optimization techniques in main README
- [ ] Phase 4: Memory layout optimizations for better cache efficiency
- [ ] Optimize pickNodePopn with per-population node lists

### Implementation Details
- Trajectory files: `/tmp/discoal_traj_<pid>_<time>_<rand>.tmp`
- Rejected trajectories cleaned up immediately
- Signal handlers ensure cleanup on exit
- 500M step safety limit prevents runaway trajectories
- All dynamic arrays now use `calloc` for zero-initialization
- Ancestry segment trees use interval-based representation with coalescing
- Reference counting enables segment sharing between nodes
- Immutable segments can be safely shared, reducing memory overhead
- Wrapper functions in `ancestryWrapper.h` provide clean API
- AVL trees built for segment lists with ≥3 segments for O(log n) lookups
- AVL indexing prevents performance degradation under high recombination
- Active material tracked via `ActiveMaterial` structure with segment list
- Active segments coalesce automatically to minimize memory overhead
- AVL tree built for active segments when count ≥10 for efficient queries
- Mutation hash table uses MUTATION_HASH_SIZE=40009 (prime > MAXMUTS)
- Binary search for hasMutation() with adaptive threshold (>10 mutations)
- Mutations sorted after placement in makeGametesMS for O(log n) lookups
- Output generation uses pre-computed presence matrix (Phase 3)
- Matrix computation reduces getAncestryCount from 17.81% to 10.22% of runtime

### Known Issues
- Test executables must be named `discoal_edited` and `discoal_legacy_backup`

## Directory Structure

**IMPORTANT**: When working with bash commands, be aware of the directory structure:
- **Root directory** (`/Users/adk/github/discoal/`): Contains source files, Makefile, and main executables
- **Testing directory** (`/Users/adk/github/discoal/testing/`): Contains all test scripts. NEVER DELETE THESE.
- **Unit test directory** (`/Users/adk/github/discoal/test/unit/`): Contains unit test files

**Note about shell sessions**: The bash tool maintains persistent state, so after `cd testing/`, subsequent commands run in that directory. To avoid confusion:
- Use absolute paths when possible
- Or explicitly return to root with `cd ..` or `cd /Users/adk/github/discoal`
- Check current directory with `pwd` if uncertain

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
- Provides 27 comprehensive test cases covering all features
- Ensures trajectory optimization maintains compatibility

### Statistical Validation Suite (NEW)
The statistical validation suite runs multiple replicates of each test case to ensure statistical equivalence between versions:

**Usage:**
```bash
cd testing/
./statistical_validation_suite.sh              # 100 replicates, auto mode
./statistical_validation_suite.sh parallel 50  # 50 replicates, parallel mode
./statistical_validation_suite.sh 200          # 200 replicates, auto mode
```

**Features:**
- Runs configurable number of replicates (default: 100) for each test case
- Extracts segregating sites counts from each replicate
- Performs Kolmogorov-Smirnov test to compare distributions
- Tests null hypothesis that outputs come from same distribution
- Supports parallel execution with GNU parallel
- Provides detailed statistical summaries (mean, SD, min, max)
- Identifies cases where distributions differ significantly

**Key statistics analyzed:**
- Number of segregating sites per replicate
- Memory usage and performance metrics
- Distribution comparisons at 0.05 significance level

This suite is essential for validating that optimizations don't introduce systematic biases or alter the statistical properties of the simulations.

### msprime Comparison Suite (NEW)
The msprime comparison suite validates discoal against the well-established msprime coalescent simulator:

**Usage:**
```bash
cd testing/
./msprime_comparison_suite.sh
```

**Features:**
- Compares discoal and msprime across 10 test scenarios
- Tests neutral models with various sample sizes and recombination rates
- Tests selection models (hard sweeps) with different strengths and ages
- Performs statistical validation using Kolmogorov-Smirnov tests
- Tracks runtime performance for both simulators
- All tests currently pass (p > 0.05) confirming statistical equivalence

**Parameter Scaling:**
The suite handles the complex parameter conversions between discoal and msprime:
- Population size: Ne=0.5 for msprime with diploid samples
- Mutation/recombination rates: Converted from locus-wide to per-bp
- Selection coefficients: Properly scaled with msprime's fitness model
- Time scaling: Special handling for sweep timing parameters

See `docs/development.rst` for detailed parameter conversion documentation.

### Unit Tests
Unit tests use the official Unity testing framework (https://github.com/ThrowTheSwitch/Unity). Tests are located in `test/unit/`:
- `test_node.c` - Tests rootedNode structure operations (3 tests)
- `test_event.c` - Tests event structure (2 tests)
- `test_node_operations.c` - Tests node operations (4 tests)
- `test_mutations.c` - Tests mutation handling (3 tests)
- `test_ancestry_segment.c` - Tests ancestry segment trees (13 tests)
- `test_active_segment.c` - Tests active material tracking (12 tests)
- `test_trajectory.c` - Tests trajectory file handling (12 tests)
- `test_coalescence_recombination.c` - Tests core simulation operations (11 tests)
- `test_memory_management.c` - Tests memory management functions (17 tests)

**Running unit tests:**
- `make run_tests` - Run all tests individually
- `make run_all_tests` - Run all tests with unified runner
- `make test_<name>` - Build and run specific test suite

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
  - Dynamic arrays for mutations (`muts`)
  - Ancestry segment tree (`ancestryRoot`) for tracking ancestral material
  - Memory management with capacity tracking
  - Population and sweep population tracking

- **event** (`discoal.h:60-68`) - Demographic events (population size changes, splits, admixture)

### Memory Management

The codebase has been optimized for memory efficiency:
- Dynamic allocation for breakpoints and mutations
- Ancestry tracked via segment trees instead of arrays (80% memory reduction)
- Segment sharing via reference counting (additional 10-16% reduction)
- Capacity tracking with `*Capacity` fields in structures
- Memory allocation functions: `initialize*()`, `ensure*Capacity()`, `cleanup*()`
- Memory-mapped files for large trajectory arrays
- Multiple executable versions exist for comparison (legacy vs optimized)

### Trajectory Storage (NEW)
- Sweep trajectories written to temporary files during generation
- Files are memory-mapped read-only after acceptance
- Rejected trajectories cleaned up immediately
- Signal handlers ensure cleanup on exit/interrupt

### Simulation Architecture

1. **Parameter parsing** in main()
2. **Initialization** of global arrays and structures
3. **Event-driven simulation** processing demographic events chronologically
4. **Coalescent process** with recombination and mutation
5. **Output generation** in ms format or tree format (with `-T` flag)

### Key Global Variables

- `nodes[]` - Active nodes in coalescent tree (dynamically allocated)
- `breakPoints[]` - Recombination breakpoints (dynamically allocated)
- `activeMaterial[]` - Tracks which sites have ancestral material
- `events[]` - Array of demographic events sorted by time
- `currentTrajectory` - Memory-mapped pointer to trajectory data
- `trajectoryFd` - File descriptor for mmap'd trajectory

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
- **master-backup branch**: Stable baseline containing legacy code for comparison testing
- **master branch**: Main development branch that can evolve independently
- **Feature branches** (e.g., `mem`, `optimization`, etc.): Active development branches
- **Development workflow**: Work on feature branches, merge to master when ready

### Legacy Testing Convention
- Legacy binary (`discoal_legacy_backup`) is automatically built from master-backup branch
- This ensures consistent baseline for functional compatibility testing
- Master-backup branch serves as the "source of truth" for original functionality
- All optimizations must maintain identical output to master-backup branch version

### Commit Guidelines
- Only commit to feature branches, never to master
- When making commits, don't add Claude as a co-author
- Run testing suite before all commits to ensure no regressions
- Document optimization benefits and performance improvements in commit messages

### Testing Requirements
- Legacy version always built from master-backup branch via `git show master-backup:filename`
- Optimized version built from current working branch
- All tests compare feature branch (optimized) vs master-backup branch (legacy)
- This guarantees that optimizations maintain compatibility with baseline functionality