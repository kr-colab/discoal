# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the discoal codebase.

## Quick Start
When opening this project:
1. Read this file to understand project structure and recent work
2. Check `git log -5` for recent commits
3. Run `git status` to see uncommitted changes
4. Review the Makefile for build commands
5. Acknowledge that you've read this

## CRITICAL: File Safety
**NEVER DELETE TEST SUITE FILES**:
- ✅ OK to delete: Test result directories (`*_TIMESTAMP/`), temp files (`.tmp`, `.out`)
- ❌ NEVER delete: Test scripts (`*.sh` in testing/), source code, config files
- Always check `git status` before cleanup

## Current Status (`tskit-integration` branch)

### Recent Major Work
1. **tskit Integration** - Full tree sequence recording with tskit C API
   - Tree sequences output with `-ts filename.trees`
   - Automatic simplification prevents fixed sites
   - Complete mutation recording with proper scaling
   - Provenance metadata for reproducibility
   - Edge squashing optimization (50-75% edge reduction)

2. **Memory Optimizations** - Achieved 80-99% memory reduction
   - Dynamic allocation replaced fixed arrays
   - Ancestry segment trees with reference counting
   - Memory-mapped sweep trajectories
   - Node recycling during simulation
   - AVL tree indexing for O(log n) lookups

3. **Performance Improvements**
   - Mutation handling: Up to 50x speedup for high θ
   - Per-population node lists: O(1) selection
   - discoal runs ~28x faster than msprime
   - High recombination scenarios now feasible

4. **Bug Fixes**
   - Zero memory leaks (verified with valgrind)
   - Fixed uninitialized memory access
   - Resolved use-after-free errors
   - Command line parsing safety improvements

### Key Implementation Details
- Tree sequences use coalescent time units (2N generations)
- Mutation rate includes 0.5 factor for time scaling
- Minimal mode: recombination parents don't get tskit nodes
- Sample node IDs updated after simplification for MS output
- All pointer arrays use `calloc` for initialization

## Project Structure

### Core Files
- `discoal_multipop.c` - Main entry point, command parsing, simulation loop
- `discoalFunctions.c` - Core simulation functions
- `alleleTraj.c` - Allele trajectory calculations for sweeps
- `tskit_interface.c` - Tree sequence recording and output
- `discoal.h` - Main header with data structures

### Key Data Structures
- **rootedNode** - Coalescent tree node with:
  - Ancestry segment tree (`ancestryRoot`)
  - Population tracking
  - tskit node ID mapping
  - Memory recycling flags
  
- **AncestrySegment** - Interval-based ancestry tracking:
  - Reference counting for sharing
  - AVL tree support for fast lookups
  - tskit node ID field

### Directories
- `/` - Source files, Makefile, executables
- `testing/` - Test scripts (NEVER DELETE)
- `test/unit/` - Unit tests (77 tests total)
- `examples/` - Example usage scripts
- `docs/` - Documentation

## Build & Test

### Build Commands
```bash
make discoal         # Build main executable
make clean          # Remove build artifacts
make run_tests      # Run all unit tests
make test_<name>    # Run specific unit test
```

### Testing
1. **Quick validation**: `cd testing && ./focused_validation_suite.sh`
2. **Full validation**: `cd testing && ./comprehensive_validation_suite.sh`
3. **Statistical tests**: `cd testing && ./statistical_validation_suite.sh`
4. **msprime comparison**: `cd testing && ./msprime_comparison_suite.sh`

## Development Workflow

1. Always run tests before commits
2. Use git branches for features
3. Don't add Claude as co-author
4. Document memory/performance improvements in commits

### Recent Changes (2024-12-20)
- Removed conditional recombination mode (-C option) as undocumented legacy code
- Fixed -ej option (ms-compatible syntax) for population joining events
- Cleaned up mutation function redundancy:
  - Removed empty stub functions (tskit_record_mutations, tskit_populate_discoal_mutations)
  - Removed redundant wrapper (tskit_place_mutations_directly)
  - Consolidated makeGametesMS by removing unnecessary wrapper
- Fixed comprehensive validation suite (-es option test removed, speedup calculation fixed)
- Cleaned up ancestry segment code:
  - Fixed critical AVL tree threshold bug (was 30, should be 3)
  - Removed unused splitResult struct and static functions
  - Removed debug code and associated variables
  - Improved code maintainability

### Active Development
- [ ] Document memory optimization techniques in README
- [ ] Phase 4: Memory layout optimizations for cache efficiency
- [ ] Investigate periodic tskit simplification

### Memory Management
- Dynamic arrays with capacity tracking
- Reference counting for segment sharing
- Memory-mapped files for trajectories
- Node recycling with mark-and-sweep GC