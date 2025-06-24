# Discoal Project Notes

## Latest Changes and Optimizations

### 2025-06-24: YAML Configuration System Debugging

**Objective**: Fix issues preventing YAML configurations from producing identical output to command line arguments.

**Issues Fixed**:

1. **Seed Array Parsing Bug**:
   - Problem: YAML parser was incorrectly parsing seed arrays `[12345, 67890]` as key-value pairs
   - Solution: Added special handling in `parseSimulationSection()` to detect numeric keys as seed array elements
   - Result: Seeds now correctly parsed from YAML arrays

2. **Parameter Initialization Order Bug**:
   - Problem: Variable initialization in `getParameters()` was happening AFTER YAML configuration was applied
   - Solution: Moved ALL variable initialization (seeds, theta, rho, etc.) to the beginning of `getParameters()`, before YAML processing
   - Result: YAML values are no longer overwritten by defaults

3. **Gene Conversion Crossover Ratio Mode**:
   - Problem: YAML crossover ratio mode wasn't matching command line `-gr` behavior
   - Solution: Modified `applyConfiguration()` to prioritize crossover ratio mode and skip setting gamma directly
   - Result: Gene conversion rate correctly calculated as `rho * gammaCoRatio` when using crossover ratio

**Validation**:
- Created comprehensive test suite `testing/yaml_validation_suite.sh` with 10 test cases
- All tests now pass with 100% success rate
- YAML configurations produce identical output to equivalent command line arguments

**Files Modified**:
- `src/core/configInterface.c`: Fixed seed parsing and crossover ratio application
- `src/core/discoal_multipop.c`: Fixed initialization order in `getParameters()`
- `testing/yaml_validation_suite.sh`: Created comprehensive test suite

### 2025-06-23: Edge Recording Refactoring

**Objective**: Refactor duplicated edge recording logic into reusable convenience functions.

**Changes Made**:

1. **Created `recordEdgesForChild()` function**:
   - Handles edge recording from parent to child based on ancestry segments
   - Supports both full ARG mode and minimal tree sequence mode
   - Eliminates duplicated segment traversal logic

2. **Created `recordCoalescentEdges()` function**:
   - Orchestrates edge recording for coalescent events (one parent, two children)
   - Handles child ordering requirements
   - Ensures parent segments are marked as recorded
   - Used in both `coalesceAtTimePopn()` and `coalesceAtTimePopnSweep()`

3. **Created `recordRecombinationEdges()` function**:
   - Handles edge recording for recombination and gene conversion events (one child, two parents)
   - Only operates in full ARG mode (skips minimal mode)
   - Consolidates logic used in 4 different locations:
     - Recombination events
     - Gene conversion events  
     - Recombination during selective sweeps
     - Gene conversion during selective sweeps

**Benefits**:
- Eliminated ~150 lines of duplicated code
- Centralized edge recording logic for easier maintenance
- Improved code readability and organization
- Maintained all existing functionality and performance optimizations

**Validation**:
- All 30 comprehensive test cases pass (100% success rate)
- No performance regressions
- Memory usage and speed improvements preserved

### 2025-06-24: Demes-C API Integration

**Objective**: Integrate demes-c API into discoal to import demographic events from YAML files.

**Changes Made**:

1. **Created demesInterface module**:
   - `demesInterface.h` and `demesInterface.c` for demes integration
   - Added `-D` flag for loading demes YAML files
   - Integrated demes-c library into build system via Makefile

2. **Implemented parameter conversions**:
   - Time scaling: generations / (2N) to match discoal's internal units
   - Migration rate scaling: demes_rate * 4N to get 4Nm
   - Population size scaling: relative to reference N (first present-day population)
   - Reference population selection: uses first population that exists at time 0

3. **Fixed critical issues**:
   - Asymmetric migration: Changed from bidirectional to unidirectional event creation
   - Population size changes: Use previous epoch size for transitions
   - Time scaling: Accounted for discoal's 2x multiplier in command line parsing
   - Population splits: Changed event type from 'S' to 'p' for compatibility

4. **Added comprehensive error handling**:
   - Exponential growth epochs (not supported)
   - Linear growth epochs (not supported)
   - Selfing and cloning rates (not supported)
   - Disconnected populations validation (infinite coalescent time)

5. **Created validation test suite**:
   - `demes_output_comparison_suite.sh` with 5 comprehensive tests
   - All tests pass - demes and command line produce identical outputs
   - Added detailed documentation about time conversion complexities

**Key Insights**:
- Discoal time units are confusing: command line uses 2N, internal uses 4N
- Population indexing matters: demes populations must map correctly to discoal indices
- Event order and timing must match exactly for reproducibility

**Benefits**:
- Users can now specify complex demography via standard demes YAML format
- Maintains backward compatibility with command line specification
- Enables interoperability with other tools using demes format

**Technical Details**:
- Time conversion: Demes generations → Internal time = generations/(2N)
- Command line time gets multiplied by 2 during parsing
- Example: 200k generations = 0.1 internal = 0.05 command line (with N=1M)