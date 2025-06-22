Development
===========

This section covers the development environment setup and testing procedures for contributing to discoal.

Recent Improvements
-------------------

Memory Optimizations
^^^^^^^^^^^^^^^^^^^^

The codebase has undergone significant memory optimizations:

* **Ancestry tracking**: Replaced fixed-size ``ancSites`` arrays with dynamic ancestry segment trees
* **Sample size limit increased**: Maximum sample size increased from 254 to 65,535
* **Memory usage reduced**: 15-99% memory savings across different simulation scenarios
* **Dynamic allocation**: All major data structures now use dynamic memory allocation

These changes eliminate the need for the ``-DBIG`` compilation flag, which is now obsolete.

Development Environment
-----------------------

We use conda to manage the development environment. This ensures all developers have consistent tools and dependencies.

Setting Up the Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Install conda or miniforge if you haven't already
2. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/kern-lab/discoal.git
      cd discoal

3. Create the conda environment:

   .. code-block:: bash

      conda env create -f environment.yml

4. Activate the environment:

   .. code-block:: bash

      conda activate discoal_dev

The environment includes:

* Python 3.12
* msprime (for comparison testing)
* Sphinx and sphinx-rtd-theme (for documentation)

Building discoal
----------------

With the development environment activated:

.. code-block:: bash

   make clean
   make discoal

For testing, you'll need both optimized and legacy versions:

.. code-block:: bash

   make test_binaries

This builds:

* ``discoal_edited``: The optimized version with all memory improvements
* ``discoal_legacy_backup``: The reference version from the master-backup branch

Testing
-------

discoal has a comprehensive testing framework to ensure code changes maintain correctness and performance.

Unit Tests
^^^^^^^^^^

discoal has a comprehensive unit testing framework using the Unity test framework. The unit tests cover all major components of the codebase.

**Running Unit Tests**

To run all unit tests:

.. code-block:: bash

   make run_tests

To run a specific test suite:

.. code-block:: bash

   make test_node         # Test node operations
   make test_event        # Test event handling
   # ... etc

**Test Coverage**

The unit test suite includes 49 tests across 7 test files:

1. **Node Operations** (``test_node.c`` - 3 tests):
   
   * Node initialization and property setting
   * Creation of new rooted nodes
   * Basic node structure validation

2. **Event Handling** (``test_event.c`` - 2 tests):
   
   * Event structure initialization
   * Event property manipulation

3. **Node Operations** (``test_node_operations.c`` - 4 tests):
   
   * Creating and destroying nodes
   * Adding and removing nodes from active set
   * Node selection by population
   * Population size tracking

4. **Node Fields** (``test_mutations.c`` - 3 tests):
   
   * Basic node creation and field validation
   * Node state tracking (parentsRecorded, isFullyRecorded, etc.)
   * Ancestry-related fields testing

5. **Ancestry Segment Trees** (``test_ancestry_segment.c`` - 13 tests):
   
   * Segment creation and validation
   * Reference counting (retain/release)
   * Shallow vs deep copying
   * Tree merging and splitting operations
   * Ancestry count queries
   * NULL safety checks

6. **Active Material Segments** (``test_active_segment.c`` - 12 tests):
   
   * Active material initialization
   * Site activity queries
   * Fixed region removal
   * Segment coalescing
   * AVL tree integration
   * Verification functions

7. **Trajectory Handling** (``test_trajectory.c`` - 12 tests):
   
   * Trajectory capacity management
   * File cleanup for rejected trajectories
   * Memory-mapped file operations
   * Large file handling
   * File persistence and cleanup
   * Concurrent trajectory management

**Note**: Two test suites (``test_coalescence_recombination.c`` and ``test_memory_management.c``) were removed as they referenced obsolete data structures that have been replaced in the tskit integration.

**Building Individual Tests**

Each test suite can be built separately:

.. code-block:: bash

   make test_ancestry_segment
   make test_node_operations
   # etc.

**Test Development**

When adding new functionality:

1. Create a new test file in ``test/unit/`` following the naming convention ``test_<component>.c``
2. Include the Unity framework headers
3. Write setUp() and tearDown() functions for test fixtures
4. Add test functions following the pattern ``test_<functionality>_<scenario>()``
5. Update the Makefile with build rules for the new test

**Test Infrastructure**

* Unity test framework is automatically downloaded to ``extern/Unity/`` on first build
* ``test_globals.c`` provides minimal global variables (seed1, seed2, currentSize) needed for tests
* Tests compile with actual discoal source code (no mocks) for integration testing

**Debugging Tests**

To debug a failing test:

.. code-block:: bash

   # Build with debug symbols
   gcc -g -O0 -I. -I./test/unit -o test_name test/unit/test_name.c test/unit/unity.c \
       discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c \
       ancestrySegmentAVL.c ancestryVerify.c activeSegment.c -lm -fcommon
   
   # Run with gdb
   gdb ./test_name

**Unity Test Framework**

The tests use the official Unity test framework (https://github.com/ThrowTheSwitch/Unity) which provides:

* Rich assertion macros (TEST_ASSERT_EQUAL, TEST_ASSERT_FLOAT_WITHIN, etc.)
* Automatic test discovery and execution
* Clear failure messages with file and line information
* Test fixtures with setUp/tearDown support

The framework files are located in ``test/unit/``:
* ``unity.h`` - Main header file
* ``unity.c`` - Implementation
* ``unity_internals.h`` - Internal definitions

**Quick Testing Reference**

Common testing commands during development:

.. code-block:: bash

   # Run all unit tests
   make run_tests
   
   # Run specific test suite
   make test_node_operations && ./test_node_operations
   
   # Clean and rebuild all tests
   make clean && make run_tests
   
   # Quick validation during development
   cd testing/ && ./focused_validation_suite.sh
   
   # Full validation before commits
   cd testing/ && ./comprehensive_validation_suite.sh
   
   # Statistical validation (if needed)
   cd testing/ && ./statistical_validation_suite.sh
   
   # Run comprehensive tests (optimized vs legacy from master-backup)
   make test_comprehensive
   
   # Run comprehensive tests (current working dir vs HEAD of branch)
   make test_comprehensive_head

**Make Targets for Comprehensive Testing**

The Makefile provides convenient targets that build the required binaries and run the comprehensive test suite:

* ``make test_comprehensive``: 
  
  * Builds ``discoal_edited`` (optimized version from current working directory)
  * Builds ``discoal_legacy_backup`` from the ``master-backup`` branch
  * Runs the comprehensive validation suite comparing these two versions
  * Use this to ensure your optimizations maintain compatibility with the original implementation

* ``make test_comprehensive_head``:
  
  * Builds ``discoal_edited`` (optimized version from current working directory)
  * Builds ``discoal_legacy_backup`` from HEAD of the current branch
  * Runs the comprehensive validation suite comparing working changes against the last commit
  * Use this to measure performance improvements of uncommitted changes

These targets automatically handle the complex process of building from different sources and are the recommended way to run comprehensive tests during development.

Comprehensive Validation Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary testing framework compares the optimized version against the legacy version to ensure identical output:

.. code-block:: bash

   cd testing/
   ./comprehensive_validation_suite.sh

This suite:

* Runs 27 test cases covering all documented features
* Compares output between optimized and legacy versions
* Profiles memory usage and performance
* Reports any differences or regressions

Test categories include:

* Basic coalescent simulations
* Recombination and gene conversion
* Multiple populations with migration
* Selection (hard/soft/partial sweeps)
* Complex demographic scenarios
* Stress tests with extreme parameters

Focused Validation Suite
^^^^^^^^^^^^^^^^^^^^^^^^

For rapid testing during development:

.. code-block:: bash

   cd testing/
   ./focused_validation_suite.sh

This runs a subset of critical tests for quick feedback.

Statistical Validation Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To ensure optimizations don't introduce statistical biases:

.. code-block:: bash

   cd testing/
   ./statistical_validation_suite.sh              # 100 replicates, auto mode
   ./statistical_validation_suite.sh parallel 50  # 50 replicates, parallel mode
   ./statistical_validation_suite.sh 200          # 200 replicates

This suite:

* Runs multiple replicates of each test case
* Extracts segregating sites statistics
* Performs Kolmogorov-Smirnov tests
* Verifies distributions are statistically equivalent

msprime Comparison Suite
^^^^^^^^^^^^^^^^^^^^^^^^

To validate discoal against the well-established msprime coalescent simulator:

.. code-block:: bash

   cd testing/
   ./msprime_comparison_suite.sh

This suite compares discoal and msprime across:

* Neutral models with and without recombination
* Various sample sizes and mutation rates
* Selection models (hard sweeps with different strengths and ages)

The comparison includes runtime performance metrics and statistical tests to ensure equivalent output distributions.

**Parameter Scaling for msprime Comparisons**

When comparing discoal with msprime, careful parameter conversion is required due to different conventions:

1. **Population Size**: discoal uses scaled parameters assuming Ne=1. For msprime, we use Ne=0.5 with diploid samples (n_samples/2) and ploidy=2 to match discoal's haploid output.

2. **Mutation Rate**: 
   
   * discoal: θ = 4 × Ne × μ × L (over whole locus)
   * msprime: mutation_rate = θ / (4 × Ne × L) (per base pair)

3. **Recombination Rate**:
   
   * discoal: ρ = 4 × Ne × r × L
   * msprime: recombination_rate = ρ / (4 × Ne × L)

4. **Selection Coefficient** (for sweeps):
   
   * discoal: α = 2 × Ne × s
   * msprime: s = α / (2 × Ne) × 2 (factor of 2 for msprime's fitness model)

5. **Sweep Timing**:
   
   * When τ > 0 in discoal, we rescale to Ne=0.25 in msprime for consistent time units
   * Allele frequencies use the original Ne to ensure valid [0,1] bounds

These scaling conventions ensure that both simulators produce statistically equivalent results, as validated by the comparison suite.

Development Workflow
--------------------

1. **Create a feature branch** from the main development branch
2. **Make changes** to the code
3. **Run unit tests** during development:

   .. code-block:: bash

      make run_tests

4. **Run focused validation tests** frequently:

   .. code-block:: bash

      cd testing/ && ./focused_validation_suite.sh

5. **Run comprehensive tests** before committing:

   .. code-block:: bash

      cd testing/ && ./comprehensive_validation_suite.sh

6. **Document performance improvements** in commit messages
7. **Submit pull request** - CI will automatically run tests

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^

discoal uses GitHub Actions for continuous integration. The CI workflow automatically runs on every push and pull request to ensure code quality.

**CI Jobs:**

1. **Unit Tests**: Runs all unit tests on multiple OS versions (Ubuntu 20.04, 22.04, latest) with both gcc and clang compilers
2. **Basic Functionality Tests**: Verifies core discoal functionality (coalescent, recombination, selection, tree sequence output)
3. **Memory Leak Check**: Uses valgrind to ensure no memory leaks in common usage scenarios

**CI Badge:**

The build status is displayed at the top of the README. Green indicates all tests are passing.

To run the same tests locally before pushing:

.. code-block:: bash

   # Run unit tests
   make run_tests
   
   # Check for memory leaks (requires valgrind)
   valgrind --leak-check=full ./discoal 10 1 100 -t 10

Code Organization
-----------------

Key source files:

* ``discoal_multipop.c``: Main program entry and command-line parsing
* ``discoalFunctions.c``: Core simulation functions
* ``alleleTraj.c``: Allele trajectory calculations for sweeps
* ``ancestrySegment.c``: Memory-efficient ancestry tracking
* ``activeSegment.c``: Active material tracking
* ``discoal.h``: Main header with data structures

Memory Optimizations
--------------------

Recent optimizations have achieved significant memory reductions:

* Dynamic allocation for all major arrays
* Segment trees for ancestry tracking (80% reduction)
* Reference counting for segment sharing (10-16% additional reduction)
* AVL tree indexing for high-recombination scenarios
* Memory-mapped files for sweep trajectories

When developing, maintain these optimizations and ensure new features don't regress memory usage.

Documentation
-------------

To build the documentation locally:

.. code-block:: bash

   cd docs/
   make html

View the built documentation:

.. code-block:: bash

   open _build/html/index.html  # macOS
   xdg-open _build/html/index.html  # Linux

Before submitting changes, ensure documentation is updated for any new features or parameter changes.