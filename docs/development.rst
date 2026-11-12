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

      git clone https://github.com/kr-colab/discoal.git
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

The unit tests have 219 tests in 11 suites:

==============================  =====  ============================================================
suite                           tests  what it tests
==============================  =====  ============================================================
``test_node.c``                 3      node creation and fields
``test_event.c``                3      demographic event records
``test_node_operations.c``      4      adding and removing nodes, choice by population, sizes
``test_mutations.c``            3      node fields
``test_ancestry_segment.c``     13     segment lists: copying, reference counts, merges, splits
``test_active_segment.c``       12     sites without an MRCA: removal of fixed regions, queries
``test_fenwick.c``              4      the Fenwick tree of recombination weights
``test_shapes.c``               64     size and migration shapes: constant, exponential, linear
``test_alleleTraj.c``           5      sweep trajectories
``test_trajectory.c``           6      trajectory files: capacity, cleanup
``test_config_interface.c``     102    the YAML parser and its checks
==============================  =====  ============================================================

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

* The Unity test framework is in ``extern/Unity/``
* ``test_globals.c`` provides minimal global variables (seed1, seed2, currentSize) needed for tests
* Tests compile with actual discoal source code (no mocks) for integration testing

**Debugging Tests**

To debug a failing test:

.. code-block:: bash

   # Build one test with debug symbols and no optimization
   make test_ancestry_segment TEST_CFLAGS='-g -O0 $(INCLUDE)'
   
   # Run with gdb
   gdb ./build/test_ancestry_segment

**Profiling**

``make discoal_profile`` builds ``build/discoal_profile``: the release
build (the same optimization flags) with debug symbols and frame pointers,
so that profilers see the code that users run and get full call stacks:

.. code-block:: bash

   make discoal_profile
   perf record --call-graph fp ./build/discoal_profile 20 1 10000 -t 10 -r 10000
   perf report
   valgrind --tool=callgrind ./build/discoal_profile 20 1 10000 -t 10 -r 1000

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
   make test_node_operations && ./build/test_node_operations

   # Clean and rebuild all tests
   make clean && make run_tests

   # Everything that CI checks for a pull request
   make check

Checks of the simulator
^^^^^^^^^^^^^^^^^^^^^^^

The ``testing/`` directory holds the checks of discoal's output: command
lines that must run, YAML and demes input that must give the same output as
the command line, and statistical comparisons with msprime (units,
demography, linkage, sweeps, gene conversion) and with SLiM (sweeps).
``testing/README.md`` lists each check, what it checks, and where it runs
(every pull request, weekly, or by hand).

``make check`` runs the checks of every pull request: the unit tests, then
each script in the same order as CI. ``CHECK_WORKERS`` sets the number of
parallel workers of the statistical checks (default 4):

.. code-block:: bash

   make check CHECK_WORKERS=16

The statistical checks need Python with msprime, tskit, numpy, scipy, demes
and stdpopsim. They use fixed seeds, so their results do not change from run
to run.

Development Workflow
--------------------

1. **Create a feature branch** from the main development branch
2. **Make changes** to the code
3. **Run unit tests** during development:

   .. code-block:: bash

      make run_tests

4. **Run all checks** before a pull request:

   .. code-block:: bash

      make check

5. **Submit pull request** - CI runs the same checks

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^

discoal uses GitHub Actions for continuous integration. The CI workflow automatically runs on every push and pull request to ensure code quality.

**CI Jobs (every push and pull request):**

1. **Unit Tests**: all unit tests on Ubuntu 22.04 and the latest Ubuntu, with gcc and clang
2. **Basic Functionality Tests**: smoke tests, and the YAML, preflight and demes suites
3. **Unit Conventions vs msprime**: the statistical checks against msprime (units, sweeps, trajectory files, recurrent sweeps, gene conversion, linkage, a fuzz test of random models)
4. **Parity with msprime on demes models**: a fixed set of core models and a random sample of the others
5. **Memory Leak Check**: valgrind on common runs

A weekly workflow runs the full parity check against msprime, without and
with recombination.

**CI Badge:**

The build status is displayed at the top of the README. Green indicates all tests are passing.

To run the same tests locally before pushing:

.. code-block:: bash

   make check

   # Check for memory leaks (requires valgrind)
   valgrind --leak-check=full ./build/discoal 10 1 100 -t 10

Code Organization
-----------------

Key source files, in ``src/core/``:

* ``discoal_multipop.c``: Main program entry and command-line parsing
* ``discoalFunctions.c``: Core simulation functions
* ``alleleTraj.c``: Allele trajectory calculations for sweeps
* ``ancestrySegment.c``: Memory-efficient ancestry tracking
* ``activeSegment.c``: Active material tracking
* ``discoal.h``: Main header with data structures

Memory Optimizations
--------------------

These structures keep memory use small:

* Dynamic allocation for all major arrays
* Segment lists for ancestry tracking, in place of per-site arrays
* Reference counting for segment sharing
* In minimal tree sequence mode, recombination splits a lineage's segment
  list in place, and a coalescence builds the parent's list in one walk
  that also gives its span
* Sweep trajectories in temporary files, memory-mapped

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