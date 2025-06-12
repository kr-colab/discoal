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

Basic unit tests for core data structures:

.. code-block:: bash

   make run_tests

This runs tests for:

* Node operations
* Event handling  
* Mutation tracking
* Node manipulation functions

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
3. **Run focused tests** frequently during development:

   .. code-block:: bash

      cd testing/ && ./focused_validation_suite.sh

4. **Run comprehensive tests** before committing:

   .. code-block:: bash

      cd testing/ && ./comprehensive_validation_suite.sh

5. **Document performance improvements** in commit messages
6. **Submit pull request** with test results

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