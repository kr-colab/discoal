Advanced Features
=================

Parameter Priors
----------------

discoal supports uniform and exponential prior distributions for parameters, useful for ABC and machine learning applications.

Uniform Priors
^^^^^^^^^^^^^^

Most parameters support uniform priors specified with ``-P`` flags:

.. list-table:: Prior Distribution Flags
   :header-rows: 1
   :widths: 20 40 40

   * - Flag
     - Parameters
     - Description
   * - ``-Pt``
     - low high
     - Prior on θ (mutation rate)
   * - ``-Pr``
     - low high
     - Prior on ρ (recombination rate)
   * - ``-Pa``
     - low high  
     - Prior on α (selection coefficient)
   * - ``-Pu``
     - low high
     - Prior on τ (time since fixation)
   * - ``-PuA``
     - low high
     - Prior on uA (recurrent adaptive mutation)
   * - ``-Px``
     - low high
     - Prior on sweep position
   * - ``-Pf``
     - low high
     - Prior on f₀ (initial frequency)
   * - ``-Pc``
     - low high
     - Prior on final sweep frequency

**Example:**

.. code-block:: bash

   # Variable mutation and selection parameters
   ./discoal 20 1000 10000 -Pt 5 50 -Pa 100 2000 -ws 0.01 -x 0.5

Exponential Prior
^^^^^^^^^^^^^^^^^

Recombination rate can follow a truncated exponential:

.. code-block:: bash

   # Exponential with mean 10, truncated at 100
   ./discoal 20 100 10000 -t 20 -Pre 10 100

Demographic Priors
^^^^^^^^^^^^^^^^^^

Population size changes can have priors:

.. code-block:: bash

   # Prior on first size change: time 0.01-0.5, size 0.1-10
   ./discoal 20 100 10000 -t 20 -Pe1 0.01 0.5 0.1 10
   
   # Prior on second size change
   ./discoal 20 100 10000 -t 20 -Pe1 0.01 0.5 0.1 10 -Pe2 0.5 2.0 0.5 5.0

Conditional Simulations
-----------------------

Simulate conditional on observing recombination in a specific region:

.. code-block:: bash

   # Condition on recombination between sites 400-600
   ./discoal 20 100 1000 -t 10 -r 20 -C 400 600

The simulator will retry until this condition is met.

Tree Output Mode
----------------

Output genealogical trees in Newick format with ``-T``:

.. code-block:: bash

   ./discoal 10 1 10000 -t 20 -r 10 -T

Output format:

.. code-block:: text

   //
   [100](0:0.5,1:0.5);
   [50]((0:0.2,2:0.2):0.3,1:0.5);
   ...

Each line shows: ``[number_of_sites]newick_tree``

Memory and Performance
----------------------

Large-scale Simulations
^^^^^^^^^^^^^^^^^^^^^^^

For very large simulations:

1. **Increase MAXSITES**: Edit ``discoal.h`` and recompile
2. **Sample size limit**: The maximum sample size is now 65,535 (previously 254)

Performance Tuning
^^^^^^^^^^^^^^^^^^

* **Time discretization**: Lower ``-i`` values speed up sweeps at potential accuracy cost
* **Memory efficiency**: Current version uses 70-99% less memory than older versions
* **Parallel runs**: Use different random seeds for embarrassingly parallel execution

Setting Random Seeds
--------------------

For reproducibility:

.. code-block:: bash

   # Specify both seeds
   ./discoal 20 10 10000 -t 20 -d 12345 67890

Seeds must be positive integers less than 2^31-1.

Debugging Features
------------------

Build with ancestry verification:

.. code-block:: bash

   make discoal_debug
   ./discoal_debug 10 1 1000 -t 10 -r 10

This enables additional checks for debugging genealogy construction.