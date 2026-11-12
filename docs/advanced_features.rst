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

Each prior needs ``low <= high`` inside the parameter's range, as the
YAML priors do: ``-Pt``, ``-Pr``, ``-Pa``, ``-Pu`` and ``-PuA`` need
``0 <= low``; ``-Px`` needs ``high <= 1``; ``-Pf`` and ``-Pc`` need
``high < 1``. ``-Pre`` needs a mean and an upper bound above 0.
discoal stops with an error otherwise. The fixed values follow the same
rules (``-t``, ``-r``, ``-g``, ``-a``, ``-uA`` and the sweep time
``>= 0``, ``-x`` in [0, 1], ``-f`` in [0, 1), ``-i`` above 0).

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

Fixed Segregating Sites and Conditional Simulation
--------------------------------------------------

discoal does not have ms's ``-s`` (a fixed number of segregating sites),
and it no longer has ``-C`` (keep only replicates with recombination in
a region); see :doc:`changes`. Both are simple to do with tskit on tree
sequence output (``-ts``): for example, place a set number of mutations
on the genealogy, or keep only the replicates whose trees change in a
region.

Sweep Trajectory From a File
----------------------------

With ``-S file``, the sweep follows a trajectory from a file instead of
one that discoal simulates, for example a trajectory from a forward
simulation or an estimate from data:

.. code-block:: text

   # time (4N0 generations before the present)  frequency
   0       0.6
   0.0005  0.476
   0.001   0.354
   ...
   0.015   0

.. code-block:: bash

   ./discoal 20 100 10000 -t 10 -r 10 -N 10000 -S trajectory.txt -x 0.3

The sweep starts (going back in time) at the first time, here the
present. The last frequency must be 0, the allele's origin. Between the
rows, discoal interpolates the frequency linearly at each of its sweep
time steps, of ``1/(i N)`` in units of 2N₀ generations (``-i``, default
40, and ``-N``). As with a simulated sweep, the sweep ends at the first
step at which the frequency is at or below 1/(2N), so rows after that
point have no effect, and frequencies above 1 - 1/(2N) are used as
1 - 1/(2N).
If the first frequency is below 1, as here, the sweep is partial, as
with ``-c``: each sampled lineage carries the allele with that
probability, and the selected SNP is in the output. The file gives the
sweep's time and strength, so ``-S`` cannot be combined with the sweep
options that set them (``-ws``, ``-wd``, ``-wn``, ``-ls``, ``-ld``,
``-ln``, ``-R``, ``-L``, ``-a``, ``-f``, ``-c``, and their priors).
``-x``, ``-Px``, ``-uA``, ``-h`` and ancient samples work as
with a simulated trajectory.

Tree Sequence Output
--------------------

``-ts file.trees`` writes the genealogies and mutations of each replicate
as a tskit tree sequence, in place of the ms-format haplotypes:

.. code-block:: bash

   ./discoal 10 3 10000 -t 20 -r 10 -ts sim.trees

With one replicate the file is ``sim.trees``; with more, discoal writes
``sim_rep1.trees``, ``sim_rep2.trees``, and so on. Load them with tskit:

.. code-block:: python

   import tskit
   ts = tskit.load("sim_rep1.trees")
   print(ts.num_trees, ts.diversity(mode="branch"))

As in msprime, discoal records only the nodes and edges where lineages
coalesce. ``-F`` (given after ``-ts``) records the full ancestral
recombination graph: every node and edge, also of recombination events.
Node times are in units of 2N0 generations. tskit can write the trees in
Newick format (``ts.trees()`` and ``tree.as_newick()``).

Memory and Performance
----------------------

Large-scale Simulations
^^^^^^^^^^^^^^^^^^^^^^^

For very large simulations:

1. **Number of sites**: at most 10^8 (``MAXSITES`` in ``src/core/discoal.h``)
2. **Sample size**: at most 65,535

Performance Tuning
^^^^^^^^^^^^^^^^^^

* **Time discretization**: Lower ``-i`` values speed up sweeps at potential accuracy cost
* **High recombination rates**: discoal records tree sequences as msprime
  does, so memory stays at tens of megabytes and run time grows slowly
  with rho
* **Parallel runs**: Use different random seeds for embarrassingly parallel execution

Setting Random Seeds
--------------------

For reproducibility:

.. code-block:: bash

   # Specify both seeds
   ./discoal 20 10 10000 -t 20 -d 12345 67890

The first seed must be in [1, 2147483562] and the second in
[1, 2147483398]; discoal stops with an error otherwise.

Debugging Features
------------------

Build with ancestry verification:

.. code-block:: bash

   make discoal_debug
   ./build/discoal_debug 10 1 1000 -t 10 -r 10

This builds the same program with debug symbols and no optimization, for
gdb and valgrind.