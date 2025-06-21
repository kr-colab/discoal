Population Structure
====================

Multiple Populations
--------------------

discoal can simulate multiple populations with migration and complex demographic histories using the ``-p`` flag.

**Basic syntax:**

.. code-block:: bash

   ./discoal totalSampleSize numReps nSites -t theta -p numPops size1 size2 ...

**Example:** Three populations with 2 samples each:

.. code-block:: bash

   ./discoal 6 2 100 -t 2 -r 2.4 -p 3 2 2 2

.. note::
   
   Populations are zero-indexed (0, 1, 2, ...). A population can have zero samples, which is useful for ghost populations.

Migration
---------

Symmetric Migration
^^^^^^^^^^^^^^^^^^^

Set all pairwise migration rates to the same value with ``-M``:

.. code-block:: bash

   # Island model with migration rate 4Nm = 0.05
   ./discoal 6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -M 0.05

Asymmetric Migration
^^^^^^^^^^^^^^^^^^^^

Set specific migration rates between population pairs with ``-m``:

.. code-block:: bash

   # Migration from pop 0 to pop 1 at rate 0.1
   # Migration from pop 1 to pop 0 at rate 0.05
   ./discoal 4 2 100 -t 2 -p 2 2 2 -m 0 1 0.1 -m 1 0 0.05

Population Splits
-----------------

Model population divergence events using ``-ed`` (forward in time, populations merge backward in time):

**Example:** ((pop0, pop1), pop2) topology

.. code-block:: bash

   # pop0 and pop1 split at time 1.0
   # Their ancestor and pop2 split at time 5.0
   ./discoal 6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -ed 1.0 0 1 -ed 5.0 1 2

Admixture Events
----------------

Model admixture using ``-ea`` where one population derives ancestry from two source populations:

**Syntax:** ``-ea time admixedPop sourcePop1 sourcePop2 proportion``

**Example:**

.. code-block:: bash

   # At time 0.02, pop0 derives 15% ancestry from pop2, 85% from pop1
   ./discoal 6 2 100 -t 2 -r 2.4 -p 3 2 2 2 -ea 0.02 0 1 2 0.85

Population Size Changes
-----------------------

Change population sizes at specific times with ``-en``:

.. code-block:: bash

   # Pop 0: bottleneck to 10% size at time 0.5
   # Pop 1: expansion to 5x size at time 0.3
   ./discoal 4 2 1000 -t 10 -p 2 2 2 -en 0.5 0 0.1 -en 0.3 1 5.0

Complex Demographic Models
--------------------------

Combine multiple demographic events:

.. code-block:: bash

   # Three population model with:
   # - Migration between sister populations
   # - Recent admixture event  
   # - Population size changes
   # - Ancient population splits
   
   ./discoal 30 10 100000 -t 100 -r 100 \
     -p 3 10 10 10 \
     -m 0 1 0.5 -m 1 0 0.5 \
     -ea 0.01 0 0 2 0.95 \
     -en 0.05 0 0.1 -en 0.2 0 1.0 \
     -en 0.05 1 0.1 -en 0.2 1 1.0 \
     -ed 0.5 0 1 \
     -ed 2.0 1 2

Ancient Samples
---------------

Sample lineages from different time points using ``-A``:

.. code-block:: bash

   # 5 modern samples from pop 0
   # 3 ancient samples from pop 0 at time 0.1
   # 2 ancient samples from pop 1 at time 0.5
   ./discoal 10 5 10000 -t 20 -p 2 5 0 \
     -A 3 0 0.1 \
     -A 2 1 0.5

Limitations
-----------

.. warning::
   
   * Selection only operates in population 0
   * During selective sweeps, migration is suspended
   * Time-varying migration rates are not currently implemented