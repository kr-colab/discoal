Selection
=========

discoal can simulate various types of natural selection, including hard sweeps, soft sweeps, partial sweeps, and recurrent hitchhiking. Selection is modeled by conditioning the coalescent on allele frequency trajectories.

Types of Sweeps
---------------

Hard Sweeps
^^^^^^^^^^^

A classical selective sweep from a single new mutation:

* ``-wd``: Deterministic trajectory
* ``-ws``: Stochastic trajectory  
* ``-wn``: Neutral fixation

**Parameters:**

* ``tau``: Time since fixation (in units of 4N generations)
* ``-a alpha``: Selection strength (2Ns)
* ``-x position``: Location of selected site (0-1)

**Example:**

.. code-block:: bash

   # Stochastic sweep fixed 0.05 time units ago
   ./discoal 20 10 10000 -t 20 -r 20 -ws 0.05 -a 1000 -x 0.5

Soft Sweeps
^^^^^^^^^^^

Selection on standing variation or recurrent mutation.

**From standing variation** (``-f``):

.. code-block:: bash

   # Sweep from initial frequency 0.1
   ./discoal 20 10 10000 -t 20 -r 20 -ws 0.05 -a 1000 -x 0.5 -f 0.1

**From recurrent mutation** (``-uA``):

.. code-block:: bash

   # Recurrent adaptive mutation rate 0.0001
   ./discoal 20 10 10000 -t 20 -r 20 -ws 0.05 -a 1000 -x 0.5 -uA 0.0001

Partial Sweeps
^^^^^^^^^^^^^^^

Sweeps that stop before fixation using ``-c``:

.. code-block:: bash

   # Partial sweep to frequency 0.8
   ./discoal 20 10 10000 -t 20 -r 20 -ws 0.05 -a 1000 -x 0.5 -c 0.8
   
   # Partial soft sweep from standing variation
   ./discoal 20 10 10000 -t 20 -r 20 -ws 0.05 -a 1000 -x 0.5 -f 0.1 -c 0.8

Linked Selection
----------------

Simulate effects of a sweep outside the sampled region:

* ``-ls``: Stochastic sweep to the left
* ``-ld``: Deterministic sweep to the left  
* ``-ln``: Neutral fixation to the left

**Example:**

.. code-block:: bash

   # Sweep to the left at genetic distance 4Nr = 100
   ./discoal 20 10 10000 -t 20 -r 20 -ls 0.05 -a 1000

The position is drawn uniformly between 0 and the specified genetic distance.

Recurrent Hitchhiking
---------------------

Multiple sweeps over time:

* ``-R rate``: Sweeps within the locus
* ``-L rate``: Sweeps to the left of the locus

Rate is per 2N individuals per generation.

**Example:**

.. code-block:: bash

   # Recurrent sweeps at rate 0.001
   ./discoal 20 10 10000 -t 20 -r 20 -R 0.001 -a 1000
   
   # Recurrent partial sweeps
   ./discoal 20 10 10000 -t 20 -r 20 -R 0.001 -a 1000 -c 0.7

Combining with Demography
-------------------------

Selection can be combined with population size changes:

.. code-block:: bash

   # Sweep during population expansion
   ./discoal 20 10 10000 -t 20 -r 20 \
     -ws 0.05 -a 1000 -x 0.5 \
     -en 0.02 0 10.0 \
     -en 0.1 0 1.0

Advanced Options
----------------

**Time discretization** (``-i``):

Controls time step size during sweeps (default 40):

.. code-block:: bash

   # Finer time steps (slower but more accurate)
   ./discoal 20 10 10000 -t 20 -ws 0.05 -a 1000 -i 400

**Effective population size during sweep** (``-N``):

.. code-block:: bash

   # Set sweep effective size to 500,000
   ./discoal 20 10 10000 -t 20 -ws 0.05 -a 1000 -N 500000

Output Considerations
---------------------

By default, the selected SNP is included in the output. To exclude it:

.. code-block:: bash

   # Hide selected SNP
   ./discoal 20 10 10000 -t 20 -ws 0.05 -a 1000 -x 0.5 -h

.. warning::
   
   Including the selected SNP will bias estimates of θ. Simulations with θ=0 will still contain the selected SNP.