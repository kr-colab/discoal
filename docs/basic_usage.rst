Basic Usage
===========

Command Line Structure
----------------------

At its most basic, a discoal command line looks like:

.. code-block:: bash

   ./discoal sampleSize numReplicates nSites -t theta

Where:

* ``sampleSize``: The number of chromosomes to sample (maximum: 65,535)
* ``numReplicates``: The number of independent samples to generate
* ``nSites``: The number of discrete sites in the sequence
* ``theta``: The population mutation rate (4N₀μ) where N₀ is the current population size and μ is the mutation rate per generation for the entire locus

Example
-------

.. code-block:: bash

   $ ./discoal 3 2 100 -t 2
   ./discoal 3 2 100 -t 2
   1665047201 686400060

   //
   segsites: 4
   positions: 0.01835 0.09557 0.46556 0.72880
   1000
   0110
   0111

   //
   segsites: 1
   positions: 0.07594
   0
   1
   0

Output Format
-------------

The output follows Hudson's ms format:

1. **Header**: Command line and random number seeds
2. **Sample blocks**: Each starting with ``//``
3. **Segregating sites**: Number of variable sites
4. **Positions**: Relative positions along the sequence (0-1)
5. **Haplotypes**: Binary strings where 0 is ancestral and 1 is derived

This format is compatible with existing ms analysis software.

Common Parameters
-----------------

Basic mutation and recombination:

.. code-block:: bash

   # With recombination rate rho = 4Nr
   ./discoal 10 5 10000 -t 10 -r 20

   # With gene conversion
   ./discoal 10 5 10000 -t 10 -r 20 -g 10 300

Population size changes:

.. code-block:: bash

   # Bottleneck: size drops to 10% at time 0.5, recovers to 80% at time 1.2
   ./discoal 10 5 10000 -t 10 -en 0.5 0 0.1 -en 1.2 0 0.8

Random Number Seeds
-------------------

By default, seeds are taken from ``/dev/urandom``, making discoal suitable for cluster computing. You can specify seeds manually:

.. code-block:: bash

   ./discoal 10 5 10000 -t 10 -d 12345 67890

Time Units
----------

discoal measures time in units of 4N₀ generations (coalescent time units), where N₀ is the current effective population size. This is the same convention used by ms.