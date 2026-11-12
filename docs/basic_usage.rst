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

All scaled parameters follow ms. N₀ is the reference size and g is one generation.

==========================================  ==========================  =====================================
Parameter                                   discoal unit                msprime equivalent
==========================================  ==========================  =====================================
``-t`` theta                                4N₀μ for the whole locus    ``mutation_rate`` = theta / (4N₀ L)
``-r`` rho                                  4N₀r for the whole locus    ``recombination_rate`` = rho / (4N₀ L)
``-g`` gamma                                4N₀g for the whole locus    ``gene_conversion_rate`` = gamma / (4N₀ L)
``-g`` mean tract length                    sites                       ``gene_conversion_tract_length``
times (``-en``, ``-ed``, ``-em``, ``-ws``)  4N₀ generations             generations = time × 4N₀
``-en`` size                                size / N₀                   ``initial_size`` = size × N₀
``-eg`` alpha                               per 4N₀ generations         ``growth_rate`` = alpha / (4N₀)
``-el`` gamma                               size / N₀ per 4N₀ g         size change per g = gamma / 4
``-M``, ``-m``, ``-em`` migRate             4N₀m                        migration rate = migRate / (4N₀)
``-a`` alpha                                2N₀s                        s = alpha / (2N₀)
==========================================  ==========================  =====================================

``-m i j`` moves lineages from population i to population j backward in time, which is the ms and msprime convention. In YAML configs, times and rates use the same units as the command line. Demes files use generations (or years) and per-generation rates; discoal converts them. Tree sequence output stores node times in units of 2N₀ generations, discoal's internal scale.

``testing/units_convention_check.py`` checks each of these conventions against msprime; ``testing/gene_conversion_check.py`` checks gene conversion.