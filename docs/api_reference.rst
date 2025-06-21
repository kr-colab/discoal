API Reference
=============

Command Line Options
--------------------

Basic Parameters
^^^^^^^^^^^^^^^^

.. option:: sampleSize

   Number of chromosomes to sample (required)

.. option:: numReplicates

   Number of independent simulations to run (required)

.. option:: nSites

   Number of discrete sites in the simulated region (required)

.. option:: -t theta

   Population mutation rate (4N₀μ)

.. option:: -r rho

   Population recombination rate (4N₀r)

.. option:: -g gamma tractLength

   Gene conversion parameters:
   
   * gamma: Population gene conversion rate (4N₀g)
   * tractLength: Mean conversion tract length in bp

.. option:: -gr ratio tractLength

   Gene conversion as ratio to crossover rate

Demographic Options
^^^^^^^^^^^^^^^^^^^

.. option:: -p numPops size1 size2 ...

   Multiple populations with specified sample sizes

.. option:: -en time popID size

   Population size change:
   
   * time: When change occurs (4N₀ generations)
   * popID: Population ID (0-indexed)
   * size: New size relative to N₀

.. option:: -ed time pop1 pop2

   Population split (forward in time):
   
   * time: When populations merge (backward)
   * pop1: Source population
   * pop2: Destination population

.. option:: -ea time admixPop pop1 pop2 prop

   Admixture event:
   
   * time: When admixture occurs
   * admixPop: Admixed population
   * pop1, pop2: Source populations
   * prop: Proportion from pop1

.. option:: -M migRate

   Set all migration rates to same value (4N₀m)

.. option:: -m pop1 pop2 migRate

   Set specific migration rate from pop1 to pop2

.. option:: -A numSamples popID time

   Ancient samples:
   
   * numSamples: Number of ancient lineages
   * popID: Population to sample from
   * time: Sampling time in past

Selection Options
^^^^^^^^^^^^^^^^^

.. option:: -ws tau

   Stochastic selective sweep:
   
   * tau: Time since fixation (4N₀ generations)

.. option:: -wd tau

   Deterministic selective sweep

.. option:: -wn tau

   Neutral fixation

.. option:: -a alpha

   Selection coefficient (2Ns)

.. option:: -x position

   Position of selected site (0-1)

.. option:: -f frequency

   Initial frequency for soft sweep

.. option:: -uA rate

   Recurrent adaptive mutation rate

.. option:: -c frequency

   Final frequency for partial sweep

.. option:: -ls tau

   Stochastic sweep to the left of locus

.. option:: -ld tau  

   Deterministic sweep to the left of locus

.. option:: -ln tau

   Neutral fixation to the left of locus

.. option:: -R rate

   Recurrent hitchhiking rate at locus

.. option:: -L rate

   Recurrent hitchhiking rate to left of locus

Prior Distributions
^^^^^^^^^^^^^^^^^^^

.. option:: -Pt low high

   Uniform prior on θ

.. option:: -Pr low high

   Uniform prior on ρ

.. option:: -Pre mean upper

   Exponential prior on ρ (truncated)

.. option:: -Pa low high

   Uniform prior on α

.. option:: -Pu low high

   Uniform prior on τ

.. option:: -PuA low high

   Uniform prior on uA

.. option:: -Px low high

   Uniform prior on sweep position

.. option:: -Pf low high

   Uniform prior on f₀

.. option:: -Pc low high

   Uniform prior on partial sweep frequency

.. option:: -Pe1 timeLow timeHigh sizeLow sizeHigh

   Prior on first demographic event

.. option:: -Pe2 timeLow timeHigh sizeLow sizeHigh

   Prior on second demographic event

Advanced Options
^^^^^^^^^^^^^^^^

.. option:: -d seed1 seed2

   Set random number seeds

.. option:: -N size

   Effective population size during sweeps (default: 1000000)

.. option:: -i scalar

   Time increment scalar for sweeps (default: 40)

.. option:: -T

   Output genealogical trees in Newick format

.. option:: -C leftBound rightBound

   Condition on recombination in specified region

.. option:: -U time

   Only record mutations more recent than time

.. option:: -h

   Hide selected SNP from output

Output Format
-------------

Standard Output
^^^^^^^^^^^^^^^

ms-compatible format::

   command_line
   seed1 seed2
   
   //
   segsites: n
   positions: pos1 pos2 ... posn
   haplotype1
   haplotype2
   ...

Tree Output
^^^^^^^^^^^

With ``-T`` flag::

   //
   [nsites1]tree1;
   [nsites2]tree2;
   ...

Where trees are in Newick format and nsites indicates how many sites have that genealogy.

Exit Codes
----------

* 0: Success
* 1: General error
* 666: Parameter error or invalid configuration