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

.. option:: -eg time popID alpha

   Exponential growth, as in ms. From ``time`` back into the past, the size
   is N(t) = N(time) exp(-alpha (t - time)).

   * time: When growth starts, looking backward (4N₀ generations)
   * popID: Population ID (0-indexed)
   * alpha: Growth rate per 4N₀ generations. Positive alpha means the
     population grew forward in time, so the past was smaller. The
     msprime ``growth_rate`` is alpha / (4N₀).

   Events at the same time are applied in command-line order, as in ms.
   ``-en`` stops growth, so ``-eg t ... -en t ...`` ends with no growth,
   and ``-en t ... -eg t ...`` grows from the new size.

.. option:: -eG time alpha

   Same as ``-eg`` for every population.

.. option:: -el time popID gamma

   Linear growth. From ``time`` back into the past, the size is
   N(t) = N(time) - gamma (t - time), relative to N₀.

   * time: When growth starts, looking backward (4N₀ generations)
   * popID: Population ID (0-indexed)
   * gamma: Change in relative size (size / N₀) per 4N₀ generations.
     Positive gamma means the past was smaller. discoal stops with an
     error if the size reaches zero before the next event.

.. option:: -eL time gamma

   Same as ``-el`` for every population.

.. option:: -ed time pop1 pop2

   Population split (forward in time):
   
   * time: When populations merge (backward)
   * pop1: Source population
   * pop2: Destination population

   ``-ej`` is an alias. The merge sets migration between pop1 and pop2 to
   zero. It does not change migration between pop1 and other populations;
   use ``-em`` for that.

.. option:: -ea time admixPop pop1 pop2 prop

   Admixture event:
   
   * time: When admixture occurs
   * admixPop: Admixed population
   * pop1, pop2: Source populations
   * prop: Proportion from pop1: backward in time, each lineage of
     admixPop moves to pop1 with probability prop and to pop2 otherwise

.. option:: -M migRate

   Set all migration rates to same value (4N₀m)

.. option:: -m pop1 pop2 migRate

   Set specific migration rate from pop1 to pop2 (4N₀m). Backward in time,
   lineages move from pop1 to pop2 at this rate, as in ms and msprime.

.. option:: -em time pop1 pop2 migRate

   From ``time`` back into the past, set the migration rate from pop1 to
   pop2 (4N₀m, same direction as ``-m``).

.. option:: -eM time migRate

   From ``time`` back into the past, set all migration rates (4N₀m).

.. option:: -A numSamples popID time

   Ancient samples:
   
   * numSamples: Number of ancient lineages
   * popID: Population to sample from
   * time: Sampling time in past, in units of 4N0 generations

   The ancient samples are part of the population's sample count: with
   ``-p 2 10 10 -A 4 1 0.3``, population 1 has 6 samples taken at time 0
   and 4 taken at time 0.3. The ancient samples of a population must not
   outnumber its samples, and the time must be ``>= 0``. Give ``-A`` more
   than once for samples at different times.

   In the output, the samples of each population come in population
   order. Within a population, the ancient samples come first, the
   earliest sampling time first, then the samples taken at time 0. In the
   example above, the rows of population 1 are 4 ancient samples, then 6
   present-day samples.

   An ancient sample taken while a sweep is in progress (in population 0)
   carries the beneficial allele with the allele's frequency at the
   sampling time; with a partial sweep (``-c``), its carriers carry the
   selected SNP.

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

   Final frequency for partial sweep. It must be above 1/(2N), one copy
   of the allele (N from ``-N``); with ``-Pc``, the high bound must be. A
   frequency at or below one copy when the sweep starts (a draw from
   ``-Pc``, or after a size change) gives a sweep that ends at once.

.. option:: -S file

   Sweep with its trajectory read from a file instead of simulated. Each
   line holds a time before the present (in units of 4N₀ generations) and
   the frequency of the selected allele, separated by white space; blank
   lines and lines that start with ``#`` are skipped. Times must increase
   down the file, every frequency but the last must be in (0, 1], and the
   last must be 0, the allele's origin. The sweep starts (going back in
   time) at the first time. discoal interpolates the frequency linearly at
   each of its sweep time steps (``-N``, ``-i``), and the sweep ends at the
   first step at which the frequency is at or below 1/(2N), as a simulated
   sweep does. Frequencies above 1 - 1/(2N) are used as 1 - 1/(2N). If the
   first frequency is below 1, the sweep is partial, as with ``-c``.
   ``-S`` cannot be combined with ``-ws``, ``-wd``, ``-wn``, ``-ls``,
   ``-ld``, ``-ln``, ``-R``, ``-L``, ``-a``, ``-f``, ``-c``, ``-Pa``,
   ``-Pu``, ``-Pf`` or ``-Pc``.

.. option:: -ls tau leftRho

   Stochastic sweep to the left of the locus. The selected site is at
   genetic distance 4Nr = 2 x leftRho from the locus; that is, leftRho
   is 2Nr. (This factor of 2 is a long-standing convention of
   discoal; ``-L`` draws the distance 4Nr directly.)

.. option:: -ld tau leftRho

   Deterministic sweep to the left of the locus

.. option:: -ln tau leftRho

   Neutral fixation to the left of the locus

.. option:: -R rate

   Recurrent stochastic sweeps at random positions inside the locus, at
   this rate per 2N individuals per generation

.. option:: -L rate

   Recurrent stochastic sweeps to the left of the locus; each sweep's
   distance (4Nr) is drawn uniformly from (0, 4Ns)

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

   Uniform prior on f₀, with ``0 <= low <= high < 1``.

.. option:: -Pc low high

   Uniform prior on partial sweep frequency, with ``0 <= low <= high < 1``
   and ``high`` above 1/(2N). With a soft sweep (``-f`` or ``-Pf``),
   every f₀ must be below every final frequency.

.. option:: -Pe1 timeLow timeHigh sizeLow sizeHigh

   A size change of population 0 with uniform priors on its time (units
   of 4N generations) and its size (relative to N), drawn for each
   replicate. The bounds must satisfy 0 <= timeLow <= timeHigh and
   0 < sizeLow <= sizeHigh. The size trajectories are checked before the
   run at the lower bounds, and again after each replicate's draws;
   discoal stops if a draw makes a population size reach zero.

.. option:: -Pe2 timeLow timeHigh sizeLow sizeHigh

   A second size change with priors, as ``-Pe1``

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

.. option:: -U time

   Only record mutations more recent than time

.. option:: -h

   Hide selected SNP from output

.. option:: -v

   Verbose. Print progress and demes conversion details to stderr. Without
   ``-v``, discoal writes to stderr only for errors and warnings. Works with
   ``-Y``.

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