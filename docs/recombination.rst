Recombination and Gene Conversion
=================================

Crossing Over (Recombination)
-----------------------------

Recombination is modeled using the ``-r`` flag, which takes the population recombination rate ρ = 4Nr as a parameter, where:

* N is the current effective population size
* r is the probability of crossover per base pair per generation

Recombination can occur between any of the discrete sites being modeled.

**Example:**

.. code-block:: bash

   # Simulate with recombination rate rho = 2.4
   ./discoal 3 2 100 -t 2 -r 2.4

Gene Conversion
---------------

Gene conversion (non-crossover recombination events) is simulated using the ``-g`` flag. This models recombination without exchange of flanking markers.

**Parameters:**

1. γ = 4N₀gL: The population gene conversion rate for the whole locus of L sites, scaled as ρ is

   * g is the probability per generation that a gene conversion tract starts at a given site

2. Mean tract length: The mean number of sites in a tract

   * Tract lengths are geometric on 1, 2, 3, ..., with this mean

**The model** is that of ms and msprime (Wiuf and Hein 2000). A tract starts at a
site and runs to the right; the sites it covers come from one parent, and the other
sites from the other parent. The chromosome goes on beyond the locus, so a tract
can also start to the left of a lineage's first site and run into its material,
and one can start in the locus and run past its end. A tract that changes no
lineage's material (it falls in a gap of the material, or covers all of it) has
no effect. ``testing/gene_conversion_check.py`` checks this model against the
two-locus theory and against msprime.

**Example:**

.. code-block:: bash

   # Gene conversion with rate 2.4 and mean tract length 10bp
   ./discoal 3 2 100 -t 2 -r 2.4 -g 2.4 10

Alternative Gene Conversion Model
---------------------------------

The ``-gr`` flag allows specifying gene conversion as a ratio to crossover events:

.. code-block:: bash

   # Gene conversion rate = 0.5 * recombination rate, tract length 50bp
   ./discoal 10 5 10000 -t 10 -r 20 -gr 0.5 50

Here the gene conversion rate is γ = ρ × conversionToCrossoverRatio, so the ratio is g/r, as for ms's ``-c f``.

Combining Recombination Models
------------------------------

Both crossing over and gene conversion can be active simultaneously:

.. code-block:: bash

   # Both crossovers and gene conversions
   ./discoal 20 10 50000 -t 50 -r 40 -g 20 100

This simulates:
* Crossover rate: ρ = 40
* Gene conversion rate: γ = 20  
* Mean conversion tract length: 100bp

Implementation Notes
--------------------

* Recombination breakpoints are tracked dynamically
* The simulator uses efficient segment-based tracking of ancestral material
* Recombination events only occur in regions with ancestral material
* During selective sweeps, recombination continues but migration is suspended