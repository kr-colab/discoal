Changes in discoal 2.0
======================

discoal 2.0 follows discoal 0.1.7 (June 2025). It writes tree sequences,
reads YAML configuration files and demes models, has new demographic
events, checks every value it reads, fixes several defects in the sweep
and recombination models, and is much faster at high recombination rates.
Some commands that worked with 0.1.x now stop with an error, and the
output for the same seeds changes. The sections below list what changed.

New Features
------------

Tree sequence output
^^^^^^^^^^^^^^^^^^^^

``-ts file.trees`` writes each replicate as a tskit tree sequence
(``file.trees``, or ``file_rep1.trees``, ``file_rep2.trees``, ... with more
than one replicate). As in msprime, it records
only the nodes and edges where lineages coalesce. ``-F`` records the full
ancestral recombination graph instead: every node and edge, also of
recombination. Node times are in units of 2N0 generations.

YAML configuration and demes models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``-Y config.yaml`` reads the parameters of a run from a YAML file. YAML
has a key for each model option of the command line, except continuous
growth (``-eg``, ``-eG``, ``-el``, ``-eL``): give these on the command
line with ``-Y``, or use a demes model. A demes model (the demes
specification of the PopSim consortium) goes in
``demography.demes_filename``; discoal converts its sizes, times,
migration rates, pulses and splits. See :doc:`yaml_configuration`,
:doc:`demes_integration` and the examples in ``config_examples/``.

New demographic events
^^^^^^^^^^^^^^^^^^^^^^

* ``-eg time popnID alpha`` and ``-eG time alpha``: exponential growth.
* ``-el time popnID gamma`` and ``-eL time gamma``: linear growth.
* ``-em time popnID1 popnID2 rate`` and ``-eM time rate``: changes of
  migration rates.

Times are in units of 4N0 generations, and alpha and gamma are rates per
4N0 generations, as in ms. Sweeps with continuous growth of population 0
are not modeled; discoal stops with an error (use ``-en`` size changes).

Other new options
^^^^^^^^^^^^^^^^^

* ``-S file``: a sweep whose allele-frequency trajectory comes from a file
  (lines of time and frequency).
* ``-v``: print progress and the details of a demes conversion to stderr.
  Without it, discoal prints only errors and warnings.

Checks of the input
^^^^^^^^^^^^^^^^^^^

discoal checks the range of every value on the command line, in YAML files
and in priors, and stops with a message that names the option. Models that
cannot finish (lineages that can never meet, a size that reaches zero) also
stop with a message.

Speed and memory
^^^^^^^^^^^^^^^^

discoal 2.0 records tree sequences as msprime does, so time and memory grow
slowly with the recombination rate. For 20 samples, 10,000 sites and
rho = 10,000, one replicate takes 0.8 s and 18 MB, against 18 s and 3 GB
with discoal 0.1.7, and 1.8 s with msprime. Writing ms output with many
segregating sites is up to 17 times faster. With more than 1000 samples,
discoal 0.1.x wrote past the end of an array at the start of each
replicate, and could crash or give wrong output; discoal 2.0 does not.

Breaking Changes
----------------

These changes can make a command that worked with discoal 0.1.x fail, or
give different output.

The random number generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

discoal 2.0 uses xoshiro256++. The same seeds (``-d``) therefore give
different output than in discoal 0.1.x. The distribution of the output is
the same, apart from the fixes below.

``-T`` and ``-U`` are removed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``-T`` wrote Newick trees, and ``-U time`` recorded mutations only back to
a given time. discoal 2.0 stops with an error for both. Write tree
sequences (``-ts``) instead: tskit writes Newick trees, and puts mutations
on the branches you choose.

``-C`` is removed
^^^^^^^^^^^^^^^^^

discoal 0.1.x had ``-C lSpot rSpot``, which kept only the replicates with
recombination in sites ``lSpot`` to ``rSpot`` during a sweep. discoal 2.0
does not have this option, and stops with "Unknown option '-C'".

To get the same result, write tree sequences (``-ts``) and use tskit to
keep only the replicates whose trees change in the region.

``-s`` is not supported
^^^^^^^^^^^^^^^^^^^^^^^

In ms, ``-s k`` gives exactly ``k`` segregating sites in each replicate.
discoal 0.1.x read ``-s`` but did not use it: the number of sites came from
``-t`` only, and discoal gave no warning. discoal 2.0 stops with an error
if you give ``-s``.

To get a fixed number of segregating sites, write tree sequences
(``-ts``) and use tskit to put that number of mutations on the trees.

Gene conversion needs a tract length of at least 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tract lengths are geometric, with the mean that you give to ``-g`` or
``-gr`` (in YAML, ``gene_conversion_tract_length``). discoal 0.1.x accepted
a mean of 0, for example ``-g 10 0``, and then drew tract lengths that
were not numbers. discoal 2.0 stops with an error if the gene conversion
rate is above 0 and the mean tract length is below 1.

With a mean of 1, discoal 0.1.x drew tracts of 0 sites, so gene conversion
had no effect. In discoal 2.0 every tract then has 1 site, so the output
of such a run changes.

Sweeps when population 0 was larger in the past
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A sweep trajectory is accepted with a probability proportional to the
size of population 0 at the allele's origin. discoal 0.1.x divided this
size by the largest size along each trajectory only. When population 0
was larger further in the past, a trajectory with a recent origin never
reached that size, so it was accepted too often. discoal 2.0 divides by
the largest size of population 0 at any time an origin can fall in. The
output of such models changes; when population 0 was not larger in the
past, it does not.

A sweep whose trajectories are almost never accepted (a size drop by a
factor of about 10^6) now stops with an error, and does not run for
ever.

A linear size change must not reach zero in the last interval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``-el`` or ``-eL`` and a positive gamma, the population gets smaller
back in time. If no later event changes its size, it reaches zero, and
lineages that are still in it coalesced at once. discoal 2.0 stops with
an error for such a model, for example ``-el 0 0 0.5`` alone. Add a size
change (``-en``) before the time at which the size reaches zero.

Recombination splits a lineage anywhere in its span
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A lineage carries chunks of material whose sites have not found their
most recent common ancestor (MRCA), with gaps between them. In Hudson's
model, as in ms and msprime, a breakpoint anywhere between the first and
the last of these sites splits the lineage. discoal 0.1.x did not split a
lineage at a breakpoint in a gap whose sites had found their MRCA. Thus
two distant sites stayed linked for too long: the probability that they
have the same MRCA was 2 to 4 times that of msprime. discoal 2.0 splits
the lineage at any breakpoint in its span, and agrees with msprime. The
output of every run with recombination (``-r``) changes.

Gene conversion (``-g``, ``-gr``) now starts a tract at any site in the
lineage's span, also in such a gap. The gene conversion model does not
change otherwise: a tract starts inside the span, and a tract that
starts to the left of the lineage's first site is not modeled.

Gene conversion in the wild-type class during a sweep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When gene conversion splits a lineage during a sweep, the parent that
does not get the sweep site is a random draw from the population: it
joins the lineage's class with that class's frequency. For a lineage in
the wild-type class, discoal 0.1.x used the frequency of the beneficial
allele instead. This gave too little diversity at the flanks of a sweep.
The output of sweeps with gene conversion changes.

One single sweep
^^^^^^^^^^^^^^^^

The sweep sampler follows one sweep at a time, but discoal 0.1.x accepted
several single-sweep options. discoal 2.0 stops with an error if you give
more than one of ``-ws``, ``-wd``, ``-wn``, ``-ls``, ``-ld``, ``-ln`` and
``-S``.
Recurrent sweeps (``-R``, ``-L``) are not affected.

Other fixes that change the output
----------------------------------

* Soft sweeps from standing variation (``-f``): the neutral part of the
  trajectory, below f0, was wrong, and left too much diversity (14% more
  segregating sites than forward simulations in SLiM).
* Recurrent sweeps (``-R``, ``-L``): a run crashed as soon as a sweep
  happened.
* Sweeps with more than one population: the events within population 0
  were chosen with the wrong total rate, and a run could stop with "tried
  to pick allele from popn 0".
* Ancient samples taken during a sweep (``-A``): they now get the
  frequency class of the sweep at their sampling time.

Known Differences from msprime
------------------------------

discoal's gene conversion model differs from msprime's in two points
(GitHub issue #143): a tract starts only inside a lineage's span, so tracts
that start to the left of its first site are not modeled, and a tract that
holds no material gives a parent with no material, where msprime makes no
event. Apart from gene conversion, discoal and msprime agree in
distribution on the stdpopsim catalog and the demes specification's models
(checked every week).
