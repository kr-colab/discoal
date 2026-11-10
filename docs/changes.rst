Changes in discoal 2.0
======================

Breaking Changes
----------------

These changes can make a command that worked with discoal 1.x fail, or
give different output.

``-C`` is removed
^^^^^^^^^^^^^^^^^

discoal 1.x had ``-C lSpot rSpot``, which kept only the replicates with
recombination in sites ``lSpot`` to ``rSpot`` during a sweep. discoal 2.0
does not have this option, and stops with "Unknown option '-C'".

To get the same result, write tree sequences (``-ts``) and use tskit to
keep only the replicates whose trees change in the region.

``-s`` is not supported
^^^^^^^^^^^^^^^^^^^^^^^

In ms, ``-s k`` gives exactly ``k`` segregating sites in each replicate.
discoal 1.x read ``-s`` but did not use it: the number of sites came from
``-t`` only, and discoal gave no warning. discoal 2.0 stops with an error
if you give ``-s``.

To get a fixed number of segregating sites, write tree sequences
(``-ts``) and use tskit to put that number of mutations on the trees.

Gene conversion needs a tract length of at least 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tract lengths are geometric, with the mean that you give to ``-g`` or
``-gr`` (in YAML, ``gene_conversion_tract_length``). discoal 1.x accepted
a mean of 0, for example ``-g 10 0``, and then drew tract lengths that
were not numbers. discoal 2.0 stops with an error if the gene conversion
rate is above 0 and the mean tract length is below 1.

With a mean of 1, discoal 1.x drew tracts of 0 sites, so gene conversion
had no effect. In discoal 2.0 every tract then has 1 site, so the output
of such a run changes.

Sweeps when population 0 was larger in the past
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A sweep trajectory is accepted with a probability proportional to the
size of population 0 at the allele's origin. discoal 1.x divided this
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
