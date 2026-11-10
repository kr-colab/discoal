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
