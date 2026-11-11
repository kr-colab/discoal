Introduction
============

**discoal** is a program aimed at generating samples under a coalescent model with recombination, step-wise changes in population size, and selection. For users familiar with Richard Hudson's ``ms``, the usage of discoal will be quite familiar, and indeed was meant to "play nice" with programs the user, or others, may have written for analyzing ``ms`` style output. 

discoal is not meant to take the place of ``ms``—indeed in comparison it is quite limited in what it can do—but instead is meant to add a few models not covered by other simulation programs. In particular, discoal can quickly generate samples from models with selective sweeps in the history of the sample along with stepwise population size changes.

Why discoal?
------------

discoal gets its name from the contraction of "discrete" and "coalescent", because it handles recombination along the chromosome, and its associated programmatic bookkeeping, by modeling a discrete number of ancestral sites. 

The current optimized version uses dynamic memory allocation and segment-based tracking for efficient memory usage, allowing simulations of large genomic regions (up to 100 million sites) with reasonable memory requirements. Memory usage now scales efficiently with the actual complexity of the simulation rather than pre-allocating for maximum possible sizes.

Key Features
------------

* **Selective sweeps**: Both hard and soft sweeps, partial sweeps, and recurrent hitchhiking
* **Population structure**: Multiple populations with migration and admixture
* **Demographic events**: Step-wise population size changes and population splits
* **Recombination**: Both crossing over and gene conversion
* **Ancient samples**: Ability to sample lineages at different time points
* **Efficient memory usage**: Dynamic allocation and segment-based ancestry tracking
* **ms compatibility**: Output format compatible with existing ``ms`` analysis pipelines

Technical Details
-----------------

discoal implements several optimizations for efficient simulation:

* **Dynamic memory allocation**: Arrays grow as needed rather than pre-allocating maximum sizes
* **Segment-based ancestry tracking**: Each lineage keeps a list of the segments of its ancestral material, not an array of sites
* **Tree sequence recording as in msprime**: Only the nodes and edges where lineages coalesce are recorded, and segments whose sites have found their common ancestor are dropped
* **Weighted recombination draws**: Each population keeps a Fenwick tree of its lineages' spans, so a recombination picks its lineage and breakpoint in O(log n)
* **Trajectory files**: Sweep trajectories are kept in temporary files, memory-mapped, and deleted when closed

Memory use stays at tens of megabytes even at high recombination rates.