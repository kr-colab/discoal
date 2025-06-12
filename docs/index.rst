.. discoal documentation master file

discoal - a coalescent simulator with selection
================================================

**discoal** is a coalescent simulator for generating genetic data under models with recombination, population structure, and selection. It is designed to be familiar to users of Richard Hudson's ``ms`` program while adding the ability to simulate selective sweeps and other forms of natural selection.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   installation
   basic_usage
   recombination
   population_structure
   selection
   advanced_features
   examples
   development
   api_reference

Features
--------

* Coalescent simulation with recombination and gene conversion
* Multiple population models with migration and admixture
* Various selection models including hard and soft sweeps
* Partial sweeps and recurrent hitchhiking
* Efficient memory usage with dynamic allocation
* Output compatible with ``ms`` format

Quick Start
-----------

.. code-block:: bash

   # Download and compile
   git clone https://github.com/kern-lab/discoal.git
   cd discoal
   make discoal

   # Basic simulation: 10 samples, 5 replicates, 10000 sites, theta=10
   ./discoal 10 5 10000 -t 10

   # With recombination (rho=10)
   ./discoal 10 5 10000 -t 10 -r 10

   # With a selective sweep
   ./discoal 10 5 10000 -t 10 -r 10 -ws 0.01 -a 1000 -x 0.5

Citing discoal
--------------

If you use discoal in your research, please cite:

Kern, A.D. and Schrider, D.R. (2016). discoal: flexible coalescent simulations with selection. *Bioinformatics*, 32(24), 3839-3841.

Support
-------

* Report issues on `GitHub <https://github.com/kern-lab/discoal/issues>`_
* Contact: adkern@uoregon.edu

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`