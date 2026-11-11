Installation
============

System Requirements
-------------------

discoal should build on most Unix-like systems including:

* Linux (various distributions)
* macOS 
* Windows with WSL (Windows Subsystem for Linux)

Prerequisites:

* C compiler (gcc or clang)
* make
* Standard C libraries including math library

The libraries that discoal uses (tskit, libyaml, libcyaml, demes-c) are in
the repository, in ``extern/``; you do not install them.

Download
--------

The source code for discoal is available from GitHub:

.. code-block:: bash

   git clone https://github.com/kr-colab/discoal.git

Or download the latest release as a zip file from the `releases page <https://github.com/kr-colab/discoal/releases>`_.

Compilation
-----------

Basic compilation is straightforward:

.. code-block:: bash

   cd discoal
   make

This builds ``build/discoal`` and links it as ``./discoal``. ``make discoal``
builds only ``build/discoal``. The build uses ``-march=native``, so the
program runs only on processors like the one that built it.

Additional build targets:

.. code-block:: bash

   # Build with debugging symbols
   make discoal_debug
   
   # Clean build artifacts
   make clean
   
   # Build and run the unit tests
   make run_tests

   # Build, run the unit tests and the checks that CI runs
   # (needs Python with msprime, tskit, numpy, scipy, demes and stdpopsim)
   make check

Compilation Options
-------------------

discoal supports up to 65,535 samples.

To change the maximum number of sites (default 100 million), edit ``MAXSITES`` in ``src/core/discoal.h`` and recompile.

Testing the Installation
------------------------

Run a simple test to verify the installation:

.. code-block:: bash

   ./discoal 10 1 1000 -t 10

This should produce output in ms format showing segregating sites and haplotypes.

Troubleshooting
---------------

**Random number seeding**: discoal uses ``/dev/urandom`` for random number generation. This makes it suitable for cluster computing where multiple jobs launch simultaneously without worrying about seed collisions.

**Memory issues**: If you encounter memory allocation errors, ensure your system has sufficient RAM.

**Compilation errors**: Ensure you have a C99-compatible compiler. On older systems, you may need to add ``-std=c99`` to the compiler flags.