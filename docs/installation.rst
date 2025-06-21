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

Download
--------

The source code for discoal is available from GitHub:

.. code-block:: bash

   git clone https://github.com/kern-lab/discoal.git

Or download the latest release as a zip file from the `releases page <https://github.com/kern-lab/discoal/releases>`_.

Compilation
-----------

Basic compilation is straightforward:

.. code-block:: bash

   cd discoal
   make discoal

This will create the ``discoal`` executable in the current directory.

Additional build targets:

.. code-block:: bash

   # Build with debugging symbols
   make discoal_debug
   
   # Clean build artifacts
   make clean
   
   # Build and run tests
   make run_tests

Compilation Options
-------------------

The standard compilation supports up to 65,535 samples. The previous limitation of 254 samples (and the need for the ``-DBIG`` flag) has been removed due to memory optimizations.

To change the maximum number of sites (default 100 million), edit ``MAXSITES`` in ``discoal.h`` and recompile.

Testing the Installation
------------------------

Run a simple test to verify the installation:

.. code-block:: bash

   ./discoal 10 1 1000 -t 10

This should produce output in ms format showing segregating sites and haplotypes.

Troubleshooting
---------------

**Random number seeding**: discoal uses ``/dev/urandom`` for random number generation. This makes it suitable for cluster computing where multiple jobs launch simultaneously without worrying about seed collisions.

**Memory issues**: If you encounter memory allocation errors, ensure your system has sufficient RAM. The optimized version significantly reduces memory usage compared to earlier versions.

**Compilation errors**: Ensure you have a C99-compatible compiler. On older systems, you may need to add ``-std=c99`` to the compiler flags.