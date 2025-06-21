Examples
========

This section provides practical examples for common use cases.

Basic Neutral Simulations
-------------------------

**Simple neutral model:**

.. code-block:: bash

   # 20 samples, 100 replicates, 50kb region, θ=40
   ./discoal 20 100 50000 -t 40

**With recombination:**

.. code-block:: bash

   # Add recombination rate ρ=40
   ./discoal 20 100 50000 -t 40 -r 40

**With gene conversion:**

.. code-block:: bash

   # Crossovers + gene conversion (rate 20, tract length 500bp)
   ./discoal 20 100 50000 -t 40 -r 40 -g 20 500

Demographic Models
------------------

**Population expansion:**

.. code-block:: bash

   # 10-fold expansion 0.1 time units ago
   ./discoal 20 100 50000 -t 40 -r 40 -en 0.1 0 0.1

**Bottleneck:**

.. code-block:: bash

   # Severe bottleneck (1% size) followed by recovery
   ./discoal 20 100 50000 -t 40 -r 40 \
     -en 0.05 0 0.01 \
     -en 0.1 0 1.0

**Exponential growth:**

.. code-block:: bash

   # Approximate exponential growth with multiple size changes
   ./discoal 20 100 50000 -t 40 -r 40 \
     -en 0.001 0 0.9 \
     -en 0.002 0 0.8 \
     -en 0.004 0 0.65 \
     -en 0.008 0 0.4 \
     -en 0.016 0 0.16

Selection Examples
------------------

**Recent hard sweep:**

.. code-block:: bash

   # Strong sweep (2Ns=1000) completed 100 generations ago
   # Assuming N=10,000, this is tau = 100/(4*10000) = 0.0025
   ./discoal 20 100 50000 -t 40 -r 40 -ws 0.0025 -a 1000 -x 0.5

**Soft sweep from standing variation:**

.. code-block:: bash

   # Selection on variant at 1% frequency
   ./discoal 20 100 50000 -t 40 -r 40 -ws 0.01 -a 500 -x 0.5 -f 0.01

**Incomplete sweep:**

.. code-block:: bash

   # Sweep to 80% frequency
   ./discoal 20 100 50000 -t 40 -r 40 -ws 0.01 -a 1000 -x 0.5 -c 0.8

**Recurrent hitchhiking:**

.. code-block:: bash

   # Sweeps occur at rate 0.0001 per generation
   ./discoal 20 100 50000 -t 40 -r 40 -R 0.0001 -a 500

Population Structure
--------------------

**Two populations with migration:**

.. code-block:: bash

   # 10 samples each, symmetric migration 4Nm=1
   ./discoal 20 100 50000 -t 40 -r 40 -p 2 10 10 -M 1.0

**Three population phylogeny:**

.. code-block:: bash

   # ((pop0,pop1),pop2) with realistic parameters
   ./discoal 30 100 50000 -t 60 -r 60 -p 3 10 10 10 \
     -ed 0.5 0 1 \
     -ed 2.0 1 2 \
     -m 0 1 0.1 -m 1 0 0.1

**Admixture model:**

.. code-block:: bash

   # Pop0 is 20% pop2, 80% pop1 ancestry (0.01 time units ago)
   ./discoal 30 100 50000 -t 60 -r 60 -p 3 10 10 10 \
     -ea 0.01 0 1 2 0.8

Complex Scenarios
-----------------

**Human-like demographic model:**

.. code-block:: bash

   # Out-of-Africa model approximation
   # Pop0=Africa, Pop1=Europe, Pop2=Asia
   ./discoal 60 100 50000 -t 100 -r 80 -p 3 20 20 20 \
     -en 0.0 1 0.2 \
     -en 0.0 2 0.3 \
     -ed 0.015 1 2 \
     -ed 0.02 2 0 \
     -en 0.02 0 0.25 \
     -m 0 1 0.5 -m 1 0 0.5 \
     -m 0 2 0.5 -m 2 0 0.5 \
     -m 1 2 0.5 -m 2 1 0.5

**Selection with demography:**

.. code-block:: bash

   # Sweep during population bottleneck
   ./discoal 50 100 100000 -t 100 -r 100 \
     -en 0.01 0 0.05 \
     -en 0.05 0 1.0 \
     -ws 0.02 -a 500 -x 0.5

**Ancient DNA with selection:**

.. code-block:: bash

   # 20 modern + 10 ancient samples
   # Sweep occurred between ancient and modern sampling
   ./discoal 30 100 50000 -t 80 -r 80 \
     -A 10 0 0.05 \
     -ws 0.02 -a 1000 -x 0.5

ABC/Machine Learning Applications
----------------------------------

**Parameter estimation setup:**

.. code-block:: bash

   # Generate training data with parameter priors
   for i in {1..10000}; do
     ./discoal 50 1 100000 \
       -Pt 10 100 \
       -Pr 10 100 \
       -Pa 100 2000 \
       -Pu 0.001 0.1 \
       -ws 0.01 -x 0.5 >> training_data.txt
   done

**Feature-rich simulations:**

.. code-block:: bash

   # Complex model with multiple features for ML
   ./discoal 100 1 100000 \
     -Pt 20 80 \
     -Pr 20 80 \
     -p 2 50 50 \
     -M 1.0 \
     -Pe1 0.01 0.5 0.1 10 \
     -Pa 100 2000 \
     -Px 0.1 0.9 \
     -Pf 0.001 0.1 \
     -ws 0.01