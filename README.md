[![CI](https://github.com/kr-colab/discoal/actions/workflows/ci.yml/badge.svg)](https://github.com/kr-colab/discoal/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/discoal/badge/?version=latest)](https://discoal.readthedocs.io/en/latest/?badge=latest)

# discoal
discoal is a coalescent simulation program capable of simulating models with recombination, selective sweeps, and demographic changes including population splits and admixture events.

**discoal 2.0** writes tree sequences (`-ts`), reads YAML configuration files and demes models (`-Y`), has growth and migration-change events (`-eg`, `-el`, `-em`), and is checked against msprime and SLiM. It is faster than msprime at high recombination rates. See [the list of changes](https://discoal.readthedocs.io/en/latest/changes.html), which includes changes that can make old commands fail or give different output.

## Documentation

Full documentation is available at [https://discoal.readthedocs.io](https://discoal.readthedocs.io)

The documentation includes:
- Installation instructions
- Detailed usage guides and examples
- Complete parameter reference
- Advanced features and workflows

Read the Docs also offers the documentation as a PDF and as an ePub.

### Additional Documentation

- [YAML Configuration](https://discoal.readthedocs.io/en/latest/yaml_configuration.html) - Using YAML files to specify simulation parameters
- [Demes Integration](https://discoal.readthedocs.io/en/latest/demes_integration.html) - Using demes-format demographic models

## Quick Start

discoal needs a C compiler (gcc or clang) and make; the libraries it uses
are in `extern/`. Clone this repository, change into its directory, and
type `make`. This builds `build/discoal` and links it as `./discoal`
(`make discoal` builds only `build/discoal`).

`./discoal --help` prints the usage:
```
$ ./discoal --help
discoal version 2.0.0
usage: discoal sampleSize numReplicates nSites -ws tau
parameters: 
	 -t theta 
	 -r rho (=zero if not specified)
	 -g conversionRate tractLengthMean (gene conversion)
	 -gr conversionToCrossoverRatio tractLengthMean (gene conversion where initiation rate = rho*conversionToCrossoverRatio)
	 -p npops sampleSize1 sampleSize2 etc.
	 -Y configFile.yaml (load parameters from a YAML configuration file; a demes model goes in demography.demes_filename)
	 -en time popnID size (changes size of popnID)
	 -eg time popnID alpha (sets popnID to exponential growth at rate alpha per 4N0 generations, as in ms; positive alpha = past was smaller)
	 -eG time alpha (-eg applied to all populations)
	 -el time popnID gamma (sets popnID to linear growth; gamma is the change in relative size per 4N0 generations; positive gamma = past was smaller)
	 -eL time gamma (-el applied to all populations)
	 -ed time popnID1 popnID2 (joins popnID1 into popnID2; -ej is an alias)
	 -ea time daughterPopnID founderPopnID1 founderPopnID2 admixProp (admixture-- back in time daughterPopnID into two founders)
	 -em time popnID1 popnID2 migRate (changes migration rate from popnID1 to popnID2 starting at time)
	 -eM time migRate (sets all off-diagonal migration rates to migRate starting at time)
	 -ws tau (sweep ended tau ago, in units of 4N0 generations; stochastic sweep)
	 -wd tau (sweep ended tau ago, in units of 4N0 generations; deterministic sweep)
	 -wn tau (sweep ended tau ago, in units of 4N0 generations; neutral sweep)
	 -ls tau leftRho (stochastic sweep some genetic distance to the left of the simulated window--at 4Nr = 2 x leftRho)
		 similarly, ld and ln simulate deterministic and neutral sweeps to the left of the window, respectively
	 -f first frequency at which selection acts on allele (F0; sweep models only)
	 -uA rate at which adaptive mutation recurs during the sweep phase (sweep models only)
	 -N sweepEffectivePopnSize (sweep models only)
	 -a alpha (=2Ns)
	 -x sweepSite (0-1)
	 -c partialSweepFinalFrequency (partial sweeps)
	 -i dt (sweep time increment scalar; default 40 -> 1/40N)
	 -S file (sweep trajectory from a file: lines of time (4N0 generations, increasing) and allele frequency, ending at frequency 0)
	 -M migRate (sets all rates to migRate)
	 -m popnID1 popnID2 migRate (sets migRate from popnID1 to popnID2)
	 -A sampleSize popnID time (ancient sample from popnID at specified time)
	 -Pt low high (prior on theta)
	 -Pr low high (prior on rho)
	 -Pre mean upperBound (prior on rho -- exponentially distributed but truncated at an upper bound)
	 -Pa low high (prior on alpha)
	 -Pu low high (prior on tau; sweep models only; still must use "-ws tau" (or a YAML sweep, not both) and "tau" will be ignored)
	 -PuA low high (prior on uA; sweep models only)
	 -Px low high (prior on sweepSite; sweep models only)
	 -Pf low high (prior on F0; sweep models only)
	 -Pc low high (prior on partialSweepFinalFreq; sweep models only)
	 -Pe1 lowTime highTime lowSize highSize (priors on first demographic move time and size)
	 -Pe2 lowTime highTime lowSize highSize (priors on second demographic move time and size)
	 -R rhhRate (recurrent hitch hiking mode at the locus; rhh is rate per 2N individuals / generation)
	 -L rhhRate (recurrent hitch hiking mode to the side of locus; leftRho is ~Unif(0,4Ns); rhh is rate per 2N individuals / generation)
	 -h (hide selected SNP in partial sweep mode)
	 -v (verbose: print progress and demes conversion details to stderr)
	 -ts filename (tree sequence output mode - saves to tskit file)
	 -F (full ARG mode - record all edges including recombination; use with -ts)
	 -d seed1 seed2 (set random number generator seeds)
	 --help (print this text)
	 --version (print the version)
```

