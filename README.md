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
A coalescent simulator with recombination, population structure and
selective sweeps.

usage: discoal sampleSize numReplicates nSites [options]
       discoal -Y config.yaml [options]

examples:
  discoal 10 5 10000 -t 10 -r 10
      neutral, with recombination
  discoal 10 5 10000 -t 10 -r 10 -ws 0.01 -a 1000 -x 0.5
      a hard sweep at the middle of the locus
  discoal 20 1 10000 -t 10 -p 2 10 10 -ed 1.0 1 0
      two populations that join
  discoal 10 3 10000 -t 10 -r 10 -ts out.trees
      tree sequences (tskit)

Each option has a short name, as in ms, and a long name; they do the same.

Times are in units of 4N0 generations, sizes in units of N0, and theta,
rho and migration rates are 4N0 times the rate (as in ms). Populations are
numbered from 0.

mutation and recombination:
  -t, --theta theta       4N0 u for the whole locus
  -r, --rho rho           4N0 r for the whole locus (default 0)
  -g, --gene-conversion gamma tractLength
                          gene conversion: 4N0 g for the whole locus, and the
                          mean tract length in sites
  -gr, --gene-conversion-ratio ratio tractLength
                          gene conversion at ratio x rho

populations and samples:
  -p, --populations npops n1 n2 ...
                          number of populations and the sample size of each
  -M, --migration migRate
                          migration between all populations at migRate
                          (give -p first)
  -m, --migration-rate pop1 pop2 migRate
                          migration from pop1 to pop2 (give -p first)
  -A, --ancient-sample n pop time
                          n of the samples of pop are ancient, taken at time

demographic events:
  -en, --size-change time pop size
                          size of pop changes to size
  -eg, --exponential-growth time pop alpha
                          exponential growth of pop at rate alpha per 4N0
                          generations (positive: smaller in the past)
  -eG, --exponential-growth-all time alpha
                          -eg for all populations
  -el, --linear-growth time pop gamma
                          linear growth of pop; gamma is the change in size
                          per 4N0 generations (positive: smaller in the past)
  -eL, --linear-growth-all time gamma
                          -el for all populations
  -ed, --join time pop1 pop2
                          pop1 joins pop2, backward in time (-ej is the same)
  -ea, --admixture time pop src1 src2 prop
                          admixture: backward in time a lineage of pop goes
                          to src1 with probability prop, else to src2
  -em, --migration-change time pop1 pop2 migRate
                          migration from pop1 to pop2 changes to migRate
  -eM, --migration-change-all time migRate
                          all migration rates change to migRate

sweeps:
  -ws, --sweep tau        stochastic sweep that ended tau ago
  -wd, --sweep-deterministic tau
                          deterministic sweep that ended tau ago
  -wn, --sweep-neutral tau
                          neutral fixation that ended tau ago
  -ls, --linked-sweep tau leftRho
                          stochastic sweep to the left of the locus, at a
                          distance of 4N0 r = 2 x leftRho
  -ld, --linked-sweep-deterministic tau leftRho
                          as -ls, deterministic
  -ln, --linked-sweep-neutral tau leftRho
                          as -ls, neutral fixation
  -S, --trajectory-file file
                          sweep whose trajectory is in file: lines of time
                          (increasing) and allele frequency, ending at 0
  -R, --recurrent-sweeps rate
                          recurrent sweeps in the locus; rate per 2N
                          individuals per generation
  -L, --recurrent-linked-sweeps rate
                          recurrent sweeps to the left of the locus, at a
                          distance drawn uniformly from (0, 4Ns)
  -a, --alpha alpha       strength of selection, 2Ns
  -x, --sweep-position position
                          position of the selected site in the locus, 0 to 1
  -f, --initial-frequency f0
                          soft sweep: selection starts at frequency f0
  -uA, --beneficial-mutation-rate rate
                          soft sweep: rate of recurrent beneficial mutation
  -c, --final-frequency frequency
                          partial sweep: final frequency of the allele
  -h, --hide-selected-snp
                          hide the selected SNP of a partial sweep
  -N, --sweep-population-size size
                          effective population size during a sweep
                          (default 1000000)
  -i, --sweep-time-step dt
                          time step of sweep trajectories, 1/(dt N)
                          (default 40)

priors (a value is drawn for each replicate):
  -Pt, --prior-theta low high
                          theta, uniform
  -Pr, --prior-rho low high
                          rho, uniform
  -Pre, --prior-rho-exponential mean upperBound
                          rho, exponential, cut at upperBound
  -Pa, --prior-alpha low high
                          alpha, uniform
  -Pu, --prior-sweep-time low high
                          tau, uniform (give a sweep option too; its tau is
                          not used)
  -PuA, --prior-beneficial-mutation-rate low high
                          uA, uniform
  -Px, --prior-sweep-position low high
                          position of the selected site, uniform
  -Pf, --prior-initial-frequency low high
                          f0, uniform
  -Pc, --prior-final-frequency low high
                          final frequency of a partial sweep, uniform
  -Pe1, --prior-size-change lowTime highTime lowSize highSize
                          a size change of population 0: time and size,
                          uniform
  -Pe2, --prior-size-change lowTime highTime lowSize highSize
                          a second such size change

input and output:
  -Y, --config file       read the model from a YAML file; a demes model
                          goes in demography.demes_filename
  -ts, --tree-sequence file
                          write tree sequences (tskit) to file; with several
                          replicates, to file_rep1, file_rep2, ...
  -F, --full-arg          with -ts (and after it): the full ARG, with the
                          nodes and edges of recombination
  -d, --seed seed1 seed2  seeds of the random number generator
  -v, --verbose           print progress and demes details to stderr
  --help                  print this text
  --version               print the version

Documentation: https://discoal.readthedocs.io
```

