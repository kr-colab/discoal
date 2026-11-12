# YAML Configuration

A YAML file can give all the parameters of a discoal run, in place of
the command line.

## Usage

```bash
discoal -Y config.yaml
discoal -Y config.yaml -d 1 2      # the seeds 1 2, not those of the file
```

- With `-Y`, the positional arguments (sample size, replicates, sites)
  come from the file. discoal stops with an error if you also give them
  on the command line.
- Other command-line options, before or after `-Y`, apply after the
  file:
  - an option that sets a value (for example `-t`, `-r`, `-d`, `-N`)
    replaces the value of the file;
  - an event option (for example `-en`, `-em`, `-A`) adds an event to
    the events of the file;
  - `-m i j rate` sets one entry of the migration matrix;
  - `-p` replaces `deme_sample_size`. The sum must still be equal to
    `sample_size`.
  - With a demes file, options that change the demography or the
    samples (`-p`, `-A`, `-m`, `-M`, all `-e` options, `-Pe1`, `-Pe2`)
    are errors. The demes file gives the demography, and
    `deme_sample_size` and `ancient_samples` the samples.
- `config_examples/` has examples that run. Each one starts with a line
  `# Equivalent to: discoal ...`, which is the same run on the command
  line. `config_examples/all_options.yaml` lists every field, but it
  does not run, because it gives some fields more than once.

## Structure

A file has up to six blocks:

```yaml
simulation:  { ... }   # required
genetics:    { ... }   # required
demography:  { ... }   # optional; default: one population
selection:   { ... }   # optional; default: no sweep
priors:      { ... }   # optional; parameters drawn for each replicate
output:      { ... }   # optional; default: ms-style output on stdout
```

Each example below is a block. To run it, put it in a file with a
`simulation` block and a `genetics` block.

All times are in units of 4N generations, as on the command line and in
ms (see "Units").

## simulation

```yaml
simulation:
  sample_size: 20          # required; chromosomes sampled, 1 to 65535
  num_replicates: 2        # required; > 0
  num_sites: 10000         # required; sites in the locus, 1 to 10^8
  seed: [12345, 67890]     # optional; two integers (-d)
```

Without `seed`, discoal reads its seeds from `/dev/urandom`.

## genetics

```yaml
genetics:
  mutation_rate: 20.0                # required; theta = 4N mu L (-t)
  recombination_rate: 15.0           # required; rho = 4N r L (-r)
  gene_conversion_rate: 10.0         # optional (-g)
  gene_conversion_tract_length: 100  # mean tract length in sites
```

- All values must be `>= 0`. A mutation or recombination rate of 0 is
  allowed.
- Gene conversion has two forms. Set one of them, not both:
  - `gene_conversion_rate`: the rate for the whole locus, 4N g L, scaled
    as rho is (`-g rate length`).
  - `crossover_ratio`: the rate is `recombination_rate x
    crossover_ratio` (`-gr ratio length`).
- With either form, `gene_conversion_tract_length` must be at least 1.
  The tract lengths are geometric, with this mean.

## demography

Use this block for more than one population, for size changes or
migration, or for a demes file.

```yaml
demography:
  deme_sample_size: [10, 5, 5]       # required; samples of each population
  effective_population_size: 10000   # optional; N (-N), default 10^6
  migration_matrix:                  # optional; 4Nm, diagonal 0
    - row: [0.0, 0.1, 0.1]
    - row: [0.1, 0.0, 0.1]
    - row: [0.1, 0.1, 0.0]
  demographic_events:                # optional
    population_size_changes:
      - {time: 0.5, size: 0.1, population: 0}
    exponential_growth_changes:
      - {time: 0.0, rate: 2.0, population: 1}
    migration_rate_changes:
      - {time: 0.2, rate: 0.5, source: 1, destination: 2}
    ancient_samples:
      - {time: 0.1, population: 2, sample_size: 2}
    population_splits:
      - {time: 2.0, derived: 1, ancestral: 0}
      - {time: 3.0, derived: 2, ancestral: 0}
```

A model with one population and a sweep that needs `-N` also uses this
block, with `deme_sample_size: [n]`.

### Samples and N

- `deme_sample_size` has one entry for each population (at most 121).
  Each entry must be `>= 0`, and the sum must be equal to
  `simulation.sample_size`.
- `effective_population_size` is N, as `-N` sets it. Only sweeps use
  N. It must be an integer. 0 is the same as no value, and gives the
  default, 10^6.

### Migration

- `migration_matrix` has one row for each population. The diagonal
  must be 0, and the other entries `>= 0`. The entry in row i, column
  j is `-m i j`: the rate at which lineages in population i move to
  population j, backward in time.
- A `migration_rate_change` is `-em time source destination rate`: from
  `time` into the past, lineages in `source` move to `destination` at
  `rate`. `time` and `rate` must be `>= 0`, and the two populations
  must be different.
- If a rate change to 0 leaves lineages in populations that can never
  meet again, discoal stops with an error that names these
  populations. To prevent this, keep a migration rate between them, or
  merge them with a `population_split`.

### Size changes

A `population_size_change` is `-en time population size`. `size` is a
fraction of N, and must be `> 0`. `time` must be `>= 0`. `population`
is required.

### Growth

- An entry of `exponential_growth_changes` is `-eg time population
  rate`: from `time` back into the past, the size of `population`
  changes exponentially at `rate` per 4N generations, from its size at
  `time`. A positive rate makes the population smaller in the past.
- An entry of `linear_growth_changes` is `-el time population rate`:
  the size changes linearly, by `rate` (a fraction of N) per 4N
  generations. A positive rate makes the population smaller in the
  past, and discoal stops with an error if the size reaches zero.
- `time` must be `>= 0`, and `rate` a finite number. The growth goes on
  until the next size change or growth entry of the population; a rate
  of 0 stops it. A size change resets the growth, as `-en` does.
- For `-eG` or `-eL` (all populations), give one entry for each
  population. In the same way, `-eM` is one `migration_rate_changes`
  entry for each pair of populations.
- A simulated sweep trajectory does not follow growth, so discoal stops
  with an error if population 0 grows at or after the start of a single
  sweep (a growth that a size change ends before the sweep is allowed),
  or if any population grows with recurrent sweeps. A sweep from a
  `trajectory_file` follows the sizes, and is not restricted.

For example, the size of population 0 grows tenfold over the last 0.23
units of time (rate ln(10) / 0.23), and is constant before that:

```yaml
    exponential_growth_changes:
      - {time: 0.0, rate: 10.0, population: 0}
    population_size_changes:
      - {time: 0.23, size: 0.1, population: 0}
```

### Splits and admixtures

- A `population_split` is a merge, backward in time: at `time`, the
  lineages of `derived` move into `ancestral`. All migration into and
  out of `derived` then stops. This is `-ed time derived ancestral`,
  with an `-em time ... 0` for each migration rate of `derived`. A
  command-line `-ed` alone keeps these rates, as in ms. `time` must be
  `> 0`, and the two populations must be different.
- An `admixture`, in a list `admixtures` in `demographic_events`, is
  `-ea time population source_1 source_2 proportion`. Backward in
  time, each lineage of `population` moves to `source_1` with
  probability `proportion`, and to `source_2` if not. All migration
  into and out of `population` then stops, as for a split. A
  command-line `-ea` alone keeps the migration rates, as in ms. Rules:
  - the three populations must be different;
  - the two sources must not have ended;
  - a population can have one admixture only;
  - `time` must be `> 0`;
  - `proportion` must be in `[0, 1]`.

  For example:

  ```yaml
      admixtures:
        - {time: 0.1, population: 0, source_1: 1, source_2: 2, proportion: 0.3}
  ```

- A split or an admixture ends a population. At that time and after
  it, a migration rate change or an ancient sample of that population
  is an error.

### Ancient samples

- An `ancient_sample` is `-A sample_size population time`: it takes
  `sample_size` of the samples of `population` at `time`, not at time
  0. With the example above, population 2 has 3 samples from time 0
  and 2 from time 0.1.
- The ancient samples of a population must not be more than its
  samples. `time` must be `>= 0`, and `sample_size` `>= 1`.
- In the output, the samples come in population order. In a
  population, the ancient samples come first, the most recent first,
  then the samples from time 0.
- An ancient sample taken during a sweep in population 0 carries the
  beneficial allele with the allele's frequency at that time.

### Events at the same time

The order of events at the same time can change the output. Events at
the same time take effect in this order:

1. `ancient_samples`;
2. `migration_rate_changes`;
3. `population_size_changes`;
4. `exponential_growth_changes`, then `linear_growth_changes`, so that
   a growth entry starts from the new size of a size change at the same
   time;
5. `population_splits`, then `admixtures`, each followed by the
   migration changes that end the population's migration;
6. size changes drawn from priors;
7. the sweep;
8. events from command-line options.

In each list, the order is the file order. On the command line, events
at the same time take effect in the order of the options. To get the
same output from the command line, give the options in the order above.

## selection

```yaml
selection:
  sweep_mode: "stochastic"        # required; see below
  selection_coefficient: 200      # alpha = 2Ns (-a)
  sweep_position: 0.5             # optional; in [0, 1] (-x), default 0.5
  fixation_time_ago: 0.05         # tau (-ws, -wd, -wn)
  initial_frequency: 0.1          # optional; soft sweep from standing variation (-f)
  final_frequency: 0.9            # optional; partial sweep (-c)
  beneficial_mutation_rate: 0.5   # optional; soft sweep from recurrent mutation (-uA)
  sweep_time_step: 40             # optional; trajectory step 1/(step x N) (-i), default 40
```

`sweep_mode` is required. It is one of:

| `sweep_mode` | sweep | option |
|---|---|---|
| `"stochastic"` | stochastic trajectory | `-ws` |
| `"deterministic"` | deterministic trajectory | `-wd` |
| `"neutral"` | neutral fixation | `-wn` |
| `"trajectory_file"` | trajectory read from a file | `-S` |

### Single sweeps

- A single sweep needs `fixation_time_ago`, or a prior on it.
  `fixation_time_ago: 0` is a sweep that ends at the present (`-ws 0`).
- A run has one single sweep at most. With a single sweep in the file,
  a sweep option on the command line (`-ws`, `-wd`, `-wn`, `-ls`,
  `-ld`, `-ln`, `-S`, `-R`, `-L`) is an error.
- `left_rho` puts the selected site outside the locus, to its left, at
  `4Nr = 2 x left_rho`. This is `-ls tau left_rho` (or `-ld`, `-ln`),
  with the same factor of 2. `left_rho: 0` puts the site at the left end
  of the locus. Do not give `sweep_position` with `left_rho`.
- With `sweep_mode: "trajectory_file"`, give `trajectory_file:
  "path"` (`-S`). The file gives the times (in 4N generations) and
  frequencies of the sweep; see `-S` in `docs/api_reference.rst`. If
  the first frequency is below 1, the sweep is partial. A relative
  path is relative to the directory in which discoal runs.
  The file gives the sweep's time and strength, so
  `selection_coefficient`, `fixation_time_ago`, `initial_frequency`,
  `final_frequency`, `left_rho`, the recurrent sweep fields, and priors
  on these are errors.

### Recurrent sweeps

```yaml
selection:
  sweep_mode: "stochastic"
  selection_coefficient: 500
  recurrent_sweep_rate: 2.0                # -R
  recurrent_sweeps_left_of_locus: false    # optional; true is -L
```

- Recurrent sweeps happen at `recurrent_sweep_rate` per 2N individuals
  per generation, at random times. Do not give `fixation_time_ago`.
- Without `recurrent_sweeps_left_of_locus`, the sweeps are at random
  positions in the locus (`-R`). With it, they are to the left of the
  locus (`-L`), at a distance (4Nr) drawn uniformly from (0, 4Ns).
- Recurrent sweeps are always stochastic, so `sweep_mode` must be
  `"stochastic"`.
- `initial_frequency`, or a prior on it, and a prior on
  `fixation_time_ago` are errors with recurrent sweeps.

### Rules

- `selection_coefficient` is required, except with
  `"trajectory_file"`. It must be `>= 0`, and 0 with `"neutral"`.
- `initial_frequency` and `final_frequency` must be in `[0, 1)`. If
  both are given, `initial_frequency` must be less than
  `final_frequency`. `final_frequency`, or the high bound of a prior on
  it, must be more than 1/(2N), one copy of the allele.
- `beneficial_mutation_rate`, `recurrent_sweep_rate`, `left_rho` and
  `fixation_time_ago` must be `>= 0`, and `sweep_time_step` `> 0`.
- The rules of the command line also apply. For example:
  - a deterministic sweep with a population size change is an error;
  - with a sweep, each size of population 0, at any time, must be at
    least one individual (N x size >= 1);
  - population 0 must not grow or shrink continuously during a sweep
    (see `docs/selection.rst`).

## priors

Each entry replaces a fixed value with a value drawn for each
replicate, as the `-P` options do. The same seeds give the same values
as the equivalent command line. Priors on sweep parameters need a
`selection` block.

```yaml
priors:
  mutation_rate: {low: 5, high: 15}               # theta (-Pt)
  recombination_rate: {low: 0, high: 50}          # rho (-Pr)
  selection_coefficient: {low: 100, high: 1000}   # alpha (-Pa)
  fixation_time_ago: {low: 0, high: 0.1}          # tau (-Pu)
  beneficial_mutation_rate: {low: 0, high: 0.01}  # uA (-PuA)
  sweep_position: {low: 0.4, high: 0.6}           # -Px
  initial_frequency: {low: 0.01, high: 0.1}       # f0 (-Pf)
  final_frequency: {low: 0.5, high: 0.9}          # -Pc
  population_size_changes:                        # -Pe1, -Pe2
    - {time_low: 0.01, time_high: 0.5, size_low: 0.1, size_high: 10, population: 0}
```

- A uniform prior draws from `[low, high]`. `low` must be less than or
  equal to `high`, and both must be in the range of the parameter.
  With a soft sweep that is also partial, every value of
  `initial_frequency` must be less than every value of
  `final_frequency`.
- `recombination_rate_exponential: {mean: 10, upper: 100}` (`-Pre`)
  draws rho from an exponential distribution with this mean, and
  replaces a value above `upper` with `upper`. Do not give it with
  `priors.recombination_rate`.
- With a prior on `fixation_time_ago`, the `selection` block needs no
  `fixation_time_ago`. On the command line, `-Pu` still needs a `-ws`
  value, which the drawn value replaces.
- A prior on `sweep_position` is an error with `left_rho` or with
  recurrent sweeps to the left of the locus.
- Each entry of `population_size_changes` is a size change of
  `population` (default 0), at a time drawn from
  `[time_low, time_high]`, to a size drawn from
  `[size_low, size_high]`. `-Pe1` and `-Pe2` are the first two entries,
  for population 0. You can give up to 64 entries. The bounds must be
  `0 <= time_low <= time_high` and `0 < size_low <= size_high`.
  With a linear size change (`-el`), a drawn value can make a
  population size go to 0 before the change ends. discoal then stops
  with an error.

## output

```yaml
output:
  output_type: "tree_sequence"         # "haplotype" (default) or "tree_sequence"
  tree_sequence_filename: "out.trees"  # required with "tree_sequence"
  unsimplified_tree_sequence: false    # optional; with "tree_sequence" (-F)
  hide_partial_snp: false              # optional (-h)
```

- `"haplotype"` writes ms-style output on stdout: for each replicate,
  `//`, `segsites:`, `positions:` and one row for each sample.
- `"tree_sequence"` writes tree sequences (`-ts`), not ms-style
  output. With one replicate, the file name is
  `tree_sequence_filename`. With more, discoal puts `_repN` before the
  last `.` of the name, for replicates N = 1, 2, ...: `out.trees`
  gives `out_rep1.trees`, `out_rep2.trees`, ... A name with no `.`
  gets `_repN` at its end.
- `tree_sequence_filename` without `output_type: "tree_sequence"` is
  an error.
- `unsimplified_tree_sequence: true` writes the full ancestral
  recombination graph.
- `hide_partial_snp: true` removes the selected site of a partial sweep
  from the output.

## Units

| quantity | unit | option |
|---|---|---|
| times | 4N generations | `-en`, `-ws`, ... |
| population sizes | fraction of N | `-en` |
| migration rates | 4Nm | `-m`, `-em` |
| selection | alpha = 2Ns | `-a` |
| mutation, recombination | 4N mu L, 4N r L | `-t`, `-r` |

For example, 40,000 generations with N = 10,000 is a time of
40,000 / (4 x 10,000) = 1.0. A demes file uses generations or years,
and discoal converts them. For all units, including the time unit of
tree sequence output, see the table in `docs/basic_usage.rst`.

## Demes files

```yaml
demography:
  deme_sample_size: [10, 10]
  demes_filename: "config_examples/demes_example.demes.yaml"
```

- A demes file gives the demography. To convert its times, sizes and
  rates, discoal uses as N the present-day size of the first deme that
  exists at time 0. See `docs/demes_integration.md`.
- This N converts the demes file only. A sweep uses
  `effective_population_size`. With a sweep, set it to the same value.
- `deme_sample_size` follows the order of the demes in the file,
  extinct demes included.
- With a demes file, do not give `migration_matrix`.
  `demographic_events` can contain only `ancient_samples`. For these:
  - `time` is in 4N generations, with the N above;
  - a sample is taken after the events of the demes file at the same
    time, as in msprime: a pulse into the deme at that time does not
    move the sample;
  - the deme must exist at that time, from its end time (included) to
    its start time (not included);
  - a deme that does not exist at time 0 can have samples if all of
    them are ancient.

  See `config_examples/ancient_demes_example.yaml`.

## Validation

`testing/yaml_validation_suite.sh` runs each example in
`config_examples/` from the YAML file and from its `# Equivalent to:`
command line, with the same seeds, and checks that the outputs are
identical. It compares ms-style output without its first line (the
command line). It compares tree sequences with tskit, without their
provenance, which has a time stamp.
