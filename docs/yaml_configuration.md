# YAML Configuration Interface

discoal accepts a YAML configuration file as an alternative to the
positional + flag-based command line. YAML is convenient for
multi-population models, sweep parameter sets, and any case where
you want the simulation specification to live next to the data
rather than in a shell history.

The YAML interface is parsed by `src/core/configInterface.c`. The
authoritative options index — including every accepted field name —
is the runnable file `config_examples/all_options.yaml`.

## Usage

```bash
discoal -Y config.yaml
```

`-Y` is mutually exclusive with the positional command line.
Mixing them is rejected at startup.

## Top-level structure

Every configuration is a mapping with five named blocks:

```yaml
simulation:   { ... }   # required
genetics:     { ... }   # required
demography:   { ... }   # optional (default: single panmictic population)
selection:    { ... }   # optional (default: neutral)
output:       { ... }   # optional (default: ms-style haplotype text)
```

The blocks below are documented in turn.

## simulation

```yaml
simulation:
  sample_size: 20            # required, total chromosomes sampled
  num_replicates: 2          # required, independent replicates
  num_sites: 10000           # required, locus length in sites
  seed: [12345, 67890]       # optional, two-integer pair for reproducibility
```

If `seed` is omitted, discoal seeds from `/dev/urandom` per the
usual command-line behavior.

`sample_size`, `num_replicates`, and `num_sites` must each be
strictly positive.

## genetics

```yaml
genetics:
  mutation_rate: 20.0                # required, theta = 4 N mu L
  recombination_rate: 15.0           # required, rho = 4 N r L
  gene_conversion_rate: 1.0e-1       # optional, gamma
  gene_conversion_tract_length: 100  # optional, mean tract length in sites
  crossover_ratio: 0.5               # optional, alternative to gene_conversion_rate
```

Gene conversion can be specified in either of two ways:

- absolute rate: set `gene_conversion_rate` and
  `gene_conversion_tract_length`. Equivalent to `-g <rate> <tract>`.
- ratio mode: set `crossover_ratio` and
  `gene_conversion_tract_length`. The conversion rate is then
  computed as `rho * crossover_ratio`. Equivalent to `-gr <ratio>
  <tract>`.

Setting both `gene_conversion_rate` and `crossover_ratio` is
rejected. All four numeric fields must be `>= 0`; zero is a legal
value for `mutation_rate` and `recombination_rate` (giving a
neutral or non-recombining model).

## demography

Used for any model with more than one population, with population
size changes, with migration, or that loads an external demes
file.

```yaml
demography:
  deme_sample_size: [10, 5, 5]            # required if demography is set
  effective_population_size: 10000        # optional, sweep effective N (matches CLI -N)
  migration_matrix:                       # optional, constant matrix; off-diagonals are 4Nm
    - row: [0.0, 1.0e-1, 1.0e-1]
    - row: [1.0e-1, 0.0, 1.0e-1]
    - row: [1.0e-1, 1.0e-1, 0.0]
  demographic_events:                     # optional
    population_size_changes:
      - { time: 0.5, size: 0.1, population: 0 }
      - { time: 1.0, size: 2.0, population: 1 }
    migration_rate_changes:
      - { time: 0.2, rate: 0.5, source: 1, destination: 2 }
    ancient_samples:
      - { time: 0.1, population: 2, sample_size: 2 }
    population_splits:
      - { time: 2.0, derived: 1, ancestral: 0 }
      - { time: 3.0, derived: 2, ancestral: 0 }
  demes_filename: "config_examples/example.demes.yaml"  # optional, see below
```

Notes:

- The elements of `deme_sample_size` must each be `>= 0` and sum
  to exactly `simulation.sample_size`. Mismatch is the most common
  configuration error in multi-population YAMLs.
- `effective_population_size` is the same value the CLI's `-N`
  flag sets — it is the effective N used during sweep simulation,
  not a generic per-population size scaler. It must be a
  non-negative integer; fractional or out-of-range values are
  rejected.
- `migration_matrix` rows must be the same length as the number of
  demes. The diagonal must be exactly zero; off-diagonal entries
  must be `>= 0`.
- A `population_split` is interpreted backward in time: at the
  given time the `derived` population merges into the `ancestral`
  one, and all migration into and out of the `derived` population
  stops. This is `-ed time derived ancestral` on the CLI plus
  `-em time i derived 0` and `-em time derived i 0` for every other
  population `i`. (A CLI `-ed` alone, as in ms, keeps migration into
  the merged population, so with constant migration, lineages can be
  stranded there; discoal then stops with an error.)
  `time` must be `> 0` and the two population indices must differ.
- For `population_size_changes`, `time` must be `>= 0` and `size`
  must be `> 0` (extinction at zero is not modeled).
- An `admixture` (in `admixtures`, for example
  `{time: 0.1, population: 0, source_1: 1, source_2: 2, proportion: 0.3}`)
  is `-ea time population source_1 source_2 proportion`:
  backward in time, each lineage of `population` moves to `source_1`
  with probability `proportion` and to `source_2` otherwise. Like a
  split, a YAML admixture ends `population`: all migration into and
  out of it stops at that time (a CLI `-ea` keeps it, as in ms). The
  three populations must differ, the sources must still exist, a
  population is admixed once, `time` must be `> 0` and `proportion`
  in `[0, 1]`.
- A `migration_rate_change` sets the rate (4Nm) at which lineages in
  `source` move to `destination` backward in time, from `time` into
  the past, until a later change. This is
  `-em time source destination rate` on the CLI. `time` and `rate`
  must be `>= 0`, and `source` and `destination` must differ. A
  change that involves a population at or after the time it merges
  into another one (a `population_split` with it as `derived`) is
  rejected: the split already stops its migration.
- A rate change to zero can leave lineages in populations that never
  exchange migrants again and never merge. discoal then stops with an
  error that names the stranded populations. To avoid this, keep a
  connecting migration rate, or add a `population_split` that merges
  the populations.
- An `ancient_sample` takes `sample_size` of the population's
  samples at `time` instead of at time 0. This is
  `-A sample_size population time` on the CLI. The ancient samples
  are part of `deme_sample_size`: with `deme_sample_size: [10, 5, 5]`
  and the entry above, population 2 has 3 samples taken at time 0 and
  2 taken at time 0.1. The ancient samples of a population must not
  outnumber its samples. `time` must be `>= 0` and `sample_size`
  `>= 1`. A sample from a population at or after the time it merges
  into another one is rejected.
- In the output, the samples of each population come in population
  order; within a population, the ancient samples come first (the
  earliest sampling time first), then the samples taken at time 0.
- An ancient sample taken while a sweep is in progress (in population
  0) carries the beneficial allele with the allele's frequency at the
  sampling time; with a partial sweep, its carriers carry the selected
  SNP.
- Events at the same time take effect in this order:
  `ancient_samples`, `migration_rate_changes`,
  `population_size_changes`, `population_splits`, `admixtures` (each
  split and admixture followed by the migration changes that stop the
  ended population's migration), then the size changes with priors,
  and within each list in file order. On the CLI, events
  at the same time take effect in command-line order, so an equivalent
  command line must list them in this order to give the same output.
- With `demes_filename`, `migration_matrix` must be omitted, and
  `demographic_events` may contain only `ancient_samples`: the demes
  file gives the demography. See "Combining with the demes format".


## selection

```yaml
selection:
  sweep_mode: "stochastic"            # one of "stochastic" / "deterministic" / "neutral" / "trajectory_file"
  selection_coefficient: 10.0         # alpha = 2Ns (-a); not with "trajectory_file"
  trajectory_file: "traj.txt"         # with sweep_mode "trajectory_file" only (-S)
  sweep_position: 0.5                 # optional, in [0, 1] (-x; default 0.5)
  fixation_time_ago: 1.0              # single sweep: tau, >= 0 (-ws/-wd/-wn)
  initial_frequency: 1.0e-1           # optional, soft sweep from standing variation (-f)
  final_frequency: 0.9                # optional, partial sweep (-c)
  beneficial_mutation_rate: 1.0e-2    # optional, soft sweep from recurrent mutation (-uA)
  left_rho: 20                        # optional, sweep outside the locus at 4Nr = 2 x left_rho (-ls/-ld/-ln)
  sweep_time_step: 40                 # optional, trajectory time step 1/(step N) (-i; default 40)
  recurrent_sweep_rate: 0.1           # optional, recurrent hitchhiking (-R)
  recurrent_sweeps_left_of_locus: true  # optional, with recurrent_sweep_rate (-L)
```

`sweep_mode` maps directly to the CLI sweep family (`-ws`, `-wd`,
`-wn`). A single sweep needs `fixation_time_ago` (or a prior on it,
see "priors"); `fixation_time_ago: 0` is a sweep that ends at the
present, as `-ws 0`. For recurrent hitchhiking, set
`recurrent_sweep_rate` instead and omit `fixation_time_ago`.

- `left_rho` puts the selected site outside the locus, at genetic
  distance `4Nr = 2 x left_rho` to its left (that is, `left_rho` is
  2Nr): `-ls tau left_rho` (or `-ld`, `-ln`, by `sweep_mode`), with
  the same long-standing factor of 2. `left_rho: 0` puts it at the
  left end of the locus. There is then no `sweep_position`.
- `recurrent_sweeps_left_of_locus: true` with `recurrent_sweep_rate`
  is `-L rate`: each sweep's distance from the locus (4Nr) is drawn
  uniformly from (0, 4Ns). Without it, recurrent sweeps are at
  random positions inside the locus (`-R rate`).
- Recurrent sweeps are always stochastic, as with `-R` and `-L`, so
  `recurrent_sweep_rate` needs `sweep_mode: "stochastic"`.
- `sweep_mode: "trajectory_file"` with `trajectory_file` is `-S`: the
  sweep follows the trajectory in the file (times in 4N0 generations
  before the present, and frequencies, ending at 0; the sweep ends at
  the first frequency at or below 1/(2N); see `-S` in
  `docs/api_reference.rst`). A first frequency below 1 makes the sweep
  partial. The file gives the sweep's time and strength, so
  `selection_coefficient`, `fixation_time_ago`, `initial_frequency`,
  `final_frequency`, `left_rho`, the recurrent sweep fields, and priors
  on these are rejected. A relative path is relative to the directory
  that discoal runs in.

Validation: `selection_coefficient` is required, except with
`sweep_mode: "trajectory_file"`. It must be `>= 0`, and exactly
zero when `sweep_mode: "neutral"`. `sweep_position` is a fraction
in `[0, 1]`. `initial_frequency` and `final_frequency`, when set,
must be in `[0, 1)`; when both are set,
`initial_frequency < final_frequency`. `final_frequency` (or the high
bound of a prior on it) must be above 1/(2N), one copy of the allele,
with N the `effective_population_size`. `beneficial_mutation_rate`,
`recurrent_sweep_rate`, `left_rho` and `fixation_time_ago` must be
`>= 0`, and `sweep_time_step` `> 0`.

## priors

Each entry replaces a fixed value with a new draw for every
replicate, as the CLI `-P` flags do. The same seeds give the same
draws as the equivalent command line.

```yaml
priors:
  mutation_rate: {low: 5, high: 15}             # theta, uniform (-Pt)
  recombination_rate: {low: 0, high: 50}        # rho, uniform (-Pr)
  # recombination_rate_exponential: {mean: 10, upper: 100}  # or this (-Pre)
  selection_coefficient: {low: 100, high: 1000} # alpha (-Pa)
  fixation_time_ago: {low: 0, high: 0.1}        # tau (-Pu)
  beneficial_mutation_rate: {low: 0, high: 0.01}  # uA (-PuA)
  sweep_position: {low: 0.4, high: 0.6}         # -Px
  initial_frequency: {low: 0.01, high: 0.1}     # f0 (-Pf)
  final_frequency: {low: 0.5, high: 0.9}        # partial sweep (-Pc)
  population_size_changes:                      # -Pe1, -Pe2
    - {time_low: 0.01, time_high: 0.5, size_low: 0.1, size_high: 10,
       population: 0}
```

Notes:

- Uniform priors draw from `[low, high]`, and need `low <= high`
  inside the parameter's range (rates `>= 0`, `sweep_position` in
  `[0, 1]`, `initial_frequency` in `[0, 1)`, `final_frequency`
  inside `(0, 1)`).
- `recombination_rate_exponential` draws rho from an exponential
  distribution with the given mean and sets a draw above `upper` to
  `upper` (a cap, not a truncation), as `-Pre`. It cannot be combined
  with `recombination_rate`.
- Priors on sweep parameters need a selection block. With a prior on
  `fixation_time_ago`, the selection block needs no
  `fixation_time_ago`; on the CLI, `-Pu` still needs a `-ws` value,
  which the draw replaces. A prior on `sweep_position` cannot be
  combined with `left_rho` or `recurrent_sweeps_left_of_locus`.
- Each `population_size_changes` entry is a size change of
  `population` (default 0) at a time uniform on
  `[time_low, time_high]` (units of 4N generations) to a size
  uniform on `[size_low, size_high]` (relative to N). Up to 64
  entries may be given; `-Pe1` and `-Pe2` are the first two, for
  population 0. The bounds must satisfy
  `0 <= time_low <= time_high` and `0 < size_low <= size_high`. The
  size trajectories are checked before the run at the lower bounds,
  and again after each replicate's draws; discoal stops if a draw
  makes a population size reach zero (for example, with a linear
  decline).

## output

```yaml
output:
  output_type: "tree_sequence"         # optional; "haplotype" (default) or "tree_sequence"
  hide_partial_snp: false              # optional; suppresses the selected SNP (matches CLI -h)
  unsimplified_tree_sequence: false    # optional; only meaningful with tree_sequence output
  tree_sequence_filename: "out.trees"  # required if output_type is "tree_sequence"
```

Notes:

- `output_type: "haplotype"` (the default) gives discoal's ms-format
  output on stdout: one `//` block per replicate with `segsites:`,
  `positions:`, and one row per haplotype. `"tree_sequence"` writes
  tree sequence files instead (the CLI's `-ts`). The schema used to
  accept `"snp_array"` and a `finite_output` field, which had no
  effect; both are now rejected with an error.
- When `output_type` is `"tree_sequence"`, discoal writes one
  `.trees` file per replicate, named by stripping the `.trees`
  suffix from `tree_sequence_filename` and appending
  `_repN.trees` for `N = 1..num_replicates`. If
  `unsimplified_tree_sequence: true` the simplification step is
  skipped and each file retains the full ARG.
- `hide_partial_snp: true` is the YAML form of the CLI `-h` flag:
  the selected SNP is excluded from the haplotype output. Useful
  when measuring θ from sweep simulations, where including the
  selected site biases the estimate.

## Units

Times in the YAML file are in **units of 4N generations**, where N
is the reference population size. This is the same convention as
the command line and as ms.

For a user thinking in generations, the conversion is

```
yaml_time = generations / (4 * N_reference)
```

So 40,000 generations with `N_reference = 10,000` gives
`yaml_time = 1.0`. This is the same value you would write on the
command line as `-en 1.0 ...`.

This convention applies to every time field in the schema:
`population_size_changes.time`, `population_splits.time`,
`migration_rate_changes.time`, `ancient_samples.time`,
`fixation_time_ago`, etc.

A demes file does not use this convention: its times are in
generations or years, and discoal converts them itself (see
"Combining with the demes format"). For all units, including the
2N time scale of tree sequence output, see the units table in
`docs/basic_usage.rst`.

Other unit conventions:

- **Population sizes** in `population_size_changes.size` and the
  matrix entries are relative ratios, not absolute counts. A size
  of `0.1` means "this population has 10% of the reference N at
  that time".
- **Migration rates** (`migration_matrix` rows and
  `migration_rate_changes.rate`) are 4Nm.
- **Selection coefficient** is 2Ns, matching `-a`.

## Combining with the demes format

`demography.demes_filename` loads an external demes-format YAML
specifying population sizes, splits, and migration windows.
discoal converts the demes graph into its internal event stream
using the present-day size of the first present-day deme as the
reference N. Demes migration rates (per-generation) are scaled to
4Nm.

When a demes file is in use, the `migration_matrix` block must be
omitted, and `demographic_events` may contain only
`ancient_samples` — the parser rejects other combinations rather
than silently overriding either. For ancient samples with a demes
model:

- `population` is an index into `deme_sample_size`, which follows
  the order of demes in the demes file, extinct demes included.
- `time` is in units of 4 N_ref generations, where N_ref is the
  present-day size of the first deme that exists at time 0:
  `time = generations / (4 * N_ref)`.
- The deme must exist at the sampling time: from its end time
  (included) back to its start time (excluded). discoal stops with
  an error otherwise.
- A deme that is extinct at time 0 can have samples if all of them
  are ancient.

See `config_examples/ancient_demes_example.yaml`.

Note that the demes importer's handling of multi-window migration
is currently coarser than the demes spec strictly intends. If you
hit unexpected results with a demes file containing migration
intervals whose participant set changes over time, reproduce the
same demography with explicit `demographic_events` and
`migration_matrix` to diagnose.

## Worked example

This file is exercised by the validation suite as
`config_examples/demographic_example.yaml`.

```yaml
simulation:
  sample_size: 20
  num_replicates: 2
  num_sites: 10000
  seed: [12345, 10102]

genetics:
  mutation_rate: 1.0e-1
  recombination_rate: 1.0e-1

demography:
  deme_sample_size: [10, 5, 5]
  effective_population_size: 1.0
  migration_matrix:
    - row: [0.0, 1.0e-1, 1.0e-1]
    - row: [1.0e-1, 0.0, 1.0e-1]
    - row: [1.0e-1, 1.0e-1, 0.0]
  demographic_events:
    population_size_changes:
      - { time: 0.5, size: 0.1, population: 0 }
      - { time: 1.0, size: 2.0, population: 1 }
    population_splits:
      - { time: 2.0, derived: 1, ancestral: 0 }
      - { time: 3.0, derived: 2, ancestral: 0 }
```

The equivalent positional command line, also exercised by the
validation suite, is

```
discoal 20 2 10000 -t 0.1 -r 0.1 -p 3 10 5 5 -M 0.1 \
        -en 0.5 0 0.1 -en 1.0 1 2.0 \
        -ed 2.0 1 0 -em 2.0 0 1 0 -em 2.0 1 0 0 -em 2.0 2 1 0 -em 2.0 1 2 0 \
        -ed 3.0 2 0 -em 3.0 0 2 0 -em 3.0 2 0 0 -em 3.0 1 2 0 -em 3.0 2 1 0 \
        -d 12345 10102
```

## Validation

The repository ships a parity-check suite at
`testing/yaml_validation_suite.sh`. For each runnable fixture
under `config_examples/`, the suite reads a `# Equivalent to:
discoal …` header, runs both the YAML form and the CLI form with
the same seeds, and confirms that the output is identical. ms
fixtures are diffed on stdout (modulo discoal's own argv echo on
line 1); tree-sequence fixtures are compared as
`tskit.TableCollection`s with provenance ignored, since discoal
stamps a timestamp into each `.trees` file's provenance table.
