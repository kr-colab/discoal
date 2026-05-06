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
  one. Equivalent to `-ed time derived ancestral` on the CLI.
  `time` must be `> 0` and the two population indices must differ.
- For `population_size_changes`, `time` must be `>= 0` and `size`
  must be `> 0` (extinction at zero is not modeled).
- `demes_filename` and explicit `demographic_events` /
  `migration_matrix` are mutually exclusive — use one or the other.

Time-varying migration *is* supported: use the CLI primitives
`-em time popID1 popID2 rate` (single pair) or
`-eM time rate` (all off-diagonal pairs), or load a demes file
with multi-window migration via `-D` / `demes_filename`.

The `demographic_events` block also accepts the field name
`ancient_samples`, which is part of the schema but is **not yet
wired into the simulator** and will exit with an error if set.

## selection

```yaml
selection:
  sweep_mode: "stochastic"            # one of "stochastic" / "deterministic" / "neutral"
  selection_coefficient: 10.0         # alpha = 2Ns
  sweep_position: 0.5                 # in [0, 1]
  fixation_time_ago: 1.0              # required unless using recurrent_sweep_rate
  initial_frequency: 1.0e-1           # optional, triggers soft sweep from standing variation
  final_frequency: 0.9                # optional, triggers partial sweep
  beneficial_mutation_rate: 1.0e-2    # optional, soft sweep from recurrent mutation (uA)
  recurrent_sweep_rate: 0.1           # optional, recurrent hitchhiking; mutually exclusive with fixation_time_ago
```

`sweep_mode` maps directly to the CLI sweep family (`-ws`, `-wd`,
`-wn`). For a single sweep, `fixation_time_ago` is required. For
recurrent hitchhiking, set `recurrent_sweep_rate` instead and omit
`fixation_time_ago`.

Validation: `selection_coefficient` must be `>= 0`, and exactly
zero when `sweep_mode: "neutral"`. `sweep_position` is a fraction
in `[0, 1]`. `initial_frequency` and `final_frequency`, when set,
must be in `[0, 1)`; when both are set,
`initial_frequency < final_frequency`. `beneficial_mutation_rate`
and `recurrent_sweep_rate` must be `>= 0`.

## output

```yaml
output:
  output_type: "tree_sequence"         # see notes; "haplotype" / "snp_array" are currently no-ops
  hide_partial_snp: false              # optional; suppresses the selected SNP (matches CLI -h)
  unsimplified_tree_sequence: false    # optional; only meaningful with tree_sequence output
  tree_sequence_filename: "out.trees"  # required if output_type is "tree_sequence"
  finite_output: false                 # accepted but currently has no effect; see notes
```

Notes:

- Only `output_type: "tree_sequence"` actually changes simulator
  behavior. Setting `"haplotype"` or `"snp_array"` is accepted by
  the parser but currently has no effect — in either case discoal
  emits its default ms-style stdout (one `//` block per replicate
  with `segsites:`, `positions:`, and per-haplotype rows).
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
- `finite_output` is in the schema but is currently a stub: the
  parser sets a global flag that nothing in the simulator
  consults. Setting it has no effect today; treat it as reserved.

## Units

discoal's YAML and CLI both expose **2N units of time**, scaled
internally by a factor of 2 to land in 4N coalescent units. This
is the same convention used by the positional CLI flags.

For a user thinking in generations, the conversion is

```
yaml_time = generations / (2 * N_reference)
```

So 40,000 generations with `N_reference = 10,000` gives
`yaml_time = 2.0`. This is the same value you would write on the
command line as `-en 2.0 ...`.

This convention applies to every time field in the schema:
`population_size_changes.time`, `population_splits.time`,
`fixation_time_ago`, etc.

Other unit conventions:

- **Population sizes** in `population_size_changes.size` and the
  matrix entries are relative ratios, not absolute counts. A size
  of `0.1` means "this population has 10% of the reference N at
  that time".
- **Migration rates** (`migration_matrix` rows) are 4Nm.
- **Selection coefficient** is 2Ns, matching `-a`.

## Combining with the demes format

`demography.demes_filename` loads an external demes-format YAML
specifying population sizes, splits, and migration windows.
discoal converts the demes graph into its internal event stream
using the present-day size of the first present-day deme as the
reference N. Demes migration rates (per-generation) are scaled to
4Nm.

When a demes file is in use, the `demographic_events` and
`migration_matrix` blocks must be omitted — the parser rejects
the combination rather than silently overriding either.

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
        -ed 2.0 1 0   -ed 3.0 2 0 \
        -d 12345 10102
```

## Limitations

The YAML interface does not currently expose:

- prior distributions (`-Pt`, `-Pa`, `-Pu`, `-Px`, etc.);
- linked-locus sweeps (`-ls`, `-ld`, `-ln`);
- conditional simulation (`-C`);
- ancient samples and time-varying migration rates (the schema
  accepts these field names but the parser rejects them at
  runtime — see the demography section).

Use the positional command-line interface for these. The CLI
reference is in `docs/basic_usage.rst`,
`docs/population_structure.rst`, `docs/selection.rst`, and
`docs/advanced_features.rst`.

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
