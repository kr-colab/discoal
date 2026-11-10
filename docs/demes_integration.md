# Demes Integration in Discoal

Discoal supports the [demes](https://popsim-consortium.github.io/demes-spec/) specification for describing demographic models. This allows you to use standardized demographic models that are compatible with other simulation tools like msprime, SLiM, and fwdpy11.

## Usage

To use a demes file with discoal:

```bash
discoal <sample_size> <num_replicates> <sequence_length> -D model.yaml [other options]
```

Or within a YAML configuration file:

```yaml
events:
  demes_file: model.yaml
```

## Supported Demes Features

### Supported
- Multiple populations with specified sizes
- Population size changes (instantaneous)
- Exponential growth epochs (translated to discoal `-eg` events). As
  with `-eg`, a sweep cannot overlap an exponential epoch of the first
  deme (population 0), and recurrent sweeps cannot be combined with any
  exponential epoch; see "Combining with Demography" in
  `docs/selection.rst`.
- Population splits and mergers
- Symmetric and asymmetric migration, including multi-window
  piecewise-constant migration matrices
- Pulse migration events
- Sample size specification per population

### Not Supported
- Selfing rates
- Cloning rates

If your demes file contains unsupported features, discoal will report an error and exit.

## Parameter Conversion

Discoal uses different time and rate units than the demes specification. The integration handles these conversions automatically:

### Time Conversion
- Demes: time in the units the graph specifies (`time_units: generations`
  or `time_units: years`).
- Discoal command line: time in units of 4N generations (as in ms).

```
cli_time = demes_time / generation_time / (4 * N_ref)
```

Where `N_ref` is the size of the first present-day deme in the demes file.
For `time_units: generations`, the demes spec requires `generation_time = 1`.
For `time_units: years`, `generation_time` is the number of years per
generation. (Internally discoal uses units of 2N generations, so the importer
stores twice `cli_time`.)

Example: with `N_ref = 10000`, a split 40000 generations ago is
`-ed 1.0 ...` on the command line.

### Population Size Conversion
- Demes: absolute population size
- Discoal: size relative to reference N
- Conversion: `discoal_size = demes_size / N`

### Migration Rate Conversion
- Demes: proportion of migrants per generation
- Discoal: scaled migration rate 4Nm
- Conversion: `discoal_rate = demes_rate * 4N`

## Example Demes Files

### Single Population with Size Changes

```yaml
description: Single population with bottleneck and recovery
time_units: generations
demes:
  - name: pop0
    epochs:
      - start_time: .inf
        end_time: 20000
        size_function: constant
        size: 10000
      - start_time: 20000
        end_time: 10000
        size_function: constant
        size: 1000      # Bottleneck
      - start_time: 10000
        end_time: 0
        size_function: constant
        size: 50000     # Recovery
```

### Two Populations with Migration

```yaml
description: Two populations with symmetric migration
time_units: generations
demes:
  - name: pop0
    epochs:
      - start_time: .inf
        end_time: 0
        size_function: constant
        size: 10000
  - name: pop1
    epochs:
      - start_time: 50000
        end_time: 0
        size_function: constant
        size: 10000
        ancestors: [pop0]
migrations:
  - demes: [pop0, pop1]
    rate: 0.0001
    start_time: 50000
    end_time: 0
```

### Population Split with Asymmetric Migration

```yaml
description: Population split with asymmetric gene flow
time_units: generations
demes:
  - name: ancestral
    epochs:
      - start_time: .inf
        end_time: 100000
        size_function: constant
        size: 20000
  - name: popA
    epochs:
      - start_time: 100000
        end_time: 0
        size_function: constant
        size: 15000
        ancestors: [ancestral]
  - name: popB
    epochs:
      - start_time: 100000
        end_time: 0
        size_function: constant
        size: 5000
        ancestors: [ancestral]
migrations:
  - source: popA
    dest: popB
    rate: 0.0002
    start_time: 100000
    end_time: 50000
  - source: popB
    dest: popA
    rate: 0.0001
    start_time: 100000
    end_time: 0
```

## Combining Demes with Other Parameters

When using a demes file, you can still specify other evolutionary parameters:

```bash
# Add mutation and recombination to a demes model
discoal 100 1000 50000 -D model.yaml -t 20 -r 15

# Add selection to a demes model
discoal 100 1000 50000 -D model.yaml -t 20 -r 15 -ws 0.01 -a 100 -x 0.5
```

## Sample Size Specification

Give the number of samples per deme in `demography.deme_sample_size` of the
`-Y` config. The list follows the order of demes in the demes file,
**extinct demes included**. Give 0 for a deme that does not exist at time 0,
or take all its samples as ancient samples (below); discoal stops with an
error if such a deme gets samples at time 0.

For a demes file with demes `ANC` (ancestral, extinct), `A` and `B`:

```yaml
demography:
  deme_sample_size: [0, 50, 50]   # ANC, A, B
  demes_filename: "model.demes.yaml"
```

### Ancient Samples

The `demographic_events` block of the `-Y` config may give
`ancient_samples` together with `demes_filename` (no other event lists).
Each entry takes `sample_size` of a deme's `deme_sample_size` samples at
`time` instead of at time 0. `population` is the deme's index in the demes
file, and `time` is in units of 4 N_ref generations
(`time = generations / (4 * N_ref)`, as in "Time Conversion"). The deme
must exist at that time: from its end time (included) back to its start
time (excluded).

For the demes file above, if A's present size is 10000 (so N_ref = 10000)
and A exists 2000 generations ago, this takes four of A's samples then:

```yaml
demography:
  deme_sample_size: [0, 50, 50]   # ANC, A, B
  demes_filename: "model.demes.yaml"
  demographic_events:
    ancient_samples:
      - {time: 0.05, population: 1, sample_size: 4}
```

`config_examples/ancient_demes_example.yaml` samples an extinct ghost deme.

## Validation and Debugging

### Checking Conversion
Discoal reports the loaded demographic model:

```
Loaded 2 populations and 5 events from demes file 'model.yaml'
```

### Common Issues

1. **Infinite coalescent times**: Ensure all populations are connected through migration or common ancestors
2. **Time conversion confusion**: discoal command-line times are in units of 4N generations. Internal times, and node times in tree sequence output, are in units of 2N generations
3. **Sample size mismatch**: Ensure `-p` sample sizes match the number of present-day demes

### Validation Script

Use the demes validation suite to compare demes and manual specifications:

```bash
testing/demes_validation_suite.sh
```

## Tips for Writing Demes Files

1. **Use descriptive names**: Help future users understand your model
2. **Add descriptions**: Document assumptions and parameter sources
3. **Validate with other tools**: Test your demes file with msprime to ensure correctness
4. **Start simple**: Build complex models incrementally
5. **Check units**: Always specify `time_units: generations`

## Converting Existing Models

### From msprime
Demes files written for msprime work directly with discoal, with two exceptions: linear size epochs (the demes-c parser rejects them) and selfing or cloning rates (not supported).

### From discoal command line
To convert a discoal command line to demes:
1. Convert times: `demes_time = discoal_cmd_time * 4N`
2. Convert sizes: `demes_size = discoal_relative_size * N`
3. Convert migration rates: `demes_rate = discoal_rate / (4N)`

### Example Conversion

Discoal command:
```bash
discoal 100 1000 50000 -t 20 -r 15 -p 2 50 50 -ej 0.05 1 0 -en 0.1 0 2.0
```

Equivalent demes (with N=10000):
```yaml
description: Two population model with merger
time_units: generations
demes:
  - name: pop0
    epochs:
      - end_time: 4000      # 0.1 * 4 * 10000
        start_size: 20000   # 2.0 * 10000
      - end_time: 0
        start_size: 10000
  - name: pop1
    ancestors: [pop0]       # merges into pop0
    start_time: 2000        # 0.05 * 4 * 10000
    epochs:
      - end_time: 0
        start_size: 10000
```

With the same seed, this demes file (loaded through a `-Y` config with
`demography.demes_filename`) gives output identical to the command line
above.

## Further Resources

- [Demes specification](https://popsim-consortium.github.io/demes-spec/)
- [Demes tutorial](https://popsim-consortium.github.io/demes-docs/latest/tutorial.html)
- [stdpopsim catalog](https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html) - Pre-built demographic models