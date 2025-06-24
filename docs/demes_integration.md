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
- Population splits and mergers
- Symmetric and asymmetric migration
- Pulse migration events
- Sample size specification per population

### Not Supported
- Exponential growth epochs (discoal only supports constant size epochs)
- Linear growth epochs
- Selfing rates
- Cloning rates

If your demes file contains unsupported features, discoal will report an error and exit.

## Parameter Conversion

Discoal uses different time and rate units than the demes specification. The integration handles these conversions automatically:

### Time Conversion
- Demes: time in generations
- Discoal: time in units of 4N generations (coalescent time)
- Conversion: `discoal_time = demes_time / (4N)`

Where N is the reference effective population size (the size of the first present-day deme).

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

When using multi-population demes models, specify sample sizes after the `-p` flag:

```bash
# 50 samples from pop0, 50 from pop1
discoal 100 1000 50000 -D two_pop_model.yaml -p 2 50 50 -t 20
```

The population order matches the order of present-day demes in the file.

## Validation and Debugging

### Checking Conversion
Discoal reports the loaded demographic model:

```
Loaded 2 populations and 5 events from demes file 'model.yaml'
```

### Common Issues

1. **Infinite coalescent times**: Ensure all populations are connected through migration or common ancestors
2. **Time conversion confusion**: Remember that discoal command line times are in units of 2N, but internal times are 4N
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
Demes files written for msprime work directly with discoal (if they don't use growth epochs).

### From discoal command line
To convert a discoal command line to demes:
1. Convert times: `demes_time = discoal_cmd_time * 2N`
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
      - start_time: .inf
        end_time: 2000      # 0.1 * 2 * 10000
        size_function: constant
        size: 10000
      - start_time: 2000
        end_time: 0
        size_function: constant
        size: 20000         # 2.0 * 10000
  - name: pop1
    epochs:
      - start_time: .inf
        end_time: 1000      # 0.05 * 2 * 10000
        size_function: constant
        size: 10000
        ancestors: [pop0]   # Merge into pop0
```

## Further Resources

- [Demes specification](https://popsim-consortium.github.io/demes-spec/)
- [Demes tutorial](https://popsim-consortium.github.io/demes-docs/latest/tutorial.html)
- [stdpopsim catalog](https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html) - Pre-built demographic models