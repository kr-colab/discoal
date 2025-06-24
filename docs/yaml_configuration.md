# YAML Configuration Interface

Discoal supports YAML configuration files as an alternative to command line arguments. This allows for more readable and maintainable simulation specifications, especially for complex demographic scenarios.

## Usage

To use a YAML configuration file, use the `-Y` flag:

```bash
discoal -Y config.yaml
```

Command line arguments can still be used alongside YAML files. Command line arguments take precedence over YAML values.

## YAML Configuration Structure

The YAML configuration file is organized into sections:

### 1. Simulation Section

Basic simulation parameters:

```yaml
simulation:
  sample_size: 20              # Total number of samples
  num_replicates: 100          # Number of independent replicates
  num_sites: 10000            # Length of simulated region
  seed: [12345, 67890]        # Random seeds (optional)
  effective_popn_size: 1000000 # Reference effective population size (default: 1000000)
```

### 2. Genetics Section

Mutation, recombination, and gene conversion parameters:

```yaml
genetics:
  mutation_rate: 10.0          # Population-scaled mutation rate (theta)
  recombination_rate: 8.0      # Population-scaled recombination rate (rho)
  
  # Gene conversion (optional)
  gene_conversion:
    rate: 2.0                  # Gene conversion initiation rate
    tract_length: 500          # Mean tract length
    
    # OR use crossover ratio mode (don't specify rate if using this)
    crossover_ratio: 0.5       # Ratio of gene conversion to crossover
    tract_length: 500          # Mean tract length
```

**Note**: When using `crossover_ratio`, the gene conversion rate is calculated as `rho * crossover_ratio`. Do not specify both `rate` and `crossover_ratio`.

### 3. Populations Section

Multi-population configuration:

```yaml
populations:
  count: 2                     # Number of populations
  sample_sizes: [10, 10]       # Samples per population
```

### 4. Events Section

Demographic and selection events:

#### Demographic Events

```yaml
events:
  demographic:
    - type: size_change
      time: 0.1                # Time in coalescent units (4N generations)
      population: 0            # Population ID
      size: 2.0                # New size relative to N
      
    - type: population_split
      time: 0.2
      source_pop: 1            # Source population (disappears)
      dest_pop: 0              # Destination population
      
    - type: migration_change
      time: 0.05
      source_pop: 0
      dest_pop: 1
      rate: 0.5                # Migration rate (4Nm)
```

#### Selection Events

```yaml
events:
  selection:
    - type: sweep
      mode: stochastic         # Options: stochastic, deterministic, neutral
      time: 0.01               # Time when allele reaches frequency 1
      selection_coeff: 100.0   # Population-scaled selection coefficient (2Ns)
      position: 0.5            # Position along sequence (0-1)
      initial_freq: 0.1        # Initial frequency (optional, for soft sweeps)
```

#### Demes File Integration

Instead of specifying demographic events directly, you can load a demes-format demographic model:

```yaml
events:
  demes_file: model.yaml       # Path to demes YAML file
```

**Note**: When using a demes file, demographic events in the YAML config are ignored.

## Complete Example

```yaml
# Simulation of a population with recent growth and a selective sweep
simulation:
  sample_size: 50
  num_replicates: 1000
  num_sites: 50000
  seed: [42, 12345]

genetics:
  mutation_rate: 20.0
  recombination_rate: 15.0
  gene_conversion:
    tract_length: 300
    crossover_ratio: 0.3

populations:
  count: 1
  sample_sizes: [50]

events:
  demographic:
    - type: size_change
      time: 0.01
      population: 0
      size: 0.1              # Bottleneck to 10% of original size
      
    - type: size_change
      time: 0.001
      population: 0
      size: 10.0             # Recent expansion to 10x original size
      
  selection:
    - type: sweep
      mode: stochastic
      time: 0.005
      selection_coeff: 500.0
      position: 0.7
```

## Parameter Conversion Notes

### Time Units
- YAML times are in coalescent units (4N generations)
- To convert from generations: `yaml_time = generations / (4N)`
- Example: 100,000 generations with N=10,000 → yaml_time = 2.5

### Population Sizes
- Sizes are relative to the reference effective population size N
- Example: `size: 2.0` means 2N individuals

### Migration Rates
- Migration rates are in units of 4Nm
- This matches the standard coalescent parameterization

### Selection Coefficients
- Selection coefficients are population-scaled: 2Ns
- s is the per-generation selection coefficient

## Validation

To ensure your YAML configuration produces expected results:

1. Use the same random seeds in YAML and command line versions
2. Run both versions and compare outputs
3. Use the validation test suite: `testing/yaml_validation_suite.sh`

## Limitations

Currently, the YAML interface does not support:
- Prior distributions for ABC sampling
- Some specialized sweep modes (partial sweeps with specific end frequencies)
- Complex migration matrices (use demes format for these)

For these features, use command line arguments or demes format files.