# YAML Configuration Examples

This directory contains example YAML configuration files for discoal.

## Basic Examples

- `basic_example.yaml` - Simple simulation with mutation and recombination
- `sweep_example.yaml` - A soft (-f), partial (-c) sweep, with the selected site hidden (-h)
- `hard_sweep_example.yaml` - A hard sweep from one new mutation (-ws)
- `deterministic_sweep_example.yaml` - The same sweep with a deterministic trajectory (-wd)
- `neutral_sweep_example.yaml` - A neutral fixation, a null model for sweeps (-wn)
- `recurrent_mutation_sweep_example.yaml` - A soft sweep from recurrent beneficial mutation (-uA)
- `recurrent_sweeps_example.yaml` - Recurrent sweeps in the locus (-R)
- `gene_conversion_example.yaml` - Gene conversion at a set rate (-g)
- `gene_conversion_ratio_example.yaml` - Gene conversion relative to crossing over (-gr)
- `demographic_example.yaml` - Three populations with migration, size changes and splits
- `demes_example.yaml` - Using an external demes file for demography
- `migration_change_example.yaml` - Migration rates that change over time
- `growth_example.yaml` - Exponential (-eg) and linear (-el) growth
- `ancient_samples_example.yaml` - Ancient samples in two populations
- `ancient_demes_example.yaml` - Ancient samples with a demes model, from an extinct deme
- `admixture_example.yaml` - An admixed population (-ea)
- `linked_sweep_example.yaml` - A sweep outside the locus (-ls) with a set trajectory time step (-i)
- `recurrent_linked_sweeps_example.yaml` - Recurrent sweeps to the left of the locus (-L)
- `priors_example.yaml` - Parameter values drawn for each replicate (-P flags)
- `trajectory_file_example.yaml` - A partial sweep with its trajectory read from `trajectory_example.traj` (-S)
- `tree_sequence_example.yaml` - Tree sequence output (-ts)

## Running Examples

To run any example:

```bash
# From the discoal root directory
./build/discoal -Y config_examples/basic_example.yaml

# Command line options after the file override its values
./build/discoal -Y config_examples/sweep_example.yaml -d 1 2
```

Each example starts with a line `# Equivalent to: discoal ...`, the
same simulation on the command line. `testing/yaml_validation_suite.sh`
runs both forms of every example and checks that the outputs are
identical.

## Creating Your Own Configuration

1. Copy one of these examples as a starting point
2. Modify the parameters for your needs
3. See [docs/yaml_configuration.md](../docs/yaml_configuration.md) for full documentation

## Tips

- Use seeds for reproducibility
- Start with small sample sizes and sequence lengths for testing
- Validate your configuration by comparing with equivalent command line arguments
- Comment your YAML files to document parameter choices
