# YAML Configuration Examples

This directory contains example YAML configuration files for discoal.

## Basic Examples

- `basic_example.yaml` - Simple simulation with mutation and recombination
- `sweep_example.yaml` - Selective sweep simulation
- `demographic_example.yaml` - Population size changes and gene conversion
- `demes_example.yaml` - Using an external demes file for demography

## Running Examples

To run any example:

```bash
# From the discoal root directory
./discoal -Y config_examples/basic_example.yaml

# With additional command line options
./discoal -Y config_examples/sweep_example.yaml -T

# Using demes with population specification
./discoal -Y config_examples/demes_example.yaml -p 2 50 50
```

## Creating Your Own Configuration

1. Copy one of these examples as a starting point
2. Modify the parameters for your needs
3. See [docs/yaml_configuration.md](../docs/yaml_configuration.md) for full documentation

## Tips

- Use seeds for reproducibility
- Start with small sample sizes and sequence lengths for testing
- Validate your configuration by comparing with equivalent command line arguments
- Comment your YAML files to document parameter choices