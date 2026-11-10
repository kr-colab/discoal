# Demes to Discoal Converter

This tool converts demes YAML demographic models to equivalent discoal command line arguments.

## Usage

```bash
python scripts/demes_to_discoal.py model.yaml [--N0 1000000] [--samples 20 20]
```

## Examples

### Two populations with asymmetric migration

```yaml
# asymmetric_migration.yaml
time_units: generations
generation_time: 1
demes:
- name: pop0
  epochs:
  - {end_time: 0, start_size: 1000000}
- name: pop1
  epochs:
  - {end_time: 0, start_size: 1000000}
migrations:
- {source: pop0, dest: pop1, rate: 2e-6}
- {source: pop1, dest: pop0, rate: 1e-6}
```

Convert to discoal:
```bash
python scripts/demes_to_discoal.py asymmetric_migration.yaml --samples 10 10
# Output: -p 2 10 10 -em 0.0 1 0 8.0 -em 0.0 0 1 4.0
```

### Population split with size change

```yaml
# pop_split.yaml
time_units: generations
generation_time: 1
demes:
- name: deme0
  epochs:
  - {end_time: 200000, start_size: 2000000}
  - {end_time: 0, start_size: 1000000}
- name: deme1
  ancestors: [deme0]
  start_time: 200000
  epochs:
  - {end_time: 0, start_size: 1000000}
```

Convert to discoal:
```bash
python scripts/demes_to_discoal.py pop_split.yaml --samples 20 20
# Output: -p 2 20 20 -en 0.05 0 2.0 -ej 0.05 1 0 -en 0.05 0 2.0
```

## Parameter Conversions

The tool follows the same conversion rules as the demes integration in discoal:

- **Time**: `generations / (4 * N0)` → discoal time units
- **Population size**: `demes_size / N0` → relative size
- **Migration rate**: `demes_rate * 4 * N0` → discoal 4Nm units
- **Reference N0**: Default is 1,000,000 (matching discoal's `EFFECTIVE_POPN_SIZE`)

## Limitations

- Admixture events (multiple ancestors) are not fully implemented yet
- Pulse migrations are not implemented yet
- Exponential growth is not supported (consistent with discoal limitations)

## Verification

You can verify the conversion by comparing outputs:

For a simple model (two demes, no extinct demes, the same present-day
size for both, and `--N0` set to the present-day size of the first
deme), the two runs below give identical output. Run them in bash, so
that `$DISCOAL_CMD` splits into words. For other models the script's
command line can differ from discoal's own import (for example in the
reference size, other present-day sizes, and splits from an extinct
deme), so compare summary statistics instead.

```bash
# Run with the demes file, through a -Y config (config.yaml) with
# sample_size 20, num_replicates 5, num_sites 10000, mutation_rate 10,
# recombination_rate 5, and
# demography: {deme_sample_size: [10, 10], demes_filename: model.yaml}
discoal -Y config.yaml -d 123 456

# Run with converted command line
DISCOAL_CMD=$(python scripts/demes_to_discoal.py model.yaml --samples 10 10 --N0 <N0>)
discoal 20 5 10000 -t 10 -r 5 $DISCOAL_CMD -d 123 456
```
