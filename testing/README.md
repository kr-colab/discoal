# discoal checks

This directory holds the checks of discoal's output. The unit tests of the
C code are in `test/unit` (`make run_tests`).

## Quick start

```bash
make check                    # the checks of every pull request
make check CHECK_WORKERS=16   # more parallel workers for the statistical checks
```

`make check` builds discoal, runs the unit tests, and then runs each check
below that CI runs for every pull request, in the same order. The
statistical checks need Python with msprime, tskit, numpy, scipy, demes and
stdpopsim. They use fixed seeds, so their results do not change from run to
run, and each one exits with 0 if every comparison passes.

## The checks

| check | what it checks | runs |
|---|---|---|
| `yaml_parser_suite.sh` | the YAML parser accepts the valid files in `debug_yaml/` and rejects the invalid ones with the right error | every PR |
| `yaml_validation_suite.sh` | each example in `config_examples/` gives the same output as the equivalent command line (ms output and tree sequences, compared with `cmp_trees.py`) | every PR |
| `preflight_validation_suite.sh` | many command lines run and give output, or stop with the expected error; regression cases of fixed defects | every PR |
| `demes_validation_suite.sh` | demes models load and run, or are rejected where they must be (through `demes_run.sh`) | every PR |
| `demes_output_comparison_suite.sh` | a demes model and the same model given on the command line give the same output | every PR |
| `units_convention_check.py` | time, size, growth and migration units against msprime | every PR |
| `sweep_snp_check.py` | the selected SNP of partial sweeps: its position and the number of its carriers | every PR |
| `trajectory_file_check.py` | sweeps from a trajectory file (`-S`) | every PR |
| `recurrent_sweep_check.py` | recurrent sweeps (`-R`, `-L`) follow the demographic model between sweeps | every PR |
| `sweep_msprime_check.py` | hard sweeps against msprime's `SweepGenicSelection` | every PR |
| `gc_crossover_check.py` | gene conversion with tracts longer than the locus equals crossover at twice its rate, in the neutral phase and in a sweep | every PR |
| `gene_conversion_check.py` | gene conversion against the two-locus theory (at the left end, in the middle and at the right end of the locus) and against msprime | every PR |
| `linkage_msprime_check.py` | linkage of two sites (same MRCA, TMRCA correlation) against msprime | every PR |
| `fuzz_msprime_check.py` | random command-line demographic models against msprime | every PR (150 models) |
| `stdpopsim_parity.py` | stdpopsim and demes models against msprime | every PR (a core set and a random sample); weekly (all models, without and with recombination) |
| `slim_sweep_check.py` | hard, soft and partial sweeps and ancient samples against forward simulations in SLiM | by hand (needs SLiM, pyslim, stdpopsim) |
| `yaml_parser_fuzz.py` | the YAML parser on many random variants of valid files: no crash, no silent acceptance of bad input | by hand |

Helpers: `cmp_trees.py` compares two tree-sequence files (provenance
ignored); `demes_run.sh` runs discoal on a demes file through a YAML
configuration; `debug_yaml/` holds the YAML parser's test files.

## Running one check

Each script runs from the repository root and uses `build/discoal`:

```bash
bash testing/preflight_validation_suite.sh
python3 testing/sweep_msprime_check.py --workers 8
python3 testing/fuzz_msprime_check.py --models 150 --reps 500 --workers 8
```

The Python checks accept `--discoal PATH` to test another build, and print
one line per comparison (`ok` or `FAIL`, with the z-score and the KS p-value)
and a summary line.

## The parity check against msprime

`stdpopsim_parity.py` needs the demes specification for its demes models:

```bash
git clone https://github.com/popsim-consortium/demes-spec
python3 testing/stdpopsim_parity.py --reps 1000 --generic 20 --seed 1 \
    --demes-dir demes-spec/examples --demes-dir demes-spec/test-cases/valid \
    --workers 16
```

Add `--rho 100 --sites 10000` to compare with recombination (on the number
of trees). Results are cached in `testing/.parity_cache` (only what changed
runs again; `--no-cache` turns it off). `--tsv` writes the per-test results,
`--jobs-tsv` the per-model status and run times. The weekly workflow
(`.github/workflows/parity-weekly.yml`) runs all models with and without
recombination.
