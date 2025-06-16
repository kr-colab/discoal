# discoal Examples

This directory contains example scripts demonstrating how to use discoal with tskit output for various analyses.

## Prerequisites

Before running these examples, ensure you have:

1. Built discoal with tskit support:
   ```bash
   cd ..
   make clean
   make discoal
   ```

2. Installed required Python packages:
   ```bash
   pip install tskit numpy matplotlib
   ```

## Examples

### 1. Sweep Diversity Analysis (`sweep_diversity_analysis.py`)

This script demonstrates a complete workflow for analyzing diversity patterns around selective sweeps:

- Runs multiple replicate simulations with selective sweeps
- Outputs tree sequences in tskit format
- Calculates windowed nucleotide diversity
- Plots the average diversity pattern showing the characteristic valley around the sweep

**Usage:**
```bash
python sweep_diversity_analysis.py --replicates 10 --samples 20 --length 5000000
```

**Options:**
- `--replicates`: Number of simulation replicates (default: 10)
- `--samples`: Number of haploid samples (default: 20)
- `--length`: Sequence length in bp (default: 5,000,000)
- `--window`: Window size for diversity calculation (default: 5,000)
- `--cores`: Number of parallel cores to use (default: 4)
- `--discoal`: Path to discoal executable (default: ../discoal)

### 2. Simple Sweep Diversity (`sweep_diversity_simple.py`)

A simplified version for quick demonstration with fewer replicates and shorter sequences.

**Usage:**
```bash
python sweep_diversity_simple.py
```

This will:
- Run 10 replicates with 20 samples
- Simulate 5 Mb sequences
- Calculate diversity in 5 kb windows
- Save a plot showing the diversity valley

## Understanding the Output

The plots show:
- **Blue line**: Mean nucleotide diversity across replicates
- **Red dashed line**: Expected neutral diversity (θ/(1+θ))
- **Green dotted line**: Position of the selective sweep
- **Blue shaded area**: 95% confidence interval (in full analysis)

You should see a characteristic valley of reduced diversity centered on the sweep position, with the width and depth depending on:
- Selection strength (2Ns)
- Time since sweep (τ)
- Recombination rate (ρ)

## Parameter Settings

The examples use moderate to strong selection:
- θ = 40 (moderate diversity)
- ρ = 1000 (high recombination)
- 2Ns = 200 (strong selection)
- τ = 0.00 (very recent/ongoing sweeps)

These parameters create clear diversity valleys that are easy to visualize.

## Extending the Examples

You can modify these scripts to:
- Test different selection strengths
- Vary sweep ages
- Analyze other statistics (e.g., Tajima's D, haplotype homozygosity)
- Compare multiple sweep scenarios
- Add demographic models

## Troubleshooting

If you get errors:
1. Check that discoal is built: `ls ../discoal`
2. Ensure tskit is installed: `python -c "import tskit; print(tskit.__version__)"`
3. For parallel processing issues, reduce cores: `--cores 1`

## Citation

If you use these examples in your research, please cite the discoal paper and acknowledge the tskit integration.