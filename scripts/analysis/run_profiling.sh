#!/bin/bash
#
# Run discoal vs mspms profiling comparison
# This script activates the correct conda environment and runs the profiling
#

set -e

echo "Activating discoal_dev conda environment..."
eval "$(conda shell.bash hook)"
conda activate discoal_dev

echo "Checking for mspms..."
if ! command -v mspms &> /dev/null; then
    echo "ERROR: mspms not found in discoal_dev environment"
    echo "Please ensure mspms is installed: conda install -c conda-forge msprime"
    exit 1
fi

echo "Running profiling comparison..."
python profile_discoal_vs_mspms.py "$@"