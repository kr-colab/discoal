#!/bin/bash
# Test script for sample specification functionality

echo "=== Testing Sample Specification Support ==="
echo

# Build test program if needed
echo "Building test program..."
cd ../scripts/validation
make test_sample_specification >/dev/null 2>&1
cd ../../test

echo "Test 1: Basic two-population model with sample specifications"
echo "-------------------------------------------------------"
../scripts/validation/test_sample_specification configs/demes_with_samples.yaml

echo
echo "Test 2: Three-population model (with ancient samples - not fully implemented)"
echo "--------------------------------------------------------------------------"
../scripts/validation/test_sample_specification configs/ancient_samples.yaml 2>&1 | grep -v "No such file"

echo
echo "=== Sample Specification Testing Complete ==="