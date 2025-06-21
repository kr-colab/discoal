#!/bin/bash
# Test suite to verify xoshiro256++ compatibility with legacy RNG

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "============================================"
echo "xoshiro256++ RNG Compatibility Test Suite"
echo "============================================"

# Function to compare outputs
compare_outputs() {
    local test_name=$1
    local file1=$2
    local file2=$3
    
    if diff -q "$file1" "$file2" > /dev/null; then
        echo -e "${GREEN}✓ $test_name: PASS${NC}"
        return 0
    else
        echo -e "${RED}✗ $test_name: FAIL${NC}"
        echo "  Outputs differ between legacy and xoshiro256++ versions"
        # Show first few differences
        echo "  First 5 differences:"
        diff "$file1" "$file2" | head -5
        return 1
    fi
}

# Build both versions
echo "Building discoal versions..."
make clean > /dev/null 2>&1
make discoal > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to build legacy discoal${NC}"
    exit 1
fi

if [ ! -f discoal_xoshiro ]; then
    make discoal_xoshiro > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to build discoal_xoshiro${NC}"
        exit 1
    fi
fi

echo -e "${GREEN}Successfully built both versions${NC}"
echo

# Test 1: Basic simulation with fixed seeds
echo "Test 1: Basic simulation with fixed seeds"
./discoal 10 1 1000 -d 12345 67890 > legacy_basic.out 2>/dev/null
./discoal_xoshiro 10 1 1000 -d 12345 67890 > xoshiro_basic.out 2>/dev/null
compare_outputs "Basic simulation" legacy_basic.out xoshiro_basic.out

# Test 2: Simulation with recombination
echo -e "\nTest 2: Simulation with recombination"
./discoal 20 1 10000 -r 100 -d 98765 43210 > legacy_recomb.out 2>/dev/null
./discoal_xoshiro 20 1 10000 -r 100 -d 98765 43210 > xoshiro_recomb.out 2>/dev/null
compare_outputs "Recombination" legacy_recomb.out xoshiro_recomb.out

# Test 3: Simulation with demographic events
echo -e "\nTest 3: Simulation with demographic events"
./discoal 15 1 5000 -en 0.1 0 0.5 -en 0.2 0 2.0 -d 11111 22222 > legacy_demog.out 2>/dev/null
./discoal_xoshiro 15 1 5000 -en 0.1 0 0.5 -en 0.2 0 2.0 -d 11111 22222 > xoshiro_demog.out 2>/dev/null
compare_outputs "Demographic events" legacy_demog.out xoshiro_demog.out

# Test 4: Multiple replicates
echo -e "\nTest 4: Multiple replicates"
./discoal 10 5 1000 -d 55555 66666 > legacy_multi.out 2>/dev/null
./discoal_xoshiro 10 5 1000 -d 55555 66666 > xoshiro_multi.out 2>/dev/null
compare_outputs "Multiple replicates" legacy_multi.out xoshiro_multi.out

# Test 5: Edge case seeds
echo -e "\nTest 5: Edge case seeds"
./discoal 5 1 500 -d 1 1 > legacy_edge1.out 2>/dev/null
./discoal_xoshiro 5 1 500 -d 1 1 > xoshiro_edge1.out 2>/dev/null
compare_outputs "Minimum seeds" legacy_edge1.out xoshiro_edge1.out

./discoal 5 1 500 -d 2147483562 2147483398 > legacy_edge2.out 2>/dev/null
./discoal_xoshiro 5 1 500 -d 2147483562 2147483398 > xoshiro_edge2.out 2>/dev/null
compare_outputs "Maximum seeds" legacy_edge2.out xoshiro_edge2.out

# Performance comparison
echo -e "\n============================================"
echo "Performance Comparison"
echo "============================================"

# Small simulation
echo -e "\nSmall simulation (10 samples, 1000 sites, 100 replicates):"
echo -n "Legacy: "
time -p ./discoal 10 100 1000 -d 12345 67890 > /dev/null 2>&1
echo -n "xoshiro256++: "
time -p ./discoal_xoshiro 10 100 1000 -d 12345 67890 > /dev/null 2>&1

# Medium simulation with recombination
echo -e "\nMedium simulation with recombination (20 samples, 10000 sites, 100 replicates, rho=100):"
echo -n "Legacy: "
time -p ./discoal 20 100 10000 -r 100 -d 12345 67890 > /dev/null 2>&1
echo -n "xoshiro256++: "
time -p ./discoal_xoshiro 20 100 10000 -r 100 -d 12345 67890 > /dev/null 2>&1

# Large simulation with high recombination
echo -e "\nLarge simulation with high recombination (20 samples, 100000 sites, 10 replicates, rho=1000):"
echo -n "Legacy: "
time -p ./discoal 20 10 100000 -r 1000 -d 12345 67890 > /dev/null 2>&1
echo -n "xoshiro256++: "
time -p ./discoal_xoshiro 20 10 100000 -r 1000 -d 12345 67890 > /dev/null 2>&1

# Clean up
rm -f legacy_*.out xoshiro_*.out

echo -e "\n============================================"
echo "Test suite complete"
echo "============================================"