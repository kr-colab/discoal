#!/bin/bash
# Statistical validation for xoshiro256++ RNG implementation

echo "============================================"
echo "xoshiro256++ Statistical Validation"
echo "============================================"

# Test that both RNGs produce valid simulations with expected properties
echo "Running statistical tests..."

# Test 1: Check that both produce valid ms-format output
echo -e "\nTest 1: Output format validation"
./discoal 10 10 1000 -t 10 -d 12345 67890 > legacy_test.ms 2>/dev/null
./discoal_xoshiro 10 10 1000 -t 10 -d 12345 67890 > xoshiro_test.ms 2>/dev/null

# Count segregating sites
legacy_segsites=$(grep "segsites:" legacy_test.ms | awk '{sum += $2} END {print sum}')
xoshiro_segsites=$(grep "segsites:" xoshiro_test.ms | awk '{sum += $2} END {print sum}')

echo "Legacy RNG: Total segregating sites = $legacy_segsites"
echo "xoshiro256++: Total segregating sites = $xoshiro_segsites"

# Test 2: Distribution of segregating sites over many replicates
echo -e "\nTest 2: Distribution test (1000 replicates)"
./discoal 10 1000 1000 -t 10 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{print $2}' | sort | uniq -c > legacy_dist.txt
./discoal_xoshiro 10 1000 1000 -t 10 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{print $2}' | sort | uniq -c > xoshiro_dist.txt

# Calculate mean segregating sites
legacy_mean=$(./discoal 10 1000 1000 -t 10 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum += $2; count++} END {print sum/count}')
xoshiro_mean=$(./discoal_xoshiro 10 1000 1000 -t 10 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum += $2; count++} END {print sum/count}')

echo "Legacy mean segregating sites: $legacy_mean"
echo "xoshiro256++ mean segregating sites: $xoshiro_mean"

# Test 3: Recombination breakpoints distribution
echo -e "\nTest 3: Recombination test"
./discoal 10 100 10000 -r 100 -d 12345 67890 > legacy_recomb.ms 2>/dev/null
./discoal_xoshiro 10 100 10000 -r 100 -d 12345 67890 > xoshiro_recomb.ms 2>/dev/null

# Count unique haplotypes as a proxy for recombination
legacy_haplos=$(grep -A 10 "segsites:" legacy_recomb.ms | grep "^[01]" | sort -u | wc -l)
xoshiro_haplos=$(grep -A 10 "segsites:" xoshiro_recomb.ms | grep "^[01]" | sort -u | wc -l)

echo "Legacy unique haplotypes: $legacy_haplos"
echo "xoshiro256++ unique haplotypes: $xoshiro_haplos"

# Test 4: RNG direct function tests
echo -e "\nTest 4: RNG function validation"
cat > test_rng_functions.c << 'EOF'
#include <stdio.h>
#include <math.h>
#include "ranlib.h"

int main() {
    // Initialize RNG
    setall(12345, 67890);
    
    // Test ranf() - should be uniform [0,1)
    double sum = 0.0;
    int n = 10000;
    for (int i = 0; i < n; i++) {
        double x = ranf();
        sum += x;
        if (x < 0.0 || x >= 1.0) {
            printf("ERROR: ranf() returned %f\n", x);
            return 1;
        }
    }
    printf("ranf() mean: %f (expected ~0.5)\n", sum/n);
    
    // Test ignuin() - should be uniform integers
    long sum_int = 0;
    for (int i = 0; i < n; i++) {
        long x = ignuin(1, 100);
        sum_int += x;
        if (x < 1 || x > 100) {
            printf("ERROR: ignuin(1,100) returned %ld\n", x);
            return 1;
        }
    }
    printf("ignuin(1,100) mean: %f (expected ~50.5)\n", (double)sum_int/n);
    
    // Test sexpo() - should be exponential with mean 1
    sum = 0.0;
    for (int i = 0; i < n; i++) {
        double x = sexpo();
        sum += x;
        if (x < 0.0) {
            printf("ERROR: sexpo() returned negative %f\n", x);
            return 1;
        }
    }
    printf("sexpo() mean: %f (expected ~1.0)\n", sum/n);
    
    // Test genexp() - exponential with different means
    sum = 0.0;
    double mean = 5.0;
    for (int i = 0; i < n; i++) {
        double x = genexp(mean);
        sum += x;
    }
    printf("genexp(5.0) mean: %f (expected ~5.0)\n", sum/n);
    
    return 0;
}
EOF

# Compile and run with both RNG implementations
echo -e "\nTesting legacy RNG functions:"
gcc -O2 -o test_rng_legacy test_rng_functions.c ranlibComplete.c -lm
./test_rng_legacy

echo -e "\nTesting xoshiro256++ RNG functions:"
gcc -O2 -DUSE_XOSHIRO256PP -o test_rng_xoshiro test_rng_functions.c xoshiro256pp_compat.c -lm
./test_rng_xoshiro

# Performance benchmark
echo -e "\n============================================"
echo "Performance Benchmark"
echo "============================================"

echo -e "\nBenchmarking high-recombination scenario (20 samples, 100000 sites, rho=1000):"
echo "Running 5 replicates each..."

# Legacy timing
echo -n "Legacy RNG times: "
for i in {1..5}; do
    /usr/bin/time -f "%e" ./discoal 20 1 100000 -r 1000 -d $i $i > /dev/null 2>&1
done | tr '\n' ' '
echo

# xoshiro timing  
echo -n "xoshiro256++ times: "
for i in {1..5}; do
    /usr/bin/time -f "%e" ./discoal_xoshiro 20 1 100000 -r 1000 -d $i $i > /dev/null 2>&1
done | tr '\n' ' '
echo

# Cleanup
rm -f legacy_*.ms xoshiro_*.ms legacy_*.txt xoshiro_*.txt
rm -f test_rng_functions.c test_rng_legacy test_rng_xoshiro

echo -e "\n============================================"
echo "Statistical validation complete"
echo "============================================"