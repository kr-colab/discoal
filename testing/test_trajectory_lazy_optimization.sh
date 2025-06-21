#!/bin/bash

# Test script to validate that lazy trajectory generation produces identical output
# Compares legacy vs trajectory-optimized versions

TEST_DIR="trajectory_lazy_test_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "=== Lazy Trajectory Generation Validation Test ==="
echo "Test directory: $TEST_DIR"
echo ""

# Test scenarios covering different trajectory usage patterns
declare -a TRAJECTORY_TEST_CASES=(
    # Format: "name:command_args"
    "basic_sweep:5 1 1000 -t 3 -r 2 -ws 0.3 -a 1000 -x 0.5"
    "stoch_sweep:10 1 2000 -t 8 -r 5 -ws 0.5 -a 500 -x 0.3"
    "neutral_sweep:8 1 1500 -t 6 -r 4 -ws 0.6 -a 0 -x 0.4"
    "partial_sweep:12 1 2500 -t 10 -r 8 -p 2 6 6 -ws 0.4 -a 800 -x 0.5 -f 0.8"
    "complex_demog:15 1 3000 -t 12 -r 10 -ws 0.2 -a 1200 -x 0.6 -en 0.1 0 0.5 -en 0.2 0 2.0"
    "variable_popsize:10 1 2000 -t 8 -r 6 -ws 0.35 -a 600 -x 0.45 -en 0.05 0 0.3 -en 0.15 0 3.0"
    "high_alpha:8 1 1500 -t 6 -r 4 -ws 0.4 -a 2000 -x 0.5"
    "low_alpha:12 1 2000 -t 9 -r 7 -ws 0.6 -a 100 -x 0.3"
    "edge_sweep_left:6 1 1200 -t 5 -r 3 -ws 0.05 -a 800 -x 0.4"
    "edge_sweep_right:6 1 1200 -t 5 -r 3 -ws 0.95 -a 800 -x 0.4"
)

TOTAL_TESTS=${#TRAJECTORY_TEST_CASES[@]}
PASSED_TESTS=0
FAILED_TESTS=0

for test_case in "${TRAJECTORY_TEST_CASES[@]}"; do
    # Parse test case
    test_name=$(echo $test_case | cut -d: -f1)
    test_args=$(echo $test_case | cut -d: -f2)
    
    echo "Testing: $test_name"
    echo "Args: $test_args"
    
    # Fixed seed for reproducible results
    SEED1=42424
    SEED2=84848
    
    # Run legacy version
    echo "  Running legacy version..."
    ./discoal_legacy $test_args -d $SEED1 $SEED2 > "$TEST_DIR/${test_name}_legacy.out" 2>&1
    legacy_exit=$?
    
    # Run trajectory optimized version  
    echo "  Running trajectory optimized version..."
    ../discoal_trajectory_optimized $test_args -d $SEED1 $SEED2 > "$TEST_DIR/${test_name}_optimized.out" 2>&1
    optimized_exit=$?
    
    # Check exit codes
    if [ $legacy_exit -ne $optimized_exit ]; then
        echo "  ❌ FAIL: Different exit codes (legacy: $legacy_exit, optimized: $optimized_exit)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        continue
    fi
    
    # Compare outputs (skip lines containing executable names)
    sed 's/discoal_legacy/discoal/g; ../discoal_trajectory_optimized/discoal/g' "$TEST_DIR/${test_name}_legacy.out" > "$TEST_DIR/${test_name}_legacy_filtered.out"
    sed 's/discoal_legacy/discoal/g; ../discoal_trajectory_optimized/discoal/g' "$TEST_DIR/${test_name}_optimized.out" > "$TEST_DIR/${test_name}_optimized_filtered.out"
    
    if diff -q "$TEST_DIR/${test_name}_legacy_filtered.out" "$TEST_DIR/${test_name}_optimized_filtered.out" > /dev/null; then
        echo "  ✅ PASS: Outputs identical (ignoring executable names)"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        # Clean up identical outputs to save space
        rm "$TEST_DIR/${test_name}_legacy.out" "$TEST_DIR/${test_name}_optimized.out"
        rm "$TEST_DIR/${test_name}_legacy_filtered.out" "$TEST_DIR/${test_name}_optimized_filtered.out"
    else
        echo "  ❌ FAIL: Outputs differ"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        # Keep failed outputs for inspection
        echo "    Legacy output saved to: $TEST_DIR/${test_name}_legacy.out"
        echo "    Optimized output saved to: $TEST_DIR/${test_name}_optimized.out"
        
        # Show first few lines of diff for quick inspection
        echo "    First 10 lines of content diff:"
        diff "$TEST_DIR/${test_name}_legacy_filtered.out" "$TEST_DIR/${test_name}_optimized_filtered.out" | head -10 | sed 's/^/    /'
        rm "$TEST_DIR/${test_name}_legacy_filtered.out" "$TEST_DIR/${test_name}_optimized_filtered.out"
    fi
    
    echo ""
done

echo "=== TRAJECTORY LAZY OPTIMIZATION VALIDATION SUMMARY ==="
echo "Total tests: $TOTAL_TESTS"
echo "Passed: $PASSED_TESTS"
echo "Failed: $FAILED_TESTS"

if [ $FAILED_TESTS -eq 0 ]; then
    echo "🎉 ALL TESTS PASSED - Lazy trajectory generation maintains identical functionality!"
    # Clean up test directory if all tests passed
    rmdir "$TEST_DIR" 2>/dev/null || echo "Test directory contains failure artifacts: $TEST_DIR"
else
    echo "⚠️  SOME TESTS FAILED - Review outputs in $TEST_DIR"
    exit 1
fi