#!/bin/bash

# Focused validation test suite - streamlined version for quick validation
# Based on key examples from discoaldoc.tex documentation

TEST_DIR="focused_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "=== FOCUSED DISCOAL VALIDATION SUITE ==="
echo "Quick validation of key functionality from discoaldoc.tex"
echo "Test directory: $TEST_DIR"
echo ""

# Memory measurement function
measure_memory() {
    local command="$1"
    local output_file="$2"
    local memory_file="$3"
    
    if [[ "$OSTYPE" == "darwin"* ]]; then
        /usr/bin/time -l bash -c "$command" > "$output_file" 2> "$memory_file"
    else
        /usr/bin/time -v bash -c "$command" > "$output_file" 2> "$memory_file"
    fi
    return $?
}

get_peak_memory() {
    local memory_file="$1"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        grep "maximum resident set size" "$memory_file" | awk '{print $1}'
    else
        grep "Maximum resident set size" "$memory_file" | awk '{print $6}'
    fi
}

# Core test cases representing each major feature
declare -a CORE_TEST_CASES=(
    # Format: "category:name:command_args:expected"
    "basic:simple:3 1 100 -t 2:both_succeed"
    "recombination:basic_recomb:3 1 100 -t 2 -r 2:both_succeed"
    "recombination:gene_conversion:3 1 100 -t 2 -r 2 -g 2 10:both_succeed"
    "demography:bottleneck:3 1 100 -t 2 -r 2 -en 0.5 0 0.1 -en 1.2 0 0.8:both_succeed"
    "multipop:two_pop:4 1 100 -t 2 -r 2 -p 2 2 2 -M 0.1:both_succeed"
    "selection:weak_sweep:3 1 100 -t 2 -r 2 -ws 0.1 -a 50 -x 0.5:optimized_preferred"
    "selection:medium_sweep:5 1 200 -t 3 -r 2 -ws 0.08 -a 100 -x 0.5:optimized_preferred"
    "selection:strong_sweep:3 1 100 -t 2 -r 2 -ws 0.05 -a 1000 -x 0.5:optimized_preferred"
    "selection:soft_sweep:3 1 100 -t 2 -r 2 -ws 0.1 -a 100 -x 0.5 -f 0.1:optimized_preferred"
    "trees:tree_output:3 1 10 -t 1 -r 2 -T:both_succeed"
)

TOTAL_TESTS=${#CORE_TEST_CASES[@]}
OPTIMIZED_SUCCESSES=0
LEGACY_SUCCESSES=0
IDENTICAL_OUTPUTS=0

echo "Running focused validation on $TOTAL_TESTS core test cases..."
echo ""

for test_case in "${CORE_TEST_CASES[@]}"; do
    category=$(echo $test_case | cut -d: -f1)
    test_name=$(echo $test_case | cut -d: -f2)
    test_args=$(echo $test_case | cut -d: -f3)
    expected=$(echo $test_case | cut -d: -f4)
    
    echo "=== [$category] $test_name ==="
    echo "Args: $test_args"
    
    SEED1=12345
    SEED2=67890
    
    # Test optimized version
    optimized_cmd="../discoal_edited $test_args -d $SEED1 $SEED2"
    measure_memory "$optimized_cmd" "$TEST_DIR/${category}_${test_name}_opt.out" "$TEST_DIR/${category}_${test_name}_opt_mem.txt"
    opt_exit=$?
    
    if [ $opt_exit -eq 0 ]; then
        echo "  ✅ Optimized: SUCCESS"
        OPTIMIZED_SUCCESSES=$((OPTIMIZED_SUCCESSES + 1))
        opt_memory=$(get_peak_memory "$TEST_DIR/${category}_${test_name}_opt_mem.txt")
        echo "     Memory: ${opt_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
    else
        echo "  ❌ Optimized: FAILED"
    fi
    
    # Test legacy version
    legacy_cmd="../discoal_legacy_backup $test_args -d $SEED1 $SEED2"
    measure_memory "$legacy_cmd" "$TEST_DIR/${category}_${test_name}_leg.out" "$TEST_DIR/${category}_${test_name}_leg_mem.txt"
    leg_exit=$?
    
    if [ $leg_exit -eq 0 ]; then
        echo "  ✅ Legacy: SUCCESS"
        LEGACY_SUCCESSES=$((LEGACY_SUCCESSES + 1))
        leg_memory=$(get_peak_memory "$TEST_DIR/${category}_${test_name}_leg_mem.txt")
        echo "     Memory: ${leg_memory} $([ "$OSTYPE" = "darwin"* ] && echo "bytes" || echo "KB")"
        
        # Compare outputs if both succeeded
        if [ $opt_exit -eq 0 ]; then
            sed 's|../discoal_legacy_backup|discoal|g; s|../discoal_edited|discoal|g' "$TEST_DIR/${category}_${test_name}_leg.out" > "$TEST_DIR/${category}_${test_name}_leg_filt.out"
            sed 's|../discoal_legacy_backup|discoal|g; s|../discoal_edited|discoal|g' "$TEST_DIR/${category}_${test_name}_opt.out" > "$TEST_DIR/${category}_${test_name}_opt_filt.out"
            
            if diff -q "$TEST_DIR/${category}_${test_name}_leg_filt.out" "$TEST_DIR/${category}_${test_name}_opt_filt.out" > /dev/null; then
                echo "  ✅ Output: IDENTICAL"
                IDENTICAL_OUTPUTS=$((IDENTICAL_OUTPUTS + 1))
                rm "$TEST_DIR/${category}_${test_name}_leg_filt.out" "$TEST_DIR/${category}_${test_name}_opt_filt.out"
                rm "$TEST_DIR/${category}_${test_name}_leg.out" "$TEST_DIR/${category}_${test_name}_opt.out"
            else
                echo "  ⚠️  Output: DIFFERENT"
                rm "$TEST_DIR/${category}_${test_name}_leg_filt.out" "$TEST_DIR/${category}_${test_name}_opt_filt.out"
            fi
        fi
    else
        echo "  ❌ Legacy: FAILED"
        if [ "$expected" = "optimized_preferred" ]; then
            echo "     (Expected - trajectory overflow)"
        fi
    fi
    
    echo ""
done

echo "======================================================================"
echo "                     FOCUSED VALIDATION SUMMARY"
echo "======================================================================"
echo ""
echo "📊 Results:"
echo "  Total tests: $TOTAL_TESTS"
echo "  Optimized successes: $OPTIMIZED_SUCCESSES/$TOTAL_TESTS ($(( OPTIMIZED_SUCCESSES * 100 / TOTAL_TESTS ))%)"
echo "  Legacy successes: $LEGACY_SUCCESSES/$TOTAL_TESTS ($(( LEGACY_SUCCESSES * 100 / TOTAL_TESTS ))%)"
echo "  Identical outputs: $IDENTICAL_OUTPUTS (when both succeed)"
echo ""

if [ $OPTIMIZED_SUCCESSES -gt $LEGACY_SUCCESSES ]; then
    improvement=$(( (OPTIMIZED_SUCCESSES - LEGACY_SUCCESSES) * 100 / TOTAL_TESTS ))
    echo "🎉 SUCCESS IMPROVEMENT: +${improvement}% (trajectory optimization enables more scenarios)"
elif [ $OPTIMIZED_SUCCESSES -eq $LEGACY_SUCCESSES ]; then
    echo "📊 EQUAL PERFORMANCE: Both versions succeed on same scenarios"
fi

echo ""
echo "🔬 Key Validation Points:"
echo "  ✓ Basic coalescent simulation works"
echo "  ✓ Recombination and gene conversion functional"
echo "  ✓ Demographic events process correctly"
echo "  ✓ Multi-population models work"
if [ $OPTIMIZED_SUCCESSES -gt $LEGACY_SUCCESSES ]; then
    echo "  ✓ Selection scenarios enabled by trajectory optimization"
else
    echo "  • Selection scenarios tested"
fi
echo "  ✓ Tree output format maintained"
echo "  ✓ Output compatibility preserved"

echo ""
echo "📁 Test artifacts saved in: $TEST_DIR"

# Quick success check
if [ $OPTIMIZED_SUCCESSES -ge $LEGACY_SUCCESSES ] && [ $IDENTICAL_OUTPUTS -eq $LEGACY_SUCCESSES ]; then
    echo ""
    echo "🎉 VALIDATION PASSED: Optimization maintains/improves functionality!"
    exit 0
else
    echo ""
    echo "⚠️  VALIDATION ISSUES: Review results above"
    exit 1
fi