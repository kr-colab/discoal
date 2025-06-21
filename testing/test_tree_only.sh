#!/bin/bash

# Test tree-only version against edited version

echo "Testing tree-only version..."

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

TREE_ONLY="../discoal_tree_only"
EDITED="../discoal_edited"

# Test cases
declare -a tests=(
    "5 1 100 -t 10"
    "5 1 100 -t 10 -r 10"
    "10 1 1000 -t 5 -r 5"
    "5 1 100 -t 10 -r 10 -g 5 20"
    "10 1 1000 -t 10 -r 10 -ws 0.1 -a 100"
    "10 1 1000 -t 10 -r 10 -p 2 5 5 -ed 0.1 0 1"
)

passed=0
failed=0

for test in "${tests[@]}"; do
    echo -n "Testing: $test ... "
    
    # Run both versions with same seed
    $TREE_ONLY $test -d 123 456 > tree_only.out 2>&1
    $EDITED $test -d 123 456 > edited.out 2>&1
    
    # Compare outputs (ignoring first line with command)
    if tail -n +2 tree_only.out | diff -q - <(tail -n +2 edited.out) > /dev/null; then
        echo -e "${GREEN}PASS${NC}"
        ((passed++))
    else
        echo -e "${RED}FAIL${NC}"
        echo "  Differences found:"
        tail -n +2 tree_only.out | diff - <(tail -n +2 edited.out) | head -10
        ((failed++))
    fi
done

echo ""
echo "Summary: $passed passed, $failed failed"

# Clean up
rm -f tree_only.out edited.out

exit $failed