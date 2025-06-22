#!/bin/bash

echo "=== Simple Growth Comparison Test ==="
echo ""

# Test 1: No growth baseline
echo "Test 1: No growth (should match)"
echo -n "Discoal: "
./discoal 10 100 1000 -t 10 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'
echo -n "MS:      "
./extern/ms 10 100 -t 10 -seed 12345 67890 0 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'

echo ""
echo "Test 2: Growth alpha=1.0"
echo -n "Discoal: "
./discoal 10 100 1000 -t 10 -eg 0.0 0 1.0 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'
echo -n "MS:      "
./extern/ms 10 100 -t 10 -eg 0.0 1 1.0 -seed 12345 67890 0 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'

echo ""
echo "Test 3: Growth alpha=2.0"
echo -n "Discoal: "
./discoal 10 100 1000 -t 10 -eg 0.0 0 2.0 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'
echo -n "MS:      "
./extern/ms 10 100 -t 10 -eg 0.0 1 2.0 -seed 12345 67890 0 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'

echo ""
echo "Test 4: Growth alpha=0.5"
echo -n "Discoal: "
./discoal 10 100 1000 -t 10 -eg 0.0 0 0.5 -d 12345 67890 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'
echo -n "MS:      "
./extern/ms 10 100 -t 10 -eg 0.0 1 0.5 -seed 12345 67890 0 2>/dev/null | grep "segsites:" | awk '{sum+=$2; count++} END {printf "%.2f\n", sum/count}'