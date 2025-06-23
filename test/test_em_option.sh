#!/bin/bash
# Test the new -em option

echo "Testing -em option for migration rate changes"
echo

# Test 1: Basic migration rate change
echo "Test 1: Basic migration rate change"
echo "Command: ./build/discoal 10 1 1000 -t 10 -p 2 5 5 -em 0.1 0 1 0.001 -em 0.2 0 1 0"
./build/discoal 10 1 1000 -t 10 -p 2 5 5 -em 0.1 0 1 0.001 -em 0.2 0 1 0 > /dev/null
if [ $? -eq 0 ]; then
    echo "✓ Test 1 passed"
else
    echo "✗ Test 1 failed"
fi
echo

# Test 2: Multiple population migration changes
echo "Test 2: Multiple population migration changes"
echo "Command: ./build/discoal 12 1 1000 -t 10 -p 3 4 4 4 -em 0.05 0 1 0.002 -em 0.1 1 2 0.001 -em 0.15 0 1 0"
./build/discoal 12 1 1000 -t 10 -p 3 4 4 4 -em 0.05 0 1 0.002 -em 0.1 1 2 0.001 -em 0.15 0 1 0 > /dev/null
if [ $? -eq 0 ]; then
    echo "✓ Test 2 passed"
else
    echo "✗ Test 2 failed"
fi
echo

# Test 3: Migration with other events
echo "Test 3: Migration with other demographic events"
echo "Command: ./build/discoal 10 1 1000 -t 10 -p 2 5 5 -en 0.1 0 0.5 -em 0.15 0 1 0.001 -ed 0.2 1 0"
./build/discoal 10 1 1000 -t 10 -p 2 5 5 -en 0.1 0 0.5 -em 0.15 0 1 0.001 -ed 0.2 1 0 > /dev/null
if [ $? -eq 0 ]; then
    echo "✓ Test 3 passed"
else
    echo "✗ Test 3 failed"
fi

echo
echo "All tests completed"