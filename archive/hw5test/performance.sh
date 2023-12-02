#!/bin/bash

# Check for correct number of arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: ./performance_test.sh <filename>"
    exit 1
fi

FILENAME=$1

# Function to run the program and extract the execution time
run_test() {
    NUM_THREADS=$1
    OUTPUT=$(./bellmanford "$FILENAME" "$NUM_THREADS")
    TIME=$(echo "$OUTPUT" | tail -n 1)
    echo "Threads: $NUM_THREADS, Time: $TIME"
}

# Test with 1, 2, and 4 threads
echo "Performance Test for $FILENAME"
run_test 1
run_test 2
run_test 4

