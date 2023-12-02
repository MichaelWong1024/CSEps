#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: ./compare_outputs.sh <filename> <num_threads> [destination_vertex]"
    exit 1
fi

FILENAME=$1
NUM_THREADS=$2
DESTINATION=${3:-''}  # Set to an empty string if not provided

# Run the C program and exclude the last line (execution time)
OUTPUT_C=$(./bellmanford "$FILENAME" "$NUM_THREADS" "$DESTINATION" | head -n -1)

# Run the Python program
OUTPUT_PY=$(python3 bellmanford.py "$FILENAME" "$NUM_THREADS" "$DESTINATION")

# Compare outputs (excluding execution time)
echo "Comparing outputs..."
echo "C Program Output:"
echo "$OUTPUT_C"
echo ""
echo "Python Program Output:"
echo "$OUTPUT_PY"
echo ""

if [ "$OUTPUT_C" == "$OUTPUT_PY" ]; then
    echo "Outputs match!"
else
    echo "Outputs differ!"
fi
