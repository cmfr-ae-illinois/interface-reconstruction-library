#!/bin/bash

# Parameters
NX=16
METHOD="Jibben"
SAMPLE_NX=32
RADIUS=2.5
SHAPE="sphere"

# Executable
EXEC="./compiled_examples/level_set_reconstruction/level_set_reconstruction"

# Weighting functions
WEIGHTS=("Wu2" "Wu4" "Wendland2" "Wendland4" "Wendland6")

# Run each weighting function
for WEIGHT in "${WEIGHTS[@]}"; do
    OUTPUT="level_set_viz/${WEIGHT}"

    mkdir -p "$OUTPUT"

    echo "Running ${WEIGHT}..."

    "$EXEC" \
        "$NX" \
        "$METHOD" \
        "$OUTPUT" \
        "$SAMPLE_NX" \
        "$RADIUS" \
        "$SHAPE" \
        "$WEIGHT"

    echo "Finished ${WEIGHT}"
    echo "--------------------------------"
done

echo "All simulations completed."