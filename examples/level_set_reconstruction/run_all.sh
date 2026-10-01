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

# Master output
MASTER_OUTPUT="level_set_viz/errors_master.txt"

# Remove old master file
rm -f "$MASTER_OUTPUT"

FIRST_RUN=true

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

    ERROR_FILE="${OUTPUT}/errors.txt"

    if [ "$FIRST_RUN" = true ]; then
        cat "$ERROR_FILE" > "$MASTER_OUTPUT"
        FIRST_RUN=false
    else
        tail -n 1 "$ERROR_FILE" >> "$MASTER_OUTPUT"
    fi

    echo "Finished ${WEIGHT}"
    echo "--------------------------------"
done

echo "All simulations completed."
echo "Master errors written to ${MASTER_OUTPUT}"