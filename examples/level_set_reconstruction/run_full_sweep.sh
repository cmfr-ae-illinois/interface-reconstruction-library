#!/bin/bash

# ============================================================
# Parameter sweeps
# ============================================================

NX_VALUES=(16 32 64 128)
METHODS=("LVIRA" "Jibben")
SHAPES=("ellipsoid")

# Radius: 1.0, 1.1, ..., 5.0
RADIUS_VALUES=$(seq 1.0 0.1 5.0)

# Fixed parameter
SAMPLE_NX=8

# Weighting functions
WEIGHTS=("Wu4" "Wendland2" "Wendland4" "Wendland6" "Wu2")
# For all:
# WEIGHTS=("Wu2" "Wu4" "Wendland2" "Wendland4" "Wendland6")

# Executable
EXEC="./compiled_examples/level_set_reconstruction/level_set_reconstruction"

# Master output
MASTER_OUTPUT="level_set_viz/errors_master_Ellipsoid.txt"

rm -f "$MASTER_OUTPUT"

FIRST_RUN=true
RUN=0

# ============================================================
# Parameter sweep
# ============================================================

for WEIGHT in "${WEIGHTS[@]}"; do
    for NX in "${NX_VALUES[@]}"; do
        for METHOD in "${METHODS[@]}"; do
            for SHAPE in "${SHAPES[@]}"; do
                for RADIUS in $RADIUS_VALUES; do

                    RUN=$((RUN + 1))

                    OUTPUT="level_set_viz/${WEIGHT}/NX${NX}/${METHOD}/${SHAPE}/R${RADIUS}"

                    mkdir -p "$OUTPUT"

                    echo "========================================"
                    echo "Run ${RUN}"
                    echo "Weight: ${WEIGHT}"
                    echo "NX:     ${NX}"
                    echo "Method: ${METHOD}"
                    echo "Shape:  ${SHAPE}"
                    echo "Radius: ${RADIUS}"
                    echo "========================================"

                    "$EXEC" \
                        "$NX" \
                        "$METHOD" \
                        "$OUTPUT" \
                        "$SAMPLE_NX" \
                        "$RADIUS" \
                        "$SHAPE" \
                        "$WEIGHT"

                    ERROR_FILE="${OUTPUT}/errors.txt"

                    # Add results to master file
                    if [ "$FIRST_RUN" = true ]; then
                        cat "$ERROR_FILE" > "$MASTER_OUTPUT"
                        FIRST_RUN=false
                    else
                        tail -n 1 "$ERROR_FILE" >> "$MASTER_OUTPUT"
                    fi

                done
            done
        done
    done
done

echo "========================================"
echo "All simulations completed."
echo "Total simulations: ${RUN}"
echo "Master errors: ${MASTER_OUTPUT}"
echo "========================================"