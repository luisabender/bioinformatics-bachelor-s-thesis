#!/bin/bash

# Define the root directories for single and paired alignments
ROOT_DIR="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/kallisto_paired"

# Output file
OUTPUT_FILE="kallisto_alignment_summary.csv"

# Write the header to the output file
echo "sample_ID,p_pseudoaligned,n_pseudoaligned, n_processed" > "$OUTPUT_FILE"

# Loop through the root directories
#for ROOT_DIR in "${ROOT_DIRS[@]}"; do
    # Find all run_info.json files in the directory
    find "$ROOT_DIR" -type f -name "run_info.json" | while read -r JSON_FILE; do
        # Extract the sample name from the file path
        SAMPLE_NAME=$(basename "$(dirname "$JSON_FILE")")

        # Extract values for p_pseudoaligned and n_pseudoaligned using grep and sed
        P_PSEUDOALIGNED=$(grep '"p_pseudoaligned"' "$JSON_FILE" | sed -E 's/.*: ([0-9.]+),?/\1/')
        N_PSEUDOALIGNED=$(grep '"n_pseudoaligned"' "$JSON_FILE" | sed -E 's/.*: ([0-9]+),?/\1/')
        N_PROCESSED=$(grep '"n_processed"' "$JSON_FILE" | sed -E 's/.*: ([0-9]+),?/\1/')

        # Append the results to the output file
        echo "$SAMPLE_NAME,$P_PSEUDOALIGNED,$N_PSEUDOALIGNED, $N_PROCESSED" >> "$OUTPUT_FILE"
    done
#done
