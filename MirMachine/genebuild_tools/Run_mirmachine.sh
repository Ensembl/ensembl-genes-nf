#!/bin/bash

# Simplified Run_mirmachine.sh
# Usage: ./Run_mirmachine.sh <path_to_csv>

# The Mirmachine pipeline is designed to run on a single species at a time.
# This is how to run the pipeline on codon for all species in the csv downloaded from rapid

# First set up the config file to ensure outputs are to the correct location
# Then run the pipeline
#
# Usage: Run_mirmachine.sh <path_to_csv>

INPUT_FILE="$1"
PIPELINE_PATH="/hps/software/users/ensembl/genebuild/ereboperezsilva/modenv/mirm/ensembl-genes-nf/MirMachine"
SCRIPT_PATH="${PIPELINE_PATH}/genebuild_tools/run_mirmachine_species.sh"

# Convert to tab-delimited if necessary
if head -n 1 "$INPUT_FILE" | grep -q ','; then
    tr ',' '\t' < "$INPUT_FILE" > "${INPUT_FILE}.tab"
    INPUT_FILE="${INPUT_FILE}.tab"
fi

# Process each line of the input file
while IFS=$'\t' read -r line; do
    species=$(echo "$line" | cut -f1 | tr -d '"')
    assembly=$(echo "$line" | cut -f7 | tr -d '"')
    if [ -n "$species" ] && [ -n "$assembly" ]; then
        echo "Processing: $species with assembly $assembly"
        sbatch --job-name="mirMachine_${assembly}" --output="logs/output_${assembly}.log" --mail-user="${USER}@ebi.ac.uk" "$SCRIPT_PATH" "${PIPELINE_PATH}/main.nf" "$species" "$assembly" 
    fi
done < "$INPUT_FILE"

echo "All jobs submitted."
