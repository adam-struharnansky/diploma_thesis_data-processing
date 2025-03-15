#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# Define paths
INDEX="genetic_data/indexes/salmon_index_mus_musculus"
READS="genetic_data/reads/beers_all_bias"
OUTPUT="genetic_data/counts/salmon_beers_all_bias"

# Loop through each subdirectory
for sample_dir in "$READS"/*/; do
    # Extract sample name
    echo "$sample_dir"
    sample_name=$(basename "$sample_dir")

    # Define paired-end FASTQ files
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    if [[ -f "$r1" && -f "$r2" ]]; then
        # Define output directory for this sample
        sample_out="${OUTPUT}/${sample_name}"
        mkdir -p "$sample_out"

        # Run Salmon quant
        salmon quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" -p 8 -o "$sample_out" --validateMappings

        echo "Finished processing $sample_name. Results stored in $sample_out"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "All samples processed."
