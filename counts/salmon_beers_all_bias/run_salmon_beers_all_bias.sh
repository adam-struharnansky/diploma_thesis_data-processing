#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# File paths setting
INDEX="genetic_data/indexes/salmon_index_mus_musculus"
READS="genetic_data/reads/beers_all_bias"
OUTPUT="genetic_data/counts/salmon_beers_all_bias"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT"

# Loop through each subdirectory
for sample_dir in "$READS"/*/; do
    # Sample name extraction
    sample_name=$(basename "$sample_dir")

    # FASTQ files definitions
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    if [[ -f "$r1" && -f "$r2" ]]; then
        # Sample output directory definition
        sample_out="${OUTPUT}/${sample_name}"
        mkdir -p "$sample_out"

        # salmon run with given parameters
        salmon quant -i "$INDEX" -l SF -1 "$r1" -2 "$r2" -p 8 -o "$sample_out" --validateMappings

        echo "Finished processing $sample_name"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "OK, all files processed!"
