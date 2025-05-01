#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# File paths setting
INDEX="genetic_data/indexes/kallisto_index_homo_sapiens/Homo_sapiens.kallisto.idx"
READS="genetic_data/reads/seqcA"
OUTPUT="genetic_data/counts/kallisto_seqcA"

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
        # Sample output directory definitions
        sample_out="${OUTPUT}/${sample_name}"
        mkdir -p "$sample_out"

        # kallisto run with given parameters
        kallisto quant -i "$INDEX" -o "$sample_out" --threads=8 "$r1" "$r2"

        echo "Finished processing $sample_name"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "OK, all files processed!"
