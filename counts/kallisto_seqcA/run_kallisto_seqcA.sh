#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

INDEX="genetic_data/indexes/kallisto_index_homo_sapiens/Homo_sapiens.kallisto.idx"
READS="genetic_data/reads/seqcA"
OUTPUT="genetic_data/counts/kallisto_seqcA"

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

        # Run Kallisto quant
        kallisto quant -i "$INDEX" -o "$sample_out" --threads=8 "$r1" "$r2"

        echo "Finished processing $sample_name. Results stored in $sample_out"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "All samples processed."
