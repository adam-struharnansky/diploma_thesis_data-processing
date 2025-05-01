#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate fastp_env

# File paths setting
BASE_DIR="genetic_data/reads/seqcA"

# Files traversal
find "$BASE_DIR" -type f -name "*_1.fastq.gz" | while read r1; do
    # Sample name extraction
    sample_dir=$(dirname "$r1")
    r1_name=$(basename "$r1")
    sample_prefix="${r1_name%_1.fastq.gz}"
    r2="${sample_dir}/${sample_prefix}_2.fastq.gz"

    # fastp run in paired-end mode
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${sample_dir}/${sample_prefix}_1.trimmed.fastq.gz" \
        -O "${sample_dir}/${sample_prefix}_2.trimmed.fastq.gz" \
        -h "${sample_dir}/${sample_prefix}_fastp.html" \
        -j "${sample_dir}/${sample_prefix}_fastp.json" \
        -w 4

    # Replacement of original files with trimmed ones
    mv "${sample_dir}/${sample_prefix}_1.trimmed.fastq.gz" "$r1"
    mv "${sample_dir}/${sample_prefix}_2.trimmed.fastq.gz" "$r2"

done

echo "OK, all files trimmed"
