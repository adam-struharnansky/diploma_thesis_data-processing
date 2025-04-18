#!/bin/bash

# Activate fastp environment
source miniconda3/etc/profile.d/conda.sh
conda activate fastp_env

# Get the base directory as input
BASE_DIR="genetic_data/reads/seqcA"

# Traverse all subdirectories for *_1.fastq.gz files
find "$BASE_DIR" -type f -name "*_1.fastq.gz" | while read r1; do
    sample_dir=$(dirname "$r1")
    r1_name=$(basename "$r1")
    sample_prefix="${r1_name%_1.fastq.gz}"
    r2="${sample_dir}/${sample_prefix}_2.fastq.gz"

    echo "Processing $sample_prefix in $sample_dir..."

    # Run fastp in paired-end mode
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${sample_dir}/${sample_prefix}_1.trimmed.fastq.gz" \
        -O "${sample_dir}/${sample_prefix}_2.trimmed.fastq.gz" \
        -h "${sample_dir}/${sample_prefix}_fastp.html" \
        -j "${sample_dir}/${sample_prefix}_fastp.json" \
        -w 4

    # Replace original files with trimmed ones
    mv "${sample_dir}/${sample_prefix}_1.trimmed.fastq.gz" "$r1"
    mv "${sample_dir}/${sample_prefix}_2.trimmed.fastq.gz" "$r2"

    echo "Finished trimming $sample_prefix"
done

echo "All paired-end FASTQ files trimmed and replaced."
