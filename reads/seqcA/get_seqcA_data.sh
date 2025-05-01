#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

# File paths setting
BASE_DIR=genetic_data/reads/seqcA

# Output directory creation if nonexistence
mkdir -p "$BASE_DIR"

# Parameter setting
echo -e "SRR896663\nSRR896665\nSRR896667\nSRR896669\nSRR896671\nSRR896673\nSRR896675\nSRR896677" > srr_list.txt

# Files download
while read SRR; do
    prefetch "$SRR" -O "$BASE_DIR/$SRR"
done < srr_list.txt

# Files conversion
find "$BASE_DIR" -type f -name "*.sra" | while read FILE; do
    # SRR number extraction
    SRR=$(basename "$FILE" .sra)
    # SRR directory creation
    mkdir -p "$BASE_DIR/$SRR"
    # SRR to FASTQ convertion
    fastq-dump --split-files --gzip -O "$BASE_DIR/$SRR" "$FILE"
    # Original file removal
    rm "$FILE"
done

# Removal of temporary files
find "$BASE_DIR" -type f -name "*.sra" -delete
rm srr_list.txt

echo "OK, all files downloaded!"
