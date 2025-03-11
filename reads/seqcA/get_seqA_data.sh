#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

# Define the base directory for storing data
BASE_DIR=~/genetic_data/reads/seqcA

# Create the SRR list
echo -e "SRR896663\nSRR896665\nSRR896667\nSRR896669\nSRR896671\nSRR896673\nSRR896675\nSRR896677" > srr_list.txt

# Create the base directory if it doesn't exist
mkdir -p "$BASE_DIR"

# Step 1: Download all files using prefetch
echo "Downloading SRA files..."
while read SRR; do
    prefetch "$SRR" -O "$BASE_DIR/$SRR"
done < srr_list.txt

# Step 2: Convert .sra files to FASTQ
echo "Converting SRA to FASTQ..."
find "$BASE_DIR" -type f -name "*.sra" | while read FILE; do
    SRR=$(basename "$FILE" .sra)  # Extract SRR number
    mkdir -p "$BASE_DIR/$SRR"  # Create SRR directory
    fastq-dump --split-files --gzip -O "$BASE_DIR/$SRR" "$FILE"  # Convert to FASTQ
    rm "$FILE"  # Remove the original SRA file after conversion
done

rm srr_list.txt

echo "All SEQC A data downloaded correctly."