#!/bin/bash

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
    DIR=$(dirname "$FILE")  # Get the directory name
    fastq-dump --split-files --gzip -O "$DIR" "$FILE"  # Convert to FASTQ
done

# Step 3: Remove the original SRA files
echo "Removing SRA files..."
find "$BASE_DIR" -type f -name "*.sra" -delete
rm srr_list.txt

echo "All SEQC A data downloaded."