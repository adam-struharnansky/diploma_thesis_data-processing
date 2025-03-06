#!/bin/bash

# Define the base directory for storing data
BASE_DIR=~/genetic_data/reads/seqcB

# Create the SRR list
echo -e "SRR896743\nSRR896745\nSRR896747\nSRR896749\nSRR896751\nSRR896753\nSRR896755\nSRR896757" > srr_list_b.txt

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