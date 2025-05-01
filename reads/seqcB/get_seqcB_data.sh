#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

# File paths setting
BASE_DIR=genetic_data/reads/seqcB

# Output directory creation if nonexistence
mkdir -p "$BASE_DIR"

# Parameter setting
echo -e "SRR896743\nSRR896745\nSRR896747\nSRR896749\nSRR896751\nSRR896753\nSRR896755\nSRR896757" > srr_list.txt

# Files download
echo "Downloading SRA files..."
while read SRR; do
    prefetch "$SRR" -O "$BASE_DIR/$SRR"
done < srr_list.txt

# Files conversion
echo "Converting SRA to FASTQ..."
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
