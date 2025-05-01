#!/bin/bash

# File paths setting
TARGET_DIR="genetic_data/reads/beers_all_bias/"

# Output directory creation if nonexistence
mkdir -p "$TARGET_DIR"

# BEERS reads download (from http://bioinf.itmat.upenn.edu/BEERS2/paper/index.html)
wget -P $TARGET_DIR https://s3.amazonaws.com/itmat.data/BEERS2/datasets/MouseLiver/fastqs/all_bias.tar.gz

# Data extraction
tar -xvzf "$TARGET_DIR/all_bias.tar.gz" -C "$TARGET_DIR"

# Directories creation and movement of data
for sample in S{1..8}; do
  mkdir -p "$TARGET_DIR/$sample"
  mv "$TARGET_DIR/data/all_bias/beers/results/${sample}_1.fastq" "$TARGET_DIR/$sample/"
  mv "$TARGET_DIR/data/all_bias/beers/results/${sample}_2.fastq" "$TARGET_DIR/$sample/"

  # Files zip
  gzip "$TARGET_DIR/$sample/${sample}_1.fastq"
  gzip "$TARGET_DIR/$sample/${sample}_2.fastq"
done

# Unnecessary files/directories removal
rm -rf "$TARGET_DIR/data"
rm "$TARGET_DIR/all_bias.tar.gz"

echo "OK, all files downloaded"
