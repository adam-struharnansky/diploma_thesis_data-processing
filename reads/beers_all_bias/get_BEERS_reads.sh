#!/bin/bash

TARGET_DIR="genetic_data/reads/beers_all_bias/"

# Download reads from BEERS http://bioinf.itmat.upenn.edu/BEERS2/paper/index.html
wget -P $TARGET_DIR https://s3.amazonaws.com/itmat.data/BEERS2/datasets/MouseLiver/fastqs/all_bias.tar.gz

# Extract data
tar -xvzf "$TARGET_DIR/all_bias.tar.gz" -C "$TARGET_DIR"

# Create directories and move paired files
for sample in S{1..8}; do
  mkdir -p "$TARGET_DIR/$sample"
  mv "$TARGET_DIR/data/all_bias/beers/results/${sample}_1.fastq" "$TARGET_DIR/$sample/"
  mv "$TARGET_DIR/data/all_bias/beers/results/${sample}_2.fastq" "$TARGET_DIR/$sample/"
done

# Remove unnecessary files / directories
rm -rf "$TARGET_DIR/data"
rm "$TARGET_DIR/all_bias.tar.gz"