#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# Set parameters
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/Mus_musculus.GRCm38.cdna.all.fa.gz"
OUTPUT_DIR="genetic_data/indexes/salmon_index_mus_musculus"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run Salmon index
salmon index -t "$TRANSCRIPTOME_FILE" -i "$OUTPUT_DIR" --gencode
