#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# File paths setting
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/Mus_musculus.GRCm38.cdna.all.fa.gz"
OUTPUT_DIR="genetic_data/indexes/salmon_index_mus_musculus"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Salmon indexing run
salmon index -t "$TRANSCRIPTOME_FILE" -i "$OUTPUT_DIR" --gencode

echo "OK, index created!"
