#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

# File paths setting
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/gencode.v19.pc_transcripts.fa.gz"
OUTPUT_DIR="genetic_data/indexes/salmon_index_homo_sapiens"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Salmon indexing run
salmon index -t "$TRANSCRIPTOME_FILE" -i "$OUTPUT_DIR" --gencode
