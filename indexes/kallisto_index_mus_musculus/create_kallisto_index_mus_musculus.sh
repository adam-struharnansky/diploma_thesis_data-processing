#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# File paths setting
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/Mus_musculus.GRCm38.cdna.all.fa.gz"
OUTPUT_DIR="genetic_data/indexes/kallisto_index_mus_musculus"
INDEX_FILE="${OUTPUT_DIR}/Mus_musculus.kallisto.idx"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Kallisto indexing run
kallisto index -i "$INDEX_FILE" "$TRANSCRIPTOME_FILE"

echo "OK, index created!"
