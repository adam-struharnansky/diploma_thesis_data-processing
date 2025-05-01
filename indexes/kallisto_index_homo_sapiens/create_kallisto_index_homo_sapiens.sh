#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# File paths setting
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/gencode.v19.pc_transcripts.fa.gz"
OUTPUT_DIR="genetic_data/indexes/kallisto_index_homo_sapiens"
INDEX_FILE="${OUTPUT_DIR}/Homo_sapiens.kallisto.idx"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Kallisto index run
kallisto index -i "$INDEX_FILE" "$TRANSCRIPTOME_FILE"

echo "OK, index created!"
