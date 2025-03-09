#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# Set parameters
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/gencode.v19.pc_transcripts.fa.gz"
OUTPUT_DIR="genetic_data/indexes/kallisto_index_homo_sapiens"
INDEX_FILE="${OUTPUT_DIR}/Homo_sapiens.kallisto.idx"

# Run Kallisto index
kallisto index -i "$INDEX_FILE" "$TRANSCRIPTOME_FILE"
