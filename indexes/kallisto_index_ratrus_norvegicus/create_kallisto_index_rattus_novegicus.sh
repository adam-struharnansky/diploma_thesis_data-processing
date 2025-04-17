#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# Set parameters
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
OUTPUT_DIR="genetic_data/indexes/kallisto_index_rattus_norvegicus"
INDEX_FILE="${OUTPUT_DIR}/Rattus_norvegicus.kallisto.idx"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run Kallisto index
kallisto index -i "$INDEX_FILE" "$TRANSCRIPTOME_FILE"
