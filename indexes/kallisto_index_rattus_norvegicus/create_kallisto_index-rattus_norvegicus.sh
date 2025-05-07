#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# File paths setting
TRANSCRIPTOME_FILE="genetic_data/transcriptomes/Rattus_norvegicus.Rnor_5.0.cdna.all.fa.gz"
OUTPUT_DIR="genetic_data/indexes/kallisto_index_rattus_norvegicus"
INDEX_FILE="${OUTPUT_DIR}/Rattus_norvegicus.kallisto.idx"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Kallisto indexing run
kallisto index -i "$INDEX_FILE" "$TRANSCRIPTOME_FILE"

echo "OK, index created!"
