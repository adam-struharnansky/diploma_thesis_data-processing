#!/bin/bash

# File paths setting
DEST_DIR="genetic_data/outputs"

# Output directory creation if nonexistence
mkdir -p "$DEST_DIR"

# RT-PCR data file download
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56457/suppl/GSE56457_SEQC_qPCR_GEOSub.txt.gz -P "$DEST_DIR"

# RT-PCR data unzip
gunzip "$DEST_DIR/GSE56457_SEQC_qPCR_GEOSub.txt.gz"

echo "OK, file downloaded!"
