#!/bin/bash

# Destination directory
DEST_DIR="genetic_data/outputs"

# Create the directory if it doesn't exist
mkdir -p "$DEST_DIR"

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56457/suppl/GSE56457_SEQC_qPCR_GEOSub.txt.gz -P "$DEST_DIR"
gunzip "$DEST_DIR/GSE56457_SEQC_qPCR_GEOSub.txt.gz"

echo "Download completed! Files are saved in: $DEST_DIR"
