#!/bin/bash

# Destination directory
DEST_DIR="genetic_data/outputs"

mkdir -p "$DEST_DIR"

wget https://s3.amazonaws.com/itmat.data/BEERS2/datasets/MouseLiver/beers.true_TPM.parquet -P "$DEST_DIR"

echo "Download completed! Files are saved in: $DEST_DIR"