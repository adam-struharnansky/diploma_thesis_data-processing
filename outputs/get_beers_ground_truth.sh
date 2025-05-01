#!/bin/bash

# File paths setting
DEST_DIR="genetic_data/outputs"

# Output directory creation if nonexistence
mkdir -p "$DEST_DIR"

# Ground truth file download
wget https://s3.amazonaws.com/itmat.data/BEERS2/datasets/MouseLiver/beers.true_TPM.parquet -P "$DEST_DIR"

echo "OK, file downloaded!"
