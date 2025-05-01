#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/annotations/

# Annotation file download
wget -P ~/genetic_data/annotations/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

echo "OK, file downloaded"
