#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/annotations/

# Annotation file download
wget -P ~/genetic_data/annotations/ https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz

echo "OK, file downloaded."
