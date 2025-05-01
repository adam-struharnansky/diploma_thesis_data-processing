#!/bin/bash

# Output directory creation if nonexistence
mkdir -p genetic_data/transcriptomes/

# Transcriptome file download
wget -P genetic_data/transcriptomes/ https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

echo "OK, all files downloaded!"
