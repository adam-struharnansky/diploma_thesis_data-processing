#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/genomes/

# Genome file download
wget -P ~/genetic_data/genomes/ https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

echo "OK, file downloaded!"