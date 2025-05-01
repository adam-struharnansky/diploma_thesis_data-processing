#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/genomes/

# Genome file download
wget -P ~/genetic_data/genomes/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

echo "OK, file downloaded!"
