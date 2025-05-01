#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/genomes/

# Genome file download
wget -P ~/genetic_data/genomes/ https://ftp.ensembl.org/pub/release-77/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.dna.toplevel.fa.gz

echo "OK, file downloaded!"