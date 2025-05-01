#!/bin/bash

# Output directory creation if nonexistence
mkdir -p genetic_data/transcriptomes/

# Transcriptome file download
wget -P genetic_data/transcriptomes/ https://ftp.ensembl.org/pub/release-77/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_5.0.cdna.all.fa.gz

echo "OK, all files downloaded!"
