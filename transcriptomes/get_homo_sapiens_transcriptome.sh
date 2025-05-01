#!/bin/bash

# Output directory creation if nonexistence
mkdir -p genetic_data/transcriptomes/

# Trancriptome file download
wget -P genetic_data/transcriptomes/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz

echo "OK, all files downloaded!"
