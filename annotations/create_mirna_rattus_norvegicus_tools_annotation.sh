#!/bin/bash

# File paths setting
INPUT_GTF="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_GTF="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.all.gtf"

# Filtering only miRNA
zgrep 'gene_biotype "miRNA"' "$INPUT_GTF" > "$OUTPUT_GTF"

echo "OK, file created."
