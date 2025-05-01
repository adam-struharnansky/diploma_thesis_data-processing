#!/bin/bash

# File paths setting
INPUT_GTF="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
TEMP_GTF="genetic_data/annotations/temp_miRNA.gtf"
OUTPUT_GTF="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf"

# Filtering only miRNA
zgrep 'gene_biotype "miRNA"' "$INPUT_GTF" > "$TEMP_GTF"

# Filtering only genes
awk '$3 == "gene"' "$TEMP_GTF" > "$OUTPUT_GTF"

# Removal of temporary file
rm "$TEMP_GTF"

echo "OK, file created."
