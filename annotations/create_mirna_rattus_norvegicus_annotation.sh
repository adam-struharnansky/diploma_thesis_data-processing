#!/bin/bash

INPUT_GTF="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
TEMP_GTF="genetic_data/annotations/temp_miRNA.gtf"
OUTPUT_GTF="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf"

echo "Extracting miRNA entries from $INPUT_GTF ..."
zgrep 'gene_biotype "miRNA"' "$INPUT_GTF" > "$TEMP_GTF"

echo "Filtering to keep only 'gene' entries ..."
awk '$3 == "gene"' "$TEMP_GTF" > "$OUTPUT_GTF"

# Cleaning up temporary file
rm "$TEMP_GTF"
