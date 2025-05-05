#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_GENE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.all.gtf"
OUTPUT_DIR="genetic_data/counts/star_featureCounts_rat"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
THREADS=4

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME with featureCounts..."

    # Strandness and annotation file choosing
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRAND_OPTION=1
        GTF_FILE="$GTF_MIRNA"
    else
        STRAND_OPTION=0
        GTF_FILE="$GTF_GENE"
    fi

    # featureCount run with given parameters
    featureCounts -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" \
        -o "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt" "$BAM_FILE"

    featureCounts -M -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" \
        -o "$OUTPUT_DIR/${SAMPLE_NAME}_best.txt" "$BAM_FILE"

    featureCounts -M -O -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" \
        -o "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt" "$BAM_FILE"

    featureCounts -M --fraction -O -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" \
        -o "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt" "$BAM_FILE"

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
