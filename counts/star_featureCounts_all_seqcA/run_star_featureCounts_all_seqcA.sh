#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_seqcA_all"
GTF_FILE="genetic_data/annotations/gencode.v19.annotation.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_featureCounts_all_seqcA"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
THREADS=4

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    # featureCount run with given parameters
    featureCounts -p -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt" "$BAM_FILE"

    featureCounts -p -M -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_best.txt" "$BAM_FILE"

    featureCounts -p -M -O -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt" "$BAM_FILE"

    featureCounts -p -M --fraction -O -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt" "$BAM_FILE"

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
