#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_beers_all_bias_all"
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_featureCounts_beers_all_bias"
THREADS=4

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Extract the sample name
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME..."

    # Run featureCounts with different multimapped read options
    featureCounts -p -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt" "$BAM_FILE"

    featureCounts -p -M -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_best.txt" "$BAM_FILE"

    featureCounts -p -M -O -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt" "$BAM_FILE"

    featureCounts -p -M --fraction -O -T "$THREADS" -a "$GTF_FILE" -o "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt" "$BAM_FILE"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"
