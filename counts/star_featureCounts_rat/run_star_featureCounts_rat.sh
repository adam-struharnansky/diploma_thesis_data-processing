#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_GENE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_featureCounts_rat"
THREADS=4

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME with featureCounts..."

    # Detect strandness and select GTF
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRAND_OPTION=1
        GTF_FILE="$GTF_MIRNA"
    else
        STRAND_OPTION=0
        GTF_FILE="$GTF_GENE"
    fi

    # Run featureCounts with various multimapping options
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

echo "All featureCounts runs complete!"
