#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_seqcA_all" 
GTF_FILE="genetic_data/annotations/gencode.v19.annotation.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_HTSeq_all_seqcA"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME..."

    # Run HTSeq-count with different multimapped read options
    htseq-count -f bam -r name -t exon -i gene_id --stranded=no \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count -f bam -r name -t exon -i gene_id --stranded=no --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count -f bam -r name -t exon -i gene_id --stranded=no --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"