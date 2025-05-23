#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_beers_all_bias_all" 
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_HTSeq_beers_all_bias"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    # Run HTSeq-count with different multimapped read options
    htseq-count -f bam -r name -t exon -i gene_id \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
