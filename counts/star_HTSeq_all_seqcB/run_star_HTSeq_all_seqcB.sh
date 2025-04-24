#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_seqcB_all"
GTF_GENE="genetic_data/annotations/gencode.v19.annotation.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_gencode.v19.annotation.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_HTSeq_all_seqcB"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "🔄 Processing $SAMPLE_NAME..."

    # Select appropriate GTF file
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        GTF_FILE="$GTF_MIRNA"
    else
        GTF_FILE="$GTF_GENE"
    fi

    # Run HTSeq-count with different multimapping handling modes
    htseq-count -f bam -r name -t exon -i gene_id --stranded=no \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count -f bam -r name -t exon -i gene_id --stranded=no --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count -f bam -r name -t exon -i gene_id --stranded=no --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All HTSeq runs complete!"
