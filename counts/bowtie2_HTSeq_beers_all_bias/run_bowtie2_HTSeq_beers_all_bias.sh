#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_beers_all_bias_all"
GTF_GENE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Mus_musculus.GRCm38.102.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_HTSeq_beers_all_bias"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)
    echo "Processing $SAMPLE_NAME..."

    # Choose GTF file based on sample type
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        GTF_FILE="$GTF_MIRNA"
    else
        GTF_FILE="$GTF_GENE"
    fi

    # Run HTSeq-count with different multimapped read options
    htseq-count -f bam -r name -t exon -i gene_id \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"
