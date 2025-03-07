#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_seqA_all" 
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf"
OUTPUT_DIR="genetic_data/counts/bowtie2_featureCounts_all_seqcA"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through name sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Extract the sample name (removes the directory path and suffix)
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME..."

    # Run HTSeq-count with different multimapped read options
    htseq-count -f bam -r name -t exon -i gene_id \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    htseq-count -f bam -r name -t exon -i gene_id --nonunique random \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_random.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"