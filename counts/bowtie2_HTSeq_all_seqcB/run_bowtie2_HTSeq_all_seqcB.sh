#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_seqcB_all" 
GTF_FILE="genetic_data/annotations/gencode.v19.annotation.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_HTSeq_all_seqcB"

# Loop through name sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Extract the sample name (removes the directory path and suffix)
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME..."

    # Run HTSeq-count with different multimapped read options
    htseq-count --paired-end -f bam -r name -t exon -i gene_id --stranded=no \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

    htseq-count --paired-end -f bam -r name -t exon -i gene_id --stranded=no --nonunique all \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

    htseq-count --paired-end -f bam -r name -t exon -i gene_id --stranded=no --nonunique fraction \
        "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"