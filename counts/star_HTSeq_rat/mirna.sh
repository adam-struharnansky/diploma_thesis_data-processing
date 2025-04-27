#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate htseq_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_HTSeq_rat"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        echo "🔁 Reprocessing miRNA sample $SAMPLE_NAME with HTSeq..."

        STRAND_OPTION="yes"
        GTF_FILE="$GTF_MIRNA"

        # Run HTSeq-count with different multimapped read handling
        htseq-count -f bam -r name -s "$STRAND_OPTION" -t exon -i gene_id \
            "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_ignore.txt"

        htseq-count -f bam -r name -s "$STRAND_OPTION" -t exon -i gene_id --nonunique all \
            "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_all.txt"

        htseq-count -f bam -r name -s "$STRAND_OPTION" -t exon -i gene_id --nonunique fraction \
            "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/${SAMPLE_NAME}_fraction.txt"

        echo "✅ Finished reprocessing $SAMPLE_NAME"
    fi
done

echo "🎯 miRNA HTSeq reruns complete!"
