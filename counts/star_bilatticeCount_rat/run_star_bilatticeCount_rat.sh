#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_GENE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_bilatticeCount_rat"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py"

# List of basic t-norms
T_NORMS=("drastic" "product" "minimum" "Łukasiewicz")

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME with basic t-norms..."

    # Detect strandness and choose correct GTF file
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRANDNESS="stranded"
        GTF_FILE="$GTF_MIRNA"
    else
        STRANDNESS="unstranded"
        GTF_FILE="$GTF_GENE"
    fi

    for T_NORM in "${T_NORMS[@]}"; do
        echo "Using t-norm: $T_NORM"

        python "$SCRIPT" \
            --annotation_file "$GTF_FILE" \
            --alignments_file "$BAM_FILE" \
            --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_${T_NORM}.txt" \
            --t_norm "$T_NORM" \
            --strandness "$STRANDNESS" \
            --seed 42 \
            --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "All basic t-norm STAR-based processing complete!"
