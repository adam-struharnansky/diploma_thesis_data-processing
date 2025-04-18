#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_rat_all" 
GTF_FILE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_bilatticeCount_rat"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# List of t-norms to evaluate
T_NORMS=("drastic" "product" "minimum" "Łukasiewicz")

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    echo "Processing $SAMPLE_NAME with different t-norms..."

    # Detect strandness
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRANDNESS="stranded"
    else
        STRANDNESS="unstranded"
    fi

    for T_NORM in "${T_NORMS[@]}"; do
        echo "Using t-norm: $T_NORM"

        python "$SCRIPT" --annotation_file "$GTF_FILE" \
                         --alignments_file "$BAM_FILE" \
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_${T_NORM}.txt" \
                         --t_norm "$T_NORM" \
                         --paired_end \
                         --strandness "$STRANDNESS" \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"
