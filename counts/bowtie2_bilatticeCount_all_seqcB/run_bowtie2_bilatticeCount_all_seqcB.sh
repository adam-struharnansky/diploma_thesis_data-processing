#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# File paths setting
INPUT_DIR="genetic_data/alignments/bowtie2_seqcB_all" 
GTF_FILE="genetic_data/annotations/gencode.v19.annotation.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_bilatticeCount_all_seqcB"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
T_NORMS=("drastic" "product" "minimum" "Łukasiewicz")

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    # Loop through all parameter settings
    for T_NORM in "${T_NORMS[@]}"; do
        # bilatticeCount run with given parameters
        python "$SCRIPT" --annotation_file "$GTF_FILE" \
                         --alignments_file "$BAM_FILE" \
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_${T_NORM}.txt" \
                         --t_norm "$T_NORM" \
                         --paired_end \
                         --strandness unstranded \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
