#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_seqcB_all" 
GTF_FILE="genetic_data/annotations/gencode.v19.annotation.gtf"
OUTPUT_DIR="genetic_data/counts/star_bilatticeCount_all_seqcB"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# Schweizer–Sklar t-norm parameters to try
P_VALUES=(-4 -2 -1 -0.5 -0.1 0.1 0.5 2 4)

format_p_for_filename() {
    local p=$1
    if (( $(echo "$p < 0" | bc -l) )); then
        echo "n$(echo ${p#-} | sed 's/\.//')"
    else
        echo "$(echo $p | sed 's/\.//')"
    fi
}

# Loop through BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)
    echo "Processing $SAMPLE_NAME with Schweizer–Sklar t-norms..."
    
    for P in "${P_VALUES[@]}"; do
        echo "Using Schweizer–Sklar t-norm with p=$P"
        P_SAFE=$(format_p_for_filename "$P")
        python "$SCRIPT" --annotation_file "$GTF_FILE" \
                         --alignments_file "$BAM_FILE" \
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_schweizer_sklar_p${P_SAFE}.txt" \
                         --t_norm "schweizer_sklar" \
                         --t_norm_param "$P" \
                         --paired_end \
                         --strandness unstranded \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"
