#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_rat_all" 
GTF_FILE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_bilatticeCount_rat"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# Schweizer–Sklar t-norm parameters to try
P_VALUES=(-4 -2 -1 -0.5 -0.1 0.1 0.5 2 4 8 12)

# Format float p-value for safe filenames
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

    # Detect strandness
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRANDNESS="stranded"
    else
        STRANDNESS="unstranded"
    fi

    for P in "${P_VALUES[@]}"; do
        echo "Using Schweizer–Sklar t-norm with p=$P"
        P_SAFE=$(format_p_for_filename "$P")

        python "$SCRIPT" --annotation_file "$GTF_FILE" \
                         --alignments_file "$BAM_FILE" \
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_schweizer_sklar_p${P_SAFE}.txt" \
                         --t_norm "schweizer_sklar" \
                         --t_norm_param "$P" \
                         --strandness "$STRANDNESS" \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "All files processed!"
