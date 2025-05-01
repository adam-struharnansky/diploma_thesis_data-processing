#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# File paths setting
INPUT_DIR="genetic_data/alignments/bowtie2_beers_all_bias_all" 
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_bilatticeCount_beers_all_bias"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
P_VALUES=(-4 -2 -1 -0.5 -0.1 0.1 0.5 2 4 8 12)

# Name formatting function
format_p_for_filename() {
    local p=$1
    if (( $(echo "$p < 0" | bc -l) )); then
        echo "n$(echo ${p#-} | sed 's/\.//')"
    else
        echo "$(echo $p | sed 's/\.//')"
    fi
}

# Loop through name-sorted BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    # Sample name extraction
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    # Loop through all parameter settings
    for P in "${P_VALUES[@]}"; do
        P_SAFE=$(format_p_for_filename "$P")
        # bilatticeCount run with given parameters
        python "$SCRIPT" --annotation_file "$GTF_FILE" \
                         --alignments_file "$BAM_FILE" \
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_schweizer_sklar_p${P_SAFE}.txt" \
                         --t_norm "schweizer_sklar" \
                         --t_norm_param "$P" \
                         --paired_end \
                         --strandness stranded \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, All files processed!"
