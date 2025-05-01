#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_beers_all_bias_all" 
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz"
OUTPUT_DIR="genetic_data/counts/star_bilatticeCount_beers_all_bias"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py" 

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
P_VALUES=(0.0 0.2 0.5 2.0 5.0 10.0 20.0)

# Name formatting function
format_p_for_filename() {
    local p=$1
    echo "$(echo $p | sed 's/\.//')"
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
                         --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_hamacher_p${P_SAFE}.txt" \
                         --t_norm "hamacher" \
                         --t_norm_param "$P" \
                         --paired_end \
                         --strandness stranded \
                         --seed 42 \
                         --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
