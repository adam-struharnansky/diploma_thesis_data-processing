#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# File paths setting
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_GENE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf"
OUTPUT_DIR="genetic_data/counts/star_bilatticeCount_rat"
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

    # Strandness, annotation file and feature choosing
    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        STRANDNESS="stranded"
        GTF_FILE="$GTF_MIRNA"
        FEATURE="gene"
    else
        STRANDNESS="unstranded"
        GTF_FILE="$GTF_GENE"
        FEATURE="exon"
    fi

    # Loop through all parameter settings
    for P in "${P_VALUES[@]}"; do
        P_SAFE=$(format_p_for_filename "$P")
        # bilatticeCount run with given parameters
        python "$SCRIPT" \
            --annotation_file "$GTF_FILE" \
            --alignments_file "$BAM_FILE" \
            --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_schweizer_sklar_p${P_SAFE}.txt" \
            --t_norm "schweizer_sklar" \
            --t_norm_param "$P" \
            --strandness "$STRANDNESS" \
            --feature "$FEATURE" \
            --seed 42 \
            --verbose
    done

    echo "Finished processing $SAMPLE_NAME"
done

echo "OK, all files processed!"
