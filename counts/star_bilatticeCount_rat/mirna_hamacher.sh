#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate bilattice_env

# Set parameters
INPUT_DIR="genetic_data/alignments/star_rat_all"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf"
OUTPUT_DIR="genetic_data/counts/star_bilatticeCount_rat"
SCRIPT="genetic_data/bilattice_count/bilattice_count.py"

# Hamacher t-norm parameters to try
P_VALUES=(0.0 0.2 0.5 2.0 5.0 10.0 20.0)

# Helper to format p-values for filenames
format_p_for_filename() {
    local p=$1
    echo "$(echo $p | sed 's/\.//')"
}

# Loop through BAM files
for BAM_FILE in "$INPUT_DIR"/*_sorted_by_name.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" _sorted_by_name.bam)

    if [[ "$SAMPLE_NAME" == *_mirna_* ]]; then
        echo "🔁 Reprocessing $SAMPLE_NAME with Hamacher t-norms (miRNA only)..."

        STRANDNESS="stranded"
        GTF_FILE="$GTF_MIRNA"
        FEATURE="gene"

        for P in "${P_VALUES[@]}"; do
            echo " → Using Hamacher t-norm with p=$P"
            P_SAFE=$(format_p_for_filename "$P")

            python "$SCRIPT" \
                --annotation_file "$GTF_FILE" \
                --alignments_file "$BAM_FILE" \
                --output_file "$OUTPUT_DIR/${SAMPLE_NAME}_hamacher_p${P_SAFE}.txt" \
                --t_norm "hamacher" \
                --t_norm_param "$P" \
                --strandness "$STRANDNESS" \
                --feature "$FEATURE" \
                --seed 42 \
                --verbose
        done

        echo "✅ Finished reprocessing $SAMPLE_NAME"
    fi
done

echo "🎯 miRNA STAR-based Hamacher reruns complete!"
