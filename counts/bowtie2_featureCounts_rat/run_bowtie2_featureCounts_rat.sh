#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_rat_all"
GTF_FILE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_featureCounts_rat"
THREADS=4

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through position sorted BAM files
for file in "$INPUT_DIR"/*_sorted.bam; do
    if [[ -f "$file" ]]; then
        filename=$(basename -- "$file")
        sample_name="${filename%_sorted.bam}"
        output_file="$OUTPUT_DIR/${sample_name}_featureCounts.txt"

        # Set strandness: 1 if mirna, otherwise 0
        if [[ "$sample_name" == *_mirna_* ]]; then
            STRAND_OPTION=1
        else
            STRAND_OPTION=0
        fi

        echo "Running featureCounts for $filename (strandness: $STRAND_OPTION)"
        featureCounts -M -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" -o "$output_file" "$file"
    fi
done

echo "Processing complete!"
