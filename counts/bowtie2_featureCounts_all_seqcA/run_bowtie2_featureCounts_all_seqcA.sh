#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate featureCounts_env

# Set parameters
INPUT_DIR="genetic_data/alignments/bowtie2_seqA_all"
GTF_FILE="genetic_data/annotations/Mus_musculus.GRCm38.102.gtf"
OUTPUT_DIR="genetic_data/counts/bowtie2_featureCounts_all_seqcA"
THREADS=4 

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through position sorted BAM files in input directory
for file in "$INPUT_DIR"/*_sorted.bam; do
    if [[ -f "$file" ]]; then
        filename=$(basename -- "$file")
        sample_name="${filename%_sorted.bam}"
        output_file="$OUTPUT_DIR/${sample_name}_featureCounts.txt"

        # featureCounts command
        echo "Running featureCounts for $filename"
        featureCounts -p -M -T "$THREADS" -a "$GTF_FILE" -o "$output_file" "$file"
    fi
done

echo "Processing complete!"