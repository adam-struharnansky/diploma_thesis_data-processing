#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate feature_counts_env

# File paths setting
INPUT_DIR="genetic_data/alignments/bowtie2_rat_all"
GTF_GENE="genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
GTF_MIRNA="genetic_data/annotations/miRNA_Rattus_norvegicus.Rnor_5.0.77.gtf.gz"
OUTPUT_DIR="genetic_data/counts/bowtie2_featureCounts_rat"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Parameter setting
THREADS=4

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through position-sorted BAM files
for file in "$INPUT_DIR"/*_sorted.bam; do
    if [[ -f "$file" ]]; then
        # Sample name extractions
        filename=$(basename -- "$file")
        sample_name="${filename%_sorted.bam}"
        output_file="$OUTPUT_DIR/${sample_name}_featureCounts.txt"

        # Strandness and annotation file choosing
        if [[ "$sample_name" == *_mirna_* ]]; then
            STRAND_OPTION=1
            GTF_FILE="$GTF_MIRNA"
        else
            STRAND_OPTION=0
            GTF_FILE="$GTF_GENE"
        fi

        # featureCount run with given parameters
        featureCounts -M -T "$THREADS" -s "$STRAND_OPTION" -a "$GTF_FILE" -o "$output_file" "$file"
        echo "Finished processing $file"
    fi
done

echo "OK, all files processed!"
