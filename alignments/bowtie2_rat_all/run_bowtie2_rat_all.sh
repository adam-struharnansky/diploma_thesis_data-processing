#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate bowtie2_env

# Directories
READS_DIR="genetic_data/reads/rat"
OUTPUT_DIR="genetic_data/alignments/bowtie2_rat_all"
INDEX="genetic_data/indexes/bowtie2_index_rattus_norvegicus/bowtie2_rattus_norvegicus"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each FASTQ file
for fastq_file in "$READS_DIR"/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)
    echo "Processing $sample_name..."

    # Run Bowtie2 alignment for single-end reads
    bowtie2 -x "$INDEX" -U "$fastq_file" --very-sensitive -k 10 -S "${OUTPUT_DIR}/${sample_name}.sam"

    # Convert SAM to BAM
    samtools view -Sb "${OUTPUT_DIR}/${sample_name}.sam" > "${OUTPUT_DIR}/${sample_name}.bam"

    # Sort BAM by name
    samtools sort -n "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted_by_name.bam"

    # Sort BAM by coordinates
    samtools sort "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Index BAM
    samtools index "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Clean up intermediate files
    rm "${OUTPUT_DIR}/${sample_name}.sam" "${OUTPUT_DIR}/${sample_name}.bam"

    echo "Finished $sample_name. Results saved in $OUTPUT_DIR"
done

echo "All samples aligned."
