#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bowtie2_env

# File paths setting
READS_DIR="genetic_data/reads/rat"
OUTPUT_DIR="genetic_data/alignments/bowtie2_rat_all"
INDEX="genetic_data/indexes/bowtie2_index_rattus_norvegicus/bowtie2_rattus_norvegicus_index"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Loop through each sample
for fastq_file in "$READS_DIR"/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    # Bowtie2 alignment run
    bowtie2 -x "$INDEX" -U "$fastq_file" --very-sensitive -k 10 -S "${OUTPUT_DIR}/${sample_name}.sam"

    # SAM to BAM conversion
    samtools view -Sb "${OUTPUT_DIR}/${sample_name}.sam" > "${OUTPUT_DIR}/${sample_name}.bam"

    # BAM by name sorting
    samtools sort -n "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted_by_name.bam"

    # BAM by coordinates sorting
    samtools sort "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted.bam"

     # BAM file indexation
    samtools index "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Unsorted BAM and SAM files removal
    rm "${OUTPUT_DIR}/${sample_name}.sam" "${OUTPUT_DIR}/${sample_name}.bam"
done

echo "OK, all samples processed."
