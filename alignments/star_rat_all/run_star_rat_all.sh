#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate star_env

# File paths setting
READS_DIR="genetic_data/reads/rat"
OUTPUT_DIR="genetic_data/alignments/star_rat_all"
INDEX="genetic_data/indexes/star_index_rattus_norvegicus"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT_DIR"

# Loop through each sample
for fastq_file in "$READS_DIR"/*.fastq.gz; do
    # Sample name extraction
    sample_name=$(basename "$fastq_file" .fastq.gz)

    # STAR alignment run
    STAR --runThreadN 8 \
         --genomeDir "$INDEX" \
         --readFilesIn "$fastq_file" \
         --readFilesCommand zcat \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix "${OUTPUT_DIR}/${sample_name}_" \
         --outFilterMultimapNmax 10 \
         --outSAMattributes NH HI AS nM \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000

    # Renaming to fit convention
    mv "${OUTPUT_DIR}/${sample_name}_Aligned.out.bam" "${OUTPUT_DIR}/${sample_name}.bam"

    # BAM by name sorting
    samtools sort -n "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted_by_name.bam"

    # BAM by coordinates sorting
    samtools sort "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted.bam"

     # BAM file indexation
    samtools index "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Unsorted BAM file removal
    rm "${OUTPUT_DIR}/${sample_name}.bam"

done

echo "OK, all samples processed."
