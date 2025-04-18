#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate star_env

# Directories
READS_DIR="genetic_data/reads/rat"
OUTPUT_DIR="genetic_data/alignments/star_rat_all"
INDEX="genetic_data/indexes/star_index_rattus_norvegicus"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each FASTQ file
for fastq_file in "$READS_DIR"/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)
    echo "Processing $sample_name..."

    # Run STAR alignment for single-end reads
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

    # Rename the BAM file
    mv "${OUTPUT_DIR}/${sample_name}_Aligned.out.bam" "${OUTPUT_DIR}/${sample_name}.bam"

    # Sort BAM by name
    samtools sort -n "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted_by_name.bam"

    # Sort BAM by position
    samtools sort "${OUTPUT_DIR}/${sample_name}.bam" -o "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Index the sorted BAM file
    samtools index "${OUTPUT_DIR}/${sample_name}_sorted.bam"

    # Remove the original unsorted BAM
    rm "${OUTPUT_DIR}/${sample_name}.bam"

    echo "Finished $sample_name. Results saved in $OUTPUT_DIR"
done

echo "All STAR alignments finished."
