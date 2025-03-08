#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate bowtie2_env

# Directories
INDEX="genetic_data/indexes/bowtie2_index_homo_sapiens/bowtie2_homo_sapiens_index"
READS="genetic_data/reads/seqcB"
OUTPUT="genetic_data/alignments/bowtie2_seqcB_all"

# Loop through each subdirectory
for sample_dir in "$READS"/*/; do
    # Extract sample name
    echo "$sample_dir"
    sample_name=$(basename "$sample_dir")

    # Define paired-end FASTQ files
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    echo "$sample_name"

    # Check if both R1 and R2 exist before running bowtie
    if [[ -f "$r1" && -f "$r2" ]]; then
        # Run Bowtie2 alignment
        bowtie2 -x "$INDEX" -1 "$r1" -2 "$r2" --very-sensitive -k 10 -S "${OUTPUT}/${sample_name}.sam"

        # Convert SAM to BAM
        samtools view -Sb "${OUTPUT}/${sample_name}.sam" > "${OUTPUT}/${sample_name}.bam"

        # Sort BAM by name
        samtools sort -n "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted_by_name.bam"

        # Sort the BAM file
        samtools sort "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted.bam"

        # Index the sorted BAM file
        samtools index "${OUTPUT}/${sample_name}_sorted.bam"

        # Remove the unsorted BAM and SAM files
        rm "${OUTPUT}/${sample_name}.sam" "${OUTPUT}/${sample_name}.bam"

        echo "Finished $sample_name. Results stored in $OUTPUT"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "All samples processed."
