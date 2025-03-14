#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate hisat2_env

# Directories
INDEX="genetic_data/indexes/hisat2_index_homo_sapiens/genome_index"
READS="genetic_data/reads/seqcB"
OUTPUT="genetic_data/alignments/hisat2_seqcB_all"

# Loop through each sample directory
for sample_dir in "$READS"/*/; do
    # Extract sample name
    sample_name=$(basename "$sample_dir")

    # Paired-end FASTQ files
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    echo "Processing: $sample_name"

    # Check if both R1 and R2 exist
    if [[ -f "$r1" && -f "$r2" ]]; then
        # Run Hisat2 alignment
        hisat2 -x "$INDEX" \
               -1 "$r1" \
               -2 "$r2" \
               -k 10 \
               --very-sensitive \
               --threads 4 \
               --summary-file "${OUTPUT}/${sample_name}_summary.txt" \
               -S "${OUTPUT}/${sample_name}.sam"

        # Convert SAM to BAM
        samtools view -bS "${OUTPUT}/${sample_name}.sam" > "${OUTPUT}/${sample_name}.bam"

        # Sort BAM by name
        samtools sort -n "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted_by_name.bam"

        # Sort BAM file
        samtools sort "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted.bam"

        # Index the sorted BAM file
        samtools index "${OUTPUT}/${sample_name}_sorted.bam"

        # Step 6: Remove unnecessary files
        rm "${OUTPUT}/${sample_name}.sam" "${OUTPUT}/${sample_name}.bam"

        echo "Finished $sample_name"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "All samples processed."
