#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bowtie2_env

# File paths setting
INDEX="genetic_data/indexes/bowtie2_index_homo_sapiens/bowtie2_homo_sapiens_index"
READS="genetic_data/reads/seqcB"
OUTPUT="genetic_data/alignments/bowtie2_seqcB_all"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT"

# Loop through each sample
for sample_dir in "$READS"/*/; do
    # Sample name extraction
    sample_name=$(basename "$sample_dir")

    # FASTQ files definitions
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    echo "$sample_name"

    # Both files must exist
    if [[ -f "$r1" && -f "$r2" ]]; then
        # Bowtie2 alignment run
        bowtie2 -x "$INDEX" -1 "$r1" -2 "$r2" --very-sensitive -k 10 -S "${OUTPUT}/${sample_name}.sam"

        # SAM to BAM conversion
        samtools view -Sb "${OUTPUT}/${sample_name}.sam" > "${OUTPUT}/${sample_name}.bam"

        # BAM by name sorting
        samtools sort -n "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted_by_name.bam"

        # BAM by coordinates sorting
        samtools sort "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted.bam"

         # BAM file indexation
        samtools index "${OUTPUT}/${sample_name}_sorted.bam"

        # Unsorted BAM and SAM files removal
        rm "${OUTPUT}/${sample_name}.sam" "${OUTPUT}/${sample_name}.bam"

    else
        echo "Error, skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "OK, all samples processed."
