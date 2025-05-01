#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate star_env

# File paths setting
INDEX="genetic_data/indexes/star_index_homo_sapiens"
READS="genetic_data/reads/seqcA"
OUTPUT="genetic_data/alignments/star_seqcA_all"

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
        # STAR alignment run
        STAR --runThreadN 8 \
             --genomeDir "$INDEX" \
             --readFilesIn "$r1" "$r2" \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix "${OUTPUT}/${sample_name}_" \
             --outFilterMultimapNmax 10 \
             --outSAMattributes NH HI AS nM \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000

        # Renaming to fit convention
        mv "${OUTPUT}/${sample_name}_Aligned.out.bam" "${OUTPUT}/${sample_name}.bam"

        # BAM by name sorting
        samtools sort -n "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted_by_name.bam"

        # BAM by coordinates sorting
        samtools sort "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted.bam"

         # BAM file indexation
        samtools index "${OUTPUT}/${sample_name}_sorted.bam"

        # Unsorted BAM file removal
        rm "${OUTPUT}/${sample_name}.bam"

    else
        echo "Error, skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "OK, all samples processed."
