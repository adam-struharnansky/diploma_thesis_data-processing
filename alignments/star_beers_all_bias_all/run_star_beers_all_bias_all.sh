#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate star_env


# Directories
INDEX="genetic_data/indexes/star_index_mus_musculus"
READS="genetic_data/reads/beers_all_bias"
OUTPUT="genetic_data/alignments/star_beers_all_bias_all"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT"

# Loop through each subdirectory
for sample_dir in "$READS"/*/; do
    # Extract sample name
    sample_name=$(basename "$sample_dir")

    # Define paired-end FASTQ files
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    echo "$sample_name"

    # Check if both R1 and R2 exist before running STAR
    if [[ -f "$r1" && -f "$r2" ]]; then
        # Run STAR alignment
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

        # Rename to fit convention
        mv "${OUTPUT}/${sample_name}_Aligned.out.bam" "${OUTPUT}/${sample_name}.bam"

        # Sort BAM by name
        samtools sort -n "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted_by_name.bam"

        # Sort BAM by position
        samtools sort "${OUTPUT}/${sample_name}.bam" -o "${OUTPUT}/${sample_name}_sorted.bam"

        # Index the sorted BAM file
        samtools index "${OUTPUT}/${sample_name}_sorted.bam"

        # Remove the original unsorted BAM file
        rm "${OUTPUT}/${sample_name}.bam"

        echo "Finished $sample_name. Results stored in $OUTPUT"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done

echo "All samples processed."
