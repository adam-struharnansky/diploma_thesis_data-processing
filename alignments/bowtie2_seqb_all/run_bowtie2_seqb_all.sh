#!/bin/bash

# Directories
INDEX="~/genetic_data/indexes/bowtie2_index_homo_sapiens/bowtie2_homo_sapiens_index"
READS="~/genetic_data/reads/seqc_b"
OUTPUT="~/genetic_data/alignments/bowtie2_seqb_all"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT"

# Loop through each subdirectory (sample)
for sample_dir in "$READS"/*/; do
    # Extract sample name
    sample_name=$(basename "$sample_dir")

    # Define paired-end FASTQ files
    r1="${sample_dir}/${sample_name}_1.fastq.gz"
    r2="${sample_dir}/${sample_name}_2.fastq.gz"

    # Check if both R1 and R2 exist before running bowtie
    if [[ -f "$r1" && -f "$r2" ]]; then
        # Output directory for this sample
        sample_out="${OUTPUT}/${sample_name}"
        mkdir -p "$sample_out"

        echo "Processing $sample_name..."

        # Run Bowtie2 alignment
        bowtie2 -x "$INDEX" -1 "$r1" -2 "$r2" --very-sensitive -k 10 -S "${sample_out}/alignment.sam"

        echo "Finished $sample_name. Results stored in $sample_out"
    else
        echo "Skipping $sample_name: Missing R1 or R2 file."
    fi
done
echo "All samples processed."
