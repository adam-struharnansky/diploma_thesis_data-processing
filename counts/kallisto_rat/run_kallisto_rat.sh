#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate kallisto_env

# File paths setting
INDEX="genetic_data/indexes/kallisto_index_rattus_norvegicus/Rattus_norvegicus.kallisto.idx"
READS="genetic_data/reads/rat"
OUTPUT="genetic_data/counts/kallisto_rat"

# Output directory creation if nonexistence
mkdir -p "$OUTPUT"

# Loop through each FASTQ file that does not contain 'mirna'
for fastq_file in "$READS"/*.fastq.gz; do
    filename=$(basename "$fastq_file")

    # Skip files with 'mirna' in the name
    if [[ "$filename" == *mirna* ]]; then
        continue
    fi

    # Remove file extension for sample name
    sample_name="${filename%.fastq.gz}"
    sample_out="${OUTPUT}/${sample_name}"

    mkdir -p "$sample_out"

    # kallisto quant for single-end data (adjust fragment length & sd as needed)
    kallisto quant -i "$INDEX" -o "$sample_out" --single -l 200 -s 30 --fr-stranded --threads=8 "$fastq_file"

    echo "Finished processing $sample_name"
done

echo "OK, all files processed!"
