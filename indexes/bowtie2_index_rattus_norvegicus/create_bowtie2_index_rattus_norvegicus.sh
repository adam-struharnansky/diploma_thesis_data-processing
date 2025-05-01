#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate bowtie2_env

# File paths setting
GENOME_DIR=genetic_data/indexes/bowtie2_index_rattus_norvegicus
GENOME_FASTA=genetic_data/genomes/Rattus_norvegicus.Rnor_5.0.dna.toplevel.fa.gz

# Output directory creation if nonexistence
mkdir -p "$GENOME_DIR"

# Parameter setting
THREADS=4

# Function to unzip files if they are compressed
unzip_if_needed() {
    local file=$1
    if [[ $file == *.gz ]]; then
        gunzip -k $file
        echo "${file%.gz}"
    else
        echo $file
    fi
}

# Function to re-zip files if they were originally compressed
rezip_if_needed() {
    local file=$1
    if [[ $file == *.gz ]]; then
        local unzipped_file="${file%.gz}"
        gzip -f $unzipped_file
        rm -f $unzipped_file
    fi
}

# Files unzip
GENOME_FASTA_UNZIPPED=$(unzip_if_needed $GENOME_FASTA)

# Output directory creation if nonexistence
mkdir -p $GENOME_DIR

# Run Bowtie2 indexing
bowtie2-build --threads $THREADS $GENOME_FASTA_UNZIPPED $GENOME_DIR/bowtie2_rattus_norvegicus_index

# Files zip
rezip_if_needed $GENOME_FASTA

echo "OK, index created!"
