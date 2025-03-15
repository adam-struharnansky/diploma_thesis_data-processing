#!/bin/bash

# Define variables
THREADS=4
GENOME_DIR=genetic_data/indexes/bowtie2_index_mus_musculus
GENOME_FASTA=genetic_data/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

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
        gzip ${file%.gz}
    fi
}

# Unzip files if needed
GENOME_FASTA_UNZIPPED=$(unzip_if_needed $GENOME_FASTA)

# Create output directory if not exists
mkdir -p $GENOME_DIR

# Run Bowtie2 indexing
bowtie2-build --threads $THREADS $GENOME_FASTA_UNZIPPED $GENOME_DIR/bowtie2_mus_musculus_index

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
