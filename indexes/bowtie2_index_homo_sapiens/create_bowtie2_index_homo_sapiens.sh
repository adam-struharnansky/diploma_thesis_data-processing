#!/bin/bash

# Define variables
THREADS=4
GENOME_DIR=genetic_data/indexes/bowtie2_index_homo_sapiens
GENOME_FASTA=genetic_data/genomes/GRCh37.p13.genome.fa.gz

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
bowtie2-build --threads $THREADS $GENOME_FASTA_UNZIPPED $GENOME_DIR/bowtie2_homo_sapiens_index

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
