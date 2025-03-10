#!/bin/bash

# Define variables
THREADS=4
GENOME_DIR=~/genetic_data/indexes/star_index_homo_sapiens
GENOME_FASTA=~/genetic_data/genomes/GRCh37.p13.genome.fa.gz
ANNOTATION_GTF=~/genetic_data/annotations/gencode.v19.annotation.gtf.gz

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
ANNOTATION_GTF_UNZIPPED=$(unzip_if_needed $ANNOTATION_GTF)

# Run STAR command
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $GENOME_FASTA_UNZIPPED \
     --sjdbGTFfile $ANNOTATION_GTF_UNZIPPED

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
rezip_if_needed $ANNOTATION_GTF
