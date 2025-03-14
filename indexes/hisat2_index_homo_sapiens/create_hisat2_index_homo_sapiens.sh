#!/bin/bash

# Define variables
THREADS=4
GENOME_DIR=~/genetic_data/indexes/hisat2_index_homo_sapiens
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

# Create output directory if not exists
mkdir -p $GENOME_DIR

# Extract splice sites and exons from the GTF
hisat2_extract_splice_sites.py $ANNOTATION_GTF_UNZIPPED > $GENOME_DIR/splice_sites.txt
hisat2_extract_exons.py $ANNOTATION_GTF_UNZIPPED > $GENOME_DIR/exons.txt

# Build the HISAT2 index
hisat2-build --ss $GENOME_DIR/splice_sites.txt \
             --exon $GENOME_DIR/exons.txt \
             $GENOME_FASTA_UNZIPPED \
             $GENOME_DIR/genome_index \
             -p $THREADS

# Remove intermediate files to keep things clean
rm $GENOME_DIR/splice_sites.txt
rm $GENOME_DIR/exons.txt

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
rezip_if_needed $ANNOTATION_GTF
