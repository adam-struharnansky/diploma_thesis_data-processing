#!/bin/bash

# Activate conda enviroment
source miniconda3/etc/profile.d/conda.sh
conda activate hisat2_env

# Limit the RAM usage
ulimit -v 50000000  

# Define variables
THREADS=4
GENOME_DIR=genetic_data/indexes/hisat2_index_mus_musculus
GENOME_FASTA=genetic_data/genomes/Mus_musculus.GRCh38.dna.primary_assembly.fa.gz
ANNOTATION_GTF=genetic_data/annotations/Mus_musculus.GRCm38.102.gtf.gz

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
             -p $THREADS \
             --offrate 3


# Remove intermediate files
rm $GENOME_DIR/splice_sites.txt
rm $GENOME_DIR/exons.txt

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
rezip_if_needed $ANNOTATION_GTF
