#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate star_env

# File paths setting
GENOME_DIR=genetic_data/indexes/star_index_rattus_norvegicus
GENOME_FASTA=genetic_data/genomes/Rattus_norvegicus.Rnor_5.0.dna.toplevel.fa.gz
ANNOTATION_GTF=genetic_data/annotations/Rattus_norvegicus.Rnor_5.0.77.gtf.gz

# Output directory creation if nonexistence
mkdir -p "$GENOME_DIR"

# Parameter setting
THREADS=4
MAX_RAM=51539607552  # 48GB in bytes

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
ANNOTATION_GTF_UNZIPPED=$(unzip_if_needed $ANNOTATION_GTF)

# Create genome index directory if it doesn't exist
mkdir -p $GENOME_DIR

# STAR indexing run
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $GENOME_FASTA_UNZIPPED \
     --sjdbGTFfile $ANNOTATION_GTF_UNZIPPED \
     --limitGenomeGenerateRAM $MAX_RAM

# Files zip
rezip_if_needed $GENOME_FASTA
rezip_if_needed $ANNOTATION_GTF

echo "OK, index created!"
