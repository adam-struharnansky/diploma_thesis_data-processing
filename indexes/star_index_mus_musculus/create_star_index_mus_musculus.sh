#!/bin/bash

# Activate conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate star_env

# Define variables
THREADS=4
MAX_RAM=51539607552  # 48GB in bytes
GENOME_DIR=genetic_data/indexes/star_index_mus_musculus
GENOME_FASTA=genetic_data/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
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
     --limitGenomeGenerateRAM $MAX_RAM

# Re-zip files if needed
rezip_if_needed $GENOME_FASTA
rezip_if_needed $ANNOTATION_GTF
