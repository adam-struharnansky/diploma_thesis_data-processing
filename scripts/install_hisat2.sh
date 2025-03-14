#!/bin/bash

# Create a new conda environment with Python 3.8
conda create -n hisat2_env -y

# Activate the new environment
source activate hisat2_env

# Install STAR and samtools
conda install -c bioconda -c conda-forge hisat2=2.2.1 samtools=1.9 -y