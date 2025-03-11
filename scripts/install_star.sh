#!/bin/bash

# Create a new conda environment with Python 3.8
conda create -n star_env -y

# Activate the new environment
source activate star_env

# Install STAR and samtools
conda install -c bioconda -c conda-forge star=2.7.10b samtools=1.9 -y