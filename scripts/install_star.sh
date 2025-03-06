#!/bin/bash

# Create a new conda environment (optional)
conda create -n star_env python=3.8 -y

# Activate the new environment
source activate star_aligner_env

# Install STAR using the bioconda channel
conda install -c bioconda star -y
