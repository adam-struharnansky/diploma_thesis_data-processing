#!/bin/bash

# Create a new conda environment
conda create -n htseq_env -y

# Activate the new environment
source activate htseq_env

# Install htseq-counts using the bioconda channel
conda install -c bioconda -c conda-forge htseq=0.12.3 -y