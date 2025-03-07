#!/bin/bash

# Create a new conda environment
conda create -n htseq_env -y

# Activate the new environment
source activate htseq_env

# Install HTSeq using the bioconda channel
conda install -c bioconda htseq=2.0.3 -y