#!/bin/bash

# Create a new conda environment (optional)
conda create -n feature_counts_env -y

# Activate the new environment
source activate feature_counts_env

# Install featureCounts using the bioconda channel (it is part of Subread package)
conda install -c bioconda subread -y
