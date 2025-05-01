#!/bin/bash

# Create a new conda environment
conda create -n feature_counts_env -y

# Conda enviroment activation
source activate feature_counts_env

# featureCounts instalation using the bioconda channel)
conda install -c bioconda subread -y
