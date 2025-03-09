#!/bin/bash

# Create a new conda environment for Salmon
conda create -n salmon_env -y

# Activate the new environment
source activate salmon_env

# Install Salmon from Bioconda
conda install -c bioconda -c conda-forge salmon=1.10.1 -y