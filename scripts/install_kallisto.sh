#!/bin/bash

# Create a new conda environment for Kallisto
conda create -n kallisto_env -y

# Activate the new environment
source activate kallisto_env

# Install Kallisto version 0.48.0 from Bioconda
conda install -c bioconda -c conda-forge kallisto=0.48.0 -y