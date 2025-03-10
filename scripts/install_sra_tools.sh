#!/bin/bash

# Create a new conda environment
conda create -n sra_tools_env -y

# Activate the environment
source activate sra_tools_env

# Install SRA tools using the bioconda channel
conda install -c bioconda sra-tools -y
