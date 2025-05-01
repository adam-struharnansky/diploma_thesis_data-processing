#!/bin/bash

# Create a new conda environment
conda create -n sra_tools_env -y

# Conda enviroment activation
source activate sra_tools_env

# SRA tools instalation using the bioconda channel
conda install -c bioconda sra-tools -y
