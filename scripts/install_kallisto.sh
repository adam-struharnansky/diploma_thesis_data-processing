#!/bin/bash

# Create a new conda environment for Kallisto
conda create -n kallisto_env -y

# Conda enviroment activation
source activate kallisto_env

# Kallisto instalation from Bioconda
conda install -c bioconda -c conda-forge kallisto=0.48.0 -y
