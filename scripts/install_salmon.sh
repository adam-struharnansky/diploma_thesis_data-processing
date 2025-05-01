#!/bin/bash

# Create a new conda environment for Salmon
conda create -n salmon_env -y

# Conda enviroment activation
source activate salmon_env

# Salmon instalation from Bioconda
conda install -c bioconda -c conda-forge salmon=1.10.1 -y
