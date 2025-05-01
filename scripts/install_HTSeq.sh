#!/bin/bash

# Create a new conda environment
conda create -n htseq_env -y

# Conda enviroment activation
source activate htseq_env

# htseq-counts instalation using the bioconda channel
conda install -c bioconda -c conda-forge htseq=0.12.3 -y
