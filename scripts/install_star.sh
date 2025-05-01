#!/bin/bash

# Create a new conda environment with Python 3.8
conda create -n star_env -y

# Conda enviroment activation
source activate star_env

# STAR and samtools instalation
conda install -c bioconda -c conda-forge star=2.7.10b samtools=1.9 -y
