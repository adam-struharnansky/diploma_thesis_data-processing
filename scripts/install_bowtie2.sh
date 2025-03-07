#!/bin/bash

# Create a new conda environment 
conda create -n bowtie2_env -c conda-forge -c bioconda -c defaults samtools=1.10 openssl=1.1.1 bowtie2=2.5.4 -y