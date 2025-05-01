#!/bin/bash

# Create a new conda environment with Python 3.10
conda create -n bilattice_env python=3.10 -y

# Conda enviroment activation
source activate bilattice_env

# Dependencies installation
conda install -c conda-forge -c bioconda \
    bottleneck=1.4.2 \
    intervaltree=3.1.0 \
    numexpr=2.10.1 \
    numpy=2.0.1 \
    numpy-base=2.0.1 \
    pandas=2.2.3 \
    python-dateutil=2.9.0post0 \
    python-tzdata=2023.3 \
    pytz=2024.1 \
    sortedcontainers=2.4.0 \
    pysam=0.23.0 \
    -y

# Installation confirmation
echo "Installed packages in bilattice_env:"
conda list -n bilattice_env
