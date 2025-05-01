#!/bin/bash

# Exit immediately if any command fails
set -e

# Miniconda installation path definition
MINICONDA_DIR="$HOME/miniconda3"

# Existing Miniconda installation removal if present
if [ -d "$MINICONDA_DIR" ]; then
    rm -rf "$MINICONDA_DIR"
fi

# Miniconda installer download
echo "Downloading Miniconda..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh

# Installer run
echo "Installing Miniconda..."
bash Miniconda3.sh -b -p "$MINICONDA_DIR"

# Conda initialization
echo "Initializing Conda..."
$MINICONDA_DIR/bin/conda init bash

# Changes activation
source ~/.bashrc

# Bioconda channels addition
echo "Adding Bioconda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Installation verification
echo "Conda version:"
conda --version

# Installer cleanup
rm -f Miniconda3.sh

echo "Setup complete!"
