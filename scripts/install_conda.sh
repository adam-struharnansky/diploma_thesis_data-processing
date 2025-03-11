#!/bin/bash

# Exit immediately if any command fails
set -e

# Define Miniconda installation path
MINICONDA_DIR="$HOME/miniconda3"

# Remove existing Miniconda installation if present
if [ -d "$MINICONDA_DIR" ]; then
    echo "Removing existing Miniconda installation..."
    rm -rf "$MINICONDA_DIR"
fi

# Download the latest Miniconda installer
echo "Downloading Miniconda..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh

# Run the installer
echo "Installing Miniconda..."
bash Miniconda3.sh -b -p "$MINICONDA_DIR"

# Initialize Conda
echo "Initializing Conda..."
$MINICONDA_DIR/bin/conda init bash

# Activate changes
echo "Activating Conda..."
source ~/.bashrc

# Add Bioconda channels
echo "Adding Bioconda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Enable automatic activation of base
conda config --set auto_activate_base true

# Verify installation
echo "Conda installation complete. Checking version..."
conda --version

# Cleanup installer
rm -f Miniconda3.sh

echo "Restarting shell to activate conda"
exec bash
