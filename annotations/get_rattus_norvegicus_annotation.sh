#!/bin/bash

# Output directory creation if nonexistence
mkdir -p ~/genetic_data/annotations/

# Annotation file download.
wget -P ~/genetic_data/annotations/ https://ftp.ensembl.org/pub/release-77/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_5.0.77.gtf.gz

echo "OK, file downloaded."
