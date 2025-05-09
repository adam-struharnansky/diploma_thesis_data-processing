#!/bin/bash

cd ..

# conda installation
./genetic_data/scripts//install_conda.sh

# Tools installation
./genetic_data/scripts//install_HTSeq.sh
./genetic_data/scripts//install_bowtie2.sh
./genetic_data/scripts//install_fastp.sh
./genetic_data/scripts//install_featureCounts.sh
./genetic_data/scripts//install_hisat2.sh
./genetic_data/scripts//install_kallisto.sh
./genetic_data/scripts//install_salmon.sh
./genetic_data/scripts//install_sra_tools.sh
./genetic_data/scripts//install_star.sh
./genetic_data/scripts//set_bilattice_count.sh

#
