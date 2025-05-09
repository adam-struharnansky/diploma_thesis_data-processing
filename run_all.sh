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

# Download the reads
./genetic_data/reads/seqcB/clean_seqcB_data.sh
./genetic_data/reads/seqcB/get_seqcB_data.sh
./genetic_data/reads/beers_all_bias/get_BEERS_reads.sh
./genetic_data/reads/rat/get_rat_data.sh
./genetic_data/reads/seqcA/clean_seqA_data.sh
./genetic_data/reads/seqcA/get_seqcA_data.sh


alignments  counts   indexes  reads  transcriptomes
annotations  genomes  outputs  