#!/bin/bash

# Conda enviroment activation
source miniconda3/etc/profile.d/conda.sh
conda activate sra_tools_env

# File paths setting
BASE_DIR=genetic_data/reads/rat

# Output directory creation if nonexistence
mkdir -p "$BASE_DIR"

# Parameter setting: format => RNATYPE|SRR1[,SRR2]|ORGAN|SEX|AGE
samples=(
  "mirna|SRR14265300|brain|female|006"
  "rna|SRR1169967,SRR1169968|brain|female|006"
  "mirna|SRR14265316|brain|male|006"
  "rna|SRR1170001,SRR1170002|brain|male|006"
  "mirna|SRR14265496|spleen|female|021"
  "rna|SRR1170371,SRR1170372|spleen|female|021"
)

# Loop through each sample
for sample in "${samples[@]}"; do
  IFS="|" read -r TYPE SRRS ORGAN SEX AGE <<< "$sample"
  OUTPUT_NAME="rn_${ORGAN}_${TYPE}_${SEX}_${AGE}"
  OUTPUT_PATH="${BASE_DIR}/${OUTPUT_NAME}.fastq.gz"

  # Simple miRNAs saving
  if [[ "$SRRS" != *","* ]]; then
    SRR=$SRRS
    fasterq-dump "$SRR" -e 4 -O "$BASE_DIR"
    gzip "${BASE_DIR}/${SRR}.fastq"
    mv "${BASE_DIR}/${SRR}.fastq.gz" "$OUTPUT_PATH"

  # RNA concatenation and saving
  else
    IFS="," read -r SRR1 SRR2 <<< "$SRRS"
    fasterq-dump "$SRR1" -e 4 -O "$BASE_DIR"
    fasterq-dump "$SRR2" -e 4 -O "$BASE_DIR"
    cat "${BASE_DIR}/${SRR1}.fastq" "${BASE_DIR}/${SRR2}.fastq" > "${BASE_DIR}/tmp_${OUTPUT_NAME}.fastq"
    gzip "${BASE_DIR}/tmp_${OUTPUT_NAME}.fastq"
    mv "${BASE_DIR}/tmp_${OUTPUT_NAME}.fastq.gz" "$OUTPUT_PATH"
    rm "${BASE_DIR}/${SRR1}.fastq" "${BASE_DIR}/${SRR2}.fastq"
  fi
done

echo "OK, all files downaloaded!"
