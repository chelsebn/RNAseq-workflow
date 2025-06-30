#!/bin/bash

# create results directory
mkdir -p /mnt/p/RNAseq_data/results/2_trimmed_output

# list of sample names (no _R1/_R2 suffix)
samples=(
  68A 68B 68C
  70A_T 70B_T 70C_T
  77A 77B 77C
  82A_T 82B_T 82C_T
  F14A F14B F14C
  FA_T FB_T FC_T
)

# loop over samples
for sample in "${samples[@]}"
do
  echo "Processing $sample"
  
  trim_galore \
    --paired \
    --fastqc \
    --illumina \
    -j 6 \
    --trim-n \
    --length 25 \
    --output_dir /mnt/p/RNAseq_data/results/2_trimmed_output \
    /mnt/p/RNAseq_data/input/${sample}_R1.fastq.gz /mnt/p/RNAseq_data/input/${sample}_R2.fastq.gz
done
