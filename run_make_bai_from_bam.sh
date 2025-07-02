#!/bin/bash

#SBATCH --job-name=make_bai_from_bam
#SBATCH -A CHELSEBN
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --output=index_bams-%j.out
#SBATCH --error=index_bams-%j.err

# load conda environment with samtools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate RNAseq_data

# move to your BAM directory
cd /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam || {
  echo "ERROR: BAM folder not found."
  exit 1
}

# loop through each BAM file
for bam in *.bam
do
  echo "Indexing $bam ..."
  samtools index "$bam"
done

echo "All BAM files indexed successfully."
