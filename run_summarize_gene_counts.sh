#!/bin/bash

#SBATCH --job-name=summarize_gene_counts    ## Name of the job.
#SBATCH -A  CHELSEBN 
#SBATCH -p standard              ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=4     ## number of cores the job needs
#SBATCH --error=summarize_gene_counts-%J.err  ## error log file
#SBATCH --output=summarize_gene_counts-%J.out ## output log file

source ~/miniconda3/etc/profile.d/conda.sh
conda activate RNAseq_data

# Move to the folder containing the BAM files
cd /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam || {
  echo "Error: cannot cd into /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam"
  exit 1
}

# Store list of BAM files (space-separated)
dirlist=$(ls -t ./*.bam | tr '\n' ' ')
echo "Processing the following BAM files:"
echo "$dirlist"

# Run featureCounts
featureCounts \
  -a /pub/chelsebn/BakerLab/RNAseq_data/annotation/* \
  -o /pub/chelsebn/BakerLab/RNAseq_data/results/4_final_counts/final_counts.txt \
  -g gene_id \
  -T 4 \
  -M \
  --fraction \
  -p \
  $dirlist
