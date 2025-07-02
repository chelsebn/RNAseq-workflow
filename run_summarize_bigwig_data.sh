#!/bin/bash

#SBATCH --job-name=summarize_bigwig_data    ## Name of the job.
#SBATCH -A  CHELSEBN 
#SBATCH -p standard              ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=4     ## number of cores the job needs
#SBATCH --error=summarize_bigwig_data-%J.err  ## error log file
#SBATCH --output=summarize_bigwig_data-%J.out ## output log file

# load conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools_env

# move to bigWig folder
cd /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam/bamcoverage || {
  echo "ERROR: bamcoverage folder not found"
  exit 1
}

# collect all bigwig files
bigwigs=$(ls *.bw | tr '\n' ' ')

# run multiBigwigSummary
multiBigwigSummary bins \
  -b $bigwigs \
  -o bw-summary.npz

echo "multiBigwigSummary completed successfully in bamcoverage."
