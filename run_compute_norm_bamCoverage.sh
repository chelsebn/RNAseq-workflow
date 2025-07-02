#!/bin/bash

#SBATCH --job-name=compute_norm_bamCoverage    ## Name of the job.
#SBATCH -A  CHELSEBN 
#SBATCH -p standard              ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=4     ## number of cores the job needs
#SBATCH --error=compute_norm_bamCoverage-%J.err  ## error log file
#SBATCH --output=compute_norm_bamCoverage-%J.out ## output log file

source ~/miniconda3/etc/profile.d/conda.sh
#conda create -n deeptools_env -c bioconda -c conda-forge python=3.9 deeptools
conda activate deeptools_env

# move to BAM folder
cd /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam || {
  echo "BAM folder missing!"
  exit 1
}

for bamfile in *.bam
do
  sample=$(basename "$bamfile" .bam)
  echo "Converting $bamfile to bigWig..."
  bamCoverage \
    --bam "$bamfile" \
    -o "bamcoverage/${sample}.bw" \
    --normalizeUsing BPM \
    --numberOfProcessors 4
done
