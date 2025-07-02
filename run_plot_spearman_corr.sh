#!/bin/bash

#SBATCH --job-name=plot_spearman_corr    ## Name of the job.
#SBATCH -A  CHELSEBN 
#SBATCH -p standard              ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=4     ## number of cores the job needs
#SBATCH --error=plot_spearman_corr-%J.err  ## error log file
#SBATCH --output=plot_spearman_corr-%J.out ## output log file

# load conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate RNASeq_data

plotCorrelation -in /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam/bamcoverage/bw-summary.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam/bamcoverage/heatmap_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix /pub/chelsebn/BakerLab/RNAseq_data/results/3_aligned_sequences/aligned_bam/bamcoverage/SpearmanCorr_readCounts.tab
