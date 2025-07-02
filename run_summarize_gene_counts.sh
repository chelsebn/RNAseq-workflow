#!/bin/bash

# Move to the folder containing the BAM files
cd results/3_aligned_sequences/aligned_bam || {
  echo "Error: cannot cd into results/3_aligned_sequences/aligned_bam"
  exit 1
}

# Store list of BAM files (space-separated)
dirlist=$(ls -t ./*.bam | tr '\n' ' ')
echo "Processing the following BAM files:"
echo "$dirlist"

# Create output directory if it does not exist
mkdir -p ../../results/4_final_counts

# Run featureCounts
featureCounts \
  -a ../../annotation/* \
  -o ../../results/4_final_counts/final_counts.txt \
  -g gene_name \
  -T 4 \
  -M \
  --fraction \
  -p \
  $dirlist
