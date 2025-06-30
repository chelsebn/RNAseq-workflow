#!/bin/bash

# define directories
genomeDir="/mnt/p/RNAseq_data/star_index"
trimmedDir="/mnt/p/RNAseq_data/results/2_trimmed_output"
alignedDir="/mnt/p/RNAseq_data/results/4_aligned_sequences"
bamDir="${alignedDir}/aligned_bam"
logDir="${alignedDir}/aligned_logs"

# make output directories if missing
mkdir -p "$alignedDir" "$bamDir" "$logDir"

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
  echo "Aligning $sample"
  
  STAR \
    --genomeDir "$genomeDir" \
    --readFilesCommand zcat \
    --readFilesIn "${trimmedDir}/${sample}_R1_001_val_1.fq.gz" "${trimmedDir}/${sample}_R2_001_val_2.fq.gz" \
    --outFilterMismatchNmax 2 \
    --runThreadN 4 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix "${alignedDir}/${sample}_"

  # Move the BAM file to BAM folder
  mv -v "${alignedDir}/${sample}_Aligned.sortedByCoord.out.bam" "$bamDir"

  # Move the logs to the log folder
  mv -v "${alignedDir}/${sample}_Log.out" "$logDir"
  mv -v "${alignedDir}/${sample}_Log.final.out" "$logDir"
  mv -v "${alignedDir}/${sample}_Log.progress.out" "$logDir"
  mv -v "${aligned

