# install needed plotting packages if not yet present
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(ggplot2)
library(reshape2)

# read your CSV
counts_df <- read.csv("/mnt/p/Data/GeneCounts_R1_R2_together/RNAData_R1_R2_together_filtered.csv", header=TRUE)

# counts matrix (columns 7 onward)
counts_mat <- counts_df[, 7:ncol(counts_df)]
rownames(counts_mat) <- counts_df$Gene

# gene lengths in base pairs
gene_length <- counts_df$Length
names(gene_length) <- counts_df$Gene

# library sizes (total mapped reads per sample)
lib_size <- colSums(counts_mat)

# manual FPKM calculation
calculateFPKM <- function(counts, gene_length, libsize) {
  fpkm <- sweep(counts, 2, libsize, FUN="/") * 1e9
  fpkm <- sweep(fpkm, 1, gene_length, FUN="/")
  return(fpkm)
}

fpkm <- calculateFPKM(counts_mat, gene_length, lib_size)

# reshape for ggplot
fpkm_df <- as.data.frame(fpkm)
fpkm_df$Gene <- rownames(fpkm_df)

fpkm_long <- melt(fpkm_df,
                  id.vars = "Gene",
                  variable.name = "Sample",
                  value.name = "FPKM")

# apply pseudocount
fpkm_long$FPKM <- fpkm_long$FPKM + 0.0001

# plot
plot <- ggplot(fpkm_long, aes(x=FPKM)) +
  geom_histogram(bins=50, color="black", fill="gray") +
  scale_x_log10() +
  facet_wrap(~ Sample, scales="free_y") +
  labs(title="FPKM Distribution per Sample (with pseudocount 0.0001)",
       x="FPKM (log10 scale)",
       y="Number of Genes") +
  theme_bw()

# display the plot
print(plot)

# save the plot
ggsave("/mnt/p/fpkm_histograms.png", plot=plot, width=12, height=8, dpi=300)

# save the FPKM table
write.csv(fpkm, "/mnt/p/fpkm_results_pseudo_zero_count.csv")
