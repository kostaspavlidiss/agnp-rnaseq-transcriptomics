# AgNP RNA-Seq Project
# Step 4: Volcano Plot

results <- read.csv(
  "results/DESeq2_all_genes.csv"
)

# Keep valid results
valid <- !is.na(results$padj) &
  is.finite(results$log2FoldChange)

results <- results[valid, ]

# Classify genes
results$Status <- "Below joint thresholds"

results$Status[
  results$padj < 0.05 &
  results$log2FoldChange > 1
] <- "Upregulated"

results$Status[
  results$padj < 0.05 &
  results$log2FoldChange < -1
] <- "Downregulated"

# Calculate negative log10 FDR
results$negLogFDR <- -log10(
  pmax(results$padj, .Machine$double.xmin)
)

# Colors
colors <- c(
  "Below joint thresholds" = "#BFC4CC",
  "Upregulated" = "#8260B0",
  "Downregulated" = "#3487A0"
)

# Save volcano plot
png(
  "figures/03_Volcano_Plot.png",
  width = 2000,
  height = 1500,
  res = 200
)

par(mar = c(5, 5, 4, 2))

plot(
  results$log2FoldChange,
  results$negLogFDR,
  pch = 16,
  cex = 0.5,
  col = adjustcolor(
    colors[results$Status],
    alpha.f = 0.65
  ),
  xlab = "log2 Fold Change (AgNP / Control)",
  ylab = "-log10 adjusted p-value",
  main = "Transcriptomic Response to Silver Nanoparticles"
)

# Significance thresholds
abline(
  v = c(-1, 1),
  lty = 2,
  col = "grey55"
)

abline(
  h = -log10(0.05),
  lty = 2,
  col = "grey55"
)

legend(
  "topright",
  legend = names(colors),
  col = colors,
  pch = 16,
  bty = "n",
  cex = 0.85
)

dev.off()

# Save annotated classification
write.csv(
  results,
  "results/Volcano_gene_classification.csv",
  row.names = FALSE
)

cat("\n=== VOLCANO PLOT SUMMARY ===\n")
print(table(results$Status))

cat("\nVolcano plot saved successfully!\n")
