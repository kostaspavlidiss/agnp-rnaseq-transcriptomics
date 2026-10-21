# AgNP RNA-Seq Mini Project
# Step 2: Sample Correlation

# Import counts
data <- read.delim(
  "data/count_table_Ag10.txt",
  check.names = FALSE
)

counts <- as.matrix(data[, -1])
rownames(counts) <- data[[1]]

# Library sizes
lib_sizes <- colSums(counts)

# Filter low-count features
keep <- rowSums(counts >= 10) >= 3

filtered <- counts[keep, , drop = FALSE]

# CPM normalization
cpm <- sweep(
  filtered,
  2,
  lib_sizes,
  "/"
) * 1e6

log_cpm <- log2(cpm + 1)

# Pearson correlations
cor_matrix <- cor(
  log_cpm,
  method = "pearson"
)

cat("\n=== SAMPLE CORRELATIONS ===\n")
print(round(cor_matrix, 3))

# Save numerical results
write.csv(
  cor_matrix,
  "results/sample_correlations.csv"
)

# Hierarchical clustering
hc <- hclust(
  as.dist(1 - cor_matrix),
  method = "average"
)

# Generate heatmap
png(
  "figures/02_Correlation_Heatmap.png",
  width = 1600,
  height = 1400,
  res = 180
)

heatmap(
  cor_matrix,
  Rowv = as.dendrogram(hc),
  Colv = as.dendrogram(hc),
  symm = TRUE,
  scale = "none",
  col = colorRampPalette(
    c("#F7FCFD", "#7FBBC7", "#12475E")
  )(100),
  margins = c(9, 9),
  main = "RNA-Seq Sample Correlation"
)

dev.off()

cat("\nCorrelation heatmap saved!\n")
