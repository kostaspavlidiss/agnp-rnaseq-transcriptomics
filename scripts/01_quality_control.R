# AgNP RNA-Seq Mini Project
# Step 1: Quality Control and PCA

# Load count table
data <- read.delim(
  "data/count_table_Ag10.txt",
  check.names = FALSE
)

genes <- data[[1]]
counts <- as.matrix(data[, -1])

# Validate counts
stopifnot(
  is.numeric(counts),
  !anyNA(counts),
  all(counts >= 0),
  all(counts == floor(counts)),
  !anyDuplicated(genes)
)

rownames(counts) <- genes

# Define biological conditions
samples <- colnames(counts)

condition <- c(
  rep("Control", 3),
  rep("AgNP", 3)
)

stopifnot(
  identical(samples, paste0("P865_20", 1:6))
)

# Basic quality control
library_sizes <- colSums(counts)
expressed <- colSums(counts > 0)

qc <- data.frame(
  Sample = samples,
  Condition = condition,
  Total_Counts = library_sizes,
  Detected_Features = expressed
)

cat("\n=== SAMPLE QUALITY CONTROL ===\n")
print(qc, row.names = FALSE)

write.csv(
  qc,
  "results/sample_qc.csv",
  row.names = FALSE
)

# Filter low-count features
keep <- rowSums(counts >= 10) >= 3

cat("\nFeatures before filtering:",
    nrow(counts), "\n")
cat("Features after filtering:",
    sum(keep), "\n")

# Counts per million (CPM)
filtered <- counts[keep, , drop = FALSE]

stopifnot(all(library_sizes > 0), nrow(filtered) > 1)

cpm <- sweep(
  filtered, 2, library_sizes, "/"
) * 1e6

log_cpm <- log2(cpm + 1)

# Principal Component Analysis
pca <- prcomp(
  t(log_cpm),
  center = TRUE,
  scale. = FALSE
)

variance <- 100 * pca$sdev^2 /
  sum(pca$sdev^2)

cat("\n=== PCA EXPLAINED VARIANCE ===\n")
print(round(variance[1:2], 2))

# Save PCA coordinates
coordinates <- data.frame(
  Sample = samples,
  Condition = condition,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

write.csv(
  coordinates,
  "results/PCA_coordinates.csv",
  row.names = FALSE
)

# Publication-style PCA figure
png(
  "figures/01_PCA_AgNP.png",
  width = 1800,
  height = 1300,
  res = 180
)

colors <- ifelse(
  condition == "Control",
  "#267D95",
  "#8661B5"
)

plot(
  pca$x[, 1:2],
  type = "n",
  main = "RNA-Seq | Silver Nanoparticle Exposure",
  xlab = paste0("PC1 (", round(variance[1], 1), "%)"),
  ylab = paste0("PC2 (", round(variance[2], 1), "%)")
)

abline(h = 0, v = 0, col = "grey85", lty = 3)

points(
  pca$x[, 1:2],
  col = colors,
  pch = 19,
  cex = 2
)

text(
  pca$x[, 1:2],
  labels = samples,
  pos = 3,
  cex = 0.75
)

legend(
  "topright",
  legend = c("Control", "AgNP-treated"),
  col = c("#267D95", "#8661B5"),
  pch = 19,
  bty = "n"
)

dev.off()

cat("\nPCA figure saved successfully!\n")
