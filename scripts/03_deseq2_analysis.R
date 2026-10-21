# AgNP RNA-Seq Mini Project
# Step 3: Differential Gene Expression
# Dataset: E-MTAB-6321

suppressPackageStartupMessages(library(DESeq2))

# 1. Load raw gene counts
data <- read.delim(
  "data/count_table_Ag10.txt",
  check.names = FALSE
)

genes <- data[[1]]
counts_matrix <- as.matrix(data[, -1])

stopifnot(
  is.numeric(counts_matrix),
  !anyNA(counts_matrix),
  all(counts_matrix >= 0),
  all(counts_matrix == round(counts_matrix)),
  !anyDuplicated(genes)
)

rownames(counts_matrix) <- genes

# 2. Experimental design
expected <- paste0("P865_20", 1:6)
stopifnot(identical(colnames(counts_matrix), expected))

condition <- factor(
  c(rep("Control", 3), rep("AgNP", 3)),
  levels = c("Control", "AgNP")
)

sample_info <- data.frame(
  condition = condition,
  row.names = colnames(counts_matrix)
)

# 3. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info,
  design = ~ condition
)

# 4. Remove low-count features
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

cat("\nFeatures retained:", nrow(dds), "\n")

# 5. Differential expression
dds <- DESeq(dds, parallel = FALSE, quiet = TRUE)

res <- results(
  dds,
  contrast = c("condition", "AgNP", "Control"),
  alpha = 0.05
)

# 6. Export results
results_df <- data.frame(
  Gene_ID = rownames(res),
  as.data.frame(res)
)

results_df <- results_df[
  order(results_df$padj, na.last = TRUE),
]

write.csv(
  results_df,
  "results/DESeq2_all_genes.csv",
  row.names = FALSE
)

# 7. Identify significant genes
valid <- !is.na(results_df$padj) &
         !is.na(results_df$log2FoldChange)

significant <- valid &
  results_df$padj < 0.05 &
  abs(results_df$log2FoldChange) > 1

up <- significant &
  results_df$log2FoldChange > 1

down <- significant &
  results_df$log2FoldChange < -1

write.csv(
  results_df[significant, ],
  "results/DESeq2_significant_genes.csv",
  row.names = FALSE
)

# 8. Print results
cat("\n=== SIZE FACTORS ===\n")
print(round(sizeFactors(dds), 3))

cat("\n=== DIFFERENTIAL EXPRESSION RESULTS ===\n")
cat("Genes tested:", sum(!is.na(res$pvalue)), "\n")
cat("FDR < 0.05:",
    sum(valid & results_df$padj < 0.05), "\n")
cat("Significant with |log2FC| > 1:",
    sum(significant), "\n")
cat("Upregulated:", sum(up), "\n")
cat("Downregulated:", sum(down), "\n")

cat("\n=== TOP 10 GENES ===\n")
print(
  head(
    results_df[, c(
      "Gene_ID",
      "log2FoldChange",
      "pvalue",
      "padj"
    )],
    10
  ),
  row.names = FALSE
)

# Save reproducibility information
saveRDS(dds, "results/DESeq2_object.rds")
capture.output(
  sessionInfo(),
  file = "results/R_session_info.txt"
)

cat("\nAnalysis completed successfully!\n")
