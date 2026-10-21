# AgNP RNA-seq Project
# Step 6: Biological Signature Analysis

suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS("results/DESeq2_object.rds")

res <- read.csv(
  "results/DESeq2_annotated_genes.csv",
  stringsAsFactors = FALSE
)

# Predefined panel: epithelial, ECM and EMT markers
panel <- c(
  "TP63", "KRT6A", "DSP", "FGFR2",
  "COL17A1", "FAT2", "S100A14",
  "TIMP3", "ITGA11", "BGN",
  "DPYSL3", "IGFBP4",
  "CDH1", "CDH2", "VIM", "FN1",
  "MMP2", "ZEB2", "TGFB1", "NOTCH3"
)

# Match gene symbols to Ensembl identifiers
idx <- match(panel, res$Gene_Symbol)
selected <- res[idx[!is.na(idx)], , drop = FALSE]

cat("\n=== MISSING PANEL GENES ===\n")
print(setdiff(panel, selected$Gene_Symbol))

# Retrieve DESeq2-normalized counts
normalized <- counts(dds, normalized = TRUE)

selected <- selected[
  selected$Gene_ID %in% rownames(normalized),
  ,
  drop = FALSE
]

mat <- log2(
  normalized[selected$Gene_ID, , drop = FALSE] + 1
)

rownames(mat) <- selected$Gene_Symbol

# Remove genes with no variation
sd_values <- apply(mat, 1, sd)
keep <- is.finite(sd_values) & sd_values > 0

mat <- mat[keep, , drop = FALSE]
selected <- selected[keep, , drop = FALSE]

stopifnot(nrow(mat) >= 2)

# Gene-wise Z-scores
z <- t(scale(t(mat)))
z <- pmax(pmin(z, 2), -2)

# Export selected genes
selected$Significant <- !is.na(selected$padj) &
  selected$padj < 0.05 &
  abs(selected$log2FoldChange) > 1

write.csv(
  selected,
  "results/Biological_signature_genes.csv",
  row.names = FALSE
)

cat("\n=== BIOLOGICAL SIGNATURE ===\n")
print(
  selected[, c(
    "Gene_Symbol",
    "log2FoldChange",
    "padj",
    "Significant"
  )],
  row.names = FALSE
)

# Group colors
group_colors <- c(
  rep("#267D95", 3),
  rep("#8661B5", 3)
)

# Generate heatmap
png(
  "figures/04_Biological_Signature.png",
  width = 1800,
  height = 2000,
  res = 180
)

heatmap(
  z,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  ColSideColors = group_colors,
  col = colorRampPalette(
    c("#267D95", "white", "#8661B5")
  )(101),
  breaks = seq(-2, 2, length.out = 102),
  margins = c(9, 11),
  labCol = c(
    "C1", "C2", "C3",
    "Ag1", "Ag2", "Ag3"
  ),
  main = "AgNP Exposure | Epithelial & ECM Signature",
  cexRow = 0.85,
  cexCol = 0.9
)

dev.off()

cat("\nBiological signature heatmap saved!\n")
