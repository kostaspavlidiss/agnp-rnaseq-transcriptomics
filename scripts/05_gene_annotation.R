# AgNP RNA-seq Project
# Step 5: Gene Annotation with HGNC

# Load DESeq2 results
res <- read.csv(
  "results/DESeq2_all_genes.csv",
  stringsAsFactors = FALSE
)

# Load official HGNC annotations
hgnc <- read.delim(
  "data/hgnc_complete_set.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required <- c(
  "ensembl_gene_id", "symbol", "name", "status"
)

stopifnot(all(required %in% names(hgnc)))

# Retain approved mappings
hgnc <- hgnc[
  hgnc$status == "Approved" &
  !is.na(hgnc$ensembl_gene_id) &
  nzchar(hgnc$ensembl_gene_id),
]

mapping <- hgnc[, c(
  "ensembl_gene_id", "symbol", "name"
)]

# Exclude ambiguous mappings
duplicated_id <- duplicated(mapping$ensembl_gene_id) |
  duplicated(mapping$ensembl_gene_id, fromLast = TRUE)

mapping <- mapping[!duplicated_id, ]

# Match Ensembl IDs
clean_id <- sub("\\.[0-9]+$", "", res$Gene_ID)

idx <- match(
  clean_id,
  mapping$ensembl_gene_id
)

res$Gene_Symbol <- mapping$symbol[idx]
res$Gene_Name <- mapping$name[idx]

# Identify significant genes
sig <- !is.na(res$padj) &
  res$padj < 0.05 &
  abs(res$log2FoldChange) > 1

# Export annotated results
write.csv(
  res,
  "results/DESeq2_annotated_genes.csv",
  row.names = FALSE
)

write.csv(
  res[sig, ],
  "results/Annotated_significant_genes.csv",
  row.names = FALSE
)

cat("\n=== ANNOTATION SUMMARY ===\n")

cat("Successfully annotated:",
    sum(!is.na(res$Gene_Symbol)),
    "of", nrow(res), "\n")

cat("Significant genes annotated:",
    sum(sig & !is.na(res$Gene_Symbol)),
    "of", sum(sig), "\n")

cat("\n=== TOP 15 DIFFERENTIAL GENES ===\n")

print(
  head(
    res[sig, c(
      "Gene_ID",
      "Gene_Symbol",
      "Gene_Name",
      "log2FoldChange",
      "padj"
    )],
    15
  ),
  row.names = FALSE
)

cat("\nGene annotation completed!\n")
