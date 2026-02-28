# RNA-seq differential expression (DESeq2) - OVERALL MODEL
# Inputs:
#   - COUNT_FILE: CSV with gene_id column + sample columns
# Outputs:
#   - ORA universe after count-filter: ORA_universe.txt (1 column)
#   - Filtered DEGs per comparison: DESeq2_<comparison>_lfc2_fdr0.01.csv

library(DESeq2)
library(dplyr)

# -----------------------------
# a) Inputs
# -----------------------------
COUNT_FILE <- "LWM_geneCountMatrix.csv"
DESIGN_COL <- "tissuetype"

# thresholds
MIN_ROW_SUM <- 10
LFC_CUTOFF  <- 2
FDR_CUTOFF  <- 0.01

# -----------------------------
# b) Read counts + build sample info from count-matrix column names
# -----------------------------
counts_data <- as.matrix(read.csv(COUNT_FILE, row.names = "gene_id", check.names = FALSE))

# clean sample names: remove "-LFM..." suffix (if present)
colnames(counts_data) <- sub("-LFM.*", "", colnames(counts_data))

# expected format after cleaning: LWM-fb-1, LWM-gut-2, LWM-sg-3
sample_id   <- colnames(counts_data)
tissue_code <- sub("^.*?-([^-]+)-[0-9]+$", "\\1", sample_id)

tissuetype <- dplyr::recode(
  tissue_code,
  fb  = "FatBody",
  gut = "Gut",
  sg  = "Labial",
  .default = NA_character_
)

# sanity check: all samples mapped
stopifnot(!any(is.na(tissuetype)))

colData <- data.frame(
  tissuetype = factor(tissuetype),
  row.names  = sample_id
)
stopifnot(all(colnames(counts_data) == rownames(colData)))

# -----------------------------
# c) Build dds + filter
# -----------------------------
counts_data_int <- round(counts_data)
storage.mode(counts_data_int) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_data_int,
  colData   = colData,
  design    = ~ tissuetype
)

# remove low-count genes
keep <- rowSums(counts(dds)) >= MIN_ROW_SUM
dds <- dds[keep, ]

# ORA universe (genes that entered the model after filtering)
write.table(
  data.frame(Gene = rownames(dds)),
  file = "ORA_universe.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# -----------------------------
# d) Model fitting
# -----------------------------
dds <- DESeq(dds)

# (optional) overall results object
res_all <- results(dds)
summary(res_all)

# -----------------------------
# e) Pairwise comparisons (filtered outputs only)
# -----------------------------
comparisons <- list(
  Labial_vs_FatBody = c(DESIGN_COL, "Labial", "FatBody"),
  Gut_vs_FatBody    = c(DESIGN_COL, "Gut", "FatBody"),
  Gut_vs_Labial     = c(DESIGN_COL, "Gut", "Labial")
)

for (nm in names(comparisons)) {
  
  res <- results(
    dds,
    contrast     = comparisons[[nm]],
    lfcThreshold = LFC_CUTOFF,
    alpha        = FDR_CUTOFF
  )
  
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  res_filt <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj < FDR_CUTOFF) %>%
    dplyr::filter(abs(log2FoldChange) > LFC_CUTOFF)
  
  message("\n=== ", nm, " ===")
  print(summary(res))
  
  write.csv(
    res_filt,
    file = paste0("DESeq2_", nm, "_lfc2_fdr0.01.csv"),
    row.names = FALSE
  )
}