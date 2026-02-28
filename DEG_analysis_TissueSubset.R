# DESeq2 differential expression with subsetting (no external metadata)
# - Derives sample metadata (tissuetype) from count-matrix column names
# - Subsets tissues and fits an independent model per subset
# - Rounds fractional counts (IsoformAnalyzer-style) to integers for DESeq2
# - Writes: (1) filtered DEG CSV per contrast, (2) ORA universe per subset
#
# Expected sample name format (after cleaning):
#   LWM-fb-1, LWM-gut-2, LWM-sg-3
# Column-name cleaning removes everything from "-LFM" to end.

library(DESeq2)
library(dplyr)

# -----------------------------
# Function
# -----------------------------
run_deseq_subset_no_meta <- function(
    count_file,
    keep_tissues,
    outdir       = "DESeq2_subset",
    design_col   = "tissuetype",
    min_row_sum  = 10,
    lfc_cutoff   = 2,
    fdr_cutoff   = 0.01
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- read counts
  counts_data <- as.matrix(read.csv(count_file, row.names = "gene_id", check.names = FALSE))
  
  # ---- clean sample names: remove "-LFM..." suffix
  colnames(counts_data) <- sub("-LFM.*", "", colnames(counts_data))
  sample_id <- colnames(counts_data)
  
  # ---- derive tissuetype from sample names
  # captures fb/gut/sg from: <anything>-<code>-<rep>
  tissue_code <- sub("^.*?-([^-]+)-[0-9]+$", "\\1", sample_id)
  
  tissuetype <- dplyr::recode(
    tissue_code,
    fb  = "FatBody",
    gut = "Gut",
    sg  = "Labial",
    .default = NA_character_
  )
  
  # fail fast if sample naming is unexpected
  if (any(is.na(tissuetype))) {
    bad <- unique(sample_id[is.na(tissuetype)])
    stop(
      "Could not map tissuetype for some samples.\n",
      "Expected names like: LWM-fb-1, LWM-gut-2, LWM-sg-3 (after -LFM cleanup).\n",
      "Unmapped examples: ", paste(head(bad, 5), collapse = ", ")
    )
  }
  
  colData <- data.frame(
    tissuetype = factor(tissuetype),
    row.names  = sample_id
  )
  
  # ---- subset samples first (independent model)
  colData_sub <- colData[colData[[design_col]] %in% keep_tissues, , drop = FALSE]
  colData_sub[[design_col]] <- droplevels(colData_sub[[design_col]])
  
  if (nrow(colData_sub) < 2) {
    stop("After subsetting, fewer than 2 samples remain. Check keep_tissues.")
  }
  
  counts_sub <- counts_data[, rownames(colData_sub), drop = FALSE]
  stopifnot(all(colnames(counts_sub) == rownames(colData_sub)))
  
  # ---- DESeq2 requires integers; IsoformAnalyzer often yields fractional counts
  counts_sub <- round(counts_sub)
  storage.mode(counts_sub) <- "integer"
  
  # ---- build dds
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = colData_sub,
    design    = as.formula(paste0("~ ", design_col))
  )
  
  # ---- low-count filter
  keep <- rowSums(counts(dds)) >= min_row_sum
  dds <- dds[keep, ]
  
  # ---- ORA universe (genes that entered the model after filtering)
  write.table(
    data.frame(Gene = rownames(dds)),
    file = file.path(outdir, paste0("Universe_", paste(keep_tissues, collapse = "_"), ".txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  # ---- fit
  dds <- DESeq(dds)
  
  # ---- pairwise comparisons within subset only
  # preserve user-specified order instead of alphabetical factor order
  tissue_levels <- keep_tissues[keep_tissues %in% levels(colData(dds)[[design_col]])]
  
  if (length(tissue_levels) < 2) {
    stop("Need at least two tissue levels after filtering to run contrasts.")
  }
  
  # generate ordered pairwise comparisons based on keep_tissues order
  comps <- list()
  for (i in seq_len(length(tissue_levels) - 1)) {
    for (j in seq((i + 1), length(tissue_levels))) {
      comps[[length(comps) + 1]] <- c(tissue_levels[i], tissue_levels[j])
    }
  }
  
  for (cmp in comps) {
    nm <- paste0(cmp[1], "_vs_", cmp[2])
    
    res <- results(
      dds,
      contrast = c(design_col, cmp[1], cmp[2]),
      lfcThreshold = lfc_cutoff,
      alpha = fdr_cutoff
    )
    
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    
    # output filtered DEGs only
    res_filt <- res_df %>%
      filter(!is.na(padj)) %>%
      filter(padj < fdr_cutoff) %>%
      filter(abs(log2FoldChange) > lfc_cutoff)
    
    write.csv(
      res_filt,
      file = file.path(
        outdir,
        paste0(
          "DESeq2_", nm,
          "_lfc", lfc_cutoff,
          "_fdr", fdr_cutoff,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  invisible(TRUE)
}

# -----------------------------
# runs
# -----------------------------
COUNT_FILE <- "LWM_geneCountMatrix.csv"

# 1) Labial vs FatBody only
run_deseq_subset_no_meta(
  count_file   = COUNT_FILE,
  keep_tissues = c("Labial", "FatBody"),
  outdir       = "DESeq2_subset")

# 2) Gut vs FatBody only
run_deseq_subset_no_meta(
  count_file   = COUNT_FILE,
  keep_tissues = c("Gut", "FatBody"),
  outdir       = "DESeq2_subset")

# 3) Gut vs Labial only
run_deseq_subset_no_meta(
  count_file   = COUNT_FILE,
  keep_tissues = c("Gut", "Labial"),
  outdir       = "DESeq2_subset")
