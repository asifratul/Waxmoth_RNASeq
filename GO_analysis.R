################### GO analysis ######################

#### Load libraries####

library(clusterProfiler)
library(GO.db)
library(GSEABase)
library(tidyverse)
library(readxl)
library(data.table)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

################# extract D. melanogaster GO terms (TERM2NAME) #####################
# get entrez to flybase
keytypes(org.Dm.eg.db)
head(keys(org.Dm.eg.db, keytype = "FLYBASE"))
k <- keys(org.Dm.eg.db, keytype = "FLYBASE")
cols <- c("ENTREZID", "FLYBASECG")
flybase_map <- AnnotationDbi::select(org.Dm.eg.db, keys = k, columns = cols, keytype = "FLYBASE")
head(flybase_map)

entrez_ids <- keys(org.Dm.eg.db, keytype = "ENTREZID")

## ---------------------------------------------------------
## 1. Function to extract TERM2NAME for given ont & level
## ---------------------------------------------------------
extract_TERM2NAME <- function(ont, level, outdir = "Eggnoc/TERM2NAME") {
  message(sprintf("Processing ontology = %s | level = %d", ont, level))
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  ggo <- groupGO(
    gene      = entrez_ids,
    OrgDb     = org.Dm.eg.db,
    ont       = ont,
    level     = level,
    readable  = FALSE
  )
  
  term2name <- ggo %>%
    filter(Count > 0) %>%
    select(ID, Description)
  
  outfile <- file.path(outdir, paste0("TERM2NAME_", ont, level, ".csv"))
  write.csv(term2name, outfile, row.names = FALSE)
}

## ---------------------------------------------------------
## 2. Run required combinations for this study
##    (MF4, BP4, CC3)
## ---------------------------------------------------------
extract_TERM2NAME("MF", 4)
extract_TERM2NAME("BP", 4)
extract_TERM2NAME("CC", 3)


###################### Extract GO from eggNOG-mapper ####################

emapper <- read_excel("GO_annotation_files/out.emapper.annotations.xlsx", sheet = 1)

Gene2emapperGO <- emapper %>%
  dplyr::select(query, GOs) %>%
  filter(!is.na(GOs)) %>%
  rowwise() %>%
  mutate(GO_list = strsplit(GOs, ",")) %>%
  unnest(cols = c(GO_list)) %>%
  filter(!grepl("-", GO_list)) %>%
  dplyr::select(Protein_ID = query, GO_ID = GO_list) %>%
  distinct()

# Final GO universe comes ONLY from eggNOG
loc_go_combined <- Gene2emapperGO
write.csv(loc_go_combined, "Gene2AllGO_Unique.csv", row.names = FALSE)

#################################################################################
############ Parent to offspring (MF4 / BP4 / CC3) ###################
#################################################################################

GO_ID <- loc_go_combined %>% dplyr::select(GO_ID) %>% distinct()
go_ids <- GO_ID$GO_ID %>% as.character()
myCollection <- GOCollection(go_ids)

TERM2GENE_DIR <- "Eggnoc/TERM2GENE"
TERM2NAME_DIR <- "Eggnoc/TERM2NAME"
dir.create(TERM2GENE_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: build TERM2GENE for a GO level file (MF4/BP4/CC3)
# -----------------------------
make_TERM2GENE <- function(term2name_csv, offspring_obj, parent_col, out_csv,
                           loc_go = loc_go_combined, collection = myCollection) {
  
  slimdf <- read.csv(term2name_csv) %>% dplyr::select(ID)
  
  gomap <- as.list(offspring_obj[slimdf$ID])
  
  # Map each parent ID to the observed child GO IDs
  slimdf$ids <- sapply(lapply(gomap, intersect, ids(collection)), paste, collapse = ";")
  slimdf$ids[slimdf$ids == "character(0)"] <- ""
  
  # Parent-to-offspring table
  P2O <- (slimdf %>% data.table())[, .(GO_ID = unlist(strsplit(ids, ";"))), by = ID]
  setnames(P2O, c(parent_col, "GO_ID"))
  
  # Add parent terms explicitly (so parents always appear as self-mapped)
  UniqueDat <- unique(P2O[, 1])
  P2O1 <- rbind(P2O, data.frame(tmp = UniqueDat, GO_ID = UniqueDat), use.names = FALSE)
  setnames(P2O1, c(parent_col, "GO_ID"))
  
  # Join to genes
  full_join(loc_go, P2O1, by = "GO_ID", relationship = "many-to-many") %>%
    na.omit() %>%
    dplyr::select(all_of(parent_col), "Protein_ID") %>%
    write.csv(out_csv, row.names = FALSE)
  
  message("✅ Wrote: ", out_csv)
}

# Run for MF4 / BP4 / CC3
make_TERM2GENE(
  term2name_csv = file.path(TERM2NAME_DIR, "TERM2NAME_MF4.csv"),
  offspring_obj = GOMFOFFSPRING,
  parent_col    = "MF.GO4",
  out_csv       = file.path(TERM2GENE_DIR, "Term2Gene_MF4.csv")
)

make_TERM2GENE(
  term2name_csv = file.path(TERM2NAME_DIR, "TERM2NAME_BP4.csv"),
  offspring_obj = GOBPOFFSPRING,
  parent_col    = "BP.GO4",
  out_csv       = file.path(TERM2GENE_DIR, "Term2Gene_BP4.csv")
)

make_TERM2GENE(
  term2name_csv = file.path(TERM2NAME_DIR, "TERM2NAME_CC3.csv"),
  offspring_obj = GOCCOFFSPRING,
  parent_col    = "CC.GO3",
  out_csv       = file.path(TERM2GENE_DIR, "Term2Gene_CC3.csv")
)

########################################################
##################### KEGG #############################

kegg_clean <- emapper %>%
  dplyr::select(query, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway)) %>%
  rowwise() %>%
  mutate(KEGG_list = strsplit(KEGG_Pathway, ",")) %>%
  unnest(cols = c(KEGG_list)) %>%
  mutate(KEGG_list = str_trim(KEGG_list)) %>%
  filter(!grepl("^map", KEGG_list, ignore.case = TRUE), !grepl("-", KEGG_list)) %>%
  dplyr::select(KEGG_ID = KEGG_list, Protein_ID = query)

write.csv(kegg_clean, file.path(TERM2GENE_DIR, "Term2Gene_KEGG.csv"), row.names = FALSE)

R.utils::setOption("clusterProfiler.download.method", "auto")
ko.pathways <- download_KEGG(species = "ko")
term2name <- ko.pathways$KEGGPATHID2NAME
write.csv(term2name, file.path(TERM2NAME_DIR, "TERM2NAME_KEGG.csv"), row.names = FALSE)

##################################### Enzyme ######################################################

kegg1 <- emapper %>%
  select(query, EC) %>%
  filter(!is.na(EC), EC != "-") %>%
  mutate(Enzymes = strsplit(EC, ",")) %>%
  unnest(cols = c(Enzymes)) %>%
  mutate(
    Enzymes = str_trim(Enzymes),
    Enzymes = na_if(Enzymes, "-"),
    Enzymes = str_replace(Enzymes, "^([0-9]+\\.[0-9]+)\\..*$", "\\1"),
    Enzymes = if_else(str_detect(Enzymes, "^\\d+\\.\\d$"), str_replace(Enzymes, "(\\d+)\\.(\\d)$", "\\1.0\\2"), Enzymes)
  ) %>%
  filter(!is.na(Enzymes), Enzymes != "") %>% distinct(query, Enzymes)

write.csv(kegg1, file.path(TERM2GENE_DIR, "Term2Gene_Enzymes.csv"), row.names = FALSE)

############################################################################################
######################## ORA analysis ######################################################

run_ORA <- function(diff_gene_file, universe_file, term_pairs, output_dir = "ORA_results", suffix = "") {
  if (suffix != "") output_dir <- file.path(output_dir, paste0("ORA_", suffix))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  UP_T <- read.delim(diff_gene_file, header = TRUE)
  Universe <- read.delim(universe_file, header = TRUE)
  
  for (term_name in names(term_pairs)) {
    TERM2GENE <- read.csv(term_pairs[[term_name]]$TERM2GENE)
    TERM2NAME <- read.csv(term_pairs[[term_name]]$TERM2NAME)
    
    enricher(
      gene = UP_T[[1]], universe = Universe[[1]], TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 5000
    ) %>% as.data.frame() %>%
      write.csv(file.path(output_dir, paste0("ORA_", term_name, ifelse(suffix=="","",paste0("_",suffix)), ".csv")), row.names = FALSE)
  }
}

term_files_list <- list(
  MF4    = list(TERM2GENE = "Eggnoc/TERM2GENE/Term2Gene_MF4.csv", TERM2NAME = "Eggnoc/TERM2NAME/TERM2NAME_MF4.csv"),
  BP4    = list(TERM2GENE = "Eggnoc/TERM2GENE/Term2Gene_BP4.csv", TERM2NAME = "Eggnoc/TERM2NAME/TERM2NAME_BP4.csv"),
  CC3    = list(TERM2GENE = "Eggnoc/TERM2GENE/Term2Gene_CC3.csv", TERM2NAME = "Eggnoc/TERM2NAME/TERM2NAME_CC3.csv"),
  KEGG   = list(TERM2GENE = "Eggnoc/TERM2GENE/Term2Gene_KEGG.csv", TERM2NAME = "Eggnoc/TERM2NAME/TERM2NAME_KEGG.csv"),
  Enzyme = list(TERM2GENE = "Eggnoc/TERM2GENE/Term2Gene_Enzymes.csv", TERM2NAME = "Eggnoc/TERM2NAME/Term2Name_Enzyme.csv")
)

# ------------------------------------------------------------------
# ORA execution
# ------------------------------------------------------------------

# Labial_vs_Fatbody
run_ORA(
  diff_gene_file = "DESeq2_subset/DESeq2_Labial_vs_FatBody_lfc2_fdr0.01.csv",
  universe_file  = "DESeq2_subset/Universe_Labial_FatBody.txt",
  term_pairs     = term_files_list,
  output_dir     = "Eggnoc/ORA_results",
  suffix         = "Labial_vs_Fatbody")

# Gut_vs_FatBody
run_ORA(
  diff_gene_file = "DESeq2_subset/DESeq2_Gut_vs_FatBody_lfc2_fdr0.01.csv",
  universe_file  = "DESeq2_subset/Universe_Gut_FatBody.txt",
  term_pairs     = term_files_list,
  output_dir     = "Eggnoc/ORA_results",
  suffix         = "Gut_vs_FatBody"
)


# Gut_vs_Labial
run_ORA(
  diff_gene_file = "DESeq2_subset/DESeq2_Gut_vs_Labial_lfc2_fdr0.01.csv",
  universe_file  = "DESeq2_subset/Universe_Gut_Labial.txt",
  term_pairs     = term_files_list,
  output_dir     = "Eggnoc/ORA_results",
  suffix         = "Gut_vs_Labial")


# If you later generate SignalP / SignalP_UP subsets,
# keep the SAME comparison names and only change the DEG/Universe filenames.
