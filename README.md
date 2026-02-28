# Waxmoth_RNASeq

This repository contains the analysis scripts used in the manuscript:

**Transcriptomic analyses of labial glands and gut tissues from two wax moths, *Achroia grisella* and *Galleria mellonella***

## Authors
**Reginald Young¹**, **Khandaker Asif Ahmed²**, **Leon Court¹**, **Rahul Rane³**, **Tom Walsh¹**, **Gunjan Pandey¹***  

¹ CSIRO Environment, Acton, ACT 2601, Australia  
² CSIRO Agriculture and Food, Acton, ACT 2601, Australia  
³ CSIRO Health and Biosecurity, Parkville, VIC 3052, Australia  

---
## Scripts

### 1. `DEG_analysis_Overall.R`
Runs DESeq2 using all samples in a single model and performs pairwise contrasts:
- Labial vs FatBody
- Gut vs FatBody
- Gut vs Labial

Outputs:
- `ORA_universe.txt`
- `DESeq2_<comparison>_lfc2_fdr0.01.csv`

---

### 2. `DEG_analysis_TissueSubset.R`
Fits DESeq2 independently for tissue subsets.

For each subset:
- Builds an independent model
- Generates a subset-specific ORA universe
- Outputs filtered DEG tables

Outputs:
- `Universe_<Tissue1>_<Tissue2>.txt`
- `DESeq2_<Tissue1>_vs_<Tissue2>_lfc2_fdr0.01.csv`

---

### 3. `GO_analysis.R`
Performs functional enrichment analysis using:

- GO Molecular Function (Level 4)
- GO Biological Process (Level 4)
- GO Cellular Component (Level 3)
- KEGG pathways
- Enzyme classes (EC)

Annotation source: eggNOG-mapper

The script:
- Generates TERM2NAME and TERM2GENE tables
- Runs ORA using DESeq2 DEG outputs
- Writes enrichment result tables per comparison

---

## Notes

- Sample metadata is inferred from count matrix column names.
- Count values are rounded to integers before DESeq2 fitting.
- ORA is performed using eggNOG-derived annotations.
