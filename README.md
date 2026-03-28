# Bulk RNAseq Differential Expression Analysis Pipeline

Pipeline for conducting differential expression analyses on bulk RNA-seq data from neurological and infectious disease studies. Developed by **Jon Sanchez-Valle** (BSC) and **Andrea Marti Sarrias** (UB).

---

## Overview

This pipeline performs end-to-end differential expression analysis using **DESeq2** and **limma/voom**, starting from raw count data and metadata downloaded from [GREIN](https://www.ilincs.org/apps/grein/) and [ARCHS4](https://maayanlab.cloud/archs4/). It covers quality control, data preprocessing, differential expression analysis, and downstream GSEA preparation across multiple disease datasets.

**Diseases covered:**
- COVID-19 (including Long COVID)
- Alzheimer's Disease (AD)
- Parkinson's Disease (PD)

---

## Datasets

| Dataset ID | Disease | Sample Type |
|---|---|---|
| GSE95587 | Alzheimer's Disease | Brain tissue |
| GSE147507 | COVID-19 | Cell lines |
| GSE174745 | COVID-19 | Dopaminergic neurons (hiPSC-derived) |
| GSE179923 | COVID-19 | hiPSC-derived cells |
| GSE152418 | COVID-19 | PBMCs |
| GSE251849 | Long COVID | PBMCs |
| GSE136666 | Parkinson's Disease | Brain tissue |
| GSE216281 | Parkinson's Disease | Brain tissue |
| GSE68719 | Parkinson's Disease | Brain tissue |

Counts and QC files are downloaded from **GREIN**; additional counts for comparative analyses are obtained from **ARCHS4**.

---

## Requirements

**R version:** ≥ 4.0 recommended

**R packages:**

```r
install.packages(c("data.table", "lsa", "dplyr", "devtools", "calibrate", "dendextend"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "PCAtools", "sva", "edgeR",
  "EnhancedVolcano", "VennDiagram", "SummarizedExperiment"
))
```

---

## Directory Structure

Before running the pipeline, the following structure is expected under a `GREIN/` directory:

```
GREIN/
├── QC/              # QC tables downloaded from GREIN (one .txt per study)
├── Counts/          # Raw count matrices downloaded from GREIN
└── Metadata/        # Sample metadata tables downloaded from GREIN
```

Output directories are created automatically by the pipeline.

---

## Usage

The script is run from the command line with a single argument specifying the analysis step:

```bash
Rscript RNAseq_analyses.R <step>
```

Steps must be run **in order**.

---

## Pipeline Steps

### 1. `check_quality_of_samples`

Performs a preliminary quality assessment of all studies downloaded from GREIN.

- Parses QC files and extracts the percentage of aligned reads per sample.
- Retains samples with **≥ 70% aligned reads**.
- Filters studies with **at least 10 passing samples**.
- Outputs summary tables:
  - `GREIN/Summary_by_study.txt` — QC metrics for all studies.
  - `GREIN/Summary_by_study_working_datasets.txt` — Studies that pass the quality threshold.

---

### 2. `create_summarized_experiment`

Builds `SummarizedExperiment` objects for each dataset, harmonising metadata across studies.

- Selects samples that pass the 70% alignment threshold.
- Standardises metadata column names (`sample`, `category`, `sex`, `age`, `tissue`, etc.).
- Labels samples consistently as `case` or `control`.
- Saves two versions of each dataset:
  - `GREIN/SummarizedExperiments/` — all samples.
  - `GREIN/SummarizedExperiments_CaseControl/` — case and control samples only.

---

### 3. `organize_data_as_needed`

Filters lowly expressed genes and prepares objects for differential expression.

- Creates `DGEList` objects (edgeR) and computes log-CPM values.
- Removes genes with log-CPM ≤ 1 in a number of samples equal to the smallest group size.
- Calculates TMM normalisation factors (`calcNormFactors`).
- Builds `DESeqDataSet` objects with dataset-specific design formulas accounting for available covariates (sex, tissue, treatment, etc.).
- Outputs saved to `GREIN/processed_summarized_experiments/`:
  - `dgeobj_*.rds` — filtered and normalised DGEList objects (for limma/voom).
  - `deseqobj_*.rds` — DESeqDataSet objects (for DESeq2).

---

### 4. `differential_expression_analyses`

Runs differential expression analysis with both **DESeq2** and **limma/voom** for each dataset.

- **DESeq2:** `case` vs. `control` contrast with Benjamini-Hochberg correction (α = 0.05).
- **limma/voom:** Mean-variance modelling with `voom`, linear model fitting (`lmFit`), empirical Bayes smoothing (`eBayes`), and `topTable` for results.
- Design matrices are dataset-specific and may include covariates such as sex, tissue, treatment, severity, or age.
- Results (log fold-changes, adjusted p-values) are written to `GREIN/Differential_expression_profiles/`.
- A summary table with the number of up- and down-regulated DEGs per dataset and method is saved to `GREIN/Summary_number_sDEGs.txt`.

---

### 5. `rankedlistofgenes`

Prepares ranked gene lists for Gene Set Enrichment Analysis (GSEA).

- Ranks genes by `log2FoldChange` (descending) from DESeq2 results.
- Writes `.rnk` files to `GREIN/Ranks/`, compatible with GSEA Preranked.
- Includes a template GSEA CLI command (commented out) using:
  - KEGG, Reactome, and GO Biological Process gene sets (MSigDB v2023.1).
  - Settings: 1000 permutations, gene set size 15–500, top 20 plots.
  - Output written to `GREIN/Pathways/`.

---

## Output Summary

| Directory | Contents |
|---|---|
| `GREIN/Aligned_reads/` | Per-study tables of alignment rates |
| `GREIN/Summary_by_study.txt` | QC summary for all studies |
| `GREIN/Summary_by_study_working_datasets.txt` | Studies passing QC filters |
| `GREIN/SummarizedExperiments/` | Full SummarizedExperiment objects |
| `GREIN/SummarizedExperiments_CaseControl/` | Case/control-only SummarizedExperiment objects |
| `GREIN/processed_summarized_experiments/` | Normalised DGEList and DESeqDataSet objects |
| `GREIN/Differential_expression_profiles/` | DEG results tables (DESeq2 and limma) |
| `GREIN/Summary_number_sDEGs.txt` | Summary of DEG counts per dataset and method |
| `GREIN/Ranks/` | Ranked gene lists for GSEA |
| `GREIN/Pathways/` | GSEA output (generated externally via GSEA CLI) |

---

## Authors

- **Jon Sanchez-Valle** — Barcelona Supercomputing Center (BSC)
- **Andrea Marti Sarrias** — Universitat de Barcelona (UB)
