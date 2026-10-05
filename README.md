<div align="center">

# 🧬 Genome-wide DNA Methylation Analysis in CVID

**Illumina EPICv2 array · Differentially methylated regions · Immunodeficiency gene panels · Sibling-pair profiling**

![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-sesame%20%7C%20DMRcate-1a81c2)
![Array](https://img.shields.io/badge/Array-EPICv2%20(hg38)-2b8cbe)
![Status](https://img.shields.io/badge/status-exploratory-orange)

</div>

---

## 📖 Overview

This repository contains an R Markdown pipeline for analysing DNA methylation in patients with **common variable immunodeficiency (CVID)** compared with controls, using the **Illumina Infinium MethylationEPIC v2.0** array.

The workflow takes raw IDAT files through quality control and normalisation, calls **differentially methylated regions (DMRs)** genome-wide, and then zooms in on curated **CVID, PID and PIDCAN gene panels**. It also includes a **sibling-pair analysis** that benchmarks two related individuals against the control distribution.

> ⚠️ **Status:** exploratory analysis. Findings are hypothesis-generating and should be validated before biological interpretation.

---

## 📑 Table of contents

- [Workflow at a glance](#-workflow-at-a-glance)
- [Quality-control thresholds](#-quality-control-thresholds)
- [Getting started](#-getting-started)
- [Input files](#-input-files)
- [Outputs](#-outputs)
- [Limitations](#-limitations)
- [Data availability](#-data-availability-and-privacy)
- [Author](#-author)

---

## 🔬 Workflow at a glance

```mermaid
flowchart TD
    A[Raw IDAT files + SampleSheet] --> B[sesame preprocessing<br/>openSesame, prep = CDB]
    B --> C[Quality control<br/>SNP heatmap, pOOBAH detection p-values]
    C --> D[Probe and sample filtering<br/>SNP / CH probes, missing values]
    D --> E[Beta to M-values]
    E --> F[DMRcate<br/>cpg.annotate + dmrcate, FDR 0.05]
    F --> G[Gene panel filtering<br/>CVID / PID / PIDCAN]
    G --> H[Excel tables]
    G --> I[Locus plots<br/>Gviz + Ensembl, PDF]
    D --> J[Sibling vs control profiling]
    J --> K[Excel tables + PDF plots]
    D --> L[Density and MDS plots]
```

| # | Step | Tools |
|---|------|-------|
| 1 | **Import** sample sheet, build IDAT prefixes | base R |
| 2 | **Preprocess** raw signal to beta values | `sesame` (`CDB`) |
| 3 | **QC**: sample identity (SNP heatmap) and detection p-values (pOOBAH) | `sesame` |
| 4 | **Filter** failed probes/samples, SNP and non-CpG probes | `DMRcate::rmSNPandCH` |
| 5 | **Transform** beta values to M-values | `DMRcate` |
| 6 | **Differential methylation** Patient vs Control (`~ Sample_Group`) | `DMRcate`, `limma` |
| 7 | **Panel filtering** of DMRs by gene lists | `tidyverse` |
| 8 | **Visualise** density, MDS, and DMR locus tracks | `ggplot2`, `Gviz`, `AnnotationHub` |
| 9 | **Sibling analysis** against control range per CpG | `tidyverse`, `ggplot2` |

---

## ✅ Quality-control thresholds

| Level | Criterion |
|-------|-----------|
| Detection | pOOBAH p-value > 0.05 → value masked as `NA` |
| Probe | Removed if detection fails in **> 20 %** of samples |
| Sample | Removed if **≥ 10 %** of probes fail detection |
| Probe | SNP-overlapping and non-CpG (CH) probes removed |
| Probe | Probes with any remaining `NA` removed |
| Replicates | EPICv2 replicate probes collapsed with `epicv2Filter = "mean"` |
| DMP calling | FDR threshold **0.05** |
| Sibling flag | \|Δβ\| ≥ **0.2** between siblings |

---

## 🚀 Getting started

### 1. Requirements

- R (≥ 4.3) and RStudio
- Internet access on first run (reference data are downloaded by `sesameDataCache()` and AnnotationHub)

### 2. Install packages

```r
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "sesame", "sesameData", "DMRcate", "limma", "GenomicRanges",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylationEPICv2manifest",
  "Gviz", "AnnotationHub", "ensembldb"
))

# CRAN
install.packages(c(
  "tidyverse", "data.table", "writexl", "ggplot2",
  "ggrepel", "patchwork"
))
```

### 3. Run

The Rmd uses paths relative to the **parent** of the data folder (e.g. `metilasyon/SampleSheet.csv`). Set the working directory to the folder that *contains* `metilasyon/`, then:

```r
rmarkdown::render("metilasyon/mett.Rmd")
```

or run the document chunk by chunk in RStudio.

> 💡 `View()` calls only work interactively. Comment them out when knitting.

---

## 📂 Input files

```
projeler/
└── metilasyon/
    ├── mett.Rmd            ← analysis (this repository)
    ├── SampleSheet.csv
    ├── CVID.csv
    ├── PID.csv
    ├── PIDCAN.csv
    └── *_Grn.idat / *_Red.idat
```

### `SampleSheet.csv`

| Column | Description |
|--------|-------------|
| `Sentrix_ID` | Array barcode |
| `Sentrix_Position` | Position on the array, e.g. `R01C01` |
| `Sample_Name` | Unique sample label (used in plots) |
| `Sample_Group` | `P` (patient) or `C` (control) |

### Gene lists

`CVID.csv`, `PID.csv` and `PIDCAN.csv` each need a **`Gene`** column with HGNC gene symbols.

---

## 📊 Outputs

| Output | Content |
|--------|---------|
| `cvid_DMR_result.xlsx` · `pid_DMR_results.xlsx` · `pidcan_DMR_result.xlsx` | DMRs overlapping each gene panel |
| `CVID_All_Genes_Plots.pdf` · `PID_All_Genes_Plots.pdf` · `PIDCAN_All_Genes_Plots.pdf` | Genomic tracks per gene and DMR: Ensembl gene models, CpGs, DMRs, group methylation heatmaps and means |
| `sibling_CVID.xlsx` · `sibling_PID.xlsx` · `sibling_PIDCAN.xlsx` | Sibling vs control comparison with `gene_summary` and `cpg_detail` sheets |
| `sibling_*_plots.pdf` | Per-gene plots: control min–max band, control mean, sibling values |

---

## 🧭 Limitations

- **Confounders:** the model contains only the group variable. Cell-type composition, age, sex and batch are not adjusted for and can be strong confounders in blood-derived samples.
- **Relatedness:** the sibling pair is part of the main comparison without modelling the relationship.
- **Sibling analysis:** flags are exploratory, and the control min–max range is sensitive to the number of controls.
- **Annotation sources:** DMR gene overlap uses the DMRcate (Ensembl/GENCODE-based) annotation, whereas the sibling analysis uses the Illumina manifest (`UCSC_RefGene_Name`). Results may differ.
- **Panel selection:** gene panels are applied after genome-wide DMR calling, so panel-level findings should be interpreted as post hoc.
- **Reproducibility:** the latest Ensembl release is pulled from AnnotationHub at run time. Record `sessionInfo()` and the Ensembl version with your results.

---

## 🔒 Data availability and privacy

Raw IDAT files, sample sheets and per-patient result tables are **not** included in this repository because they contain patient-level data. Data access is subject to ethical approval and patient consent. For data-related requests, please contact the author.

---

## 👩‍🔬 Author

---

## 📄 License

Add a license of your choice (e.g. MIT for the code) before making the repository public.
