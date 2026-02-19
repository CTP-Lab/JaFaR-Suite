# JaFaR-Suite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18681624.svg)](https://doi.org/10.5281/zenodo.18681624)

R/Shiny suite to convert MZmine outputs ↔ W4M tables, run univariate/multivariate stats, and export to GNPS/MetGem for molecular networking.

**JaFaR-Suite** is a collection of **R/Shiny apps** designed to bridge **MZmine ↔ Workflow4Metabolomics (W4M/Galaxy)**, perform **univariate & multivariate statistics**, and generate **GNPS / SIRIUS / MetGem-ready files** for molecular networking and downstream annotation.

> Not affiliated with MZmine, Workflow4Metabolomics (W4M), GNPS, or MetGem.

---

## What this suite does (high-level)

- Converts **MZmine exports → W4M triplet** (`dataMatrix`, `sampleMetadata`, `variableMetadata`)
- Processes and filters features **inside W4M** (format check, fold-change vs blanks, quality metrics, batch correction, RSD filtering)
- Converts the updated **W4M triplet → a MZmine-like CSV matrix** and **filters the original MGF** accordingly
- Runs **univariate and multivariate statistics** from the **post-W4M triplet**
- Produces a final **filtered CSV + filtered MGF** right before annotation / networking (GNPS / MetGem / SIRIUS / other annotation tools)

---

## Typical workflow (recommended)

```text
(1) MZmine exports
    - quant table (CSV)
    - MGF (MS/MS)

        │
        ▼
(2) App: JaFaR_MZmine_to_W4M_bidirectional  [direction: MZmine → W4M]
    Output: W4M triplet
      - dataMatrix.tsv
      - sampleMetadata.tsv
      - variableMetadata.tsv

        │
        ▼
(3) W4M/Galaxy step (outside this repo)
    The W4M triplet is processed in W4M/Galaxy in the following order:

    1) Check format
    2) Intensity check to compute mean fold-change (samples vs blanks)
    3) Filter features based on a chosen mean fold-change threshold
    4) Quality metrics (pre-correction)
    5) Batch correction (e.g., QC-based drift correction)
    6) Quality metrics (post-correction)
    7) Filter features (e.g. RSD(QC) > 30%)

    Output: updated W4M triplet (post-W4M)
      - dataMatrix.tsv
      - sampleMetadata.tsv
      - variableMetadata.tsv

        │
        ├───────────────────────────────────────────────┐
        │                                               │
        ▼                                               ▼
(4) App: JaFaR_MZmine_to_W4M_bidirectional       (5) Stats apps (from post-W4M triplet)
    [direction: W4M → MZmine]                         - Univariate
    Inputs: post-W4M triplet + original MGF            - Multivariate
    Outputs:                                           Outputs include an imputed matrix CSV
    - filtered MGF (based on post-W4M features)         (used later before annotation)
    - "non-annotated" CSV matrix (MZmine-like)

        │                                               │
        └───────────────────────────────┬───────────────┘
                                        ▼
(6) App: JaFaR_CSV_and_MGF_filter_just_before_annotation
    Inputs:
    - imputed matrix CSV (from multivariate stats)
    - non-annotated CSV matrix (from W4M → MZmine conversion)
    - filtered MGF (from W4M → MZmine conversion)
    Outputs:
    - final CSV + final MGF (synchronized) for GNPS / MetGem / SIRIUS / other annotation tools
```

---

## Apps included

1) **JaFaR_MZmine_to_W4M_bidirectional**  
Folder: `apps/JaFaR_MZmine_to_W4M_bidirectional/`  
Used twice in the workflow:
- MZmine → W4M: build W4M triplet tables from MZmine exports  
- W4M → MZmine: convert post-W4M triplet to a MZmine-like CSV and filter the original MGF  

2) **JaFaR_Untargeted_Metabolomics_Univariate_Pipeline**  
Folder: `apps/JaFaR_Untargeted_Metabolomics_Univariate_Pipeline/`  
- Univariate statistics from the post-W4M triplet  
- Exports tables + plots to an output folder  

3) **JaFaR_Untargeted_Metabolomics_Multivariat_Pipelinee**  
Folder: `apps/JaFaR_Untargeted_Metabolomics_Multivariate_Pipeline/`  
- Multivariate statistics/diagnostics from the post-W4M triplet  
- Produces an imputed matrix CSV used in the final pre-annotation filtering step  

4) **JaFaR_CSV_and_MGF_filter_just_before_annotation**  
Folder: `apps/JaFaR_CSV_and_MGF_filter_just_before_annotation/`  
Final synchronization step just before annotation:
- input (A) imputed matrix CSV (multivariate)  
- input (B) non-annotated CSV matrix (W4M → MZmine conversion)  
- input (C) filtered MGF (W4M → MZmine conversion)  
- output: final filtered CSV + final filtered MGF (feature IDs aligned)  

---

## Requirements

- R (recommended >= 4.2)  
- Packages are installed automatically using `install_deps.R`  

Tip (Windows): install **Rtools** if you need compilation for some packages.

---

## Installation

1) Download or clone the repository:

https://github.com/CTP-Lab/JaFaR-Suite

2) Install dependencies (from R):

- Open R (or RStudio)
- Set your working directory to the repository folder
- Run:

```r
source("install_deps.R")
install_deps()
```

---

## Run the apps

- Open R (or RStudio)
- Set your working directory to the repository folder
- Run:

```r
source("scripts/run_app.R")
```

Then run any app with:

```r
run_app("apps/JaFaR_MZmine_to_W4M_bidirectional")
```

Other apps:

```r
run_app("apps/JaFaR_Untargeted_Metabolomics_Univariate_Pipeline")
run_app("apps/JaFaR_Untargeted_Metabolomics_Multivariate_Pipeline")
run_app("apps/JaFaR_CSV_and_MGF_filter_just_before_annotation")
```

---

## Inputs / Outputs (high level)

- MZmine → W4M  
Input: MZmine quant table + sample metadata  
Output: W4M triplet tables (`dataMatrix.tsv`, `sampleMetadata.tsv`, `variableMetadata.tsv`)

- W4M processing  
Input: W4M triplet  
Output: updated W4M triplet (post-W4M)

- W4M → MZmine (+MGF filtering)  
Input: post-W4M triplet + original MGF  
Output: filtered MGF + non-annotated CSV matrix

- Univariate stats  
Input: post-W4M triplet  
Output: univariate outputs (tables/plots)

- Multivariate stats  
Input: post-W4M triplet  
Output: multivariate outputs + imputed matrix CSV

- Pre-annotation filter  
Input: imputed matrix CSV + non-annotated CSV + filtered MGF  
Output: final synchronized CSV + final MGF

---

## How to cite

Fauquet, J. (2026). JaFaR-Suite (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18681624

See also:
- `CITATION.cff`
- `CITATION.bib`

---

## Releasing + DOI (Zenodo)

See `RELEASING.md`.

Summary:
- enable the repo in Zenodo (GitHub integration)
- create a GitHub Release (`v1.0.0`, `v1.0.1`, …)
- Zenodo archives the release and mints a DOI

---

## License

MIT License (see `LICENSE`).

---

## Support / Issues

Please open a GitHub Issue with:
- which app(s) you used and in which direction (MZmine→W4M or W4M→MZmine)
- your input file types (MZmine export format, W4M tables, MGF)
- `sessionInfo()`
- and a minimal reproducible example if possible
