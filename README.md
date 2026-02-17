# W4MineR-Suite
R/Shiny suite to convert MZmine outputs ↔ W4M tables, run univariate/multivariate stats, and export to GNPS/MetGem for molecular networking.

# W4MineR Suite

**W4MineR Suite** is a collection of **R/Shiny apps** designed to bridge **MZmine ↔ Workflow4Metabolomics (W4M/Galaxy)**, perform **univariate & multivariate statistics**, and generate **GNPS / SIRIUS / MetGem-ready files** for molecular networking and downstream annotation.

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
(2) App: W4MineR_MZmine_to_W4M_bidirectional  [direction: MZmine → W4M]
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
    3) Filter features based on the chosen mean fold-change threshold
    4) Quality metrics (pre-correction)
    5) Batch correction (e.g., QC-based drift correction)
    6) Quality metrics (post-correction)
    7) Filter features with RSD(QC) > 30%

    Output: updated W4M triplet (post-W4M)
      - dataMatrix.tsv
      - sampleMetadata.tsv
      - variableMetadata.tsv

        │
        ├───────────────────────────────────────────────┐
        │                                               │
        ▼                                               ▼
(4) App: W4MineR_MZmine_to_W4M_bidirectional       (5) Stats apps (from post-W4M triplet)
    [direction: W4M → MZmine]                         - Univariate
    Inputs: post-W4M triplet + original MGF            - Multivariate
    Outputs:                                           Outputs include an imputed matrix CSV
    - filtered MGF (based on post-W4M features)         (used later before annotation)
    - "non-annotated" CSV matrix (MZmine-like)

        │                                               │
        └───────────────────────────────┬───────────────┘
                                        ▼
(6) App: W4MineR_CSV_and_MGF_filter_just_before_annotation
    Inputs:
    - imputed matrix CSV (from multivariate stats)
    - non-annotated CSV matrix (from W4M → MZmine conversion)
    - filtered MGF (from W4M → MZmine conversion)
    Outputs:
    - final CSV + final MGF (synchronized) for GNPS / MetGem / SIRIUS / other annotation tools
