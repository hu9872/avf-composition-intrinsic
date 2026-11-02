# Where the data come from, and what was done to them

This project **reanalyses public data**. No new sequencing was performed. If you reuse
anything here, the original depositors must be credited.

## Source datasets

| Accession | File used | Role in this study |
|---|---|---|
| **GSE305423** | `GSE305423_RAW.tar` (10x per-sample triplets) | Primary single-cell atlas: 66,168 cells, 12 post-operative AVF samples (6 failed, 6 matured) plus pre-access veins |
| **GSE220796** | `GSE220796_rVazquez-Padron_counts.txt.gz` | Independent bulk cohort (paired vein / juxta-anastomotic AVF); used for the compartment-signature check |
| **GSE119296** | `GSE119296_Vazquez-Padron_RNAseq_080718_counts.txt.gz` | Positive control: vein → AVF program shift must be detectable, otherwise the pipeline is insensitive |
| GSE314264 | `GSE314264_Final_counts-AVF_timepoints.csv.gz` | Exploratory only; not used in the final analysis |
| GSE244317, GSE39488 | — | Screened for eligibility; excluded (no matching design) |

## How to obtain them

```bash
python code/00_download.py      # FTP from NCBI GEO, ~23 GB including intermediates
python code/01_prepare.py       # unpack 10x triplets, build AnnData + sample ledger
```

Both scripts are idempotent: existing files are skipped, interrupted downloads resume from a
`.part` file.

## Sample ledger

The GSM → group mapping is hard-coded in `code/01_prepare.py` and was verified against GEO
sample titles on 2026-10-31. Two early-AVF samples (`GSM9176937`, `GSM9176938`) come from the
same patient (suffix 1A / 1B) and are treated as one donor throughout.

## What is *not* available

* **Empty-droplet (raw) matrices** — GEO provides filtered matrices only. This is why
  SoupX/decontX ambient correction could not be run; the manuscript states this explicitly and
  substitutes a marker-source screen plus a rule-sensitivity analysis
  (`code/05g_gene_origin_all.py`, `code/31_rule_floor_sensitivity.py`).
* **Patient-level clinical covariates** — not deposited. Age range (28–72 years) is all that is
  available; no adjustment for comorbidity or dialysis vintage is possible.

## Licence and attribution

The source data remain under the terms of their GEO records. The two original publications are:

1. Martinez L, et al. The single-cell landscape of the human vein after arteriovenous fistula
   creation and implications for maturation failure. *Kidney Int* 2026;109:160–177.
   https://doi.org/10.1016/j.kint.2025.08.024
2. Martinez L, et al. The transcriptomics of the human vein transformation after arteriovenous
   fistula anastomosis uncovers layer-specific remodeling and hallmarks of maturation failure.
   *Kidney Int Rep* 2023;8:837–850. https://doi.org/10.1016/j.ekir.2023.01.008
