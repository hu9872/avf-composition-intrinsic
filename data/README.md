# Data are not stored in this repository

The raw and intermediate data total ~23 GB, which does not belong in git. Everything is
regenerated from public sources by two commands.

## Get the data

```bash
python code/00_download.py     # NCBI GEO FTP → data/raw/   (~23 GB, resumable)
python code/01_prepare.py      # unpack 10x triplets → data/processed/01_raw_concat.h5ad
```

## What the download covers

| Accession | Role |
|---|---|
| GSE305423 | Primary single-cell atlas (this study's main input) |
| GSE220796 | Independent bulk cohort (signature check) |
| GSE119296 | Positive control |
| GSE314264 | Exploratory; not used in the final analysis |
| GSE244317, GSE39488 | Screened and excluded |

Details, sample ledger notes and licence terms: `docs/DATA_SOURCES.md`.

## Not obtainable

Empty-droplet (raw) matrices are not deposited for GSE305423, so SoupX/decontX cannot be run.
See the Limitations section of the manuscript.
