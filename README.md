# Analysis code for: apparent transcriptional signatures of AVF maturation failure are dominated by cell composition

Pre-specified reanalysis of public human arteriovenous fistula (AVF) single-cell RNA-sequencing
data. The study asks a simple question with a methodological answer: when a bulk or
compartment-level comparison reports a "transcriptional signature" of maturation failure, how
much of it is a change in *which cells are present* rather than a change in *what those cells
transcribe*?

Everything needed to reproduce the numbers, figures and tables in the manuscript is here.

---

## What is in this repository

| Path | Contents |
|---|---|
| `code/` | Analysis pipeline, figure builders, document builders, and the numeric-provenance assertion script |
| `results/` | Every result table cited in the article (large intermediate marker tables excluded — see `docs/RESULTS_INDEX.md`) |
| `figures/` | Main figures 1–8, supplementary figures S1–S3 and the graphical abstract (PDF + 300 dpi PNG) |
| `manuscript/` | Manuscript source and the two in-house checking scripts output it is validated against |
| `docs/` | How to reproduce, where the data come from, which file backs which claim |
| `data/` | **Not included** (23 GB). Download instructions in `data/README.md` |

## Quick start

```bash
# 1. environment
conda env create -f environment.yml
conda activate avf-composition-intrinsic

# 2. point the pipeline at your checkout (scripts were developed with an absolute root)
python code/98_relocate_paths.py --to "$PWD"

# 3. get the public data (≈23 GB; needs network and disk)
python code/00_download.py
python code/01_prepare.py

# 4. re-derive results and figures
make all

# 5. check that every number in the manuscript still matches results/
python code/90_verify_numbers.py
```

`make all` runs the pipeline in dependency order. A full run takes several hours on a
16-core workstation with 64 GB RAM; the trajectory and communication steps dominate.

## Reproducibility statements

* **Random seed** — fixed at `2611` for all stochastic steps (Harmony, Leiden clustering,
  Scrublet, downsampling, permutation tests).
* **Unit of replication** — the donor, never the cell. All inference at donor level
  (n = 6 failed, 6 matured).
* **Pinned environment** — `environment.yml` / `requirements.txt`; versions frozen
  2026-11-01.
* **Numeric provenance** — `code/90_verify_numbers.py` parses the manuscript, recomputes each
  quoted value from `results/`, and fails loudly on any mismatch. `results/provenance.csv`
  records, for every claim, the source file, the producing script and the figure.
* **Code / text overlap pre-screen** — `code/95_text_overlap_check.py` compares the manuscript
  against the two papers that originally published these datasets, and includes a positive
  control so that a null result is interpretable.

## What this study does *not* claim

The negative findings are bounded, not universal. At n = 6 per group the design has 80% power
only for standardised effects of about |d| = 1.8 (`results/s30_tost.csv`). "No detectable
difference" is therefore reported as *exclusion of large effects*, not as equivalence. The
same applies to the cDC2 observation, which is reported as a hypothesis.

Ambient-RNA correction with SoupX/decontX was not possible: GEO provides filtered matrices
only, with no empty-droplet profiles. The manuscript states this and substitutes a
marker-source screen plus a rule-sensitivity analysis.

## Data sources

| Accession | Use | Original publication |
|---|---|---|
| [GSE305423](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305423) | Single-cell atlas, 66,168 cells | Martinez et al., *Kidney Int* 2026;109:160–177 ([doi](https://doi.org/10.1016/j.kint.2025.08.024)) |
| [GSE220796](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220796) | Paired bulk vein/AVF, signature check | Martinez et al., *Kidney Int Rep* 2023;8:837–850 ([doi](https://doi.org/10.1016/j.ekir.2023.01.008)) |
| [GSE119296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119296) | Positive control, vein → AVF | see GEO record |

We reanalyse these data; we did not generate them. Attribution to the original depositors is
required if you reuse anything derived here.

## Licence

Code and derived tables: MIT (see `LICENSE`). Source data remain under the terms of their
original GEO depositions. Text and figures from the two original publications are not
redistributed.

## Citation

See `CITATION.cff`. Until the article is accepted, please cite the archived release:

> Tang, Y., & Hu, D. (2026). Analysis code for: apparent transcriptional signatures of
> arteriovenous fistula maturation failure are dominated by cell composition (v1.0.1).
> Zenodo. https://doi.org/10.5281/zenodo.23887196

The concept DOI https://doi.org/10.5281/zenodo.23887195 always resolves to the most recent
version; cite the version DOI above when the exact files matter.
