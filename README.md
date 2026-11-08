# Analysis code for: Cell composition accounts for much of the transcriptomic difference between failed and matured human arteriovenous fistulas

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

The manuscript text is not part of this repository and not part of any release. The deposit is
the analysis code; the prose is submitted to the journal separately. A working copy may sit in
`manuscript/` on a local machine, where `.gitignore` keeps it out of version control and
`.gitattributes` keeps it out of every archive.

**Figure files and article figure numbers.** The file stems predate the figure order used in
the submitted article, so they do not match it one to one:

| File stem | Article figure |
|---|---|
| `Fig1_atlas` | Figure 1 |
| `Fig2_subclusters` | Figure 2 |
| `Fig3_composition` | Figure 3 |
| `Fig4_trajectory` | Figure 4 |
| `Fig6_compartments` | Figure 5 |
| `Fig7_myeloid_cdc2` | Figure 6 |
| `Fig8_signal_decomposition` | Figure 7 |
| `Fig5_communication` | Figure 8 |
| `docs/` | How to reproduce, where the data come from, which file backs which claim |
| `supplement/` | `Data_S1.xlsx` (13 worksheets, S1–S13), `Document_S1.md` (Supplemental methods, figures and tables), and the per-step CSV exports in `csv/` — see `supplement/README.md` for how the three S-numbering series differ |
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

# 5. re-derive every number from results/ (add MS_PATH to also match the prose)
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
* **Numeric provenance** — `code/90_verify_numbers.py` recomputes every quoted value from
  `results/` and, when a manuscript is available, asserts that the value appears in the text;
  it fails loudly on any mismatch. Because the prose is not shipped, a fresh clone runs in
  recompute-only mode and reports each claim as `RECOMPUTED`. `results/provenance.csv` records,
  for every claim, the source file, the producing script and the figure.
* **Code / text overlap pre-screen** — `code/95_text_overlap_check.py` compares a local
  manuscript (`--manuscript`, default: the most recent `manuscript/Manuscript_*.md`) against
  the two papers that originally published these datasets, and includes a positive control so
  that a null result is interpretable. It exits cleanly when no manuscript is present.

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

> Tang, Y., & Hu, D. (2026). Analysis code for: Cell composition accounts for much of the
> transcriptomic difference between failed and matured human arteriovenous fistulas (v1.0.1).
> Zenodo. https://doi.org/10.5281/zenodo.23887196

The concept DOI https://doi.org/10.5281/zenodo.23887195 always resolves to the most recent
version; cite the version DOI above when the exact files matter.
