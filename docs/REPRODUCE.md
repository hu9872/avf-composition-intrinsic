# Step-by-step reproduction

Target: from a clean machine to every number in the manuscript, with no manual editing.

## 0. Prerequisites

- Linux / macOS / Windows (WSL recommended on Windows)
- 16 cores, **64 GB RAM**, 60 GB free disk (23 GB data + intermediates)
- conda or mamba
- Network access to `ftp.ncbi.nlm.nih.gov`

## 1. Environment

```bash
git clone <REPOSITORY URL> avf && cd avf
conda env create -f environment.yml
conda activate avf-composition-intrinsic
```

Expected: solver finishes in 3–10 min; `python -c "import scanpy; print(scanpy.__version__)"`
prints `1.12.5`.

If conda's solver struggles with `pydeseq2`, install it separately:

```bash
pip install pydeseq2==0.5.4
```

## 2. Point the pipeline at your checkout

The scripts were developed with an absolute Windows root (`D:/AVF_GSE_2611`). Rewrite it once:

```bash
python code/98_relocate_paths.py --to "$PWD" --dry-run   # preview
python code/98_relocate_paths.py --to "$PWD"             # apply
```

Expected: 40–46 files changed. This is a pure literal substitution; nothing else is touched.
To undo, run the same command with `--to "D:/AVF_GSE_2611"`.

## 3. Data

```bash
python code/00_download.py     # ~23 GB, 30–90 min depending on link
python code/01_prepare.py      # unpacks 10x triplets → data/processed/01_prepared.h5ad
```

Expected after step 3: `data/raw/` holds the GEO supplementary files;
`data/processed/` holds the AnnData objects and `results/s2_cells_per_sample.csv` shows 12
post-operative samples.

## 4. Analysis

```bash
make analysis     # QC → annotation → composition → trajectories → communication → ECM → bulk
```

| Step | Script | Output | Runtime |
|---|---|---|---|
| QC + doublets | `02_qc.py` | `results/s2_qc_*` | 10 min |
| Annotation | `03_annotate.py` | `results/s3_*` | 25 min |
| Fibroblast subtypes | `04_fibro.py` | `results/s4_*` | 20 min |
| Composition tests | `05b_composition_pseudobulk.py` | `results/s5_composition_tests.csv` (Table 1) | 15 min |
| Compartment pseudobulk | `05c_compartment_pseudobulk.py` | `results/s5c_*` | 40 min |
| Myeloid / cDC2 | `05d_myeloid.py`, `05e_cdc2.py` | `results/s5d_*`, `s5e_*` (Table 2) | 45 min |
| Decomposition | `05f_decomposition.py` | `results/s5f_decomposition_retained.csv` (Table 3) | 25 min |
| Trajectory | `05_trajectory.py` | `results/s7_*` | 90 min |
| Communication | `06_communication.py`, `06c_pathways.py` | `results/s8c_*` | 60 min |
| Matrisome | `07_ecm.py` | `results/s9_*` | 20 min |
| Bulk cohorts | `08_bulk.py`, `08b_vein_to_avf.py` | `results/s10_*`, `s10b_*` | 15 min |

## 5. Robustness analyses

```bash
make robustness
```

| Script | Question it answers | Output |
|---|---|---|
| `30_power_tost.py` | What can n = 6 detect? Is any null an equivalence? | `s30_power_curve.csv`, `s30_tost.csv` |
| `31_rule_floor_sensitivity.py` | Does the marker-source rule or expression floor change the gene list? | `s31_*` |
| `32_gating_perturbation.py` | Does the cDC2 gate (±20%) change the conclusion? | `s32_*` |
| `33_deconvolution.py` | Can bulk data be deconvolved with this reference? | `s33_*` |
| `34_compartment_signatures.py` | Do compartment signatures replicate in bulk? | `s34_*` |
| `35_extra_checks.py` | Donor-level tests for downsampled detection, housekeeping controls | `s35_extra_checks.csv` |

## 6. Figures and tables

```bash
make figures tables
```

Outputs land in `figures/` and `supplement/`. Figures 1–8 and S1 are drawn by `18_figures_final.py`
(the definitive versions; 300 dpi PNG plus vector PDF). Figures S2 and S3 are produced earlier, by
`make robustness` (`30_power_tost.py`, `34_compartment_signatures.py`), and are not overwritten by
`make figures`. The graphical abstract is regenerated at 1200 × 1200 px, 300 dpi by
`40_graphical_abstract.py`.

## 7. Verification (the step that matters)

```bash
python code/90_verify_numbers.py
```

Expected last line:

```
checked 48 numeric claims (144 table cells) against results/ — 0 failures
```

Any failure prints `FAIL <id> <claim>: text <a> vs source <b>` and exits non-zero. This script
is the reason the manuscript's numbers can be trusted: it does not read the manuscript for
values, it recomputes them from `results/` and searches for the exact formatted string.

Text overlap pre-screen (does not need an iThenticate account):

```bash
python code/95_text_overlap_check.py
```

## 8. Common problems

| Symptom | Cause | Fix |
|---|---|---|
| `FileNotFoundError: data/raw/...` | `00_download.py` not run, or run before relocation | rerun step 2 then step 3 |
| `KeyError` on `obs['celltype']` | `03_annotate.py` not run | run `make analysis` from the top |
| Memory kill during `05_trajectory.py` | diffusion maps on 66 k cells | ≥64 GB RAM, or subset to the fibroblast compartment as in `04_fibro.py` |
| `leidenalg` import error | igraph/leidenalg mismatch | `pip install igraph==1.0.0 leidenalg==0.12.0` |
| Figures differ slightly in layout | matplotlib version | pin to `3.11.2` |

## 9. What is deterministic and what is not

Deterministic given the pinned environment: composition tests, pseudobulk DE, power/TOST,
threshold sensitivity, decomposition, figure assembly.

Seeded but hardware-sensitive: Harmony integration, Leiden clustering, diffusion pseudotime,
CellRank. Seeds are fixed at `2611`; expect label numbering to be stable but not bit-identical
across BLAS implementations. If cluster numbering shifts, `results/s3_cluster_labels.json` maps
numbers to identities — always compare by identity, not by index.
