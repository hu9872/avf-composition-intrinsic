# Supplementary material

**Article:** Cell composition accounts for much of the transcriptomic difference between failed and matured human arteriovenous fistulas

This folder holds the supplementary files that accompany the manuscript, plus the individual CSV
versions of every supplementary table.

| File | Content | Built by |
|---|---|---|
| `Document_S1.docx` | Supplemental methods S1–S13, figure legends S1–S3, tables S1–S12 | `code/45_build_document_s1.py` (source `Document_S1.md`) |
| `Data_S1.xlsx` | Source data behind the figure panels, 13 sheets | `code/44_make_supp_tables.py` |
| `csv/T_S*.csv` | Tables S1–S12 of Document S1, one file per table | `code/44_make_supp_tables.py` |
| `csv/S*.csv` | Intermediate exports from individual pipeline steps, including unfiltered versions where the workbook holds only significant rows. The number is the step order, not a worksheet number — see the note below | `code/20_supp_tables.py` |

**Two entries in that table are build products, not archived binaries.** `Document_S1.docx` does not
exist in this archive as a file: it is generated from `Document_S1.md` by
`code/45_build_document_s1.py`. Likewise two CSV exports (`csv/S3_pseudobulk_DE.csv`,
`csv/S5_compartment_DE.csv`) exceed the archive size limit and are written by `code/20_supp_tables.py`.
Run `make tables` to recreate all three; the files the journal itself needs are submitted separately.

**Three numbering series, deliberately distinct.** This archive carries three independent S-number
series. (i) *Supplemental methods S1–S13* — the method sections of Document S1. (ii) *Data S1
worksheets S1–S13* — the source-data workbook, named after the analysis step that produced each sheet.
(iii) *`csv/S*.csv`* — intermediate exports, numbered by pipeline step order (S4b, S6b and S7b are
sub-steps). They are not aligned with one another: for example `csv/S6_myeloid_markers.csv` is a
marker table from step 6, while Data S1 sheet `S6_decomposition` holds the composition–intrinsic
decomposition, and Document S1 § S6 is the cell-of-origin screen. Read the file names, not the
numbers, when crossing between the three.

The three supplementary figures are not duplicated here. They live in `../figures/`
(`FigS1_compartment_genes`, `FigS2_power`, `FigS3_compartment_signatures`; PDF and 300 dpi PNG) and are
embedded into `Document_S1.docx` at build time from those files.

## Document S1

`Supplemental methods` covers, in order: data assembly (S1); quality control, integration and
compartment annotation (S2); mural and fibroblast subclustering (S3); myeloid subclustering and gated
annotation (S4); depth control within cDC2 (S5); cell-of-origin screening and the expression floor
(S6); the composition–intrinsic decomposition (S7); trajectory analysis (S8); ligand–receptor analysis
(S9); ECM programmes and the bulk cohorts (S10); the power and equivalence analysis (S11); compartment
signature scores in an independent bulk cohort (S12); software and reproducibility (S13).

Every value in the methods was read from the analysis scripts or from `results/`, not reconstructed
from memory, so the pipeline can be re-run without the original author. Two tables sit inside the text
(state marker panels; the ten gating rules in priority order) and are unnumbered on purpose — they
belong to the method description, not to the data.

Tables S1–S12 map to figures as follows: S1 → Figure 1; S2, S3 → Figure 2; S4 → Figure 4C;
S5 → Figures 3, 4, 6 and S2; S6 → Figure 6A; S7 → Figure 7; S8 → Figure 8B; S9, S10 → Figure 6 and
Figure S1; S11 → Figure S3; S12 → donor characteristics quoted in the Limitations.

## Data S1

Thirteen sheets of source data, named after the analysis step that produced them:

| Sheet | Rows | Content |
|---|---|---|
| `S1_compartment_counts` | 21 | Cells per donor and compartment |
| `S2_myeloid_counts` | 20 | Myeloid state counts per donor |
| `S3_mural_fractions` | 12 | Mural/fibroblast state fractions per donor |
| `S4_composition_tests` | 6 | Donor-level composition tests |
| `S5_fibroblast_DE` | 12 | Differentially expressed genes in the fibroblast lineage |
| `S6_decomposition` | 7 | Composition–intrinsic decomposition of the seven myeloid genes |
| `S7_gene_origin` | 30 | Compartment-resolved genes with cell-of-origin screening outcome |
| `S8_ligand_receptor_axes` | 9 | Pre-specified ligand–receptor axes |
| `S9_ecm_donor` | 12 | Matrisome programme scores per donor |
| `S10_bulk_DE` | 45,799 | Whole-tissue bulk differential expression (GSE220796) |
| `S11_power_TOST` | 14 | Detectable effects and equivalence tests |
| `S12_rule_sensitivity` | 30 | Cell-of-origin rule and expression floor sensitivity |
| `S13_signature_scores` | 8 | Compartment signature scores in GSE220796 |

**Unfiltered tables.** `csv/S3_pseudobulk_DE.csv` and `csv/S5_compartment_DE.csv` contain every tested
gene, not only significant ones.

**Code and intermediate objects.** Analysis scripts (`code/`) and all result tables (`results/`) are
provided so that every number in the manuscript can be traced to its source file;
`code/90_verify_numbers.py` re-derives the quoted values and fails if any of them drifts.
