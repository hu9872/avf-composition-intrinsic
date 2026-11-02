# Supplementary material

**Article:** Apparent transcriptional signatures of arteriovenous fistula maturation failure are dominated by cell composition

This folder holds the two supplementary files that accompany the manuscript, plus the individual CSV
versions of every supplementary table.

| File | Content | Built by |
|---|---|---|
| `Document_S1.docx` | Supplemental methods S1–S13, figure legends S1–S3, tables S1–S12 | `code/45_build_document_s1.py` (source `Document_S1.md`) |
| `Data_S1.xlsx` | Source data behind the figure panels, 13 sheets | `code/44_make_supp_tables.py` |
| `csv/S*.csv`, `csv/T_S*.csv` | The same tables as individual files, including unfiltered versions where the workbook holds only significant rows | `code/20_supp_tables.py`, `code/44_make_supp_tables.py` |

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
| `S7_decomposition` | 7 | Composition–intrinsic decomposition of the seven myeloid genes |
| `S8_gene_origin` | 30 | Compartment-resolved genes with cell-of-origin screening outcome |
| `S9_ligand_receptor_axes` | 9 | Pre-specified ligand–receptor axes |
| `S10_ecm_donor` | 12 | Matrisome programme scores per donor |
| `S11_bulk_DE` | 45,799 | Whole-tissue bulk differential expression (GSE220796) |
| `S12_power_TOST` | 14 | Detectable effects and equivalence tests |
| `S13_rule_sensitivity` | 30 | Cell-of-origin rule and expression floor sensitivity |
| `S14_signature_scores` | 8 | Compartment signature scores in GSE220796 |

**Unfiltered tables.** `csv/S3_pseudobulk_DE.csv` and `csv/S5_compartment_DE.csv` contain every tested
gene, not only significant ones.

**Code and intermediate objects.** Analysis scripts (`code/`) and all result tables (`results/`) are
provided so that every number in the manuscript can be traced to its source file;
`code/90_verify_numbers.py` re-derives the quoted values and fails if any of them drifts.
