# Document S1

**Supplemental information for:** Apparent transcriptional signatures of arteriovenous fistula maturation failure are dominated by cell composition

**Contents**

- Supplemental methods S1–S13
- Figure legends for Figures S1–S3
- Tables S1–S12

---

## Supplemental methods

### S1. Data assembly

GSE305423 was downloaded as a tar archive and its MD5 checksum compared with the value published with the series (899,010,560 bytes; md5 ef25e38feb218883dd0fc320e91ff3a9). The 21 libraries come from 20 patients: 6 pre-access cephalic veins, 3 early AVF resections from 2 patients, 6 matured and 6 failed juxta-anastomotic AVFs. The two early libraries from the same patient were merged into one donor, so the object holds 20 donors and not 21; the third early library belongs to a second patient. All three early libraries appear in the atlas but are excluded from outcome comparisons, which rest on the 12 post-operative donors. The pre-access compartment (2,409 cells from 6 donors) supplies the root for pseudotime and takes no part in any group test.

GSE220796 was obtained as a count matrix with sample annotation. Of 41 libraries, 20 patients had both a basilic vein and a juxta-anastomotic AVF sample; one further failed AVF had no paired vein and was excluded, leaving 9 matured and 11 failed pairs. GSE119296 has 19 vein and 19 AVF libraries, but only one vein library could be matched to its AVF counterpart by identifier. That contrast is between patient groups rather than within patients, so the dataset is used only as a positive control for instrument sensitivity.

### S2. Quality control, integration and compartment annotation

Filters were applied in this order: at least 200 detected genes, at least 500 UMIs, mitochondrial reads below 15%, and fewer than 6,000 detected genes. This retained 74,285 of 97,070 cells (76.5%). Doublets were then called with Scanpy's Scrublet implementation run per library (batch_key = sample), giving 497 calls that were removed, leaving 73,788 cells.

Normalisation was to 10,000 counts per cell followed by natural-log transformation. Highly variable genes were selected with the Seurat flavour on the donor batch (2,000 genes), principal components were computed with zero-centring disabled (50 components), and libraries were integrated with Harmony on donor (maximum 20 iterations). The neighbour graph used 15 neighbours and 40 principal components.

Leiden clustering at resolution 0.8 gave 40 clusters. For each cluster, rank-sum markers were computed against all other cells (Wilcoxon test, log2FC above 0.5, FDR below 0.05) and the top 50 retained. Eight canonical panels were defined: fibroblast (DCN, LUM, COL1A1, COL3A1, PDGFRA), endothelial (PECAM1, VWF, CDH5, KDR), smooth muscle (MYH11, ACTA2, TAGLN, CNN1), pericyte (RGS5, KCNJ8, PDGFRB), myeloid (CD68, CD14, LYZ, C1QA), T (CD3D, CD3E, IL7R, CD2), NK (NKG7, GNLY, KLRD1) and B/plasma (MS4A1, CD79A, IGHG1). A cluster was assigned to the panel contributing the largest number of its top 50 markers, normalised by panel length; clusters with fewer than three supporting markers were labelled "Mixed" and removed.

Erythrocyte contamination was scored as the mean of HBB, HBA1 and HBA2. Contaminated clusters sat at 3.13–7.28 log-normalised units and clean clusters at 0.04–0.36, so the 1.0 cut-off falls in an empty interval. Together the two rules removed 7,620 cells (10.3%), leaving 66,168 in eight compartments.

### S3. Mural and fibroblast subclustering

The 30,009 cells of the fibroblast, smooth muscle and pericyte compartments were re-embedded on 1,500 highly variable genes (Seurat flavour, donor batch) and 40 principal components, integrated with Harmony on donor (20 iterations) and graphed with 15 neighbours and 30 components. Leiden was run at 0.6 and at 1.0. Resolution 0.6 merged contractile and synthetic smooth muscle into one cluster; 1.0 separated them without splitting the fibroblast states, and was used. This gave 24 subclusters, collapsed into six states by the panels below, taking the panel with the highest mean standardised expression per subcluster. Marker support (the number of panel genes among the subcluster's top 50 rank-sum markers) is reported per subcluster in Table S2.

| State | Panel genes |
|---|---|
| Resident adventitial-like fibroblast | PI16, ABCA8, MFAP5, GSN, TPPP3, HTRA3, SERPINF1 |
| Stromal fibroblast | DCN, LUM, COL1A1, COL1A2, COL3A1, BGN, ASPN, C7, CFD |
| Myofibroblast | POSTN, COL8A1, FN1, CTHRC1, TGFBI, SPP1, COMP, TIMP1, COL5A1, COL10A1 |
| Pericyte | RGS5, KCNJ8, ABCC9, NOTCH3, HIGD1B, COX4I2 |
| Contractile SMC | MYH11, CNN1, MYOCD, LMOD1, DES, ACTG2 |
| Synthetic SMC | TAGLN, MYL9, TPM2, ELN, COL4A1, COL4A2 |

Programme scores used scanpy.tl.score_genes (control set 50 genes, 25 bins). Seven programmes were scored: contractile, synthetic myofibroblast, adventitial, TGF-β response, interferon response, matrisome and crosslinking.

### S4. Myeloid subclustering and gated annotation

Immunoglobulin genes (symbols beginning with IG) and haemoglobin genes (HBB, HBA1, HBA2, ALAS2) were removed before embedding. The remaining cells were embedded on 1,500 highly variable genes, 40 principal components and Harmony on donor, graphed with 15 neighbours and 30 components, and clustered at Leiden 0.9, giving 22 clusters. A marker-support rule of the kind used at atlas level left 76% of cells unassigned and recovered no cDC2 cluster, because rank-sum markers in this compartment are dominated by ribosomal (RPL/RPS) and cytoskeletal transcripts.

Annotation was therefore done by a gated rule on cluster-mean expression, applied in priority order, with the first satisfied rule deciding:

| Priority | State | Rule (cluster mean, log-normalised units) |
|---|---|---|
| 1 | Mast cell | TPSAB1 ≥ 0.5 or CPA3 ≥ 0.5 |
| 2 | Neutrophil | S100A8 ≥ 1.5 with CD14 < 0.5 and C1QA < 0.5 |
| 3 | cDC2 | CD1C ≥ 0.4, or CLEC10A ≥ 0.6 with C1QA < 1.6 |
| 4 | cDC1 | CLEC9A ≥ 0.3 or XCR1 ≥ 0.2 |
| 5 | pDC | LILRA4 ≥ 0.15 with CLEC4C ≥ 0.08 |
| 6 | Classical monocyte | C1QA ≥ 1.0 with CD14 ≥ 1.5, C1QA < 2.2 and FCN1 ≥ 0.15 |
| 7 | Macrophage | C1QA ≥ 1.0 (and rule 6 not met) |
| 8 | Classical monocyte | CD14 ≥ 1.0 |
| 9 | FCGR3A monocyte | FCGR3A ≥ 1.5 |
| 10 | Low quality | none of the above |

The thresholds are permissive relative to published cDC2 marker levels; their function is to prevent assignment to a lineage whose defining gene is absent, not to define a lineage by expression magnitude. Applying them gave 11,937 macrophages, 1,578 classical monocytes, 1,158 cDC2, 73 pDC, 10 neutrophils and 89 low-quality cells. No cluster satisfied the cDC1 rule: CLEC9A was detected in 6.4% of myeloid cells at a mean of 0.12 counts and XCR1 in 0.4%, with no cluster mean above 0.09 for either gene. No mast-cell cluster survived because mast cells were partitioned at atlas level, before this subclustering.

Two sensitivity analyses were applied to the gating. First, a cluster-independent definition classified a cell as cDC2 when at least two of CD1C, CLEC10A, FCER1A and CD1E were detected, without reference to clustering. Second, every numeric threshold in the rule above was multiplied by 0.8 or by 1.2 and the whole annotation recomputed, then cDC2 counts per 1,000 tissue cells were retested (Table S7).

### S5. Depth control within cDC2

Failed-group cDC2 cells had a median of 14,474 UMIs against 17,070 in matured fistulas, a 15% difference that lowers the detection probability of every transcript and would on its own produce an apparent loss of markers. Each cell was therefore downsampled to exactly 5,000 UMIs by multinomial resampling of its count vector, and detection rate and counts per 10,000 recomputed on the resampled matrix. Cells with fewer than 5,000 UMIs were excluded; 1,086 of 1,158 cDC2 cells qualified, from 20 donors. The donor-level comparison used the 12 post-operative donors (6 failed, 6 matured).

Seven housekeeping genes served as positive controls. Standardised differences (failed minus matured, counts per 10,000) were RPS18 0.17, EEF1A1 0.12, GAPDH 0.31, TMSB10 −0.33, B2M −0.66, ACTB −0.97 and RPL13A 1.04. Four of the seven had p above 0.5; RPL13A (p = 0.11) and ACTB (p = 0.14) did not reach significance, so equalisation was partial but without a detectable residual difference in the controls. This is reported in the main text and constrains how the CD1C difference can be read.

A CD1C-low, C1QA-positive state was defined by re-embedding the dendritic-like cells (cDC2, pDC and the CD1C-negative monocyte fraction) on 800 highly variable genes, 25 principal components, Harmony on donor (15 iterations), 15 neighbours on 20 components and Leiden 0.8. Its frequency was expressed as a centred log-ratio per donor.

### S6. Cell-of-origin screening and expression floor

For each compartment-level differentially expressed gene, the mean log-normalised expression was computed in all eight compartments. A gene was retained if its mean in the target compartment reached half of the maximum across compartments, or 0.15 units, whichever was more permissive. Retained genes were then dichotomised by a uniform 0.10 floor: genes above it were interpreted mechanistically, genes below it were reported and carried through the decomposition arithmetic without mechanistic weight. The origin rule was subsequently perturbed to 25%, 50% and 80% of the maximum, applied both with "or" and with "and" against the 0.15 floor, and the floor was varied at 0.05, 0.10 and 0.15 (Tables S9 and S10).

### S7. Composition–intrinsic decomposition

For each compartment-level gene the observed log2 fold change was written as the sum of two terms:

log2FC(observed) ≈ log2(f_failed / f_mature) + log2(rate_failed / rate_mature)

where f is the share of the compartment occupied by the gene's dominant state and rate the mean expression of that gene within that state. The first term is the composition term, the second the cell-intrinsic term, and the residual measures how far the two-term model is from the observation. Observed values are ratios of pseudobulk mean rates rather than shrunken DESeq2 estimates, so they differ marginally from the gene-level fold changes; both are reported. The identity holds when one state dominates expression of the gene, which is why a broadly expressed gene (YEATS4) fits poorly.

### S8. Trajectory

PAGA was run on the six-state graph built from the 15-neighbour, 30-component graph described above, and the layout stored for plotting. Diffusion maps used 15 components; diffusion pseudotime used 10. The root was the medoid of the pre-access resident adventitial-like cells in diffusion space: of 795 candidate cells, the one with the smallest Euclidean distance to the compartment centroid. Taking instead the extreme of the first diffusion component placed the root at the smooth-muscle end, because that extreme is determined by a few outlier cells, and inverted the axis; the resulting state ordering was checked against the marker-based expectation before that option was rejected.

The transition matrix for fate absorption was the connectivities graph, oriented so that pseudotime increases along every retained edge. Cells with no forward neighbour were made absorbing. Three terminal sets were the upper pseudotime quartile of the myofibroblast, contractile SMC and resident adventitial-like states. A leakage of 1 × 10⁻³ was applied so that the linear system stays non-singular, with the consequence that absorption probabilities sum to slightly below 1 and are interpreted as relative fate bias. Absorption was solved exactly and averaged per donor.

RNA velocity was not estimated because the deposited matrices are filtered (cellranger filtered_feature_bc_matrix), which removes the intronic reads a reliable velocity estimate requires.

### S9. Ligand–receptor analysis

Thirty-four pairs were specified in advance from CellPhoneDB v2: four TGF-β, four PDGF, three Notch, five chemokine, four interleukin and TNF, seven matrix–integrin, and seven further pairs spanning ephrin, angiopoietin, VEGF, Wnt, serpin and HGF signalling. A pair was scored as the product of the sending state's ligand mean and the receiving state's receptor mean, restricted to cells in which both partners were detected in at least 10% of the population (the CellPhoneDB default). Scores were averaged per donor and compared by Welch's t-test.

The donor-level design admits only 924 distinct partitions of 12 donors into two groups of 6, so the smallest attainable two-sided p is about 0.0011 and no interaction can survive FDR control across 744 tests. The pathway-level analysis collapses the 744 interactions into nine axes and is the pre-specified inferential unit.

### S10. ECM programmes and bulk cohorts

Matrisome genes were the 372 genes annotated to GO:0062023, cross-checked against the core matrisome of Naba et al. Non-negative matrix factorisation was run on the log-normalised, z-scored matrisome matrix for k = 3 to 7, with NNDSVD-a initialisation and 800 iterations. Reconstruction error fell from 1,696 at k = 3 to 1,576 at k = 7, improving by less than 2% per additional factor beyond k = 5; k = 5 also reproduced the categories expected from the literature (interstitial proteoglycan, basement membrane, microfibrillar, fibrillar collagen, fibronectin/versican) and was retained. Zonation was defined as the media score minus the adventitia score over layer-specific gene sets. It carries no spatial information because the dataset has no coordinates.

Bulk counts were analysed with PyDESeq2 (design ~outcome). For GSE119296, Ensembl identifiers were mapped to HGNC symbols and only unambiguous one-to-one mappings retained. Programme scores were the mean of gene-wise z-scores computed within cohort, so a score of zero corresponds to the cohort mean.

### S11. Power and equivalence analysis

Power was computed for a two-sample t-test with equal group sizes by non-central t distribution, at a two-sided α of 0.05. The detectable standardised effect at 80% power was solved numerically for n from 3 to 40; it is 1.80 at n = 6 per group, 1.20 at n = 12 and 0.79 at n = 26 (Table S5). Equivalence was assessed by two one-sided tests against a smallest effect of interest of d = 0.80, which corresponds to 0.80 times the pooled standard deviation of each endpoint; the pooled standard deviation was recovered from the reported difference and Cohen's d, and the 90% confidence interval reported alongside. All 14 endpoints tested failed to establish equivalence at this margin, which is the expected outcome for a design of this size and is why no null result in this paper is described as demonstrating equivalence.

### S12. Compartment signature scores in an independent bulk cohort

Deconvolution of the bulk cohort against a same-tissue single-cell reference was attempted first and was not identifiable: with 11,219 shared genes and unconstrained non-negative least squares, weights collapsed onto a single compartment and the median coefficient of determination was 0.31. Compartment abundance was therefore probed with signature scores. For each of the eight compartments, the 50 genes with the highest specificity in the study's own reference profile were selected, specificity being log2 CPM in the compartment minus the mean over the other compartments, restricted to genes above 10 CPM and to positive specificity. Each sample was scored as the mean gene-wise z-score across samples, and a housekeeping-normalised variant was computed using the mean of RPL13A, RPLP0, B2M, TBP, GAPDH, ACTB and PPIA. Scores were compared between outcomes by Welch's t-test. Sample-level global properties (mean log2 CPM, number of detected genes, library size) did not differ between groups (all p > 0.9), which excludes a library-quality explanation for the observed differences.

### S13. Software and reproducibility

Python 3.13.14; scanpy 1.12.5; anndata 0.13.4; numpy 2.5.4; pandas 3.1.0; scipy 1.18.1; scikit-learn 1.9.1; statsmodels 0.15.0; PyDESeq2 0.5.4; harmonypy 0.0.9; leidenalg 0.12.0; igraph 1.0.0; matplotlib 3.11.2. Random seed 2611 throughout. Scripts are numbered in execution order and results are written one table per step, so that each number in the main text traces to a file. The environment is pinned in environment.yml in the code archive.

---

## Figure legends

[[FIG:FigS1_compartment_genes.png]]

**Figure S1. All 30 compartment-resolved genes.** Each gene is plotted with its log2 fold change (failed versus matured) and its mean expression in the target compartment, coloured by retention after cell-of-origin screening and by whether it exceeds the 0.10 expression floor. Genes below the floor are shown open. Related to Figure 6.

[[FIG:FigS2_power.png]]

**Figure S2. What the design can detect.** (A) Power as a function of standardised effect size for 6, 12 and 26 donors per group, computed for a two-sided two-sample t-test at α = 0.05. The dashed line marks 80% power at n = 6, reached at d = 1.80. (B) Difference between failed and matured fistulas for 14 endpoints, with 90% confidence intervals, against equivalence bounds of d = 0.80 (shaded). No endpoint met the equivalence criterion. Related to Figures 3, 4 and 6.

[[FIG:FigS3_compartment_signatures.png]]

**Figure S3. Compartment signature scores in an independent bulk cohort.** Per-sample signature score (mean gene-wise z) for eight compartments in GSE220796 (11 failed, 9 matured); horizontal bars are group means. The myeloid and endothelial signatures were lower in failed fistulas and the SMC signature higher. Related to Figure 7 and to the Limitations of study.

---

## Tables

**Table S1. Cells per donor and compartment.** Cell counts after quality control, doublet removal and compartment assignment, for all 20 donors. Related to Figure 1.

[[TABLE:T_S1]]

**Table S2. Mural and fibroblast subclusters and their marker support.** Cluster identity, cell number, classification margin, and the number of panel genes among each subcluster's top 50 rank-sum markers. Related to Figure 2.

[[TABLE:T_S2]]

**Table S3. Differentially expressed genes in the fibroblast lineage.** Donor-level pseudobulk, failed versus matured, FDR below 0.05 and absolute log2 fold change above 0.585. Related to Figure 2C.

[[TABLE:T_S3]]

**Table S4. Fate absorption probabilities by donor.** Multi-fate absorption for the three terminal states, with donor-level tests. Related to Figure 4C.

[[TABLE:T_S4]]

**Table S5. Detectable effects and equivalence tests.** For each of 14 endpoints: difference, Cohen's d, Welch p, the equivalence bound corresponding to a smallest effect of interest of d = 0.80, the 90% confidence interval, the two one-sided test p value, whether equivalence was established, and the smallest standardised effect the data exclude. At n = 6 per group the design reaches 80% power only at d = 1.80. Related to Figures 3, 4, 6 and S2.

[[TABLE:T_S5]]

**Table S6. Compartment-level differential expression.** Number of genes tested and number reaching threshold in each compartment, with the number up and down. Related to Figure 6A.

[[TABLE:T_S6]]

**Table S7. Gating threshold perturbation.** cDC2 annotation recomputed with every numeric threshold in the gating rule multiplied by 0.8, 1.0 or 1.2, and the donor-level comparison of cDC2 per 1,000 tissue cells repeated. Related to Figure 7.

[[TABLE:T_S7]]

**Table S8. CD1C detection within cDC2 after depth equalisation.** Per-donor detection rate for CD1C and three further cDC2 markers, plus three housekeeping controls, after every cell was downsampled to 5,000 UMIs. Related to Figure 8B.

[[TABLE:T_S8]]

**Table S9. Cell-of-origin rule sensitivity.** All 30 compartment-resolved genes under four variants of the origin rule. Cell means are log-normalised units. Related to Figure 6 and Figure S1.

[[TABLE:T_S9]]

**Table S10. Expression floor sensitivity.** Retained genes classified as interpretable or not at floors of 0.05, 0.10 and 0.15 log-normalised units. Related to Figure S1.

[[TABLE:T_S10]]

**Table S11. Compartment signature scores in GSE220796.** Mean signature score per compartment by outcome, with donor-level tests. Related to Figure S3.

[[TABLE:T_S11]]

**Table S12. Donor characteristics.** Sampling group, sex, age, myeloid cell number and total cells for all 20 donors. Related to the Limitations of study.

[[TABLE:T_S12]]

**Data S1. Source data for Figures 1–8 and S1–S3.** Per-donor cell counts by compartment and myeloid subtype, mural state fractions, composition tests, fibroblast pseudobulk results, decomposition terms, gene-of-origin screen, ligand–receptor axis tests, ECM programme scores, bulk differential expression, power and equivalence output, rule sensitivity and signature scores, one sheet per analysis step.
