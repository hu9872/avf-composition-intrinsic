# Apparent transcriptional signatures of arteriovenous fistula maturation failure are dominated by cell composition

**Authors**

Yuewu Tang^1, Daijun Hu^1,^*,^†

^1 Department of Nephrology, Chongqing University Three Gorges Hospital, 165 Xincheng Road, Wanzhou District, Chongqing 404000, China

^† Lead contact
^\* Correspondence: Daijun Hu, 283624940@qq.com

**Summary**

Maturation failure limits autologous arteriovenous fistulas (AVF), and bulk studies attribute it to myofibroblast expansion and a shift in the extracellular matrix programme. Such comparisons mix cell abundance with per-cell transcription. We reanalysed public human AVF single-cell data (66,168 cells; six failed and six matured fistulas), testing three pre-specified hypotheses with the donor as the unit of replication. Myofibroblast abundance, fate bias and all five matrisome programmes were unchanged. Compartment-level differences were largely explained by ambient RNA and cell composition: myeloid cells were 1.8-fold more abundant in failed fistulas, whereas cDC2 numbers were unchanged, so cDC2 transcripts were diluted rather than lost. A composition–intrinsic decomposition attributed a median 51% of the remaining signal to composition. In an independent bulk cohort, however, myeloid signature scores were not higher in failed fistulas. Transcriptional signatures of AVF failure largely reflect tissue composition, and gene-level biomarkers should be read against it.

**Keywords**

arteriovenous fistula; vascular access; maturation failure; single-cell RNA sequencing; cell composition; ambient RNA

## Introduction

Autologous arteriovenous fistulas (AVF) remain the preferred access for haemodialysis, but maturation often fails: many never reach a usable diameter or flow, and many that do need later intervention for stenosis [[1]],[[2]]. Below this sits venous neointimal hyperplasia, built from proliferating, migrating cells and excess extracellular matrix (ECM), set against the flow-driven remodelling that decides whether the vein dilates or narrows [[3]],[[4]].

Work on the lesion has rested on two ideas. One is cellular: adventitial fibroblasts activate, migrate inward and take on a myofibroblast phenotype, the "outside-in" sequence [[5]],[[6]],[[7]]. The other is molecular: a failing fistula is thought to differ from a maturing one in the magnitude of ECM programme activation, so that failure is repair that overshoots. Both point therapeutically towards anti-fibrotic and anti-myofibroblast treatment.

The human evidence is indirect. Bulk transcriptomics of paired pre-access vein and postoperative AVF showed strong matrisome activation after anastomosis, with proteoglycans and collagens distributed by wall layer [[8]]; later analyses of the same material ranked candidate genes and built classifiers [[9]]. The limitation of these designs is rarely stated: a bulk or compartment-level comparison cannot separate a change in how many cells of a state are present from a change in what those cells transcribe. A transcript confined to a rare population therefore looks down-regulated whenever that population is diluted, even if every surviving cell expresses it normally. Composition varies several-fold between donors, and dissociation releases message from abundant bystander cells [[10]],[[11]].

Single-cell RNA sequencing of human AVF tissue is now available for 20 patients, sampled at access creation, early after anastomosis and again at second-stage surgery [[12]]. Its authors read it as support for both ideas: persistent inflammation and an exaggerated healing response in failed fistulas, proinflammatory and adhesive macrophages with an activated Myddosome programme, activated vasa vasorum endothelium, hyperactive fibroblasts and myofibroblasts, and a macrophage-derived osteopontin signal on the surrounding wall, with morphometry and immunofluorescence documenting the wall expansion that accompanies it. They also noted that the global transcriptional difference between outcomes is subtle. That reading concerns what cells are doing. The same libraries allow a narrower question, how much of each difference is how many cells and how much is what those cells transcribe, and it gives a different answer.

We asked whether a synthetic myofibroblast state is more abundant, and a resident adventitial-like state less abundant, in failed fistulas (i); whether fibroblast-lineage cells sit further along an activation trajectory (ii); and whether ECM programmes shift in composition rather than uniformly increase (iii). All three were pre-specified and tested with the donor as the unit of replication.

None of the three held. Two further steps were needed to see why: compartment-level pseudobulk localised what remained, and decomposition split each difference into a composition term and a cell-intrinsic one. A negative result is worth only as much as the instrument behind it [[13]], so we quantified what the design can and cannot detect, tested the robustness of every threshold that could drive the result, and applied the same instrument to a contrast in which large remodelling is known to occur.

## Results

### An eight-compartment atlas and six mural states

Quality control and doublet removal left 73,788 cells; 7,620 (10.3%) lacked a coherent identity or carried erythrocyte contamination, leaving 66,168 cells in eight compartments (Figure 1A; Table S1). Fibroblasts were the largest (24,106 cells, 36.4%), myeloid second (14,845) and endothelium third (9,904); every sampling stage contributed to every compartment (Figure 1B–C).

Subclustering the 30,009 mural and fibroblast cells gave six states (Figure 2A): resident adventitial-like fibroblasts (10,298), myofibroblasts (7,809), stromal fibroblasts (4,855), contractile smooth muscle cells (SMC; 3,762), pericytes (2,730) and synthetic SMC (555). Programme scores returned the expected arrangement, except that stromal fibroblasts and not myofibroblasts carried the highest matrisome score (Figure 2B). Both outcome groups appeared in every state.

### Myofibroblast abundance is unchanged, and the design bounds what can be excluded

The myofibroblast fraction was nearly identical (ΔCLR +0.12, d = 0.12, p = 0.84), as was the resident adventitial-like fraction (+0.21, d = 0.19, p = 0.75) (Figure 3; Table 1). Pericytes differed most, lower in failed fistulas (−0.94, d = −1.19, Welch p = 0.069, Mann–Whitney p = 0.041), without surviving correction (FDR = 0.41).

The donor had to be the unit of analysis. Between-donor spread exceeded any group effect, the resident fraction running from 4.5% to 72.7%, and the synthetic SMC state, which at cell level looked confined to failed fistulas (527 versus 28 cells), came from one donor contributing 517 of them; only 4 of 6 failed and 3 of 6 matured donors contained any, and the donor-level test was null (+0.69 CLR, d = 0.46, p = 0.46).

A null test at this size is a statement about resolution, not about equivalence, so we quantified both. With six donors per group and a two-sided test at α = 0.05, 80% power is reached only at |d| ≈ 1.80 (Figure S2A; Table S5). Dual one-sided testing against a smallest effect of interest of d = 0.80 did not establish equivalence for any of the 14 endpoints examined (all TOST p > 0.09; Figure S2B; Table S5), and 90% confidence intervals bounded the myofibroblast difference to between −0.90 and +1.14 CLR units, excluding standardised effects above about 1.2. These data therefore exclude a large redistribution of mural states; they do not exclude a modest one.

### The fibroblast transcriptome is close to unchanged

Donor-level pseudobulk within the fibroblast lineage returned 12 differentially expressed genes from 24,157 tested, and within states none for resident or stromal fibroblasts and one each for myofibroblasts (MT1X), pericytes (NR4A2) and contractile SMC (GREM2) (Figure 2C; Table S3). The genes that motivated the study, POSTN, COL8A1, FN1, TIMP1 and PI16, were below threshold in every state.

The twelve genes are informative for what they are not. None belongs to the matrisome; several mark lymphocytes or lymphatic endothelium (XCL2 −7.14, CCL21 −3.03), others low-abundance regulators (CDKN2A +2.02, CD200R1 −2.54). A coordinated activation programme would move collagens, glycoproteins and crosslinking enzymes together. That pattern is missing.

### The activation trajectory is a gradient, and its use is not displaced

Absorption into the myofibroblast fate was 0.117 in failed and 0.171 in matured fistulas (d = −0.75, p = 0.23, FDR = 0.49); into contractile SMC 0.278 versus 0.228 (d = 0.42, p = 0.49); into resident adventitial-like fibroblasts 0.352 versus 0.280 (d = 0.44, p = 0.47) (Figure 4C; Table S4). No fate survived correction, and the myofibroblast point estimate ran opposite to the hypothesis. The 90% confidence interval for the myofibroblast displacement (−0.129 to +0.021) excludes shifts larger than about 1.8 times the between-donor standard deviation, so a gross redistribution is ruled out and a modest one is not.

Pseudotime ordered the states as expected, though not as stages. Median percentile was 0.18 in resident adventitial-like fibroblasts, 0.48 in stromal fibroblasts, 0.56 in myofibroblasts, and 0.83 and 0.84 in contractile and synthetic SMC (Figure 4B, 4D). The two intermediate states share a band rather than following one another, which fits their sharing most programmes (Figure 2B). PAGA gave the same architecture: one connected fibroblast module joined to the SMC module through myofibroblasts and stromal cells, with pericytes attached to contractile SMC (Figure 4A).

### Matrisome programmes show no compositional shift

Non-negative matrix factorisation of the 372 matrisome genes gave five programmes: interstitial proteoglycan (SFRP2, DCN, LUM, MGP, SERPINF1), basement membrane (COL4A2, COL6A2, TIMP3, IGFBP7), microfibrillar (FBN1, MFAP5, TNXB, CLU), fibrillar collagen (COL1A1, COL3A1, COL5A1/2, SERPINH1) and fibronectin/versican (FN1, VCAN, CCN2, ELN). None differed between groups (p = 0.29–0.99), and the zonation index was flat (Δ = 0.015, d = 0.01, p = 0.98) (Figure 6B).

Interpretation of these five nulls carries the same caveat as the abundance nulls. The narrowest 90% confidence interval among the five programmes excluded standardised effects above 1.05, and the widest only above 1.69 (Table S5), so a coordinated programme shift of moderate size would not have been detected.

### Compartment-level differences are largely ambient or compositional

Whole-tissue bulk of the same patients had reported hundreds of differentially expressed genes [[8]], so we asked which compartment carries them. The same model per compartment produced a steep gradient (Figure 6A; Table S6): 16 genes in myeloid cells, 7 in endothelium, 2 in mast cells, 2 in T cells and 0–1 per fibroblast state. Read plainly, that places the biology in myeloid and endothelial cells.

The gene-of-origin screen says otherwise. Nine of the 16 myeloid genes were not myeloid. MAL, the most strongly down-regulated at −5.91 log2, had a mean of 0.22 in T cells and 0.00 in every myeloid state; its apparent fall was made of ambient T-cell message, plentiful because T cells are 10.5% of the tissue. IL32 behaved the same way (NK 1.81, T 1.63, myeloid 0.07), as did six further myeloid genes and two endothelial ones. IL1RL1, nominally up 4.34-fold, had a mean of 1.23 in mast cells against 0.09 in endothelium.

Five endothelial and seven myeloid genes passed the origin rule, but a fold change on a transcript expressed at 0.01–0.06 units is mostly sampling noise. Under the 0.10 floor, two endothelial genes remain interpretable (ACE −1.70, mean 0.19; RND1 −1.69, mean 0.18) and three myeloid genes (CD1C −1.90, mean 0.10; CCL3L1 +2.12, mean 2.01; SLC16A10 +1.39, mean 0.68) (Figure S1).

### The myeloid compartment expands while cDC2 numbers are unchanged

Four of the seven attributable myeloid genes are cDC2 markers (CD1C, APOBEC3A, TIFAB, CHN2), which points to loss of the CD1C^+ population. Re-clustering the myeloid compartment gave 11,937 macrophages, 1,578 classical monocytes, 1,158 cDC2, 73 pDC, 10 neutrophils and 89 low-quality cells (Figure 7A–B), with no cDC1 state (no cluster mean above 0.09 for CLEC9A or XCR1) and no separate mast-cell state, mast cells having been partitioned at atlas level.

cDC2 abundance was the same per 1,000 tissue cells (17.7 versus 16.3; d = 0.12, p = 0.84) and in absolute count (65.7 versus 49.7 per donor; d = 0.42, p = 0.49) (Figure 7C; Table 2). The marker-based definition agreed (24.5 versus 24.9; d = −0.03, p = 0.96), as did dropping the two donors with fewest myeloid cells (d = 0.15, p = 0.81) and rescaling every gating threshold by 0.8 or 1.2 (d = 0.59, p = 0.33 and d = 0.12, p = 0.84; Table S7).

What differed was the compartment around them. Myeloid cells were 1.8-fold more abundant in failed fistulas (262 versus 146 per 1,000; d = 1.37, p = 0.044; absolute count d = 1.52, p = 0.026; Figure 7E), so the cDC2 share of that compartment fell from 12.4% to 6.2% (d = −1.06, p = 0.12) (Figure 7D). The population was diluted rather than depleted.

### Within cDC2 the CD1C programme is attenuated, and the compartment signal is dilution

Raw detection of CD1C within cDC2 was lower in failed fistulas (44.7% versus 61.4%), but those cells held 15% fewer UMIs (median 14,474 versus 17,070), which lowers detection of every gene. After downsampling, detection was 31.2% versus 48.1% (d = −1.00, p = 0.13) and counts 1.80 versus 3.23 per 10,000 (d = −1.01, p = 0.13) (Figure 8B; Table S8). CLEC10A moved the same way (d = −1.12) and FCER1A more weakly (d = −0.40), while CD1E was flat (d = −0.16). Depth equalisation was partial: four of seven housekeeping genes showed no difference (all p > 0.5), while RPL13A counts remained higher in failed fistulas (d = 1.04, p = 0.11) and ACTB lower (d = −0.97, p = 0.14), neither significantly. A CD1C-low, C1QA-positive state was more frequent in failed fistulas (0.68 versus 0.10 CLR, d = 0.75, p = 0.23). The direction fits partial attenuation of antigen presentation: all four markers except CD1E point the same way, the magnitude rests on CD1C and CLEC10A, and nothing reaches significance at this sample size.

If the population is intact, the compartment signal has to be dilution. CD1C (observed −1.95) split into a composition term of −1.00 and an intrinsic term of −0.90, predicting −1.90 with a residual of −0.05; APOBEC3A, TIFAB and CHN2 fitted similarly, composition dominated for CCL3L1 (+1.46 of +1.98) and accounted for all of SLC16A10 (+1.46 against +1.25), while YEATS4 fitted badly (residual −1.25) (Figure 8A; Table 3). Across the seven genes the median absolute residual was 0.24 log2 units against observed effects of 0.7–2.7, and composition supplied a median 51% of the change.

### Threshold choices do not drive these conclusions

Two thresholds could have manufactured them: the rule that a compartment-level gene must reach half its maximum across compartments, and the 0.10 expression floor below which no mechanism is claimed. Both were perturbed.

Tightening the origin rule to require 80% of the maximum and a mean above 0.15 left 5 of the 30 genes; relaxing it to 25% of the maximum left 21; the primary rule left 19 (Table S9). The myeloid result survived in the sense that matters: five of the seven attributable myeloid genes are cDC2 markers regardless of which rule is used, and CD1C (mean 0.101, maximum elsewhere 0.066) passes every origin variant, since its signal is higher in myeloid cells than anywhere else. What is not floor-independent is the interpretation of the CD1C fold change: at a floor of 0.15, CD1C and CHN2 drop below the interpretable threshold, leaving CCL3L1 and SLC16A10 in myeloid cells and ACE and RND1 in endothelium (Table S10). The counting result on which the dilution argument rests does not depend on any expression floor, because it uses cell numbers, not transcript levels.

### Fibronectin–integrin is the only axis with nominal evidence of elevation

Of 744 interactions, none reached FDR < 0.05 under exact donor permutation, consistent with a resolution limit of 924 partitions. At pathway level, fibronectin–integrin was higher in failed fistulas (+0.40 log2, d = 1.99, p = 0.014, FDR = 0.13 across nine axes), and stayed higher when restricted to signals received by fibroblasts (+0.38, p = 0.014) (Figure 5A). TGF-β (−0.003, p = 0.98), collagen–integrin/DDR (+0.01, p = 0.96) and PDGF (−0.16, p = 0.47) were unchanged. The largest individual effects were FN1 signals to myeloid and myofibroblast receivers (+0.51 to +0.99), including endothelial→myofibroblast FN1–ITGAV (d = 3.16), and weaker TNF–TNFRSF1A signalling to pericytes (−1.57) (Figure 5B). Some myeloid-directed signal is expected from abundance alone; restricting to fibroblast receivers keeps the effect, and FN1 was up in independent bulk (+0.56, FDR = 0.23). This axis is reported as hypothesis-generating: it survives correction only at pathway level.

### Bulk comparison and a positive control

In GSE220796 (11 failed, 9 matured) 594 genes met threshold. Among pre-specified genes COL8A1 was up (+0.90, FDR = 0.021) and TIMP1 down (−0.75, FDR = 0.048), while FN1 (+0.56), LOXL1 (+0.73), PI16 (−0.80), POSTN (−0.34), ACTA2 (+0.18) and MYH11 (+0.10) were not (Figure 6C). No programme score differed (p = 0.29–0.83); the largest was a lower interferon score in failed fistulas (−0.48, p = 0.29). Gene-level differences without a programme-level shift are what led to the compartment and decomposition analyses. These samples are the paired vein–fistula bulk series reported earlier [[8, 35]].

To check that the instruments are not insensitive, the same scoring was applied to GSE119296 [[37]], comparing pre-access vein with post-anastomosis AVF. There 4,867 genes were differentially expressed and every ECM programme moved as expected: synthetic myofibroblast +1.10, TGF-β +0.67, matrisome +0.61, crosslinking +0.60 (all p ≤ 0.015), with contractile (−0.90) and medial identity (−0.82) falling and interferon (+0.31, p = 0.23) flat. The instruments detect this remodelling and almost nothing when failed and matured fistulas are compared.

### An independent bulk cohort does not reproduce the myeloid difference

The myeloid expansion is the positive compositional finding of this study and rests on 12 donors. We therefore asked whether the same tissue, from an independent cohort, shows it. Deconvolution was not identifiable with a same-tissue reference: non-negative least squares on the 11,219 shared genes collapsed onto a single compartment (median fit R² = 0.31). We used compartment signature scores instead, the mean z-score of each compartment's 50 most specific genes, which test direction rather than proportion.

In GSE220796 the myeloid signature was lower, not higher, in failed fistulas (d = −1.50, p = 0.004; housekeeping-normalised d = −1.18, p = 0.021), as was the endothelial signature (d = −1.72, p = 0.002), while the SMC signature was higher (d = 0.84, p = 0.12) (Figure S3; Table S11). Global sample properties were balanced between groups (mean log2 CPM, detected genes and library size all p > 0.9), so this is not a library-quality artefact. Three differences between cohorts could produce it: GSE220796 samples juxta-anastomotic tissue after anastomosis whereas the single-cell failed samples come from second-stage surgery; signature scores are relative to the cellular transcript pool, so an increase in matrix or medial transcripts lowers every other score; and the two cohorts share no patients. We report the discrepancy rather than resolve it, and we do not claim external validation of the myeloid expansion.

## Limitations of study

The central limitation is resolution. Six fistulas per group exclude standardised effects above roughly 1.2–1.8 depending on the endpoint, and no endpoint reached equivalence at a smallest effect of interest of d = 0.80 (Table S5). Every null reported here should be read as "a large shift is excluded", never as "there is no difference".

The second is cohort structure. Only one public human AVF single-cell dataset exists, so no independent single-cell replication was possible, and the independent bulk cohort disagreed with the single-cell cohort on the direction of the myeloid difference (Figure S3). The specific compositional difference we report is therefore cohort-specific; the decomposition framework that produced it is not.

The third is phenotype. Failure is clinically heterogeneous, spanning early thrombosis, failure to dilate and later stenosis [[28]], and the deposited metadata do not separate these, so "failed" is one category. Time from surgery, dialysis status, comorbidity and access interventions are not recorded, and neither is angiotensin-converting enzyme inhibitor or angiotensin-receptor blocker use, which is directly relevant to the endothelial ACE difference.

The fourth is ambient RNA. The deposited matrices are filtered, so empty-droplet-based soup estimation, and therefore formal correction with tools such as SoupX or decontX [[11]], could not be applied. Contamination was handled by screening genes for cell of origin, which discards affected genes rather than recovering them; the threshold sensitivity in Table S9 shows that the retained gene list is rule-dependent even though the cDC2 conclusion is not.

The fifth is annotation. Myeloid annotation required a gated rule because rank-sum markers were dominated by ribosomal and cytoskeletal genes, and the rule is a judgement; counts were stable under threshold rescaling (Table S7) but a different marker set could move individual cells between states. No cDC1 cells were detected, which may reflect biology or the limits of a 3′ assay on frozen tissue.

The sixth is spatial. Wall-layer identity was inferred from expression, not measured, and no public human AVF spatial dataset exists. A thin neointima cannot be separated from the media and adventitia that surround it, and the bulk signal behind the excess-matrix model came from a single post-anastomosis time point [[8]] whereas our failed samples are established lesions.

Finally, all inferences are observational and correlative. Nothing here establishes that composition changes cause failure rather than accompany it, and the sex imbalance between groups (five female and one male in the failed group, four and two in the matured group) is too small a sample to adjust for.

## Discussion

None of the three pre-specified hypotheses survived. Myofibroblast abundance was not raised in failed fistulas, the fibroblast lineage was not displaced towards a myofibroblast fate (if anything the reverse, d = −0.75), and none of five matrisome programmes moved. What was detectable at compartment level turned out to be mostly compositional. In this cohort a failed fistula differs from a matured one more in how many of which cells it holds than in what those cells are doing.

That is not the reading the atlas was published with. Martinez and colleagues attributed failure to persistent inflammation and an exaggerated healing response, with proinflammatory and adhesive macrophages, activated vasa vasorum endothelium and hyperactive fibroblasts and myofibroblasts, and proposed macrophage-derived osteopontin as the paracrine driver [[12]]. Analysing the same 21 libraries, we see the same cell-level differences and read them differently. The myeloid compartment is larger in failed fistulas, as both accounts show; in our decomposition the expansion is in cell number rather than per-cell transcription. Fibroblasts and myofibroblasts are present in the failed wall; their matrisome output per cell is not higher. The immunofluorescence and in situ hybridisation in that study place specific proteins in specific wall compartments, which these data cannot do, and nothing here contests those images. What is at issue is the step from a cell-rich, protein-rich lesion to the inference that the same cells are transcribing more. The readings differ in the order of the question: whether a difference is assigned to a cell type and then interpreted, or split into a composition term and an intrinsic one first. The atlas reports that the global transcriptional difference between outcomes is subtle, which is where this analysis starts.

The trajectory changes how abundance should be read. Stromal fibroblasts and myofibroblasts sat in the same band of pseudotime (median percentiles 0.48 and 0.56, against 0.18 for resident fibroblasts and 0.83 for smooth muscle), and PAGA resolved one connected module rather than two separable states. The lineage behaves as a gradient with an early and a late pole, so "myofibroblast" names a region of that gradient rather than an attractor a cell enters or does not. Any statistic counting cells on one side of a cut-off, abundance, fate probability or a signature score, depends on where the cut-off falls, which may be one reason such readouts travel poorly between studies. A shift in the distribution without a change in the proportion beyond any threshold would also escape this design.

None of this contradicts the histology. Failing fistulas do contain neointimal hyperplasia and matrix. What it contradicts is one reading of it, that failure means more myofibroblast-driven matrix production than successful maturation. Architecture and cell number are different measurements: a matrix-rich, cell-poor lesion can be built by few long-lived myofibroblasts, and the matrix remains after those cells go quiet. Timing matters too, since our failed samples came from second-stage surgery, when the lesion is established, whereas the bulk signal behind the excess-matrix model came from a single post-anastomosis time point [[8]]; and juxta-anastomotic tissue contains the whole wall, so a thin neointima is diluted by media and adventitia.

The compartment analysis first appeared to locate the biology: 16 genes in myeloid cells, 7 in endothelium, almost none in fibroblasts. Two checks took that apart. Screening for cell of origin showed 9 of 16 myeloid and 2 of 7 endothelial genes to be ambient transcripts from elsewhere. MAL, nominally the most down-regulated myeloid gene at −5.91 log2, is a T-cell transcript found in a compartment with no T cells in it. The scale is not marginal: dissociating vascular tissue releases message from whatever is most abundant, and where T cells are a tenth of all cells the ambient pool carries a great deal of lymphocyte RNA [[10]],[[11]].

The cDC2 result is the one most likely to have been believed. The surviving myeloid genes were dominated by cDC2 markers, which invites the conclusion that the CD1C^+ population collapses: attractive, consistent with the antigen-presentation framing of neointimal lesions [[29]], and wrong. Counting showed the population unchanged in absolute terms under two definitions, after the smallest donors were dropped and after gating thresholds were rescaled. The apparent loss came from a 1.8-fold expansion of the myeloid compartment, overwhelmingly macrophages, which halved the cDC2 share without touching the count. Markers of rare populations are the most fragile readout in any bulk comparison, because apparent abundance is the product of how many such cells exist and how large the compartment around them has become.

That mechanism is internal to one cohort, and the independent bulk cohort did not reproduce the expansion it depends on. Two readings are compatible with both datasets. Either the failure phenotype differs between the cohorts, since one samples tissue at second-stage surgery and the other after anastomosis; or the myeloid expansion is a feature of this particular series of 12 donors. Distinguishing them would need a second single-cell cohort with outcome annotation sampled at the same stage, which the public data do not provide. What does not depend on the point is the arithmetic: whatever the compartment size, marker-based abundance is its product with per-cell expression, and the two terms must be separated before a biological claim is made.

What is left as cell-intrinsic signal is small and insecure. CD1C detection within cDC2 was 17 percentage points lower at matched depth (d = −1.00, p = 0.13), with CLEC10A pointing the same way; we treat it as a hypothesis. On the endothelial side the screen and the floor leave ACE (−1.70) and RND1 (−1.69), two genes in six donors. Lower endothelial ACE is plausible in a shunt, where local angiotensin II handling shapes inward remodelling [[30]], and RND1 is a Rho-family GTPase that promotes actin disassembly and loss of adhesion [[31]]. A confound: ACE inhibitors and angiotensin-receptor blockers are common in CKD stage 5 and unrecorded, so the ACE difference could be treatment. The pericyte difference (d = −1.19) fits detachment and microvascular instability [[32]] but was significant only on rank testing.

Fibronectin–integrin stands apart because it was pre-specified, internally consistent and not attributable to composition: higher at pathway level (+0.40 log2, d = 1.99, p = 0.014), higher when restricted to fibroblast receivers, carried mainly by FN1–ITGAV/ITGB1 from stromal and resident fibroblasts to myeloid cells and myofibroblasts, and matched by FN1 up-regulation in independent bulk (+0.56, FDR = 0.23). Integrin engagement is required for myofibroblast adhesion and for retaining inflammatory cells [[33]], so a rise here without movement in TGF-β, PDGF or collagen–integrin would put the difference in matrix–integrin engagement rather than matrix quantity. Two caveats: FN1 is an acute-phase reactant and rises with the systemic inflammation common in ESKD, and the effect survives correction only at pathway level (FDR = 0.13). If it is real it is testable: block integrin engagement rather than TGF-β in an AVF model.

The transferable point is about inference. A substantial literature reports bulk signatures and classifiers for maturation failure [[9]]. Many will be proxies for composition, which does not make them clinically useless, since a classifier can predict outcome from a compositional signature, but it does make their biology uncertain. Reporting the composition of the samples compared, screening candidates for cell of origin, decomposing differences where composition differs, and stating the expression level a claim rests on would prevent most of the errors described here. The same discipline applies to our own positive findings: the myeloid expansion, the strongest compositional signal in this cohort, failed its external check, and we report it as cohort-specific rather than as a property of AVF failure. The point is not hypothetical. A single-nucleus study of vein taken at access creation, rather than at second-stage surgery, reports endothelial interferon and TNF-NFκB programmes and suppressed smooth-muscle bioenergetics in fistulas that later failed [[34]]. Same disease, different window, different reading: the sampling design belongs in the claim.

## Resource availability

### Lead contact

Requests for further information and resources should be directed to and will be fulfilled by the lead contact, Daijun Hu (283624940@qq.com).

### Materials availability

This study did not generate new materials.

### Data and code availability

- Single-cell and bulk data are publicly available at NCBI Gene Expression Omnibus under accession numbers GSE305423 [[35]], GSE220796 [[36]] and GSE119296 [[37]]. Accession numbers are listed in the key resources table. No new primary data were generated.
- All analysis code, parameter files, random seeds and environment specifications are available at https://github.com/hu9872/avf-composition-intrinsic.git and have been archived at Zenodo under the concept DOI https://doi.org/10.5281/zenodo.23887195 [[38]], which always resolves to the most recent released version; the tagged release analysed here is listed in the key resources table.
- Any additional information required to reanalyse the data reported in this paper is available from the lead contact on request.

## Acknowledgments

We thank the investigators who generated and deposited the publicly available datasets reanalysed in this study (GSE305423, GSE220796 and GSE119296) for making their data and metadata openly available. This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

## Author contributions

Yuewu Tang: conceptualisation, methodology, software, formal analysis, investigation, data curation, writing – original draft, visualisation. Daijun Hu: conceptualisation, methodology, validation, supervision, project administration, writing – review and editing. All authors approved the final manuscript and agree to be accountable for all aspects of the work.

## Declaration of interests

The authors declare no competing interests.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the authors used a large language model assistant to draft and refactor analysis code and to edit the manuscript for clarity and language. After using this tool, the authors reviewed and edited the content and verified every reported value against the analysis outputs, and take full responsibility for the content of the published article. No text, image or dataset was generated by the tool, and no result in this paper was produced or interpreted by it without independent verification by the authors.

## Supplemental information

Document S1. Figures S1–S3, Tables S1–S12, and supplemental methods.
Data S1. Per-donor composition, pseudobulk results, decomposition terms and sensitivity analyses, as separate sheets.

## Figure titles and legends

**Figure 1. Single-cell atlas of human AVF tissue.** (A) UMAP of 66,168 cells coloured by compartment. (B) The same embedding coloured by sampling stage. (C) Compartment composition per sampling stage. n = 20 donors; cells are shown for display only and are not treated as replicates.

**Figure 2. Mural and fibroblast states.** (A) UMAP of 30,009 mural and fibroblast cells coloured by state. (B) Mean programme score per state, z-scored across states. (C) Donor-level pseudobulk in the fibroblast lineage (24,157 genes tested); 12 genes reached FDR < 0.05 with an absolute log2 fold change above 0.585 (red). Dashed line, FDR = 0.05; dotted lines, 1.5-fold change. Tests are two-sided Welch t-tests on donor-level values with Benjamini–Hochberg correction; n = 6 failed and 6 matured donors.

**Figure 3. State composition.** (A) Per-donor percentage of each state in matured (green) and failed (red) fistulas; horizontal bars are medians. (B) Difference in centred log-ratio abundance with 95% confidence intervals, Cohen's d and Welch p. No state survived correction across the six states. Symbols: \*p < 0.05, \*\*p < 0.01, \*\*\*p < 0.001; ns, not significant.

**Figure 4. Trajectory and fate bias.** (A) PAGA graph; node area scales with state size, edge width with connectivity. (B) UMAP coloured by pseudotime percentile, rooted at the medoid of the pre-access resident adventitial-like compartment. (C) Difference in absorption probability (failed − matured) for the three terminal fates, with 95% confidence intervals and Cohen's d. (D) Marker dynamics along pseudotime in 20 equal-count bins; lines are LOESS-smoothed means, shaded bands are ±1 SEM across bins.

**Figure 5. Ligand–receptor analysis.** (A) Nine pre-specified axes, donor-level difference with Cohen's d and Welch p. (B) The 12 most significant individual interactions, shown for completeness; none survived correction across 744 tests. Permutation p values are exact over all 924 donor partitions; the smallest attainable value is 0.0011.

**Figure 6. Where the difference lies.** (A) Differentially expressed genes per compartment before and after cell-of-origin screening. (B) Five matrisome programmes and the zonation index; none significant. (C) Pre-specified genes in independent whole-tissue bulk (GSE220796); COL8A1 and TIMP1 reached FDR < 0.05. Bulk tests used PyDESeq2 with outcome as the only factor; n = 11 failed and 9 matured samples.

**Figure 7. The cDC2 compartment is numerically preserved in failed AVF.** (A) UMAP of 14,845 myeloid cells coloured by gated state. (B) Mean marker expression per state; red rows (CLEC9A, XCR1) sit at ambient level. (C) cDC2 per 1,000 tissue cells (d = 0.12, p = 0.84). (D) cDC2 as a percentage of the myeloid compartment (d = −1.06, p = 0.12). (E) Myeloid cells per 1,000 tissue cells (d = 1.37, p = 0.044). Points are donors; horizontal bars are means; p values are two-sided Welch t-tests with n = 6 per group.

**Figure 8. The myeloid pseudobulk signal reflects composition and depth, not cDC2 loss.** (A) Observed compartment log2 fold change against the value predicted from composition plus cell-intrinsic expression (Table 3); the diagonal is perfect agreement. (B) CD1C detection within cDC2 after downsampling every cell to 5,000 UMIs (d = −1.00, p = 0.13). (C) Median UMIs per cDC2 cell before correction. Points are donors; 12 donors (6 failed, 6 matured) contributed cDC2 cells to the downsampled comparison.

## Tables

**Table 1. Mural and fibroblast state composition, failed versus matured (donor level, n = 6 + 6).**

| State | Failed, mean CLR | Matured, mean CLR | Difference (95% CI) | Cohen's d | Welch p | Mann–Whitney p | FDR |
|---|---|---|---|---|---|---|---|
| Myofibroblast | 0.074 | −0.047 | 0.121 (−0.983 to 1.225) | 0.12 | 0.836 | 0.818 | 0.836 |
| AF resident | 1.748 | 1.535 | 0.213 (−1.054 to 1.479) | 0.19 | 0.749 | 0.818 | 0.836 |
| AF stromal | 0.795 | 1.221 | −0.426 (−1.027 to 0.174) | −0.80 | 0.195 | 0.589 | 0.586 |
| SMC contractile | 0.738 | 0.396 | 0.342 (−0.696 to 1.379) | 0.37 | 0.541 | 0.937 | 0.811 |
| SMC synthetic | −3.228 | −3.923 | 0.694 (−1.026 to 2.415) | 0.46 | 0.455 | 0.485 | 0.811 |
| Pericyte | −0.126 | 0.817 | −0.943 (−1.844 to −0.043) | −1.19 | 0.069 | 0.041 | 0.411 |

CLR, centred log-ratio; AF, adventitial fibroblast; SMC, smooth muscle cell; CI, confidence interval; FDR, false discovery rate across the six states (Benjamini–Hochberg). Proportions were transformed after a 0.5 pseudo-count. Tests are two-sided.

**Table 2. Myeloid and cDC2 abundance, failed versus matured (donor level, n = 6 + 6).**

| Metric | Failed, mean | Matured, mean | Difference (95% CI) | Cohen's d | Welch p | FDR |
|---|---|---|---|---|---|---|
| Myeloid cells per 1,000 tissue cells | 262.4 | 145.5 | 116.9 (20.7 to 213.1) | 1.38 | 0.044 | 0.175 |
| Myeloid cells per donor | 1015.3 | 525.3 | 490.0 (124.4 to 855.6) | 1.52 | 0.026 | 0.175 |
| cDC2 per 1,000 tissue cells | 17.7 | 16.3 | 1.3 (−11.2 to 13.8) | 0.12 | 0.842 | 0.962 |
| cDC2 cells per donor | 65.7 | 49.7 | 16.0 (−27.3 to 59.3) | 0.42 | 0.486 | 0.777 |
| cDC2, marker rule, per 1,000 | 24.5 | 24.9 | −0.4 (−16.3 to 15.5) | −0.03 | 0.962 | 0.962 |
| cDC2 as % of myeloid compartment | 6.2 | 12.4 | −6.2 (−12.9 to 0.4) | −1.06 | 0.120 | 0.283 |

cDC2, conventional type 2 dendritic cell; CI, confidence interval. Two cDC2 definitions are shown: a gated rule on cluster-mean expression (rows 3–6) and a cluster-independent marker rule requiring two of CD1C, CLEC10A, FCER1A and CD1E (row 5). FDR is Benjamini–Hochberg across the eight metrics. Tests are two-sided.

**Table 3. Decomposition of myeloid pseudobulk differences into composition and cell-intrinsic terms (log2 units).**

| Gene | Dominant state | Observed log2FC | Composition term | Intrinsic term | Predicted | Residual | % from composition |
|---|---|---|---|---|---|---|---|
| CD1C | cDC2 | −1.95 | −1.00 | −0.90 | −1.90 | −0.05 | 51 |
| APOBEC3A | cDC2 | −2.68 | −1.00 | −1.61 | −2.61 | −0.07 | 37 |
| TIFAB | cDC2 | −2.43 | −1.00 | −1.65 | −2.65 | 0.22 | 41 |
| CHN2 | cDC2 | −1.61 | −1.00 | −1.24 | −2.24 | 0.63 | 62 |
| CCL3L1 | Classical monocyte | 1.98 | 1.46 | 0.88 | 2.33 | −0.35 | 74 |
| SLC16A10 | Classical monocyte | 1.25 | 1.46 | 0.03 | 1.49 | −0.24 | 117 |
| YEATS4 | pDC | 0.72 | −0.87 | 2.84 | 1.97 | −1.25 | −120 |

log2FC, log2 fold change. Predicted = composition + intrinsic; residual = observed − predicted. All seven genes passed cell-of-origin screening; CD1C, CCL3L1 and SLC16A10 exceed the 0.10 log-normalised expression floor. Median |residual| = 0.24; median composition contribution = 51% of the observed effect. Percentages outside 0–100 indicate genes whose observed change is smaller than, or opposite to, the composition term.

## STAR Methods

### Resource availability

#### Lead contact

Further information and requests should be directed to and will be fulfilled by the lead contact, Daijun Hu (283624940@qq.com).

#### Materials availability

This study did not generate new materials.

#### Data and code availability

Public datasets are listed in the key resources table and in Resource availability above. Analysis code is deposited at https://github.com/hu9872/avf-composition-intrinsic.git and archived at Zenodo under the concept DOI https://doi.org/10.5281/zenodo.23887195, which always resolves to the most recent released version. No accession numbers are restricted.

### Experimental model and study participant details

This is a secondary analysis of publicly deposited, de-identified human data; no new human or animal subjects were involved and no ethical approval was required. GSE305423 contains 20 patients with chronic kidney disease stage 5 or end-stage kidney disease undergoing access surgery [[12]]; six pre-access cephalic veins, three early AVF resections from two patients, and six matured plus six failed juxta-anastomotic AVFs. Age at sampling ranged from 28 to 72 years (Table S12). Sex was recorded for every donor: the failed group comprised five female and one male donor and the matured group four female and two male donors; the pre-access group comprised four female and two male donors. With six donors per group the design has no power to test sex effects, and none was attempted; the imbalance is reported in Limitations of study rather than adjusted for. GSE220796 and GSE119296 are bulk cohorts of vein and AVF tissue; donor-level covariates beyond outcome are not deposited, and sex is therefore not modelled in the bulk analyses.

### Method details

#### Data sources

Three public datasets were retrieved from GEO and analysed locally. GSE305423 (10x Genomics 3′ scRNA-seq) holds 21 libraries from 20 patients [[12]]; the two early libraries from one patient were merged into one donor and kept out of outcome comparisons. GSE220796 provides paired bulk RNA-seq of basilic vein and juxta-anastomotic AVF [[8]]; 20 of its 41 samples formed complete pairs (9 matured, 11 failed). GSE119296 (19 vein and 19 AVF bulk libraries, one shared identifier) has no usable pairing and served only as a positive control. No public human AVF spatial dataset exists, and no spatial claim is made.

#### Preprocessing and annotation

Cells passed a gene, UMI, mitochondrial and complexity filter (74,285 of 97,070); 497 doublets were removed with Scrublet [[14]],[[15]]. Counts were log-normalised and integrated with Harmony by donor [[16]], and Leiden clusters were assigned to the compartment whose canonical panel supplied most of their top 50 rank-sum markers, normalised by panel length. Clusters with fewer than three supporting markers, or with mean haemoglobin above 1.0, were dropped, leaving 66,168 cells in eight compartments (Table S1).

#### Mural and fibroblast subclustering

The 30,009 mural and fibroblast cells were re-clustered at Leiden resolution 1.0, the resolution that separated mural states without fragmenting fibroblasts. Six states were named from their highest-scoring canonical panel; panels, marker support and per-state scores are in Table S2.

#### Myeloid subclustering and gated annotation

The 14,845 myeloid cells were re-clustered after immunoglobulin and haemoglobin genes were removed (Leiden 0.9). A marker-support rule left 76% of cells unassigned and found no cDC2, rank-sum markers being dominated by ribosomal and cytoskeletal genes. Clusters were therefore annotated by a gated rule on cluster-mean expression, with thresholds applied in priority order and deliberately permissive against published marker levels [[21]], to avoid assigning a cluster to a lineage whose defining gene is absent (Table S2). No cDC1 state emerged. A cluster-independent definition (at least two of CD1C, CLEC10A, FCER1A, CD1E) and rescaling of every threshold by 0.8 or 1.2 served as sensitivity analyses (Table S7).

#### Cell-of-origin screening and expression floor

Each compartment-level gene was screened for cell of origin, being retained if its mean in the target compartment reached half the maximum across compartments, or 0.15 log-normalised units. Origin and quantifiability are different questions, so a uniform 0.10 floor was applied afterwards: genes above it were interpreted, those below were reported and carried through the arithmetic without mechanistic weight. The origin rule was perturbed to 25%, 50% and 80% of the maximum, alone and combined with a 0.15 floor, and the floor was varied at 0.05, 0.10 and 0.15 (Tables S9 and S10).

#### Depth control

Within cDC2, every cell with at least 5,000 UMIs was downsampled to exactly 5,000 by multinomial resampling, leaving 1,086 of 1,158 cells from 20 donors and 12 post-operative donors for the outcome test; seven housekeeping genes served as positive controls (Table S8).

#### Composition–intrinsic decomposition

Each log2 fold change was decomposed as log2FC(observed) ≈ log2(f_failed/f_matured) + log2(rate_failed/rate_matured), where f is the share of the compartment held by the gene's dominant state and rate its expression within that state. The residual, observed minus predicted, measures agreement. The decomposition assumes that one state dominates the gene's signal within the compartment; genes for which that assumption fails carry a large residual and are flagged (Table 3).

#### Trajectory and fate bias

PAGA described state connectivity [[22]]. Diffusion pseudotime was rooted at the medoid of the pre-access resident adventitial-like cells [[23]]; rooting at the extreme of the first diffusion component, whose position is set by outliers, placed the origin at the smooth-muscle end and was not used. Because a single outlying cell stretches the absolute scale to 1.0 while more than 99% of cells lie below 0.1, pseudotime is reported as a percentile. Fate bias was a multi-fate absorption probability on a Markov chain oriented along increasing pseudotime, the directed formulation of CellRank [[24]]; unoriented, absorption into any single terminal set is trivially 1. Terminal fates were the upper pseudotime quartile of each state, with 1 × 10⁻³ leakage to keep the system non-singular. RNA velocity was not estimated, as only filtered matrices are deposited.

#### Ligand–receptor analysis

A pre-specified panel of 34 ligand–receptor pairs (CellPhoneDB v2) was scored as the product of sender ligand and receiver receptor means [[25]]; donor-averaged scores were tested by exact permutation over all 924 partitions and aggregated into nine pre-specified axes, the smallest attainable p (about 0.0011) precluding FDR control across 744 interactions.

#### Matrisome programmes

Matrisome genes came from GO:0062023 (372 genes), cross-checked against a curated core matrisome [[26]], and decomposed by non-negative matrix factorisation [[27]]; error fell by less than 2% per factor beyond k = 5, and k = 5 recovered the expected categories. Wall-layer identity was scored for intima, media and adventitia, and a zonation index defined as media minus adventitia. These samples have no spatial coordinates, so zoning is reconstructed from expression and carries no spatial claim.

#### Bulk analyses

Bulk cohorts were analysed with PyDESeq2 (failed versus matured, 11 versus 9) [[20]], programme scores being the mean of gene-wise z-scores. Genes were retained at FDR < 0.05 and an absolute log2 fold change above 0.585.

#### Compartment signature scores

Because deconvolution with a same-tissue reference was not identifiable (non-negative least squares on 11,219 shared genes collapsed onto one compartment, median fit R² = 0.31), compartment abundance was probed with signature scores: for each compartment, the 50 genes with the highest specificity in the study's own reference profile (log2 CPM in the compartment minus the mean over other compartments, restricted to genes above 10 CPM), scored as the mean gene-wise z-score across bulk samples, with a housekeeping-normalised variant using the mean of RPL13A, RPLP0, B2M, TBP, GAPDH, ACTB and PPIA. Scores were compared between outcomes by Welch's t-test (Table S11).

### Quantification and statistical analysis

The donor was the unit of replication throughout; no test was run at cell level [[17]]. Proportions were centred-log-ratio transformed after a 0.5 pseudo-count [[18]] and compared between groups by Welch's t-test and Mann–Whitney U, with Cohen's d, 95% confidence intervals and Benjamini–Hochberg correction within each family of tests. Pseudobulk counts were summed per donor and state (at least 10 cells) and tested with PyDESeq2, design ~outcome, on post-operative libraries only [[19]],[[20]]. Tests are two-sided and α = 0.05 unless stated. No correction was applied across families; the fibronectin–integrin result is therefore reported with its within-family FDR (0.13) and treated as hypothesis-generating.

Because six donors per group detect only large effects, non-significant results are reported with effect size and interval and read as absence of a large effect, not as equivalence [[13]]. Power was computed for a two-sample t-test at α = 0.05 by non-central t distribution; 80% power requires |d| = 1.80 at n = 6 per group and |d| = 0.79 at n = 26 (Table S5). Equivalence was tested by two one-sided tests against a smallest effect of interest of d = 0.80, with the 90% confidence interval reported alongside.

Software: Python 3.13.14, scanpy 1.12.5, anndata 0.13.4, numpy 2.5.4, pandas 3.1.0, scipy 1.18.1, statsmodels 0.15.0, PyDESeq2 0.5.4, matplotlib 3.11.2. The random seed was 2611 for all stochastic steps. Full parameter values, marker panels and per-step cell counts are given in Document S1.

### Additional resources

The analysis plan, including the three pre-specified hypotheses and the decision rule for the appended cDC2 hypothesis, was written before the corresponding analyses were run. The study was not prospectively registered in a public registry; no registration number exists. The code repository doubles as the record of analysis order.

### Key resources table

| Reagent or resource | Source | Identifier |
|---|---|---|
| Deposited data | | |
| Human AVF single-cell RNA-seq, 21 libraries, 20 patients | Martinez et al. [[12]] | GEO: GSE305423 |
| Human vein and AVF bulk RNA-seq, paired | Martinez et al. [[8]] | GEO: GSE220796 |
| Human vein and AVF bulk RNA-seq | Vazquez-Padron laboratory | GEO: GSE119296 |
| HGNC gene symbols | HGNC | https://www.genenames.org |
| Software and algorithms | | |
| scanpy 1.12.5 | Wolf et al. [[14]] | RRID:SCR_018139 |
| anndata 0.13.4 | scverse | https://anndata.readthedocs.io |
| Scrublet | Wolock et al. [[15]] | RRID:SCR_018098 |
| Harmony | Korsunsky et al. [[16]] | https://github.com/immunogenomics/harmony |
| PAGA | Wolf et al. [[22]] | RRID:SCR_018138 |
| Diffusion pseudotime | Haghverdi et al. [[23]] | RRID:SCR_018138 |
| CellRank | Lange et al. [[24]] | https://cellrank.readthedocs.io |
| CellPhoneDB v2 | Efremova et al. [[25]] | RRID:SCR_017054 |
| PyDESeq2 0.5.4 | Muzellec et al. [[20]] | https://github.com/scverse/pydeseq2 |
| DESeq2 (reference implementation) | Love et al. [[19]] | RRID:SCR_000154 |
| scikit-learn NMF | Pedregosa et al. | RRID:SCR_002577 |
| SciPy 1.18.1 | Virtanen et al. | RRID:SCR_008058 |
| NumPy 2.5.4 | Harris et al. | RRID:SCR_008633 |
| pandas 3.1.0 | McKinney | RRID:SCR_018214 |
| Matplotlib 3.11.2 | Hunter | RRID:SCR_008624 |
| Analysis code, release v1.0.1 (Zenodo concept DOI 10.5281/zenodo.23887195) | This paper | https://github.com/hu9872/avf-composition-intrinsic.git; https://doi.org/10.5281/zenodo.23887195 |

## Resource identification initiative

RRIDs are given in the key resources table for all software with an available identifier. Datasets are identified by GEO accession, which is the persistent identifier used by the repository; no RRID is assigned to data. No antibodies, cell lines or organisms were used.

## References

1. Al-Jaishi, A.A., Oliver, M.J., Thomas, S.M., Lok, C.E., Zhang, J.C., Garg, A.X., et al. (2014). Patency rates of the arteriovenous fistula for hemodialysis: a systematic review and meta-analysis. Am. J. Kidney Dis. 63, 464–478. https://doi.org/10.1053/j.ajkd.2013.08.023
2. Lok, C.E., Huber, T.S., Lee, T., Shenoy, S., Yevzlin, A.S., Abreo, K., et al. (2020). KDOQI clinical practice guideline for vascular access: 2019 update. Am. J. Kidney Dis. 75, S1–S164. https://doi.org/10.1053/j.ajkd.2019.12.001
3. Roy-Chaudhury, P., Sukhatme, V.P., and Cheung, A.K. (2006). Hemodialysis vascular access dysfunction: a cellular and molecular viewpoint. J. Am. Soc. Nephrol. 17, 1112–1127. https://doi.org/10.1681/ASN.2005050615
4. Lee, T., Chauhan, V., Krishnamoorthy, M., Wang, Y., Arend, L., Mistry, M.J., et al. (2011). Severe venous neointimal hyperplasia prior to dialysis access surgery. Nephrol. Dial. Transplant. 26, 2264–2270. https://doi.org/10.1093/ndt/gfq733
5. Stenmark, K.R., Yeager, M.E., El Kasmi, K.C., Nozik-Grayck, E., Gerasimovskaya, E.V., Li, M., et al. (2013). The adventitia: essential regulator of vascular wall structure and function. Annu. Rev. Physiol. 75, 23–47. https://doi.org/10.1146/annurev-physiol-030212-183802
6. Maiellaro, K., and Taylor, W.R. (2007). The role of the adventitia in vascular inflammation. Cardiovasc. Res. 75, 640–648. https://doi.org/10.1016/j.cardiores.2007.06.023
7. Coen, M., Gabbiani, G., and Bochaton-Piallat, M.L. (2011). Myofibroblast-mediated adventitial remodeling: an underestimated player in arterial pathology. Arterioscler. Thromb. Vasc. Biol. 31, 2391–2396. https://doi.org/10.1161/ATVBAHA.111.231548
8. Martinez, L., Rojas, M.G., Tabbara, M., Pereira-Simon, S., Santos Falcon, N., Rauf, M.A., et al. (2023). The transcriptomics of the human vein transformation after arteriovenous fistula anastomosis uncovers layer-specific remodeling and hallmarks of maturation failure. Kidney Int. Rep. 8, 837–850. https://doi.org/10.1016/j.ekir.2023.01.008
9. Li, Y., Chen, Y., Cui, W., Wang, J., Chen, X., Zhang, C., et al. (2023). Identification and validation of hub genes in the stenosis of arteriovenous fistula. J. Pers. Med. 13, 207. https://doi.org/10.3390/jpm13020207
10. van den Brink, S.C., Sage, F., Vértesy, Á., Spanjaard, B., Peterson-Maduro, J., Baron, C.S., et al. (2017). Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nat. Methods 14, 935–936. https://doi.org/10.1038/nmeth.4437
11. Young, M.D., and Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9, giaa151. https://doi.org/10.1093/gigascience/giaa151
12. Martinez, L., Stoyell-Conti, F.F., Tabbara, M., Rojas, M.G., Pereira-Simon, S., Santos Falcon, N., et al. (2026). The single-cell landscape of the human vein after arteriovenous fistula creation and implications for maturation failure. Kidney Int. 109, 160–177. https://doi.org/10.1016/j.kint.2025.08.024
13. Altman, D.G., and Bland, J.M. (1995). Absence of evidence is not evidence of absence. BMJ 311, 485. https://doi.org/10.1136/bmj.311.7003.485
14. Wolf, F.A., Angerer, P., and Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 19, 15. https://doi.org/10.1186/s13059-017-1382-0
15. Wolock, S.L., Lopez, R., and Klein, A.M. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. Cell Syst. 8, 281–291.e9. https://doi.org/10.1016/j.cels.2018.11.005
16. Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat. Methods 16, 1289–1296. https://doi.org/10.1038/s41592-019-0619-0
17. Squair, J.W., Gautier, M., Kathe, C., Anderson, M.A., James, N.D., Hutson, T.H., et al. (2021). Confronting false discoveries in single-cell differential expression. Nat. Commun. 12, 5692. https://doi.org/10.1038/s41467-021-25960-2
18. Gloor, G.B., Macklaim, J.M., Pawlowsky-Glahn, V., and Egozcue, J.J. (2017). Microbiome datasets are compositional: and this is not optional. Front. Microbiol. 8, 2224. https://doi.org/10.3389/fmicb.2017.02224
19. Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550. https://doi.org/10.1186/s13059-014-0550-8
20. Muzellec, B., Teleńczuk, M., Cabeli, V., and Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics 39, btad547. https://doi.org/10.1093/bioinformatics/btad547
21. Villani, A.C., Satija, R., Reynolds, G., Sarkizova, S., Shekhar, K., Fletcher, J., et al. (2017). Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors. Science 356, eaah4573. https://doi.org/10.1126/science.aah4573
22. Wolf, F.A., Hamey, F.K., Plass, M., Solana, J., Dahlin, J.S., Göttgens, B., et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biol. 20, 59. https://doi.org/10.1186/s13059-019-1663-x
23. Haghverdi, L., Büttner, M., Wolf, F.A., Buettner, F., and Theis, F.J. (2016). Diffusion pseudotime robustly reconstructs lineage branching. Nat. Methods 13, 845–848. https://doi.org/10.1038/nmeth.3971
24. Lange, M., Bergen, V., Klein, M., Setty, M., Reuter, B., Bakhti, M., et al. (2022). CellRank for directed single-cell fate mapping. Nat. Methods 19, 159–170. https://doi.org/10.1038/s41592-021-01346-6
25. Efremova, M., Vento-Tormo, M., Teichmann, S.A., and Vento-Tormo, R. (2020). CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. Nat. Protoc. 15, 1484–1506. https://doi.org/10.1038/s41596-020-0292-x
26. Naba, A., Clauser, K.R., Hoersch, S., Liu, H., Carr, S.A., and Hynes, R.O. (2012). The matrisome: in silico definition and in vivo characterization by proteomics of normal and tumor extracellular matrices. Mol. Cell. Proteomics 11, M111.014647. https://doi.org/10.1074/mcp.M111.014647
27. Lee, D.D., and Seung, H.S. (1999). Learning the parts of objects by non-negative matrix factorization. Nature 401, 788–791. https://doi.org/10.1038/44565
28. Beathard, G.A., Arnold, P., Jackson, J., and Litchfield, T. (2003). Aggressive treatment of early fistula failure. Kidney Int. 64, 1487–1494. https://doi.org/10.1046/j.1523-1755.2003.00210.x
29. Zernecke, A. (2015). Dendritic cells in atherosclerosis: evidence in mice and humans. Arterioscler. Thromb. Vasc. Biol. 35, 763–770. https://doi.org/10.1161/ATVBAHA.114.303566
30. Fleming, I. (2006). Signaling by the angiotensin-converting enzyme. Circ. Res. 98, 887–896. https://doi.org/10.1161/01.RES.0000217340.40936.53
31. Nobes, C.D., Lauritzen, I., Mattei, M.G., Paris, S., Hall, A., and Chardin, P. (1998). A new member of the Rho family, Rnd1, promotes disassembly of actin filament structures and loss of cell adhesion. J. Cell Biol. 141, 187–197. https://doi.org/10.1083/jcb.141.1.187
32. Armulik, A., Genové, G., and Betsholtz, C. (2011). Pericytes: developmental, physiological, and pathological perspectives, problems, and promises. Dev. Cell 21, 193–215. https://doi.org/10.1016/j.devcel.2011.07.001
33. Humphries, J.D., Chastney, M.R., Askari, J.A., and Humphries, M.J. (2019). Signal transduction via integrin adhesion complexes. Curr. Opin. Cell Biol. 56, 14–21. https://doi.org/10.1016/j.ceb.2018.08.004
34. Shah, N.A., Thomas, R., Thomas, S.D., Au, A.Y.M., Davidson, T., O'Keeffe, C.J., et al. (2026). Seq-ing the truth: baseline snRNA-seq links endothelial interferon/TNF-NFκB signalling and smooth muscle bioenergetic suppression to arteriovenous fistula maturation. Kidney360. Published online July 17, 2026. https://doi.org/10.34067/KID.0000001283
35. NCBI Gene Expression Omnibus (2025). GSE305423: single-cell RNA sequencing of human arteriovenous fistula tissue. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305423
36. NCBI Gene Expression Omnibus (2023). GSE220796: paired bulk RNA sequencing of human vein and arteriovenous fistula. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220796
37. NCBI Gene Expression Omnibus (2019). GSE119296: bulk RNA sequencing of human vein and arteriovenous fistula. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119296
38. Tang, Y., and Hu, D. (2026). Analysis code for: apparent transcriptional signatures of arteriovenous fistula maturation failure are dominated by cell composition. Zenodo. https://doi.org/10.5281/zenodo.23887195
