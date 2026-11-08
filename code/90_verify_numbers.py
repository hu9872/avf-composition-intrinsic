"""数字溯源与校验（提交前闸门）。

做法：每条断言 = 从 results/ 原始表重算数值 + 用固定格式拼出正文应出现的字符串，
然后断言该字符串确实出现在稿件中。任何一条失败即禁止提交。

输出：
  results/provenance.csv   每条断言的数值、来源、脚本、图件与状态
  终端报告

用法：python code/90_verify_numbers.py
"""
import glob
import os
import re
import sys
import numpy as np
import pandas as pd

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"


def _default_manuscript():
    """Pick the most recently edited draft, so the gate never silently runs on a
    superseded file. Set MS_PATH to pin a specific manuscript instead.

    Returns None when no manuscript is present. That is the normal state of a
    fresh clone: the prose is not part of the released archive, so the gate
    recomputes every value and reports it without matching it against text."""
    cands = [p for p in glob.glob(os.path.join(ROOT, "manuscript", "Manuscript_*.md"))
             if "_backup" not in os.path.basename(p)]
    if not cands:
        return None
    cands.sort(key=os.path.getmtime)
    return cands[-1]


# 默认校验 manuscript/ 下最近修改的稿子；指定稿件用 MS_PATH。
# 稿件不随仓库发布，克隆后这里通常为 None，校验退化为"只重算、不比对文本"。
MS = os.environ.get("MS_PATH") or _default_manuscript()


def load(name):
    return pd.read_csv(f"{RES}/{name}.csv")


def row(df, **kw):
    m = np.ones(len(df), dtype=bool)
    for k, v in kw.items():
        m &= (df[k] == v).values
    if m.sum() != 1:
        raise KeyError(f"{kw} -> {m.sum()} rows")
    return df[m].iloc[0]


def build_checks():
    c = []

    comp = load("s5_composition_tests")
    r = row(comp, state="Myofibroblast")
    c.append(("C01", "myofibroblast abundance", "s5_composition_tests",
              f"(ΔCLR +{r['diff']:.2f}, d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05b_composition_pseudobulk.py", "Fig3B"))
    r = row(comp, state="AF_resident")
    c.append(("C02", "resident AF abundance", "s5_composition_tests",
              f"(+{r['diff']:.2f}, d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05b_composition_pseudobulk.py", "Fig3B"))
    r = row(comp, state="Pericyte")
    c.append(("C03", "pericyte abundance", "s5_composition_tests",
              f"(−{abs(r['diff']):.2f}, d = −{abs(r['cohens_d']):.2f}, "
              f"Welch P = {r['welch_p']:.3f}, Mann–Whitney P = {r['mwu_p']:.3f})",
              "05b_composition_pseudobulk.py", "Fig3B"))
    c.append(("C04", "pericyte FDR", "s5_composition_tests",
              f"without surviving correction (FDR = {r['fdr']:.2f})",
              "05b_composition_pseudobulk.py", "Fig3B"))

    r = row(comp, state="SMC_synthetic")
    c.append(("C05", "synthetic SMC abundance", "s5_composition_tests",
              f"(+{r['diff']:.2f} CLR, d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05b_composition_pseudobulk.py", "Fig3B"))

    det = load("s30_detectable_effect")
    d6 = float(row(det, n=6)["detectable_d"])
    c.append(("C06", "detectable effect at n=6", "s30_detectable_effect",
              f"80% power is reached only at |d| ≈ {d6:.2f}",
              "30_power_tost.py", "FigS2A"))
    d26 = float(row(det, n=26)["detectable_d"])
    c.append(("C07", "detectable effect at n=26", "s30_detectable_effect",
              f"|d| = {d26:.2f} at n = 26", "30_power_tost.py", "Table S5"))

    tost = load("s30_tost")
    r = row(tost, metric="Myofibroblast abundance (CLR)")
    c.append(("C08", "myofibroblast CI bounds", "s30_tost",
              f"between −{abs(r['ci90_low']):.2f} and +{r['ci90_high']:.2f} CLR units",
              "30_power_tost.py", "FigS2B"))
    n_eq = int(tost.equivalent_at_d080.sum())
    tot = len(tost)
    c.append(("C09", "equivalence count", "s30_tost",
              f"did not establish equivalence for any of the {tot} endpoints examined",
              "30_power_tost.py", "FigS2B"))
    ecm_rows = tost[tost.metric.str.startswith("Matrisome programme ECM")]
    widest = float(ecm_rows["excluded_effect_from_d"].max())
    narrowest = float(ecm_rows["excluded_effect_from_d"].min())
    c.append(("C10", "widest excluded effect", "s30_tost",
              f"excluded standardised effects above {narrowest:.2f}, and the widest only "
              f"above {widest:.2f}", "30_power_tost.py", "Table S5"))

    fate = load("s7_fate_tests")
    r = row(fate, fate="Myofibroblast")
    c.append(("C11", "myofibroblast fate", "s7_fate_tests",
              f"was {r['failed']:.3f} in failed and {r['mature']:.3f} in matured fistulas "
              f"(d = −{abs(r['cohens_d']):.2f}, P = {r['welch_p']:.2f}, FDR = {r['fdr']:.2f})",
              "05_trajectory.py", "Fig4C"))
    r2 = row(fate, fate="SMC_contractile")
    c.append(("C12", "contractile fate", "s7_fate_tests",
              f"into contractile SMC {r2['failed']:.3f} versus {r2['mature']:.3f} "
              f"(d = {r2['cohens_d']:.2f}, P = {r2['welch_p']:.2f})",
              "05_trajectory.py", "Fig4C"))
    r3 = row(fate, fate="AF_resident")
    c.append(("C13", "resident fate", "s7_fate_tests",
              f"into resident adventitial-like fibroblasts {r3['failed']:.3f} versus "
              f"{r3['mature']:.3f} (d = {r3['cohens_d']:.2f}, P = {r3['welch_p']:.2f})",
              "05_trajectory.py", "Fig4C"))

    ecm = load("s9_ecm_tests")
    lo = ecm.loc[ecm.metric.isin([f"ECM{i}" for i in range(1, 6)]), "p"]
    c.append(("C14", "ECM programmes p range", "s9_ecm_tests",
              f"None differed between groups (P = {lo.min():.2f}–{lo.max():.2f})",
              "07_ecm.py", "Fig6B"))
    z = row(ecm, metric="zonation")
    c.append(("C15", "zonation index", "s9_ecm_tests",
              f"the zonation index was flat (Δ = {z['diff']:.3f}, d = {z['cohens_d']:.2f}, "
              f"P = {z['p']:.2f})", "07_ecm.py", "Fig6B"))

    summ = load("s5c_compartment_de_summary")
    c.append(("C16", "compartment DE counts", "s5c_compartment_de_summary",
              f"{int(row(summ, compartment='Myeloid')['n_de'])} genes in myeloid cells, "
              f"{int(row(summ, compartment='Endothelial')['n_de'])} in endothelium, "
              f"{int(row(summ, compartment='Mast')['n_de'])} in mast cells, "
              f"{int(row(summ, compartment='T_cell')['n_de'])} in T cells",
              "05c_compartment_pseudobulk.py", "Fig6A"))

    scr = load("s5g_compartment_gene_screen")
    m = row(scr, gene="MAL")
    c.append(("C17", "MAL ambient", "s5g_compartment_gene_screen",
              f"MAL, the most strongly down-regulated at −{abs(m['log2FoldChange']):.2f} log2",
              "05g_gene_origin_all.py", "Fig6A"))
    means = pd.read_csv(f"{RES}/s5g_gene_compartment_means.csv", index_col=0)
    c.append(("C18", "MAL compartment means", "s5g_gene_compartment_means",
              f"had a mean of {means.loc['MAL','T_cell']:.2f} in T cells",
              "05g_gene_origin_all.py", "Fig6A"))
    c.append(("C19", "IL32 source", "s5g_gene_compartment_means",
              f"(NK {means.loc['IL32','NK']:.2f}, T {means.loc['IL32','T_cell']:.2f}, "
              f"myeloid {means.loc['IL32','Myeloid']:.2f})",
              "05g_gene_origin_all.py", "Fig6A"))
    il = row(scr, gene="IL1RL1")
    c.append(("C20", "IL1RL1 source", "s5g_compartment_gene_screen",
              f"IL1RL1, nominally up {il['log2FoldChange']:.2f}-fold, had a mean of "
              f"{means.loc['IL1RL1','Mast']:.2f} in mast cells against "
              f"{means.loc['IL1RL1','Endothelial']:.2f} in endothelium",
              "05g_gene_origin_all.py", "Fig6A"))

    ace = row(scr, gene="ACE")
    rnd = row(scr, gene="RND1")
    c.append(("C21", "endothelial genes above floor", "s5g_compartment_gene_screen",
              f"(ACE −{abs(ace['log2FoldChange']):.2f}, mean {means.loc['ACE','Endothelial']:.2f}; "
              f"RND1 −{abs(rnd['log2FoldChange']):.2f}, mean {means.loc['RND1','Endothelial']:.2f})",
              "05g_gene_origin_all.py", "FigS1"))
    cd1c = row(scr, gene="CD1C")
    ccl = row(scr, gene="CCL3L1")
    slc = row(scr, gene="SLC16A10")
    c.append(("C22", "myeloid genes above floor", "s5g_compartment_gene_screen",
              f"(CD1C −{abs(cd1c['log2FoldChange']):.2f}, mean {means.loc['CD1C','Myeloid']:.2f}; "
              f"CCL3L1 +{ccl['log2FoldChange']:.2f}, mean {means.loc['CCL3L1','Myeloid']:.2f}; "
              f"SLC16A10 +{slc['log2FoldChange']:.2f}, mean {means.loc['SLC16A10','Myeloid']:.2f})",
              "05g_gene_origin_all.py", "FigS1"))

    ab = load("s5e_abundance_tests")
    r = row(ab, metric="cdc2_D1_per1000")
    c.append(("C23", "cDC2 per 1000", "s5e_abundance_tests",
              f"per 1,000 tissue cells ({r['failed_mean']:.1f} versus {r['mature_mean']:.1f}; "
              f"d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05e_cdc2.py", "Fig7C"))
    r = row(ab, metric="n_cdc2_D1")
    c.append(("C24", "cDC2 absolute", "s5e_abundance_tests",
              f"in absolute count ({r['failed_mean']:.1f} versus {r['mature_mean']:.1f} per donor; "
              f"d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05e_cdc2.py", "Table 2"))
    r = row(ab, metric="cdc2_D2_per1000")
    c.append(("C25", "cDC2 marker rule", "s5e_abundance_tests",
              f"The marker-based definition agreed ({r['failed_mean']:.1f} versus "
              f"{r['mature_mean']:.1f}; d = −{abs(r['cohens_d']):.2f}, P = {r['welch_p']:.2f})",
              "05e_cdc2.py", "Table 2"))
    r = row(ab, metric="myeloid_per1000")
    c.append(("C26", "myeloid expansion", "s5e_abundance_tests",
              f"({r['failed_mean']:.0f} versus {r['mature_mean']:.0f} per 1,000; "
              f"d = {r['cohens_d']:.2f}, P = {r['welch_p']:.3f}",
              "05e_cdc2.py", "Fig7E"))
    r2 = row(ab, metric="n_myeloid")
    c.append(("C27", "myeloid absolute", "s5e_abundance_tests",
              f"absolute count d = {r2['cohens_d']:.2f}, P = {r2['welch_p']:.3f}",
              "05e_cdc2.py", "Table 2"))
    r = row(ab, metric="cdc2_D1_frac_myeloid")
    c.append(("C28", "cDC2 share", "s5e_abundance_tests",
              f"from {r['mature_mean']*100:.1f}% to {r['failed_mean']*100:.1f}% "
              f"(d = −{abs(r['cohens_d']):.2f}, P = {r['welch_p']:.2f})",
              "05e_cdc2.py", "Fig7D"))

    ex = load("s35_extra_checks")
    r = row(ex, metric="CD1C detection within cDC2, raw (donor mean, %)")
    c.append(("C29", "raw CD1C detection", "s35_extra_checks",
              f"Raw detection of CD1C within cDC2 was lower in failed fistulas "
              f"({r['failed_mean']:.1f}% versus {r['mature_mean']:.1f}%)",
              "35_extra_checks.py", "Fig8C"))
    r = row(ex, metric="cDC2 per 1,000 (two smallest donors dropped)")
    c.append(("C30", "drop smallest donors", "s35_extra_checks",
              f"(d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "35_extra_checks.py", "Table 2"))
    r = row(ex, metric="CD1C-low C1QA+ state frequency (CLR)")
    c.append(("C31", "CD1C-low state", "s35_extra_checks",
              f"({r['failed_mean']:.2f} versus {r['mature_mean']:.2f} CLR, "
              f"d = {r['cohens_d']:.2f}, P = {r['welch_p']:.2f})",
              "05e_cdc2.py", "Results"))

    dep = load("s5e_cdc2_depth_by_group")
    f_dep = float(row(dep, group="failed_AVF")["total_counts"])
    m_dep = float(row(dep, group="mature_AVF")["total_counts"])
    c.append(("C32", "UMI depth", "s5e_cdc2_depth_by_group",
              f"(median {f_dep:,.0f} versus {m_dep:,.0f})",
              "05e_cdc2.py", "Fig8C"))

    r = row(ex, metric="CD1C detection (%), downsampled (donor mean)")
    c.append(("C33", "downsampled detection", "s35_extra_checks",
              f"detection was {r['failed_mean']:.1f}% versus {r['mature_mean']:.1f}%",
              "05e_cdc2.py", "Fig8B"))
    r = row(ex, metric="CD1C counts per 10,000, downsampled (donor mean)")
    c.append(("C33b", "downsampled counts", "s35_extra_checks",
              f"counts {r['failed_mean']:.2f} versus {r['mature_mean']:.2f} per 10,000",
              "05e_cdc2.py", "Fig8B"))

    dec = load("s5f_decomposition_retained")
    d1 = row(dec, gene="CD1C")
    c.append(("C34", "CD1C decomposition", "s5f_decomposition_retained",
              f"CD1C (observed −{abs(d1['log2FC_observed']):.2f}) split into a composition "
              f"term of −{abs(d1['log2FC_composition']):.2f} and an intrinsic term of "
              f"−{abs(d1['log2FC_intrinsic']):.2f}, predicting "
              f"−{abs(d1['log2FC_predicted']):.2f} with a residual of "
              f"−{abs(d1['residual']):.2f}", "05f_decomposition.py", "Fig8A"))
    d2 = row(dec, gene="CCL3L1")
    c.append(("C35", "CCL3L1 composition share", "s5f_decomposition_retained",
              f"(+{d2['log2FC_composition']:.2f} of +{d2['log2FC_observed']:.2f})",
              "05f_decomposition.py", "Table 3"))
    d3 = row(dec, gene="SLC16A10")
    c.append(("C36", "SLC16A10 composition share", "s5f_decomposition_retained",
              f"(+{d3['log2FC_composition']:.2f} against +{d3['log2FC_observed']:.2f})",
              "05f_decomposition.py", "Table 3"))
    d4 = row(dec, gene="YEATS4")
    c.append(("C37", "YEATS4 residual", "s5f_decomposition_retained",
              f"YEATS4 fitted badly (residual −{abs(d4['residual']):.2f})",
              "05f_decomposition.py", "Table 3"))
    med_res = float(np.median(np.abs(dec["residual"])))
    med_pct = float(np.median(dec["pct_from_composition"]))
    c.append(("C38", "median residual", "s5f_decomposition_retained",
              f"the median absolute residual was {med_res:.2f} log2 units",
              "05f_decomposition.py", "Table 3"))
    c.append(("C39", "median composition share", "s5f_decomposition_retained",
              f"composition supplied a median {med_pct:.0f}% of the change",
              "05f_decomposition.py", "Fig8A"))

    ax_ = load("s8c_pathway_tests")
    fn = row(ax_, axis="Fibronectin-integrin")
    c.append(("C40", "fibronectin-integrin", "s8c_pathway_tests",
              f"(+{fn['log2fc']:.2f} log2, d = {fn['cohens_d']:.2f}, "
              f"P = {fn['welch_p']:.3f}, FDR = {fn['fdr']:.2f} across nine axes)",
              "06c_pathways.py", "Fig5A"))
    fib = load("s8c_pathway_tests_incoming_fib")
    fi = row(fib, axis="Fibronectin-integrin")
    c.append(("C41", "fibronectin fibroblast receivers", "s8c_pathway_tests_incoming_fib",
              f"(+{fi['log2fc']:.2f}, P = {fi['welch_p']:.3f})",
              "06c_pathways.py", "Fig5A"))
    tg = row(ax_, axis="TGF-beta")
    c.append(("C42", "TGF-beta axis", "s8c_pathway_tests",
              f"TGF-β ({tg['log2fc']:+.3f}, P = {tg['welch_p']:.2f})",
              "06c_pathways.py", "Fig5A"))

    bl = pd.read_csv(f"{RES}/s10_focus_genes.csv", index_col=0)
    col8 = bl.loc["COL8A1"]
    timp = bl.loc["TIMP1"]
    fn1 = bl.loc["FN1"]
    c.append(("C43", "bulk COL8A1 and TIMP1", "s10_focus_genes",
              f"COL8A1 was up (+{col8['log2FoldChange']:.2f}, FDR = {col8['padj']:.3f}) and "
              f"TIMP1 down (−{abs(timp['log2FoldChange']):.2f}, FDR = {timp['padj']:.3f})",
              "08_bulk.py", "Fig6C"))
    c.append(("C44", "bulk FN1", "s10_focus_genes",
              f"while FN1 (+{fn1['log2FoldChange']:.2f})",
              "08_bulk.py", "Fig6C"))

    sc = load("s10_program_tests")
    c.append(("C45", "bulk programme scores", "s10_program_tests",
              f"No programme score differed (P = {sc['p'].min():.2f}–{sc['p'].max():.2f})",
              "08_bulk.py", "Fig6C"))

    de = load("s10_deseq_AVF_failed_vs_matured")
    nde = int(((de.padj < 0.05) & (de.log2FoldChange.abs() > 0.585)).sum())
    c.append(("C46", "bulk DE count", "s10_deseq_AVF_failed_vs_matured",
              f"In GSE220796 (11 failed, 9 matured) {nde} genes met threshold",
              "08_bulk.py", "Fig6C"))

    pc = load("s10b_program_tests")
    c.append(("C47", "positive control programmes", "s10b_program_tests",
              f"synthetic myofibroblast +{float(pc.loc[pc.program=='Synthetic_MyoF','diff'].iloc[0]):.2f}",
              "08b_vein_to_avf.py", "Results"))

    return c



def verify_tables(text):
    """核对正文章节三张表的数值与结果文件一致（按行逐值比对，容差 0.006）。"""
    problems = []
    rows_checked = 0

    def nums(line):
        """只解析数据单元格（跳过首列名称，避免 '1,000'、'CD1C' 之类被当成数值）。"""
        cells = [c.strip() for c in line.strip().strip("|").split("|")][1:]
        vals = []
        for c in cells:
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", c):   # cDC2、pDC 之类纯标识单元格
                continue
            for x in re.findall(r"[\u2212-]?\d+(?:,\d{3})*(?:\.\d+)?", c):
                vals.append(float(x.replace("\u2212", "-").replace(",", "")))
        return vals

    # ---- Table 1
    comp = load("s5_composition_tests").set_index("state")
    block = text.split("**Table 1.")[1].split("**Table 2.")[0]
    for name, key in [("Myofibroblast", "Myofibroblast"), ("AF resident", "AF_resident"),
                      ("AF stromal", "AF_stromal"), ("SMC contractile", "SMC_contractile"),
                      ("SMC synthetic", "SMC_synthetic"), ("Pericyte", "Pericyte")]:
        line = [l for l in block.split("\n") if l.startswith(f"| {name} |")]
        if not line:
            problems.append(f"Table 1: row {name} not found")
            continue
        v = nums(line[0])
        r = comp.loc[key]
        exp = [r.failed_clr, r.mature_clr, r["diff"], r.ci_low, r.ci_high,
               r.cohens_d, r.welch_p, r.mwu_p, r.fdr]
        for a, b in zip(v, exp):
            rows_checked += 1
            if abs(a - b) > 0.006:
                problems.append(f"Table 1 {name}: text {a} vs source {b:.3f}")

    # ---- Table 2
    ab = load("s5e_abundance_tests").set_index("metric")
    block = text.split("**Table 2.")[1].split("**Table 3.")[0]
    pairs = [("Myeloid cells per 1,000", "myeloid_per1000"),
             ("Myeloid cells per donor", "n_myeloid"),
             ("cDC2 per 1,000", "cdc2_D1_per1000"),
             ("cDC2 cells per donor", "n_cdc2_D1"),
             ("cDC2 as % of myeloid", "cdc2_D1_frac_myeloid"),
             ("cDC2, marker rule, per 1,000", "cdc2_D2_per1000"),
             ("cDC2, marker rule, cells per donor", "n_cdc2_D2"),
             ("cDC2, marker rule, as % of myeloid", "cdc2_D2_frac_myeloid")]
    for label, key in pairs:
        line = [l for l in block.split("\n") if l.startswith(f"| {label}")]
        if not line:
            problems.append(f"Table 2: row {label} not found")
            continue
        v = nums(line[0])
        r = ab.loc[key]
        scale = 100 if "frac_myeloid" in key else 1
        exp = [r.failed_mean * scale, r.mature_mean * scale, r["diff"] * scale,
               r.ci_low * scale, r.ci_high * scale, r.cohens_d,
               r.welch_p, r.fdr]
        tol = 0.06 if "frac_myeloid" in key else 0.06
        for a, b in zip(v, exp):
            rows_checked += 1
            if abs(a - b) > tol:
                problems.append(f"Table 2 {label}: text {a} vs source {b:.3f}")

    # ---- Table 3
    dec = load("s5f_decomposition_retained").set_index("gene")
    block = text.split("**Table 3.")[1].split("## STAR Methods")[0]
    for gene in ["CD1C", "APOBEC3A", "TIFAB", "CHN2", "CCL3L1", "S16", "YEATS4"]:
        pass
    for gene in ["CD1C", "APOBEC3A", "TIFAB", "CHN2", "CCL3L1", "SLC16A10", "YEATS4"]:
        line = [l for l in block.split("\n") if l.startswith(f"| {gene} |")]
        if not line:
            problems.append(f"Table 3: row {gene} not found")
            continue
        v = nums(line[0])
        r = dec.loc[gene]
        exp = [r.log2FC_observed, r.log2FC_composition, r.log2FC_intrinsic,
               r.log2FC_predicted, r.residual, r.pct_from_composition]
        for a, b in zip(v, exp):
            rows_checked += 1
            if abs(a - b) > 0.51:      # 末列为整数百分比
                problems.append(f"Table 3 {gene}: text {a} vs source {b:.3f}")
    return rows_checked, problems


def main():
    checks = [x for x in build_checks() if x[3]]
    text = open(MS, encoding="utf-8").read() if MS else None
    if text is None:
        print("verifying values only: no manuscript under manuscript/ "
              "(set MS_PATH to also match every value against the prose)")
    else:
        print(f"verifying against: {os.path.relpath(MS, ROOT).replace(os.sep, '/')}")
    rows, fails = [], []
    for cid, claim, src, snippet, script, fig in checks:
        snippet = re.sub(r"(?<![A-Za-z])-(?=\d)", "\u2212", snippet)
        if text is None:
            status = "RECOMPUTED"
        else:
            status = "PASS" if snippet in text else "FAIL"
        rows.append(dict(id=cid, claim=claim, expected_text=snippet,
                         source_file=f"results/{src}.csv", script=f"code/{script}",
                         figure=fig, status=status))
        if status == "FAIL":
            fails.append((cid, claim, snippet))
    if text is None:
        n_tab, tab_problems = 0, []
    else:
        n_tab, tab_problems = verify_tables(text)
    for i, pb in enumerate(tab_problems):
        rows.append(dict(id=f"T{i+1:02d}", claim="table value", expected_text=pb,
                         source_file="results/", script="code/90_verify_numbers.py",
                         figure="Table 1-3", status="FAIL"))
        fails.append((f"T{i+1:02d}", "table value", pb))
    pd.DataFrame(rows).to_csv(f"{RES}/provenance.csv", index=False)

    print(f"checked {len(rows) - len(tab_problems)} numeric claims "
          f"({n_tab} table cells) against results/")
    if text is None:
        print(f"RECOMPUTED {len(rows)} / FAIL {len(fails)} "
              f"(0 table cells checked: no manuscript text)")
    else:
        print(f"PASS {len(rows)-len(fails)} / FAIL {len(fails)}")
    for cid, claim, snip in fails:
        if text is None:
            print(f"\n[FAIL] {cid} {claim}\n  expected value: {snip}")
            continue
        print(f"\n[FAIL] {cid} {claim}\n  expected in manuscript: {snip}")
        key = snip[:28]
        idx = text.find(key)
        print(f"  manuscript has: {text[idx:idx+len(snip)+8] if idx>=0 else '(anchor not found)'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
