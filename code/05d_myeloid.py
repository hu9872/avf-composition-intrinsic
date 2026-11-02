"""S5d v2: 髓系区室亚群细分 —— 检验 CD1C+ cDC2 是否为整群丢失。

上一版注释失败的原因已定位：髓系区室以巨噬细胞为主体，Wilcoxon 秩和前 40 名标记
被核糖体/肌动蛋白等管家基因占据，面板命中数上不去，导致 76% 细胞被判 Mixed、
cDC2 归零。本版改用两套互相独立、互不依赖的判定：

  判读 1（聚类层）：逐细胞 score_genes + 簇内多数投票，面板只用谱系特异基因
  判读 2（无注释层）：cDC2 签名分数阈值直接定义细胞，完全绕过聚类

四个互斥假设的正交检验：
  A 群体层：供体水平 cDC2 占比（以髓系为分母 + 以样本全部细胞为分母）
  B 转录层：仅在 cDC2 细胞内部比较 CD1C 检出率与表达量（含 CD1C 无关身份签名）
  C 状态层：DC 区室再聚类，检索 CD1C-low cDC2 亚状态
  D 表达层：cDC2 内部 pseudobulk 差异表达
  E 归一化：排除解离/采样偏倚
"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats
from statsmodels.stats.multitest import multipletests

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
GROUPS = ["mature_AVF", "failed_AVF"]

# 只用谱系特异基因；刻意剔除 CD74/HLA-DP/CST3/ITGAX 等共享基因，
# 否则单核-巨噬的强信号会把 DC 面板整体抬高
PANELS = {
    "Mono_classical": ["CD14", "FCN1", "S100A8", "S100A9", "S100A12", "CCR2", "VCAN",
                       "MNDA", "EREG", "CSF3R"],
    "Mono_nonclassical": ["FCGR3A", "MS4A7", "CDKN1C", "CX3CR1", "LILRB1", "LST1"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2", "SEPP1",
                   "SELENOP", "MAF", "NR1H3"],
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "CADM1", "WDFY4"],
    "cDC2": ["CD1C", "CLEC10A", "FCER1A", "CD1E", "CD1D"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "PLD4", "TCF4", "IRF7"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"],
    "Neutrophil": ["FCGR3B", "CXCR2", "FPR1"],
}
# CD1C 无关的 cDC2 身份基因：用于在 CD1C 沉默时仍能判断群体是否还在
CDC2_FREE = ["CLEC10A", "FCER1A", "CD1E", "CD1D"]
DC_LIKE = ["cDC1", "cDC2", "pDC"]


def clr(mat, pseudo=0.5):
    x = mat + pseudo
    lg = np.log(x)
    return lg.sub(lg.mean(axis=1), axis=0)


def donor_test(df, col, group_col="group", g1="failed_AVF", g2="mature_AVF", min_n=4):
    sub = df[[col, group_col]].dropna()
    a = sub.loc[sub[group_col] == g1, col].values.astype(float)
    b = sub.loc[sub[group_col] == g2, col].values.astype(float)
    if len(a) < min_n or len(b) < min_n:
        return None
    t, p = stats.ttest_ind(a, b, equal_var=False)
    try:
        u, pu = stats.mannwhitneyu(a, b)
    except ValueError:
        pu = np.nan
    sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    d = (a.mean() - b.mean()) / sp if sp > 0 else np.nan
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return dict(failed_mean=a.mean(), mature_mean=b.mean(),
                diff=a.mean() - b.mean(), cohens_d=d, welch_p=p, mwu_p=pu,
                ci_low=(a.mean() - b.mean()) - 1.96 * se,
                ci_high=(a.mean() - b.mean()) + 1.96 * se,
                n_failed=len(a), n_mature=len(b),
                values_failed=";".join(f"{v:.3f}" for v in sorted(a)),
                values_mature=";".join(f"{v:.3f}" for v in sorted(b)))


def main():
    adata = ad.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    n_all = adata.obs.groupby("donor", observed=True).size().rename("n_all_cells")

    mye = adata[adata.obs["celltype"] == "Myeloid"].copy()
    mye.X = mye.layers["lognorm"].copy()
    print("myeloid cells:", mye.shape)

    bad = mye.var_names.str.match(r"^(IG[HKL][VDJ]|IGHG|IGKC|IGLC|HBB|HBA1|HBA2|ALAS2)")
    mye = mye[:, ~bad].copy()
    print("after gene filter:", mye.shape)

    sc.pp.highly_variable_genes(mye, n_top_genes=1500, batch_key="donor", flavor="seurat")
    hv = mye[:, mye.var.highly_variable].copy()
    sc.pp.pca(hv, n_comps=40, zero_center=False, random_state=SEED)
    import harmonypy  # noqa: F401
    sc.external.pp.harmony_integrate(hv, "donor", random_state=SEED, max_iter_harmony=20)
    mye.obsm["X_pca"] = hv.obsm["X_pca_harmony"]
    sc.pp.neighbors(mye, n_neighbors=15, n_pcs=30, random_state=SEED)
    for r, key in [(0.5, "mye_0.5"), (0.9, "mye_0.9")]:
        sc.tl.leiden(mye, resolution=r, key_added=key, random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(mye, random_state=SEED)
    print("clusters @0.9:", mye.obs["mye_0.9"].nunique())

    # ------------------------------------------------ 判读 1：逐细胞打分 + 多数投票
    print("\n=== 判读 1：score_genes ===")
    used = {}
    for k, gl in PANELS.items():
        g = [x for x in gl if x in mye.var_names]
        used[k] = g
        sc.tl.score_genes(mye, g, score_name=f"sg_{k}", ctrl_size=50, n_bins=25,
                          use_raw=False, random_state=SEED)
    sgc = [f"sg_{k}" for k in PANELS]
    S = mye.obs[sgc].copy()
    S.columns = list(PANELS)
    # 簇内 z-score（跨簇标准化）：消除不同面板基因数/表达量级差异
    cl_mean = S.groupby(mye.obs["mye_0.9"], observed=True).mean()
    cl_z = (cl_mean - cl_mean.mean()) / (cl_mean.std() + 1e-9)
    cl_z.to_csv(f"{RES}/s5d_score_zscores.csv")

    # 逐细胞 argmax -> 簇内多数投票
    cell_call = S.idxmax(1)
    mye.obs["cell_call"] = cell_call.values
    vote = (pd.crosstab(mye.obs["mye_0.9"], mye.obs["cell_call"], normalize="index"))
    vote.to_csv(f"{RES}/s5d_cluster_vote.csv")
    labels, flags, rows = {}, {}, []
    for cl in cl_z.index:
        zrow = cl_z.loc[cl].sort_values(ascending=False)
        majority = vote.loc[cl].idxmax()
        agree = float(vote.loc[cl].max())
        # 簇标签取「z 值最高」与「簇内多数投票」的交集；冲突时以 z 值为主并记录
        lab = zrow.index[0]
        if majority != lab and abs(zrow.iloc[0] - zrow[majority]) < 0.5:
            lab = majority
        labels[cl] = lab
        flags[cl] = dict(assigned=lab, z_max=round(float(zrow.iloc[0]), 2),
                         z_runner=round(float(zrow.iloc[1]), 2),
                         z_by_panel={k: round(float(v), 2) for k, v in zrow.items()},
                         majority_call=majority, majority_agreement=round(agree, 3),
                         n=int((mye.obs["mye_0.9"] == cl).sum()))
        rows.append(f"cluster {cl:>2}: {lab:16s} z1={zrow.index[0]}:{zrow.iloc[0]:5.2f} "
                    f"z2={zrow.index[1]}:{zrow.iloc[1]:5.2f} "
                    f"vote={majority}({agree:.0%}) n={flags[cl]['n']}")
    open(f"{RES}/s5d_annotation_check.txt", "w", encoding="utf8").write("\n".join(rows))
    print("\n".join(rows))

    mye.obs["mye_subtype"] = mye.obs["mye_0.9"].map(labels).astype(str)
    print("\nsubtype counts:")
    print(mye.obs["mye_subtype"].value_counts().to_string())

    # ------------------------------------------------ 判读 2：无注释的 cDC2 量化
    # 以 CLEC10A/FCER1A/CD1E/CD1D（不含 CD1C）的签名分数定义 cDC2，
    # 与 CD1C 本身解耦，用于判别「群体丢失」vs「基因沉默」
    sc.tl.score_genes(mye, [g for g in CDC2_FREE if g in mye.var_names],
                      score_name="sg_cDC2_free", ctrl_size=50, n_bins=25,
                      use_raw=False, random_state=SEED)
    free_score = mye.obs["sg_cDC2_free"].values
    # 阈值：混合分布双峰之间的谷底（在 0 附近）；用 KDE 找局部极小
    hist, edges = np.histogram(free_score, bins=120, range=(-1, 2))
    mid = (edges[:-1] + edges[1:]) / 2
    lo = np.where((mid > -0.3) & (mid < 0.8))[0]
    thr = float(mid[lo[np.argmin(hist[lo])]]) if len(lo) else 0.0
    mye.obs["is_cdc2_free"] = (free_score > thr).astype(int)
    print(f"\n判读 2: cDC2 无关签名阈值 = {thr:.3f}，命中 {int(mye.obs.is_cdc2_free.sum())} 细胞")
    pd.DataFrame({"score_mid": mid, "count": hist}).to_csv(
        f"{RES}/s5d_cdc2_score_hist.csv", index=False)

    # ------------------------------------------------ 检验 E：绝对丰度
    cnt = (mye.obs.groupby(["donor", "group", "mye_subtype"], observed=True)
           .size().unstack(fill_value=0))
    cnt.to_csv(f"{RES}/s5d_subtype_counts.csv")
    per_k = cnt.div(n_all.reindex(cnt.index.get_level_values("donor")).values, axis=0) * 1000
    per_k.to_csv(f"{RES}/s5d_subtype_per1000.csv")

    # 无注释法的每千细胞计数
    free_cnt = (mye.obs.groupby(["donor", "group"], observed=True)["is_cdc2_free"]
                .sum().rename("n_cdc2_free").reset_index())
    free_cnt["per1000"] = (free_cnt.n_cdc2_free /
                           n_all.reindex(free_cnt.donor).values * 1000)
    free_cnt.to_csv(f"{RES}/s5d_cdc2_free_per1000.csv", index=False)

    # ------------------------------------------------ 检验 A：组成（CLR，供体为单位）
    post = cnt.loc[cnt.index.get_level_values("group").isin(GROUPS)]
    cl = clr(post.astype(float))
    cl.to_csv(f"{RES}/s5d_subtype_clr.csv")

    rows = []
    for s in post.columns:
        r = donor_test(cl.reset_index(), s)
        if r is None:
            continue
        r["subtype"] = s
        rows.append(r)
    comp = pd.DataFrame(rows).set_index("subtype")
    comp["fdr"] = multipletests(comp.welch_p, method="fdr_bh")[1]
    comp.to_csv(f"{RES}/s5d_composition_tests.csv")
    print("\n=== A. 供体水平组成（CLR，分母=髓系） ===")
    print(comp[["failed_mean", "mature_mean", "diff", "cohens_d", "welch_p",
                "mwu_p", "fdr"]].round(3).to_string())

    rows = []
    for s in per_k.columns:
        r = donor_test(per_k.reset_index(), s)
        if r is None:
            continue
        r["subtype"] = s
        rows.append(r)
    pk = pd.DataFrame(rows).set_index("subtype")
    pk["fdr"] = multipletests(pk.welch_p, method="fdr_bh")[1]
    pk.to_csv(f"{RES}/s5d_per1000_tests.csv")
    print("\n=== E. 每千细胞绝对丰度（分母=样本全部细胞） ===")
    print(pk[["failed_mean", "mature_mean", "diff", "cohens_d", "welch_p", "fdr"]]
          .round(3).to_string())

    # 无注释法 cDC2 丰度检验（最关键的一条）
    r_free = donor_test(free_cnt, "per1000")
    r_free_n = donor_test(free_cnt, "n_cdc2_free")
    pd.DataFrame([dict(metric="cDC2_free_per1000", **r_free),
                  dict(metric="cDC2_free_count", **r_free_n)]).to_csv(
        f"{RES}/s5d_cdc2_free_test.csv", index=False)
    print("\n=== A2. 无注释法 cDC2 丰度（供体水平） ===")
    for nm, r in [("per1000", r_free), ("count", r_free_n)]:
        print(f"  {nm}: failed={r['failed_mean']:.2f} mature={r['mature_mean']:.2f} "
              f"d={r['cohens_d']:.2f} welch_p={r['welch_p']:.4f} mwu_p={r['mwu_p']:.4f}")
        print(f"      failed values: {r['values_failed']}")
        print(f"      mature values: {r['values_mature']}")

    mye_tot = mye.obs.groupby(["donor", "group"], observed=True).size().rename("n").reset_index()
    mye_tot = pd.merge(mye_tot, n_all.reset_index(), on="donor")
    mye_tot["per1000"] = mye_tot["n"] / mye_tot["n_all_cells"] * 1000
    mye_tot.to_csv(f"{RES}/s5d_myeloid_total.csv", index=False)
    r_tot = donor_test(mye_tot, "n")
    r_totp = donor_test(mye_tot, "per1000")
    print("\n髓系总量(细胞数): d=%.2f p=%.4f | 每千细胞: d=%.2f p=%.4f" %
          (r_tot["cohens_d"], r_tot["welch_p"], r_totp["cohens_d"], r_totp["welch_p"]))
    pd.DataFrame([dict(metric="myeloid_count", **r_tot),
                  dict(metric="myeloid_per1000", **r_totp)]).to_csv(
        f"{RES}/s5d_myeloid_total_test.csv", index=False)

    # ------------------------------------------------ 检验 B：cDC2 内部 CD1C 水平
    is_cdc2 = (mye.obs["mye_subtype"] == "cDC2").values
    cdc2 = mye[is_cdc2].copy()
    print("\ncDC2 cells (cluster-based):", cdc2.n_obs,
          "| donors:", cdc2.obs["donor"].nunique())

    if cdc2.n_obs > 100:
        gb = ["CD1C", "CLEC10A", "FCER1A", "CD1E", "CD1D"]
        gb = [g for g in gb if g in cdc2.var_names]
        Cc = np.asarray(cdc2[:, gb].layers["counts"].todense())
        Ln = np.asarray(cdc2[:, gb].layers["lognorm"].todense())
        det = pd.DataFrame((Cc > 0).astype(float),
                           columns=[f"det_{g}" for g in gb], index=cdc2.obs_names)
        det["donor"] = cdc2.obs["donor"].values
        det["group"] = cdc2.obs["group"].values
        det["n_cells"] = 1
        det_rate = det.groupby(["donor", "group"], observed=True).mean()
        det_rate.to_csv(f"{RES}/s5d_cdc2_detection_rate.csv")

        expr = pd.DataFrame(Ln, columns=[f"exp_{g}" for g in gb], index=cdc2.obs_names)
        expr["donor"] = cdc2.obs["donor"].values
        expr["group"] = cdc2.obs["group"].values
        expr_mean = expr.groupby(["donor", "group"], observed=True).mean()
        expr_mean.to_csv(f"{RES}/s5d_cdc2_mean_expr.csv")

        ncells = det.groupby(["donor", "group"], observed=True)["n_cells"].sum()

        rows = []
        for g in gb:
            for metric, tbl, col in [("detection", det_rate, f"det_{g}"),
                                     ("mean_lognorm", expr_mean, f"exp_{g}")]:
                r = donor_test(tbl.reset_index().merge(
                    ncells.rename("n").reset_index(), on=["donor", "group"]), col)
                if r is None:
                    continue
                r.update(gene=g, metric=metric)
                rows.append(r)
        bt = pd.DataFrame(rows)
        bt["fdr"] = multipletests(bt.welch_p, method="fdr_bh")[1]
        bt.to_csv(f"{RES}/s5d_cdc2_cd1c_test.csv", index=False)
        print("\n=== B. cDC2 内部 CD1C 检出率与表达（供体水平） ===")
        print(bt[["gene", "metric", "failed_mean", "mature_mean", "diff",
                  "cohens_d", "welch_p", "fdr"]].round(4).to_string(index=False))

        # CD1C 无关身份签名：若两组相当 -> 群体在、只是 CD1C 下调（H2/H3）
        efm = expr_mean[[f"exp_{g}" for g in ["CLEC10A", "FCER1A", "CD1E", "CD1D"]
                         if f"exp_{g}" in expr_mean.columns]].copy()
        efm["cd1c_free_identity"] = efm.mean(1)
        efm.to_csv(f"{RES}/s5d_cdc2_identity_signature.csv")
        r = donor_test(efm.reset_index(), "cd1c_free_identity")
        print("\ncDC2 身份签名（CLEC10A/FCER1A/CD1E/CD1D，不含 CD1C）: "
              f"failed={r['failed_mean']:.3f} mature={r['mature_mean']:.3f} "
              f"d={r['cohens_d']:.2f} p={r['welch_p']:.4f}")

        # cDC2 细胞中 CD1C+ 的比例
        posf = pd.DataFrame({"donor": cdc2.obs["donor"].values,
                             "group": cdc2.obs["group"].values,
                             "cd1c_pos": (np.asarray(cdc2[:, ["CD1C"]].layers["counts"]
                                                     .todense()).ravel() > 0).astype(float)})
        posr = posf.groupby(["donor", "group"], observed=True)["cd1c_pos"].mean()
        posr.to_csv(f"{RES}/s5d_cdc2_cd1c_positive_fraction.csv")
        r = donor_test(posr.reset_index(), "cd1c_pos")
        print(f"cDC2 内 CD1C+ 细胞比例: failed={r['failed_mean']:.3f} "
              f"mature={r['mature_mean']:.3f} d={r['cohens_d']:.2f} p={r['welch_p']:.4f}")

    # ------------------------------------------------ 检验 C：DC 区室再聚类
    is_dc = mye.obs["mye_subtype"].isin(DC_LIKE).values
    dc = mye[is_dc].copy()
    print("\nDC-like cells:", dc.n_obs, "| by subtype:",
          dc.obs["mye_subtype"].value_counts().to_dict())
    if dc.n_obs > 200:
        sc.pp.highly_variable_genes(dc, n_top_genes=800, batch_key="donor", flavor="seurat")
        dhv = dc[:, dc.var.highly_variable].copy()
        sc.pp.pca(dhv, n_comps=25, zero_center=False, random_state=SEED)
        sc.external.pp.harmony_integrate(dhv, "donor", random_state=SEED, max_iter_harmony=15)
        dc.obsm["X_pca"] = dhv.obsm["X_pca_harmony"]
        sc.pp.neighbors(dc, n_neighbors=15, n_pcs=20, random_state=SEED)
        sc.tl.leiden(dc, resolution=0.8, key_added="dc_0.8", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
        sc.tl.rank_genes_groups(dc, "dc_0.8", method="wilcoxon", use_raw=False)
        sc.get.rank_genes_groups_df(dc, None).to_csv(f"{RES}/s5d_dc_markers.csv", index=False)

        Cc = np.asarray(dc[:, ["CD1C"]].layers["counts"].todense()).ravel()
        dc.obs["cd1c_pos"] = (Cc > 0).astype(int)
        posf = dc.obs.groupby("dc_0.8", observed=True)["cd1c_pos"].mean()
        cellof = pd.crosstab(dc.obs["dc_0.8"], dc.obs["mye_subtype"])
        summary = pd.DataFrame({"n": dc.obs["dc_0.8"].value_counts(),
                                "cd1c_pos_frac": posf})
        summary = summary.join(cellof.add_prefix("n_"))
        summary.to_csv(f"{RES}/s5d_dc_cluster_summary.csv")
        print(summary.round(3).to_string())

        cnt_dc = (dc.obs.groupby(["donor", "group", "dc_0.8"], observed=True)
                  .size().unstack(fill_value=0))
        cnt_dc.to_csv(f"{RES}/s5d_dc_cluster_counts.csv")
        cl_dc = clr(cnt_dc.astype(float))
        rows = []
        for s in cnt_dc.columns:
            r = donor_test(cl_dc.reset_index(), s)
            if r is None:
                continue
            r["dc_cluster"] = s
            rows.append(r)
        if rows:
            dd = pd.DataFrame(rows).set_index("dc_cluster")
            dd["fdr"] = multipletests(dd.welch_p, method="fdr_bh")[1]
            dd.to_csv(f"{RES}/s5d_dc_cluster_tests.csv")
            print("\n=== C. DC 亚簇供体水平组成 ===")
            print(dd[["failed_mean", "mature_mean", "diff", "cohens_d",
                      "welch_p", "fdr"]].round(3).to_string())

    # ------------------------------------------------ 检验 D：cDC2 内部 pseudobulk
    post_mask = mye.obs["group"].isin(GROUPS) & (mye.obs["mye_subtype"] == "cDC2")
    if post_mask.sum() > 300:
        sub = mye[post_mask.values].copy()
        X = sub.layers["counts"].tocsr()
        obs = sub.obs
        recs, meta = [], []
        for (donor, group), idx in obs.groupby(["donor", "group"], observed=True).indices.items():
            if len(idx) < 10:
                print(f"  skip {donor} ({len(idx)} cells)")
                continue
            recs.append(np.asarray(X[np.isin(sub.obs_names, obs.index[idx])].sum(0)).ravel())
            meta.append(dict(sample=donor, donor=donor, group=group, n=len(idx)))
        if len(meta) >= 8 and len({m["group"] for m in meta}) == 2:
            pb = ad.AnnData(X=np.vstack(recs).astype(np.float32),
                            obs=pd.DataFrame(meta).set_index("sample"),
                            var=pd.DataFrame(index=mye.var_names))
            pb.obs["group"] = pd.Categorical(pb.obs["group"], categories=GROUPS)
            print("\ncDC2 pseudobulk: donors =", len(meta),
                  pb.obs.group.value_counts().to_dict())
            try:
                from pydeseq2.dds import DeseqDataSet
                from pydeseq2.ds import DeseqStats
                dds = DeseqDataSet(counts=pb.X, metadata=pb.obs, design="~group",
                                   refit_cooks=True, quiet=True)
                dds.deseq2()
                st = DeseqStats(dds, contrast=["group", "failed_AVF", "mature_AVF"], quiet=True)
                st.summary()
                d = st.results_df.copy()
                d.index.name = "gene"
                d.to_csv(f"{RES}/s5d_cdc2_pseudobulk_de.csv")
                sig = d[(d.padj < 0.05) & (d.log2FoldChange.abs() > 0.585)]
                print("cDC2 内部 DE 基因数:", len(sig))
                focus = [g for g in ["CD1C", "CLEC10A", "FCER1A", "CD1E", "CD1D"]
                         if g in d.index]
                print(d.loc[focus, ["log2FoldChange", "padj"]].round(4).to_string())
            except Exception as e:
                print("DESeq2 失败:", type(e).__name__, str(e)[:150])

    # ------------------------------------------------ 供体明细（防单供体驱动）
    det = pd.crosstab([mye.obs["donor"], mye.obs["group"]], mye.obs["mye_subtype"])
    det["__total__"] = det.sum(1)
    det.to_csv(f"{RES}/s5d_donor_by_subtype.csv")
    print("\n=== 供体 x 亚型细胞数 ===")
    print(det.to_string())

    mye.write_h5ad(os.path.join(WORK, "05d_myeloid.h5ad"), compression="gzip")
    json.dump({"labels": labels, "decision": {k: {kk: (vv if not isinstance(vv, dict)
                                                      else vv)
                                                  for kk, vv in v.items()}
                                              for k, v in flags.items()},
               "cdc2_free_threshold": thr,
               "n_myeloid": int(mye.n_obs)},
              open(f"{RES}/s5d_annotation.json", "w"), indent=2, ensure_ascii=False,
              default=float)
    print("\nwritten:", os.path.join(WORK, "05d_myeloid.h5ad"))


if __name__ == "__main__":
    main()
