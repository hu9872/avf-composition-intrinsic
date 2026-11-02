"""S5e v2: cDC2 群体存续性裁决。

修正 05d 的三处缺陷：
  (1) z 值 argmax 无绝对门槛 -> CLEC9A 实际最高仅 0.111 却判出 3,634 个 cDC1。改用门控注释。
  (2) 降采样时把「细胞总 UMI」误算成「面板基因计数和」，导致全部细胞被判深度不足。
      改为以细胞真实 total_counts 为分母做多项分布重抽样。
  (3) GMM 等密度点落在分布外。改为「标记规则」定义敏感性分析，阈值可解释。

cDC2 的两种独立定义：
  D1 门控聚类注释（cluster 水平，marker 实测校验）
  D2 标记规则（单细胞水平：{CD1C, CLEC10A, FCER1A, CD1E} 中 ≥2 个检出）
两者一致性 + 各自独立做丰度检验。
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
TARGET_DEPTH = 5000
CTRL_GENES = ["RPL13A", "RPS18", "ACTB", "B2M", "TMSB10", "EEF1A1", "GAPDH"]
CDC2_MARKERS = ["CD1C", "CLEC10A", "FCER1A", "CD1E"]


def clr(mat, pseudo=0.5):
    x = mat + pseudo
    lg = np.log(x)
    return lg.sub(lg.mean(axis=1), axis=0)


def fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def donor_test(df, col, group_col="group", g1="failed_AVF", g2="mature_AVF", min_n=4):
    sub = df[[col, group_col]].dropna()
    a = sub.loc[sub[group_col] == g1, col].values.astype(float)
    b = sub.loc[sub[group_col] == g2, col].values.astype(float)
    if len(a) < min_n or len(b) < min_n:
        return None
    if (len(a) == len(b) and np.allclose(a, b)) or (np.nanstd(a) == 0 and np.nanstd(b) == 0):
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


def gated_annotation(cm_row):
    g = cm_row
    if g.get("TPSAB1", 0) >= 0.5 or g.get("CPA3", 0) >= 0.5:
        return "Mast"
    if g.get("S100A8", 0) >= 1.5 and g.get("CD14", 0) < 0.5 and g.get("C1QA", 0) < 0.5:
        return "Neutrophil"
    if g.get("CD1C", 0) >= 0.4 or (g.get("CLEC10A", 0) >= 0.6 and g.get("C1QA", 0) < 1.6):
        return "cDC2"
    if g.get("CLEC9A", 0) >= 0.3 or g.get("XCR1", 0) >= 0.2:
        return "cDC1"
    if g.get("LILRA4", 0) >= 0.15 and g.get("CLEC4C", 0) >= 0.08:
        return "pDC"
    if g.get("C1QA", 0) >= 1.0:
        if g.get("CD14", 0) >= 1.5 and g.get("C1QA", 0) < 2.2 and g.get("FCN1", 0) >= 0.15:
            return "Mono_classical"
        return "Macrophage"
    if g.get("CD14", 0) >= 1.0:
        return "Mono_classical"
    if g.get("FCGR3A", 0) >= 1.5:
        return "Mono_FCGR3A"
    return "Low_quality"


def main():
    a = ad.read_h5ad(os.path.join(WORK, "05d_myeloid.h5ad"))
    cm = pd.read_csv(f"{RES}/s5d_cluster_meanexpr.csv", index_col=0)
    cm.index = cm.index.astype(str)

    # ---------------------------------------------------------- 1. 门控注释
    print("=== 1. 绝对表达门控注释 ===")
    lab = {cl: gated_annotation(cm.loc[cl]) for cl in cm.index}
    a.obs["mye_subtype"] = a.obs["mye_0.9"].astype(str).map(lab).fillna("Low_quality")
    vc = a.obs["mye_subtype"].value_counts()
    print(vc.to_string())
    print("  cDC1 命中数:", int(vc.get("cDC1", 0)), "(CLEC9A/XCR1 未过门槛，判定为未检出)")

    check = []
    for s in vc.index:
        m = (a.obs["mye_subtype"] == s).values
        row = {"subtype": s, "n": int(m.sum())}
        for g in ["CD1C", "CLEC10A", "FCER1A", "C1QA", "CD14", "FCGR3A",
                  "CLEC9A", "XCR1", "LILRA4", "TPSAB1"]:
            if g in a.var_names:
                X = np.asarray(a[m, [g]].layers["lognorm"].todense()).ravel()
                row[g] = round(float(X.mean()), 3)
        check.append(row)
    ck = pd.DataFrame(check).set_index("subtype").sort_values("n", ascending=False)
    ck.to_csv(f"{RES}/s5e_subtype_marker_check.csv")
    print(ck.to_string())

    absent = {}
    for g in ["CLEC9A", "XCR1", "LILRA4", "CLEC4C"]:
        if g in a.var_names:
            X = np.asarray(a[:, [g]].layers["counts"].todense()).ravel()
            absent[g] = dict(pct_positive=round(float((X > 0).mean() * 100), 2),
                             mean_counts=round(float(X.mean()), 3))
    print("\n关键标记阳性率（全部髓系）:", absent)
    json.dump(absent, open(f"{RES}/s5e_absent_markers.json", "w"), indent=2)

    # ---------------------------------------------------------- 2. 两种 cDC2 定义
    print("\n=== 2. cDC2 的两种独立定义 ===")
    d1 = (a.obs["mye_subtype"] == "cDC2").values
    Mk = np.vstack([np.asarray(a[:, [g]].layers["counts"].todense()).ravel() > 0
                    for g in CDC2_MARKERS if g in a.var_names]).T
    npos = Mk.sum(1)
    for k in [1, 2, 3]:
        print(f"  ≥{k}/4 标记阳性: {int((npos >= k).sum())} 细胞")
    d2 = npos >= 2
    a.obs["is_cDC2_cluster"] = d1.astype(int)
    a.obs["is_cDC2_marker"] = d2.astype(int)
    a.obs["cdc2_marker_n"] = npos

    ov = pd.crosstab(pd.Series(d1, name="D1_cluster"),
                     pd.Series(d2, name="D2_marker"))
    ov.to_csv(f"{RES}/s5e_cdc2_definition_overlap.csv")
    print("\n  D1(聚类) × D2(标记规则) 交叉表:")
    print(ov.to_string())
    inter = int((d1 & d2).sum())
    print(f"  交集 {inter}；D1 内被 D2 支持的比例 {inter / max(d1.sum(),1):.1%}")

    pur = []
    for g in ["CD1C", "CLEC10A", "FCER1A", "CD1E", "C1QA", "CD14", "LYZ"]:
        if g in a.var_names:
            C = np.asarray(a[:, [g]].layers["counts"].todense()).ravel()
            pur.append(dict(gene=g,
                            pct_pos_D1=round(float((C[d1] > 0).mean() * 100), 1),
                            pct_pos_D2=round(float((C[d2] > 0).mean() * 100), 1),
                            pct_pos_rest=round(float((C[~(d1 | d2)] > 0).mean() * 100), 1)))
    pu = pd.DataFrame(pur)
    pu.to_csv(f"{RES}/s5e_cdc2_purity.csv", index=False)
    print("\n  纯度核验（阳性率 %）:")
    print(pu.to_string(index=False))

    # ---------------------------------------------------------- 3. 群体层丰度
    print("\n=== 3. 群体层丰度（供体为单位） ===")
    full = ad.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"), backed="r")
    n_all = full.obs.groupby("donor", observed=True).size().rename("n_all_cells")
    full.file.close()

    obs = a.obs
    grp = obs.groupby(["donor", "group"], observed=True)
    comp = pd.DataFrame({"n_myeloid": grp.size(),
                         "n_cdc2_D1": grp["is_cDC2_cluster"].sum(),
                         "n_cdc2_D2": grp["is_cDC2_marker"].sum()}).reset_index()
    ntot = n_all.reindex(comp.donor).values
    comp["cdc2_D1_frac_myeloid"] = comp.n_cdc2_D1 / comp.n_myeloid
    comp["cdc2_D2_frac_myeloid"] = comp.n_cdc2_D2 / comp.n_myeloid
    comp["cdc2_D1_per1000"] = comp.n_cdc2_D1 / ntot * 1000
    comp["cdc2_D2_per1000"] = comp.n_cdc2_D2 / ntot * 1000
    comp["myeloid_per1000"] = comp.n_myeloid / ntot * 1000
    comp.to_csv(f"{RES}/s5e_donor_abundance.csv", index=False)
    post = comp[comp.group.isin(GROUPS)].copy()

    rows = []
    for metric in ["n_cdc2_D1", "cdc2_D1_frac_myeloid", "cdc2_D1_per1000",
                   "n_cdc2_D2", "cdc2_D2_frac_myeloid", "cdc2_D2_per1000",
                   "n_myeloid", "myeloid_per1000"]:
        r = donor_test(post, metric)
        if r:
            r["metric"] = metric
            rows.append(r)
    ab = pd.DataFrame(rows).set_index("metric")
    ab["fdr"] = fdr(ab.welch_p)
    ab.to_csv(f"{RES}/s5e_abundance_tests.csv")
    print(ab[["failed_mean", "mature_mean", "diff", "cohens_d", "welch_p",
              "mwu_p", "fdr"]].round(4).to_string())
    print("\n逐供体明细:")
    print(post[["donor", "group", "n_myeloid", "n_cdc2_D1", "cdc2_D1_frac_myeloid",
                "cdc2_D1_per1000", "n_cdc2_D2"]].round(3).to_string(index=False))

    # 各亚型供体水平组成（CLR，分母=髓系）
    cnt = pd.crosstab([obs.donor, obs.group], obs.mye_subtype)
    cnt.to_csv(f"{RES}/s5e_subtype_counts.csv")
    cntp = cnt.loc[cnt.index.get_level_values("group").isin(GROUPS)]
    cl = clr(cntp.astype(float))
    rows = []
    for s in cntp.columns:
        r = donor_test(cl.reset_index(), s)
        if r:
            r["subtype"] = s
            rows.append(r)
    ctest = pd.DataFrame(rows).set_index("subtype")
    ctest["fdr"] = fdr(ctest.welch_p)
    ctest.to_csv(f"{RES}/s5e_subtype_composition_tests.csv")
    print("\n=== 3b. 各亚型供体水平组成（CLR，分母=髓系） ===")
    print(ctest[["failed_mean", "mature_mean", "diff", "cohens_d", "welch_p",
                 "fdr"]].round(3).to_string())

    # ---------------------------------------------------------- 4. 深度控制的 CD1C 检验
    print("\n=== 4. 等深度降采样（5,000 UMI）后的 CD1C 检验 ===")
    cdc2 = a[obs["is_cDC2_cluster"].values == 1].copy()
    print(f"  cDC2(D1) 细胞 {cdc2.n_obs}")
    total_full = cdc2.obs["total_counts"].values.astype(float)
    keep = total_full >= TARGET_DEPTH
    print(f"  深度 ≥{TARGET_DEPTH}: {keep.sum()} / {cdc2.n_obs}")
    cdc2 = cdc2[keep].copy()
    total_full = cdc2.obs["total_counts"].values.astype(float)

    genes4 = [g for g in CDC2_MARKERS + ["CD1D"] + CTRL_GENES if g in cdc2.var_names]
    C = cdc2[:, genes4].layers["counts"].tocsr()
    rng = np.random.default_rng(SEED)
    det_rows, cpm_rows = [], []
    for i in range(cdc2.n_obs):
        v = C[i].toarray().ravel().astype(float)
        p = v / total_full[i]                       # 真实总深度为分母
        p_other = max(0.0, 1.0 - p.sum())
        newc = rng.multinomial(TARGET_DEPTH, np.append(p, p_other))[:-1]
        det_rows.append((newc > 0).astype(float))
        cpm_rows.append(newc / TARGET_DEPTH * 1e4)

    det = pd.DataFrame(np.array(det_rows), columns=[f"det_{g}" for g in genes4],
                       index=cdc2.obs_names)
    cpm = pd.DataFrame(np.array(cpm_rows), columns=[f"cpm_{g}" for g in genes4],
                       index=cdc2.obs_names)
    det["donor"] = cdc2.obs["donor"].values
    det["group"] = cdc2.obs["group"].values
    cpm["donor"] = cdc2.obs["donor"].values
    cpm["group"] = cdc2.obs["group"].values
    byd = det.groupby(["donor", "group"], observed=True)[
        [c for c in det.columns if c.startswith("det_")]].mean()
    byd.to_csv(f"{RES}/s5e_cdc2_downsampled_detection.csv")
    byc = cpm.groupby(["donor", "group"], observed=True)[
        [c for c in cpm.columns if c.startswith("cpm_")]].mean()
    byc.to_csv(f"{RES}/s5e_cdc2_downsampled_cpm.csv")

    ncell = cdc2.obs.groupby(["donor", "group"], observed=True).size().rename("n").reset_index()
    valid = ncell[ncell.n >= 10]
    print(f"  进入统计的供体 {len(valid)} 个 (每供体 cDC2 ≥10 细胞)")
    rows = []
    for g in genes4:
        r = donor_test(byd.reset_index().merge(valid, on=["donor", "group"]), f"det_{g}")
        if r:
            r.update(gene=g, metric="detection@5000UMI")
            rows.append(r)
        r2 = donor_test(byc.reset_index().merge(valid, on=["donor", "group"]), f"cpm_{g}")
        if r2:
            r2.update(gene=g, metric="CPM@5000UMI")
            rows.append(r2)
    dt = pd.DataFrame(rows)
    dt["fdr"] = fdr(dt.welch_p)
    dt.to_csv(f"{RES}/s5e_cdc2_depth_controlled.csv", index=False)
    print("\n  等深度检验:")
    print(dt[["gene", "metric", "failed_mean", "mature_mean", "diff",
              "cohens_d", "welch_p", "fdr"]].round(4).to_string(index=False))
    print("\n  阳性对照（管家基因；若无差异则深度已被拉平）:")
    print(dt[dt.gene.isin(CTRL_GENES)][["gene", "metric", "failed_mean",
                                        "mature_mean", "cohens_d",
                                        "welch_p"]].round(4).to_string(index=False))

    dep = cdc2.obs.groupby("group")[["total_counts", "n_genes_by_counts"]].median()
    dep.to_csv(f"{RES}/s5e_cdc2_depth_by_group.csv")
    print("\n  降采样前 cDC2 深度中位数:")
    print(dep.round(0).to_string())

    # 未降采样的原始结果（供对照，说明深度校正前后差异）
    print("\n=== 4b. 未校正深度的原始值（对照） ===")
    raw = []
    Cd = np.asarray(cdc2[:, genes4].layers["counts"].todense())
    Ld = np.asarray(cdc2[:, genes4].layers["lognorm"].todense())
    for j, g in enumerate(genes4):
        tmp = pd.DataFrame({"donor": cdc2.obs["donor"].values,
                            "group": cdc2.obs["group"].values,
                            "det": (Cd[:, j] > 0).astype(float),
                            "exp": Ld[:, j]})
        for m, col in [("detection_raw", "det"), ("lognorm_raw", "exp")]:
            bym = tmp.groupby(["donor", "group"], observed=True)[col].mean().reset_index()
            r = donor_test(bym.merge(valid, on=["donor", "group"]), col)
            if r:
                r.update(gene=g, metric=m)
                raw.append(r)
    rw = pd.DataFrame(raw)
    rw["fdr"] = fdr(rw.welch_p)
    rw.to_csv(f"{RES}/s5e_cdc2_raw_expression.csv", index=False)
    print(rw[rw.gene.isin(CDC2_MARKERS)][["gene", "metric", "failed_mean",
                                          "mature_mean", "cohens_d",
                                          "welch_p"]].round(4).to_string(index=False))

    # ---------------------------------------------------------- 5. cDC2 亚状态
    print("\n=== 5. cDC2 亚状态再聚类 ===")
    sub = a[obs["is_cDC2_cluster"].values == 1].copy()
    sub.X = sub.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=600, batch_key="donor", flavor="seurat")
    shv = sub[:, sub.var.highly_variable].copy()
    sc.pp.pca(shv, n_comps=20, zero_center=False, random_state=SEED)
    sc.external.pp.harmony_integrate(shv, "donor", random_state=SEED, max_iter_harmony=15)
    sub.obsm["X_pca"] = shv.obsm["X_pca_harmony"]
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=15, random_state=SEED)
    for r, k in [(0.4, "c2_state"), (0.8, "c2_state_hi")]:
        sc.tl.leiden(sub, resolution=r, key_added=k, random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(sub, "c2_state", method="wilcoxon", use_raw=False)
    sc.get.rank_genes_groups_df(sub, None).to_csv(f"{RES}/s5e_cdc2_state_markers.csv",
                                                  index=False)

    Cc = np.asarray(sub[:, ["CD1C"]].layers["counts"].todense()).ravel()
    sub.obs["cd1c_pos"] = (Cc > 0).astype(int)
    summ = pd.DataFrame({"n": sub.obs["c2_state"].value_counts(),
                         "cd1c_pos_frac": sub.obs.groupby("c2_state", observed=True)["cd1c_pos"].mean()})
    for g in ["CD1C", "CLEC10A", "FCER1A", "C1QA", "CD14", "HLA-DRA", "CCR7", "FSCN1"]:
        if g in sub.var_names:
            X = np.asarray(sub[:, [g]].layers["lognorm"].todense()).ravel()
            summ[g] = pd.Series(X, index=sub.obs_names).groupby(
                sub.obs["c2_state"], observed=True).mean().round(3)
    summ.to_csv(f"{RES}/s5e_cdc2_state_summary.csv")
    print(summ.round(3).to_string())

    cnts = pd.crosstab([sub.obs["donor"], sub.obs["group"]], sub.obs["c2_state"])
    cnts.to_csv(f"{RES}/s5e_cdc2_state_counts.csv")
    cnts_post = cnts.loc[cnts.index.get_level_values("group").isin(GROUPS)]
    cl2 = clr(cnts_post.astype(float))
    rows = []
    for s_ in cnts.columns:
        r = donor_test(cl2.reset_index(), s_)
        if r:
            r["state"] = s_
            rows.append(r)
    st = pd.DataFrame(rows).set_index("state")
    st["fdr"] = fdr(st.welch_p)
    st.to_csv(f"{RES}/s5e_cdc2_state_tests.csv")
    print("\n  亚状态供体水平组成:")
    print(st[["failed_mean", "mature_mean", "diff", "cohens_d",
              "welch_p", "mwu_p", "fdr"]].round(3).to_string())
    print("\n  供体 x 状态细胞数:")
    print(cnts.to_string())

    a.write_h5ad(os.path.join(WORK, "05d_myeloid.h5ad"), compression="gzip")
    sub.write_h5ad(os.path.join(WORK, "05e_cdc2.h5ad"), compression="gzip")
    print("\ndone.")


if __name__ == "__main__":
    main()
