"""S7: PAGA + 扩散拟时序 + 多命运定向吸收概率（CellRank 思路的稀疏实现）。

为什么必须定向：kNN 图是无向连通的，若只设单一吸收终态，任何细胞最终被吸收的
概率恒为 1（无信息量）。正确做法是
  1) 以扩散拟时序给每条边定向（只允许 pt 上升方向的转移）；
  2) 同时设定多个终末命运（肌成纤维 / 收缩型 SMC / 外膜驻留），
     求各命运的相对吸收概率，三者之和为 1。
统计单位仍为供体。
"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
FATES = ["Myofibroblast", "SMC_contractile", "AF_resident"]
GENES = ["POSTN", "COL8A1", "ACTA2", "DCN", "PI16", "SERPINE1", "FN1", "MYH11"]


def directed_transition(conn, pt):
    """按拟时序定向的转移矩阵；无前向邻居者置自环（吸收）。"""
    P = conn.tocsr().astype(np.float64).copy()
    rows, cols = P.nonzero()
    data = np.asarray(P[rows, cols]).ravel()
    fwd = pt[cols] > pt[rows]
    P = sp.coo_matrix((np.where(fwd, data, 0.0), (rows, cols)),
                      shape=P.shape).tocsr()
    P.eliminate_zeros()
    deg = np.asarray(P.sum(1)).ravel()
    stuck = deg == 0
    if stuck.any():
        P = P + sp.diags(stuck.astype(np.float64))
        deg = np.asarray(P.sum(1)).ravel()
    return sp.diags(1.0 / deg) @ P


def absorption(P, term_mask, fate_mask, leakage=1e-3):
    """多命运吸收概率：终末集合整体吸收，fate_mask 指定本次求解的目标命运。

    leakage：给转移矩阵掺入极小的自留概率，保证 (I - P_oo) 非奇异
    （定向图里可能存在无出路细胞，纯解析解会因子分解失败）。
    """
    n = P.shape[0]
    # 只缩放、不补自环：使每行和 = 1-eps < 1（严格次随机），
    # 保证 (I - P_oo) 非奇异；丢失的 eps 质量解释为"未抵达任何命运"。
    P = (1.0 - leakage) * P

    term = np.where(term_mask)[0]
    other = np.where(~term_mask)[0]
    P = P.tolil()
    P[term, :] = 0.0
    for i in term:
        P[i, i] = 1.0
    P = P.tocsr()

    h = np.zeros(n)
    h[term] = fate_mask[term].astype(float)
    if len(other) == 0:
        return h
    A = sp.identity(len(other), format="csr") - P[other][:, other]
    b = P[other][:, term] @ fate_mask[term].astype(float)
    try:
        h[other] = spsolve(A.tocsc(), b)
    except RuntimeError:
        from scipy.sparse.linalg import bicgstab
        h[other], _ = bicgstab(A.tocsc(), b, rtol=1e-8, maxiter=5000)
    return np.clip(h, 0, 1)


def main():
    sub = sc.read_h5ad(os.path.join(WORK, "04_fibro.h5ad"))
    print(sub.shape, sub.obs["subtype"].value_counts().to_dict())

    sc.tl.paga(sub, groups="subtype")
    sc.pl.paga(sub, plot=False)

    sc.tl.diffmap(sub, n_comps=15, random_state=SEED)
    # 根节点规则（写入 Methods）：候选 = 移植前（pre_access）的 AF_resident 细胞；
    # 取其在扩散空间中的 medoid（到候选集质心欧氏距离最小的细胞）作为根。
    # 不用 argmin(DC1)：DC1 的极值端由少数离群细胞决定，会把根放到 SMC 端，
    # 导致拟时序方向与生物学方向相反（AF_resident 反而最晚）。
    root_mask = ((sub.obs["group"] == "pre_access") &
                 (sub.obs["subtype"] == "AF_resident")).values
    if root_mask.sum() == 0:
        root_mask = (sub.obs["group"] == "pre_access").values
    dm = sub.obsm["X_diffmap"]
    cand = np.where(root_mask)[0]
    cent = dm[cand].mean(0)
    root = cand[np.argmin(((dm[cand] - cent) ** 2).sum(1))]
    sub.uns["iroot"] = int(root)
    sc.tl.dpt(sub, n_dcs=10)
    print("root:", sub.obs_names[root],
          sub.obs.iloc[root][["sample", "subtype"]].to_dict())

    pt = sub.obs["dpt_pseudotime"].values.astype(float)
    finite = np.isfinite(pt)
    if (~finite).any():
        print("non-finite pseudotime cells:", int((~finite).sum()),
              "-> 置为该命运的末端（最大 pt）")
        pt[~finite] = np.nanmax(pt[finite]) + 1e-6
    sub.obs["dpt_pseudotime"] = pt
    # DPT 被单个离群细胞拉伸（max=1 而绝大多数 <0.1），绝对尺度不宜直接解读；
    # 同时保存百分位，供图示与文字描述使用。
    sub.obs["dpt_percentile"] = pd.Series(pt).rank(pct=True).values

    diag = (pd.DataFrame({"subtype": sub.obs["subtype"].astype(str).values,
                          "pt": pt,
                          "pct": sub.obs["dpt_percentile"].values})
            .groupby("subtype").agg(n=("pt", "size"),
                                    pt_mean=("pt", "mean"),
                                    pt_median=("pt", "median"),
                                    pct_median=("pct", "median"))
            .sort_values("pct_median"))
    print("\n== pseudotime by state (sorted by median percentile) ==")
    print(diag.round(4).to_string())
    diag.to_csv(f"{RES}/s7_pseudotime_by_state.csv")

    P = directed_transition(sub.obsp["connectivities"], pt)

    # 终末细胞 = 各命运亚型中拟时序位于上四分位的细胞（而非整个亚型，
    # 否则 73% 的细胞都成了终态，矩阵近奇异且失去意义）
    term_mask = np.zeros(sub.n_obs, dtype=bool)
    fate_masks = {}
    for f in FATES:
        fm = (sub.obs["subtype"] == f).values
        if fm.sum() == 0:
            fate_masks[f] = fm
            continue
        thr = np.quantile(pt[fm], 0.75)
        t = fm & (pt >= thr)
        fate_masks[f] = t
        term_mask |= t
        print(f"  fate {f:16s} subtype n={int(fm.sum()):6d} -> terminal n={int(t.sum()):5d}")
    print("terminal cells total:", int(term_mask.sum()), "/", sub.n_obs)

    for f in FATES:
        sub.obs[f"abs_{f}"] = absorption(P, term_mask, fate_masks[f])

    cols = [f"abs_{f}" for f in FATES]
    tot = sub.obs[cols].sum(1)
    print("absorption sum (should be ~1): mean %.3f, min %.3f, max %.3f"
          % (tot.mean(), tot.min(), tot.max()))
    sub.obs[cols].to_csv(f"{RES}/s7_absorption_cells.csv")

    post = sub.obs[sub.obs["group"].isin(["mature_AVF", "failed_AVF"])]
    summ = (post.groupby(["donor", "group"], observed=True)
            .agg(n=("dpt_pseudotime", "size"),
                 dpt_mean=("dpt_pseudotime", "mean"),
                 **{f"abs_{f}": (f"abs_{f}", "mean") for f in FATES}).reset_index())
    summ.to_csv(f"{RES}/s7_donor_fate.csv", index=False)
    print(summ.round(4).to_string(index=False))

    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    rows = []
    for f in FATES:
        a = summ.loc[summ.group == "failed_AVF", f"abs_{f}"].values
        b = summ.loc[summ.group == "mature_AVF", f"abs_{f}"].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        u, pu = stats.mannwhitneyu(a, b)
        d = (a.mean() - b.mean()) / (np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) + 1e-12)
        rows.append(dict(fate=f, failed=a.mean(), mature=b.mean(),
                         diff=a.mean() - b.mean(), cohens_d=d,
                         welch_p=p, mwu_p=pu, n_failed=len(a), n_mature=len(b)))
    td = pd.DataFrame(rows)
    td["fdr"] = multipletests(td.welch_p, method="fdr_bh")[1]
    td.to_csv(f"{RES}/s7_fate_tests.csv", index=False)
    print(td.round(4).to_string(index=False))
    json.dump(td.to_dict("records"), open(f"{RES}/s7_fate_test.json", "w"), indent=2, default=float)

    tr = []
    for g in GENES:
        if g not in sub.var_names:
            continue
        x = np.asarray(sub[:, g].X.todense()).ravel()
        ok = np.isfinite(pt)
        bins = pd.qcut(pt[ok], 20, duplicates="drop")
        df = pd.DataFrame({"bin": bins, "x": x[ok]}).groupby("bin", observed=True).mean()
        df["gene"] = g
        tr.append(df.reset_index())
    pd.concat(tr).to_csv(f"{RES}/s7_gene_trend.csv", index=False)

    sub.write_h5ad(os.path.join(WORK, "05_traj.h5ad"))
    print("saved 05_traj.h5ad")


if __name__ == "__main__":
    main()
