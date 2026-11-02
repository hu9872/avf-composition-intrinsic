"""S8: 配体-受体通讯（CellPhoneDB 式发送者 x 接收者打分 + 供体水平精确置换检验）。

设计要点（对应 Methods）：
1. 不以细胞为重复，统计单位为供体；
2. 打分采用 CellPhoneDB 的"表达细胞均值之积"，并要求两个 partner 在相应细胞状态
   的表达比例 >= EXPR_PROP；
3. 组间检验：Welch t + 对供体分组标签做精确置换（6 vs 6，C(12,6)=924 种划分），
   报告 BH-FDR；
4. 仅纳入至少一端为成纤维/肌成纤维细胞的互作，聚焦本研究的科学问题。
"""
import os, sys, json, itertools
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
EXPR_PROP = 0.10
GROUPS = ["mature_AVF", "failed_AVF"]

# CellPhoneDB v2 核心条目节选（覆盖 TGF-β / PDGF / ECM-integrin / 趋化 / Notch 五轴）
LR = [
    ("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"), ("TGFB2", "TGFBR1"), ("TGFB3", "TGFBR2"),
    ("PDGFB", "PDGFRB"), ("PDGFA", "PDGFRA"), ("PDGFC", "PDGFRA"), ("PDGFD", "PDGFRB"),
    ("JAG1", "NOTCH3"), ("JAG1", "NOTCH1"), ("DLL4", "NOTCH1"),
    ("CXCL12", "CXCR4"), ("CCL2", "CCR2"), ("CCL5", "CCR5"), ("CXCL8", "CXCR1"),
    ("IL6", "IL6R"), ("IL1B", "IL1R1"), ("TNF", "TNFRSF1A"),
    ("FN1", "ITGB1"), ("FN1", "ITGAV"), ("COL1A1", "ITGB1"), ("COL1A1", "DDR1"),
    ("COL1A2", "DDR1"), ("COL3A1", "DDR1"), ("POSTN", "ITGAV"), ("SPP1", "CD44"),
    ("TIMP1", "CD63"), ("EFNB2", "EPHB4"), ("ANGPT2", "TEK"), ("VEGFA", "KDR"),
    ("WNT5A", "FZD4"), ("SERPINE1", "LRP1"), ("HGF", "MET"), ("IL11", "IL11RA"),
]
FIBRO_HINT = ("Fib", "Myofibroblast", "Pericyte", "SMC")


def build_states():
    """合并主群注释与成纤维亚群注释，得到统一 state 标签。"""
    allc = sc.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    sub = sc.read_h5ad(os.path.join(WORK, "04_fibro.h5ad"))
    m = sub.obs["subtype"].astype(str)
    allc.obs["state"] = allc.obs["celltype"].astype(str)
    hit = allc.obs_names.isin(m.index)
    allc.obs.loc[hit, "state"] = "Fib_" + m.reindex(allc.obs_names[hit]).astype(str).values
    cnt = allc.obs["state"].value_counts()
    keep_states = cnt[cnt >= 30].index.tolist()
    allc = allc[allc.obs["state"].isin(keep_states)].copy()
    print("states kept:", len(keep_states), "cells:", allc.n_obs)
    print(cnt.head(30).to_string())
    return allc


def state_gene_means(allc, present):
    """返回 (均值矩阵 dict[donor]->DataFrame states x genes, 表达比例 dict)。"""
    X = allc[:, present].layers["lognorm"]
    X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float32)
    df = pd.DataFrame(X, columns=present)
    df["donor"] = allc.obs["donor"].values
    df["state"] = allc.obs["state"].values
    df["group"] = allc.obs["group"].values
    mean = df.groupby(["donor", "state"], observed=True)[present].mean()
    frac = (df[present] > 0).groupby([df["donor"], df["state"]], observed=True).mean()
    return mean, frac, df


def exact_partitions(labels, n1):
    """枚举全部 C(n, n1) 种分组划分（n=12 时 924 种）。"""
    idx = range(len(labels))
    return list(itertools.combinations(idx, n1))


def main():
    allc = build_states()
    genes = sorted({g for p in LR for g in p})
    present = [g for g in genes if g in allc.var_names]
    missing = [g for g in genes if g not in allc.var_names]
    lr_ok = [(l, r) for l, r in LR if l in present and r in present]
    print(f"LR pairs usable: {len(lr_ok)}/{len(LR)}; missing genes: {missing}")

    mean, frac, _ = state_gene_means(allc, present)

    donors = sorted(mean.index.get_level_values(0).unique())
    dgroup = {}
    for d in donors:
        g = allc.obs.loc[allc.obs.donor == d, "group"]
        dgroup[d] = g.iloc[0]
    post = [d for d in donors if dgroup[d] in GROUPS]
    print("post-op donors:", len(post), {d: dgroup[d] for d in post})

    # donor -> (states x genes) 均值/比例矩阵
    M, F = {}, {}
    for d in post:
        M[d] = mean.loc[d]
        F[d] = frac.loc[d]
    states = sorted(set().union(*[set(M[d].index) for d in post]))
    states = [s for s in states if all(s in M[d].index for d in post)]
    fib_states = [s for s in states if any(h in s for h in FIBRO_HINT)]
    print("states used:", len(states), "| fibro-related:", fib_states)

    # 每个供体在每个互作上的得分
    recs = []
    for d in post:
        Md = M[d].loc[states]
        Fd = F[d].loc[states]
        for s in states:
            for r in states:
                for l, rr in lr_ok:
                    if Fd.loc[s, l] < EXPR_PROP or Fd.loc[r, rr] < EXPR_PROP:
                        continue
                    recs.append((d, dgroup[d], s, r, l, rr,
                                 float(Md.loc[s, l] * Md.loc[r, rr])))
    sc_df = pd.DataFrame(recs, columns=["donor", "group", "sender", "receiver",
                                        "ligand", "receptor", "score"])
    print("raw interaction records:", sc_df.shape)
    sc_df.to_csv(f"{RES}/s8_lr_donor_scores.csv", index=False)

    # 只保留至少一端为成纤维/肌成纤维的互作
    focus = sc_df[sc_df.sender.isin(fib_states) | sc_df.receiver.isin(fib_states)]
    key = ["sender", "receiver", "ligand", "receptor"]
    piv = focus.pivot_table(index=key, columns="donor", values="score", aggfunc="mean")
    piv = piv.dropna(axis=0, thresh=int(0.8 * len(post)))
    print("focused interactions with sufficient donors:", piv.shape)

    nA = sum(dgroup[d] == "failed_AVF" for d in piv.columns)
    parts = exact_partitions(list(piv.columns), nA)
    vals = piv.values
    obs_diff = np.nanmean(vals[:, [i for i, c in enumerate(piv.columns)
                                   if dgroup[c] == "failed_AVF"]], axis=1) - \
               np.nanmean(vals[:, [i for i, c in enumerate(piv.columns)
                                   if dgroup[c] == "mature_AVF"]], axis=1)
    null = np.empty((len(parts), len(obs_diff)))
    for j, comb in enumerate(parts):
        a = list(comb)
        b = [i for i in range(len(piv.columns)) if i not in comb]
        null[j] = np.nanmean(vals[:, a], axis=1) - np.nanmean(vals[:, b], axis=1)
    perm_p = (np.abs(null) >= np.abs(obs_diff)).mean(0)

    rows = []
    for i, k in enumerate(piv.index):
        a = vals[i, [j for j, c in enumerate(piv.columns) if dgroup[c] == "failed_AVF"]]
        b = vals[i, [j for j, c in enumerate(piv.columns) if dgroup[c] == "mature_AVF"]]
        t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        rows.append(dict(sender=k[0], receiver=k[1], ligand=k[2], receptor=k[3],
                         failed_mean=np.nanmean(a), mature_mean=np.nanmean(b),
                         log2fc=np.log2((np.nanmean(a) + 1e-6) / (np.nanmean(b) + 1e-6)),
                         cohens_d=(np.nanmean(a) - np.nanmean(b)) /
                                  (np.sqrt((np.nanvar(a, ddof=1) + np.nanvar(b, ddof=1)) / 2) + 1e-12),
                         welch_p=p, perm_p=float(perm_p[i]),
                         n_failed=int((~np.isnan(a)).sum()), n_mature=int((~np.isnan(b)).sum())))
    td = pd.DataFrame(rows)
    td["fdr"] = multipletests(td.perm_p, method="fdr_bh")[1]
    td = td.sort_values("perm_p")
    td.to_csv(f"{RES}/s8_lr_tests.csv", index=False)
    print(td.head(15).round(4).to_string(index=False))

    # 聚焦：指向肌成纤维细胞的入向信号，供体水平汇总
    myo = [s for s in states if "Myofibroblast" in s]
    if myo:
        inc = focus[focus.receiver.isin(myo)]
        agg = inc.groupby(["donor", "group"], observed=True)["score"].sum().reset_index()
        agg.to_csv(f"{RES}/s8_incoming_to_myof.csv", index=False)
        a = agg.loc[agg.group == "failed_AVF", "score"].values
        b = agg.loc[agg.group == "mature_AVF", "score"].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        print(f"total incoming signal to myofibroblast: failed {a.mean():.4f} vs "
              f"mature {b.mean():.4f}, Welch p={p:.4g}")
        json.dump(dict(failed=float(a.mean()), mature=float(b.mean()), welch_p=float(p),
                       n_failed=int(len(a)), n_mature=int(len(b))),
                  open(f"{RES}/s8_incoming_test.json", "w"), indent=2)

    json.dump(dict(n_lr_pairs=len(lr_ok), missing_genes=missing, expr_prop=EXPR_PROP,
                   n_tested=int(len(td)), n_fdr05=int((td.fdr < 0.05).sum()),
                   n_perm_partitions=len(parts), focus_states=fib_states),
              open(f"{RES}/s8_summary.json", "w"), indent=2)
    allc.write_h5ad(os.path.join(WORK, "06_states_all.h5ad"))
    print("saved 06_states_all.h5ad | FDR<0.05 pairs:", int((td.fdr < 0.05).sum()))


if __name__ == "__main__":
    main()
