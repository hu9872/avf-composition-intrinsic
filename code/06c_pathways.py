"""S8c: 配体-受体信号的通路级聚合（解决 744 个单互作检验的多重校正瓶颈）。

动机：精确置换（6 vs 6，C(12,6)=924 种划分）的最小可达 p 为 0.0011，
在 744 个单互作检验上做 BH 校正不可能达到 FDR<0.05。这不是"无差异"的证据，
而是检验单元的粒度问题。因此把互作归并为 9 条信号轴，在供体水平检验轴均值。
"""
import os, itertools, json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

RES = r"D:/AVF_GSE_2611/results"
AXES = {
    "TGF-beta": ["TGFB1", "TGFB2", "TGFB3"],
    "PDGF": ["PDGFB", "PDGFA", "PDGFC", "PDGFD"],
    "Fibronectin-integrin": ["FN1"],
    "Collagen-integrin/DDR": ["COL1A1", "COL1A2", "COL3A1"],
    "Periostin-integrin": ["POSTN"],
    "Chemokine": ["CXCL12", "CCL2", "CCL5", "CXCL8"],
    "Inflammatory": ["IL6", "IL1B", "TNF", "IL11"],
    "Notch": ["JAG1", "DLL4"],
    "Angiogenic": ["ANGPT2", "VEGFA", "EFNB2", "HGF", "WNT5A", "SERPINE1"],
}


def main():
    d = pd.read_csv(f"{RES}/s8_lr_donor_scores.csv")
    d["axis"] = "Other"
    for ax, lig in AXES.items():
        d.loc[d.ligand.isin(lig), "axis"] = ax
    d = d[d.axis != "Other"].copy()

    d["group"] = np.where(d.donor.str.startswith("F"), "failed_AVF", "mature_AVF")
    donors = sorted(d.donor.unique())
    dg = d.groupby("donor")["group"].first()

    rows = []
    for ax, g in d.groupby("axis"):
        piv = g.pivot_table(index=["sender", "receiver", "ligand", "receptor"],
                            columns="donor", values="score", aggfunc="mean")
        # 供体水平轴均值 = 该轴内所有互作得分的均值
        per_donor = piv.mean(0)
        a = per_donor[[x for x in per_donor.index if dg[x] == "failed_AVF"]].values
        b = per_donor[[x for x in per_donor.index if dg[x] == "mature_AVF"]].values
        if len(a) < 3 or len(b) < 3:
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        u, pu = stats.mannwhitneyu(a, b)
        dv = (a.mean() - b.mean()) / (np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) + 1e-12)
        rows.append(dict(axis=ax, failed=a.mean(), mature=b.mean(),
                         log2fc=np.log2((a.mean() + 1e-6) / (b.mean() + 1e-6)),
                         diff=a.mean() - b.mean(), cohens_d=dv,
                         welch_p=p, mwu_p=pu,
                         n_failed=len(a), n_mature=len(b), n_interactions=piv.shape[0]))
    td = pd.DataFrame(rows).sort_values("welch_p")
    td["fdr"] = multipletests(td.welch_p, method="fdr_bh")[1]
    td.round(4).to_csv(f"{RES}/s8c_pathway_tests.csv", index=False)
    print(td.round(4).to_string(index=False))

    # 反向补充：仅"指向成纤维/肌成纤维"的入向信号
    fib = d.receiver.str.contains("Fib_")
    rows2 = []
    for ax in AXES:
        g = d[fib & (d.axis == ax)]
        if g.empty:
            continue
        per = g.groupby(["donor"])["score"].mean()
        a = per[[x for x in per.index if dg[x] == "failed_AVF"]].values
        b = per[[x for x in per.index if dg[x] == "mature_AVF"]].values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows2.append(dict(axis=ax, failed=a.mean(), mature=b.mean(),
                          log2fc=np.log2((a.mean() + 1e-6) / (b.mean() + 1e-6)),
                          welch_p=p, n_failed=len(a), n_mature=len(b)))
    t2 = pd.DataFrame(rows2).sort_values("welch_p")
    t2["fdr"] = multipletests(t2.welch_p, method="fdr_bh")[1]
    t2.round(4).to_csv(f"{RES}/s8c_pathway_tests_incoming_fib.csv", index=False)
    print("\n--- 仅指向成纤维谱系的入向信号 ---")
    print(t2.round(4).to_string(index=False))

    json.dump(dict(n_axes=len(td), n_sig=int((td.fdr < 0.05).sum()),
                   min_perm_resolution=1 / 924),
              open(f"{RES}/s8c_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
