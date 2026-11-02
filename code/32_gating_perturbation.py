"""S32: cDC2 门控注释的阈值扰动敏感性 + 供体人口学汇总。

回答审稿质疑：
  Q「门控阈值是人为设定的，换阈值 cDC2 计数结论是否翻转？」
做法：把 05e 的门控阈值整体缩放 0.8 / 1.0 / 1.2 三档，重算 22 个髓系 cluster 的
标签与每供体 cDC2 计数，再按同一供体级 Welch 检验比较 failed vs matured。

同时导出供体级 sex / age（SGBA 政策要求报告性别对结果的影响）。

输入：
  data/processed/05d_myeloid.h5ad   （只读 obs，backed 模式）
  results/s5d_cluster_meanexpr.csv  （cluster x gene 均值）
  results/s5d_myeloid_total.csv     （每供体组织细胞总数）
输出：
  results/s32_gating_perturbation.csv
  results/s32_gating_cluster_labels.csv
  results/s32_donor_demographics.csv
  results/s32_summary.json
"""
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"
WORK = f"{ROOT}/data/processed"

BASE = dict(mast_tpsab1=0.5, mast_cpa3=0.5, neut_s100a8=1.5,
            cdc2_cd1c=0.4, cdc2_clec10a=0.6, cdc2_c1qa_max=1.6,
            cdc1_clec9a=0.3, pdc_lilra4=0.15, pdc_clec4c=0.08,
            mac_c1qa=1.0, mono_cd14=1.5, mono_c1qa_max=2.2,
            mono_fcn1=0.15, mono_cd14_min=1.0, mono_fcgr3a=1.5)


def gate(g, s=1.0):
    """05e_cdc2.py 的门控规则，全部阈值乘以缩放因子 s。"""
    b = {k: v * s for k, v in BASE.items()}
    if g.get("TPSAB1", 0) >= b["mast_tpsab1"] or g.get("CPA3", 0) >= b["mast_cpa3"]:
        return "Mast"
    if g.get("S100A8", 0) >= b["neut_s100a8"] and g.get("CD14", 0) < 0.5 and g.get("C1QA", 0) < 0.5:
        return "Neutrophil"
    if g.get("CD1C", 0) >= b["cdc2_cd1c"] or (g.get("CLEC10A", 0) >= b["cdc2_clec10a"]
                                             and g.get("C1QA", 0) < b["cdc2_c1qa_max"]):
        return "cDC2"
    if g.get("CLEC9A", 0) >= b["cdc1_clec9a"] or g.get("XCR1", 0) >= 0.2:
        return "cDC1"
    if g.get("LILRA4", 0) >= b["pdc_lilra4"] and g.get("CLEC4C", 0) >= b["pdc_clec4c"]:
        return "pDC"
    if g.get("C1QA", 0) >= b["mac_c1qa"]:
        if g.get("CD14", 0) >= b["mono_cd14"] and g.get("C1QA", 0) < b["mono_c1qa_max"] \
                and g.get("FCN1", 0) >= b["mono_fcn1"]:
            return "Mono_classical"
        return "Macrophage"
    if g.get("CD14", 0) >= b["mono_cd14_min"]:
        return "Mono_classical"
    if g.get("FCGR3A", 0) >= b["mono_fcgr3a"]:
        return "Mono_FCGR3A"
    return "Low_quality"


def cohen_d(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sp if sp > 0 else np.nan


def welch(a, b):
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(p)


def main():
    cm = pd.read_csv(f"{RES}/s5d_cluster_meanexpr.csv", index_col=0)
    cm.index = cm.index.astype(str)
    tot = pd.read_csv(f"{RES}/s5d_myeloid_total.csv")[["donor", "group", "n_all_cells"]]

    a = ad.read_h5ad(f"{WORK}/05d_myeloid.h5ad", backed="r")
    obs = a.obs[["donor", "group", "mye_0.9", "sex", "age"]].copy()
    obs["mye_0.9"] = obs["mye_0.9"].astype(str)
    obs["donor"] = obs["donor"].astype(str)
    a.file.close()

    # 供体人口学
    demo = obs.groupby("donor", as_index=False).agg(
        group=("group", "first"), sex=("sex", "first"), age=("age", "first"),
        n_myeloid=("donor", "size"))
    demo = demo.merge(tot[["donor", "n_all_cells"]], on="donor", how="left")
    demo.to_csv(f"{RES}/s32_donor_demographics.csv", index=False)
    print(demo.to_string(index=False))

    # 各档标签
    lab_by_s = {}
    rows = []
    for s in [0.8, 1.0, 1.2]:
        lab = {cl: gate(cm.loc[cl], s) for cl in cm.index}
        lab_by_s[s] = lab
        obs[f"lab_{s}"] = obs["mye_0.9"].map(lab).fillna("Low_quality")
        cnt = obs.groupby(["donor", f"lab_{s}"]).size().rename("n").reset_index()
        cdc2 = cnt[cnt[f"lab_{s}"] == "cDC2"].set_index("donor")["n"].reindex(
            tot.donor).fillna(0)
        df = tot.assign(n_cdc2=cdc2.values)
        df["per1000"] = df.n_cdc2 / df.n_all_cells * 1000
        df["share"] = df.n_cdc2 / df.merge(
            obs.groupby("donor").size().rename("nm").reset_index(), on="donor").nm
        f = df[df.group == "failed_AVF"].per1000.values
        m = df[df.group == "mature_AVF"].per1000.values
        rows.append(dict(scale=s, n_cdc2_total=int(df.n_cdc2.sum()),
                         n_clusters_cdc2=int(sum(1 for v in lab.values() if v == "cDC2")),
                         failed_mean_per1000=round(float(f.mean()), 3),
                         mature_mean_per1000=round(float(m.mean()), 3),
                         cohens_d=round(float(cohen_d(f, m)), 3),
                         welch_p=round(welch(f, m), 4)))
        df.to_csv(f"{RES}/s32_per1000_scale_{s}.csv", index=False)

    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/s32_gating_perturbation.csv", index=False)
    print("\n=== 门控阈值扰动 ===")
    print(out.to_string(index=False))

    labdf = pd.DataFrame(lab_by_s)
    labdf.index.name = "cluster"
    labdf.to_csv(f"{RES}/s32_gating_cluster_labels.csv")

    sex_tab = demo.pivot_table(index="group", columns="sex", values="donor",
                               aggfunc="count", fill_value=0)
    print("\n=== 性别分布（供体数）===")
    print(sex_tab.to_string())

    summary = dict(
        scales={str(r.scale): dict(n_cdc2=int(r.n_cdc2_total), d=r.cohens_d,
                                   p=r.welch_p)
                for r in out.itertuples()},
        conclusion_stable=bool(all(r.welch_p > 0.05 for r in out.itertuples())),
        sex_by_group=sex_tab.to_dict(),
        age_range=[float(demo.age.min()), float(demo.age.max())],
        n_donors=int(demo.donor.nunique()),
    )
    with open(f"{RES}/s32_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, default=str)
    print("\n" + json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
