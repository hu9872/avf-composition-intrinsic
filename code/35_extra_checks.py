"""S35: 额外稳健性数值（供正文引用，避免使用无法从结果文件复现的数字）。

内容：
 1. 剔除髓系细胞最少的两名供体后，cDC2 与髓系每千细胞比较；
 2. cDC2 内 CD1C 原始检出率（供体均值与合并口径两种）；
 3. 髓系亚型检验的对照（两类 cDC2 定义）。
输出：results/s35_extra_checks.csv
"""
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"


def cohen_d(a, b):
    sp = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1))
                 / (len(a) + len(b) - 2))
    return (np.mean(a) - np.mean(b)) / sp if sp > 0 else np.nan


def main():
    rows = []
    # ---------------------------------------------------- 1. 剔除最小供体
    a = pd.read_csv(f"{RES}/s5e_donor_abundance.csv")
    b = a[a.group.isin(["failed_AVF", "mature_AVF"])].sort_values("n_myeloid")
    drop = b.head(2).donor.tolist()
    c = b[~b.donor.isin(drop)]
    for col, label in [("cdc2_D1_per1000", "cDC2 per 1,000 (two smallest donors dropped)"),
                       ("myeloid_per1000", "Myeloid per 1,000 (two smallest donors dropped)")]:
        f = c.loc[c.group == "failed_AVF", col].values
        m = c.loc[c.group == "mature_AVF", col].values
        rows.append(dict(metric=label, failed_mean=round(float(f.mean()), 3),
                         mature_mean=round(float(m.mean()), 3),
                         cohens_d=round(float(cohen_d(f, m)), 3),
                         welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4),
                         n_failed=len(f), n_mature=len(m),
                         note="dropped " + ";".join(drop)))

    # -------------------------------------------- 2. CD1C 原始检出率
    r = pd.read_csv(f"{RES}/s5d_cdc2_detection_rate.csv")
    for grp, lab in [("failed_AVF", "failed"), ("mature_AVF", "mature_AVF")]:
        pass
    f = r.loc[r.group == "failed_AVF", "det_CD1C"].values * 100
    m = r.loc[r.group == "mature_AVF", "det_CD1C"].values * 100
    rows.append(dict(metric="CD1C detection within cDC2, raw (donor mean, %)",
                     failed_mean=round(float(f.mean()), 1),
                     mature_mean=round(float(m.mean()), 1),
                     cohens_d=round(float(cohen_d(f, m)), 3),
                     welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4),
                     n_failed=len(f), n_mature=len(m), note="before depth control"))

    adata = ad.read_h5ad(f"{ROOT}/data/processed/05e_cdc2.h5ad")
    obs = adata.obs
    pooled = obs.groupby("group")["cd1c_pos"].mean() * 100
    rows.append(dict(metric="CD1C detection within cDC2, raw (pooled cells, %)",
                     failed_mean=round(float(pooled["failed_AVF"]), 1),
                     mature_mean=round(float(pooled["mature_AVF"]), 1),
                     cohens_d=np.nan, welch_p=np.nan,
                     n_failed=int((obs.group == "failed_AVF").sum()),
                     n_mature=int((obs.group == "mature_AVF").sum()),
                     note="cell-level pooling, descriptive only"))

    # --------------------------------------- 3. 状态 3（CD1C-low C1QA+）
    t = pd.read_csv(f"{RES}/s5e_cdc2_state_tests.csv")
    s3 = t[t.state == 3].iloc[0]
    rows.append(dict(metric="CD1C-low C1QA+ state frequency (CLR)",
                     failed_mean=round(float(s3.failed_mean), 3),
                     mature_mean=round(float(s3.mature_mean), 3),
                     cohens_d=round(float(s3.cohens_d), 3),
                     welch_p=round(float(s3.welch_p), 4),
                     n_failed=int(s3.n_failed), n_mature=int(s3.n_mature),
                     note="state 3 of the cDC2 re-embedding"))

    # ------------------------- 4. 降采样后 CD1C 检出率与 CPM 的供体级检验
    for fn, pref, unit in [("s5e_cdc2_downsampled_detection", "det_", "detection (%)"),
                           ("s5e_cdc2_downsampled_cpm", "cpm_", "counts per 10,000")]:
        dd = pd.read_csv(f"{RES}/{fn}.csv")
        dd = dd[dd.group.isin(["failed_AVF", "mature_AVF"])]
        for g in ["CD1C", "CLEC10A", "FCER1A", "CD1E", "RPL13A", "ACTB"]:
            col = pref + g
            if col not in dd.columns:
                continue
            f = dd.loc[dd.group == "failed_AVF", col].values * (100 if pref == "det_" else 1)
            m = dd.loc[dd.group == "mature_AVF", col].values * (100 if pref == "det_" else 1)
            rows.append(dict(metric=f"{g} {unit}, downsampled (donor mean)",
                             failed_mean=round(float(f.mean()), 2),
                             mature_mean=round(float(m.mean()), 2),
                             cohens_d=round(float(cohen_d(f, m)), 3),
                             welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4),
                             n_failed=len(f), n_mature=len(m), note=fn))

    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/s35_extra_checks.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
