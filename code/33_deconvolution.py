"""S33: 独立 bulk 队列的组成反卷积（GSE220796）。

目的：为单细胞层面「失败组髓系扩增」提供外部方向一致性检验。
公开数据中没有人 AVF 的第二套 scRNA 队列，因此改在同一组织的独立 bulk 队列上
用非负最小二乘（NNLS）估计区室比例。

参考谱：本研究图谱（03_annotated.h5ad，66,168 细胞，8 区室）的区室平均 counts，
       按区室总计数归一为 CPM 后作为签名矩阵。
目标：GSE220796 的 juxta-anastomotic AVF 样本（11 failed / 9 matured），CPM。

诚实边界：参考谱与本队列同组织但不同患者、不同平台（10x 3' vs bulk polyA），
NNLS 假设签名完整、无批次效应。结果只作方向一致性参考，不作定量结论；
若方向不一致，如实报告为「未获外部支持」。

输出：
  results/s33_reference_cpm.csv
  results/s33_deconvolution.csv
  results/s33_deconvolution_tests.csv
  figures/FigS3_deconvolution.pdf / .png
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"
RAW = f"{ROOT}/data/raw"
FIG = f"{ROOT}/figures"
WORK = f"{ROOT}/data/processed"


def build_reference():
    cache = f"{RES}/s33_reference_cpm.csv"
    if os.path.exists(cache):
        return pd.read_csv(cache, index_col=0)
    a = ad.read_h5ad(f"{WORK}/03_annotated.h5ad", backed="r")
    ct = a.obs["celltype"].astype(str).values
    X = a.layers["counts"]
    prof = []
    for c in sorted(set(ct)):
        m = np.where(ct == c)[0]
        sub = X[m]
        s = np.asarray(sub.sum(axis=0)).ravel().astype(float)
        prof.append(pd.Series(s / s.sum() * 1e6, index=a.var_names, name=c))
    a.file.close()
    ref = pd.concat(prof, axis=1)
    ref.to_csv(cache)
    return ref


def main():
    meta = pd.read_csv(f"{RES}/s10_meta_GSE220796.csv")
    cnt = pd.read_csv(f"{RAW}/GSE220796_rVazquez-Padron_counts.txt.gz", sep="\t", index_col=0)
    cnt.index = [str(x).strip('"') for x in cnt.index]
    cnt = cnt.groupby(level=0).sum()
    meta = meta[meta.code.isin(cnt.columns)]
    n_tissue = meta.groupby("patient").tissue.nunique()
    paired = n_tissue[n_tissue == 2].index
    avf = meta[meta.tissue.str.contains("Juxta") & meta.patient.isin(paired)].set_index("code")

    ref = build_reference()
    genes = [g for g in ref.index if g in cnt.index]
    print("shared genes:", len(genes), "| compartments:", list(ref.columns))

    y_all = cnt.loc[genes, avf.index]
    bulk_cpm = y_all / y_all.sum(axis=0) * 1e6
    keep = (bulk_cpm.median(axis=1) > 5.0) & (ref.loc[genes].max(axis=1) > 1.0)
    X = ref.loc[genes][keep.values].values
    Y = bulk_cpm.loc[keep].values
    print("genes used:", int(keep.sum()))

    rows, r2 = [], []
    for j, smp in enumerate(bulk_cpm.columns):
        y = Y[:, j]
        w, _ = nnls(X, y)
        tot = w.sum()
        rows.append(dict(sample=smp, **{c: v / tot for c, v in zip(ref.columns, w)}))
        fit = X @ w
        r2.append(1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    P = pd.DataFrame(rows).set_index("sample")
    P["fit_r2"] = r2
    P["outcome"] = avf.loc[P.index, "outcome"].values
    P["patient"] = avf.loc[P.index, "patient"].values
    P.to_csv(f"{RES}/s33_deconvolution.csv")
    print("\nmedian fit R2:", round(float(np.median(r2)), 3))
    print(P.round(3).to_string())

    out = []
    for c in ref.columns:
        f = P.loc[P.outcome == "Failed", c].values
        m = P.loc[P.outcome == "Matured", c].values
        sp = np.sqrt(((len(f) - 1) * np.var(f, ddof=1) + (len(m) - 1) * np.var(m, ddof=1))
                     / (len(f) + len(m) - 2))
        out.append(dict(compartment=c,
                        failed_mean=round(float(f.mean()), 4),
                        mature_mean=round(float(m.mean()), 4),
                        diff=round(float(f.mean() - m.mean()), 4),
                        ratio=round(float(f.mean() / m.mean()), 3) if m.mean() > 0 else np.nan,
                        cohens_d=round(float((f.mean() - m.mean()) / sp), 3) if sp > 0 else np.nan,
                        welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4),
                        n_failed=len(f), n_matured=len(m)))
    T = pd.DataFrame(out).sort_values("welch_p")
    T.to_csv(f"{RES}/s33_deconvolution_tests.csv", index=False)
    print("\n=== 反卷积比例比较 ===")
    print(T.to_string(index=False))

    comps = ["Myeloid", "Fibroblast", "Endothelial", "T_cell"]
    fig, ax = plt.subplots(figsize=(4.6, 3.2), dpi=300)
    for i, c in enumerate(comps):
        for j, (oc, col) in enumerate([("Matured", "#2E86C1"), ("Failed", "#C0392B")]):
            v = P.loc[P.outcome == oc, c].values
            ax.scatter(np.full(len(v), i + (j - 0.5) * 0.28), v * 100, s=15,
                       color=col, alpha=0.85, linewidth=0,
                       label=oc if i == 0 else None)
        ax.plot([i - 0.24, i + 0.24],
                [P.loc[P.outcome == "Matured", c].mean() * 100,
                 P.loc[P.outcome == "Failed", c].mean() * 100],
                color="0.2", lw=1.2, zorder=4)
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels([c.replace("_", " ") for c in comps], fontsize=7.5)
    ax.set_ylabel("Deconvolved proportion (%)", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_title("GSE220796 bulk, NNLS (11 failed, 9 matured)", fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(f"{FIG}/FigS3_deconvolution.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/FigS3_deconvolution.png", bbox_inches="tight", dpi=300)
    print("saved FigS3")


if __name__ == "__main__":
    main()
