"""S34: 独立 bulk 队列的区室签名评分（外部方向一致性检验）。

背景：NNLS 反卷积（33_deconvolution.py）在同一组织参考谱下不可辨识
（拟合 R2 中位 0.31，权重塌缩到单一区室），因此改用更稳健的签名评分：
  1) 由本研究图谱的区室 CPM 参考谱（s33_reference_cpm.csv）推导每区室 top-50 特异基因；
  2) 在 GSE220796 的 AVF bulk 样本上，对每个签名取基因 z 分数均值；
  3) 比较 failed (11) vs matured (9)。
签名评分检验的是「该区室特征基因的整体丰度是否更高」，不是细胞比例，
因此只能回答方向是否一致，不能替代比例估计。

输出：
  results/s34_compartment_signatures_genelist.csv
  results/s34_compartment_scores.csv
  results/s34_compartment_score_tests.csv
  figures/FigS3_compartment_signatures.pdf/.png
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"
RAW = f"{ROOT}/data/raw"
FIG = f"{ROOT}/figures"
TOP_N = 50


def top_markers(ref, n=TOP_N):
    lg = np.log2(ref + 1)
    genes = {}
    for c in ref.columns:
        others = [x for x in ref.columns if x != c]
        spec = lg[c] - lg[others].mean(axis=1)
        ok = (ref[c] > 10) & (spec > 0.5)
        genes[c] = spec[ok].sort_values(ascending=False).head(n).index.tolist()
    return genes


def main():
    ref = pd.read_csv(f"{RES}/s33_reference_cpm.csv", index_col=0)
    markers = top_markers(ref)
    pd.DataFrame([dict(compartment=c, rank=i + 1, gene=g)
                  for c, gl in markers.items() for i, g in enumerate(gl)]
                 ).to_csv(f"{RES}/s34_compartment_signatures_genelist.csv", index=False)

    meta = pd.read_csv(f"{RES}/s10_meta_GSE220796.csv")
    cnt = pd.read_csv(f"{RAW}/GSE220796_rVazquez-Padron_counts.txt.gz", sep="\t", index_col=0)
    cnt.index = [str(x).strip('"') for x in cnt.index]
    cnt = cnt.groupby(level=0).sum()
    meta = meta[meta.code.isin(cnt.columns)]
    n_tissue = meta.groupby("patient").tissue.nunique()
    paired = n_tissue[n_tissue == 2].index
    avf = meta[meta.tissue.str.contains("Juxta") & meta.patient.isin(paired)].set_index("code")

    lg = np.log2(cnt / cnt.sum(axis=0) * 1e6 + 1)
    lg = lg[avf.index]

    rows = []
    for c, gl in markers.items():
        g = [x for x in gl if x in lg.index]
        if len(g) < 10:
            continue
        z = lg.loc[g].apply(lambda r: (r - r.mean()) / (r.std() + 1e-9), axis=1)
        s = z.mean(axis=0)
        for smp, v in s.items():
            rows.append(dict(sample=smp, compartment=c, score=float(v),
                             outcome=avf.loc[smp, "outcome"], n_genes=len(g)))
    S = pd.DataFrame(rows)
    S.to_csv(f"{RES}/s34_compartment_scores.csv", index=False)

    out = []
    for c, d in S.groupby("compartment"):
        f = d.loc[d.outcome == "Failed", "score"].values
        m = d.loc[d.outcome == "Matured", "score"].values
        sp = np.sqrt(((len(f) - 1) * np.var(f, ddof=1) + (len(m) - 1) * np.var(m, ddof=1))
                     / (len(f) + len(m) - 2))
        out.append(dict(compartment=c, n_genes=int(d.n_genes.iloc[0]),
                        failed_mean=round(float(f.mean()), 3),
                        mature_mean=round(float(m.mean()), 3),
                        diff=round(float(f.mean() - m.mean()), 3),
                        cohens_d=round(float((f.mean() - m.mean()) / sp), 3),
                        welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4),
                        n_failed=len(f), n_matured=len(m)))
    T = pd.DataFrame(out).sort_values("welch_p")
    T.to_csv(f"{RES}/s34_compartment_score_tests.csv", index=False)
    print(T.to_string(index=False))

    # 管家基因标化变体：除以 RPL13A/RPLP0/B2M/TBP/GAPDH/ACTB/PPIA 均值的 log2 CPM
    hk = [g for g in ["RPL13A", "RPLP0", "B2M", "TBP", "GAPDH", "ACTB", "PPIA"]
          if g in lg.index]
    hk_mean = lg.loc[hk].mean(axis=0)
    hrows = []
    for c, gl in markers.items():
        g = [x for x in gl if x in lg.index]
        if len(g) < 10:
            continue
        s_norm = lg.loc[g].mean(axis=0) - hk_mean
        f = s_norm[avf.outcome == "Failed"].values
        m = s_norm[avf.outcome == "Matured"].values
        sp = np.sqrt(((len(f) - 1) * np.var(f, ddof=1) + (len(m) - 1) * np.var(m, ddof=1))
                     / (len(f) + len(m) - 2))
        hrows.append(dict(compartment=c, failed_mean=round(float(f.mean()), 3),
                          mature_mean=round(float(m.mean()), 3),
                          diff=round(float(f.mean() - m.mean()), 3),
                          cohens_d=round(float((f.mean() - m.mean()) / sp), 3),
                          welch_p=round(float(stats.ttest_ind(f, m, equal_var=False).pvalue), 4)))
    pd.DataFrame(hrows).sort_values("welch_p").to_csv(
        f"{RES}/s34_housekeeping_normalised.csv", index=False)
    print("\n=== 管家基因标化 ===")
    print(pd.DataFrame(hrows).sort_values("welch_p").to_string(index=False))

    order = ["Myeloid", "Fibroblast", "Endothelial", "SMC", "T_cell", "NK", "B_cell", "Mast"]
    order = [c for c in order if c in set(S.compartment)]
    fig, ax = plt.subplots(figsize=(5.6, 3.2), dpi=300)
    for i, c in enumerate(order):
        for j, (oc, col) in enumerate([("Matured", "#2E86C1"), ("Failed", "#C0392B")]):
            v = S[(S.compartment == c) & (S.outcome == oc)].score.values
            ax.scatter(np.full(len(v), i + (j - 0.5) * 0.28), v, s=15, color=col,
                       alpha=0.85, linewidth=0, label=oc if i == 0 else None)
        ax.plot([i - 0.24, i + 0.24],
                [S[(S.compartment == c) & (S.outcome == "Matured")].score.mean(),
                 S[(S.compartment == c) & (S.outcome == "Failed")].score.mean()],
                color="0.2", lw=1.2, zorder=4)
    ax.axhline(0, color="0.7", lw=0.7, zorder=0)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([c.replace("_", " ") for c in order], fontsize=7, rotation=30,
                       ha="right")
    ax.set_ylabel("Compartment signature score (mean z)", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=2)
    ax.set_title("GSE220796 bulk, top-50 signatures (11 failed, 9 matured)",
                 fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(f"{FIG}/FigS3_compartment_signatures.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/FigS3_compartment_signatures.png", bbox_inches="tight", dpi=300)
    print("saved FigS3_compartment_signatures")


if __name__ == "__main__":
    main()
