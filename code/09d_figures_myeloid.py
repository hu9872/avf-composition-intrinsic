"""S9d: 髓系/cDC2 模块出图（Fig 7 与 Fig 8）。"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import anndata as ad

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
FIG = r"D:/AVF_GSE_2611/figures"
os.makedirs(FIG, exist_ok=True)

FAIL, MAT = "#C44E52", "#55A868"
GREY = "#9E9E9E"
mpl.rcParams.update({"font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "font.family": "DejaVu Sans", "pdf.fonttype": 42})

CL_COLORS = {
    "Macrophage": "#B7C9E2",
    "Mono_classical": "#8FB3D9",
    "cDC2": "#E08A3C",
    "pDC": "#7A5195",
    "Neutrophil": "#BBBBBB",
    "Low_quality": "#DDDDDD",
}


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  saved", name)


def fig7():
    a = ad.read_h5ad(os.path.join(WORK, "05d_myeloid.h5ad"))
    fig = plt.figure(figsize=(12.4, 3.5))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[1.15, 1.25, 1, 1],
                  wspace=0.48, left=0.06, right=0.985, top=0.79, bottom=0.19)

    # A. UMAP
    ax = fig.add_subplot(gs[0, 0])
    U = a.obsm["X_umap"]
    order = ["Low_quality", "Neutrophil", "Macrophage", "Mono_classical", "pDC", "cDC2"]
    for s in order:
        m = (a.obs["mye_subtype"] == s).values
        if m.sum() == 0:
            continue
        ax.scatter(U[m, 0], U[m, 1], s=1.1, c=CL_COLORS[s], linewidths=0,
                   rasterized=True, label=f"{s} ({int(m.sum()):,})")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("UMAP1", fontsize=7); ax.set_ylabel("UMAP2", fontsize=7)
    ax.set_title("Myeloid compartment\n(14,845 cells)", fontsize=8, loc="left", pad=6)
    ax.legend(markerscale=5, fontsize=5.6, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.10), ncol=2, handletextpad=0.2,
              columnspacing=0.8)

    # B. marker heatmap（z-score across subtypes）+ 未检出标记
    ax = fig.add_subplot(gs[0, 1])
    cm = pd.read_csv(f"{RES}/s5d_cluster_meanexpr.csv", index_col=0)
    sub_lab = pd.read_csv(f"{RES}/s5e_subtype_marker_check.csv", index_col=0)
    subs = [s for s in ["Macrophage", "Mono_classical", "cDC2", "pDC"]
            if s in a.obs["mye_subtype"].values]
    genes = ["CD14", "FCN1", "S100A8", "FCGR3A", "C1QA", "CD163", "MRC1",
             "CD1C", "CLEC10A", "FCER1A", "CD1E", "LILRA4", "CLEC4C",
             "CLEC9A", "XCR1"]
    genes = [g for g in genes if g in sub_lab.columns]
    M = sub_lab.loc[subs, genes].T.astype(float)
    im = ax.imshow(M.values, aspect="auto", cmap="Blues", vmin=0, vmax=2.2)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs, rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            ax.text(j, i, "0" if v < 0.005 else f"{v:.2f}", ha="center", va="center",
                    fontsize=4.4, color="white" if v > 1.40 else "#222222")
    # 未检出标记：绝对表达 <0.10 视为未检出
    nd = [i for i, g in enumerate(genes) if M.loc[g].max() < 0.10]
    nd = [i for i, g in enumerate(genes) if M.loc[g].max() < 0.10]
    if nd:
        ax.axhline(min(nd) - 0.5, color="#B22222", lw=0.9)
        for i in nd:
            ax.get_yticklabels()[i].set_color("#B22222")
    ax.set_title("Subtype markers, mean log-norm.\n(red = not detected)", fontsize=8,
                 loc="left", pad=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.ax.tick_params(labelsize=5.5); cb.outline.set_linewidth(0.4)

    # C. cDC2 绝对丰度（每千细胞）
    comp = pd.read_csv(f"{RES}/s5e_donor_abundance.csv")
    post = comp[comp.group.isin(["failed_AVF", "mature_AVF"])].copy()
    ax = fig.add_subplot(gs[0, 2])
    series = [("failed_AVF", FAIL, "Failed"), ("mature_AVF", MAT, "Matured")]
    for i, (g, c, lab) in enumerate(series):
        v = post.loc[post.group == g, "cdc2_D1_per1000"].values
        x = np.random.default_rng(1).normal(i, 0.045, len(v))
        ax.scatter(x, v, s=16, c=c, edgecolor="white", linewidth=0.4, zorder=3)
        ax.hlines(np.mean(v), i - 0.22, i + 0.22, color=c, lw=1.6, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([s[2] for s in series], fontsize=7)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("cDC2 per 1,000 cells", fontsize=7)
    ax.set_title("cDC2 abundance unchanged\n(d = 0.12, p = 0.84)", fontsize=8, loc="left")
    ax.tick_params(labelsize=6.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    # D. cDC2 占髓系比例
    ax = fig.add_subplot(gs[0, 3])
    for i, (g, c, lab) in enumerate(series):
        v = post.loc[post.group == g, "cdc2_D1_frac_myeloid"].values * 100
        x = np.random.default_rng(2).normal(i, 0.045, len(v))
        ax.scatter(x, v, s=16, c=c, edgecolor="white", linewidth=0.4, zorder=3)
        ax.hlines(np.mean(v), i - 0.22, i + 0.22, color=c, lw=1.6, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([s[2] for s in series], fontsize=7)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("cDC2 (% of myeloid)", fontsize=7)
    ax.set_title("Diluted by myeloid expansion\n(d = −1.06, p = 0.12)", fontsize=8, loc="left")
    ax.tick_params(labelsize=6.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    fig.suptitle("Figure 7  The cDC2 compartment is numerically preserved in failed AVF",
                 fontsize=9, x=0.06, ha="left", y=0.97)
    save(fig, "Fig7_myeloid_cdc2")


def fig8():
    dec = pd.read_csv(f"{RES}/s5f_decomposition_retained.csv")
    dec = dec.reindex(dec.log2FC_observed.abs().sort_values().index)
    det = pd.read_csv(f"{RES}/s5e_cdc2_downsampled_detection.csv")

    fig = plt.figure(figsize=(11.8, 3.6))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.5, 1, 1],
                  wspace=0.42, left=0.06, right=0.985, top=0.80, bottom=0.20)

    # A. 分解：观测 vs 构成 vs 内部
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(dec))
    ax.barh(y + 0.18, dec.log2FC_observed, height=0.34, color="#444444",
            label="Observed (pseudobulk)", zorder=3)
    ax.barh(y - 0.18, dec.log2FC_predicted, height=0.34, color="#C9A227",
            label="Predicted (composition + intrinsic)", zorder=3, alpha=0.9)
    ax.scatter(dec.log2FC_composition, y + 0.18, marker="|", s=90, c=FAIL, lw=1.6,
               zorder=4, label="Composition component")
    ax.set_yticks(y); ax.set_yticklabels(dec.gene, fontsize=6.8)
    ax.axvline(0, c="k", lw=0.6)
    ax.set_xlabel("log2 fold-change (failed / matured)", fontsize=7)
    ax.set_title("Compartment pseudobulk decomposes into\ncomposition × intrinsic expression",
                 fontsize=8, loc="left")
    ax.legend(fontsize=5.6, frameon=False, loc="lower left")
    ax.tick_params(labelsize=6.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    # B. 等深度 CD1C 检出率
    ax = fig.add_subplot(gs[0, 1])
    d = det[det.group.isin(["failed_AVF", "mature_AVF"])]
    series = [("failed_AVF", FAIL, "Failed"), ("mature_AVF", MAT, "Matured")]
    for i, (g, c, lab) in enumerate(series):
        v = d.loc[d.group == g, "det_CD1C"].values * 100
        x = np.random.default_rng(3).normal(i, 0.045, len(v))
        ax.scatter(x, v, s=18, c=c, edgecolor="white", linewidth=0.4, zorder=3)
        ax.hlines(np.mean(v), i - 0.22, i + 0.22, color=c, lw=1.6, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([s[2] for s in series], fontsize=7)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("CD1C+ cells within cDC2 (%)", fontsize=7)
    ax.set_title("Within cDC2, at equal depth\n(d = −0.84, p = 0.23)", fontsize=8, loc="left")
    ax.tick_params(labelsize=6.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    # C. 单个 cDC2 细胞测序深度
    ax = fig.add_subplot(gs[0, 2])
    dep = pd.read_csv(f"{RES}/s5e_cdc2_depth_by_group.csv")
    dep = dep[dep.group.isin(["failed_AVF", "mature_AVF"])]
    cols = [FAIL if g == "failed_AVF" else MAT for g in dep.group]
    ax.bar(["Failed", "Matured"], dep.set_index("group").loc[
        ["failed_AVF", "mature_AVF"], "total_counts"].values, color=cols, width=0.55)
    for i, g in enumerate(["failed_AVF", "mature_AVF"]):
        v = dep.set_index("group").loc[g, "total_counts"]
        ax.text(i, v, f"{v/1000:.1f}k", ha="center", va="bottom", fontsize=6.5)
    ax.set_ylabel("Median UMIs per cDC2 cell", fontsize=7)
    ax.set_title("Depth confounder corrected\nby downsampling to 5,000 UMIs",
                 fontsize=8, loc="left")
    ax.tick_params(labelsize=6.5)
    ax.set_ylim(0, dep.total_counts.max() * 1.22)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    fig.suptitle("Figure 8  The myeloid pseudobulk signal reflects composition and depth, "
                 "not cDC2 loss", fontsize=9, x=0.06, ha="left", y=0.97)
    save(fig, "Fig8_signal_decomposition")


if __name__ == "__main__":
    fig7()
    fig8()
