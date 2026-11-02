"""终版图表：Fig1-Fig8 + FigS1。
设计约束：面板字母、无文字重叠、统一配色与字号、PDF 矢量 + 300 dpi PNG。
所有数值直接读 results/ 下的真实结果表。
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

RES = r"D:/AVF_GSE_2611/results"
FIG = r"D:/AVF_GSE_2611/figures"
os.makedirs(FIG, exist_ok=True)


def d_p(a, b):
    """Cohen's d (pooled SD) and Welch p for failed vs matured, computed from the
    values actually drawn in the panel so that annotation and points cannot drift."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    return f"d = {(a.mean() - b.mean()) / sp:.2f}, p = {stats.ttest_ind(a, b, equal_var=False).pvalue:.2f}"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# ---------------------------------------------------------------- palette
GROUP_COL = {"pre_access": "#8C8C8C", "early_AVF": "#E0A33E",
             "mature_AVF": "#2E7D5B", "failed_AVF": "#B03A2E"}
GROUP_LAB = {"pre_access": "Pre-access vein", "early_AVF": "Early AVF",
             "mature_AVF": "Matured AVF", "failed_AVF": "Failed AVF"}
GROUP_ORDER = ["pre_access", "early_AVF", "mature_AVF", "failed_AVF"]

CT_COL = {"Fibroblast": "#4C9F70", "Myeloid": "#8172B2", "Endothelial": "#E07B39",
          "T_cell": "#5B8DB8", "SMC": "#C44E52", "NK": "#7F7F7F",
          "Mast": "#D4A6C8", "B_cell": "#2F6C8F"}
CT_LAB = {"Fibroblast": "Fibroblast", "Myeloid": "Myeloid", "Endothelial": "Endothelial",
          "T_cell": "T cell", "SMC": "SMC", "NK": "NK cell", "Mast": "Mast cell",
          "B_cell": "B cell"}

ST_COL = {"AF_resident": "#2E7D5B", "AF_stromal": "#7FB77E", "Myofibroblast": "#C0392B",
          "Pericyte": "#E08A3C", "SMC_contractile": "#4C72B0", "SMC_synthetic": "#8172B2"}
ST_LAB = {"AF_resident": "Resident AF", "AF_stromal": "Stromal AF",
          "Myofibroblast": "Myofibroblast", "Pericyte": "Pericyte",
          "SMC_contractile": "Contractile SMC", "SMC_synthetic": "Synthetic SMC"}
ST_ORDER = ["AF_resident", "AF_stromal", "Myofibroblast", "Pericyte",
            "SMC_contractile", "SMC_synthetic"]

MY_COL = {"Macrophage": "#8172B2", "Mono_classical": "#4C72B0", "cDC2": "#C0392B",
          "pDC": "#E0A33E", "Neutrophil": "#55A868", "Low_quality": "#BFBFBF"}
MY_LAB = {"Macrophage": "Macrophage", "Mono_classical": "Classical monocyte",
          "cDC2": "cDC2", "pDC": "pDC", "Neutrophil": "Neutrophil",
          "Low_quality": "Low quality"}
MY_ORDER = ["Macrophage", "Mono_classical", "cDC2", "pDC", "Neutrophil", "Low_quality"]

RED, GREEN, GREY = "#B03A2E", "#2E7D5B", "#9E9E9E"


def panel(ax, letter, dx=-0.055, dy=1.045):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left")


def umap_frame(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)


def repel_labels(ax, xs, ys, labels, colors, fontsize=6.5):
    """在显示坐标中做各向异性排斥，避免标签互相遮挡。"""
    fig = ax.figure
    fig.canvas.draw()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    w = np.array([max(26, 3.4 * fontsize * len(l)) for l in labels])
    h = np.full(len(labels), 3.8 * fontsize)
    for _ in range(400):
        moved = False
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                dx = pts[j, 0] - pts[i, 0]
                dy = pts[j, 1] - pts[i, 1]
                ox = (w[i] + w[j]) / 2 + 6
                oy = (h[i] + h[j]) / 2 + 4
                if abs(dx) < ox and abs(dy) < oy:
                    sx = (ox - abs(dx)) / 2 * np.sign(dx if dx else 1)
                    sy = (oy - abs(dy)) / 2 * np.sign(dy if dy else 1)
                    pts[i] -= [sx * 0.6, sy * 0.6]
                    pts[j] += [sx * 0.6, sy * 0.6]
                    moved = True
        if not moved:
            break
    inv = ax.transData.inverted()
    for (x, y), lab, c in zip(inv.transform(pts), labels, colors):
        ax.annotate(lab, xy=(x, y), fontsize=fontsize, color="#1A1A1A", fontweight="bold",
                    ha="center", va="center", zorder=20,
                    bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=c, lw=0.5, alpha=0.88))


def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf")
    fig.savefig(f"{FIG}/{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", name)


# ================================================================ Fig 1
atlas = pd.read_csv(f"{RES}/plot_atlas.csv")
fig = plt.figure(figsize=(13.2, 4.0))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.12], wspace=0.30,
              left=0.03, right=0.99, top=0.90, bottom=0.14)

ax = fig.add_subplot(gs[0, 0])
order = ["Fibroblast", "Myeloid", "Endothelial", "SMC", "NK", "Mast", "T_cell", "B_cell"]
d = atlas[atlas.celltype.isin(order)]
for ct in order:
    s = d[d.celltype == ct]
    ax.scatter(s.U1, s.U2, s=0.5, c=CT_COL[ct], linewidths=0, rasterized=True, alpha=0.85)
cent = d.groupby("celltype")[["U1", "U2"]].median()
repel_labels(ax, cent.loc[order, "U1"].values, cent.loc[order, "U2"].values,
             [CT_LAB[c] for c in order], [CT_COL[c] for c in order])
umap_frame(ax)
ax.set_title("66,168 cells, 8 compartments", pad=4)
panel(ax, "A")

ax = fig.add_subplot(gs[0, 1])
for g in GROUP_ORDER:
    s = atlas[atlas.group == g]
    ax.scatter(s.U1, s.U2, s=0.35, c=GROUP_COL[g], linewidths=0, rasterized=True, alpha=0.75)
umap_frame(ax)
ax.set_title("Sampling stage", pad=4)
handles = [Line2D([], [], marker="o", ls="", ms=3.5, color=GROUP_COL[g],
                  label=GROUP_LAB[g]) for g in GROUP_ORDER]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=4,
          frameon=False, handletextpad=0.2, columnspacing=0.9)
panel(ax, "B")

ax = fig.add_subplot(gs[0, 2])
tab = pd.crosstab(atlas.group, atlas.celltype)[order].reindex(GROUP_ORDER)
frac = tab.div(tab.sum(1), axis=0)
bottom = np.zeros(len(frac))
xpos = np.arange(len(frac))
for ct in order:
    v = frac[ct].values
    ax.bar(xpos, v, bottom=bottom, width=0.66, color=CT_COL[ct], edgecolor="white", linewidth=0.4)
    bottom += v
ax.set_xticks(xpos)
ax.set_xticklabels([GROUP_LAB[g].replace(" ", "\n") for g in GROUP_ORDER], fontsize=6.5)
ax.set_ylabel("Fraction of cells")
ax.set_ylim(0, 1.0)
ax.set_title("Compartment composition by stage", pad=4)
ax.spines[["top", "right"]].set_visible(False)
handles = [Patch(facecolor=CT_COL[c], edgecolor="none", label=CT_LAB[c]) for c in order]
ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
          frameon=False, handlelength=0.9, handletextpad=0.4, labelspacing=0.35)
panel(ax, "C", dx=-0.12)
save(fig, "Fig1_atlas")

# ================================================================ Fig 2
fib = pd.read_csv(f"{RES}/plot_fibro.csv")
de = pd.read_csv(f"{RES}/s5_pseudobulk_de.csv")
pr = pd.read_csv(f"{RES}/plot_fibro_programmes.csv", index_col=0)

fig = plt.figure(figsize=(11.2, 8.4))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.35, 1.35, 1.0], height_ratios=[1.0, 0.95],
              wspace=0.34, hspace=0.34, left=0.055, right=0.975, top=0.95, bottom=0.075)

ax = fig.add_subplot(gs[0, 0:2])
for st in ["SMC_synthetic", "Pericyte", "SMC_contractile", "AF_stromal", "Myofibroblast", "AF_resident"]:
    s = fib[fib.subtype == st]
    ax.scatter(s.U1, s.U2, s=0.4, c=ST_COL[st], linewidths=0, rasterized=True, alpha=0.8)
cent = fib.groupby("subtype")[["U1", "U2"]].median()
repel_labels(ax, cent.loc[ST_ORDER, "U1"].values, cent.loc[ST_ORDER, "U2"].values,
             [ST_LAB[s] for s in ST_ORDER], [ST_COL[s] for s in ST_ORDER])
umap_frame(ax)
ax.set_title("30,009 mural and fibroblast cells, six states", pad=4)
panel(ax, "A")

ax = fig.add_subplot(gs[0, 2])
pcols = ["pr_Contractile", "pr_Synthetic_MyoF", "pr_Adventitial", "pr_TGFb", "pr_IFN",
         "pr_Matrisome", "pr_Crosslink"]
plab = ["Contractile", "Synthetic\nmyofib.", "Adventitial", "TGF-b", "IFN",
        "Matrisome", "Crosslink"]
M = pr.loc[ST_ORDER, pcols].values.astype(float)
Mz = (M - M.mean(0)) / M.std(0)
vmax = np.abs(Mz).max()
im = ax.imshow(Mz, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_yticks(range(len(ST_ORDER)))
ax.set_yticklabels([ST_LAB[s] for s in ST_ORDER], fontsize=6.2)
ax.set_xticks(range(len(pcols)))
ax.set_xticklabels(plab, rotation=45, ha="right", fontsize=6.2)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("z score", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
ax.set_title("Programme score", pad=4)
panel(ax, "B", dx=-0.30)

ax = fig.add_subplot(gs[1, :])
lin = de[de.level == "lineage"].copy()
lin["nl"] = -np.log10(lin.padj.clip(lower=1e-300))
sig = (lin.padj < 0.05) & (lin.log2FoldChange.abs() > 0.585)
ax.scatter(lin.loc[~sig, "log2FoldChange"], lin.loc[~sig, "nl"], s=2.5, c="#CFCFCF",
           linewidths=0, rasterized=True, label="Not significant")
ax.scatter(lin.loc[sig, "log2FoldChange"], lin.loc[sig, "nl"], s=14, c=RED,
           edgecolor="white", linewidth=0.3, label=f"FDR < 0.05 (n = {int(sig.sum())})")
ax.axhline(-np.log10(0.05), ls="--", lw=0.6, c="#666666")
for xv in (-0.585, 0.585):
    ax.axvline(xv, ls=":", lw=0.6, c="#666666")
ax.set_xlabel("log$_2$ fold change (failed vs matured)")
ax.set_ylabel("-log$_{10}$ FDR")
ax.set_xlim(-9.5, 5.5)
ax.set_ylim(-0.05, 2.4)
lab = lin[sig].sort_values("padj").head(4)
for _, r in lab.iterrows():
    ax.annotate(r.gene, (r.log2FoldChange, r.nl), textcoords="offset points",
                xytext=(6, 4), fontsize=6.3, color="#333333")
ax.legend(frameon=False, loc="upper left", handletextpad=0.3)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Donor-level pseudobulk within the fibroblast lineage (24,157 genes)", pad=4)
panel(ax, "C", dx=-0.035)
save(fig, "Fig2_subclusters")

# ================================================================ Fig 3
comp = pd.read_csv(f"{RES}/s5_composition_tests.csv")
fr = pd.read_csv(f"{RES}/s5_subtype_fraction_plot.csv")
fr = fr[fr.group.isin(["failed_AVF", "mature_AVF"])]

fig = plt.figure(figsize=(11.0, 4.1))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.55, 1.0], wspace=0.30,
              left=0.065, right=0.985, top=0.90, bottom=0.20)

ax = fig.add_subplot(gs[0, 0])
ordr = comp.sort_values("diff").state.tolist()
rng = np.random.default_rng(7)
for i, st in enumerate(ordr):
    for k, (g, col) in enumerate([("mature_AVF", GROUP_COL["mature_AVF"]),
                                  ("failed_AVF", GROUP_COL["failed_AVF"])]):
        v = fr[fr.group == g][st].values * 100
        xj = i + (k - 0.5) * 0.30 + rng.uniform(-0.055, 0.055, len(v))
        ax.scatter(xj, v, s=13, c=col, edgecolor="white", linewidth=0.35, zorder=3)
        ax.hlines(np.median(v), i + (k - 0.5) * 0.30 - 0.10, i + (k - 0.5) * 0.30 + 0.10,
                  color=col, lw=1.4, zorder=4)
ax.set_xticks(range(len(ordr)))
ax.set_xticklabels([ST_LAB[s].replace(" ", "\n") for s in ordr], fontsize=6.5)
ax.set_ylabel("% of mural and fibroblast cells")
ax.set_ylim(-2, 76)
ax.spines[["top", "right"]].set_visible(False)
handles = [Line2D([], [], marker="o", ls="", ms=4, color=GROUP_COL[g], label=GROUP_LAB[g])
           for g in ("mature_AVF", "failed_AVF")]
ax.legend(handles=handles, frameon=False, loc="upper right", ncol=2, handletextpad=0.2)
ax.set_title("State composition, one point per donor (n = 6 + 6)", pad=4)
panel(ax, "A")

ax = fig.add_subplot(gs[0, 1])
c2 = comp.set_index("state").loc[ordr]
ypos = np.arange(len(ordr))
for i, (st, r) in enumerate(c2.iterrows()):
    col = RED if r["diff"] < 0 else "#4C72B0"
    ax.hlines(i, r.ci_low, r.ci_high, color=col, lw=1.1)
    ax.plot([r.ci_low, r.ci_low], [i - 0.13, i + 0.13], color=col, lw=0.9)
    ax.plot([r.ci_high, r.ci_high], [i - 0.13, i + 0.13], color=col, lw=0.9)
    ax.plot(r["diff"], i, "o", ms=4.2, color=col, zorder=4)
    ax.text(2.55, i, f"d = {r.cohens_d:+.2f}, p = {r.welch_p:.2f}", fontsize=6.2,
            va="center", ha="left")
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(ypos)
ax.set_yticklabels([ST_LAB[s] for s in ordr], fontsize=6.8)
ax.set_xlim(-2.2, 4.1)
ax.set_xticks([-2, -1, 0, 1, 2])
ax.set_xlabel("Difference in CLR abundance (failed - matured)")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Effect size with 95% CI", pad=4)
panel(ax, "B", dx=-0.42)
save(fig, "Fig3_composition")

# ================================================================ Fig 4
ft = pd.read_csv(f"{RES}/s7_fate_tests.csv")
fd = pd.read_csv(f"{RES}/s7_donor_fate.csv")
gt = pd.read_csv(f"{RES}/s7_gene_trend.csv")
tu = pd.read_csv(f"{RES}/plot_traj_umap.csv")
paga = pd.read_csv(f"{RES}/plot_paga.csv", index_col=0)
ppos = pd.read_csv(f"{RES}/plot_paga_pos.csv", index_col=0)
sizes = tu.subtype.value_counts()

fig = plt.figure(figsize=(11.2, 7.9))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0],
              wspace=0.30, hspace=0.34, left=0.06, right=0.925, top=0.945, bottom=0.075)

ax = fig.add_subplot(gs[0, 0])
nn = paga.index.tolist()
for i, a_ in enumerate(nn):
    for j, b_ in enumerate(nn):
        if j <= i:
            continue
        w = paga.loc[a_, b_]
        if w >= 0.05:
            ax.plot([ppos.loc[a_, "x"], ppos.loc[b_, "x"]], [ppos.loc[a_, "y"], ppos.loc[b_, "y"]],
                    color="#7A7A7A", lw=0.4 + 5.0 * w, alpha=0.75, zorder=1)
for s in nn:
    ax.scatter(ppos.loc[s, "x"], ppos.loc[s, "y"], s=40 + 1.7 * np.sqrt(sizes[s]),
               color=ST_COL[s], edgecolor="white", linewidth=0.8, zorder=3)
off = {"AF_resident": (-0.02, 0.26), "AF_stromal": (0.05, 0.24), "Myofibroblast": (0.0, -0.30),
       "Pericyte": (0.14, 0.22), "SMC_contractile": (0.06, -0.28), "SMC_synthetic": (0.02, 0.24)}
for s in nn:
    dx, dy = off[s]
    ax.text(ppos.loc[s, "x"] + dx, ppos.loc[s, "y"] + dy, ST_LAB[s], fontsize=6.3,
            ha="center", va="center")
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("PAGA topology", pad=4)
panel(ax, "A", dx=-0.02)

ax = fig.add_subplot(gs[0, 1])
o = np.argsort(tu.dpt_pseudotime.values)
pt = tu.dpt_pseudotime.values
rp = pt.argsort().argsort() / (len(pt) - 1) * 100
sc = ax.scatter(tu.U1.values[o], tu.U2.values[o], s=0.4, c=rp[o],
                cmap="viridis", linewidths=0, rasterized=True, vmin=0, vmax=100)
umap_frame(ax)
cb = fig.colorbar(sc, ax=ax, fraction=0.042, pad=0.02)
cb.set_label("Pseudotime percentile", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
ax.set_title("Pseudotime from resident adventitial fibroblast root", pad=4)
panel(ax, "B", dx=-0.02)

ax = fig.add_subplot(gs[1, 0])
flab = {"Myofibroblast": "Myofibroblast", "SMC_contractile": "Contractile SMC",
        "AF_resident": "Resident AF"}
yy = np.arange(len(ft))[::-1]
for i, (_, r) in enumerate(ft.iterrows()):
    f = fd[fd.group == "failed_AVF"]["abs_" + r.fate].values
    m = fd[fd.group == "mature_AVF"]["abs_" + r.fate].values
    d_ = r["diff"]
    se = np.sqrt(f.var(ddof=1) / len(f) + m.var(ddof=1) / len(m))
    dfw = (f.var(ddof=1) / len(f) + m.var(ddof=1) / len(m)) ** 2 / (
        (f.var(ddof=1) / len(f)) ** 2 / (len(f) - 1) + (m.var(ddof=1) / len(m)) ** 2 / (len(m) - 1))
    tcrit = stats.t.ppf(0.975, dfw)
    y = yy[i]
    ax.hlines(y, d_ - tcrit * se, d_ + tcrit * se, color="#4C72B0", lw=1.1)
    ax.plot(d_, y, "o", ms=4.4, color="#4C72B0")
    ax.text(-0.205, y + 0.34, f"d = {r.cohens_d:+.2f}, p = {r.welch_p:.2f}", fontsize=6.3,
            va="center", ha="left", color="#333333")
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels([flab[f] for f in ft.fate], fontsize=6.8)
ax.set_xlim(-0.215, 0.325)
ax.set_ylim(-0.55, 2.62)
ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2, 0.3])
ax.set_xlabel("Difference in absorption probability (failed - matured)")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Fate bias, donor level", pad=4)
panel(ax, "C", dx=-0.30)

ax = fig.add_subplot(gs[1, 1])
genes = ["PI16", "DCN", "POSTN", "COL8A1", "FN1", "ACTA2", "MYH11"]
gene_col = {"PI16": "#2E7D6F", "DCN": "#4C9A8F", "POSTN": "#C0504D",
            "COL8A1": "#D98C4A", "FN1": "#8E5AA8", "ACTA2": "#3B6FB6", "MYH11": "#1F3F66"}


def bin_edges(s):
    """bin 形如 '(-0.001, 0.00267]'，取其右端点作为该分位的拟时序上界。"""
    return np.array([float(str(b).split(",")[1].rstrip("]")) for b in s])


for k, g in enumerate(genes):
    s = gt[gt.gene == g].copy()
    s["pt_hi"] = bin_edges(s["bin"])
    s = s.sort_values("pt_hi")            # 按拟时序排序，不能用表达量排序
    xp = np.arange(len(s)) * (100.0 / len(s)) + 50.0 / len(s)   # 等频分箱的中点百分位
    ax.plot(xp, s.x.values, lw=1.2, color=gene_col[g], label=g)
ax.set_xlabel("Pseudotime percentile")
ax.set_xlim(-3, 103)
ax.set_xlabel("Pseudotime percentile")
ax.set_ylabel("Mean expression (log-normalised)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), handlelength=1.2)
ax.set_title("Marker dynamics along pseudotime", pad=4)
panel(ax, "D", dx=-0.22)
save(fig, "Fig4_trajectory")

# ================================================================ Fig 5
ax_tab = pd.read_csv(f"{RES}/s8c_pathway_tests.csv").sort_values("log2fc")
lr = pd.read_csv(f"{RES}/s8_lr_tests.csv").sort_values("welch_p").head(12)

fig = plt.figure(figsize=(12.6, 4.6))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.35], wspace=0.62,
              left=0.13, right=0.985, top=0.90, bottom=0.14)

ax = fig.add_subplot(gs[0, 0])
yy = np.arange(len(ax_tab))
for i, (_, r) in enumerate(ax_tab.iterrows()):
    col = RED if r.log2fc > 0 else "#4C72B0"
    ax.barh(i, r.log2fc, height=0.62, color=col, edgecolor="white", linewidth=0.4)
    ax.text(0.92, i, f"d = {r.cohens_d:+.2f}, p = {r.welch_p:.3f}",
            va="center", ha="left", fontsize=6.0)
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels(ax_tab.axis, fontsize=6.6)
ax.set_xlabel("log$_2$ ratio of axis score (failed / matured)")
ax.set_xlim(-0.95, 1.75)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Nine pre-specified signalling axes", pad=4)
panel(ax, "A", dx=-0.34)

ax = fig.add_subplot(gs[0, 1])
lbl = [f"{r.sender.replace('Fib_', '').replace('_', ' ')} \u2192 {r.receiver.replace('Fib_', '').replace('_', ' ')}  "
       f"{r.ligand}\u2013{r.receptor}" for _, r in lr.iterrows()]
yy = np.arange(len(lr))[::-1]
for i, (_, r) in enumerate(lr.iterrows()):
    col = RED if r.log2fc > 0 else "#4C72B0"
    ax.hlines(yy[i], 0, r.log2fc, color=col, lw=0.9, alpha=0.8)
    ax.plot(r.log2fc, yy[i], "o", ms=4.0, color=col)
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels(lbl, fontsize=6.0)
ax.set_xlabel("log$_2$ ratio of interaction score (failed / matured)")
ax.set_xlim(-2.0, 1.8)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Twelve most significant individual interactions (exploratory)", pad=4)
panel(ax, "B", dx=-0.62)
save(fig, "Fig5_communication")

# ================================================================ Fig 6
cc = pd.read_csv(f"{RES}/s5g_compartment_counts_before_after.csv")
ecm = pd.read_csv(f"{RES}/s9_ecm_tests.csv")
bulk = pd.read_csv(f"{RES}/s10_focus_genes.csv", index_col=0)

fig = plt.figure(figsize=(13.6, 4.5))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.0], wspace=0.52,
              left=0.10, right=0.985, top=0.90, bottom=0.16)

ax = fig.add_subplot(gs[0, 0])
cc = cc.sort_values("n_de_before")
yy = np.arange(len(cc))
ax.barh(yy + 0.18, cc.n_de_before, height=0.34, color="#BDBDBD", label="Before screening")
ax.barh(yy - 0.18, cc.n_de_after, height=0.34, color="#4C72B0", label="After screening")
for i, r in enumerate(cc.itertuples()):
    ax.text(r.n_de_before + 0.3, i + 0.18, str(r.n_de_before), va="center", fontsize=6.0)
    ax.text(r.n_de_after + 0.3, i - 0.18, str(r.n_de_after), va="center", fontsize=6.0)
ax.set_yticks(yy)
ax.set_yticklabels([s.replace("Fib_", "").replace("_", " ") for s in cc.compartment], fontsize=6.5)
ax.set_xlabel("Differentially expressed genes")
ax.set_xlim(0, 19)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="lower right", handlelength=1.0, labelspacing=0.3)
ax.set_title("Genes per compartment, before and after\ncell-of-origin screening", pad=4)
panel(ax, "A", dx=-0.34)

ax = fig.add_subplot(gs[0, 1])
keep = ["ECM1", "ECM2", "ECM3", "ECM4", "ECM5", "zonation"]
ecm = ecm[ecm.metric.isin(keep)].set_index("metric").loc[keep]
elab = ["ECM1 interstitial", "ECM2 basement\nmembrane", "ECM3 microfibril",
        "ECM4 fibrillar\ncollagen", "ECM5 fibronectin/\nversican", "Zonation index"]
yy = np.arange(len(ecm))[::-1]
for i, (_, r) in enumerate(ecm.iterrows()):
    col = RED if r["diff"] > 0 else "#4C72B0"
    ax.barh(yy[i], r["diff"], height=0.6, color=col, edgecolor="white", linewidth=0.4)
    ax.text(0.082, yy[i], f"p = {r.p:.2f}", fontsize=6.1, va="center")
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels(elab, fontsize=6.2)
ax.set_xlim(-0.03, 0.115)
ax.set_xlabel("Difference in programme usage (failed - matured)")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Matrisome programmes and computed\nzonation", pad=4)
panel(ax, "B", dx=-0.42)

ax = fig.add_subplot(gs[0, 2])
genes = ["COL8A1", "LOXL1", "FN1", "MFAP5", "COL1A1", "ACTA2", "MYH11", "DCN",
         "COL3A1", "POSTN", "PI16", "TIMP1"]
b = bulk.loc[[g for g in genes if g in bulk.index]]
b = b.sort_values("log2FoldChange")
yy = np.arange(len(b))
for i, (g, r) in enumerate(b.iterrows()):
    col = RED if r.log2FoldChange > 0 else "#4C72B0"
    ax.barh(i, r.log2FoldChange, height=0.6, color=col, edgecolor="white", linewidth=0.4)
    fdr = r.padj
    star = "*" if fdr < 0.05 else ""
    ax.text(1.18, i, f"FDR = {fdr:.2f}{star}", fontsize=5.9, va="center",
            fontweight="bold" if fdr < 0.05 else "normal")
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels(b.index, fontsize=6.6)
ax.set_xlim(-1.35, 2.0)
ax.set_xlabel("log$_2$ fold change (failed vs matured)")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Whole-tissue bulk, GSE220796\n(11 failed vs 9 matured)", pad=4)
panel(ax, "C", dx=-0.36)
save(fig, "Fig6_compartments")

# ================================================================ Fig 7
my = pd.read_csv(f"{RES}/plot_myeloid.csv")
mk = pd.read_csv(f"{RES}/s5e_subtype_marker_check.csv").set_index("subtype")
ab = pd.read_csv(f"{RES}/s5e_abundance_tests.csv").set_index("metric")
don = pd.read_csv(f"{RES}/s5e_donor_abundance.csv")

fig = plt.figure(figsize=(11.2, 8.0))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.25, 1.0, 0.95], height_ratios=[1.0, 1.0],
              wspace=0.42, hspace=0.42, left=0.06, right=0.975, top=0.945, bottom=0.075)

ax = fig.add_subplot(gs[0, 0:2])
for st in ["Low_quality", "Neutrophil", "pDC", "cDC2", "Mono_classical", "Macrophage"]:
    s = my[my.mye_subtype == st]
    ax.scatter(s.U1, s.U2, s=0.5, c=MY_COL[st], linewidths=0, rasterized=True, alpha=0.75)
cent = my.groupby("mye_subtype")[["U1", "U2"]].median()
repel_labels(ax, cent.loc[MY_ORDER, "U1"].values, cent.loc[MY_ORDER, "U2"].values,
             [MY_LAB[s] for s in MY_ORDER], [MY_COL[s] for s in MY_ORDER])
umap_frame(ax)
ax.set_title("14,845 myeloid cells, gated annotation", pad=4)
panel(ax, "A")

ax = fig.add_subplot(gs[0, 2])
cols = ["Macrophage", "Mono_classical", "cDC2", "pDC", "Neutrophil", "Low_quality"]
rows = ["CD14", "FCGR3A", "C1QA", "CD1C", "CLEC10A", "FCER1A", "LILRA4", "CLEC9A", "XCR1", "TPSAB1"]
M = mk.loc[cols, rows].T.values.astype(float)
im = ax.imshow(np.clip(M, 0, 2.5), cmap="Blues", vmin=0, vmax=2.5, aspect="auto")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([MY_LAB[c].replace(" ", "\n") for c in cols], fontsize=5.8, rotation=45, ha="right")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(rows, fontsize=6.0)
for t, g in zip(ax.get_yticklabels(), rows):
    if g in ("CLEC9A", "XCR1"):
        t.set_color(RED)
        t.set_fontweight("bold")
for i in range(len(rows)):
    for j in range(len(cols)):
        ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=4.9,
                color="white" if M[i, j] > 1.4 else "#222222")
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("Mean log-normalised", fontsize=6.2)
cb.ax.tick_params(labelsize=5.8)
ax.set_title("State markers (red: at ambient level)", pad=4)
panel(ax, "B", dx=-0.42)

ax = fig.add_subplot(gs[1, 0])
for g, lab in [("failed_AVF", "Failed"), ("mature_AVF", "Matured")]:
    v = don[don.group == g]["cdc2_D1_per1000"].values
    ax.scatter(np.full(len(v), 0 if g == "failed_AVF" else 1) + np.random.default_rng(3).uniform(-0.07, 0.07, len(v)),
               v, s=16, c=GROUP_COL[g], edgecolor="white", linewidth=0.4, zorder=3)
    ax.hlines(np.mean(v), (0 if g == "failed_AVF" else 1) - 0.16, (0 if g == "failed_AVF" else 1) + 0.16,
              color=GROUP_COL[g], lw=1.4, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Failed", "Matured"])
ax.set_xlim(-0.45, 1.45)
ax.set_ylabel("cDC2 per 1,000 tissue cells")
ax.spines[["top", "right"]].set_visible(False)
_a = don.loc[don.group == "failed_AVF", "cdc2_D1_per1000"].values
_b = don.loc[don.group == "mature_AVF", "cdc2_D1_per1000"].values
ax.set_title(f"cDC2 abundance ({d_p(_a, _b)})", pad=4, fontsize=7.8)
panel(ax, "C", dx=-0.30)

ax = fig.add_subplot(gs[1, 1])
for g, lab in [("failed_AVF", "Failed"), ("mature_AVF", "Matured")]:
    v = don[don.group == g]["cdc2_D1_frac_myeloid"].values * 100
    ax.scatter(np.full(len(v), 0 if g == "failed_AVF" else 1) + np.random.default_rng(5).uniform(-0.07, 0.07, len(v)),
               v, s=16, c=GROUP_COL[g], edgecolor="white", linewidth=0.4, zorder=3)
    ax.hlines(np.mean(v), (0 if g == "failed_AVF" else 1) - 0.16, (0 if g == "failed_AVF" else 1) + 0.16,
              color=GROUP_COL[g], lw=1.4, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Failed", "Matured"])
ax.set_xlim(-0.45, 1.45)
ax.set_ylabel("cDC2, % of myeloid cells")
ax.spines[["top", "right"]].set_visible(False)
_a = don.loc[don.group == "failed_AVF", "cdc2_D1_frac_myeloid"].values
_b = don.loc[don.group == "mature_AVF", "cdc2_D1_frac_myeloid"].values
ax.set_title(f"Dilution within myeloid ({d_p(_a, _b)})", pad=4, fontsize=7.8)
panel(ax, "D", dx=-0.30)

ax = fig.add_subplot(gs[1, 2])
for g, lab in [("failed_AVF", "Failed"), ("mature_AVF", "Matured")]:
    v = don[don.group == g]["myeloid_per1000"].values
    ax.scatter(np.full(len(v), 0 if g == "failed_AVF" else 1) + np.random.default_rng(11).uniform(-0.07, 0.07, len(v)),
               v, s=16, c=GROUP_COL[g], edgecolor="white", linewidth=0.4, zorder=3)
    ax.hlines(np.mean(v), (0 if g == "failed_AVF" else 1) - 0.16, (0 if g == "failed_AVF" else 1) + 0.16,
              color=GROUP_COL[g], lw=1.4, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Failed", "Matured"])
ax.set_xlim(-0.45, 1.45)
ax.set_ylabel("Myeloid cells per 1,000")
ax.spines[["top", "right"]].set_visible(False)
_a = don.loc[don.group == "failed_AVF", "myeloid_per1000"].values
_b = don.loc[don.group == "mature_AVF", "myeloid_per1000"].values
ax.set_title(f"Myeloid expansion ({d_p(_a, _b)})", pad=4, fontsize=7.8)
panel(ax, "E", dx=-0.30)
save(fig, "Fig7_myeloid_cdc2")

# ================================================================ Fig 8
dc = pd.read_csv(f"{RES}/s5f_decomposition_retained.csv")
det = pd.read_csv(f"{RES}/s5e_cdc2_downsampled_detection.csv")
dp = pd.read_csv(f"{RES}/s5e_cdc2_depth_by_group.csv").set_index("group")

fig = plt.figure(figsize=(12.4, 4.1))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.35, 0.95, 0.62], wspace=0.42,
              left=0.075, right=0.985, top=0.90, bottom=0.18)

ax = fig.add_subplot(gs[0, 0])
dc = dc.sort_values("log2FC_observed")
yy = np.arange(len(dc))
for i, (_, r) in enumerate(dc.iterrows()):
    ax.hlines(i, r.log2FC_predicted, r.log2FC_observed, color="#BBBBBB", lw=1.2, zorder=1)
    ax.plot(r.log2FC_observed, i, "o", ms=4.6, color="#333333", zorder=3,
            label="Observed" if i == 0 else None)
    ax.plot(r.log2FC_predicted, i, "X", ms=4.6, color="#E0A33E", zorder=3,
            label="Predicted (composition + intrinsic)" if i == 0 else None)
ax.axvline(0, color="#444444", lw=0.7)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.gene} ({r.subtype.replace('Mono_classical', 'mono').replace('cDC2', 'cDC2')})"
                    for _, r in dc.iterrows()], fontsize=6.3)
ax.set_xlabel("log$_2$ fold change (failed vs matured)")
ax.set_xlim(-3.2, 2.9)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="lower right", handletextpad=0.3)
ax.set_title("Compartment signal decomposed", pad=4)
panel(ax, "A", dx=-0.42)

ax = fig.add_subplot(gs[0, 1])
for g, lab in [("failed_AVF", "Failed"), ("mature_AVF", "Matured")]:
    v = det[det.group == g]["det_CD1C"].values * 100
    ax.scatter(np.full(len(v), 0 if g == "failed_AVF" else 1) + np.random.default_rng(9).uniform(-0.08, 0.08, len(v)),
               v, s=17, c=GROUP_COL[g], edgecolor="white", linewidth=0.4, zorder=3)
    ax.hlines(np.mean(v), (0 if g == "failed_AVF" else 1) - 0.17, (0 if g == "failed_AVF" else 1) + 0.17,
              color=GROUP_COL[g], lw=1.4, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Failed", "Matured"])
ax.set_xlim(-0.45, 1.45)
ax.set_ylabel("CD1C$^+$ cells within cDC2 (%)")
ax.spines[["top", "right"]].set_visible(False)
_a = det.loc[det.group == "failed_AVF", "det_CD1C"].values * 100
_b = det.loc[det.group == "mature_AVF", "det_CD1C"].values * 100
ax.set_title(f"At matched depth, 5,000 UMIs\n({d_p(_a, _b)})", pad=4, fontsize=7.8)
panel(ax, "B", dx=-0.32)

ax = fig.add_subplot(gs[0, 2])
vals = [dp.loc["failed_AVF", "total_counts"], dp.loc["mature_AVF", "total_counts"]]
bars = ax.bar([0, 1], vals, width=0.6, color=[GROUP_COL["failed_AVF"], GROUP_COL["mature_AVF"]],
              edgecolor="white", linewidth=0.5)
for x, v in zip([0, 1], vals):
    ax.text(x, v + 350, f"{v/1000:.1f}k", ha="center", fontsize=6.6)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Failed", "Matured"])
ax.set_ylim(0, 20000)
ax.set_ylabel("Median UMIs per cDC2 cell")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Depth imbalance before\ncorrection", pad=4, fontsize=7.8)
panel(ax, "C", dx=-0.44)
save(fig, "Fig8_signal_decomposition")

# ================================================================ Fig S1
scr = pd.read_csv(f"{RES}/s5g_compartment_gene_screen.csv")
scr["low"] = scr.focal_mean < 0.10
scr = scr.sort_values(["compartment", "log2FoldChange"])
fig, ax = plt.subplots(figsize=(7.4, 7.2))
yy = np.arange(len(scr))
for i, (_, r) in enumerate(scr.iterrows()):
    if r.retained and not r.low:
        fc, ec, lab = "#4C72B0", "white", "Retained, above expression floor"
    elif r.retained and r.low:
        fc, ec, lab = "#9DC3E6", "#4C72B0", "Retained, below expression floor"
    else:
        fc, ec, lab = "#D9D9D9", "#7F7F7F", "Removed as ambient/non-target"
    ax.barh(i, r.log2FoldChange, height=0.66, color=fc, edgecolor=ec, linewidth=0.6)
ax.axvline(0, color="#444444", lw=0.8)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.gene} ({r.compartment.replace('Fib_', '').replace('_', ' ')})"
                    for _, r in scr.iterrows()], fontsize=6.6)
ax.set_xlabel("log$_2$ fold change (failed vs matured)")
ax.set_xlim(-6.6, 5.2)
ax.spines[["top", "right"]].set_visible(False)
handles = [Patch(facecolor="#4C72B0", edgecolor="white", label="Retained, above floor"),
           Patch(facecolor="#9DC3E6", edgecolor="#4C72B0", label="Retained, below floor"),
           Patch(facecolor="#D9D9D9", edgecolor="#7F7F7F", label="Removed (ambient / non-target)")]
ax.legend(handles=handles, frameon=False, loc="lower right", fontsize=6.3)
ax.set_title("All compartment-resolved differentially expressed genes (n = %d)" % len(scr), pad=6)
save(fig, "FigS1_compartment_genes")
