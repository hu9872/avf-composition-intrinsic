"""S14: 主图生成（图 1–6）。所有图严格基于 results/ 与 data/processed/ 的真实输出。"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
FIG = r"D:/AVF_GSE_2611/figures"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "savefig.bbox": "tight", "figure.dpi": 300})
PAL = {"pre_access": "#4C72B0", "early_AVF": "#8C8C8C",
       "mature_AVF": "#55A868", "failed_AVF": "#C44E52"}
GRP_ORDER = ["pre_access", "early_AVF", "mature_AVF", "failed_AVF"]
FAIL, MAT = "#C44E52", "#55A868"


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
    print("saved", name)


def fig1(a):
    fig, ax = plt.subplots(1, 3, figsize=(12.6, 3.7),
                           gridspec_kw={"width_ratios": [1, 1, 1.15]})
    sc.pl.umap(a, color="celltype", ax=ax[0], show=False, legend_loc="on data",
               legend_fontsize=5.5, frameon=False, title="Cell types", size=3)
    sc.pl.umap(a, color="group", ax=ax[1], show=False, palette=PAL, frameon=False,
               title="Sampling stage", legend_fontsize=6, size=3)
    comp = (a.obs.groupby(["group", "celltype"], observed=True).size().unstack(fill_value=0))
    comp = comp.div(comp.sum(axis=1), axis=0).reindex(
        [g for g in GRP_ORDER if g in comp.index])
    comp.plot(kind="bar", stacked=True, ax=ax[2], width=0.75, legend=False,
              colormap="tab20", edgecolor="none")
    ax[2].set_ylabel("Fraction of cells"); ax[2].set_xlabel("")
    ax[2].tick_params(axis="x", rotation=25, labelsize=7)
    ax[2].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=5, frameon=False)
    for _ax in ax[:2]:
        _ax.set_xlabel("UMAP-1"); _ax.set_ylabel("UMAP-2")
    save(fig, "Fig1_atlas")


def fig2(sub, pb):
    fig = plt.figure(figsize=(12.6, 6.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05])
    ax0 = fig.add_subplot(gs[0, 0])
    sc.pl.umap(sub, color="subtype", ax=ax0, show=False, legend_loc="on data",
               legend_fontsize=5, frameon=False, title="Mural/fibroblast states", size=3)
    ax1 = fig.add_subplot(gs[0, 1])
    sc.pl.umap(sub, color="group", ax=ax1, show=False, palette=PAL, frameon=False,
               title="Sampling stage", legend_fontsize=5, size=3)
    ax2 = fig.add_subplot(gs[0, 2])
    cols = [c for c in sub.obs.columns if c.startswith("pr_")]
    df = sub.obs.groupby("subtype", observed=True)[cols].mean()
    df.columns = [c.replace("pr_", "") for c in df.columns]
    z = (df - df.mean()) / (df.std() + 1e-9)
    im = ax2.imshow(z, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax2.set_xticks(range(df.shape[1]))
    ax2.set_xticklabels(df.columns, rotation=45, ha="right", fontsize=6)
    ax2.set_yticks(range(df.shape[0])); ax2.set_yticklabels(df.index, fontsize=6)
    ax2.set_title("Programme score (z across states)"); fig.colorbar(im, ax=ax2, fraction=0.045)

    ax3 = fig.add_subplot(gs[1, :])
    if pb is not None and len(pb):
        t = pb.copy()
        if "unit" in t.columns:
            pref = [u for u in t["unit"].unique() if "Myofibroblast" in str(u)] or ["lineage"]
            t = t[t["unit"] == pref[0]].copy()
            ax3.set_title(f"Donor-level pseudobulk DE, {pref[0]} (failed vs matured) — "
                          f"{int(((t.padj<0.05)&(t.log2FoldChange.abs()>0.585)).sum())} genes at FDR<0.05")
        t = t.dropna(subset=["padj"]).drop_duplicates("gene").set_index("gene")
        t["-log10padj"] = -np.log10(t["padj"].clip(lower=1e-300))
        sig = (t.padj < 0.05) & (t.log2FoldChange.abs() > 0.585)
        ax3.scatter(t.log2FoldChange, t["-log10padj"], s=3, c="#B0B0B0", linewidths=0)
        ax3.scatter(t.loc[sig, "log2FoldChange"], t.loc[sig, "-log10padj"],
                    s=8, c=FAIL, linewidths=0)
        for g in ["COL8A1", "POSTN", "TIMP1", "PI16", "LOXL1", "FN1"]:
            if g in t.index:
                r = t.loc[g]
                ax3.annotate(g, (r.log2FoldChange, -np.log10(max(r.padj, 1e-300))), fontsize=6)
        ax3.axhline(-np.log10(0.05), ls="--", lw=0.5, c="k")
        ax3.axvline(0, lw=0.5, c="k")
        ax3.set_xlabel("log2 fold-change (failed vs matured)"); ax3.set_ylabel("-log10 FDR")
        if not ax3.get_title():
            ax3.set_title("Donor-level pseudobulk differential expression")
    save(fig, "Fig2_subclusters")


def fig3(sub):
    post = sub.obs[sub.obs.group.isin(["mature_AVF", "failed_AVF"])]
    counts = post.groupby(["donor", "group", "subtype"], observed=True).size().unstack(fill_value=0)
    frac = counts.div(counts.sum(axis=1), axis=0)
    frac.to_csv(f"{RES}/s5_subtype_fraction_plot.csv")
    states = [c for c in frac.columns]
    order = ["AF_resident", "AF_stromal", "Myofibroblast", "Pericyte",
             "SMC_contractile", "SMC_synthetic"]
    states = [s for s in order if s in states] + [s for s in states if s not in order]

    fig, ax = plt.subplots(1, 2, figsize=(11.2, 3.6),
                           gridspec_kw={"width_ratios": [1.45, 1]})
    data, pos, colors = [], [], []
    for i, s in enumerate(states):
        for j, (g, c) in enumerate([("mature_AVF", MAT), ("failed_AVF", FAIL)]):
            idx = frac.index.get_level_values("group") == g
            data.append(frac.loc[idx, s].values); pos.append(i * 3 + j + 1); colors.append(c)
    bp = ax[0].boxplot(data, positions=pos, widths=0.72, patch_artist=True, showfliers=False)
    for b, c in zip(bp["boxes"], colors):
        b.set_facecolor(c); b.set_linewidth(0.5); b.set_alpha(0.85)
    for m in bp["medians"]:
        m.set_color("white"); m.set_linewidth(1)
    for i in range(len(states)):
        for j, c in enumerate([MAT, FAIL]):
            v = data[i * 2 + j]
            ax[0].scatter(np.full(len(v), i * 3 + j + 1) + np.random.default_rng(1).normal(0, .05, len(v)),
                          v, s=6, c="k", zorder=3, linewidths=0)
    ax[0].set_xticks([i * 3 + 1.5 for i in range(len(states))])
    ax[0].set_xticklabels(states, rotation=25, ha="right", fontsize=7)
    ax[0].set_ylabel("Fraction of mural/fibroblast cells")
    ax[0].legend([Rectangle((0, 0), 1, 1, fc=MAT), Rectangle((0, 0), 1, 1, fc=FAIL)],
                 ["Matured (n=6)", "Failed (n=6)"], fontsize=6, frameon=False, loc="upper right")
    ax[0].set_title("State composition per donor")

    t = pd.read_csv(f"{RES}/s5_composition_tests.csv")
    t = t.sort_values("diff")
    y = np.arange(len(t))
    ax[1].barh(y, t["diff"], color=np.where(t["diff"] > 0, FAIL, MAT), height=0.6)
    ax[1].errorbar(t["diff"], y, xerr=[t["diff"] - t["ci_low"], t["ci_high"] - t["diff"]],
                   fmt="none", ecolor="k", lw=0.7, capsize=2)
    for yy, (v, pv, dv) in enumerate(zip(t["diff"], t["welch_p"], t.cohens_d)):
        ax[1].text(t["ci_high"].iloc[yy] + 0.05, yy, f"d={dv:+.2f}, p={pv:.2f}",
                   va="center", fontsize=5.5)
    ax[1].set_yticks(y)
    ax[1].set_yticklabels([s.replace("SMC_", "SMC-") for s in t.state], fontsize=7)
    ax[1].axvline(0, c="k", lw=0.6)
    ax[1].set_xlabel("Δ CLR abundance (failed − matured)")
    ax[1].set_title("Effect size with 95% CI")
    ax[1].set_xlim(t["ci_low"].min() - 0.4, t["ci_high"].max() + 1.1)
    fig.subplots_adjust(wspace=0.42)
    save(fig, "Fig3_composition")


def fig4(sub):
    fig = plt.figure(figsize=(13.0, 3.7))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.15, 1.15], wspace=0.5)
    ax0 = fig.add_subplot(gs[0, 0])
    sc.pl.paga(sub, ax=ax0, show=False, node_size_scale=0.9, fontsize=4.5,
               edge_width_scale=0.6, title="PAGA")
    ax1 = fig.add_subplot(gs[0, 1])
    sc.pl.umap(sub, color="dpt_pseudotime", ax=ax1, show=False, frameon=False,
               title="Diffusion pseudotime", size=3)
    ax2 = fig.add_subplot(gs[0, 2])
    t = pd.read_csv(f"{RES}/s7_fate_tests.csv")
    y = np.arange(len(t))
    ax2.barh(y, t["diff"], color=np.where(t["diff"] > 0, FAIL, MAT), height=0.55)
    for yy, r in t.iterrows():
        ax2.text(0, yy, f"  d={r.cohens_d:+.2f}, p={r.welch_p:.2f}", va="center", fontsize=6)
    ax2.set_yticks(y); ax2.set_yticklabels([f.replace("SMC_", "SMC-") for f in t.fate], fontsize=7)
    ax2.axvline(0, c="k", lw=0.6)
    ax2.set_xlabel("Δ absorption probability (failed − matured)")
    ax2.set_title("Fate bias (donor level)")
    ax2.set_xlim(min(t["diff"].min() - 0.01, -0.02), t["diff"].max() + 0.09)

    ax3 = fig.add_subplot(gs[0, 3])
    tr = pd.read_csv(f"{RES}/s7_gene_trend.csv")
    for g, grp in tr.groupby("gene"):
        grp = grp.sort_values("bin")
        ax3.plot(range(len(grp)), grp["x"], lw=1.1, label=g)
    ax3.legend(fontsize=5, frameon=False, ncol=2)
    ax3.set_xlabel("Pseudotime bin"); ax3.set_ylabel("Mean expression")
    ax3.set_title("Marker dynamics along pseudotime")
    save(fig, "Fig4_trajectory")


def fig5():
    p = f"{RES}/s8c_pathway_tests.csv"
    if not os.path.exists(p):
        return
    t = pd.read_csv(p).sort_values("log2fc").reset_index(drop=True)
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 3.6),
                           gridspec_kw={"width_ratios": [1, 1], "wspace": 0.42})
    y = np.arange(len(t))
    ax[0].barh(y, t["log2fc"], color=np.where(t["log2fc"] > 0, FAIL, MAT), height=0.6)
    for yy, r in t.iterrows():
        ax[0].text(0, yy, f"  d={r.cohens_d:+.2f}, p={r.welch_p:.3f}", va="center", fontsize=5.5)
    ax[0].set_yticks(y); ax[0].set_yticklabels(t.axis, fontsize=7)
    ax[0].axvline(0, c="k", lw=0.6)
    ax[0].set_xlabel("log2 ratio of axis score (failed / matured)")
    ax[0].set_title("Ligand–receptor axes (donor level)")
    ax[0].set_xlim(t["log2fc"].min() - 0.05, t["log2fc"].max() + 0.55)

    d = pd.read_csv(f"{RES}/s8_lr_tests.csv")
    d = d.reindex(d.perm_p.sort_values().index).head(14).iloc[::-1]
    lbl = (d.sender.str.replace("Fib_", "") + "→" + d.receiver.str.replace("Fib_", "") +
           "  " + d.ligand + "–" + d.receptor)
    ax[1].barh(lbl, d.log2fc, color=np.where(d.log2fc > 0, FAIL, MAT), height=0.62)
    for yy, (pv, dv) in enumerate(zip(d.perm_p, d.cohens_d)):
        ax[1].text(0, yy, f"  d={dv:+.2f}", va="center", fontsize=5)
    ax[1].axvline(0, c="k", lw=0.6)
    ax[1].tick_params(axis="y", labelsize=5)
    ax[1].set_xlabel("log2 ratio of interaction score (failed / matured)")
    ax[1].set_title("Top individual interactions (exploratory)")
    save(fig, "Fig5_communication")


def fig6():
    p = f"{RES}/s5g_compartment_counts_before_after.csv"
    if not os.path.exists(p):
        p = f"{RES}/s5c_compartment_de_summary.csv"
        t = pd.read_csv(p).sort_values("n_de")
        t = t.rename(columns={"n_de": "n_de_before"}).assign(n_de_after=t.n_de)
    else:
        t = pd.read_csv(p).sort_values("n_de_before")
    fig = plt.figure(figsize=(13.4, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.05], wspace=0.46)
    ax0 = fig.add_subplot(gs[0, 0])
    y = np.arange(len(t))
    ax0.barh(y + 0.19, t["n_de_before"], color="#9E9E9E", height=0.36,
             label="before screening", zorder=3)
    ax0.barh(y - 0.19, t["n_de_after"], color="#3B6FA0", height=0.36,
             label="retained after screening", zorder=3)
    for yy, r in t.reset_index(drop=True).iterrows():
        ax0.text(r.n_de_before + 0.2, yy + 0.19, f"{int(r.n_de_before)}",
                 va="center", fontsize=5.5, color="#555555")
        ax0.text(r.n_de_after + 0.2, yy - 0.19, f"{int(r.n_de_after)}",
                 va="center", fontsize=5.5, color="#3B6FA0")
    ax0.set_yticks(y); ax0.set_yticklabels(t.compartment, fontsize=6.3)
    ax0.set_xlabel("DE genes (FDR<0.05, |log2FC|>0.585)", fontsize=7)
    ax0.set_title("Where does the difference lie?", fontsize=8, loc="left")
    ax0.set_xlim(0, t.n_de_before.max() * 1.22)
    ax0.legend(fontsize=5.6, frameon=False, loc="lower right")
    ax0.tick_params(labelsize=6.3)
    for sp in ["top", "right"]:
        ax0.spines[sp].set_visible(False)

    ax1 = fig.add_subplot(gs[0, 1])
    e = pd.read_csv(f"{RES}/s9_ecm_tests.csv").sort_values("p")
    yy = np.arange(len(e))
    ax1.barh(yy, e["diff"], color=np.where(e["diff"] > 0, FAIL, MAT), height=0.6)
    for i, r in e.reset_index(drop=True).iterrows():
        ax1.text(0, i, f"  p={r.p:.2f}", va="center", fontsize=5.5)
    ax1.set_yticks(yy); ax1.set_yticklabels(e.metric, fontsize=6.5)
    ax1.axvline(0, c="k", lw=0.6)
    ax1.set_xlabel("Δ score (failed − matured)")
    ax1.set_title("ECM programmes / zonation (all n.s.)")
    ax1.set_xlim(e["diff"].min() - 0.01, e["diff"].max() + 0.05)

    ax2 = fig.add_subplot(gs[0, 2])
    f = f"{RES}/s10_focus_genes.csv"
    if os.path.exists(f):
        fr = pd.read_csv(f, index_col=0).sort_values("log2FoldChange")
        yy2 = np.arange(len(fr))
        cols = [FAIL if v > 0 else MAT for v in fr.log2FoldChange]
        cols = [c if p < 0.05 else "#B0B0B0" for c, p in zip(cols, fr.padj)]
        ax2.barh(yy2, fr.log2FoldChange, color=cols, height=0.62)
        for i, (v, p_) in enumerate(zip(fr.log2FoldChange, fr.padj)):
            ax2.text(v, i, f" FDR={p_:.2f}", va="center", fontsize=5,
                     ha="left" if v > 0 else "right")
        ax2.set_yticks(yy2); ax2.set_yticklabels(fr.index, fontsize=6.5)
        ax2.axvline(0, c="k", lw=0.6)
        ax2.set_xlabel("log2FC (failed vs matured)")
        ax2.set_title("Whole-tissue bulk (GSE220796)")
        ax2.set_xlim(fr.log2FoldChange.min() - 0.2, fr.log2FoldChange.max() + 0.55)
    save(fig, "Fig6_compartments")


def figS1_compartment_genes():
    p = f"{RES}/s5c_compartment_de_all.csv"
    if not os.path.exists(p):
        return
    d = pd.read_csv(p)
    d = d[(d.padj < 0.05) & (d.log2FoldChange.abs() > 0.585)]
    keep = d.compartment.isin(["Myeloid", "Endothelial", "Mast", "T_cell",
                               "Fib_Myofibroblast", "Fib_Pericyte", "Fib_SMC_contractile"])
    d = d[keep].sort_values("log2FoldChange")
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    cols = [FAIL if v > 0 else MAT for v in d.log2FoldChange]
    ax.barh(np.arange(len(d)), d.log2FoldChange, color=cols, height=0.7)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels([f"{g} ({c.replace('Fib_','')})" for g, c in
                        zip(d.gene, d.compartment)], fontsize=5.2)
    ax.axvline(0, c="k", lw=0.6)
    ax.set_xlabel("log2FC (failed vs matured)")
    ax.set_title("All compartment-resolved DE genes")
    save(fig, "FigS1_compartment_genes")


if __name__ == "__main__":
    a = sc.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    fig1(a)
    del a
    sub = sc.read_h5ad(os.path.join(WORK, "05_traj.h5ad"))
    pb = None
    pbp = f"{RES}/s5_pseudobulk_de.csv"
    if os.path.exists(pbp):
        pb = pd.read_csv(pbp)
    fig2(sub, pb)
    fig3(sub)
    fig4(sub)
    fig5()
    fig6()
    figS1_compartment_genes()
    print("all figures done")
