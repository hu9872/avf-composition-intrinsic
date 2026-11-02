"""S3: 主群注释（跨 cluster 相对富集打分 + 秩和标记一致性核验 + 混合群/污染群剔除）。

为什么不用绝对打分：组织解离样本存在大量环境 RNA（ambient collagen / hemoglobin），
绝对表达量高的群（尤其是髓系）会被误判为成纤维细胞。改为对每个基因先做
"跨 cluster 标准化"，再按面板取均值，衡量的是相对富集而非绝对丰度。
"""
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
os.makedirs(RES, exist_ok=True)
SEED = 2611
RES_KEY = "leiden_0.8"
MARGIN_MIN = 0.25   # top1 - top2 的 z 差值下限，低于此判为 Mixed

PANELS = {
    "Endothelial": ["PECAM1", "VWF", "CDH5", "EGFL7", "PLVAP", "KDR", "CLDN5", "AQP1",
                    "ENG", "ACKR1", "RAMP2", "ADGRL4", "CD34", "ADGRF5", "BMX", "ROBO4"],
    "Lymphatic_EC": ["PROX1", "LYVE1", "CCL21", "PDPN", "MMRN1", "FLT4", "TFF3"],
    "SMC": ["MYH11", "CNN1", "ACTA2", "TAGLN", "MYOCD", "LMOD1", "ACTG2", "MYL9",
            "TPM2", "DES", "CALD1", "PLN"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
                   "COL6A3", "COL12A1", "COL14A1", "DCN", "LUM", "BGN", "ASPN", "POSTN",
                   "CCDC80", "FBLN1", "MFAP5", "PI16", "PDGFRA", "THY1", "SFRP2", "CFD",
                   "C7", "CFH", "ABI3BP", "CD248", "C1R", "C1S"],
    "Pericyte": ["RGS5", "KCNJ8", "ABCC9", "NOTCH3", "HIGD1B", "PDGFRB", "COX4I2"],
    "Myeloid": ["CD68", "CD14", "LYZ", "CSF1R", "AIF1", "CD163", "MRC1", "C1QA", "C1QB",
                "C1QC", "ITGAM", "APOE", "CCL18", "CD1C", "CD1E", "LILRB1", "FCGR3A",
                "FCN1", "S100A8", "S100A9"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "IL7R", "TRAC", "CD2", "BCL11B", "CD52", "CCR7",
               "LTB", "TRBC1", "TRBC2"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1", "CST7", "KLRF1", "GZMB", "GZMA"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "CD37", "BANK1"],
    "Plasma": ["JCHAIN", "MZB1", "IGHG1", "DERL3", "TNFRSF17", "IGKC"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "CTSG", "HPGD", "GATA2"],
    "Adipocyte": ["ADIPOQ", "PLIN1", "LEP", "PLIN4", "CIDEC"],
    "Neural": ["S100B", "PLP1", "MPZ", "SOX10"],
    "Erythrocyte": ["HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "AHSP", "CA1", "BNIP3L"],
}
N_TOP_MARKERS = 50
MIN_LFC = 0.5
MIN_SUPPORT = 3
HB_MAX = 1.0   # 血红蛋白均值阈值：正常群 0.04–0.36，污染群 >3，间隔极大


def cluster_score(sub, keys, panels, min_mean=0.05, min_frac_clusters=0.20):
    """返回 (z 分数表, cluster x gene 均值表)。

    只保留在多数 cluster 中可检出的基因：近乎全零的基因（如脂肪细胞 panel 在
    血管组织中）标准差趋近 0，会把噪声放大成极高的 z，导致误判。
    """
    genes = sorted({g for v in panels.values() for g in v if g in sub.var_names})
    X = sub[:, genes].X
    X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float32)
    df = pd.DataFrame(X, columns=genes)
    df["__cl"] = sub.obs[keys].values
    cm = df.groupby("__cl", observed=True).mean()
    detected = ((cm > 0.1).mean(0) >= min_frac_clusters) & (cm.mean(0) >= min_mean)
    keep = [g for g in genes if detected.get(g, False)]
    cm_k = cm[keep]
    z = (cm_k - cm_k.mean()) / (cm_k.std() + 1e-6)
    out = {}
    for name, gl in panels.items():
        g = [x for x in gl if x in z.columns]
        if len(g) >= 2:
            out[name] = z[g].mean(1)
    return pd.DataFrame(out), cm


def marker_support(top, cl, panels, n_top=N_TOP_MARKERS, min_lfc=MIN_LFC, min_padj=0.05):
    """以 Wilcoxon 秩和检验的前 N 个标记基因落在各 panel 的个数（按面板大小归一）。"""
    t = top[top["group"] == cl]
    t = t[(t.logfoldchanges > min_lfc) & (t.pvals_adj < min_padj)].head(n_top)
    names = list(t["names"].astype(str))
    sup = {name: len(set(names) & set(gl)) for name, gl in panels.items()}
    norm = {k: v / np.sqrt(len(panels[k])) for k, v in sup.items()}
    return sup, norm, names


def main():
    adata = ad.read_h5ad(os.path.join(WORK, "02_qc.h5ad"))
    print("cells:", adata.shape, "| clusters:", adata.obs[RES_KEY].nunique())

    sc.tl.rank_genes_groups(adata, RES_KEY, method="wilcoxon", use_raw=False, pts=True)
    top = sc.get.rank_genes_groups_df(adata, None)
    top.to_csv(os.path.join(RES, "s3_cluster_markers.csv"), index=False)

    zsc, cm = cluster_score(adata, RES_KEY, PANELS)
    zsc.to_csv(os.path.join(RES, "s3_cluster_zscores.csv"))
    cm.to_csv(os.path.join(RES, "s3_cluster_meanexpr.csv"))

    labels, flags, rows = {}, {}, []
    # 红细胞污染：以血红蛋白基因实测表达判定（比 marker panel 更直接）
    hb = [g for g in ["HBB", "HBA1", "HBA2"] if g in adata.var_names]
    Xh = adata[:, hb].X
    Xh = np.asarray(Xh.todense() if hasattr(Xh, "todense") else Xh, dtype=np.float32)
    hb_mean = pd.Series(Xh.mean(1), index=adata.obs_names).groupby(
        adata.obs[RES_KEY].values, observed=True).mean()
    hb_mean.to_csv(os.path.join(RES, "s3_cluster_haemoglobin.csv"))

    for cl in zsc.index:
        s = zsc.loc[cl].sort_values(ascending=False)
        margin = float(s.iloc[0] - s.iloc[1])
        sup, norm, names = marker_support(top, cl, PANELS)
        ns = pd.Series(norm).sort_values(ascending=False)
        best = ns.index[0]
        n_best = int(sup[best])
        # 判定规则（写入 Methods）：
        #   0) 血红蛋白均值 > HB_MAX -> Erythrocyte_contam（红细胞/环境血红蛋白污染）
        #   主判据 = 该 cluster 相对其余细胞的秩和标记基因落在各 panel 的个数
        #            （除以 sqrt(面板长度) 归一，避免大面板占优）
        #   1) 最佳 panel 命中 < MIN_SUPPORT -> Mixed（低置信，不参与谱系分析）
        #   2) 否则采用该 panel
        if float(hb_mean.get(cl, 0)) > HB_MAX:
            lab = "Erythrocyte_contam"
        elif n_best < MIN_SUPPORT:
            lab = "Mixed"
        elif best == "Erythrocyte":
            lab = "Erythrocyte_contam"
        else:
            lab = best
        labels[cl] = lab
        flags[cl] = dict(assigned=lab, marker_support=dict(pd.Series(sup).sort_values(
                             ascending=False).head(5).astype(int)),
                         best_panel=best, n_best_markers=n_best,
                         hb_mean=round(float(hb_mean.get(cl, 0)), 3),
                         z_margin=round(margin, 3),
                         z_top1=f"{s.index[0]}:{s.iloc[0]:.2f}",
                         z_top2=f"{s.index[1]}:{s.iloc[1]:.2f}",
                         top_markers=",".join(names[:12]),
                         n=int((adata.obs[RES_KEY] == cl).sum()))
        rows.append(f"cluster {cl:>2}: {lab:16s} support={n_best:2d}/{best} "
                    f"hb={float(hb_mean.get(cl,0)):5.2f} "
                    f"zmargin={margin:5.2f} ({s.index[0]} {s.iloc[0]:.2f}/{s.index[1]} {s.iloc[1]:.2f})  "
                    f"top={','.join(names[:6])}")

    open(os.path.join(RES, "s3_annotation_check.txt"), "w").write("\n".join(rows))
    print("\n".join(rows))
    pd.DataFrame(flags).T.to_csv(os.path.join(RES, "s3_cluster_decision.csv"))

    adata.obs["celltype"] = adata.obs[RES_KEY].map(labels).astype(str)
    adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"])
    print("\n--- 注释分布（含污染/混合） ---")
    print(adata.obs["celltype"].value_counts().to_string())

    # 剔除污染群与混合群，仅保留可用于谱系分析的高置信群
    BAD = {"Erythrocyte_contam", "Mixed"}
    n0 = adata.n_obs
    adata = adata[~adata.obs["celltype"].isin(BAD)].copy()
    print(f"\nremoved low-confidence clusters: {n0 - adata.n_obs} cells "
          f"({100*(n0-adata.n_obs)/n0:.1f}%); remain {adata.n_obs}")
    adata.obs["celltype"] = adata.obs["celltype"].cat.remove_unused_categories()

    ct = (adata.obs.groupby(["donor", "group", "celltype"], observed=True)
          .size().unstack(fill_value=0))
    ct.to_csv(os.path.join(RES, "s3_counts_donor_celltype.csv"))

    adata.write_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    json.dump({"labels": labels, "decision": flags,
               "rule": f"rank-based marker support, top{N_TOP_MARKERS} genes, "
                       f"logFC>{MIN_LFC}, support>={MIN_SUPPORT}, size-normalised",
               "excluded": sorted(BAD)},
              open(os.path.join(RES, "s3_cluster_labels.json"), "w"), indent=2,
              ensure_ascii=False, default=int)
    print("saved 03_annotated.h5ad")


if __name__ == "__main__":
    main()
