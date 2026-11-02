"""S4: 成纤维谱系（成纤维/周细胞/SMC）再聚类、程序打分与亚群定义。"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.dirname(__file__))
from programs import PROGRAMS, score

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
LINEAGE = ["Fibroblast", "Pericyte", "SMC"]

SUB_MARKERS = {
    # 外膜驻留型（adventitial resident）：PI16/ABCA8 高，低收缩、低基质分泌
    "AF_resident": ["PI16", "ABCA8", "MFAP5", "GSN", "TPPP3", "HTRA3", "SERPINF1"],
    # 基质型成纤维：经典胶原/蛋白聚糖高
    "AF_stromal": ["DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "BGN", "ASPN", "C7", "CFD"],
    # 肌成纤维：收缩应力 + ECM + 组织重塑
    "Myofibroblast": ["POSTN", "COL8A1", "FN1", "CTHRC1", "TGFBI", "SPP1", "COMP",
                      "TIMP1", "COL5A1", "COL10A1"],
    # 周细胞样
    "Pericyte": ["RGS5", "KCNJ8", "ABCC9", "NOTCH3", "HIGD1B", "COX4I2"],
    # 收缩型 SMC
    "SMC_contractile": ["MYH11", "CNN1", "MYOCD", "LMOD1", "DES", "ACTG2"],
    # 合成型/去分化 SMC
    "SMC_synthetic": ["TAGLN", "MYL9", "TPM2", "ELN", "COL4A1", "COL4A2"],
}


def main():
    adata = ad.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    sub = adata[adata.obs["celltype"].isin(LINEAGE)].copy()
    sub.X = sub.layers["lognorm"].copy()
    print("lineage cells:", sub.shape)
    print(sub.obs["celltype"].value_counts().to_string())

    sc.pp.highly_variable_genes(sub, n_top_genes=1500, batch_key="donor", flavor="seurat")
    hv = sub[:, sub.var.highly_variable].copy()
    sc.pp.pca(hv, n_comps=40, zero_center=False, random_state=SEED)
    try:
        import harmonypy  # noqa: F401
        sc.external.pp.harmony_integrate(hv, "donor", random_state=SEED, max_iter_harmony=20)
        method = "harmony"
    except Exception as e:
        raise RuntimeError("Harmony 失败(%s)，不使用 ComBat 兜底" % type(e).__name__)
    sub.obsm["X_pca"] = hv.obsm["X_pca_harmony"]
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.leiden(sub, resolution=0.6, key_added="sub_0.6", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.leiden(sub, resolution=1.0, key_added="sub_1.0", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(sub, random_state=SEED)

    rep = score(sub, PROGRAMS, prefix="pr_")
    print("score coverage:", rep)
    sc.tl.rank_genes_groups(sub, "sub_0.6", method="wilcoxon", use_raw=False)
    sc.get.rank_genes_groups_df(sub, None).to_csv(f"{RES}/s4_subcluster_markers.csv", index=False)

    # 亚群注释：marker 平均表达（z-score 后按群取均值）+ 程序得分联合判定
    keys = "sub_0.6"
    genes_all = sorted({g for v in SUB_MARKERS.values() for g in v if g in sub.var_names})
    mat = sub[:, genes_all].layers["lognorm"]
    dense = pd.DataFrame(np.asarray(mat.todense(), dtype=np.float32),
                         index=sub.obs_names, columns=genes_all)
    z = (dense - dense.mean()) / (dense.std().replace(0, np.nan))
    z = z.fillna(0.0)
    z["__cl"] = sub.obs[keys].values
    cl_mean = z.groupby("__cl", observed=True).mean()

    pr_cols = [c for c in sub.obs.columns if c.startswith("pr_")]
    pr_mean = sub.obs.groupby(keys, observed=True)[pr_cols].mean()

    sup_df = sc.get.rank_genes_groups_df(sub, None)
    ann = {}
    for cl in sub.obs[keys].cat.categories:
        if cl not in cl_mean.index:
            continue
        mk = {k: float(np.nanmean([cl_mean.loc[cl, g] for g in v if g in cl_mean.columns]))
              for k, v in SUB_MARKERS.items()}
        order = sorted(mk, key=mk.get, reverse=True)
        best, margin = order[0], mk[order[0]] - mk[order[1]]
        # 秩和标记支持度交叉核验（与主群注释同一逻辑）
        t = sup_df[(sup_df["group"] == cl) & (sup_df.logfoldchanges > 0.5) &
                   (sup_df.pvals_adj < 0.05)].head(50)
        nm = set(t["names"].astype(str))
        support = {k: len(nm & set(v)) for k, v in SUB_MARKERS.items()}
        pr = pr_mean.loc[cl]
        # 谱系内的每个细胞必然属于某一壁层/成纤维状态，"Unknown" 在此无生物学意义；
        # 且占比大的状态（如肌成纤维）在跨 cluster 标准化中 z 会被压缩，
        # 因此不以绝对阈值设 Unknown，而按 argmax 赋值，并记录 margin 供审阅。
        ann[cl] = dict(label=best,
                       n=int((sub.obs[keys] == cl).sum()),
                       margin=round(float(margin), 3),
                       support_best=int(support[best]),
                       support_cross=",".join(f"{k}:{v}" for k, v in
                                              sorted(support.items(), key=lambda x: -x[1])[:3]),
                       **{k: round(v, 3) for k, v in mk.items()},
                       **{k: round(float(v), 3) for k, v in pr.items()})
    pd.DataFrame(ann).T.to_csv(f"{RES}/s4_subcluster_annotation.csv")
    print(pd.DataFrame(ann).T[["label", "n", "margin", "support_best"]].to_string())

    sub.obs["subtype"] = sub.obs[keys].map({k: v["label"] for k, v in ann.items()}).astype(str)
    sub.obs["subtype"] = pd.Categorical(sub.obs["subtype"])
    print(sub.obs["subtype"].value_counts().to_string())

    sub.uns["batch_method"] = method
    sub.write_h5ad(os.path.join(WORK, "04_fibro.h5ad"))
    json.dump(ann, open(f"{RES}/s4_subcluster_annotation.json", "w"), indent=2)
    print("saved 04_fibro.h5ad")


if __name__ == "__main__":
    main()
