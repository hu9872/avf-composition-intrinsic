"""S2: 质控、标准化、去批次、聚类。内存友好（稀疏 PCA，不做 dense scaling）。"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1
WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
os.makedirs(RES, exist_ok=True)
SEED = 2611

MIN_GENES, MIN_COUNTS, MAX_MT, MAX_GENES = 200, 500, 15.0, 6000


def main():
    adata = ad.read_h5ad(os.path.join(WORK, "01_raw_concat.h5ad"))
    print("raw:", adata.shape)

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    qc0 = adata.obs[["sample", "group", "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    qc0.to_csv(os.path.join(RES, "s2_qc_before.csv"))

    keep = ((adata.obs.n_genes_by_counts >= MIN_GENES) &
            (adata.obs.total_counts >= MIN_COUNTS) &
            (adata.obs.pct_counts_mt < MAX_MT) &
            (adata.obs.n_genes_by_counts < MAX_GENES)).values
    print("filter pass: %d / %d (%.1f%%)" % (keep.sum(), len(keep), 100 * keep.mean()))
    adata = adata[keep].copy()
    adata.layers["counts"] = adata.X.copy()  # 保留原始计数，供 pseudobulk 使用

    # 双细胞检测（scanpy 内置 scrublet，按样本分别运行）
    # 不再静默跳过：去双细胞是审稿必查项，失败必须暴露
    sc.pp.scrublet(adata, batch_key="sample", random_state=SEED, verbose=False)
    n_dbl = int(adata.obs["predicted_doublet"].sum())
    print("predicted doublets: %d (%.2f%%)" % (n_dbl, 100 * n_dbl / adata.n_obs))
    adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
    print("after doublet removal:", adata.shape)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="donor", flavor="seurat")
    hv = adata[:, adata.var.highly_variable].copy()
    sc.pp.pca(hv, n_comps=50, zero_center=False, random_state=SEED)

    # 去批次：以供体为批次单元做 Harmony。
    # 注意：不允许退回 ComBat——97k x 36.6k 稠密化约 14 GB，必然 OOM。
    try:
        import harmonypy  # noqa: F401
        sc.external.pp.harmony_integrate(hv, "donor", random_state=SEED, max_iter_harmony=20)
        method = "harmony"
    except Exception as e:
        raise RuntimeError(
            "Harmony 去批次失败(%s: %s)。请勿改用 ComBat（内存不可行）。"
            % (type(e).__name__, str(e)[:120])
        )

    adata.obsm["X_pca"] = hv.obsm["X_pca_harmony"]
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, random_state=SEED)
    for r in (0.4, 0.8, 1.2):
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=SEED)

    adata.uns["batch_method"] = method
    adata.write_h5ad(os.path.join(WORK, "02_qc.h5ad"))
    json.dump({"batch_method": method, "n_cells": int(adata.n_obs),
               "n_genes": int(adata.n_vars),
               "params": dict(MIN_GENES=MIN_GENES, MIN_COUNTS=MIN_COUNTS,
                              MAX_MT=MAX_MT, MAX_GENES=MAX_GENES)},
              open(os.path.join(RES, "s2_qc_summary.json"), "w"), indent=2)

    per = adata.obs.groupby("sample").size().rename("cells")
    per.to_csv(os.path.join(RES, "s2_cells_per_sample.csv"))
    print("done:", adata.shape, "batch:", method)
    print(per.to_string())


if __name__ == "__main__":
    main()
