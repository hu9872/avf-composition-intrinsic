"""导出绘图所需的轻量表格（UMAP 坐标 + 元数据），避免绘图脚本反复读大 h5ad。"""
import anndata as ad
import pandas as pd

DATA = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"


def dump(path, cols, out, umap_key="X_umap"):
    a = ad.read_h5ad(path, backed="r")
    df = a.obs[cols].copy()
    xy = a.obsm[umap_key]
    df["U1"] = xy[:, 0]
    df["U2"] = xy[:, 1]
    df.to_csv(f"{RES}/{out}", index=False)
    print(out, df.shape)


dump(f"{DATA}/03_annotated.h5ad", ["celltype", "group", "donor"], "plot_atlas.csv")
dump(f"{DATA}/04_fibro.h5ad", ["subtype", "group", "donor"], "plot_fibro.csv")
dump(f"{DATA}/05d_myeloid.h5ad", ["mye_subtype", "group", "donor", "is_cDC2_cluster", "is_cdc2_free"],
     "plot_myeloid.csv")

# 六个状态的程序得分均值（供 Fig2B 热图）
a = ad.read_h5ad(f"{DATA}/04_fibro.h5ad", backed="r")
pcols = [c for c in a.obs.columns if c.startswith("pr_")]
pm = a.obs.groupby("subtype", observed=True)[pcols].mean()
pm.to_csv(f"{RES}/plot_fibro_programmes.csv")
print("plot_fibro_programmes.csv", pm.shape)
