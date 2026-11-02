"""导出 PAGA 连接矩阵/布局与拟时序，供绘图脚本使用。
注：本对象的 uns['paga']['groups'] 存的是分组键名（'subtype'），
行/列顺序由 uns['subtype_sizes'] 对应的分类顺序给出。
"""
import anndata as ad
import numpy as np
import pandas as pd

RES = r"D:/AVF_GSE_2611/results"
a = ad.read_h5ad(r"D:/AVF_GSE_2611/data/processed/05_traj.h5ad", backed="r")
paga = a.uns["paga"]
conn = np.asarray(paga["connectivities"].toarray(), dtype=float)
pos = np.asarray(paga["pos"], dtype=float)
cats = list(a.obs["subtype"].cat.categories)
sizes = np.asarray(a.uns["subtype_sizes"])
assert len(cats) == conn.shape[0], (len(cats), conn.shape)
print(dict(zip(cats, sizes)))

df = pd.DataFrame(conn, index=cats, columns=cats)
df.to_csv(f"{RES}/plot_paga.csv")
pd.DataFrame(pos, index=cats, columns=["x", "y"]).to_csv(f"{RES}/plot_paga_pos.csv")
print(df.round(3).to_string())
print("size order matches categories:", sizes.tolist())

obs = a.obs[["subtype", "dpt_pseudotime", "group", "donor"]].copy()
xy = a.obsm["X_umap"]
obs["U1"], obs["U2"] = xy[:, 0], xy[:, 1]
obs.to_csv(f"{RES}/plot_traj_umap.csv", index=False)
print("saved", obs.shape, "dpt", round(obs.dpt_pseudotime.min(), 3), round(obs.dpt_pseudotime.max(), 3))
