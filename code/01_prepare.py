"""S1: 解包 GSE305423，逐样本读入 10x 三件套，构建统一 AnnData 与元数据台账。"""
import os, re, tarfile, gzip, shutil
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

RAW = r"D:/AVF_GSE_2611/data/raw"
WORK = r"D:/AVF_GSE_2611/data/processed"
os.makedirs(WORK, exist_ok=True)

# 已核对的 GSM -> 分组映射（依据 GEO 样本标题，2026-10-31 核对）
GROUP = {}
for i, g in enumerate(range(9176931, 9176937), 1):
    GROUP[f"GSM{g}"] = ("pre_access", f"pre_access_{i}", f"P{i}")
GROUP["GSM9176937"] = ("early_AVF", "early_AVF_1A", "E1")
GROUP["GSM9176938"] = ("early_AVF", "early_AVF_1B", "E1")   # 与 1A 同一患者
GROUP["GSM9176939"] = ("early_AVF", "early_AVF_2", "E2")
for i, g in enumerate(range(9176940, 9176946), 1):
    GROUP[f"GSM{g}"] = ("mature_AVF", f"mature_AVF_{i}", f"M{i}")
for i, g in enumerate(range(9176946, 9176952), 1):
    GROUP[f"GSM{g}"] = ("failed_AVF", f"failed_AVF_{i}", f"F{i}")

SEX = {  # 来自 GEO characteristics 字段
    "GSM9176931": ("Male", 71), "GSM9176932": ("Male", 56), "GSM9176933": ("Female", 69),
    "GSM9176934": ("Female", 72), "GSM9176935": ("Female", 62), "GSM9176936": ("Female", 59),
    "GSM9176937": ("Female", 64), "GSM9176938": ("Female", 64), "GSM9176939": ("Male", 59),
    "GSM9176940": ("Male", 40), "GSM9176941": ("Female", 28), "GSM9176942": ("Female", 64),
    "GSM9176943": ("Female", 68), "GSM9176944": ("Female", 59), "GSM9176945": ("Male", 35),
    "GSM9176946": ("Female", 42), "GSM9176947": ("Male", 49), "GSM9176948": ("Female", 65),
    "GSM9176949": ("Female", 56), "GSM9176950": ("Female", 55), "GSM9176951": ("Female", 47),
}


def ensure_extracted():
    tar = os.path.join(RAW, "GSE305423_RAW.tar")
    d = os.path.join(RAW, "GSE305423")
    if os.path.isdir(d) and len(os.listdir(d)) >= 63:
        return d
    os.makedirs(d, exist_ok=True)
    with tarfile.open(tar) as t:
        t.extractall(d)
    return d


def read_sample(d, gsm):
    """按 GSM 前缀聚合 barcodes/features/matrix。"""
    files = [f for f in os.listdir(d) if f.startswith(gsm)]
    b = [f for f in files if "barcodes" in f][0]
    f_ = [f for f in files if "features" in f][0]
    m = [f for f in files if "matrix" in f][0]
    M = sio.mmread(os.path.join(d, m)).T.tocsr()  # 基因 × 细胞 -> 细胞 × 基因
    with gzip.open(os.path.join(d, b), "rt") as fh:
        bc = [x.strip() for x in fh]
    feats = pd.read_csv(os.path.join(d, f_), sep="\t", header=None)
    genes = feats[1].astype(str).values if feats.shape[1] > 1 else feats[0].astype(str).values
    a = ad.AnnData(X=M.astype(np.float32))
    a.obs_names = [f"{gsm}_{x}" for x in bc]
    a.var_names = genes
    a.var_names_make_unique()
    return a


if __name__ == "__main__":
    d = ensure_extracted()
    print("extracted files:", len(os.listdir(d)))
    adatas, rows = [], []
    for gsm, (grp, label, donor) in sorted(GROUP.items()):
        a = read_sample(d, gsm)
        a.obs["gsm"] = gsm
        a.obs["group"] = grp
        a.obs["sample"] = label
        a.obs["donor"] = donor
        s, age = SEX[gsm]
        a.obs["sex"] = s
        a.obs["age"] = age
        adatas.append(a)
        rows.append(dict(gsm=gsm, sample=label, group=grp, donor=donor,
                         sex=s, age=age, n_cells=a.n_obs, n_genes=a.n_vars))
        print(f"{gsm} {label:16s} {a.n_obs:6d} cells")
    meta = pd.DataFrame(rows)
    meta.to_csv(os.path.join(WORK, "sample_manifest.csv"), index=False)

    # 基因取交集，保证 var 一致
    common = set(adatas[0].var_names)
    for a in adatas[1:]:
        common &= set(a.var_names)
    common = sorted(common)
    print("common genes:", len(common))
    for a in adatas:
        a._inplace_subset_var([g in common for g in a.var_names])
    adata = ad.concat(adatas, join="inner", merge="same", index_unique=None)
    adata.obs = adata.obs.astype(str)
    adata.obs["age"] = adata.obs["age"].astype(int)
    adata.write_h5ad(os.path.join(WORK, "01_raw_concat.h5ad"))
    print("saved:", adata.shape)
