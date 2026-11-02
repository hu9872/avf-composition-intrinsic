"""S8b（补充）: LIANA rank_aggregate 分组独立运行，用于主通讯结论的稳健性佐证。

设计：matured 与 failed 两组各自独立建对象（每 donor x state 最多 N_SUB 个细胞），
分别跑 rank_aggregate，再比较两组 rank 靠前的互作，避免"把两组混一起算rank"的混杂。
"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
N_SUB = 120
GROUPS = ["mature_AVF", "failed_AVF"]


def subsample(a, n=N_SUB, seed=SEED):
    rng = np.random.default_rng(seed)
    keep = []
    idx_by = a.obs.groupby(["donor", "state"], observed=True).indices
    for k, idx in idx_by.items():
        if len(idx) > n:
            keep.extend(rng.choice(idx, n, replace=False).tolist())
        else:
            keep.extend(idx.tolist())
    keep = sorted(set(keep))
    return a[keep].copy()


def main():
    import liana as li
    a = ad.read_h5ad(os.path.join(WORK, "06_states_all.h5ad"))
    out = {}
    for g in GROUPS:
        b = a[a.obs.group == g].copy()
        b = subsample(b)
        b.obs["state"] = b.obs["state"].astype("category")
        b.obs["state"] = b.obs["state"].cat.remove_unused_categories()
        li.mt.rank_aggregate(b, groupby="state", resource_name="consensus",
                             expr_prop=0.1, use_raw=False, verbose=False,
                             n_perms=100, seed=SEED, key_added="liana_res")
        r = b.uns["liana_res"].copy()
        r["group"] = g
        out[g] = r
        r.to_csv(f"{RES}/s8b_liana_{g}.csv", index=False)
        print(g, "cells:", b.n_obs, "interactions:", r.shape[0])

    if len(out) == 2:
        key = ["source", "target", "ligand_complex", "receptor_complex"]
        f = out["failed_AVF"].sort_values("magnitude_rank").head(300)
        m = out["mature_AVF"].sort_values("magnitude_rank").head(300)
        fs = set(map(tuple, f[key].astype(str).values))
        ms = set(map(tuple, m[key].astype(str).values))
        only_f = fs - ms
        only_m = ms - fs
        pd.DataFrame(sorted(only_f), columns=key).to_csv(f"{RES}/s8b_only_failed.csv", index=False)
        pd.DataFrame(sorted(only_m), columns=key).to_csv(f"{RES}/s8b_only_mature.csv", index=False)
        print(f"top300 overlap: shared={len(fs & ms)}, failed-only={len(only_f)}, mature-only={len(only_m)}")
        json.dump(dict(shared=len(fs & ms), failed_only=len(only_f), mature_only=len(only_m),
                       n_sub_per_state=N_SUB),
                  open(f"{RES}/s8b_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
