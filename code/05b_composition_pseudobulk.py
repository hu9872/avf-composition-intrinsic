"""S5+S6: 供体水平组成检验 + pseudobulk 差异表达。"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
os.makedirs(RES, exist_ok=True)
SEED = 2611
GROUPS = ["mature_AVF", "failed_AVF"]


def clr(mat, pseudo=0.5):
    x = mat + pseudo
    lg = np.log(x)
    return lg.sub(lg.mean(axis=1), axis=0)


def composition(sub):
    post = sub.obs[sub.obs.group.isin(GROUPS)]
    counts = (post.groupby(["donor", "group", "subtype"], observed=True)
              .size().unstack(fill_value=0))
    frac = counts.div(counts.sum(1), axis=0)
    frac.to_csv(f"{RES}/s5_subtype_fraction.csv")
    cl = clr(counts.astype(float))

    rows = []
    for s in counts.columns:
        a = cl.loc[cl.index.get_level_values("group") == "failed_AVF", s].values
        b = cl.loc[cl.index.get_level_values("group") == "mature_AVF", s].values
        if len(a) < 3 or len(b) < 3:
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        u, pu = stats.mannwhitneyu(a, b)
        sp = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        d = (a.mean() - b.mean()) / sp if sp > 0 else np.nan
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        rows.append(dict(state=s, failed_clr=a.mean(), mature_clr=b.mean(),
                         diff=a.mean() - b.mean(), cohens_d=d, welch_p=p, mwu_p=pu,
                         ci_low=(a.mean() - b.mean()) - 1.96 * se,
                         ci_high=(a.mean() - b.mean()) + 1.96 * se,
                         n_failed=len(a), n_mature=len(b)))
    td = pd.DataFrame(rows).sort_values("welch_p")
    from statsmodels.stats.multitest import multipletests
    if len(td):
        td["fdr"] = multipletests(td.welch_p, method="fdr_bh")[1]
    td.to_csv(f"{RES}/s5_composition_tests.csv", index=False)
    print(td.round(3).to_string(index=False))
    return counts


def pseudobulk(sub, counts_meta):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    X = sub.layers["counts"].tocsr()
    obs = sub.obs
    post = obs[obs.group.isin(GROUPS)]
    pos = {name: i for i, name in enumerate(obs.index)}

    out = []
    for level, key in [("lineage", None), ("state", "subtype")]:
        if key is None:
            grp = post.assign(unit="lineage").groupby(["donor", "group", "unit"], observed=True).indices
        else:
            grp = post.groupby(["donor", "group", key], observed=True).indices
        recs, meta = [], []
        for (donor, group, unit), idx in grp.items():
            if len(idx) < 10:
                continue
            rows = np.array([pos[n] for n in post.index[idx]], dtype=np.int64)
            vec = np.asarray(X[rows].sum(0)).ravel()
            recs.append(vec); meta.append(dict(sample=f"{donor}|{unit}", donor=donor,
                                               group=group, unit=unit, n=len(idx)))
        if len(recs) < 6:
            continue
        mat = pd.DataFrame(np.vstack(recs), index=[m["sample"] for m in meta],
                           columns=sub.var_names)
        md = pd.DataFrame(meta).set_index("sample")
        for unit in md.unit.unique():
            sel = md[md.unit == unit].index
            m2 = md.loc[sel]
            if m2.group.value_counts().min() < 3:
                continue
            c = mat.loc[sel].round().astype(int)
            c = c.loc[:, c.sum(0) > 10]
            dds = DeseqDataSet(counts=c, metadata=m2[["group"]],
                               design="~group", refit_cooks=True, quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["group", "failed_AVF", "mature_AVF"], quiet=True)
            st.summary()
            r = st.results_df.copy()
            r["unit"] = unit; r["level"] = level
            r["n_failed"] = int((m2.group == "failed_AVF").sum())
            r["n_mature"] = int((m2.group == "mature_AVF").sum())
            out.append(r.reset_index().rename(columns={"index": "gene"}))
            n_sig = int(((r.padj < 0.05) & (r.log2FoldChange.abs() > 0.585)).sum())
            print(f"  [{level}] {unit}: {m2.shape[0]} samples, {n_sig} DE genes")
    if out:
        de = pd.concat(out, ignore_index=True)
        de.to_csv(f"{RES}/s5_pseudobulk_de.csv", index=False)
        print("total DE rows:", de.shape)
        return de
    return None


def main():
    sub = sc.read_h5ad(os.path.join(WORK, "04_fibro.h5ad"))
    composition(sub)
    pseudobulk(sub, None)


if __name__ == "__main__":
    main()
