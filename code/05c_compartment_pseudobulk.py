"""S5c: 主要细胞类型分别做供体水平 pseudobulk，定位"差异到底在哪个区室"。

动机：整组织 bulk（GSE220796）有 594 个差异基因，而成纤维谱系 pseudobulk 仅 12 个。
若差异主要落在免疫/内皮区室而非成纤维本身，则"内瘘失败=成纤维驱动"的通行假设需要修正。
这是本研究的决定性分析之一。
"""
import os, sys, json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
GROUPS = ["mature_AVF", "failed_AVF"]
MIN_CELLS = 10
MIN_DONORS = 3


def main():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    a = sc.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    sub = sc.read_h5ad(os.path.join(WORK, "04_fibro.h5ad"))

    # 成纤维细胞用亚型标签，其余用主群标签
    st = a.obs["celltype"].astype(str).copy()
    m = sub.obs["subtype"].astype(str)
    hit = a.obs_names.isin(m.index)
    st.loc[hit] = "Fib_" + m.reindex(a.obs_names[hit]).astype(str).values
    a.obs["compartment"] = st.values

    a = a[a.obs["group"].isin(GROUPS)].copy()
    post = a.obs
    X = a.layers["counts"].tocsr()
    pos = {n: i for i, n in enumerate(a.obs_names)}

    comps = [c for c in post["compartment"].unique()]
    rows, de_all = [], []
    for comp in sorted(comps):
        obs_c = post[post["compartment"] == comp]
        grp = obs_c.groupby(["donor", "group"], observed=True).indices
        recs, meta = [], []
        for (donor, group), idx in grp.items():
            if len(idx) < MIN_CELLS:
                continue
            r = np.array([pos[n] for n in obs_c.index[idx]], dtype=np.int64)
            recs.append(np.asarray(X[r].sum(0)).ravel())
            meta.append(dict(sample=f"{donor}", donor=donor, group=group, n=len(idx)))
        if len(recs) < 2 * MIN_DONORS:
            continue
        md = pd.DataFrame(meta).set_index("sample")
        if md.group.value_counts().min() < MIN_DONORS:
            rows.append(dict(compartment=comp, n_failed=int((md.group == "failed_AVF").sum()),
                             n_mature=int((md.group == "mature_AVF").sum()),
                             note="insufficient donors", n_de=0))
            continue
        c = pd.DataFrame(np.vstack(recs), index=md.index, columns=a.var_names)
        c = c.loc[:, c.sum(0) > 10].round().astype(int)
        cfg = md[["group"]].copy()
        cfg["group"] = pd.Categorical(cfg["group"], categories=["mature_AVF", "failed_AVF"])
        dds = DeseqDataSet(counts=c, metadata=cfg, design="~group",
                           refit_cooks=True, quiet=True)
        dds.deseq2()
        stt = DeseqStats(dds, contrast=["group", "failed_AVF", "mature_AVF"], quiet=True)
        stt.summary()
        r = stt.results_df.copy()
        sig = (r.padj < 0.05) & (r.log2FoldChange.abs() > 0.585)
        n_de = int(sig.sum())
        top = (r[sig].reindex(r[sig].log2FoldChange.abs().sort_values(ascending=False).index)
               .head(12))
        rows.append(dict(compartment=comp,
                         n_failed=int((md.group == "failed_AVF").sum()),
                         n_mature=int((md.group == "mature_AVF").sum()),
                         n_cells=int(md.n.sum()), n_tested=int(r.shape[0]),
                         n_de=n_de, n_up=int(((r.padj < 0.05) & (r.log2FoldChange > 0.585)).sum()),
                         n_down=int(((r.padj < 0.05) & (r.log2FoldChange < -0.585)).sum()),
                         top_genes=";".join(f"{g}({v:+.2f})" for g, v in
                                            top.log2FoldChange.items()),
                         note=""))
        rr = r.copy(); rr["compartment"] = comp
        de_all.append(rr.reset_index().rename(columns={"index": "gene"}))
        print(f"{comp:20s} {rows[-1]['n_failed']}v{rows[-1]['n_mature']} "
              f"cells={rows[-1]['n_cells']:6d} DE={n_de:5d}")

    td = pd.DataFrame(rows).sort_values("n_de", ascending=False)
    td.to_csv(f"{RES}/s5c_compartment_de_summary.csv", index=False)
    print("\n", td[["compartment", "n_failed", "n_mature", "n_cells", "n_de",
                     "n_up", "n_down"]].to_string(index=False))
    if de_all:
        pd.concat(de_all, ignore_index=True).to_csv(f"{RES}/s5c_compartment_de_all.csv",
                                                    index=False)
    json.dump(td.to_dict("records"), open(f"{RES}/s5c_summary.json", "w"),
              indent=2, default=int)
    print("\nsaved s5c_*")


if __name__ == "__main__":
    main()
