"""S10b（正向对照）: GSE119296 术前静脉 (V) vs 术后 AVF (T)。

目的：验证我们的程序评分体系对"已知会发生的静脉->内瘘重塑"确实敏感
      （正性对照），而非仅对结局分组敏感。

诚实声明（写入 Methods）：该数据集样本编号在 V/T 两组间仅 1 例重合，
无法按患者配对，故按非配对两组比较处理；仅供方向性对照，不用于结局推断。
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from programs import PROGRAMS, MEDIA, ADVENTITIAL

RAW = r"D:/AVF_GSE_2611/data/raw"
RES = r"D:/AVF_GSE_2611/results"
FN = "GSE119296_Vazquez-Padron_RNAseq_080718_counts.txt.gz"


def load_hgnc():
    """Ensembl ID -> 基因符号（GSE119296 提供的是 Ensembl ID）。"""
    p = r"D:/AVF_GSE_2611/data/ref/hgnc_complete_set.txt"
    d = pd.read_csv(p, sep="\t", usecols=["symbol", "ensembl_gene_id"], low_memory=False)
    d = d.dropna(subset=["ensembl_gene_id"])
    return dict(zip(d.ensembl_gene_id.astype(str), d.symbol.astype(str)))


def to_symbols(df):
    m = load_hgnc()
    base = [x.split(".")[0] for x in df.index.astype(str)]
    sym = [m.get(b) for b in base]
    df = df.copy()
    df.index = [s if isinstance(s, str) else b for s, b in zip(sym, base)]
    df = df[~df.index.str.startswith("ENSG")]
    return df.groupby(level=0).sum()


def zscore_score(mat, genes):
    g = [x for x in genes if x in mat.index]
    if len(g) < 3:
        return None, 0
    z = mat.loc[g].apply(lambda r: (r - r.mean()) / (r.std() + 1e-9), axis=1)
    return z.mean(0), len(g)


def main():
    meta = pd.read_csv(f"{RES}/meta_GSE119296.csv")
    meta["tissue"] = np.where(meta.title.str.startswith("V"), "Vein", "AVF")

    cnt = pd.read_csv(os.path.join(RAW, FN), sep="\t", index_col=0)
    cnt.index = [str(x).strip('"') for x in cnt.index]
    cnt = cnt.groupby(level=0).sum()
    print("原始基因数 (Ensembl):", cnt.shape[0])
    cnt = to_symbols(cnt)
    print("映射后基因数 (symbol):", cnt.shape[0])

    code2 = meta.set_index("gsm")
    common = [c for c in cnt.columns if c in code2.index] or \
             [c for c in cnt.columns if c in meta.title.values]
    if common and common[0] in code2.index:
        X = cnt[common]
        tis = code2.loc[common, "tissue"]
    else:
        t2 = meta.set_index("title")
        X = cnt[[c for c in cnt.columns if c in t2.index]]
        tis = t2.loc[X.columns, "tissue"]
    print("samples:", X.shape, tis.value_counts().to_dict())

    # 非配对 DESeq2
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    counts = X.T.round().astype(int)
    counts = counts.loc[:, counts.sum(0) > 10]
    coldata = pd.DataFrame({"tissue": tis.values}, index=X.columns)
    dds = DeseqDataSet(counts=counts, metadata=coldata, design="~tissue",
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["tissue", "AVF", "Vein"], quiet=True)
    st.summary()
    res = st.results_df.copy()
    res.to_csv(f"{RES}/s10b_deseq_AVF_vs_vein.csv")
    n_sig = int(((res.padj < 0.05) & (res.log2FoldChange.abs() > 0.585)).sum())
    print("DE (FDR<0.05, |log2FC|>0.585):", n_sig)

    rows = []
    prog = dict(PROGRAMS); prog.setdefault("Media", MEDIA)
    for name, genes in prog.items():
        s, ng = zscore_score(X, genes)
        if s is None:
            continue
        for smp, val in s.items():
            rows.append(dict(sample=smp, program=name, score=val, tissue=tis.loc[smp]))
    sc = pd.DataFrame(rows)
    if sc.empty:
        raise RuntimeError("程序打分为空：请检查基因符号映射（Ensembl->symbol）")
    sc.to_csv(f"{RES}/s10b_program_scores.csv", index=False)

    out = []
    for name, g in sc.groupby("program"):
        a = g[g.tissue == "AVF"].score.values
        b = g[g.tissue == "Vein"].score.values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        out.append(dict(program=name, AVF=a.mean(), Vein=b.mean(),
                        diff=a.mean() - b.mean(), p=p, n_AVF=len(a), n_Vein=len(b)))
    o = pd.DataFrame(out).sort_values("p")
    o.to_csv(f"{RES}/s10b_program_tests.csv", index=False)
    print(o.round(3).to_string(index=False))

    focus = ["POSTN", "COL8A1", "FN1", "COL1A1", "COL3A1", "ACTA2", "MYH11", "PI16", "LOX"]
    fr = res.reindex([g for g in focus if g in res.index])[
        ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].round(4)
    fr.to_csv(f"{RES}/s10b_focus_genes.csv")
    print(fr.to_string())

    json.dump(dict(n_AVF=int((tis == "AVF").sum()), n_Vein=int((tis == "Vein").sum()),
                   n_de=n_sig, paired=False,
                   note="V/T 编号仅 1 例重合，按非配对处理"),
              open(f"{RES}/s10b_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
