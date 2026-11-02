"""S10: 配对 bulk 验证（GSE220796 主验证；GSE119296 复现重塑签名）。"""
import os, re, sys, json
import numpy as np
import pandas as pd
import urllib.request, socket

sys.path.insert(0, os.path.dirname(__file__))
from programs import PROGRAMS, MATRISOME_CORE, SYNTHETIC_MYOF, ADVENTITIAL, MEDIA

socket.setdefaulttimeout(40)
RAW = r"D:/AVF_GSE_2611/data/raw"
RES = r"D:/AVF_GSE_2611/results"
WORK = r"D:/AVF_GSE_2611/data/processed"
os.makedirs(RES, exist_ok=True)


def geo_titles(acc):
    import time
    cache = os.path.join(RES, f"meta_{acc}.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache)
    last = None
    for _ in range(6):
        try:
            h = urllib.request.urlopen(
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief",
                timeout=60).read().decode("utf8", "ignore")
            gsm, rows = None, []
            for line in h.split("\n"):
                if line.startswith("^SAMPLE"):
                    gsm = line.split("=")[1].strip()
                elif line.startswith("!Sample_title"):
                    rows.append((gsm, line.split("=", 1)[1].strip()))
            df = pd.DataFrame(rows, columns=["gsm", "title"])
            df.to_csv(cache, index=False)
            return df
        except Exception as e:
            last = e
            time.sleep(4)
    raise RuntimeError(f"GEO metadata 获取失败: {last}")


def load_counts(fn):
    d = pd.read_csv(os.path.join(RAW, fn), sep="\t", index_col=0)
    d.index = [str(x).strip('"') for x in d.index]
    d = d.groupby(level=0).sum()
    return d


def zscore_score(mat, genes):
    g = [x for x in genes if x in mat.index]
    if len(g) < 3:
        return None, 0
    z = mat.loc[g].apply(lambda r: (r - r.mean()) / (r.std() + 1e-9), axis=1)
    return z.mean(0), len(g)


def main():
    meta = geo_titles("GSE220796")
    pat = meta["title"].str.extract(r"Patient\s*(\d+)\s*,\s*([^,]+),\s*(Matured|Failed),\s*(\S+)")
    meta[["patient", "tissue", "outcome", "code"]] = pat
    meta.to_csv(f"{RES}/s10_meta_GSE220796.csv", index=False)

    cnt = load_counts("GSE220796_rVazquez-Padron_counts.txt.gz")
    code2meta = meta.set_index("code")
    common = [c for c in cnt.columns if c in code2meta.index]
    cnt = cnt[common]
    print("samples matched:", len(common), "/", cnt.shape[1])

    avf = code2meta.loc[common].reset_index().rename(columns={"index": "code"})
    avf = avf[avf.tissue.str.contains("Juxta")].copy()
    # 仅保留有配对术前静脉的患者，与原文配对设计一致（排除孤立的 Patient 17）
    paired_pat = code2meta.loc[common].reset_index().rename(columns={"index": "code"})
    paired_pat = paired_pat.groupby("patient").tissue.nunique()
    avf = avf[avf.patient.isin(paired_pat[paired_pat == 2].index)].copy()
    avf = avf.set_index("code")
    X = cnt[avf.index]
    print("AVF samples:", X.shape, avf.outcome.value_counts().to_dict())

    # 1) DESeq2：failed vs matured（AVF 样本）
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    counts = X.T.round().astype(int)
    counts = counts.loc[:, counts.sum(0) > 10]
    coldata = pd.DataFrame({"outcome": avf.outcome.values}, index=avf.index)
    dds = DeseqDataSet(counts=counts, metadata=coldata, design="~outcome", refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["outcome", "Failed", "Matured"], quiet=True)
    st.summary()
    res = st.results_df.copy()
    res.to_csv(f"{RES}/s10_deseq_AVF_failed_vs_matured.csv")
    print("DE genes FDR<0.05 & |log2FC|>0.585:",
          int(((res.padj < 0.05) & (res.log2FoldChange.abs() > 0.585)).sum()))

    # 2) 程序评分比较
    rows = []
    prog_list = dict(PROGRAMS)
    prog_list.setdefault("Media", MEDIA)
    for name, genes in prog_list.items():
        s, ng = zscore_score(X, genes)
        if s is None:
            continue
        for smp, val in s.items():
            rows.append(dict(sample=smp, program=name, score=val,
                             outcome=avf.loc[smp, "outcome"], patient=avf.loc[smp, "patient"]))
    sc_df = pd.DataFrame(rows)
    sc_df.to_csv(f"{RES}/s10_program_scores.csv", index=False)
    from scipy import stats
    summ = []
    for name, g in sc_df.groupby("program"):
        a = g[g.outcome == "Failed"].score.values
        b = g[g.outcome == "Matured"].score.values
        t, p = stats.ttest_ind(a, b, equal_var=False)
        summ.append(dict(program=name, failed=a.mean(), matured=b.mean(),
                         diff=a.mean() - b.mean(), p=p,
                         n_failed=len(a), n_matured=len(b)))
    summ = pd.DataFrame(summ).sort_values("p")
    summ.to_csv(f"{RES}/s10_program_tests.csv", index=False)
    print(summ.to_string(index=False))

    # 3) 配对设计：AVF vs 术前静脉（确认重塑签名，非结局比较）
    both = code2meta.loc[common].reset_index().rename(columns={"index": "code"})
    pairs = both.pivot_table(index="patient", columns="tissue", values="code", aggfunc="first")
    print("paired patients (vein+AVF):", pairs.dropna().shape[0])

    # 4) 关键基因方向核对（单细胞结论的外部检验）
    focus = ["POSTN", "COL8A1", "FN1", "COL1A1", "COL3A1", "PI16", "MFAP5", "SFRP2",
             "ACTA2", "MYH11", "LOX", "LOXL1", "TIMP1", "DCN", "LUM"]
    fr = res.reindex([g for g in focus if g in res.index])[
        ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].round(4)
    fr.to_csv(f"{RES}/s10_focus_genes.csv")
    print("\n--- 关键基因 bulk 方向 ---")
    print(fr.to_string())

    json.dump(dict(n_avf=int(X.shape[1]),
                   n_failed=int((avf.outcome == "Failed").sum()),
                   n_matured=int((avf.outcome == "Matured").sum()),
                   n_pairs=int(pairs.dropna().shape[0])),
              open(f"{RES}/s10_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
