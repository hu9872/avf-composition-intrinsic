"""S9: ECM 程序（NMF）与壁层分区计算重建 + 邻域生态位。"""
import os, sys, json, urllib.request, socket
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import NMF

sys.path.insert(0, os.path.dirname(__file__))
from programs import (MATRISOME_CORE, COLLAGENS, GLYCOPROTEINS, PROTEOGLYCANS,
                      CROSSLINKERS, INTIMA, MEDIA, ADVENTITIA)

socket.setdefaulttimeout(45)
WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
SEED = 2611
os.makedirs(RES, exist_ok=True)


def enrichr_matrisome():
    """从 Enrichr 拉取 GO:0062023（collagen-containing ECM）基因集；失败则返回 None。"""
    try:
        u = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
             "&libraryName=GO_Cellular_Component_2023")
        txt = urllib.request.urlopen(u, timeout=90).read().decode("utf8", "ignore")
        for line in txt.split("\n"):
            f = line.rstrip("\n").split("\t")
            if f and "collagen-containing extracellular matrix" in f[0].lower():
                genes = [g.strip() for g in f[1:] if g.strip()]
                open(f"{RES}/s9_matrisome_geneset.txt", "w").write("\n".join(genes))
                return genes
    except Exception as e:
        print("Enrichr 不可用:", type(e).__name__, str(e)[:80])
    return None


def main():
    sub = sc.read_h5ad(os.path.join(WORK, "04_fibro.h5ad"))
    sub.X = sub.layers["lognorm"].copy()
    print("cells:", sub.shape)

    gs = enrichr_matrisome()
    if gs:
        genes = [g for g in dict.fromkeys(gs) if g in sub.var_names]
        src = f"GO:0062023 ({len(genes)} genes matched)"
    else:
        genes = [g for g in MATRISOME_CORE if g in sub.var_names]
        src = f"curated core matrisome ({len(genes)} genes)"
    print("matrisome source:", src)

    X = sub[:, genes].X
    X = np.asarray(X.todense() if sp.issparse(X) else X, dtype=np.float32)

    # 选择 k：重构误差 + 稳定性
    errs = {}
    for k in range(3, 8):
        m = NMF(n_components=k, init="nndsvda", random_state=SEED, max_iter=500)
        W = m.fit_transform(X)
        errs[k] = float(m.reconstruction_err_)
    print("reconstruction error:", {k: round(v, 2) for k, v in errs.items()})
    k = 5
    nmf = NMF(n_components=k, init="nndsvda", random_state=SEED, max_iter=800)
    W = nmf.fit_transform(X)
    usage = pd.DataFrame(W, index=sub.obs_names, columns=[f"ECM{k+1}" for k in range(k)])
    usage.to_csv(f"{RES}/s9_nmf_usage.csv")

    # 每个程序的 top 基因与类别归属
    rows = []
    for i in range(k):
        h = pd.Series(nmf.components_[i], index=genes).sort_values(ascending=False)
        top = h.head(15)
        cats = {"collagen": sum(g in COLLAGENS for g in top.index),
                "glycoprotein": sum(g in GLYCOPROTEINS for g in top.index),
                "proteoglycan": sum(g in PROTEOGLYCANS for g in top.index),
                "crosslink": sum(g in CROSSLINKERS for g in top.index)}
        name = max(cats, key=cats.get)
        rows.append(dict(program=f"ECM{i+1}", dominant=name, **cats,
                         top_genes=";".join(top.index)))
    prog = pd.DataFrame(rows)
    prog.to_csv(f"{RES}/s9_nmf_programs.csv", index=False)
    print(prog.to_string(index=False))

    sub.obs[[f"ECM{i+1}" for i in range(k)]] = usage.values
    sub.obs["ECM_dominant"] = usage.idxmax(1).values

    # 壁层身份评分
    for name, gs_ in [("Intima", INTIMA), ("Media_id", MEDIA), ("Adventitia_id", ADVENTITIA)]:
        g = [x for x in gs_ if x in sub.var_names]
        if len(g) >= 3:
            sc.tl.score_genes(sub, g, score_name=f"id_{name}", use_raw=False)
    sub.obs["zonation"] = sub.obs["id_Media_id"] - sub.obs["id_Adventitia_id"]

    # 邻域生态位：kNN 邻居的平均 ECM 使用度
    A = sub.obsp["connectivities"].tocsr()
    A = A.multiply(1.0 / np.maximum(np.asarray(A.sum(1)), 1e-9))
    for i in range(k):
        sub.obs[f"nb_ECM{i+1}"] = np.asarray(A @ usage[f"ECM{i+1}"].values).ravel()

    # 组间比较（供体为单位，仅术后样本）
    post = sub.obs[sub.obs.group.isin(["mature_AVF", "failed_AVF"])]
    cols = [f"ECM{i+1}" for i in range(k)] + ["zonation"] + [c for c in sub.obs.columns if c.startswith("nb_ECM")]
    agg = post.groupby(["donor", "group"], observed=True)[cols].mean().reset_index()
    agg.to_csv(f"{RES}/s9_donor_ecm.csv", index=False)

    from scipy import stats
    tests = []
    for c in cols:
        a = agg.loc[agg.group == "failed_AVF", c].values
        b = agg.loc[agg.group == "mature_AVF", c].values
        if len(a) > 1 and len(b) > 1:
            t, p = stats.ttest_ind(a, b, equal_var=False)
            d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2 + 1e-12)
            tests.append(dict(metric=c, failed=a.mean(), mature=b.mean(), diff=a.mean() - b.mean(),
                              cohens_d=d, p=p))
    td = pd.DataFrame(tests).sort_values("p")
    td.to_csv(f"{RES}/s9_ecm_tests.csv", index=False)
    print(td.round(3).to_string(index=False))

    # ECM 主导类型构成
    comp = (post.groupby(["donor", "group", "ECM_dominant"], observed=True).size()
            .unstack(fill_value=0))
    comp.to_csv(f"{RES}/s9_ecm_dominant_counts.csv")

    sub.write_h5ad(os.path.join(WORK, "06_ecm.h5ad"))
    json.dump(dict(source=src, k=k, recon_err=errs), open(f"{RES}/s9_summary.json", "w"), indent=2)
    print("saved 06_ecm.h5ad")


if __name__ == "__main__":
    main()
