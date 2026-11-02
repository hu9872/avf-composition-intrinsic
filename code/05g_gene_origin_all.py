"""S5g: 对全部区室做基因来源筛查，输出「过滤前/后」差异基因数。

规则（Methods 2.8）：某基因在目标区室的平均表达 >= 其他区室最大值的 50%，
或 >= 0.15 lognorm，方判定为该区室来源。否则视为环境 RNA / 污染，剔除。
"""
import os
import numpy as np
import pandas as pd
import anndata as ad

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"

PARENT = {  # 成纤维亚状态归到成纤维母区室
    "Fib_AF_resident": "Fibroblast", "Fib_AF_stromal": "Fibroblast",
    "Fib_Myofibroblast": "Fibroblast", "Fib_Pericyte": "Fibroblast",
    "Fib_SMC_contractile": "Fibroblast",
}
MAIN = ["Myeloid", "Endothelial", "T_cell", "NK", "B_cell", "Mast", "Fibroblast"]


def main():
    a = ad.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    de = pd.read_csv(f"{RES}/s5c_compartment_de_all.csv")
    sig = de[(de.padj < 0.05) & (de.log2FoldChange.abs() > 0.585)].copy()
    print("总差异基因记录:", len(sig))

    genes = sorted(set(sig.gene) & set(a.var_names))
    print("矩阵中可查到:", len(genes))

    # 各主区室平均表达
    L = {}
    for ct in MAIN:
        m = (a.obs.celltype == ct).values
        if m.sum() == 0:
            continue
        X = np.asarray(a[m, genes].layers["lognorm"].todense())
        L[ct] = X.mean(0)
    M = pd.DataFrame(L, index=genes)
    M.to_csv(f"{RES}/s5g_gene_compartment_means.csv")

    rows = []
    for _, r in sig.iterrows():
        g = r.gene
        if g not in M.index:
            continue
        focal = PARENT.get(r.compartment, r.compartment)
        if focal not in M.columns:
            continue
        v = M.loc[g, focal]
        others = M.loc[g, [c for c in M.columns if c != focal]].max()
        keep = bool(v >= 0.5 * others or v >= 0.15)
        src = M.loc[g].idxmax()
        rows.append(dict(compartment=r.compartment, gene=g, log2FoldChange=r.log2FoldChange,
                         padj=r.padj, focal_mean=round(float(v), 3),
                         max_other_mean=round(float(others), 3),
                         argmax_compartment=src, retained=keep))
    sc = pd.DataFrame(rows)
    sc.to_csv(f"{RES}/s5g_compartment_gene_screen.csv", index=False)

    summ = (sc.groupby("compartment")
            .agg(n_de_before=("retained", "size"), n_de_after=("retained", "sum"))
            .reset_index())
    summ["n_removed"] = summ.n_de_before - summ.n_de_after
    # 补上无差异基因的区室
    allc = sig.compartment.unique()
    summ = summ.set_index("compartment").reindex(
        sorted(set(allc))).fillna(0).astype(int).reset_index()
    summ = summ.sort_values("n_de_before", ascending=False)
    summ.to_csv(f"{RES}/s5g_compartment_counts_before_after.csv", index=False)
    print("\n=== 过滤前 / 后差异基因数 ===")
    print(summ.to_string(index=False))

    print("\n=== 被剔除的基因及其来源 ===")
    rm = sc[~sc.retained].sort_values(["compartment", "padj"])
    print(rm[["compartment", "gene", "log2FoldChange", "focal_mean",
              "max_other_mean", "argmax_compartment"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
