"""S31: 细胞来源筛查规则与表达下限的敏感性分析。

回答两个可预见的审稿质疑：
  Q1「来源筛查的 50% 规则是否武断？」→ 0.25 / 0.50 / 0.80 三档 + AND/OR 两种组合；
  Q2「0.10 的表达下限是否武断？」→ 0.05 / 0.10 / 0.15 三档，逐基因判定是否可解释。

输入：results/s5g_compartment_gene_screen.csv（30 个区室层面差异基因的全区间均值）
输出：
  results/s31_rule_sensitivity.csv
  results/s31_rule_gene_matrix.csv
  results/s31_floor_sensitivity.csv
  results/s31_summary.json
"""
import json
import numpy as np
import pandas as pd

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"

KEY_GENES = ["CD1C", "APOBEC3A", "TIFAB", "CHN2", "YEATS4", "CCL3L1", "SLC16A10",
             "ACE", "RND1", "MT1X", "NR4A2", "GREM2"]


def main():
    d = pd.read_csv(f"{RES}/s5g_compartment_gene_screen.csv")
    d["max_all"] = np.maximum(d.focal_mean, d.max_other_mean)

    rules = {
        "permissive_r0.25_or_floor": lambda x: (x.focal_mean >= 0.25 * x.max_all) | (x.focal_mean >= 0.15),
        "primary_r0.50_or_floor": lambda x: (x.focal_mean >= 0.50 * x.max_all) | (x.focal_mean >= 0.15),
        "strict_r0.50_and_floor": lambda x: (x.focal_mean >= 0.50 * x.max_all) & (x.focal_mean >= 0.15),
        "very_strict_r0.80_and_floor": lambda x: (x.focal_mean >= 0.80 * x.max_all) & (x.focal_mean >= 0.15),
    }

    mat = d[["compartment", "gene", "log2FoldChange", "padj", "focal_mean",
             "max_other_mean"]].copy()
    rows = []
    for name, fn in rules.items():
        keep = fn(d).values
        mat[name] = keep
        rows.append(dict(rule=name, n_retained=int(keep.sum()), n_dropped=int((~keep).sum())))
    for comp in sorted(d.compartment.unique()):
        m = d.compartment == comp
        r = dict(rule=f"by_compartment:{comp}")
        for name in rules:
            r[name] = int(mat.loc[m, name].sum())
        rows.append(r)
    rule_df = pd.DataFrame(rows)
    rule_df.to_csv(f"{RES}/s31_rule_sensitivity.csv", index=False)
    mat.to_csv(f"{RES}/s31_rule_gene_matrix.csv", index=False)

    print("=== 来源筛查规则敏感性 ===")
    print(rule_df.to_string(index=False))

    # --------------------------------------------------- 表达下限敏感性
    floors = [0.05, 0.10, 0.15]
    frows = []
    for _, r in mat.iterrows():
        if not r["primary_r0.50_or_floor"]:
            continue
        row = dict(compartment=r.compartment, gene=r.gene,
                   log2FoldChange=round(r.log2FoldChange, 3),
                   focal_mean=round(r.focal_mean, 4),
                   max_other_mean=round(r.max_other_mean, 4))
        for f in floors:
            row[f"interpretable_ge_{f}"] = bool(r.focal_mean >= f)
        frows.append(row)
    fdf = pd.DataFrame(frows).sort_values(["compartment", "focal_mean"],
                                          ascending=[True, False])
    fdf.to_csv(f"{RES}/s31_floor_sensitivity.csv", index=False)

    summary = []
    for f in floors:
        col = f"interpretable_ge_{f}"
        sub = fdf[fdf[col]]
        summary.append(dict(
            floor=f, n_interpretable=int(len(sub)),
            compartments=";".join(f"{k}:{v}" for k, v in
                                  sub.groupby("compartment").size().items()),
            myeloid_genes=";".join(sub.loc[sub.compartment == "Myeloid", "gene"]),
            endothelial_genes=";".join(sub.loc[sub.compartment == "Endothelial", "gene"])))
    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{RES}/s31_floor_summary.csv", index=False)
    print("\n=== 表达下限敏感性 ===")
    print(sdf.to_string(index=False))

    # 关键基因的三档判定
    key = fdf[fdf.gene.isin(KEY_GENES)]
    print("\n=== 关键基因 ===")
    print(key.to_string(index=False))

    out = dict(
        n_genes_tested=int(len(d)),
        retained_primary=int(mat["primary_r0.50_or_floor"].sum()),
        retained_permissive=int(mat["permissive_r0.25_or_floor"].sum()),
        retained_strict=int(mat["strict_r0.50_and_floor"].sum()),
        retained_very_strict=int(mat["very_strict_r0.80_and_floor"].sum()),
        floor_counts={str(f): int(sdf.loc[sdf.floor == f, "n_interpretable"].iloc[0])
                      for f in floors},
        stable_across_rules=sorted(
            set(mat.loc[mat["very_strict_r0.80_and_floor"], "gene"]) &
            set(mat.loc[mat["permissive_r0.25_or_floor"], "gene"])),
        unstable_genes=sorted(
            set(mat.loc[mat["permissive_r0.25_or_floor"], "gene"]) -
            set(mat.loc[mat["very_strict_r0.80_and_floor"], "gene"])),
    )
    with open(f"{RES}/s31_summary.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n=== 摘要 ===")
    print(json.dumps({k: v for k, v in out.items() if k != "stable_across_rules"},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
