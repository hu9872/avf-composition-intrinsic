"""S5f: 髓系差异信号的定量归因分解。

问题：05c 的髓系 pseudobulk 得到 16 个差异基因（CD1C log2FC -1.9, padj 0.024）。
但 cDC2 群体在失败/成熟间绝对丰度相同（每千细胞 d=0.12），
说明该信号可能来自「髓系构成变化导致的比例稀释」，而非细胞自身转录改变。

本脚本做两件事：
  1. 剔除污染假信号：检查每个「髓系差异基因」是否真的在髓系表达
     （MAL 在 T 细胞 0.22、髓系 0.00，属环境 RNA 污染的典型形态）
  2. 对保留的基因做乘法分解：
        pseudobulk 表达 = 来源亚型占比 × 该亚型细胞内表达率
     验证 log2FC_observed ≈ log2(f_failed/f_mature) + log2(rate_failed/rate_mature)
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

WORK = r"D:/AVF_GSE_2611/data/processed"
RES = r"D:/AVF_GSE_2611/results"
GROUPS = ["mature_AVF", "failed_AVF"]

DE_GENES = ["MAL", "AC135782.1", "DNASE1L3", "APOBEC3A", "ADGRG1", "GP1BB", "TIFAB",
            "CCL3L1", "IL32", "CD1C", "PTGER3", "CBFA2T3", "CHN2", "IL21R",
            "SLC16A10", "YEATS4"]


def main():
    full = ad.read_h5ad(os.path.join(WORK, "03_annotated.h5ad"))
    mye = ad.read_h5ad(os.path.join(WORK, "05d_myeloid.h5ad"))

    # ---------------------------------------------- 1. 基因来源核查（是否真为髓系）
    print("=== 1. 差异基因的细胞来源核查 ===")
    genes = [g for g in DE_GENES if g in full.var_names]
    rows = []
    for g in genes:
        L = np.asarray(full[:, [g]].layers["lognorm"].todense()).ravel()
        C = np.asarray(full[:, [g]].layers["counts"].todense()).ravel()
        rec = {"gene": g, "pct_positive": round(float((C > 0).mean() * 100), 2)}
        for ct in ["Myeloid", "T_cell", "Endothelial", "Fibroblast", "NK"]:
            m = (full.obs.celltype == ct).values
            rec[f"{ct}_mean"] = round(float(L[m].mean()), 3)
        # 判定髓系归属：髓系均值是否不低于其他谱系最高值的 50%
        others = max(rec["T_cell_mean"], rec["Endothelial_mean"],
                     rec["Fibroblast_mean"], rec["NK_mean"])
        rec["myeloid_enriched"] = bool(rec["Myeloid_mean"] >= 0.5 * others or
                                       rec["Myeloid_mean"] >= 0.15)
        rows.append(rec)
    src = pd.DataFrame(rows).set_index("gene")
    src.to_csv(f"{RES}/s5f_gene_source.csv")
    print(src.to_string())
    keep = src[src.myeloid_enriched].index.tolist()
    print(f"\n  确认为髓系来源的基因 ({len(keep)}): {keep}")
    print(f"  判定为非髓系/污染 ({len(genes) - len(keep)}): "
          f"{[g for g in genes if g not in keep]}")

    # ---------------------------------------------- 2. 乘法分解
    print("\n=== 2. 伪批量信号的乘法分解 ===")
    obs = mye.obs
    n_mye = obs.groupby(["donor", "group"], observed=True).size().rename("n_mye")

    out = []

    def decompose(g, subtype):
        """pseudobulk_CPM ≈ f_subtype × rate_subtype；返回两组各成分与预测/实测 log2FC。"""
        m = (obs.mye_subtype == subtype).values
        C = np.asarray(mye[:, [g]].layers["counts"].todense()).ravel()
        tot = obs["total_counts"].values
        d = pd.DataFrame({"donor": obs.donor.values, "group": obs.group.values,
                          "cnt": C, "tot": tot, "in_sub": m.astype(float)})
        agg = d.groupby(["donor", "group"], observed=True).apply(
            lambda x: pd.Series({
                "n_sub": x.in_sub.sum(),
                "sub_cnt": x.loc[x.in_sub == 1, "cnt"].sum(),
                "sub_tot": x.loc[x.in_sub == 1, "tot"].sum(),
                "all_cnt": x.cnt.sum(),
                "all_tot": x.tot.sum(),
            }), include_groups=False).reset_index()
        agg = agg.merge(n_mye.reset_index(), on=["donor", "group"])
        agg["f"] = agg.n_sub / agg.n_mye                      # 亚型占髓系比例
        agg["rate"] = np.where(agg.sub_tot > 0,
                               agg.sub_cnt / agg.sub_tot * 1e4, 0)   # 亚型内 CPM
        agg["pb_cpm"] = agg.all_cnt / agg.all_tot * 1e4       # 实测伪批量 CPM
        a1 = agg[agg.group == "failed_AVF"]
        a2 = agg[agg.group == "mature_AVF"]
        if len(a1) < 4 or len(a2) < 4:
            return None
        eps = 1e-9
        obs_fc = np.log2((a1.pb_cpm.mean() + eps) / (a2.pb_cpm.mean() + eps))
        comp_fc = np.log2((a1.f.mean() + eps) / (a2.f.mean() + eps))
        expr_fc = np.log2((a1.rate.mean() + eps) / (a2.rate.mean() + eps))
        return dict(gene=g, subtype=subtype, log2FC_observed=obs_fc,
                    log2FC_composition=comp_fc, log2FC_intrinsic=expr_fc,
                    log2FC_predicted=comp_fc + expr_fc,
                    failed_f=a1.f.mean(), mature_f=a2.f.mean(),
                    failed_rate=a1.rate.mean(), mature_rate=a2.rate.mean(),
                    failed_pb=a1.pb_cpm.mean(), mature_pb=a2.pb_cpm.mean(), agg=agg)

    # 为每个基因选「主要表达亚型」
    sub_mean = {}
    for s in ["cDC2", "Macrophage", "Mono_classical", "pDC"]:
        msk = (obs.mye_subtype == s).values
        for g in keep:
            L = np.asarray(mye[msk, [g]].layers["lognorm"].todense()).ravel()
            sub_mean[(g, s)] = float(L.mean())
    sm = pd.Series(sub_mean).unstack()
    sm.to_csv(f"{RES}/s5f_subtype_expression.csv")
    print("\n  各基因在亚型内的平均表达:")
    print(sm.round(3).to_string())

    rows = []
    for g in keep:
        s = sm.loc[g].idxmax()
        r = decompose(g, s)
        if r is None:
            continue
        rr = {k: v for k, v in r.items() if k != "agg"}
        rr["dominant_subtype"] = s
        rows.append(rr)
    dec = pd.DataFrame(rows)
    dec["residual"] = dec.log2FC_observed - dec.log2FC_predicted
    dec["pct_from_composition"] = np.where(
        dec.log2FC_observed.abs() > 1e-6,
        100 * dec.log2FC_composition / dec.log2FC_observed, np.nan)
    dec.to_csv(f"{RES}/s5f_decomposition.csv", index=False)
    print("\n  分解结果（log2FC）:")
    print(dec[["gene", "dominant_subtype", "log2FC_observed", "log2FC_composition",
               "log2FC_intrinsic", "log2FC_predicted", "residual",
               "pct_from_composition"]].round(2).to_string(index=False))

    # ---------------------------------------------- 3. 环境 RNA 混杂检验
    print("\n=== 3. 环境 RNA 混杂检验（MAL 与 T 细胞含量） ===")
    tc = full.obs[full.obs.celltype == "T_cell"].groupby(
        ["donor", "group"], observed=True).size().rename("n_t")
    tot_c = full.obs.groupby(["donor", "group"], observed=True).size().rename("n_all")
    tdf = pd.concat([tc, tot_c], axis=1).reset_index()
    tdf["t_frac"] = tdf.n_t / tdf.n_all

    mal = None
    if "MAL" in mye.var_names:
        Cm = np.asarray(mye[:, ["MAL"]].layers["counts"].todense()).ravel()
        mal = pd.DataFrame({"donor": obs.donor.values, "group": obs.group.values,
                            "mal": Cm, "tot": obs.total_counts.values})
        mal = mal.groupby(["donor", "group"], observed=True).apply(
            lambda x: pd.Series({"mal_cpm": x.mal.sum() / x.tot.sum() * 1e4,
                                 "mal_pct_pos": (x.mal > 0).mean() * 100}),
            include_groups=False).reset_index()
        mg = mal.merge(tdf, on=["donor", "group"])
        mg.to_csv(f"{RES}/s5f_mal_vs_tcell.csv", index=False)
        print(mg.round(4).to_string(index=False))
        for col in ["mal_cpm", "mal_pct_pos"]:
            r, p = stats.spearmanr(mg.t_frac, mg[col])
            print(f"  Spearman(T 细胞占比, {col}) = {r:.3f}, p = {p:.4f}")
        fa = mg[mg.group == "failed_AVF"].mal_cpm
        ma = mg[mg.group == "mature_AVF"].mal_cpm
        print(f"  MAL 伪批量 CPM: failed {fa.mean():.2f} vs mature {ma.mean():.2f}, "
              f"log2FC {np.log2((fa.mean()+1e-9)/(ma.mean()+1e-9)):.2f}")

    # ---------------------------------------------- 4. 敏感性：剔除最小供体
    print("\n=== 4. 敏感性分析：剔除髓系细胞量最少的供体 ===")
    comp = pd.read_csv(f"{RES}/s5e_donor_abundance.csv")
    small = comp.nsmallest(2, "n_myeloid")
    print("  剔除:", small[["donor", "group", "n_myeloid"]].to_dict("records"))
    kept = comp[~comp.donor.isin(small.donor)]
    post = kept[kept.group.isin(GROUPS)]
    rows = []
    for metric in ["n_cdc2_D1", "cdc2_D1_frac_myeloid", "cdc2_D1_per1000",
                   "n_myeloid", "myeloid_per1000"]:
        a1 = post.loc[post.group == "failed_AVF", metric].values.astype(float)
        a2 = post.loc[post.group == "mature_AVF", metric].values.astype(float)
        if len(a1) < 3 or len(a2) < 3:
            continue
        t, p = stats.ttest_ind(a1, a2, equal_var=False)
        sp = np.sqrt((a1.var(ddof=1) + a2.var(ddof=1)) / 2)
        d = (a1.mean() - a2.mean()) / sp if sp > 0 else np.nan
        rows.append(dict(metric=metric, failed_mean=a1.mean(), mature_mean=a2.mean(),
                         cohens_d=d, welch_p=p))
    sens = pd.DataFrame(rows)
    sens.to_csv(f"{RES}/s5f_sensitivity.csv", index=False)
    print(sens.round(4).to_string(index=False))


def temp_safe(keep):
    return []


if __name__ == "__main__":
    main()
