"""生成 Document S1 的补充表（T_S1..T_S12）与 Data S1 工作簿。

所有表直接从 results/ 的原始输出裁剪，不做手工录入，保证与正文数字同源。
输出：supplement/csv/T_S*.csv 与 supplement/Data_S1.xlsx
"""
import os
import numpy as np
import pandas as pd

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"
OUT = f"{ROOT}/supplement/csv"
os.makedirs(OUT, exist_ok=True)


def save(df, name):
    df.to_csv(f"{OUT}/{name}.csv", index=False)
    print(name, df.shape)


def main():
    # ---------------------------------------------------------------- S1
    s1 = pd.read_csv(f"{RES}/s3_counts_donor_celltype.csv")
    cols = ["donor", "group", "Fibroblast", "Myeloid", "Endothelial", "T_cell",
            "NK", "Mast", "B_cell", "SMC"]
    s1 = s1[cols].copy()
    s1["total"] = s1[cols[2:]].sum(axis=1)
    tot = {c: "" for c in cols}
    tot["donor"] = "All"
    tot["group"] = ""
    for c in cols[2:]:
        tot[c] = int(s1[c].sum())
    tot["total"] = int(s1["total"].sum())
    s1 = pd.concat([s1, pd.DataFrame([tot])], ignore_index=True)
    save(s1, "T_S1")

    # ---------------------------------------------------------------- S2
    s2 = pd.read_csv(f"{RES}/s4_subcluster_annotation.csv")
    s2 = s2.rename(columns={"Unnamed: 0": "cluster"})
    s2 = s2[["cluster", "label", "n", "margin", "support_best", "support_cross"]]
    save(s2, "T_S2")

    # ---------------------------------------------------------------- S3
    s3 = pd.read_csv(f"{RES}/s5_pseudobulk_de.csv")
    lin = s3[(s3.unit == "lineage") & (s3.padj < 0.05) & (s3.log2FoldChange.abs() > 0.585)]
    lin = lin.sort_values("padj")[["gene", "baseMean", "log2FoldChange", "lfcSE",
                                   "pvalue", "padj"]]
    lin = lin.round(4)
    save(lin, "T_S3")

    # ---------------------------------------------------------------- S4
    s4 = pd.read_csv(f"{RES}/s7_fate_tests.csv").round(4)
    save(s4, "T_S4")

    # ---------------------------------------------------------------- S5
    t = pd.read_csv(f"{RES}/s30_tost.csv")
    s5 = t[["metric", "diff", "cohens_d", "welch_p", "equivalence_bound",
            "ci90_low", "ci90_high", "tost_p", "equivalent_at_d080",
            "excluded_effect_from_d"]].copy()
    s5.columns = ["Endpoint", "Difference (failed - matured)", "Cohen's d", "Welch p",
                  "Equivalence bound (SESOI d=0.8)", "90% CI low", "90% CI high",
                  "TOST p", "Equivalent at d=0.8", "Smallest excluded |d|"]
    save(s5.round(4), "T_S5")

    # ---------------------------------------------------------------- S6
    s6 = pd.read_csv(f"{RES}/s5c_compartment_de_summary.csv")
    s6 = s6[["compartment", "n_failed", "n_mature", "n_cells", "n_tested",
             "n_de", "n_up", "n_down"]]
    save(s6, "T_S6")

    # ---------------------------------------------------------------- S7
    pert = pd.read_csv(f"{RES}/s32_gating_perturbation.csv").round(4)
    pert.columns = ["Threshold scale", "cDC2 cells (total)", "cDC2 clusters",
                    "Failed, per 1,000", "Matured, per 1,000", "Cohen's d", "Welch p"]
    save(pert, "T_S7")

    # ---------------------------------------------------------------- S8
    s8 = pd.read_csv(f"{RES}/s5e_cdc2_downsampled_detection.csv")
    keep = ["donor", "group", "det_CD1C", "det_CLEC10A", "det_FCER1A", "det_CD1E",
            "det_RPL13A", "det_ACTB", "det_B2M"]
    s8 = s8[[c for c in keep if c in s8.columns]].round(4)
    s8.columns = [c.replace("det_", "Detection, ") for c in s8.columns]
    save(s8, "T_S8")

    # ---------------------------------------------------------------- S9
    m = pd.read_csv(f"{RES}/s31_rule_gene_matrix.csv")
    s9 = m[["compartment", "gene", "log2FoldChange", "padj", "focal_mean",
            "max_other_mean", "permissive_r0.25_or_floor", "primary_r0.50_or_floor",
            "strict_r0.50_and_floor", "very_strict_r0.80_and_floor"]].copy()
    s9.columns = ["Compartment", "Gene", "log2FC", "FDR", "Mean in target",
                  "Maximum elsewhere", "Permissive (25%)", "Primary (50%)",
                  "Strict (50% + floor)", "Very strict (80% + floor)"]
    save(s9.round(4), "T_S9")

    # ---------------------------------------------------------------- S10
    s10 = pd.read_csv(f"{RES}/s31_floor_sensitivity.csv")
    s10 = s10[["compartment", "gene", "log2FoldChange", "focal_mean",
               "interpretable_ge_0.05", "interpretable_ge_0.1", "interpretable_ge_0.15"]]
    s10.columns = ["Compartment", "Gene", "log2FC", "Mean in target",
                   "Interpretable at 0.05", "at 0.10", "at 0.15"]
    save(s10.round(4), "T_S10")

    # ---------------------------------------------------------------- S11
    s11 = pd.read_csv(f"{RES}/s34_compartment_score_tests.csv").round(4)
    s11.columns = ["Compartment", "Signature genes", "Failed, mean z", "Matured, mean z",
                   "Difference", "Cohen's d", "Welch p", "n failed", "n matured"]
    save(s11, "T_S11")

    # ---------------------------------------------------------------- S12
    s12 = pd.read_csv(f"{RES}/s32_donor_demographics.csv")
    s12 = s12[["donor", "group", "sex", "age", "n_myeloid", "n_all_cells"]].copy()
    s12.columns = ["Donor", "Sampling group", "Sex", "Age (years)",
                   "Myeloid cells", "All cells"]
    save(s12, "T_S12")

    # ------------------------------------------------------------ Data S1
    with pd.ExcelWriter(f"{ROOT}/supplement/Data_S1.xlsx", engine="openpyxl") as xl:
        s1.to_excel(xl, sheet_name="S1_compartment_counts", index=False)
        pd.read_csv(f"{RES}/s5e_subtype_counts.csv").to_excel(
            xl, sheet_name="S2_myeloid_counts", index=False)
        pd.read_csv(f"{RES}/s5_subtype_fraction.csv").round(6).to_excel(
            xl, sheet_name="S3_mural_fractions", index=False)
        pd.read_csv(f"{RES}/s5_composition_tests.csv").round(6).to_excel(
            xl, sheet_name="S4_composition_tests", index=False)
        lin.to_excel(xl, sheet_name="S5_fibroblast_DE", index=False)
        pd.read_csv(f"{RES}/s5f_decomposition_retained.csv").round(6).to_excel(
            xl, sheet_name="S6_decomposition", index=False)
        pd.read_csv(f"{RES}/s5g_compartment_gene_screen.csv").round(6).to_excel(
            xl, sheet_name="S7_gene_origin", index=False)
        pd.read_csv(f"{RES}/s8c_pathway_tests.csv").round(6).to_excel(
            xl, sheet_name="S8_ligand_receptor_axes", index=False)
        pd.read_csv(f"{RES}/s9_donor_ecm.csv").round(6).to_excel(
            xl, sheet_name="S9_ecm_donor", index=False)
        pd.read_csv(f"{RES}/s10_deseq_AVF_failed_vs_matured.csv").round(6).to_excel(
            xl, sheet_name="S10_bulk_DE", index=False)
        t.to_excel(xl, sheet_name="S11_power_TOST", index=False)
        s9.to_excel(xl, sheet_name="S12_rule_sensitivity", index=False)
        s11.to_excel(xl, sheet_name="S13_signature_scores", index=False)
    print("Data_S1.xlsx written")


if __name__ == "__main__":
    main()
