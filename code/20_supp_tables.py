"""生成补充材料：Table S1-S7 的独立 CSV（supplement/csv/）。

工作簿版本已由 44_make_supp_tables.py 统一输出为 supplement/Data_S1.xlsx，
本脚本只负责 CSV，避免产生第二套同名表格。
"""
import os

import pandas as pd

RES = r"D:/AVF_GSE_2611/results"
SUP = r"D:/AVF_GSE_2611/supplement"
CSVD = f"{SUP}/csv"
os.makedirs(CSVD, exist_ok=True)

S1 = pd.read_csv(f"{RES}/s3_cluster_decision.csv").rename(columns={"Unnamed: 0": "leiden_cluster"})
S1 = S1[["leiden_cluster", "assigned", "n", "marker_support", "best_panel",
         "n_best_markers", "hb_mean", "z_margin", "top_markers"]]
S1.columns = ["Leiden cluster", "Assigned compartment", "Cells", "Rank-sum support (n)",
              "Best panel", "Markers in best panel (n)", "Mean haemoglobin (log-norm)",
              "Margin (z, best - second)", "Top markers"]

S2 = pd.read_csv(f"{RES}/s4_subcluster_annotation.csv").rename(columns={"Unnamed: 0": "leiden_cluster"})
S2 = S2[["leiden_cluster", "label", "n", "margin", "support_best", "support_cross",
         "AF_resident", "AF_stromal", "Myofibroblast", "Pericyte",
         "SMC_contractile", "SMC_synthetic"]]
S2.columns = ["Leiden cluster", "Assigned state", "Cells", "Margin (best - second)",
              "Rank-sum support (n)", "Support across panels", "Resident AF score",
              "Stromal AF score", "Myofibroblast score", "Pericyte score",
              "Contractile SMC score", "Synthetic SMC score"]

S3 = pd.read_csv(f"{RES}/s5_pseudobulk_de.csv")
S3 = S3[["level", "unit", "gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
         "n_failed", "n_mature"]]
S3.columns = ["Level", "Unit", "Gene", "Base mean", "log2FC", "SE", "p", "FDR",
              "Donors failed", "Donors matured"]

S4a = pd.read_csv(f"{RES}/s7_fate_tests.csv")
S4a.columns = ["Fate", "Failed", "Matured", "Difference", "Cohen's d", "Welch p",
               "Mann-Whitney p", "n failed", "n matured", "FDR"]
S4b = pd.read_csv(f"{RES}/s7_donor_fate.csv")
S4b.columns = ["Donor", "Group", "Cells", "Mean pseudotime", "Myofibroblast",
               "Contractile SMC", "Resident AF"]

S5 = pd.read_csv(f"{RES}/s5c_compartment_de_all.csv")
S5 = S5[["compartment", "gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]]
S5.columns = ["Compartment", "Gene", "Base mean", "log2FC", "SE", "p", "FDR"]

S6a = pd.read_csv(f"{RES}/s5e_subtype_marker_check.csv")
S6a.columns = ["State", "Cells"] + [f"{c} (mean log-norm)" for c in S6a.columns[2:]]
S6b = pd.read_csv(f"{RES}/s5e_subtype_counts.csv")
S6b.columns = ["Donor", "Group"] + list(S6b.columns[2:])

S7 = pd.read_csv(f"{RES}/s5e_cdc2_depth_controlled.csv")
S7 = S7[["gene", "metric", "failed_mean", "mature_mean", "diff", "cohens_d",
         "welch_p", "mwu_p", "fdr", "n_failed", "n_mature"]]
S7.columns = ["Gene", "Metric", "Failed", "Matured", "Difference", "Cohen's d",
              "Welch p", "Mann-Whitney p", "FDR", "n failed", "n matured"]

S7b = pd.read_csv(f"{RES}/s5e_cdc2_depth_by_group.csv")
S7b.columns = ["Group", "Median UMIs", "Median genes"]

sig3 = S3[(S3.FDR < 0.05) & (S3.log2FC.abs() > 0.585)]
sig5 = S5[(S5.FDR < 0.05) & (S5.log2FC.abs() > 0.585)]
print("significant rows: S3", len(sig3), "S5", len(sig5))

tables = {
    "S1_cluster_annotation": S1,
    "S2_subcluster_annotation": S2,
    "S3_pseudobulk_DE_sig": sig3,
    "S4_fate_bias": S4a,
    "S5_compartment_DE_sig": sig5,
    "S6_myeloid_markers": S6a,
    "S7_depth_controlled": S7,
}
extra = {"S4b_donor_fate": S4b, "S6b_myeloid_counts": S6b, "S7b_group_depth": S7b}

for name, df in {**tables, **extra}.items():
    df.to_csv(f"{CSVD}/{name}.csv", index=False)
    print(name, df.shape)

print("written:", f"{CSVD}/S*.csv")
