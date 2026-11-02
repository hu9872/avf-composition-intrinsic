"""共享基因集与评分工具。所有基因集来源标注于注释，便于 Methods 交代。"""

# 收缩型 SMC 程序（经典 contractile machinery）
CONTRACTILE = ["MYH11", "CNN1", "ACTA2", "TAGLN", "DES", "MYOCD", "LMOD1", "SYNPO2"]

# 合成型 / 肌成纤维细胞程序（ECM 分泌 + 收缩应力纤维）
SYNTHETIC_MYOF = ["POSTN", "FN1", "COL1A1", "COL3A1", "COL8A1", "TIMP1", "SPARC",
                  "CTGF", "CCN2", "LOX", "TGFB1", "EDN1"]

# 外膜驻留成纤维细胞（adventitial resident）
ADVENTITIAL = ["PI16", "MFAP5", "SFRP2", "CFD", "DCN", "LUM", "GSN", "ABCA8"]

# TGF-β 响应（MSigDB HALLMARK_TGF_BETA_SIGNALING 核心成员节选）
TGFB_RESPONSE = ["SERPINE1", "JUNB", "SKIL", "SMAD7", "SMAD3", "CTGF", "CCN2",
                 "THBS1", "COL1A2", "TGFBR1", "TGFBR2", "LTBP1"]

# 干扰素响应
IFN_RESPONSE = ["IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1", "HLA-DRA", "CXCL9", "CXCL10"]

# 壁层身份
INTIMA = ["PECAM1", "VWF", "CDH5", "EGFL7", "PLVAP"]
MEDIA = ["MYH11", "CNN1", "ACTA2", "TAGLN", "MYOCD"]
ADVENTITIA = ["DCN", "LUM", "PI16", "MFAP5", "SFRP2", "COL1A1", "COL3A1"]

# 核心 matrisome（Naba et al. 2012 核心基质体代表成员，按四类划分）
COLLAGENS = ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
             "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL8A2", "COL12A1", "COL14A1",
             "COL15A1", "COL18A1"]
GLYCOPROTEINS = ["FN1", "FBN1", "FBN2", "LTBP1", "LTBP2", "LTBP4", "THBS1", "THBS2",
                 "POSTN", "TNC", "TNXB", "SPARC", "SPARCL1", "MFAP2", "MFAP5", "ELN",
                 "EMILIN1", "MMRN2", "VWA1"]
PROTEOGLYCANS = ["DCN", "LUM", "BGN", "FMOD", "PRELP", "OGN", "ASPN", "VCAN", "ACAN",
                 "HSPG2", "GPC1", "SDC1", "SDC4", "PODN", "OMD"]
CROSSLINKERS = ["LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4", "TGM2", "PLOD1", "PLOD2", "PLOD3"]

MATRISOME_CORE = [g for s in (COLLAGENS, GLYCOPROTEINS, PROTEOGLYCANS, CROSSLINKERS) for g in s]

PROGRAMS = {
    "Contractile": CONTRACTILE,
    "Synthetic_MyoF": SYNTHETIC_MYOF,
    "Adventitial": ADVENTITIAL,
    "TGFb": TGFB_RESPONSE,
    "IFN": IFN_RESPONSE,
    "Matrisome": MATRISOME_CORE,
    "Crosslink": CROSSLINKERS,
}


def score(adata, gene_sets, prefix="pr_", use_raw=False):
    """对给定基因集打分，缺失基因自动跳过并返回报告。"""
    import scanpy as sc
    report = {}
    for name, genes in gene_sets.items():
        g = [x for x in genes if x in adata.var_names]
        report[name] = (len(g), len(genes))
        if len(g) >= 3:
            sc.tl.score_genes(adata, g, score_name=prefix + name, use_raw=use_raw)
    return report


def present(adata, genes):
    return [g for g in genes if g in adata.var_names]
