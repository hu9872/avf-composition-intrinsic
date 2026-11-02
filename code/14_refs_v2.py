"""扩充版参考文献核验：每条按标题在 Europe PMC 精确检索，输出真实书目字段。
用途：为 manuscript_v4 提供 28-32 条可核验文献；无法检索到的一律标记并剔除。
"""
import json, urllib.parse, urllib.request, socket, time
import pandas as pd

socket.setdefaulttimeout(45)
RES = r"D:/AVF_GSE_2611/results"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# key: (检索标题, 在正文中的用途)
QUERIES = {
    # ---- 背景：疾病负担与病理 ----
    "patency": ("Patency rates of the arteriovenous fistula for hemodialysis: a systematic review and meta-analysis",
                "AVF 通畅率与失败负担"),
    "kdoqi": ("KDOQI Clinical Practice Guideline for Vascular Access: 2019 Update",
              "成熟失败的定义/临床标准"),
    "roy06": ("Hemodialysis vascular access dysfunction: a cellular and molecular viewpoint",
              "新生内膜增生与透析通路功能障碍的细胞机制"),
    "lee2011": ("Severe venous neointimal hyperplasia prior to dialysis access surgery",
                "术前已存在静脉内膜增生"),
    "beathard": ("Aggressive treatment of early fistula failure",
                 "早期失败临床处理与表型异质性"),
    "stenmark": ("The adventitia: essential regulator of vascular wall structure and function",
                 "外膜成纤维细胞/outside-in 假说"),
    "maiellaro": ("The role of the adventitia in vascular inflammation",
                  "外膜炎症与细胞迁移"),
    "coen": ("Myofibroblast-mediated adventitial remodeling: an underestimated player in arterial pathology",
             "外膜肌成纤维细胞来源"),
    # ---- 数据来源与既有转录组 ----
    "avf_bulk": ("The Transcriptomics of the Human Vein Transformation After Arteriovenous Fistula Anastomosis Uncovers Layer-Specific Remodeling",
                 "配对 bulk 转录组、层特异性重塑"),
    "avf_scRNA": ("The single cell landscape of the human vein after arteriovenous fistula creation and implications for maturation failure",
                  "本研究所用单细胞队列"),
    "hubgene2023": ("Identification and Validation of Hub Genes in the Stenosis of Arteriovenous Fistula",
                    "bulk 派生候选基因/机器学习分类器"),
    "alpl": ("Preoperative Serum Alkaline Phosphatase Predicts Arteriovenous Fistula Maturation Failure",
             "基于 bulk/临床标志物的预测模型"),
    # ---- 方法 ----
    "scanpy": ("SCANPY: large-scale single-cell gene expression data analysis", "Scanpy"),
    "harmony": ("Fast, sensitive and accurate integration of single-cell data with Harmony", "Harmony"),
    "scrublet": ("Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data", "双细胞去除"),
    "deseq2": ("Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2", "DESeq2 模型"),
    "pydeseq2": ("PyDESeq2: a python package for bulk RNA-seq differential expression analysis", "PyDESeq2 实现"),
    "paga": ("PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells", "PAGA"),
    "dpt": ("Diffusion pseudotime robustly reconstructs lineage branching", "扩散拟时序"),
    "cellrank": ("CellRank for directed single-cell fate mapping", "多命运吸收概率框架"),
    "cellphonedb": ("CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes", "配体-受体打分逻辑"),
    "naba": ("The matrisome: in silico definition and in vivo characterization by proteomics of normal and neoplastic extracellular matrix", "matrisome 基因集"),
    "nmf": ("Learning the parts of objects by non-negative matrix factorization", "NMF"),
    "squair": ("Confronting false discoveries in single-cell differential expression", "以供体为统计单位/pseudobulk"),
    "murphy": ("A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis", "pseudobulk 优于细胞水平检验"),
    "sccoda": ("scCODA is a Bayesian model for compositional single-cell data analysis", "细胞组成分析"),
    "clr": ("Microbiome Datasets Are Compositional: And This Is Not Optional", "CLR 变换"),
    "soupX": ("SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data", "环境 RNA 污染"),
    "decontx": ("DecontX: A Bayesian Method for Decontamination of Single-Cell RNA-Seq Data", "环境 RNA 校正方法"),
    # ---- 生物学解释 ----
    "villani": ("Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors", "cDC2 定义与标记"),
    "zernecke": ("Dendritic cells in atherosclerosis: evidence in mice and humans", "血管壁树突状细胞"),
    "armulik": ("Pericytes: developmental, physiological, and pathological perspectives, problems, and promises", "周细胞与微血管稳定"),
    "fleming": ("Signaling by the angiotensin-converting enzyme", "内皮 ACE 与血管重塑"),
    "rnd1": ("A new member of the Rho family, Rnd1, promotes disassembly of actin filament structures and loss of cell adhesion", "RND1 与细胞黏附/骨架"),
    "integrin": ("Signal transduction via integrin adhesion complexes", "整合素黏附信号"),
    "fn_myof": ("Fibronectin-integrin interactions in myofibroblast differentiation", "纤连蛋白-整合素与肌成纤维细胞分化"),
}


def fetch(q):
    u = (f"{BASE}?query={urllib.parse.quote('TITLE:' + chr(34) + q + chr(34))}"
         f"&format=json&pageSize=3&resultType=core")
    for attempt in range(3):
        try:
            d = json.loads(urllib.request.urlopen(u).read())
            return d["resultList"]["result"]
        except Exception:
            time.sleep(2)
    return []


rows = []
for key, (q, purpose) in QUERIES.items():
    res = fetch(q)
    if not res:
        rows.append(dict(key=key, status="NOT_FOUND", query=q[:90], purpose=purpose))
        print(f"{key:12s} NOT_FOUND  {q[:60]}")
        continue
    r = res[0]
    ji = r.get("journalInfo", {}) or {}
    rows.append(dict(
        key=key, status="OK", purpose=purpose,
        title=r.get("title", "").strip(),
        authors=(r.get("authorString", "") or "").strip(),
        journal=(ji.get("journal", {}) or {}).get("title", ""),
        journal_abbr=(ji.get("journal", {}) or {}).get("medlineAbbreviation", ""),
        year=r.get("pubYear", ""),
        volume=ji.get("volume", ""), issue=ji.get("issue", ""),
        pages=r.get("pageInfo", ""),
        doi=r.get("doi", ""), pmid=r.get("pmid", ""),
        pubtype=";".join((r.get("pubTypeList", {}) or {}).get("pubType", []) or []),
    ))
    x = rows[-1]
    print(f"{key:12s} OK  {x['journal_abbr'][:28]:28s} {x['year']} {x['volume']}:{x['pages']} PMID {x['pmid']}")

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/refs_verified_v2.csv", index=False)
print("\nTOTAL:", len(df), "NOT_FOUND:", df[df.status != "OK"].key.tolist())
