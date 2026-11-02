"""参考文献核验：按标题在 Europe PMC 精确检索，输出可用于 Vancouver 格式的字段。"""
import json, urllib.parse, urllib.request, socket, time, sys
import pandas as pd

socket.setdefaulttimeout(45)
RES = r"D:/AVF_GSE_2611/results"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

QUERIES = {
    "avf_scRNA": "The single cell landscape of the human vein after arteriovenous fistula creation and implications for maturation failure",
    "avf_bulk": "The Transcriptomics of the Human Vein Transformation After Arteriovenous Fistula Anastomosis Uncovers Layer-Specific Remodeling",
    "rbms1": "RBMS1 Contributes to Early Arteriovenous Fistula Neointimal Remodeling",
    "cells_time": "Sequential inflammatory and matrisome programs drive remodeling of the mouse carotid-jugular arteriovenous fistula",
    "hubgene2023": "Identification and Validation of Hub Genes in the Stenosis of Arteriovenous Fistula",
    "alpl": "Preoperative Serum Alkaline Phosphatase Predicts Arteriovenous Fistula Maturation Failure",
    "naba": "The matrisome: in silico definition and in vivo characterization by proteomics of normal and neoplastic extracellular matrix",
    "squair": "Confronting false discoveries in single-cell differential expression",
    "paga": "PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells",
    "dpt": "Diffusion pseudotime robustly reconstructs lineage branching",
    "cellphonedb": "CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes",
    "liana": "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data",
    "deseq2": "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2",
    "harmony": "Fast, sensitive and accurate integration of single-cell data with Harmony",
    "scCODA": "scCODA is a Bayesian model for compositional single-cell data analysis",
    "nmf": "Learning the parts of objects by non-negative matrix factorization",
    "outsidein": "Adventitial fibroblasts and the outside-in hypothesis of vascular remodeling",
}


def fetch(q):
    u = (f"{BASE}?query={urllib.parse.quote('TITLE:' + chr(34) + q + chr(34))}"
         f"&format=json&pageSize=3&resultType=core")
    try:
        d = json.loads(urllib.request.urlopen(u).read())
        return d["resultList"]["result"]
    except Exception as e:
        return []


rows = []
for key, q in QUERIES.items():
    res = fetch(q)
    if not res:
        rows.append(dict(key=key, status="NOT_FOUND", query=q[:80]))
        continue
    r = res[0]
    rows.append(dict(
        key=key, status="OK",
        title=r.get("title", "")[:150],
        authors=(r.get("authorString", "") or "")[:120],
        journal=(r.get("journalInfo", {}) or {}).get("journal", {}).get("title", ""),
        year=r.get("pubYear", ""),
        volume=(r.get("journalInfo", {}) or {}).get("volume", ""),
        issue=(r.get("journalInfo", {}) or {}).get("issue", ""),
        pages=r.get("pageInfo", ""),
        doi=r.get("doi", ""), pmid=r.get("pmid", ""),
    ))
    print(f"{key:12s} {rows[-1]['journal'][:35]:35s} {rows[-1]['year']} "
          f"{rows[-1]['volume']}:{rows[-1]['pages']}  PMID {rows[-1]['pmid']}")

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/refs_verified.csv", index=False)
print("\nNOT_FOUND:", df[df.status != "OK"].key.tolist())
