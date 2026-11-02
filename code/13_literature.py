"""S13: 系统性文献对比矩阵。数据源 Europe PMC REST（可核验、可复现）。"""
import json, urllib.parse, urllib.request, socket, time, re
import pandas as pd

socket.setdefaulttimeout(45)
OUT = r"D:/AVF_GSE_2611/results"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

QUERIES = [
    '("arteriovenous fistula" OR "vascular access" OR "venous neointimal hyperplasia") AND '
    '(transcriptomic* OR bioinformatic* OR "single cell" OR "single-cell" OR "RNA-seq" OR WGCNA OR "gene expression profiling")',
    '("arteriovenous fistula") AND ("machine learning" OR "hub gene" OR "differentially expressed")',
]

FIELDS = ["tissue/dataset", "n", "data type", "core method", "main claim", "validation"]


def search(q, pages=3):
    out, cursor = [], "*"
    for _ in range(pages):
        u = f"{BASE}?query={urllib.parse.quote(q)}&format=json&pageSize=100&resultType=core&cursorMark={urllib.parse.quote(cursor)}"
        d = json.loads(urllib.request.urlopen(u).read())
        out += d["resultList"]["result"]
        nc = d.get("nextCursorMark")
        if not nc or nc == cursor:
            break
        cursor = nc
        time.sleep(0.4)
    return out


def main():
    seen, rows = set(), []
    for q in QUERIES:
        for r in search(q):
            key = r.get("doi") or r.get("id")
            if key in seen:
                continue
            seen.add(key)
            yr = r.get("pubYear", "")
            try:
                yr = int(yr)
            except Exception:
                continue
            if yr < 2015:
                continue
            title = re.sub(r"<[^>]+>", "", r.get("title", "") or "")
            ab = re.sub(r"<[^>]+>", "", r.get("abstractText", "") or "")
            # 只保留与内瘘/血管通路直接相关者
            blob = (title + " " + ab).lower()
            if not any(k in blob for k in ["arteriovenous fistula", "vascular access", "neointimal", "avf"]):
                continue
            rows.append(dict(
                year=yr, journal=r.get("journalTitle") or r.get("bookOrReportDetails", {}).get("publisher", ""),
                title=title[:180], pmid=r.get("pmid", ""), doi=r.get("doi", ""),
                has_geo="GSE" in ab,
                geo=";".join(sorted(set(re.findall(r"GSE\d{5,6}", ab))))[:60],
                is_review=(r.get("pubType", "") or "").lower().startswith("review"),
                abstract=ab[:1200],
            ))
    df = pd.DataFrame(rows).sort_values("year", ascending=False)
    df.to_csv(f"{OUT}/s13_literature_raw.csv", index=False)
    print("records:", len(df))
    print(df.groupby("year").size().to_string())
    print("\n--- 含 GEO 编号者 ---")
    print(df[df.has_geo][["year", "journal", "title", "geo"]].to_string(index=False)[:4000])


if __name__ == "__main__":
    main()
