"""按 PMID 抓取最终文献清单，输出 Vancouver 格式与 key->条目 映射。
manuscript_v4 中用 {key} 占位，17_number_refs.py 会按首次出现顺序编号并替换。
"""
import json, urllib.parse, urllib.request, socket, time
import pandas as pd

socket.setdefaulttimeout(45)
RES = r"D:/AVF_GSE_2611/results"
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

KEYS = {
    "patency": "24183112",
    "kdoqi": "32778223",
    "roy06": "16565259",
    "lee2011": "21220751",
    "beathard": "12969170",
    "stenmark": "23216413",
    "maiellaro": "17662969",
    "coen": "21868702",
    "avf_bulk": "37069981",
    "hubgene2023": "36836441",
    "avf_scRNA": "40992573",
    "scanpy": "29409532",
    "harmony": "31740819",
    "scrublet": "30954476",
    "deseq2": "25516281",
    "pydeseq2": "37669147",
    "squair": "34584091",
    "clr": "29187837",
    "soupX": "33367645",
    "vandenbrink": "28960196",
    "paga": "30890159",
    "dpt": "27571553",
    "cellrank": "35027767",
    "cellphonedb": "32103204",
    "naba": "22159717",
    "nmf": "10548103",
    "villani": "28428369",
    "zernecke": "25675999",
    "armulik": "21839917",
    "fleming": "16614314",
    "rnd1": "9531558",
    "integrin": "30195153",
    "altman": "7647644",
}


def fetch_pmid(pmid):
    q = urllib.parse.quote(f"EXT_ID:{pmid} AND SRC:MED")
    u = f"{BASE}?query={q}&format=json&pageSize=1&resultType=core"
    for _ in range(3):
        try:
            d = json.loads(urllib.request.urlopen(u).read())
            r = d["resultList"]["result"]
            return r[0] if r else None
        except Exception:
            time.sleep(2)
    return None


def fmt_authors(s, n=6):
    parts = [p.strip() for p in (s or "").rstrip(".").split(",") if p.strip()]
    if not parts:
        return ""
    if len(parts) > n:
        return ", ".join(parts[:n]) + ", et al"
    return ", ".join(parts)


def clean_title(t):
    t = (t or "").strip()
    if not t.endswith((".", "?", "!")):
        t += "."
    return t


out, rows = {}, []
for key, pmid in KEYS.items():
    r = fetch_pmid(pmid)
    if r is None:
        print("MISS", key, pmid)
        continue
    ji = r.get("journalInfo", {}) or {}
    abbr = (ji.get("journal", {}) or {}).get("medlineAbbreviation", "") or \
           (ji.get("journal", {}) or {}).get("title", "")
    authors = fmt_authors(r.get("authorString", ""))
    title = clean_title(r.get("title", ""))
    year = r.get("pubYear", "")
    vol = ji.get("volume", "")
    iss = ji.get("issue", "")
    pages = r.get("pageInfo", "")
    doi = r.get("doi", "")
    volpart = vol + (f"({iss})" if iss else "")
    loc = f"{year};{volpart}:{pages}" if vol else f"{year}:{pages}"
    entry = f"{authors}. {title} {abbr}. {loc}."
    if doi:
        entry += f" doi:{doi}."
    entry += f" PMID: {pmid}."
    out[key] = entry
    rows.append(dict(key=key, pmid=pmid, journal=abbr, year=year, volume=vol,
                     issue=iss, pages=pages, doi=doi, title=r.get("title", ""),
                     authors=authors, entry=entry))
    print(f"{key:12s} {entry[:120]}")

pd.DataFrame(rows).to_csv(f"{RES}/refs_final.csv", index=False)
with open(r"D:/AVF_GSE_2611/manuscript/refs_map.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\nN =", len(out))
