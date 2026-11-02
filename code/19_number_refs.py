"""按正文首次出现顺序为 {key} 占位符编号，生成参考文献表并统计字数。"""
import json
import re
from collections import OrderedDict

import sys
MS = sys.argv[1] if len(sys.argv) > 1 else r"D:/AVF_GSE_2611/manuscript/manuscript_v4.md"
OUT = sys.argv[2] if len(sys.argv) > 2 else r"D:/AVF_GSE_2611/manuscript/manuscript_v4_final.md"
REFS = json.load(open(r"D:/AVF_GSE_2611/manuscript/refs_map.json", encoding="utf-8"))

text = open(MS, encoding="utf-8").read()

# 修正源数据中的排版瑕疵
for k in REFS:
    REFS[k] = REFS[k].replace("change and  dispersion", "change and dispersion")

order = OrderedDict()
pat = re.compile(r"(?<!\{)\{(?!\{)([A-Za-z0-9_,]+)\}")


def repl(m):
    keys = [k.strip() for k in m.group(1).split(",")]
    nums = []
    for k in keys:
        if k not in REFS:
            raise KeyError(f"unknown reference key: {k}")
        if k not in order:
            order[k] = len(order) + 1
        nums.append(order[k])
    nums = sorted(set(nums))
    return "[" + ",".join(str(n) for n in nums) + "]"


body = pat.sub(repl, text)

unused = [k for k in REFS if k not in order]
if unused:
    print("WARNING unused refs:", unused)

lines = [f"{n}. {REFS[k]}" for k, n in order.items()]
body = body.replace("{{REFERENCES}}", "\n".join(lines))

open(OUT, "w", encoding="utf-8").write(body)
print("references numbered:", len(order), "->", OUT)

# ---------------- 字数统计
def words(s):
    s = re.sub(r"\[[0-9,]+\]", "", s)
    s = re.sub(r"[#*|>]", " ", s)
    return len(re.findall(r"\b[A-Za-z0-9][A-Za-z0-9'’\-]*\b", s))


abstract = body.split("## Abstract")[1].split("**Keywords:**")[0]
kw = body.split("**Keywords:**")[1].split("---")[0]
intro = body.split("## 1. Introduction")[1]
results_start = intro.split("## 4. Discussion")[0]
body_main = body.split("## 1. Introduction")[1].split("## 5. Declarations")[0]
disc = body.split("## 4. Discussion")[1].split("## 5. Declarations")[0]

print("Abstract words:", words(abstract))
print("Keywords:", kw.strip())
for name, seg in [("Introduction", body.split("## 1. Introduction")[1].split("## 2. Methods")[0]),
                  ("Methods", body.split("## 2. Methods")[1].split("## 3. Results")[0]),
                  ("Results", body.split("## 3. Results")[1].split("## 4. Discussion")[0]),
                  ("Discussion", disc)]:
    print(f"{name} words:", words(seg))
print("Body total (Intro-Discussion):", words(body_main))
