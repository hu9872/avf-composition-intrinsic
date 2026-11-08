"""95_text_overlap_check.py

投稿前的"重复度预检"（text-overlap pre-screen）。

为什么需要这个脚本
------------------
iThenticate 需要机构账号，且结果不可预期。本脚本在提交前做三道可自查的检查：

1. 与"同一批数据的原始论文"比对（ref 8 = Kidney Int Rep 2023；ref 12 = Kidney Int 2026）。
   这两篇是本文数据集 GSE305423 / GSE220796 的原始出处，也是最容易被编辑判定为
   "重复发表 / 自我抄袭"的来源，必须在投稿前量化。
2. 与自己的前序版本比对（本人未发表前稿导出的纯文本），量化文本回收
   （text recycling）比例，用于在 cover letter 中如实声明"本稿为本人前稿的改写"。
3. 输出所有 ≥MIN_RUN 词的连续重合片段及其位置，便于逐句改写或补引注。

方法（与 iThenticate 的近似性）
------------------------------
iThenticate 用的是"连续短语匹配 + 排除引文/参考文献"的启发式。本脚本用
n-gram（默认 8 词）shingle 集合的包含度，再对命中的 shingle 反查最长的连续重合段。
数值不等于 iThenticate 的 similarity index，但可用于定位高风险句子。

用法
----
    python code/95_text_overlap_check.py                 # 默认参数
    python code/95_text_overlap_check.py --n 10          # 更严格的 10 词窗口
    python code/95_text_overlap_check.py --no-self       # 跳过自比

输出
----
    results/s95_overlap_summary.csv      每个来源的总体重合度
    results/s95_overlap_passages.csv     所有 ≥MIN_RUN 词的重合片段
    plan/12_查重预检报告.md               人读报告（由本脚本生成）
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "_ref"
RES = ROOT / "results"
PLAN = ROOT / "plan"


def _manuscript_path():
    """最近修改的稿件；manuscript/ 下没有稿件时返回 None（稿件不随仓库发布）。"""
    cands = [p for p in ROOT.glob("manuscript/Manuscript_*.md")
             if "_backup" not in p.name]
    if not cands:
        return None
    cands.sort(key=lambda p: p.stat().st_mtime)
    return cands[-1]


# 低于此词数的连续重合不报告（iThenticate 亦忽略 5-8 词以下的通用短语）
MIN_RUN = 12

STOP_ONLY = {
    "the", "of", "and", "in", "to", "a", "was", "were", "is", "are", "we", "with",
    "that", "for", "as", "by", "on", "from", "at", "be", "this", "these", "it",
}


# --------------------------------------------------------------------------
# 文本处理
# --------------------------------------------------------------------------
def strip_manuscript(md: str) -> tuple[str, list[tuple[int, str]]]:
    """去掉参考文献、表格、图片占位，保留段落级偏移以便定位。"""
    # 截断参考文献之后的内容
    ref_cut = re.search(r"^## References", md, re.M)
    body = md[: ref_cut.start()] if ref_cut else md

    # 去掉 markdown 表格行与分隔线
    body = "\n".join(
        l for l in body.split("\n")
        if not l.strip().startswith("|") and not re.fullmatch(r"\s*-{3,}\s*", l)
    )
    # 去掉图片占位、脚注式标记
    body = re.sub(r"\[\[FIG:[^\]]+\]\]", " ", body)
    body = re.sub(r"`[^`]*`", " ", body)
    body = re.sub(r"\[\[[\d,\s]+\]\]", " ", body)      # 引用标记不参与比对（iThenticate 亦排除）

    # 段落切分（保留标题作为段落）
    paras: list[tuple[int, str]] = []
    cursor = 0
    for chunk in re.split(r"\n\s*\n", body):
        start = body.find(chunk, cursor)
        cursor = max(cursor, start + len(chunk))
        clean = re.sub(r"^#+\s*", "", chunk.strip())
        clean = re.sub(r"\*\*(.+?)\*\*", r"\1", clean)
        clean = clean.replace("**", "").strip()
        if len(clean.split()) >= 4:
            paras.append((start, clean))
    text = "\n\n".join(p for _, p in paras)
    return text, paras


def tokens(text: str) -> list[str]:
    """保留字母/数字，数字统一为 <num>，避免"66,168 cells"这类数字串伪造匹配。"""
    text = text.lower()
    text = text.replace("\u2212", "-").replace("’", "'")
    raw = re.findall(r"[a-z0-9']+", text)
    out = []
    for t in raw:
        if t.replace(",", "").replace(".", "").isdigit():
            out.append("<num>")
        else:
            out.append(t.strip("'"))
    return [t for t in out if t]


def shingles(toks: list[str], n: int) -> dict[tuple[str, ...], list[int]]:
    idx: dict[tuple[str, ...], list[int]] = {}
    for i in range(len(toks) - n + 1):
        idx.setdefault(tuple(toks[i:i + n]), []).append(i)
    return idx


# --------------------------------------------------------------------------
# 比对
# --------------------------------------------------------------------------
@dataclass
class SourceResult:
    source: str
    label: str
    src_words: int
    matched_shingles: int
    manuscript_shingles: int
    containment_pct: float
    longest_run_words: int
    note: str


def load_sources() -> list[tuple[str, str, str]]:
    """返回 (tag, label, text)。"""
    out: list[tuple[str, str, str]] = []
    r12 = REF / "ref12_fulltext.txt"
    r08 = REF / "ref08_fulltext.txt"
    if r12.exists():
        out.append(("ref12", "Martinez 2026 Kidney Int（本文数据集 GSE305423 原始单细胞论文）",
                    r12.read_text(encoding="utf-8")))
    if r08.exists():
        t = r08.read_text(encoding="utf-8")
        label = ("Martinez 2023 Kidney Int Rep（GSE220796 原始 bulk 论文）"
                 if len(t) > 20000 else
                 "Martinez 2023 Kidney Int Rep（GSE220796 原始 bulk 论文；仅取到 PMC 摘要）")
        out.append(("ref08", label, t))
    # 摘要单独作为高敏感来源：摘要的重合最容易被编辑抓住
    for tag, name in [("abs_10.1016_j.kint.2025.08.024", "ref12 摘要（Elsevier/PubMed 官方摘要）"),
                      ("abs_10.1016_j.ekir.2023.01.008", "ref8 摘要（Elsevier/PubMed 官方摘要）")]:
        p = REF / f"{tag}.txt"
        if p.exists():
            out.append((tag, name, p.read_text(encoding="utf-8")))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="shingle 词数，默认 8")
    ap.add_argument("--min-run", type=int, default=MIN_RUN, help="报告的最小连续重合词数")
    ap.add_argument("--no-self", action="store_true", help="不与自己前稿比对")
    ap.add_argument("--manuscript", default=None,
                    help="稿件路径；缺省取 manuscript/ 下最近修改的 Manuscript_*.md")
    args = ap.parse_args()

    md_path = Path(args.manuscript) if args.manuscript else _manuscript_path()
    if md_path is None or not md_path.exists():
        print("跳过早前自比对：manuscript/ 下没有稿件（稿件不随仓库发布）")
        return
    print(f"稿件：{md_path.relative_to(ROOT) if md_path.is_relative_to(ROOT) else md_path}")
    md = md_path.read_text(encoding="utf-8")
    text, paras = strip_manuscript(md)
    mtoks = tokens(text)
    midx = shingles(mtoks, args.n)
    print(f"正文（不含参考文献/表格）：{len(mtoks)} 词，{len(paras)} 段，"
          f"{len(midx)} 个 {args.n}-gram")

    sources = load_sources()
    if not args.no_self:
        v7 = ROOT / "_v7_dump.txt"
        if v7.exists():
            sources.append(("self_v7", "本人未发表前稿（文本回收自查）",
                            v7.read_text(encoding="utf-8")))

    rows: list[SourceResult] = []
    passages: list[dict] = []

    for tag, label, src in sources:
        stoks = tokens(src)
        sidx = shingles(stoks, args.n)
        hits = set(midx) & set(sidx)
        containment = 100.0 * len(hits) / max(1, len(midx))

        # 逐段定位：哪些段落被该来源覆盖
        para_hits: list[tuple[int, float]] = []
        for pi, (_, p) in enumerate(paras):
            pt = tokens(p)
            if len(pt) < args.n:
                continue
            ps = shingles(pt, args.n)
            inter = len(set(ps) & set(sidx))
            if inter:
                para_hits.append((pi, 100.0 * inter / len(ps)))

        # 连续段：把"命中"标在 shingle 起点上，扫描最长连续 True 段
        run_len = 0
        i = 0
        while i < len(mtoks) - args.n + 1:
            if tuple(mtoks[i:i + args.n]) in sidx:
                j = i
                while j < len(mtoks) - args.n + 1 and tuple(mtoks[j:j + args.n]) in sidx:
                    j += 1
                run_len = j - i + args.n - 1          # j-i 个连续 shingle 覆盖 j-i+n-1 词
                passages.append(dict(
                    source=tag,
                    manuscript_words=f"{i + 1}-{j + args.n - 1}",
                    run_words=run_len,
                    passage=" ".join(mtoks[i:j + args.n - 1]),
                ))
                i = j
            else:
                i += 1
        longest = max([p["run_words"] for p in passages if p["source"] == tag], default=0)

        # 只保留最长的一条，其余按 run 降序截断
        rows.append(SourceResult(
            source=tag, label=label, src_words=len(stoks),
            matched_shingles=len(hits), manuscript_shingles=len(midx),
            containment_pct=round(containment, 2),
            longest_run_words=longest,
            note=f"{len(para_hits)} 段含 ≥{args.n} 词重合"
                 + (f"；最高单段 {max(x[1] for x in para_hits):.0f}%" if para_hits else ""),
        ))
        print(f"  {tag:34s} 重合 {containment:5.2f}%  最长连续 {longest} 词")

    # 每个来源只留前 6 条最长片段，去重
    passages.sort(key=lambda d: (-d["run_words"], d["source"]))
    seen: set[str] = set()
    keep: list[dict] = []
    for p in passages:
        key = p["passage"][:80]
        if key in seen:
            continue
        seen.add(key)
        keep.append(p)
        if len(keep) >= 40:
            break

    # 阳性对照：从 ref 12 原文取 30 词原样拼入，验证检测器确实能命中
    control = dict(target="ref12", injected_words=30, detected_words=0)
    r12 = REF / "ref12_fulltext.txt"
    if r12.exists():
        probe = tokens(r12.read_text(encoding="utf-8"))[40:70]
        sidx = shingles(tokens(r12.read_text(encoding="utf-8")), args.n)
        ctrl_toks = mtoks + probe
        run = 0
        i = 0
        while i < len(ctrl_toks) - args.n + 1:
            if tuple(ctrl_toks[i:i + args.n]) in sidx:
                j = i
                while j < len(ctrl_toks) - args.n + 1 and tuple(ctrl_toks[j:j + args.n]) in sidx:
                    j += 1
                run = max(run, j - i + args.n - 1)
                i = j
            else:
                i += 1
        control["detected_words"] = run
    control["passed"] = control["detected_words"] >= 25
    print(f"\n阳性对照：注入 30 词原文 → 检出 {control['detected_words']} 词连续重合 "
          f"（{'通过' if control['passed'] else '失败，检测器需检查'}）")

    RES.mkdir(exist_ok=True)
    import csv
    with open(RES / "s95_overlap_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    with open(RES / "s95_overlap_passages.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["source", "manuscript_words", "run_words", "passage"])
        w.writeheader()
        w.writerows(keep)

    meta = dict(n=args.n, min_run=args.min_run, manuscript_words=len(mtoks),
                paragraphs=len(paras), positive_control=control,
                sources=[asdict(r) for r in rows])
    (RES / "s95_overlap_meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2),
                                              encoding="utf-8")
    print(f"\n写出 results/s95_overlap_summary.csv、s95_overlap_passages.csv "
          f"（{len(keep)} 条片段）")


if __name__ == "__main__":
    main()
