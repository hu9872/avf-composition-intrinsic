"""96_fill_repository_doi.py

把稿件与配套文件里的两处占位符替换为真实的仓库 URL 与 Zenodo DOI，并重建 DOCX。

设计要点
--------
1. **只做字符串替换，不改任何科学内容** —— 因此每次替换都打印出改动的文件与次数，
   便于人工复核。
2. **幂等**：重复运行为安全操作；已填过的位置不会被二次改写（替换后不再含占位符）。
3. **自检**：结束时扫描整个投稿包，若仍存在占位符则非零退出，避免带占位符投稿。
4. **顺带重建 DOCX**：正文 DOCX 由 `41_build_iscience_docx.py` 从 Markdown 生成，
   因此必须重建，否则 DOCX 与 MD 不一致（这是最容易被忽略的一步）。

用法
----
    # 1) 先预览
    python code/96_fill_repository_doi.py \
        --repo https://github.com/yourname/avf-composition-intrinsic \
        --doi 10.5281/zenodo.1234567 --dry-run

    # 2) 确认后执行（会顺带重建 DOCX）
    python code/96_fill_repository_doi.py \
        --repo https://github.com/yourname/avf-composition-intrinsic \
        --doi 10.5281/zenodo.1234567

    # 3) 只填 URL、DOI 稍后（Zenodo 归档完成后）
    python code/96_fill_repository_doi.py --repo ... --doi PENDING

注意
----
* Zenodo 会给两个 DOI：**Concept DOI**（随版本更新不变，推荐用于正文"代码长期可用"）
  和 **Version DOI**（对应 v1.0.0 这一次归档）。正文用 Concept DOI，参考文献条目用
  Version DOI 亦可，但全文必须统一——本脚本默认两者都填同一个值，若需区分请分两次运行。
* 若 DOI 尚未拿到，用 `--doi PENDING --no-rebuild`，避免把 "PENDING" 写进 DOCX 后忘记。
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MS = ROOT / "manuscript"

TARGETS = [
    MS / "Manuscript_iScience_v1.md",
    ROOT / "CITATION.cff",
    ROOT / "README.md",
    MS / "submission_checklist_iScience.md",
    ROOT / ".zenodo.json",
]

# 占位符 → 替换值。顺序有意义：长模式先替换。
def rules(repo: str, doi: str | None) -> list[tuple[str, str]]:
    out = [
        ("[repository URL, to be inserted on deposition]", repo),
        ("[repository URL]", repo),
        ("https://github.com/<GITHUB-USER>/<REPO-NAME>", repo),
        ("<REPOSITORY URL>", repo),
    ]
    if doi:                                    # 只有拿到 DOI 才替换 DOI 占位
        doi_url = doi if doi.startswith("http") else f"https://doi.org/{doi}"
        out += [
            ("[DOI, to be inserted on deposition]", doi_url),
            ("[DOI to be assigned]", doi_url),
            ("https://doi.org/[DOI]", doi_url),
            ("[DOI]", doi_url),
            ("https://doi.org/<ZENODO-DOI>", doi_url),
            ("<ZENODO-DOI>", doi_url),
        ]
    return out


REPO_PAT = re.compile(r"repository URL|<GITHUB-USER>|<REPOSITORY URL>|\[repository URL\]")
DOI_PAT = re.compile(r"DOI, to be inserted|DOI to be assigned|<ZENODO-DOI>|https://doi\.org/\[DOI\]")
PLACEHOLDER_PAT = re.compile(REPO_PAT.pattern + "|" + DOI_PAT.pattern)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="GitHub 仓库 URL")
    ap.add_argument("--doi", default=None,
                    help="Zenodo DOI（Concept DOI 优先）。省略则只替换仓库 URL，DOI 占位保留")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-rebuild", action="store_true", help="不重建 DOCX")
    a = ap.parse_args()

    if not a.repo.startswith("http"):
        sys.exit("--repo 必须是完整 URL，例如 https://github.com/you/repo")
    if a.doi is not None and a.doi in ("PENDING", ""):
        sys.exit("DOI 未拿到时不要传 --doi，直接省略该参数即可（只填仓库 URL）")

    total = 0
    for path in TARGETS:
        if not path.exists():
            print(f"  跳过（不存在）：{path.relative_to(ROOT).as_posix()}")
            continue
        text = path.read_text(encoding="utf-8")
        before = text
        n = 0
        for old, new in rules(a.repo, a.doi):
            c = text.count(old)
            if c:
                text = text.replace(old, new)
                n += c
        if text != before:
            print(f"  {path.relative_to(ROOT).as_posix():46s} {n:2d} 处")
            total += n
            if not a.dry_run:
                path.write_text(text, encoding="utf-8")
        else:
            print(f"  {path.relative_to(ROOT).as_posix():46s}  0 处（已填或无占位）")

    print(f"\n{'将替换' if a.dry_run else '已替换'} {total} 处")

    if a.dry_run:
        print("（dry-run，未写盘）")
        return

    # 重建 DOCX，避免 MD 与 DOCX 不一致
    if not a.no_rebuild and a.doi:
        builder = ROOT / "code" / "41_build_iscience_docx.py"
        print(f"\n重建正文 DOCX：{builder.name}")
        r = subprocess.run([sys.executable, str(builder)], cwd=ROOT)
        if r.returncode != 0:
            sys.exit("DOCX 重建失败，请先修好再投稿")

    # 终检：投稿包里不许再有占位符
    leftovers: list[str] = []
    for p in list(MS.rglob("*.md")) + list(MS.rglob("*.docx")) + [ROOT / "README.md", ROOT / "CITATION.cff"]:
        if not p.exists() or p.name.endswith(".new.docx"):   # 兜底另存的临时文件，不算投稿件
            continue
        if p.suffix == ".docx":
            try:
                import docx
                t = "\n".join(x.text for x in docx.Document(p).paragraphs)
            except Exception:
                continue
        else:
            t = p.read_text(encoding="utf-8")
        if PLACEHOLDER_PAT.search(t):
            kind = ("DOI" if (a.doi is None and DOI_PAT.search(t)
                              and not REPO_PAT.search(t)) else "URL/DOI")
            if p.suffix == ".docx" and a.no_rebuild:
                print(f"  [提示] {p.name} 未重建，内含旧占位符——正式投稿前必须重跑 "
                      f"code/41_build_iscience_docx.py")
                continue
            leftovers.append((p.relative_to(ROOT).as_posix(), kind))

    print()
    hard = [f for f, kind in leftovers if kind != "DOI"]
    if hard:
        print("仍有仓库 URL 占位符，请处理：")
        for f in hard:
            print("  ✗", f)
        sys.exit(1)
    if leftovers:
        print("仓库 URL 已填；以下文件仍保留 DOI 占位（拿到 Zenodo DOI 后，"
              "带 --doi 再跑一次即可）：")
        for f, _ in leftovers:
            print("  •", f)
    else:
        print("✓ 投稿包内已无 repository URL / DOI 占位符")
    print(f"  仓库：{a.repo}")
    if a.doi:
        print(f"  DOI ：https://doi.org/{a.doi}")


if __name__ == "__main__":
    main()
