"""把 Markdown 稿件构建为投稿用 DOCX（体例由调用方给定）。

处理：标题/作者/Summary/各级标题/段落/项目符号/markdown 表格；
      [[1,2]] -> 上标引用；**bold**、*italic*、`code`；
      正文 12 pt、双倍行距、连续行号；表格为无底纹可编辑 Word 表格。
"""
import os
import re
import sys
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

import glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _fix_props(path):
    """Replace the docx template's placeholder properties with measured ones.

    python-docx starts from a 2013-vintage template: created 2013-12-23, every
    counter zero, AppVersion 14, and a blank-page thumbnail.  See code/docx_props.py.
    """
    try:
        import docx_props
        info = docx_props.sanitize(path, **docx_props.defaults_for(path))
        print("  props: pages=%(pages)s words=%(words)s lines=%(lines)s "
              "totalmin=%(total_time)s rev=%(revision)s" % info)
    except Exception as exc:                                  # noqa: BLE001
        print("[警告] 文档属性未处理:", exc)


def _latest_manuscript():
    """manuscript/ 下最近修改的 Manuscript_*.md，供未显式指定 MS_SRC 时使用。"""
    cands = [p for p in glob.glob(os.path.join(ROOT, "manuscript", "Manuscript_*.md"))
             if "_backup" not in os.path.basename(p)]
    if not cands:
        raise SystemExit("no manuscript under manuscript/; set MS_SRC")
    cands.sort(key=os.path.getmtime)
    return cands[-1]


def _title_from(path):
    """取稿件首个一级标题作为文档标题。"""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# "):
                return line[2:].strip()
    return ""


SRC = os.environ.get("MS_SRC") or _latest_manuscript()
OUT = os.environ.get("MS_OUT") or os.path.splitext(SRC)[0] + ".docx"
DOC_TITLE = os.environ.get("MS_TITLE") or ""
DOC_SUBJECT = os.environ.get("MS_SUBJECT") or "Research article"


def add_line_numbers(section):
    sectPr = section._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    sectPr.append(ln)


def add_page_numbers(section):
    """Centred page number in the footer; Cell Press asks for numbered pages."""
    par = section.footer.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = par.add_run()
    run.font.size = Pt(10)
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = "PAGE"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    for el in (begin, instr, end):
        run._r.append(el)


STAR = "\ue000"


def clean_inline(text):
    # `\*` marks a literal asterisk — the significance key in a figure legend,
    # e.g. `\*P < 0.05, \*\*P < 0.01, \*\*\*P < 0.001`. Restoring it here would
    # let the emphasis pass below read `**P < 0.01` as bold and swallow the
    # symbol, so park it on a private-use code point and put it back only once
    # the runs have been written.
    return text.replace("\\*", STAR)


TOKEN = re.compile(r"(\[\[[\d,\s]+\]\]|\^[0-9†*]+|\*\*.+?\*\*|`[^`]+`|(?<!\*)\*[^*]+\*(?!\*))")


def add_runs(p, text, size=12, bold=False, italic=False, sup_size=9):
    text = clean_inline(text)
    first_run = len(p.runs)
    for part in TOKEN.split(text):
        if not part:
            continue
        if part.startswith("[[") and part.endswith("]]"):
            nums = [x.strip() for x in part[2:-2].split(",") if x.strip()]
            r = p.add_run(",".join(nums))
            r.font.superscript = True
            r.font.size = Pt(sup_size)
            r.font.name = "Times New Roman"
        elif part.startswith("^") and len(part) > 1:
            r = p.add_run(part[1:])
            r.font.superscript = True
            r.font.size = Pt(size - 3)
            r.font.name = "Times New Roman"
        elif part.startswith("**") and part.endswith("**"):
            r = p.add_run(part[2:-2])
            r.bold = True
            r.font.size = Pt(size)
            r.font.name = "Times New Roman"
        elif part.startswith("`") and part.endswith("`"):
            r = p.add_run(part[1:-1])
            r.font.name = "Consolas"
            r.font.size = Pt(size - 1)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            r = p.add_run(part[1:-1])
            r.italic = True
            r.font.size = Pt(size)
            r.font.name = "Times New Roman"
        else:
            r = p.add_run(part)
            r.bold = bold
            r.italic = italic
            r.font.size = Pt(size)
            r.font.name = "Times New Roman"

    for r in p.runs[first_run:]:
        if STAR in r.text:
            r.text = r.text.replace(STAR, "*")


def flush_table(doc, rows):
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    # Fixed widths: the first column carries the state names, so it needs room to
    # keep words whole; the rest share what is left of the 16 cm text block.
    t.autofit = False
    t.allow_autofit = False
    usable, first = Cm(16.0), Cm(3.3) if ncol > 4 else Cm(4.0)
    rest = Cm((16.0 - first.cm) / (ncol - 1))
    for j, wd in enumerate([first] + [rest] * (ncol - 1)):
        t.columns[j].width = wd
        for c_ in t.columns[j].cells:
            c_.width = wd
    for i, row in enumerate(rows):
        for j, cell in enumerate(row):
            c = t.cell(i, j)
            c.text = ""
            p = c.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_runs(p, cell, size=9)
            if i == 0:
                for r in p.runs:
                    r.bold = True
    doc.add_paragraph()


def build():
    md = open(SRC, encoding="utf-8").read()
    lines = md.split("\n")

    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.5)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    add_line_numbers(sec) if os.environ.get("MS_LINE_NUMBERS", "1") == "1" else None
    add_page_numbers(sec)

    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    st.paragraph_format.space_after = Pt(0)
    st.paragraph_format.space_before = Pt(0)
    st.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

    for name, size in [("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 12)]:
        h = doc.styles[name]
        h.font.name = "Times New Roman"
        h.font.size = Pt(size)
        h.font.bold = True
        h.font.color.rgb = RGBColor(0, 0, 0)
        h.element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        h.paragraph_format.space_before = Pt(12)
        h.paragraph_format.space_after = Pt(6)
        h.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE

    tbl = []
    in_tbl = False
    first_h1 = True

    for raw in lines:
        line = raw.rstrip()
        if line.startswith("|"):
            if re.match(r"^\|[\s:\-|]+\|$", line):
                continue
            tbl.append([c.strip() for c in line.strip().strip("|").split("|")])
            in_tbl = True
            continue
        else:
            if in_tbl:
                flush_table(doc, tbl)
                tbl, in_tbl = [], False

        if not line.strip() or line.strip() == "---":
            continue

        if line.startswith("# "):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(p, line[2:], size=14, bold=True)
            first_h1 = False
        elif line.startswith("## "):
            h = doc.add_heading(level=1)
            add_runs(h, line[3:], size=14, bold=True)
        elif line.startswith("### "):
            h = doc.add_heading(level=2)
            add_runs(h, line[4:], size=12, bold=True)
        elif line.startswith("#### "):
            h = doc.add_heading(level=3)
            add_runs(h, line[5:], size=12, bold=True)
        elif line.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.8)
            add_runs(p, line[2:], size=11)
        elif re.match(r"^[-*]\s", line):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, re.sub(r"^[-*]\s", "", line), size=12)
        elif re.match(r"^\d+\.\s", line) and "References" not in "":
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.75)
            p.paragraph_format.first_line_indent = Cm(-0.75)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
            p.paragraph_format.space_after = Pt(6)
            add_runs(p, line, size=12)
        else:
            p = doc.add_paragraph()
            add_runs(p, line, size=12)

    if in_tbl:
        flush_table(doc, tbl)

    cp = doc.core_properties
    cp.author = cp.last_modified_by = "Daijun Hu"
    cp.title = DOC_TITLE or _title_from(SRC)
    cp.subject = DOC_SUBJECT
    cp.comments = ""

    try:
        doc.save(OUT)
        print("saved", OUT)
        _fix_props(OUT)
    except PermissionError:
        # Word / 预览器占用文件时另存，避免流水线中断
        alt = os.path.splitext(OUT)[0] + ".new.docx"
        doc.save(alt)
        print(f"[警告] {OUT} 被占用（可能正在预览器中打开），已另存为 {alt}")
        print("       关闭该文档后重跑本脚本即可覆盖正式文件名。")
        _fix_props(alt)


if __name__ == "__main__":
    build()
