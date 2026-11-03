"""构建 Document S1（补充材料合并 PDF 的 Word 源）：方法 + 图注 + 表 S1–S12。

[[TABLE:T_Sx]] 指令处插入 supplement/csv/T_Sx.csv，数值统一四舍五入到 4 位有效数字，
避免长小数把表格撑破。
输出：supplement/Document_S1.docx
"""
import os
import re
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = r"D:/AVF_GSE_2611"
SRC = f"{ROOT}/supplement/Document_S1.md"
CSV = f"{ROOT}/supplement/csv"
OUT = f"{ROOT}/supplement/Document_S1.docx"

TOKEN = re.compile(r"(\*\*.+?\*\*|`[^`]+`|(?<!\*)\*[^*]+\*(?!\*))")


def add_runs(p, text, size=10.5, bold=False):
    for part in TOKEN.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            r = p.add_run(part[2:-2]); r.bold = True
        elif part.startswith("`") and part.endswith("`"):
            r = p.add_run(part[1:-1]); r.font.name = "Consolas"; r.font.size = Pt(size - 1)
            continue
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.italic = True
        else:
            r = p.add_run(part)
            r.bold = bold
        r.font.size = Pt(size)
        r.font.name = "Times New Roman"


def fmt(v):
    if isinstance(v, (np.floating, float)):
        if not np.isfinite(v):
            return ""
        a = abs(v)
        if a == 0:
            return "0"
        if a >= 1000:
            return f"{v:,.1f}"
        if a >= 1:
            return f"{v:.3f}".rstrip("0").rstrip(".")
        return f"{v:.4f}".rstrip("0").rstrip(".")
    return str(v)


def insert_csv_table(doc, name):
    df = pd.read_csv(f"{CSV}/{name}.csv")
    small = df.shape[1] <= 6
    size = 7.2 if df.shape[1] <= 5 else (6.4 if df.shape[1] <= 8 else 5.8)
    t = doc.add_table(rows=len(df) + 1, cols=df.shape[1])
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for j, c in enumerate(df.columns):
        cell = t.cell(0, j); cell.text = ""
        p = cell.paragraphs[0]; p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        r = p.add_run(str(c)); r.bold = True; r.font.size = Pt(size)
        r.font.name = "Times New Roman"
    for i in range(len(df)):
        for j, c in enumerate(df.columns):
            cell = t.cell(i + 1, j); cell.text = ""
            p = cell.paragraphs[0]; p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            r = p.add_run(fmt(df.iloc[i, j])); r.font.size = Pt(size)
            r.font.name = "Times New Roman"
    doc.add_paragraph()


def add_page_numbers(section):
    """Centred page number in the footer, so reviewers can cite supplement pages."""
    par = section.footer.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = par.add_run()
    run.font.size = Pt(9)
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = "PAGE"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    for el in (begin, instr, end):
        run._r.append(el)


def build():
    md = open(SRC, encoding="utf-8").read()
    lines = md.split("\n")
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.0)
    sec.top_margin = sec.bottom_margin = Cm(2.0)
    add_page_numbers(sec)

    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    st.paragraph_format.line_spacing = 1.15
    st.paragraph_format.space_after = Pt(4)

    for name, size in [("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 11)]:
        h = doc.styles[name]
        h.font.name = "Times New Roman"; h.font.size = Pt(size); h.font.bold = True
        h.font.color.rgb = None
        h.paragraph_format.space_before = Pt(10)
        h.paragraph_format.space_after = Pt(4)

    tbl = []
    in_tbl = False
    for raw in lines:
        line = raw.rstrip()
        m = re.match(r"^\[\[TABLE:(T_S\d+)\]\]$", line.strip())
        if m:
            insert_csv_table(doc, m.group(1))
            continue
        fm = re.match(r"^\[\[FIG:(.+?)\]\]$", line.strip())
        if fm:
            path = f"{ROOT}/figures/{fm.group(1)}"
            if os.path.exists(path):
                doc.add_picture(path, width=Cm(15.5))
                doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            continue
        if line.startswith("|"):
            if re.match(r"^\|[\s:\-|]+\|$", line):
                continue
            tbl.append([c.strip() for c in line.strip().strip("|").split("|")])
            in_tbl = True
            continue
        else:
            if in_tbl:
                insert_md_table(doc, tbl)
                tbl, in_tbl = [], False
        if not line.strip() or line.strip() == "---":
            continue
        if line.startswith("### "):
            h = doc.add_heading(level=2); add_runs(h, line[4:], size=12, bold=True)
        elif line.startswith("## "):
            h = doc.add_heading(level=1); add_runs(h, line[3:], size=14, bold=True)
        elif line.startswith("# "):
            h = doc.add_heading(level=0); add_runs(h, line[2:], size=16, bold=True)
        elif re.match(r"^[-*]\s", line):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, re.sub(r"^[-*]\s", "", line), size=10.5)
        else:
            p = doc.add_paragraph()
            add_runs(p, line, size=10.5)

    if in_tbl:
        insert_md_table(doc, tbl)
    doc.save(OUT)
    print("saved", OUT)


def insert_md_table(doc, rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        for j, c in enumerate(row):
            cell = t.cell(i, j); cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_runs(p, c, size=8)
            if i == 0:
                for r in p.runs:
                    r.bold = True
    doc.add_paragraph()


if __name__ == "__main__":
    build()
