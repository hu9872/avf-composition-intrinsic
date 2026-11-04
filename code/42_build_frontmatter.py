"""构建投稿前置材料：Title page、Highlights、eTOC blurb、Cover letter。

输出（均在 manuscript/ 下）：
  Title_Page_iScience.docx
  Highlights_iScience.docx
  eTOC_blurb.docx
  Cover_Letter_iScience.docx
"""
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
import re

ROOT = r"D:/AVF_GSE_2611"
OUT = f"{ROOT}/manuscript"

TITLE = ("Apparent transcriptional signatures of arteriovenous fistula maturation "
         "failure are dominated by cell composition")

HIGHLIGHTS = [
    "Pre-specified single-cell reanalysis of human AVF maturation failure",
    "Myofibroblast abundance and all five matrisome programmes are unchanged",
    "cDC2 transcripts fall because myeloid cells expand, not because cDC2 are lost",
    "Composition accounts for much of the apparent transcriptional difference",
]

ETOC = ("Tang et al. reanalyse public single-cell data from human arteriovenous fistulas "
        "and show that transcriptional signatures of maturation failure largely reflect "
        "tissue composition rather than cell-intrinsic reprogramming, so gene-level "
        "biomarkers should be interpreted against the cellular make-up of the sample.")


def manuscript_metrics():
    """Word counts taken from the Markdown source, so the title page cannot drift.

    Returns (summary_words, main_text_words) or None if the source is missing.
    """
    try:
        text = open(f"{ROOT}/manuscript/Manuscript_iScience_v1.md", encoding="utf-8").read()
    except OSError:
        return None
    text = re.sub(r"\[\[\d+\]\]", "", text)

    def block(start, end):
        a = text.find(start)
        b = text.find(end, a + len(start)) if a >= 0 else -1
        return text[a + len(start):b] if b > a else ""

    return (len(block("**Summary**", "**Keywords**").split()),
            len(block("## Introduction", "## Resource availability").split()))


def new_doc(spacing=1.5, size=12):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.5)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    st.paragraph_format.line_spacing = spacing
    st.paragraph_format.space_after = Pt(6)
    cp = doc.core_properties
    cp.author = cp.last_modified_by = "Daijun Hu"
    cp.title = TITLE
    cp.subject = "iScience research article"
    cp.keywords = ("arteriovenous fistula; vascular access; single-cell RNA sequencing; "
                   "cell composition")
    cp.comments = ""
    return doc


def p(doc, text, size=12, bold=False, align=None, space_after=6, italic=False):
    par = doc.add_paragraph()
    if align == "center":
        par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par.paragraph_format.space_after = Pt(space_after)
    r = par.add_run(text)
    r.bold = bold
    r.italic = italic
    r.font.size = Pt(size)
    r.font.name = "Times New Roman"
    return par


def title_page():
    doc = new_doc()
    p(doc, "Title page", 11, italic=True)
    p(doc, TITLE, 14, bold=True)
    p(doc, "")
    p(doc, "Authors", 12, bold=True)
    p(doc, "Yuewu Tang1, Daijun Hu1,*,†")
    p(doc, "1 Department of Nephrology, Chongqing University Three Gorges Hospital, "
           "165 Xincheng Road, Wanzhou District, Chongqing 404000, China")
    p(doc, "† Lead contact")
    p(doc, "* Correspondence: Daijun Hu, 283624940@qq.com")
    p(doc, "")
    p(doc, "Article type", 12, bold=True)
    p(doc, "Research article (iScience)")
    p(doc, "Keywords", 12, bold=True)
    p(doc, "arteriovenous fistula; vascular access; maturation failure; "
           "single-cell RNA sequencing; cell composition; ambient RNA")
    p(doc, "Manuscript metrics", 12, bold=True)
    m = manuscript_metrics()
    if m:
        p(doc, f"Summary, {m[0]} words. Main text (Introduction to Discussion, including "
               f"Limitations of study), approximately {int(round(m[1], -1)):,} words. "
               "Display items, 8 main figures and 3 main tables; STAR Methods is presented "
               "in full in the main document. Supplemental information, Figures S1–S3, "
               "Tables S1–S12, Data S1.")
    else:
        p(doc, "Display items, 8 main figures and 3 main tables; STAR Methods is presented "
               "in full in the main document. Supplemental information, Figures S1–S3, "
               "Tables S1–S12, Data S1.")
    p(doc, "Reporting", 12, bold=True)
    p(doc, "This is a secondary analysis of publicly deposited, de-identified human data. "
           "No new human participants or animals were involved and no ethical approval was "
           "required. The study was not prospectively registered; the analysis plan, including "
           "the three pre-specified hypotheses and the decision rule for the appended cDC2 "
           "hypothesis, was written before the corresponding analyses were run. Generative AI "
           "was not used to generate results; its role in drafting and refactoring analysis "
           "code and in language editing is declared in the manuscript.")
    doc.save(f"{OUT}/Title_Page_iScience.docx")
    print("saved title page")


def highlights():
    doc = new_doc(spacing=1.5)
    p(doc, "Highlights", 14, bold=True)
    p(doc, TITLE, 11, italic=True)
    p(doc, "")
    for h in HIGHLIGHTS:
        assert len(h) <= 85, (len(h), h)
        par = doc.add_paragraph(style="List Bullet")
        r = par.add_run(h)
        r.font.size = Pt(12)
        r.font.name = "Times New Roman"
        par.paragraph_format.space_after = Pt(6)
    doc.save(f"{OUT}/Highlights_iScience.docx")
    print("saved highlights")


def etoc():
    doc = new_doc(spacing=1.5)
    p(doc, "eTOC blurb", 14, bold=True)
    p(doc, TITLE, 11, italic=True)
    p(doc, "")
    p(doc, ETOC, 12)
    doc.save(f"{OUT}/eTOC_blurb.docx")
    print("saved eTOC")


COVER = [
    ("4 November 2026", False),
    ("", False),
    ("The Editors", False),
    ("iScience", False),
    ("Cell Press", False),
    ("", False),
    ("Dear Editors,", False),
    ("", False),
    (f"Please consider our manuscript, \u201c{TITLE}\u201d, for publication in iScience "
     "as a research article.", False),
    ("", False),
    ("Maturation failure affects a large share of autologous arteriovenous fistulas, and "
     "the field's therapeutic hypotheses rest on bulk transcriptomic comparisons that "
     "report hundreds of differentially expressed genes. Our manuscript asks a narrower "
     "question of the same tissue: how much of that difference is how many cells, and how "
     "much is what those cells transcribe. Using public single-cell data from 20 patients "
     "(66,168 cells, six failed and six matured fistulas), we tested three pre-specified "
     "hypotheses with the donor as the unit of replication. None held. The compartment "
     "differences that remained were traced, gene by gene, to ambient RNA from other cell "
     "types and to cell composition: the myeloid compartment was 1.8-fold larger in failed "
     "fistulas while cDC2 numbers were unchanged, so cDC2 transcripts were diluted rather "
     "than lost. A composition-intrinsic decomposition attributed a median 51% of the "
     "remaining signal to composition.", False),
    ("", False),
    ("Three features may suit iScience in particular. First, the result is largely "
     "negative and corrective: a widely used readout, marker abundance in a rare "
     "population, is shown to be the product of two terms that most designs do not "
     "separate. Second, we did not stop at the negative. We quantified what the design can "
     "and cannot detect (80% power only at |d| = 1.80, and no endpoint met an equivalence "
     "criterion), perturbed every threshold that could have produced the result, and "
     "applied the same instruments to a contrast in which large remodelling is known to "
     "occur, where they perform as expected. Third, when an independent bulk cohort failed "
     "to reproduce the one positive compositional finding, we report that discrepancy "
     "rather than omitting it, and we describe the finding as cohort-specific.", False),
    ("", False),
    ("We believe this combination, a falsification of a mechanistic assumption, a reusable "
     "decomposition framework, and an explicit boundary on what the data can support, is "
     "the kind of study iScience has stated it welcomes, including replication studies and "
     "negative results. The work is fully reproducible: code, parameters, seeds and every "
     "result table are archived with a DOI, and all data are public.", False),
    ("", False),
    ("This manuscript is a secondary analysis of publicly deposited data (GSE305423, "
     "GSE220796, GSE119296). It reports no new patients, generates no new materials, and "
     "does not duplicate the primary publication: the primary atlas asked what the cells "
     "in a failing fistula are doing, whereas we ask how much of the observed difference "
     "is compositional. We do not reproduce or reuse any figure or passage from that work. "
     "The manuscript is not under consideration elsewhere, all authors have approved the "
     "submitted version, and we declare no competing interests.", False),
    ("", False),
    ("We suggest the following categories for review: vascular biology, single-cell "
     "genomics, and nephrology. We would be glad to suggest referees on request and to "
     "provide any additional material the editors require.", False),
    ("", False),
    ("Yours sincerely,", False),
    ("", False),
    ("Daijun Hu, on behalf of all authors", False),
    ("Department of Nephrology, Chongqing University Three Gorges Hospital, "
     "165 Xincheng Road, Wanzhou District, Chongqing 404000, China", False),
    ("283624940@qq.com", False),
]


def cover_letter():
    doc = new_doc(spacing=1.15)
    for text, _ in COVER:
        par = doc.add_paragraph()
        par.paragraph_format.space_after = Pt(6)
        r = par.add_run(text)
        r.font.size = Pt(11)
        r.font.name = "Times New Roman"
    doc.save(f"{OUT}/Cover_Letter_iScience.docx")
    print("saved cover letter")


if __name__ == "__main__":
    title_page()
    highlights()
    etoc()
    cover_letter()
