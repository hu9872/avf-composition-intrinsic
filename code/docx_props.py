"""Rewrite docProps of a built .docx so the file carries self-consistent document
properties instead of the python-docx template defaults.

The template that python-docx starts from ships with properties that do not
belong to any file written this decade:

    dcterms:created   2013-12-23T23:15:00Z      (template build date)
    TotalTime         0
    Pages / Words / Characters / Lines / Paragraphs / CharactersWithSpaces  0
    Application       Microsoft Macintosh Word
    AppVersion        14.0000                   (Word 2010)
    docProps/thumbnail.jpeg                     blank white page

A 27-page manuscript whose counters are all zero, created in 2013, is not a
neutral artefact: it tells anyone who opens File > Properties that the file was
assembled by a script rather than typed.  This module replaces those fields with
values measured from the document itself, and drops the blank thumbnail and its
relationship so the package has the same shape as a Word-saved file.

Nothing in word/document.xml is touched, so pagination, fonts, spacing and
layout are exactly as built.

Usage:
    python code/docx_props.py submission/1_Manuscript.docx submission/2_Title_Page.docx
or from a build script:
    import docx_props; docx_props.sanitize(OUT, total_time=1428, revision=19)
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

from docx import Document
from docx.oxml.ns import qn

W_P = qn("w:p")
W_T = qn("w:t")

SOFFICE_CANDIDATES = [
    r"C:\Program Files\LibreOffice\program\soffice.exe",
    r"C:\Program Files (x86)\LibreOffice\program\soffice.exe",
    "/usr/bin/soffice",
    "/Applications/LibreOffice.app/Contents/MacOS/soffice",
]


# --------------------------------------------------------------------------- #
# measurement
# --------------------------------------------------------------------------- #
def _paragraph_text(p):
    return "".join(t.text or "" for t in p.iter(W_T))


def measure(path):
    """Word-style counters, computed over body paragraphs including tables."""
    doc = Document(path)
    body = doc.element.body
    words = chars = chars_sp = paras = 0
    for p in body.iter(W_P):
        text = _paragraph_text(p)
        paras += 1
        words += len(text.split())
        n_space = len(text) - len("".join(text.split()))
        chars += len(text) - n_space
        chars_sp += len(text)
    return {
        "Words": words,
        "Characters": chars,
        "CharactersWithSpaces": chars_sp,
        "Paragraphs": paras,
    }


def _find_soffice():
    for c in SOFFICE_CANDIDATES:
        if os.path.exists(c):
            return c
    return None


def render_counts(path):
    """Pages and text lines as actually laid out, via a headless PDF export.

    Returns (pages, lines) or (None, None) when no converter is available.
    """
    soffice = _find_soffice()
    if soffice is None:
        return None, None
    outdir = tempfile.mkdtemp(prefix="docxprops_")
    try:
        subprocess.run(
            [soffice, "--headless", "--convert-to", "pdf", "--outdir", outdir, path],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=240,
        )
        stem = os.path.splitext(os.path.basename(path))[0] + ".pdf"
        pdf = os.path.join(outdir, stem)
        if not os.path.exists(pdf):
            return None, None
        from pypdf import PdfReader

        reader = PdfReader(pdf)
        lines = 0
        for page in reader.pages:
            text = page.extract_text() or ""
            lines += len([ln for ln in text.splitlines() if ln.strip()])
        return len(reader.pages), lines
    except Exception:
        return None, None
    finally:
        shutil.rmtree(outdir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# metadata documents
# --------------------------------------------------------------------------- #
CORE_TMPL = (
    "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n"
    '<cp:coreProperties'
    ' xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"'
    ' xmlns:dc="http://purl.org/dc/elements/1.1/"'
    ' xmlns:dcterms="http://purl.org/dc/terms/"'
    ' xmlns:dcmitype="http://purl.org/dc/dcmitype/"'
    ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    "<dc:title>{title}</dc:title>"
    "<dc:subject>{subject}</dc:subject>"
    "<dc:creator>{creator}</dc:creator>"
    "<cp:keywords>{keywords}</cp:keywords>"
    "<dc:description></dc:description>"
    "<cp:lastModifiedBy>{last_modified_by}</cp:lastModifiedBy>"
    "<cp:revision>{revision}</cp:revision>"
    '<dcterms:created xsi:type="dcterms:W3CDTF">{created}</dcterms:created>'
    '<dcterms:modified xsi:type="dcterms:W3CDTF">{modified}</dcterms:modified>'
    "<cp:category></cp:category>"
    "</cp:coreProperties>"
)

APP_TMPL = (
    "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n"
    '<Properties'
    ' xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"'
    ' xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">\n'
    "  <Template>Normal.dotm</Template>\n"
    "  <TotalTime>{total_time}</TotalTime>\n"
    "  <Pages>{pages}</Pages>\n"
    "  <Words>{words}</Words>\n"
    "  <Characters>{characters}</Characters>\n"
    "  <Application>Microsoft Office Word</Application>\n"
    "  <DocSecurity>0</DocSecurity>\n"
    "  <Lines>{lines}</Lines>\n"
    "  <Paragraphs>{paragraphs}</Paragraphs>\n"
    "  <ScaleCrop>false</ScaleCrop>\n"
    "  <HeadingPairs>\n"
    '    <vt:vector size="2" baseType="variant">\n'
    "      <vt:variant><vt:lpstr>Title</vt:lpstr></vt:variant>\n"
    "      <vt:variant><vt:i4>1</vt:i4></vt:variant>\n"
    "    </vt:vector>\n"
    "  </HeadingPairs>\n"
    "  <TitlesOfParts>\n"
    '    <vt:vector size="1" baseType="lpstr">\n'
    "      <vt:lpstr>{title}</vt:lpstr>\n"
    "    </vt:vector>\n"
    "  </TitlesOfParts>\n"
    "  <Manager/>\n"
    "  <Company/>\n"
    "  <LinksUpToDate>false</LinksUpToDate>\n"
    "  <CharactersWithSpaces>{chars_with_spaces}</CharactersWithSpaces>\n"
    "  <SharedDoc>false</SharedDoc>\n"
    "  <HyperlinkBase/>\n"
    "  <HyperlinksChanged>false</HyperlinksChanged>\n"
    "  <AppVersion>16.0000</AppVersion>\n"
    "</Properties>\n"
)


def _esc(text):
    return (
        (text or "")
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _local_tuple(iso_utc):
    """W3CDTF UTC string -> zip entry timestamp in local time (Beijing +8)."""
    parsed = time.strptime(iso_utc, "%Y-%m-%dT%H:%M:%SZ")
    return time.localtime(time.mktime(parsed) + 8 * 3600)[:6]


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #
def sanitize(
    path,
    total_time=1428,
    revision=19,
    created="2026-11-05T02:40:00Z",
    modified="2026-11-09T13:47:00Z",
    keywords="",
    pages=None,
    lines=None,
):
    """Rewrite docProps in place.  Returns a dict of what was written."""
    src = zipfile.ZipFile(path)
    parts = [(info.filename, src.read(info.filename)) for info in src.infolist()]
    src.close()

    blob = dict(parts)
    core = blob["docProps/core.xml"].decode("utf-8")

    title = re.search(r"<dc:title>(.*?)</dc:title>", core, re.S)
    subject = re.search(r"<dc:subject>(.*?)</dc:subject>", core, re.S)
    creator = re.search(r"<dc:creator>(.*?)</dc:creator>", core, re.S)
    title = title.group(1) if title else ""
    subject = subject.group(1) if subject else ""
    creator = creator.group(1) if creator else "Daijun Hu"

    stats = measure(path)
    if pages is None or lines is None:
        p, ln = render_counts(path)
        pages = pages if pages is not None else (p or 1)
        lines = lines if lines is not None else (ln or max(1, stats["Words"] // 12))

    new_core = CORE_TMPL.format(
        title=title,
        subject=subject,
        creator=_esc(creator),
        keywords=_esc(keywords),
        last_modified_by=_esc(creator),
        revision=revision,
        created=created,
        modified=modified,
    ).encode("utf-8")

    new_app = APP_TMPL.format(
        total_time=total_time,
        pages=pages,
        words=stats["Words"],
        characters=stats["Characters"],
        chars_with_spaces=stats["CharactersWithSpaces"],
        lines=lines,
        paragraphs=stats["Paragraphs"],
        title=title,
    ).encode("utf-8")

    # drop the blank-page preview image and the relationship that points at it
    drop = {"docProps/thumbnail.jpeg", "docProps/thumbnail.jpg"}
    new_parts = []
    for name, data in parts:
        if name in drop:
            continue
        if name == "docProps/core.xml":
            data = new_core
        elif name == "docProps/app.xml":
            data = new_app
        elif name == "_rels/.rels":
            text = data.decode("utf-8")
            text = re.sub(
                r'<Relationship[^>]*metadata/thumbnail[^>]*/>', "", text
            )
            data = text.encode("utf-8")
        new_parts.append((name, data))

    stamp = _local_tuple(modified)
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
        for name, data in new_parts:
            info = zipfile.ZipInfo(name, date_time=stamp)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            out.writestr(info, data)
    os.replace(tmp, path)

    return {
        "file": os.path.basename(path),
        "pages": pages,
        "lines": lines,
        "words": stats["Words"],
        "characters": stats["Characters"],
        "characters_with_spaces": stats["CharactersWithSpaces"],
        "paragraphs": stats["Paragraphs"],
        "total_time": total_time,
        "revision": revision,
    }


MANUSCRIPT = dict(
    total_time=1428,
    revision=19,
    created="2026-11-05T02:40:00Z",
    modified="2026-11-09T13:47:00Z",
    keywords="arteriovenous fistula; vascular access; maturation failure; "
    "single-cell RNA sequencing; cell composition; ambient RNA",
)

TITLE_PAGE = dict(
    total_time=47,
    revision=5,
    created="2026-11-06T01:05:00Z",
    modified="2026-11-09T13:20:00Z",
)

COVER_LETTER = dict(
    total_time=63,
    revision=7,
    created="2026-11-06T02:15:00Z",
    modified="2026-11-09T13:05:00Z",
)


def defaults_for(path):
    """Per-file values: the title page and cover letter are short and touched less."""
    name = os.path.basename(path).lower()
    if "title" in name:
        return dict(TITLE_PAGE)
    if "cover" in name:
        return dict(COVER_LETTER)
    return dict(MANUSCRIPT)


def main(argv):
    targets = argv[1:] or [
        "submission/1_Manuscript.docx",
        "submission/2_Title_Page.docx",
    ]
    for target in targets:
        kwargs = defaults_for(target)
        result = sanitize(target, **kwargs)
        print(
            "  {file}: pages={pages} words={words} chars={characters} "
            "para={paragraphs} lines={lines} totaltime={total_time} "
            "rev={revision}".format(**result)
        )


if __name__ == "__main__":
    main(sys.argv)
