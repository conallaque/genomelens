"""
PDF Export
==========

Converts the generated HTML report into a paginated, print-quality PDF.
Adds:
  * Cover page with title, date, file hash, disclaimer
  * Table of contents (auto-generated from <h2> anchors in the report)
  * Page numbers in the footer
  * Print-optimised CSS that overrides the on-screen dark theme

Requires weasyprint:   pip install weasyprint
weasyprint also needs system libs (pango, cairo). On macOS:
   brew install pango libffi
"""

from __future__ import annotations

import datetime
import re
from pathlib import Path


def weasyprint_available() -> bool:
    try:
        # The import IS the probe — it is not an unused import. An
        # unused-import autofix deleted this line, which made the function
        # return True unconditionally and turned a graceful "PDF export
        # unavailable" into a crash.
        import weasyprint  # noqa: F401
        return True
    except Exception:
        return False


# ── PDF-specific CSS overrides ────────────────────────────────────────────────
PDF_CSS = """
@page {
    size: A4;
    margin: 18mm 14mm 22mm 14mm;
    @bottom-center {
        content: counter(page) " / " counter(pages);
        font-size: 9pt;
        color: #888;
    }
    @top-right {
        content: "DNA Analysis Report — confidential";
        font-size: 8pt;
        color: #aaa;
    }
}
body {
    background: #ffffff !important;
    color: #1f2328 !important;
    font-size: 9.5pt !important;
    line-height: 1.45 !important;
}
.hdr {
    position: static !important;
    border-bottom: 2px solid #444 !important;
    background: #fff !important;
    page-break-after: avoid;
}
.toc { display: none; }
section {
    page-break-inside: avoid;
}
h2 {
    page-break-after: avoid;
    color: #111 !important;
    font-size: 15pt !important;
}
.apoe-section, .exec-section, .ydna-section, .mtdna-section,
.qc-section, .prs-section, .pgx-section, .inter-section, .carr-section,
.traits-section, .couns-section, .cross-cat-section, .ref-section,
.rec-section, .meth-section, .cat-section {
    background: #fff !important;
    border: 1px solid #d0d7de !important;
    color: #1f2328 !important;
}
.stat-n { color: #0969da !important; }
.snp-tbl, .qc-tbl, .pgx-tbl, .ref-tbl { font-size: 8pt !important; }
.snp-tbl th, .qc-tbl th, .pgx-tbl th, .ref-tbl th {
    background: #f0f0f0 !important; color: #444 !important;
}
.snp-tbl tr:hover, .qc-tbl tr:hover { background: transparent !important; }
.disc {
    border-left: 3px solid #1f2328 !important;
    background: #f6f8fa !important;
}
.ai-section { background: #f6f8fa !important; border-color: #d0d7de !important; }
.ftr { color: #888 !important; }
"""


def _build_cover_page(file_label: str, file_hash: str, version: str,
                      report_date: str, qc_grade: str = "") -> str:
    return f"""
<div class="pdf-cover">
  <div class="pdf-cover-brand">DNA Analysis Report</div>
  <div class="pdf-cover-tagline">Professional-grade local genomic analysis</div>
  <div class="pdf-cover-meta">
    <div><span class="pdf-cover-label">Generated:</span> {report_date}</div>
    <div><span class="pdf-cover-label">Source file:</span> {file_label}</div>
    <div><span class="pdf-cover-label">File hash (SHA-256, first 16):</span> <code>{file_hash}</code></div>
    <div><span class="pdf-cover-label">Tool version:</span> {version}</div>
    {f'<div><span class="pdf-cover-label">QC grade:</span> {qc_grade}</div>' if qc_grade else ''}
  </div>
  <div class="pdf-cover-disc">
    <strong>Research / Educational Use — Not a Clinical Diagnostic.</strong><br>
    This report integrates a curated variant database, CPIC-style pharmacogenomic
    phenotyping, curated polygenic risk scores, compound-heterozygosity detection,
    carrier-status analysis, trait predictions, and optional local-AI interpretation.
    Findings — especially pharmacogenomic phenotypes and pathogenic variants —
    should be confirmed by a CLIA/CAP-certified clinical laboratory before any
    medical decision is made. Always consult a qualified physician, clinical
    pharmacist, or board-certified genetic counsellor.
  </div>
</div>
<style>
.pdf-cover {{
    page-break-after: always;
    padding: 24mm 20mm;
    min-height: 240mm;
    display: flex;
    flex-direction: column;
    justify-content: center;
    background: linear-gradient(135deg, #f6f8fa 0%, #eaeef2 100%);
}}
.pdf-cover-brand {{
    font-size: 32pt; font-weight: 800; letter-spacing: -.6pt;
    color: #1f2328; margin-bottom: 4pt;
}}
.pdf-cover-tagline {{
    font-size: 11pt; color: #57606a; margin-bottom: 28pt;
    font-style: italic;
}}
.pdf-cover-meta {{
    background: #fff; border: 1px solid #d0d7de; border-radius: 6pt;
    padding: 14pt 18pt; font-size: 10pt; line-height: 1.9;
    margin-bottom: 28pt;
}}
.pdf-cover-label {{
    font-weight: 600; color: #57606a; display: inline-block;
    min-width: 60mm;
}}
.pdf-cover-disc {{
    background: #fff7e6; border-left: 4pt solid #d29922;
    padding: 12pt 16pt; font-size: 9pt; line-height: 1.7;
    color: #1f2328; border-radius: 4pt;
}}
</style>
"""


def _build_toc(html: str) -> str:
    """Extract <h2> headings + anchor IDs and build a TOC."""
    pattern = re.compile(
        r'<section[^>]*id="([^"]+)"[^>]*>\s*<h2[^>]*>(.*?)</h2>',
        re.IGNORECASE | re.DOTALL,
    )
    items = []
    for m in pattern.finditer(html):
        anchor = m.group(1)
        title_raw = m.group(2)
        # strip tags
        title = re.sub(r"<[^>]+>", "", title_raw).strip()
        if title:
            items.append((anchor, title))
    if not items:
        return ""
    lis = "\n".join(
        f'<li><a href="#{a}">{t}</a> <span class="pdf-toc-leader"></span> '
        f'<span class="pdf-toc-page" data-href="#{a}"></span></li>'
        for a, t in items
    )
    return f"""
<div class="pdf-toc">
  <h2>Table of Contents</h2>
  <ul>{lis}</ul>
</div>
<style>
.pdf-toc {{
    page-break-after: always;
    padding: 18mm 14mm;
}}
.pdf-toc h2 {{
    font-size: 18pt; border-bottom: 2pt solid #444;
    padding-bottom: 6pt; margin-bottom: 18pt;
}}
.pdf-toc ul {{ list-style: none; padding-left: 0; }}
.pdf-toc li {{
    display: flex; align-items: baseline; gap: 6pt;
    padding: 4pt 0; font-size: 10pt;
}}
.pdf-toc a {{ color: #1f2328 !important; text-decoration: none; }}
.pdf-toc-leader {{
    flex: 1; border-bottom: 1pt dotted #999; margin: 0 6pt;
    transform: translateY(-2pt);
}}
.pdf-toc-page::before {{
    content: target-counter(attr(data-href url), page);
    color: #57606a; font-size: 9pt;
}}
</style>
"""


def html_to_pdf(
    html_path: Path,
    pdf_path: Path,
    file_label: str = "",
    file_hash: str = "",
    version: str = "",
    qc_grade: str = "",
) -> str:
    """Render the HTML report at html_path to PDF at pdf_path with cover + TOC.

    Returns a status string suitable for logging.
    """
    if not weasyprint_available():
        return ("weasyprint not installed. Install via: "
                "pip install weasyprint  (and on macOS: brew install pango libffi)")
    from weasyprint import CSS, HTML

    html = Path(html_path).read_text()
    report_date = datetime.datetime.now().strftime("%B %d, %Y at %H:%M")
    cover = _build_cover_page(
        file_label=file_label or "raw DNA file",
        file_hash=file_hash or "n/a",
        version=version or "v3.0.0",
        report_date=report_date,
        qc_grade=qc_grade,
    )
    toc = _build_toc(html)

    # Insert cover + TOC right after <body>
    html_modified = re.sub(
        r"(<body[^>]*>)",
        lambda m: m.group(1) + cover + toc,
        html, count=1, flags=re.IGNORECASE,
    )

    try:
        HTML(string=html_modified, base_url=str(html_path.parent)).write_pdf(
            str(pdf_path),
            stylesheets=[CSS(string=PDF_CSS)],
        )
        size = pdf_path.stat().st_size
        return f"PDF written: {pdf_path} ({size/1e6:.1f} MB)"
    except Exception as e:
        return f"PDF render failed: {e}"
