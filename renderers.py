"""
Renderers — every build_*_html function that assembles the final HTML report.

V8 extracted this module out of the monolithic analyze.py (which was 4 876
lines pre-decomposition). The split is:

    analyze.py    — pipeline orchestration, tier1/tier2 lookup, parsing
    renderers.py  — every HTML / report-rendering function (this file)
    cli.py        — argparse parser + console-script entrypoint
    snp_registry  — single source of truth for SNP metadata
    tests/        — pytest suite (golden snapshots lock byte-equivalence)

The extraction is **deliberately verbatim** — every function below preserves
the exact behaviour it had inside analyze.py, with no whitespace or naming
changes. The V6 module golden snapshots verify byte-level equivalence; if a
snapshot drifts, the extraction was not behaviour-preserving and should be
investigated as a bug, not bumped.

External symbols this module references from `analyze.py` (`build_category_map`,
the `REPORT_VERSION` / `CATEGORY_ORDER` / `APOE_INFO` constants, the
`level_class` reference helper, the local-ancestry / ROH SVG renderers) are
imported lazily inside the functions that need them, so the import graph
stays one-way: ``analyze → renderers`` at runtime, never the reverse at
module-load time.
"""

from __future__ import annotations

import datetime
import math
import re
from pathlib import Path

# Late import — analyze.py finishes loading before any of these renderer
# functions are called at runtime, so circular-import issues do not arise.
from analyze import (
    APOE_INFO,
    CATEGORY_ORDER,
    REPORT_VERSION,
    build_category_map,
)

try:
    from references import level_class
except ImportError:  # references module optional
    def level_class(_) -> str:
        return ""

try:
    from roh import render_ideogram_svg
except ImportError:
    render_ideogram_svg = None  # type: ignore[assignment]

try:
    from local_ancestry import render_chromosome_painting_svg
except ImportError:
    render_chromosome_painting_svg = None  # type: ignore[assignment]


def md_to_html(text: str) -> str:
    """Convert minimal markdown to HTML for embedding in the report."""
    # Code spans
    text = re.sub(r"`(.+?)`", r"<code>\1</code>", text)
    # Bold
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    # Italic
    text = re.sub(r"\*(.+?)\*", r"<em>\1</em>", text)
    # H4 headers (## or ###)
    text = re.sub(
        r"^#{1,4}\s+(.+)$", r"<h4>\1</h4>", text, flags=re.MULTILINE
    )
    # Numbered list items
    text = re.sub(
        r"^\d+\.\s+(.+)$", r"<li>\1</li>", text, flags=re.MULTILINE
    )
    # Bullet list items
    text = re.sub(
        r"^[•\-\*]\s+(.+)$", r"<li>\1</li>", text, flags=re.MULTILINE
    )
    # Wrap consecutive <li> groups in <ul>
    text = re.sub(r"((<li>.*?</li>\n?)+)", r"<ul>\1</ul>", text, flags=re.DOTALL)
    # Paragraphs: lines that aren't already HTML
    lines = text.split("\n")
    out = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("<"):
            out.append(stripped)
        else:
            out.append(f"<p>{stripped}</p>")
    return "\n".join(out)


def _cat_id(cat: str) -> str:
    """Slugify a category name to match the id used on its <details> section."""
    return (
        cat.lower()
        .replace(" ", "-")
        .replace("/", "-")
        .replace("'", "")
        .replace("&", "")
        .replace(".", "")
    )


def _patch_ai_section_in_html(
    html: str, cat_id: str, new_inner_html: str
) -> tuple[str, bool]:
    """Replace the AI Interpretation content for a single category section in
    an existing report. Returns (new_html, success). If the category's section
    doesn't yet contain an ai-section block (e.g. the original run crashed
    before any AI text was rendered), a new one is appended inside the
    <details> block.
    """
    details_re = re.compile(
        r'(<details\b[^>]*\bid="' + re.escape(cat_id) + r'"[^>]*>)'
        r'(.*?)'
        r'(</details>)',
        re.DOTALL,
    )
    m = details_re.search(html)
    if not m:
        return html, False

    opening, body, closing = m.group(1), m.group(2), m.group(3)
    ai_re = re.compile(
        r'(<div class="ai-content">)(.*?)(</div>)',
        re.DOTALL,
    )
    if ai_re.search(body):
        new_body = ai_re.sub(
            lambda mm: mm.group(1) + new_inner_html + mm.group(3),
            body,
            count=1,
        )
    else:
        new_body = body + (
            '<div class="ai-section">'
            '<div class="ai-title">AI Interpretation</div>'
            f'<div class="ai-content">{new_inner_html}</div>'
            '</div>'
        )

    return html[:m.start()] + opening + new_body + closing + html[m.end():], True


def significance_badge(sig: str) -> str:
    classes = {
        "high": "badge-high",
        "moderate": "badge-moderate",
        "low": "badge-low",
    }
    cls = classes.get(sig, "badge-low")
    return f'<span class="badge {cls}">{sig.capitalize()}</span>'


def risk_indicator(copies: int) -> str:
    if copies == 0:
        return '<span class="risk-none">● None</span>'
    if copies == 1:
        return '<span class="risk-one">▲ 1 copy</span>'
    return '<span class="risk-two">▲▲ 2 copies</span>'


# ── Y-DNA HTML builder ────────────────────────────────────────────────────────

# Ancient-DNA reference summaries for major Y-DNA haplogroups. Used as a
# supplemental box in the Y-DNA section.
Y_ANCIENT_DNA_REFS: dict[str, str] = {
    "R1b": (
        "R1b is the dominant Y-DNA haplogroup in modern Western Europe. Ancient "
        "DNA evidence places R1b prominently in the <strong>Yamnaya culture</strong> "
        "(~5,000–4,500 years ago, Pontic-Caspian steppe) and the subsequent "
        "<strong>Bell Beaker culture</strong> (~4,500–4,000 years ago, Western Europe). "
        "In Britain, the Beaker migration replaced ~90% of the prior Neolithic Y-DNA "
        "in just a few centuries — one of the most dramatic genetic turnovers known."
    ),
    "R1a": (
        "R1a is associated with the <strong>Corded Ware culture</strong> (~5,000–4,000 "
        "years ago, Northern/Eastern Europe) and the Indo-Iranian expansion into South "
        "Asia. Ancient DNA from the <strong>Sintashta culture</strong> (~4,000 years ago, "
        "southern Urals) is dominated by R1a-Z93, the South Asian R1a branch."
    ),
    "I1": (
        "I1 represents the Mesolithic European hunter-gatherer Y-DNA legacy that "
        "survived through Neolithic and Bronze Age replacements, particularly in "
        "Scandinavia. Ancient I1 has been found in Mesolithic Scandinavian samples "
        "and is enriched in modern Nordic populations."
    ),
    "I2": (
        "I2 is the other major Mesolithic European Y-DNA lineage, particularly I2a, "
        "which is concentrated in the Western Balkans today. Cheddar Man (~10,000 "
        "years ago, Mesolithic Britain) carried I2a — making him a direct ancestral "
        "line representative for modern I2 carriers."
    ),
    "G2a": (
        "G2a was the dominant Y-DNA of <strong>Neolithic European farmers</strong> "
        "expanding from Anatolia ~8,000 years ago. Ötzi the Iceman (~5,300 years ago, "
        "Alpine Copper Age) carried G2a2b. G2a was largely replaced by R1b/R1a in "
        "the Bronze Age but persists at moderate frequencies in Sardinia, the "
        "Caucasus, and the Mediterranean."
    ),
    "N": (
        "N is the dominant Y-DNA across the Uralic-speaking peoples of northern "
        "Eurasia — Finns, Saami, Estonians, Hungarians, Yakuts, Nenets. N1c "
        "(N-M178) is found in ~60% of Finnish men. The expansion is associated "
        "with Bronze Age and Iron Age migrations across the northern forest belt."
    ),
    "Q": (
        "Q is the principal Y-DNA of the peopling of the Americas, carried across "
        "Beringia ~15,000–20,000 years ago by ancestors of Native Americans. "
        "Q-M242 and downstream Q-M3 dominate Native American Y-DNA. Q is also "
        "found in Central Asian (Kazakhs, some Mongolian groups) and small "
        "European populations (Ashkenazi Q1b)."
    ),
    "O": (
        "O is the dominant Y-DNA of East Asia and Southeast Asia. Ancient DNA "
        "shows O lineages expanding with the Neolithic agricultural transitions "
        "in the Yellow River and Yangtze River basins (~8,000–6,000 years ago) "
        "and subsequent expansions throughout Southeast Asia and Oceania."
    ),
    "T": (
        "T is rare in Europe today but has been identified in <strong>Neolithic Linear "
        "Pottery (LBK) and Cardial samples</strong> from central Europe and Iberia "
        "(~7,000 years ago), suggesting Neolithic farmer association. T is more "
        "common today in East Africa (Cushitic-speakers, Ethiopian Jews), the "
        "Levant, and parts of the Mediterranean."
    ),
}


def build_ydna_ancient_block(terminal_hg: str, path: list[dict]) -> str:
    """Return a small ancient-DNA reference box for the terminal haplogroup,
    falling back through ancestors if no entry matches the leaf."""
    # First try the terminal haplogroup, then walk up the path
    candidates = [terminal_hg] + [n.get("haplogroup", "") for n in reversed(path or [])]
    seen = set()
    for hg in candidates:
        base = hg.split("-")[0] if hg else ""
        if base in seen:
            continue
        seen.add(base)
        if base in Y_ANCIENT_DNA_REFS:
            return (
                '<div class="ydna-ancient">'
                '<div class="ydna-ancient-title">Ancient-DNA Comparisons</div>'
                f"<p>{Y_ANCIENT_DNA_REFS[base]}</p>"
                "</div>"
            )
    return ""


def build_ydna_html(y_result: dict) -> str:
    """Render the Y-DNA haplogroup section for the HTML report."""
    status = y_result.get("status", "no_y_data")
    path = y_result.get("path", [])
    terminal = y_result.get("terminal_haplogroup", "Unknown")
    migration = y_result.get("terminal_migration") or ""
    further = y_result.get("further_testing") or ""
    chip_gaps = y_result.get("chip_gaps", [])
    not_tested = y_result.get("not_tested_branches", [])
    y_count = y_result.get("y_snp_count", 0)

    # ── Breadcrumb path visualization ──
    crumbs = ""
    for i, node in enumerate(path):
        confirmed = node.get("snp_status") == "confirmed"
        gt_txt = node.get("found_genotype", "")
        tooltip = (
            f"{node['snp_name']}: {gt_txt} ({'confirmed' if confirmed else 'inferred — not on chip'})"
        )
        cls = "crumb-ok" if confirmed else "crumb-gap"
        arrow = "" if i == 0 else '<span class="crumb-arrow">›</span>'
        crumbs += (
            f'{arrow}<span class="{cls}" title="{tooltip}">'
            f"{node['haplogroup']}"
            f'<span class="crumb-snp">{node["snp_name"]}</span>'
            f"</span>"
        )

    # ── Per-node marker table ──
    rows = ""
    for node in path:
        confirmed = node.get("snp_status") == "confirmed"
        gt = node.get("found_genotype", "–")
        snp_name = node["snp_name"]
        rsids = node.get("rsids", [])
        rsid_links = ", ".join(
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{r}" target="_blank" rel="noopener">{r}</a>'
            for r in rsids
        ) if rsids else "position-based / not assigned"
        pos = node.get("pos")
        pos_txt = f"chrY:{pos:,}" if pos else "—"
        status_badge = (
            '<span class="badge badge-low">Confirmed</span>'
            if confirmed else
            '<span class="badge badge-moderate">Inferred (not on chip)</span>'
        )
        rows += (
            f"<tr>"
            f'<td class="gt-cell">{node["haplogroup"]}</td>'
            f'<td class="rsid-cell">{snp_name}</td>'
            f"<td>{rsid_links}</td>"
            f"<td>{pos_txt}</td>"
            f'<td class="gt-cell">{gt}</td>'
            f"<td>{status_badge}</td>"
            f'<td class="sum-cell">{node.get("description", "")}</td>'
            f"</tr>\n"
        )

    marker_table = f"""
<div class="tbl-wrap" style="margin-top:16px">
<table class="snp-tbl">
<thead><tr>
  <th>Haplogroup</th><th>Marker</th><th>rsID(s)</th>
  <th>Chr Y Position</th><th>Genotype</th><th>Status</th><th>Description</th>
</tr></thead>
<tbody>{rows}</tbody>
</table>
</div>""" if rows else ""

    # ── Chip gap notice ──
    gap_html = ""
    if chip_gaps:
        gap_list = ", ".join(f"<code>{g}</code>" for g in chip_gaps)
        gap_html = (
            f'<div class="ydna-gap">'
            f"<strong>Resolution limit:</strong> The following downstream markers "
            f"were not detected on your chip: {gap_list}. "
            f"Further refinement requires additional testing (see below)."
            f"</div>"
        )

    # ── Not-tested branches ──
    branch_html = ""
    if not_tested:
        items = "".join(
            f"<li><strong>{b['haplogroup']}</strong> ({b['snp_name']}): "
            f"{b['description']}. {b.get('further','')}</li>"
            for b in not_tested
        )
        branch_html = (
            f'<div class="ydna-branches">'
            f"<strong>Unresolved downstream branches</strong> "
            f"(markers not on chip):<ul>{items}</ul></div>"
        )

    # ── Migration narrative ──
    migration_html = ""
    if migration:
        migration_html = (
            f'<div class="ydna-migration">'
            f'<div class="ydna-migration-title">Migration History</div>'
            f"<p>{migration}</p>"
            f"</div>"
        )

    # ── Ancient-DNA reference box ──
    ancient_html = build_ydna_ancient_block(terminal, path)

    # ── Further testing ──
    further_html = ""
    if further:
        # wrap newlines as paragraphs / bullet points
        further_fmt = re.sub(r"\n+•", "<br>•", further)
        further_html = (
            f'<div class="ydna-further">'
            f'<strong>Further testing (FTDNA Big Y-700 / dedicated panels):</strong>'
            f"<p>{further_fmt}</p>"
            f"</div>"
        )

    # ── Status badge for the section header ──
    if status == "no_y_data":
        header_badge = '<span class="ydna-badge ydna-badge-gray">No Y Data</span>'
    elif status == "resolved":
        header_badge = '<span class="ydna-badge ydna-badge-green">Resolved</span>'
    elif status == "not_k":
        header_badge = '<span class="ydna-badge ydna-badge-blue">Pre-K Haplogroup</span>'
    else:
        header_badge = '<span class="ydna-badge ydna-badge-amber">Partial (chip gaps)</span>'

    y_count_note = f" &nbsp;·&nbsp; {y_count:,} Y-chromosome SNPs in file" if y_count else ""
    conf_badge = _confidence_badge(y_result.get("confidence"), y_result.get("confidence_note"))
    coverage_note = ""
    if y_result.get("n_markers_on_path"):
        coverage_note = (
            f' &nbsp;&middot;&nbsp; {y_result.get("n_markers_confirmed", 0)}/'
            f'{y_result.get("n_markers_on_path")} path markers confirmed'
        )

    return f"""
<section class="ydna-section" id="y-haplogroup">
<h2>Y-DNA Haplogroup {header_badge}</h2>

<div class="ydna-path-wrap">
  <div class="ydna-path-label">Haplogroup path</div>
  <div class="ydna-crumbs">{crumbs if crumbs else "<em>Could not determine path</em>"}</div>
  <div class="ydna-path-note">
    Terminal haplogroup: <strong>{terminal}</strong>{y_count_note}{coverage_note}
    &nbsp;&middot;&nbsp;
    <span class="crumb-ok" style="padding:2px 6px">&#9679; confirmed SNP</span>
    <span class="crumb-gap" style="padding:2px 6px">&#9675; inferred (not on chip)</span>
  </div>
  {conf_badge}
</div>

{marker_table}
{gap_html}
{branch_html}
{migration_html}
{ancient_html}
{further_html}
<div class="ydna-disclaimer">
<strong>Informational use only.</strong> Y-haplogroup assignment from a
genotyping array is limited by which markers the chip carries; terminal
(most-specific) subclades are frequently inferred across untyped positions
rather than directly confirmed. This is a genealogical/ancestral indicator,
not a medical or identity test.
</div>
</section>
"""


# ── mtDNA HTML builder ────────────────────────────────────────────────────────

def build_mtdna_html(mt_result: dict) -> str:
    """Render the mtDNA haplogroup section for the HTML report."""
    status = mt_result.get("status", "no_data")
    haplogroup = mt_result.get("haplogroup", "Unknown")
    confidence = mt_result.get("confidence", "low")
    matched = mt_result.get("matched_markers", [])
    migration = mt_result.get("migration", "")
    ancient_dna = mt_result.get("ancient_dna", "")
    further = mt_result.get("further_testing", "")
    n_mt_snps = mt_result.get("mt_snp_count", 0)
    path = mt_result.get("path", [])

    if status == "no_data":
        return (
            '<section class="mtdna-section" id="mt-haplogroup">'
            '<h2>mtDNA Haplogroup <span class="mtdna-badge mtdna-badge-gray">No mtDNA Data</span></h2>'
            "<p>No mitochondrial-DNA SNPs were detected on this chip in sufficient "
            "numbers to call a haplogroup. Most consumer autosomal arrays include "
            "only a handful of mtDNA markers, which is rarely enough for confident "
            "mtDNA classification. For maternal-lineage analysis, FTDNA mtFull "
            "(full mitochondrial sequence) or 23andMe (which uses a curated mtDNA "
            "subset and reports a haplogroup directly) is the appropriate test.</p>"
            "</section>"
        )

    conf_badge_cls = "mtdna-badge-pink" if confidence == "high" else "mtdna-badge-gray"
    conf_label = f"{haplogroup} ({confidence} confidence)"

    matched_rows = ""
    for m in matched[:15]:
        matched_rows += (
            f"<tr>"
            f'<td class="rsid-cell">'
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{m.get("rsid","")}" '
            f'target="_blank" rel="noopener">{m.get("rsid","")}</a></td>'
            f"<td>{m.get('haplogroup_marker','')}</td>"
            f'<td class="gt-cell">{m.get("genotype","")}</td>'
            f"<td>{m.get('description','')}</td>"
            f"</tr>\n"
        )
    matched_table = ""
    if matched_rows:
        matched_table = (
            '<div class="tbl-wrap" style="margin-top:14px"><table class="snp-tbl">'
            "<thead><tr>"
            "<th>rsID</th><th>Marker</th><th>Genotype</th><th>Haplogroup association</th>"
            "</tr></thead>"
            f"<tbody>{matched_rows}</tbody></table></div>"
        )

    migration_html = ""
    if migration:
        migration_html = (
            '<div class="mtdna-migration">'
            '<div class="mtdna-migration-title">Maternal-Line Migration History</div>'
            f"<p>{migration}</p>"
            "</div>"
        )

    ancient_html = ""
    if ancient_dna:
        ancient_html = (
            '<div class="mtdna-migration" style="background:linear-gradient(135deg,rgba(63,185,80,.06),rgba(88,166,255,.05));border-color:rgba(63,185,80,.25)">'
            '<div class="mtdna-migration-title" style="color:var(--grn)">Ancient-DNA Comparisons</div>'
            f"<p>{ancient_dna}</p>"
            "</div>"
        )

    further_html = ""
    if further:
        further_html = (
            '<div class="ydna-further">'
            "<strong>Further testing (FTDNA mtFull, full mitochondrial sequencing):</strong>"
            f"<p>{further}</p>"
            "</div>"
        )

    mt_count_note = f" &nbsp;·&nbsp; {n_mt_snps:,} mtDNA SNPs available on this chip" if n_mt_snps else ""

    # ── Maternal lineage chain ──
    # A haplogroup is the end of a chain of branch points, not a bare label.
    # The Y-DNA section has always shown that chain; this shows the maternal
    # one on the same terms, including which links are confirmed by a marker
    # on the chip and which are inferred because no marker for them was typed.
    crumbs = ""
    for i, node in enumerate(path):
        confirmed = node.get("snp_status") == "confirmed"
        ev = "; ".join(node.get("evidence") or []) or "no marker for this branch on this chip"
        tip = f"{node.get('snp_name','')}: {ev}"
        cls = "crumb-ok" if confirmed else "crumb-gap"
        arrow = "" if i == 0 else '<span class="crumb-arrow">&rsaquo;</span>'
        snp = (f'<span class="crumb-snp">{_esc(node["snp_name"])}</span>'
               if node.get("markers") is not False and node.get("snp_name")
               and node["snp_name"] != node["haplogroup"] else "")
        crumbs += (f'{arrow}<span class="{cls}" title="{_esc(tip)}">'
                   f'{_esc(node["haplogroup"])}{snp}</span>')
    n_conf = mt_result.get("n_confirmed_branch_points", 0)
    n_bp = mt_result.get("n_branch_points", 0)
    gaps = mt_result.get("chip_gaps") or []
    contradictions = mt_result.get("contradictions") or []
    mt_crumbs = ""
    if crumbs:
        gap_note = ""
        if gaps:
            gap_note = (f'<div class="mtdna-evidence" style="margin-top:5px">'
                        f'Greyed links ({", ".join(_esc(g) for g in gaps)}) are '
                        f'inferred rather than observed — no marker defining '
                        f'them is typed on this chip, so they are implied by a '
                        f'confirmed branch further down.</div>')
        contra_note = ""
        if contradictions:
            c = contradictions[0]
            contra_note = (f'<div class="mtdna-evidence" style="margin-top:5px;'
                           f'color:#b03a2e"><strong>Conflicting markers at '
                           f'{_esc(c["at"])}:</strong> '
                           f'{", ".join(_esc(b) for b in c["branches"])} are '
                           f'mutually exclusive branches and cannot both be '
                           f'yours. This usually means a mis-called genotype '
                           f'rather than an interesting result.</div>')
        mt_crumbs = (f'<div class="crumbs" style="margin:8px 0">{crumbs}</div>'
                     f'<div class="mtdna-evidence">Your maternal line, root to '
                     f'branch tip &mdash; {n_conf} of {n_bp} branch points '
                     f'confirmed by a marker on this chip.</div>'
                     f'{gap_note}{contra_note}')

    return f"""
<section class="mtdna-section" id="mt-haplogroup">
<h2>mtDNA Haplogroup <span class="mtdna-badge {conf_badge_cls}">{conf_label}</span></h2>

<div class="mtdna-result">
  <div class="mtdna-call">{haplogroup}</div>
  {mt_crumbs}
  <div class="mtdna-evidence">
    Based on {mt_result.get('n_markers_derived', 0)} derived /
    {mt_result.get('n_markers_matched', len(matched))} matched markers
    (of {mt_result.get('n_markers_expected', '—')} on this panel){mt_count_note}.
    Note: mtDNA haplogroup calls from autosomal-chip data are approximate — full
    mitochondrial sequencing remains the gold standard.
  </div>
  {_confidence_badge(confidence, mt_result.get('confidence_note',''))}
</div>

{matched_table}
{migration_html}
{ancient_html}
{further_html}
</section>
"""


# ── Professional-grade section renderers ──────────────────────────────────────

def _esc(s: str) -> str:
    """Minimal HTML escape for safety in injected strings."""
    if s is None:
        return ""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;"))


_CONF_LABELS = {
    "high": "High confidence",
    "moderate": "Moderate confidence",
    "low": "Low confidence",
    "none": "Indeterminate",
    "n/a": "Not applicable",
}


def _confidence_badge(confidence: str | None, note: str = "") -> str:
    """Render a consistent confidence pill (+ optional explanatory note) used
    across every score section. Returns '' when no confidence is supplied."""
    if not confidence:
        return ""
    key = str(confidence).lower()
    label = _CONF_LABELS.get(key, f"{confidence} confidence")
    note_html = f'<span class="conf-note">{_esc(note)}</span>' if note else ""
    return (
        f'<div class="conf-row">'
        f'<span class="conf-badge conf-{_esc(key).replace("/", "")}">{_esc(label)}</span>'
        f'{note_html}'
        f'</div>'
    )


def build_qc_html(qc: dict | None) -> str:
    if not qc:
        return ""
    domain_rows = ""
    for domain, info in qc["domain_coverage"].items():
        bar_pct = info["pct"]
        bar_class = "qc-bar-full" if bar_pct == 100 else (
            "qc-bar-good" if bar_pct >= 75 else (
                "qc-bar-fair" if bar_pct >= 50 else "qc-bar-poor"))
        called = info["called"]
        total = info["total"]
        var_chips = "".join(
            f'<span class="qc-var-{"on" if v["called"] else "off"}" '
            f'title="{_esc(v["rsid"])} — {_esc(v["gene"])}">{_esc(v["gene"])}</span>'
            for v in info["variants"]
        )
        domain_rows += (
            f'<tr>'
            f'<td class="qc-domain-name">{_esc(domain)}</td>'
            f'<td class="qc-cov-cell">'
            f'  <div class="qc-bar-wrap"><div class="qc-bar {bar_class}" '
            f'style="width:{bar_pct}%"></div></div>'
            f'  <span class="qc-bar-label">{called}/{total} ({bar_pct}%)</span>'
            f'</td>'
            f'<td class="qc-chips">{var_chips}</td>'
            f'</tr>\n'
        )

    return f"""
<section class="qc-section" id="quality-control">
<h2>Data Quality & Callability</h2>
<div class="qc-grade-card">
  <div class="qc-grade {qc['grade_class']}">{_esc(qc['grade'])}</div>
  <div class="qc-grade-text">
    <div class="qc-grade-note">{_esc(qc['grade_note'])}</div>
    <div class="qc-grade-stats">
      Total SNPs: <strong>{qc['total_snps']:,}</strong> &nbsp;·&nbsp;
      Callability: <strong>{qc['callability_pct']}%</strong> &nbsp;·&nbsp;
      {f"No-calls: <strong>{qc['no_call_count']:,}</strong> &nbsp;·&nbsp;" if qc.get('no_call_count') is not None else ''}
      Chip: <strong>{_esc(qc['file_format'])}</strong> &nbsp;·&nbsp;
      Inferred sex: <strong>{_esc(qc['inferred_sex'])}</strong>
      {f"&nbsp;·&nbsp; Avg domain callability: <strong>{qc['average_domain_callability']}%</strong>" if qc.get('average_domain_callability') is not None else ''}
    </div>
    <div class="qc-grade-stats">
      Autosomal: {qc['autosomal_count']:,} &nbsp;·&nbsp;
      X: {qc['x_count']:,} &nbsp;·&nbsp;
      Y: {qc['y_count']:,} &nbsp;·&nbsp;
      mtDNA: {qc['mt_count']:,}
    </div>
    <div class="qc-grade-stats">
      Database matches: <strong>{qc['tier1_match_count']}</strong> /
      {qc['db_size']} ({qc['tier1_match_pct']}%) &nbsp;·&nbsp;
      File hash: <code>{_esc(qc['file_hash'])}</code>
    </div>
  </div>
</div>
<h3 class="qc-h3">Per-Domain Variant Coverage</h3>
<div class="tbl-wrap">
<table class="qc-tbl"><thead><tr>
<th>Clinical Domain</th><th>Coverage</th><th>Variants</th>
</tr></thead><tbody>{domain_rows}</tbody></table>
</div>
</section>
"""


def build_prs_html(prs: dict | None) -> str:
    if not prs or not prs.get("panels"):
        return ""

    panel_cards = ""
    for name, panel in prs["panels"].items():
        result = panel["result"]
        status = result.get("status")
        if status == "not_applicable":
            panel_cards += (
                f'<div class="prs-card prs-card-na">'
                f'<div class="prs-name">{_esc(name)}</div>'
                f'<div class="prs-na">{_esc(result["reason"])}</div>'
                f'</div>'
            )
            continue
        if status == "insufficient_data":
            panel_cards += (
                f'<div class="prs-card prs-card-na">'
                f'<div class="prs-name">{_esc(name)}</div>'
                f'{_confidence_badge(result.get("confidence", "none"))}'
                f'<div class="prs-na">{_esc(result["reason"])}</div>'
                f'</div>'
            )
            continue
        tier = result["tier"]
        tier_class = result["tier_class"]
        pct = result["percentile"]
        z = result["z_score"]
        callability = result["callability"]
        confidence = result.get("confidence")
        conf_badge = _confidence_badge(confidence, result.get("confidence_note"))
        low_warn = ""
        if confidence == "low":
            low_warn = (
                '<div class="prs-lowconf">⚠️ Low coverage — this percentile is '
                'unreliable and is shown for transparency only.</div>'
            )
        n_used = result.get("n_used", len(result["used"]))
        n_exp = result.get("n_expected", len(result["used"]) + len(result["missing"]))
        details = (
            f'<div class="prs-details">'
            f'<div><span class="prs-lab">Z-score:</span> <strong>{z:+.2f}</strong></div>'
            f'<div><span class="prs-lab">Percentile:</span> <strong>{pct:.0f}th</strong></div>'
            f'<div><span class="prs-lab">Variants typed:</span> '
            f'{n_used}/{n_exp} '
            f'({callability}%)</div>'
            f'</div>'
        )
        # percentile track visual (0–100)
        track = (
            f'<div class="prs-track">'
            f'<div class="prs-track-fill" style="width:{pct}%"></div>'
            f'<div class="prs-track-pointer" style="left:{pct}%"></div>'
            f'</div>'
        )
        # Per-variant breakdown: which SNPs actually drove this score, with the
        # scored genotype, dosage, and published per-allele effect. All values
        # come straight from the panel definition + the sample's genotype —
        # nothing synthesised. Previously collapsed to a bare count.
        variants_detail = ""
        used = result.get("used", [])
        if used:
            vrows = ""
            for u in sorted(used, key=lambda x: -abs(x.get("log_or", 0.0))):
                lor = u.get("log_or", 0.0)
                orr = math.exp(lor)
                af = u.get("af")
                af_txt = f"{af:.0%}" if isinstance(af, (int, float)) else "—"
                vrows += (
                    f'<tr><td class="rsid-cell">'
                    f'<a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(u.get("rsid",""))}" '
                    f'target="_blank" rel="noopener">{_esc(u.get("rsid",""))}</a></td>'
                    f'<td><strong>{_esc(u.get("gene",""))}</strong></td>'
                    f'<td>{_esc(u.get("effect_allele",""))}</td>'
                    f'<td class="gt-cell">{_esc(u.get("genotype",""))}</td>'
                    f'<td>{u.get("dosage",0)}</td>'
                    f'<td>{orr:.2f}</td>'
                    f'<td>{af_txt}</td></tr>'
                )
            missing = result.get("missing", [])
            missing_note = ""
            if missing:
                mnames = ", ".join(
                    _esc(m.get("rsid", "")) for m in missing[:20] if m.get("rsid")
                )
                extra = f" +{len(missing) - 20} more" if len(missing) > 20 else ""
                missing_note = (
                    f'<div class="prs-missing">{len(missing)} panel variant(s) '
                    f"not typed on this chip: {mnames}{extra}. These are excluded "
                    f"from the score entirely (both the raw score and its expected "
                    f"mean/variance are computed over typed variants only), so a "
                    f"score built on fewer variants is a less complete estimate — "
                    f"see the coverage/confidence indicator above.</div>"
                )
            variants_detail = (
                f'<details class="prs-variants"><summary>Contributing variants '
                f"({len(used)})</summary>"
                f'<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>'
                f"<th>rsID</th><th>Gene</th><th>Effect allele</th>"
                f"<th>Your genotype</th><th>Copies</th><th>Per-allele OR</th>"
                f"<th>Effect-allele freq</th></tr></thead><tbody>{vrows}</tbody>"
                f"</table></div>{missing_note}</details>"
            )
        panel_cards += f"""
<div class="prs-card">
  <div class="prs-card-head">
    <div class="prs-name">{_esc(name)} <span class="prs-short">{_esc(panel.get("trait_short",""))}</span></div>
    <div class="prs-tier {tier_class}">{_esc(tier)}</div>
  </div>
  <div class="prs-desc">{_esc(panel["description"])}</div>
  {track}
  {conf_badge}
  {low_warn}
  {details}
  {variants_detail}
  <div class="prs-context">
    <strong>Population lifetime risk:</strong> {_esc(panel["population_lifetime_risk"])}<br>
    <strong>High-tier implication:</strong> {_esc(panel["high_tier_implication"])}
  </div>
  <div class="prs-ref"><em>Source: {_esc(panel["reference"])}</em></div>
</div>
"""

    headline = ""
    if prs.get("headline_findings"):
        items = "".join(
            f'<li><strong>{_esc(n)}</strong> — '
            f'<span class="prs-tier {p["result"]["tier_class"]}">'
            f'{_esc(p["result"]["tier"])}</span> '
            f'({p["result"]["percentile"]:.0f}th percentile)</li>'
            for n, p in prs["headline_findings"]
        )
        headline = (
            f'<div class="prs-headline">'
            f'<strong>Elevated / High polygenic findings:</strong>'
            f'<ul>{items}</ul></div>'
        )

    return f"""
<section class="prs-section" id="polygenic-risk-scores">
<h2>Polygenic Risk Scores <span class="pro-pill">PRO</span></h2>
<p class="prs-intro">
Curated-variant polygenic risk scores derived from genome-wide significant,
replicated GWAS hits. Each score sums log(OR) across called variants, is
standardised against an expected Hardy-Weinberg distribution (European
reference allele frequencies), and reported as a Z-score and percentile.
Tier classification: Low (&lt;5th) · Below Average (5–20th) · Average
(20–80th) · Elevated (80–95th) · High (≥95th).
</p>
<div class="prs-caveat">
<strong>Important caveats:</strong> these are <em>curated-variant</em> scores,
not full clinical-grade PGS (which use 10⁴–10⁶ SNPs). Effect magnitudes are
likely under-estimates of full PGS. Reference distribution is European —
ancestry-mismatched interpretation is less reliable. Polygenic risk is one
input alongside family history, lifestyle, and established clinical risk
factors; not diagnostic.
</div>
{headline}
<div class="prs-grid">{panel_cards}</div>
</section>
"""


def _build_novelty_panels_html(panels: list[dict]) -> str:
    if not panels:
        return ""
    section_blocks = ""
    for panel in panels:
        rows = ""
        for g in panel["genes"]:
            if g["tested"]:
                status = (
                    f'<span class="pgx-pheno pheno-nm">{_esc(g["genotype"])}</span>'
                )
            else:
                status = '<span class="pgx-na">Not tested on this chip</span>'
            rsid_cell = _esc(g["rsid"]) if g.get("rsid") else "—"
            rows += (
                f'<tr>'
                f'<td><strong>{_esc(g["gene"])}</strong></td>'
                f'<td class="rsid-cell">{rsid_cell}</td>'
                f'<td>{status}</td>'
                f'<td>{_esc(g["note"])}</td>'
                f'</tr>'
            )
        section_blocks += f"""
<div class="pgx-card">
  <div class="pgx-card-head">
    <div class="pgx-gene-name">{_esc(panel["section"])}</div>
    <div class="pgx-gene-long">{_esc(panel["description"])}</div>
  </div>
  <div class="tbl-wrap" style="margin-top:8px">
    <table class="pgx-tbl">
      <thead><tr><th>Gene</th><th>rsID</th><th>Genotype</th><th>Note</th></tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
</div>
"""
    return f"""
<div class="pgx-novelty" style="margin-top:24px">
  <h3>Exploratory Gene Panels <span style="font-size:0.75em;opacity:0.7">(for fun — non-clinical)</span></h3>
  <div class="pgx-caveat">
    <strong>Not a clinical finding.</strong> These panels are exploratory.
    Tested genotypes are shown without effect-size interpretation; rows
    marked "Not tested on this chip" lack a probe in the source data.
    Some symbols submitted to earlier versions of these panels were dropped
    because they were not HGNC-approved and had no validated assay.
    Do not use for medical decisions.
  </div>
  <div class="pgx-grid">{section_blocks}</div>
</div>
"""


def build_pgx_html(pgx: dict | None) -> str:
    if not pgx:
        return ""

    # Per-gene phenotype cards
    gene_cards = ""
    for gene_name, result in pgx["per_gene"].items():
        # Variant detail table
        var_rows = ""
        for v in result["variant_calls"]:
            if v["called"]:
                applied = v["applied_impact"]
                applied_txt = f"{applied:+.1f}" if applied != 0 else "0"
                var_rows += (
                    f'<tr>'
                    f'<td class="rsid-cell">{_esc(v["rsid"])}</td>'
                    f'<td><strong>{_esc(v["star_allele"])}</strong></td>'
                    f'<td>{_esc(v["name"])}</td>'
                    f'<td class="gt-cell">{_esc(v.get("genotype",""))}</td>'
                    f'<td>{v.get("dosage", "—")}</td>'
                    f'<td class="pgx-impact">{applied_txt}</td>'
                    f'</tr>'
                )
            else:
                var_rows += (
                    f'<tr class="pgx-uncalled">'
                    f'<td class="rsid-cell">{_esc(v["rsid"])}</td>'
                    f'<td>{_esc(v["star_allele"])}</td>'
                    f'<td>{_esc(v["name"])}</td>'
                    f'<td colspan="3" class="pgx-na">not on chip</td>'
                    f'</tr>'
                )
        if result.get("indeterminate"):
            # No defining variants typed — do NOT imply a phenotype/activity.
            score_row = (
                f'<div class="pgx-pheno {result["phenotype_class"]}">'
                f'{_esc(result["phenotype"])}</div>'
            )
        elif result.get("is_binary"):
            score_row = (
                f'<div class="pgx-score">'
                f'<span class="pgx-pheno {result["phenotype_class"]}">'
                f'{_esc(result["phenotype"])}</span>'
                f'</div>'
            )
        else:
            activity = result["activity_score"]
            baseline = result["baseline_activity"]
            score_row = (
                f'<div class="pgx-score">'
                f'<span class="pgx-activity-label">Activity score:</span> '
                f'<span class="pgx-activity">{activity}</span> '
                f'<span class="pgx-activity-base">(baseline {baseline})</span>'
                f'</div>'
                f'<div class="pgx-pheno {result["phenotype_class"]}">'
                f'{_esc(result["phenotype"])}</div>'
            )
        conf_badge = _confidence_badge(result.get("confidence"), result.get("confidence_note"))
        drug_rows = ""
        for drug in result["drug_recs"]:
            recs_inline = ""
            for k, v in drug.items():
                if k == "drug":
                    continue
                cls = ""
                if k == result["phenotype_code"]:
                    cls = "pgx-rec-active"
                recs_inline += (
                    f'<div class="pgx-rec-line {cls}">'
                    f'<span class="pgx-rec-pheno">{_esc(k)}</span> '
                    f'<span class="pgx-rec-text">{_esc(v)}</span>'
                    f'</div>'
                )
            drug_rows += (
                f'<div class="pgx-drug-row">'
                f'<div class="pgx-drug-name">{_esc(drug["drug"])}</div>'
                f'{recs_inline}'
                f'</div>'
            )
        um_caveat = ""
        if result.get("um_caveat"):
            um_caveat = f'<div class="pgx-caveat">{_esc(result["um_caveat"])}</div>'

        gene_cards += f"""
<div class="pgx-card">
  <div class="pgx-card-head">
    <div class="pgx-gene-name">{_esc(gene_name)}</div>
    <div class="pgx-gene-long">{_esc(result["long_name"])}</div>
  </div>
  {score_row}
  {conf_badge}
  <div class="pgx-callability">
    Defining variants typed: {result['callable_variants']}/{result['total_variants']}
    ({result['callability_pct']}%)
  </div>
  <details class="pgx-details">
    <summary>Variant detail</summary>
    <div class="tbl-wrap" style="margin-top:8px">
      <table class="pgx-tbl">
        <thead><tr><th>rsID</th><th>Star</th><th>Variant</th><th>GT</th>
        <th>Dose</th><th>Δ activity</th></tr></thead>
        <tbody>{var_rows}</tbody>
      </table>
    </div>
  </details>
  <div class="pgx-drugs">
    <div class="pgx-drugs-title">Drug-by-drug recommendations</div>
    {drug_rows}
  </div>
  {um_caveat}
  <div class="pgx-guideline"><em>{_esc(result["cpic_guideline"])}</em></div>
</div>
"""

    actionable = ""
    if pgx.get("actionable_findings"):
        items = ""
        for a in pgx["actionable_findings"]:
            items += (
                f'<li class="actionable-item">'
                f'<div class="actionable-head">'
                f'<span class="pgx-pheno {a["phenotype_class"]}">{_esc(a["phenotype"])}</span> '
                f'<span class="actionable-gene">{_esc(a["gene"])}</span> &middot; '
                f'<span class="actionable-drug">{_esc(a["drug"])}</span></div>'
                f'<div class="actionable-rec">{_esc(a["recommendation"])}</div>'
                f'</li>'
            )
        actionable = (
            f'<div class="pgx-actionable">'
            f'<h3>Actionable Drug-Gene Findings ({len(pgx["actionable_findings"])})</h3>'
            f'<ul class="actionable-list">{items}</ul>'
            f'</div>'
        )

    return f"""
<section class="pgx-section" id="pharmacogenomics">
<h2>Pharmacogenomic Phenotypes <span class="pro-pill">PRO</span></h2>
<p class="pgx-intro">
Genotype-to-phenotype calls for {pgx['n_genes_tested']} clinically-actionable
genes, using the <strong>CPIC activity-score model</strong>. Each gene's
metabolizer phenotype (Poor / Intermediate / Normal / Rapid / Ultra-rapid)
drives drug-specific dosing recommendations grounded in published CPIC
guidelines.
</p>
<div class="pgx-caveat">
<strong>Caveats:</strong> Star-allele calling from SNP arrays alone is
approximate — true clinical PGx panels also detect gene duplications/
deletions (most relevant for CYP2D6 ultra-rapid metabolizers) and rare
alleles. HLA-B*57:01 calls from proxy SNPs MUST be confirmed by direct
HLA typing before any abacavir prescription. Phenotype-to-drug mapping
follows CPIC guideline structure but is not a substitute for prescriber
judgment.
</div>
{actionable}
<div class="pgx-grid">{gene_cards}</div>
{_build_drug_database_html(pgx.get("database_findings") or [])}
{_build_novelty_panels_html(pgx.get("novelty_panels") or [])}
</section>
"""


def _build_drug_database_html(findings: list[dict]) -> str:
    if not findings:
        return ""
    rows = ""
    for f in findings:
        snp_bits = ", ".join(
            f'{_esc(m["rsid"])}={_esc(m["genotype"])}' for m in f.get("matched_snps", [])
        )
        genes = ", ".join(_esc(g) for g in f.get("genes", []))
        phenos = ", ".join(_esc(p) for p in f.get("phenotypes", []))
        pheno_html = f'<div class="db-drug-phenos"><em>Relevant phenotypes:</em> {phenos}</div>' if phenos else ""
        rows += (
            f'<li class="db-drug-item">'
            f'<div class="db-drug-head">'
            f'<span class="db-drug-name">{_esc(f["drug"])}</span> '
            f'<span class="db-drug-gene">{genes}</span> '
            f'<span class="db-drug-match">{f["n_matched"]}/{f["n_markers"]} markers</span>'
            f'</div>'
            f'<div class="db-drug-snps">{snp_bits}</div>'
            f'{pheno_html}'
            f'<div class="db-drug-rec">{_esc(f.get("recommendation", ""))}</div>'
            f'</li>'
        )
    return (
        f'<div class="pgx-drug-database">'
        f'<h3>Drug-Database Matches ({len(findings)})</h3>'
        f'<p class="db-intro">Drugs in the curated database where your genotype '
        f'overlaps the listed pharmacogenomic SNP markers.</p>'
        f'<ul class="db-drug-list">{rows}</ul>'
        f'</div>'
    )


def build_interactions_html(inter: dict | None) -> str:
    if not inter or not inter.get("findings"):
        return ""
    items = ""
    for f in inter["findings"]:
        sev_class = (
            "sev-high" if f["severity"] == "high" else
            "sev-mod" if f["severity"] == "moderate" else
            "sev-low"
        )
        var_chips = "".join(
            f'<code class="inter-var">{_esc(v)}</code>'
            for v in f.get("variants", [])
        )
        items += f"""
<div class="inter-finding">
  <div class="inter-head">
    <span class="inter-sev {sev_class}">{_esc(f["severity"].upper())}</span>
    <span class="inter-title">{_esc(f["title"])}</span>
  </div>
  <div class="inter-vars">{var_chips}</div>
  <div class="inter-body">
    <div class="inter-interp"><strong>Interpretation:</strong> {_esc(f["interpretation"])}</div>
    <div class="inter-action"><strong>Action:</strong> {_esc(f["action"])}</div>
  </div>
</div>
"""
    return f"""
<section class="inter-section" id="variant-interactions">
<h2>Variant Interactions & Compound Heterozygosity <span class="pro-pill">PRO</span></h2>
<p class="inter-intro">
Multi-variant patterns that change risk qualitatively beyond what any
single variant implies — the clinical interpretations most often missed
when SNPs are reported one at a time.
</p>
<div class="inter-stats">
{inter['high_severity_count']} high-severity &middot;
{inter['moderate_severity_count']} moderate-severity findings
</div>
<div class="inter-list">{items}</div>
</section>
"""


def build_carrier_html(carr: dict | None) -> str:
    if not carr:
        return ""

    def render_group(items, title, css_class):
        if not items:
            return ""
        rows = ""
        for c in items:
            implication = _esc(
                c.get("carrier_implication")
                if items is carr["carriers"]
                else (c.get("affected_implication") if items is carr["affected"] else "—")
            )
            # Standard carrier-report context the producer already computes:
            # ethnicity-specific carrier frequency and any chip caveat.
            freq = c.get("carrier_frequency")
            if freq:
                implication += (
                    f'<div class="carr-freq"><em>Carrier frequency:</em> {_esc(freq)}</div>'
                )
            caveat = c.get("chip_caveat")
            if caveat:
                implication += f'<div class="carr-caveat">{_esc(caveat)}</div>'
            rsid = c.get("rsid", "")
            rsid_cell = (
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(rsid)}" '
                f'target="_blank" rel="noopener">{_esc(rsid)}</a>'
                if rsid
                else "—"
            )
            rows += f"""
<tr>
<td class="rsid-cell">{rsid_cell}</td>
<td><strong>{_esc(c["gene"])}</strong></td>
<td>{_esc(c["variant"])}</td>
<td>{_esc(c["disease"])}</td>
<td class="gt-cell">{_esc(c.get("genotype",""))}</td>
<td class="gt-cell">{_esc(c.get("pathogenic_allele",""))}</td>
<td>{_esc(c["inheritance"])}</td>
<td class="sum-cell">{implication}</td>
</tr>
"""
        return f"""
<h3 class="carr-h3 {css_class}">{_esc(title)} ({len(items)})</h3>
<div class="tbl-wrap">
<table class="snp-tbl">
<thead><tr><th>rsID</th><th>Gene</th><th>Variant</th><th>Condition</th><th>Genotype</th>
<th>Pathogenic allele</th><th>Inheritance</th><th>Implication</th></tr></thead>
<tbody>{rows}</tbody></table>
</div>
"""

    affected_html = render_group(carr["affected"], "AFFECTED (homozygous / fully penetrant)", "carr-aff")
    carrier_html = render_group(carr["carriers"], "Carriers (heterozygous)", "carr-car")
    not_carrier_html = ""
    if carr["not_carriers"]:
        names = ", ".join(f'{_esc(c["gene"])} {_esc(c["variant"])}' for c in carr["not_carriers"][:8])
        more = f" + {len(carr['not_carriers'])-8} more" if len(carr["not_carriers"]) > 8 else ""
        not_carrier_html = (
            f'<div class="carr-nc"><strong>Non-carriers for tested conditions:</strong> '
            f'{names}{more}</div>'
        )
    untested_html = ""
    if carr["untested"]:
        names = ", ".join(f'{_esc(c["gene"])} {_esc(c["variant"])}' for c in carr["untested"][:6])
        more = f" + {len(carr['untested'])-6} more" if len(carr['untested']) > 6 else ""
        untested_html = (
            f'<div class="carr-ut"><strong>Not on chip (cannot determine):</strong> '
            f'{names}{more}. Comprehensive carrier screening covers ~250+ '
            f'recessive conditions via clinical-grade sequencing.</div>'
        )

    return f"""
<section class="carr-section" id="carrier-status">
<h2>Carrier Status & Pathogenic Variants <span class="pro-pill">PRO</span></h2>
<p class="carr-intro">
Chip-detectable carrier status for conditions of family-planning and
personal-health relevance. <strong>This is NOT a substitute for clinical
carrier-screening panels</strong> (which cover ~250+ recessive conditions
using sequencing). It is the subset detectable from this chip's variants.
</p>
{affected_html}
{carrier_html}
{not_carrier_html}
{untested_html}
</section>
"""


def build_traits_html(tr: dict | None) -> str:
    if not tr or not tr.get("predictions"):
        return ""
    cards = ""
    for t in tr["predictions"]:
        conf = t.get("confidence", "n/a")
        conf_class = f"conf-{conf}"
        not_tested_class = " trait-na" if conf == "n/a" else ""
        detail_html = ""
        if t.get("detail"):
            detail_html = f'<div class="trait-detail">{_esc(t["detail"])}</div>'
        cards += f"""
<div class="trait-card{not_tested_class}">
  <div class="trait-name">{_esc(t["trait"])}</div>
  <div class="trait-result">{_esc(t["result"])}</div>
  <div class="trait-evidence">{_esc(t.get("evidence",""))}</div>
  {detail_html}
  <div class="trait-conf {conf_class}">Confidence: {_esc(conf)}</div>
</div>
"""
    return f"""
<section class="traits-section" id="trait-predictions">
<h2>Trait Predictions <span class="pro-pill">PRO</span></h2>
<p class="traits-intro">
Concrete genotype-based predictions for observable traits — not disease
risk. Confidence reflects both genotype clarity and the strength of the
genotype→phenotype link in the literature.
</p>
<div class="traits-grid">{cards}</div>
</section>
"""


def build_counseling_html(c: dict | None) -> str:
    if not c or not c.get("triggers"):
        return ""
    items = ""
    for t in c["triggers"]:
        items += f"""
<div class="couns-finding">
  <div class="couns-trigger">{_esc(t["trigger"])}</div>
  <div class="couns-urgency"><strong>Urgency:</strong> {_esc(t["urgency"])}</div>
  <div class="couns-specialist"><strong>Specialist:</strong> {_esc(t["specialist"])}</div>
  <div class="couns-reason">{_esc(t["reason"])}</div>
</div>
"""
    return f"""
<section class="couns-section" id="counseling-triggers">
<h2>Professional Consultation Triggers <span class="pro-pill">PRO</span></h2>
<p class="couns-intro">
Findings in this report that, per ACMG / NSGC / NCCN guidance, justify
discussion with a board-certified genetic counsellor or condition-specific
specialist. {c['n_triggers']} triggers identified
({c['urgent_count']} time-sensitive).
</p>
<div class="couns-list">{items}</div>
</section>
"""


def build_references_html(refs: list[dict] | None) -> str:
    if not refs:
        return ""
    rows = ""
    cur_cat = None
    for r in refs:
        if r["category"] != cur_cat:
            cur_cat = r["category"]
            rows += f'<tr class="ref-cat-row"><td colspan="5"><strong>{_esc(cur_cat)}</strong></td></tr>'
        pmid_links = " ".join(
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}" target="_blank" '
            f'rel="noopener">PMID:{p}</a>'
            for p in r.get("pmids", [])
        ) or "<span class='ref-no-pmid'>—</span>"
        guidelines = ""
        if r.get("guidelines"):
            guidelines = "<br>" + "; ".join(_esc(g) for g in r["guidelines"])
        # ClinVar clinical-significance assertion, curated per reference but
        # previously dropped. Shown as a labelled tag in the evidence cell.
        clinvar = ""
        if r.get("clinvar"):
            clinvar = f'<br><span class="ref-clinvar">ClinVar: {_esc(r["clinvar"])}</span>'
        rows += f"""
<tr>
<td class="rsid-cell"><a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(r["rsid"])}" target="_blank" rel="noopener">{_esc(r["rsid"])}</a></td>
<td><strong>{_esc(r["gene"])}</strong> {_esc(r["variant_name"])}</td>
<td><span class="ref-level {level_class(r["evidence_level"])}">{_esc(r["evidence_level"])}</span></td>
<td class="sum-cell">{_esc(r["evidence_summary"])}{clinvar}{guidelines}</td>
<td class="ref-pmids">{pmid_links}</td>
</tr>
"""
    return f"""
<section class="ref-section" id="references">
<h2>References & Evidence <span class="pro-pill">PRO</span></h2>
<p class="ref-intro">
Per-variant evidence levels and primary literature for the catalogued
findings in this report. Levels: <strong>CPIC A/B/C</strong> (pharmacogenomic
clinical practice guidelines), <strong>Clinical-Validated</strong> (single-
variant disease association at clinical-grade evidence), <strong>GWAS-A</strong>
(genome-wide significant, replicated in meta-analysis),
<strong>Population-Validated</strong> (population-genetic studies),
<strong>Catalogued</strong> (peer-reviewed GWAS or clinical literature; not
in this report's curated reference catalog).
</p>
<div class="tbl-wrap">
<table class="ref-tbl">
<thead><tr>
<th>rsID</th><th>Gene / Variant</th><th>Evidence Level</th>
<th>Summary</th><th>PMIDs</th>
</tr></thead>
<tbody>{rows}</tbody>
</table>
</div>
</section>
"""


# ── V3 section renderers ──────────────────────────────────────────────────────

def build_imputation_html(imp: dict | None) -> str:
    if not imp:
        return ""
    if not imp.get("available"):
        return (
            '<section class="impute-section" id="imputation">'
            '<h2>Imputation Status <span class="pro-pill">V3</span></h2>'
            f'<div class="impute-warn">Imputation not run — {_esc(imp.get("reason",""))}</div>'
            '</section>'
        )
    n_chip = imp.get("n_chip", 0)
    n_imp = imp.get("n_imputed", 0)
    pct_imp = 100 * n_imp / max(n_chip + n_imp, 1)
    return f"""
<section class="impute-section" id="imputation">
<h2>Statistical Imputation <span class="pro-pill">V3</span></h2>
<p class="impute-intro">
Genotypes were imputed using <strong>Beagle 5.4</strong> against the
<strong>1000 Genomes Phase 3</strong> reference panel. Imputation infers
likely genotypes at sites not directly typed by the chip, dramatically
expanding the variant set available for downstream polygenic scoring.
Each imputed variant carries a Beagle DR2 (R²) quality score; the report
filters at DR2 ≥ {_esc(imp.get('min_r2', 0.3))}.
</p>
<div class="impute-stats">
  <div class="impute-stat"><div class="impute-n">{n_chip:,}</div><div class="impute-l">Chip-typed variants</div></div>
  <div class="impute-stat impute-stat-imp"><div class="impute-n">{n_imp:,}</div><div class="impute-l">Imputed variants (DR2≥{imp.get('min_r2',0.3)})</div></div>
  <div class="impute-stat"><div class="impute-n">{pct_imp:.0f}%</div><div class="impute-l">of dataset from imputation</div></div>
  <div class="impute-stat"><div class="impute-n">{('cache' if imp.get('from_cache') else 'fresh')}</div><div class="impute-l">Source</div></div>
</div>
{f'<div class="impute-warn">Per-chromosome failures: {len(imp.get("failures",[]))}.</div>' if imp.get('failures') else ''}
</section>
"""


def build_expanded_pgs_html(epgs: dict | None) -> str:
    if not epgs:
        return ""
    if not epgs.get("available"):
        return (
            '<section class="epgs-section" id="expanded-pgs">'
            '<h2>Expanded Polygenic Risk Scores (PGS Catalog) <span class="pro-pill">V3</span></h2>'
            f'<div class="epgs-warn">{_esc(epgs.get("reason","Not available"))}</div>'
            '</section>'
        )
    cards = ""
    for slug, panel in epgs.get("panels", {}).items():
        result = panel.get("result", {})
        status = result.get("status")
        if status == "computed":
            tier = result["tier"]
            tier_cls = result["tier_class"]
            pct = result["percentile"]
            cov = result.get("coverage", {})
            confidence = result.get("confidence")
            conf_badge = _confidence_badge(confidence, result.get("confidence_note"))
            low_warn = ""
            if confidence == "low":
                low_warn = (
                    f'<div class="prs-lowconf">⚠️ Only '
                    f'{cov.get("pct_callable",0)}% of the scoring file is covered — '
                    'this percentile is unreliable and shown for transparency only.</div>'
                )
            cards += f"""
<div class="prs-card">
  <div class="prs-card-head">
    <div class="prs-name">{_esc(panel.get('label', slug))}
      <span class="prs-short">{_esc(panel.get('short',''))}</span></div>
    <div class="prs-tier {tier_cls}">{_esc(tier)}</div>
  </div>
  <div class="prs-track">
    <div class="prs-track-pointer" style="left:{pct}%"></div>
  </div>
  {conf_badge}
  {low_warn}
  <div class="prs-details">
    <div><span class="prs-lab">Percentile:</span> <strong>{pct:.0f}th</strong></div>
    <div><span class="prs-lab">Z-score:</span> <strong>{result.get('z_score', 0):+.2f}</strong></div>
    <div><span class="prs-lab">Coverage:</span>
      <strong>{cov.get('pct_callable',0)}%</strong>
      ({cov.get('chip',0)} chip + {cov.get('imputed',0)} imputed / {cov.get('total',0)})</div>
    {f'<div><span class="prs-lab">Low-r² imputed used:</span> {cov.get("low_r2",0):,}</div>' if cov.get('low_r2') else ''}
    {f'<div><span class="prs-lab">Not covered:</span> {cov.get("missing",0):,}</div>' if cov.get('missing') else ''}
  </div>
  <div class="prs-context">
    <strong>Lifetime context:</strong> {_esc(panel.get('population_lifetime_risk','—'))}
  </div>
  <div class="prs-ref"><em>{_esc(panel.get('pgs_id',''))}</em></div>
</div>
"""
        elif status in ("insufficient_data", "not_downloaded", "not_applicable", "error"):
            reason = result.get("reason", "")
            cards += (
                f'<div class="prs-card prs-card-na">'
                f'<div class="prs-name">{_esc(panel.get("label", slug))}</div>'
                f'{_confidence_badge(result.get("confidence", "none"))}'
                f'<div class="prs-na">{_esc(reason)}</div></div>'
            )

    # Section-level headline: elevated/high panels (computed by the module but
    # previously only surfaced by the curated-PRS renderer, not here).
    headline_html = ""
    hf = epgs.get("headline_findings") or []
    if hf:
        items = "".join(
            f'<li><strong>{_esc(p.get("label", slug))}</strong> — '
            f'<span class="prs-tier {p["result"]["tier_class"]}">'
            f'{_esc(p["result"]["tier"])}</span> '
            f'({p["result"]["percentile"]:.0f}th percentile)</li>'
            for slug, p in hf
            if p.get("result", {}).get("tier_class")
        )
        if items:
            headline_html = (
                f'<div class="prs-headline"><strong>Elevated / High PGS '
                f"Catalog findings:</strong><ul>{items}</ul></div>"
            )
    return f"""
<section class="epgs-section" id="expanded-pgs">
<h2>Expanded Polygenic Risk Scores <span class="pro-pill">V3</span> <span class="pro-pill">PGS Catalog</span></h2>
<p class="epgs-intro">
Polygenic risk scores using full <strong>PGS Catalog</strong> scoring files —
typically hundreds of thousands to millions of weighted SNPs per condition.
These are the same scores used in clinical-grade polygenic risk assessment.
Coverage depends on chip density and imputation: <em>chip + imputed / total</em>.
</p>
<div class="prs-caveat">
<strong>Informational use only — not diagnostic.</strong> A genotyping array
typically covers only a small fraction of a PGS Catalog scoring file, so any
percentile shown here is an approximation whose reliability is bounded by the
<em>coverage</em> stat on each card. Scores flagged <em>low confidence</em>
were computed on too few variants to interpret.
<strong>A percentile is a rank, not a probability.</strong> For low-base-rate
conditions — e.g. schizophrenia (~1% lifetime) — even a high percentile implies
only a small <em>absolute</em> risk: being in the top few percent of a ~1%
condition still leaves the absolute risk in the low single digits. Read every
percentile against the "Lifetime context" figure on its card. Polygenic risk is
one input alongside family history, lifestyle, and clinical risk factors;
confirm anything actionable with a clinician.
</div>
{headline_html}
<div class="prs-grid">{cards}</div>
</section>
"""


def _build_crosscheck_html(cc: dict | None) -> str:
    """Render the Y-DNA / mtDNA geographic cross-check box."""
    if not cc:
        return ""
    verdict = cc.get("verdict", "uninformative")
    palette = {
        "concordant":    ("#3fb950", "✓", "Lineages agree"),
        "plausible":     ("#d29922", "≈", "Broadly compatible"),
        "discordant":    ("#f85149", "⚠", "Lineage conflict"),
        "uninformative": ("#8b949e", "•", "No strong constraint"),
    }
    color, glyph, label = palette.get(verdict, palette["uninformative"])

    def _line_card(line, which):
        if not line:
            return ""
        dist = line.get("dist", {})
        chips = "".join(
            f'<span style="display:inline-block;padding:1px 7px;margin:2px;border-radius:10px;'
            f'background:var(--bg4);font-size:.8em">{_esc(sp)} {w*100:.0f}%</span>'
            for sp, w in sorted(dist.items(), key=lambda kv: -kv[1]) if w >= 0.05
        )
        conf = line.get("confidence", "")
        return (
            f'<div style="flex:1;min-width:240px;background:var(--bg2);border:1px solid var(--bdr);'
            f'border-radius:8px;padding:12px">'
            f'<div style="font-size:.85em;color:var(--muted);text-transform:uppercase;letter-spacing:.04em">{_esc(which)}</div>'
            f'<div style="font-size:1.35em;font-weight:700;margin:2px 0">{_esc(line.get("haplogroup","?"))}'
            f'<span style="font-size:.55em;color:var(--muted);font-weight:400"> · {_esc(conf)} conf</span></div>'
            f'<div style="color:var(--muted);font-size:.9em;margin-bottom:6px">{_esc(line.get("region",""))}</div>'
            f'<div>{chips}</div>'
            f'</div>'
        )

    cards = _line_card(cc.get("paternal"), "Paternal line · Y-DNA") + \
            _line_card(cc.get("maternal"), "Maternal line · mtDNA")
    explanations = "".join(
        f'<li>{ex}</li>' for ex in cc.get("explanations", [])
    )
    return f"""
<div class="anc-crosscheck" style="border:1.5px solid {color};border-radius:10px;
     padding:16px;margin:18px 0;background:linear-gradient(180deg,var(--bg3),var(--bg2))">
  <div style="display:flex;align-items:center;gap:10px;margin-bottom:4px">
    <span style="font-size:1.5em;color:{color}">{glyph}</span>
    <div>
      <div style="font-weight:700;font-size:1.1em">Deep-Lineage Cross-Check — {_esc(label)}</div>
      <div style="color:var(--muted);font-size:.9em">{_esc(cc.get("summary",""))}</div>
    </div>
  </div>
  <div style="display:flex;gap:12px;flex-wrap:wrap;margin:12px 0">{cards}</div>
  <ul style="margin:6px 0 0 18px;line-height:1.5">{explanations}</ul>
  <div style="color:var(--dim);font-size:.82em;margin-top:10px">
    Your Y-DNA (strict paternal line) and mtDNA (strict maternal line) each trace
    a <em>single</em> ancestor back tens of thousands of years — they aren't the
    same as a whole-genome admixture estimate, but they must be geographically
    <em>compatible</em> with it. When they aren't, the small autosomal panel is
    the suspect, not the haplogroup.
  </div>
</div>
"""


def _bar(value: float, max_value: float, color: str) -> str:
    """Simple inline SVG bar."""
    pct = max(0.0, min(100.0, 100.0 * value / max_value)) if max_value else 0
    return (
        f'<div style="height:10px;border-radius:5px;background:var(--bg3,#eef1f4);overflow:hidden">'
        f'<div style="height:100%;width:{pct:.1f}%;background:{color}"></div></div>')


def build_holistic_synthesis_html(hs: dict | None) -> str:
    """Cross-panel synthesis — genome leverage score + insights + priorities.
    Rendered FIRST after the exec summary so readers see the meta-view."""
    if not hs or not hs.get("available"):
        return ""

    gl = hs.get("genome_leverage") or {}
    score = gl.get("score", 50)
    tier = gl.get("tier", "—")
    tier_color = {"Very favorable": "#1a7f37", "Favorable": "#2a9d8f",
                  "Balanced": "#d29922", "Actionable risk": "#b3261e"}.get(tier, "#41505f")

    reasons_up = "".join(f'<div>+ {_esc(u)}</div>' for u in gl.get("reasons_up", []))
    reasons_down = "".join(f'<div>− {_esc(d)}</div>' for d in gl.get("reasons_down", []))
    reasons_col = ""
    if reasons_up:
        reasons_col += (f'<div style="color:#1a7f37;font-size:.86em;line-height:1.6">'
                        f'<div style="font-weight:700">Drivers up</div>{reasons_up}</div>')
    if reasons_down:
        reasons_col += (f'<div style="color:#b3261e;font-size:.86em;line-height:1.6;margin-top:8px">'
                        f'<div style="font-weight:700">Drivers down</div>{reasons_down}</div>')

    leverage_hero = f"""
<div style="display:flex;gap:20px;align-items:center;flex-wrap:wrap;
     background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1.5px solid {tier_color};
     border-radius:14px;padding:20px 24px;margin:14px 0">
  <div style="text-align:center;min-width:140px">
    <div style="font-size:.75em;color:#5b6673;text-transform:uppercase;letter-spacing:.05em">
      Genome Leverage Score</div>
    <div style="font-size:3em;font-weight:800;color:{tier_color};line-height:1">{score}</div>
    <div style="font-size:.85em;color:{tier_color};font-weight:700">{_esc(tier)}</div>
  </div>
  <div style="flex:1;min-width:280px">
    <div style="line-height:1.6;color:#33404d">{_esc(gl.get('narrative',''))}</div>
    <div style="display:flex;gap:20px;margin-top:12px;flex-wrap:wrap">{reasons_col}</div>
  </div>
</div>"""

    # Insights
    insights_html = ""
    for i in hs.get("insights", []):
        impact_color = {"actionable": "#b3261e", "caution": "#d29922",
                        "informational": "#2b5f8e"}.get(i.get("impact"), "#41505f")
        ev = "".join(
            f'<span style="background:#eef4fb;color:#12467a;border:1px solid #dbe3ec;'
            f'border-radius:14px;padding:2px 8px;margin:0 3px 0 0;font-size:.72em">'
            f'{_esc(e.get("module",""))}</span>'
            for e in i.get("evidence", [])
        )
        insights_html += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {impact_color};
     border-radius:8px;padding:12px 14px;margin:8px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;align-items:baseline;flex-wrap:wrap">
    <span style="font-weight:700;color:#12467a">{_esc(i['title'])}</span>
    <span style="font-size:.78em;color:{impact_color};font-weight:600">
      {_esc(i.get('impact','—'))} · sev {i.get('severity',1)}</span>
  </div>
  <div style="line-height:1.55;color:#33404d;margin:6px 0">{_esc(i['explanation'])}</div>
  <div style="line-height:1.55;color:#2b5f8e;font-size:.9em">
    <strong>Action:</strong> {_esc(i['action'])}</div>
  <div style="margin-top:6px">{ev}</div>
</div>"""

    # Priority actions
    priority_html = ""
    for a in hs.get("priority_actions", []):
        priority_html += f"""
<li style="margin:8px 0;line-height:1.55">
  <strong style="color:#12467a">{a['priority']}. {_esc(a['title'])}</strong>
  <div style="color:#4a5560;font-size:.9em;margin-top:3px">
    <em>Why:</em> {_esc(a['why'])}</div>
  <div style="color:#2b5f8e;font-size:.9em">
    <em>Action:</em> {_esc(a['action'])}</div></li>"""

    return f"""
<section class="holistic-section" id="holistic-synthesis">
<h2>Holistic Synthesis — Cross-Panel Patterns <span class="pro-pill">V6.13</span></h2>
<p class="anc-intro">
The synthesis layer: patterns that only exist when you look across panels.
No single module can catch these — they emerge from the interaction of
genotype, labs, ancestry, neurochemistry, and behaviour. A Genome Leverage
Score summarises the overall pattern; specific cross-panel insights follow,
each cited to the underlying findings that triggered it.
</p>
{leverage_hero}
<h3 style="margin:18px 0 6px">Cross-panel insights ({hs['n_insights']})</h3>
{insights_html or '<p style="color:var(--muted)">No cross-panel patterns detected.</p>'}
{f'<h3 style="margin:18px 0 6px">Prioritised actions</h3><ol style="line-height:1.55;padding-left:22px">{priority_html}</ol>' if priority_html else ''}
<div class="anc-caveat" style="margin-top:14px">
This is a heuristic synthesis layer; each individual insight is grounded in the
underlying modules but the composite is judgment, not clinical guidance.
</div>
</section>
"""


def _prepend_ai_interpretation(section_html: str, ai_text: str | None) -> str:
    """Inject a local-AI interpretation callout immediately under a section's
    <h2> heading. Used for 'AI on all tiers' — every module section gets its own
    AI interpretation when AI is enabled. No-op if there's no section or no text."""
    if not section_html or not ai_text:
        return section_html
    block = (
        '<div class="ai-interp" style="border:1px solid #cdd9e5;'
        'border-left:4px solid #6f42c1;background:linear-gradient(135deg,#f7f4fc,#f1f5fb);'
        'border-radius:8px;padding:12px 16px;margin:10px 0">'
        '<div style="font-weight:700;color:#6f42c1;font-size:.9em;margin-bottom:4px">'
        '🧠 AI interpretation</div>'
        f'<div style="line-height:1.6">{md_to_html(ai_text)}</div></div>')
    idx = section_html.find("</h2>")
    if idx == -1:
        return block + section_html
    return section_html[:idx + 5] + block + section_html[idx + 5:]


def build_life_stage_playbook_html(lsp: dict | None) -> str:
    """Decade-by-decade life-stage playbook, current decade highlighted."""
    if not lsp or not lsp.get("available"):
        return ""

    decades_html = ""
    for d in lsp.get("decades", []):
        cur = d.get("is_current")
        border = "#12467a" if cur else "#e3e7ec"
        bg = "linear-gradient(135deg,#eef4fb,#f4f8fc)" if cur else "#fff"
        here = ('<span style="background:#12467a;color:#fff;font-size:.7em;'
                'font-weight:700;padding:2px 8px;border-radius:10px;margin-left:8px">'
                'YOU ARE HERE</span>') if cur else ""
        base_items = "".join(f"<li style='margin:3px 0;line-height:1.5'>{_esc(x)}</li>"
                             for x in d.get("base", []))
        genome_items = ""
        for gi in d.get("genome_items", []):
            genome_items += (
                f'<li style="margin:4px 0;line-height:1.5">{_esc(gi["text"])}'
                f'<span style="color:#9aa4b0;font-size:.8em"> — {_esc(gi["source"])}</span></li>')
        genome_block = (f'<div style="font-weight:600;font-size:.85em;color:#12467a;'
                        f'margin-top:8px">Tuned to your genome:</div>'
                        f'<ul style="margin:4px 0 0 18px">{genome_items}</ul>'
                        ) if genome_items else ""
        decades_html += f"""
<div style="border:1.5px solid {border};border-radius:10px;padding:14px 18px;
     margin:10px 0;background:{bg};break-inside:avoid">
  <div style="font-size:1.1em;font-weight:800;color:#12467a">{_esc(d['label'])}{here}</div>
  <div style="color:#5b6673;font-style:italic;font-size:.9em;margin:2px 0 8px">{_esc(d['theme'])}</div>
  <ul style="margin:0 0 0 18px">{base_items}</ul>
  {genome_block}
</div>"""

    lev = ""
    if lsp.get("leverage_tier"):
        lev = (f'<div style="color:#5b6673;font-size:.9em;margin-top:4px">'
               f'Genome Leverage: <strong>{lsp["leverage_score"]}/100 '
               f'({_esc(lsp["leverage_tier"])})</strong> — the higher this is, the '
               f'more these decade choices actually move your outcome.</div>')

    return f"""
<section class="life-stage-section" id="life-stage-playbook">
<h2>Life-Stage Playbook <span class="pro-pill">V6.18</span></h2>
<p class="anc-intro">
Decade-by-decade priorities synthesised from your whole genome — an
evidence-based preventive-medicine baseline for each decade, modulated by your
specific genetic risk profile (every genome-driven item is tagged to its source).
{_esc(lsp.get('note',''))}
</p>
{lev}
{decades_html}
<div class="anc-caveat" style="margin-top:12px">
General preventive-medicine framework personalised by genotype — not a
substitute for a primary-care physician's individualised screening schedule.
Screening ages follow current US guidelines and may differ by family history
and evolving recommendations.
</div>
</section>
"""


def build_environmental_optimization_html(eo: dict | None) -> str:
    """Environmental optimization — light timing, exercise modality, vitamin-D
    seasonality."""
    if not eo or not eo.get("available"):
        return ""

    def _protocol_list(items):
        return "".join(f"<li style='margin:4px 0;line-height:1.5'>{_esc(x)}</li>"
                       for x in (items or []))

    # Circadian
    circ_html = ""
    c = eo.get("circadian")
    if c:
        circ_html = f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #7c5cbf;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff">
  <div style="font-weight:700">🌅 Circadian light timing — <span style="color:#7c5cbf">{_esc(c['lean'])}</span></div>
  <div style="color:#4a5560;line-height:1.55;font-size:.9em;margin:4px 0">{_esc(c['detail'])}</div>
  <ul style="margin:6px 0 0 18px">{_protocol_list(c['protocol'])}</ul>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:4px">{_esc(c['gene'])} {_esc(c['rsid'])} {_esc(c['genotype'])} · 📖 {_esc(c['citation'])}</div>
</div>"""

    # Exercise
    ex_html = ""
    e = eo.get("exercise")
    if e:
        basis = "".join(f'<li style="margin:2px 0">{_esc(b["call"])} '
                        f'<span style="color:#9aa4b0">({_esc(b["genotype"])})</span></li>'
                        for b in e.get("basis", []))
        ex_html = f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #2a9d8f;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff">
  <div style="font-weight:700">🏋️ Exercise-modality fit — <span style="color:#2a9d8f">{_esc(e['lean'])}</span></div>
  <div style="display:flex;gap:14px;margin:6px 0;flex-wrap:wrap">
    <div style="padding:5px 12px;background:#fff6e5;border-radius:6px">Power score <strong>{e['power_score']}</strong></div>
    <div style="padding:5px 12px;background:#eef4fb;border-radius:6px">Endurance score <strong>{e['endurance_score']}</strong></div>
  </div>
  <ul style="margin:4px 0 6px 18px">{basis}</ul>
  <div style="font-weight:600;font-size:.88em;color:#33404d">Training emphasis:</div>
  <ul style="margin:4px 0 0 18px">{_protocol_list(e['emphasis'])}</ul>
  <div style="color:#8a94a3;font-size:.84em;font-style:italic;margin-top:6px">{_esc(e['caveat'])}</div>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:3px">📖 {_esc(e['citation'])}</div>
</div>"""

    # Vitamin D
    vd_html = ""
    v = eo.get("vitamin_d")
    if v:
        vars_txt = ""
        if v.get("variants"):
            vars_txt = ("<div style='font-size:.8em;color:#8a94a3;margin-top:4px'>"
                        + " · ".join(f"{_esc(x['gene'])} {_esc(x['genotype'])}"
                                     for x in v["variants"]) + "</div>")
        vd_html = f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #d29922;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff">
  <div style="font-weight:700">☀️ Vitamin-D seasonality — <span style="color:#8a6d3b">{_esc(v['tendency'])}</span></div>
  <div style="color:#4a5560;font-size:.88em;margin:3px 0">Supplement months (from latitude {v['latitude_assumed']:.0f}°): <strong>{_esc(v['supplement_months'])}</strong></div>
  <ul style="margin:6px 0 0 18px">{_protocol_list(v['protocol'])}</ul>
  {vars_txt}
  <div style="color:#8a94a3;font-size:.84em;font-style:italic;margin-top:6px">{_esc(v['caveat'])}</div>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:3px">📖 {_esc(v['citation'])}</div>
</div>"""

    return f"""
<section class="env-optimization-section" id="environmental-optimization">
<h2>Environmental Optimization <span class="pro-pill">V6.17</span></h2>
<p class="anc-intro">
Actionable behavioural protocols from your genotype — the "what to actually do
differently" layer. Light timing tuned to your chronotype, training emphasis
tuned to your muscle/power genetics, and vitamin-D supplementation tuned to
your pathway variants × latitude (assumed {eo['latitude_assumed']:.0f}°N —
pass your real latitude for precision).
</p>
{circ_html}
{ex_html}
{vd_html}
<div class="anc-caveat" style="margin-top:12px">{_esc(eo.get('note',''))}</div>
</section>
"""


def build_polygenic_traits_html(pt: dict | None) -> str:
    """Trait genetics — genotype-level single-variant calls + explicit no-score
    handling for polygenic/fraught traits."""
    if not pt or not pt.get("available"):
        return ""

    conf_color = {"high": "#1a7f37", "moderate": "#2b5f8e", "low": "#8a94a3"}

    domain_html = ""
    for cat in pt.get("categories", []):
        cards = ""
        for f in pt["by_category"].get(cat, []):
            cc = conf_color.get(f["confidence"], "#8a94a3")
            cards += f"""
<div style="border:1px solid #e3e7ec;border-radius:8px;padding:11px 13px;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;align-items:baseline;flex-wrap:wrap">
    <span style="font-weight:700">{_esc(f['trait'])}</span>
    <span style="font-size:.75em;color:{cc}">{_esc(f['confidence'])}</span>
  </div>
  <div style="font-weight:600;color:#12467a;margin:3px 0">{_esc(f['call'])}</div>
  <div style="color:#4a5560;line-height:1.5;font-size:.88em">{_esc(f['detail'])}</div>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:3px">
    genotype {_esc(f['genotype'])} · {_esc(f['gene'])} {_esc(f['rsid'])} · 📖 {_esc(f['citation'])}</div>
</div>"""
        domain_html += (f'<h3 style="margin:16px 0 6px">{_esc(cat)}</h3>'
                        f'<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:10px">{cards}</div>')

    # No-score polygenic section (structural honesty)
    notes_html = ""
    for n in pt.get("polygenic_notes", []):
        notes_html += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #8a94a3;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fafbfc">
  <div style="font-weight:700;color:#41505f">{_esc(n['trait'])} — <span style="color:#b3261e">no score given</span></div>
  <div style="color:#4a5560;line-height:1.55;font-size:.9em;margin-top:4px">{_esc(n['why_no_number'])}</div>
  <div style="color:#2b5f8e;line-height:1.55;font-size:.9em;margin-top:4px"><strong>Honest take:</strong> {_esc(n['honest_statement'])}</div>
</div>"""

    return f"""
<section class="polygenic-traits-section" id="polygenic-traits">
<h2>Trait Genetics <span class="pro-pill">V6.16</span></h2>
<p class="anc-intro">
Genotype-level trait calls for traits where a <em>single variant</em> genuinely
carries signal — taste, smell, hair, eye colour, chronotype. {pt['n_findings']}
calls below. Deliberately <strong>no polygenic percentiles</strong>: no trait
PGS scoring files are available locally, and inventing a percentile from a few
SNPs would be dishonest.
</p>
{domain_html}
<h3 style="margin:18px 0 6px">Why some traits get no number</h3>
<p style="color:#5b6673;font-size:.9em">Height is too polygenic to call from a
chip; cognition and personality PGS explain only a few percent of variance and
are socially fraught. For these, an explanation replaces a fabricated score —
a hard design limit, not a caveat.</p>
{notes_html}
<div class="anc-caveat" style="margin-top:12px">{_esc(pt.get('methodology_note',''))}</div>
</section>
"""


def build_clinical_variants_html(cv: dict | None) -> str:
    """Phase-2 ClinVar rare/pathogenic screen. Deliberately conservative: strong
    framing, star-graded confidence, and a MANDATORY negative-result statement
    (shown whether findings are present or not) so an empty list is never read
    as 'clear'."""
    if not cv:
        return ""
    if not cv.get("available"):
        # Still render a short informational block so the section isn't silently
        # missing (e.g. ClinVar table not downloaded).
        reason = _esc(cv.get("reason", "not available"))
        return f"""
<section class="clinvar-section" id="clinical-variants">
<h2>Clinical Variants (ClinVar) <span class="pro-pill">V6.22</span> <span class="pro-pill">Phase 2</span></h2>
<div class="anc-caveat">Rare/pathogenic-variant screening not run: {reason}</div>
<div class="anc-caveat" style="margin-top:8px">{_esc(cv.get('negative_disclaimer',''))}</div>
</section>"""

    def _stars(n):
        return "★" * n + "☆" * (4 - n)

    cat_meta = [
        ("actionable", "🚨 Actionable (ACMG secondary findings)", "#b3261e"),
        ("affected", "Affected-consistent (recessive homozygous)", "#b3261e"),
        ("possible_compound_het", "Possible compound heterozygote (phase unknown)", "#d29922"),
        ("dominant_risk", "Dominant-gene risk variant", "#d29922"),
        ("xlinked", "X-linked variant", "#d29922"),
        ("carrier", "Carrier (recessive heterozygous)", "#2a9d8f"),
        ("uncertain_inheritance", "Pathogenic; inheritance not determined", "#5b6673"),
    ]
    buckets = cv.get("buckets", {})
    groups_html = ""
    for key, label, color in cat_meta:
        items = buckets.get(key) or []
        if not items:
            continue
        rows = ""
        for f in items:
            rows += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {color};border-radius:6px;
     padding:10px 13px;margin:7px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{_esc(f['gene'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(f['condition'])}</span></span>
    <span style="font-size:.8em;color:{color};font-weight:700" title="ClinVar review confidence">{_stars(f['stars'])}</span>
  </div>
  <div style="color:#4a5560;line-height:1.5;font-size:.9em;margin-top:3px">{_esc(f['interpretation'])}</div>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:3px">
    {_esc(f['chrom'])}:{f['pos']} {_esc(f['ref'])}&gt;{_esc(f['alt'])} · {_esc(f['significance'])} ·
    {_esc(f['zygosity'])} · inheritance {_esc(f.get('inheritance','unknown'))}</div>
</div>"""
        groups_html += f'<h3 style="margin:14px 0 4px;color:{color}">{label} ({len(items)})</h3>{rows}'

    vus_html = ""
    if cv.get("vus_in_acmg"):
        vus_rows = "".join(
            f"<li>{_esc(f['gene'])} — {_esc(f['condition'])} "
            f"<span style='color:#9aa4b0'>({_esc(f['chrom'])}:{f['pos']}, {_stars(f['stars'])})</span></li>"
            for f in cv["vus_in_acmg"])
        vus_html = f"""
<h3 style="margin:14px 0 4px;color:#5b6673">Variants of uncertain significance in actionable genes ({cv['n_vus_in_acmg']})</h3>
<p style="color:#8a94a3;font-size:.88em">Uncertain = <strong>not actionable</strong>. Listed only
because they fall in ACMG genes and may be reclassified in future ClinVar releases. Do not act on these.</p>
<ul style="margin:4px 0 0 18px;font-size:.9em">{vus_rows}</ul>"""

    none_found = cv["n_plp"] == 0
    headline = (
        '<div style="background:linear-gradient(135deg,#f2f9f4,#eef4fb);border:1px solid #cfe3d6;'
        'border-radius:10px;padding:14px 18px;margin:10px 0">'
        '<strong>No ClinVar-classified pathogenic variants matched.</strong> '
        'That is reassuring but limited — read the statement below.</div>'
        if none_found else
        f'<div style="background:linear-gradient(135deg,#fdf3f3,#f7f4fc);border:1px solid #f0cfcf;'
        f'border-radius:10px;padding:14px 18px;margin:10px 0">'
        f'<strong>{cv["n_plp"]}</strong> pathogenic / likely-pathogenic variant(s) matched — '
        f'{cv["n_actionable"]} actionable, {cv["n_carrier"]} carrier, {cv["n_affected"]} '
        f'affected-consistent. <strong>These are screening signals, not diagnoses.</strong></div>')

    return f"""
<section class="clinvar-section" id="clinical-variants">
<h2>Clinical Variants (ClinVar) <span class="pro-pill">V6.22</span> <span class="pro-pill">Phase 2</span></h2>
<p class="anc-intro">
Whole-genome/exome screen against <strong>ClinVar</strong> for known
pathogenic / likely-pathogenic variants — the sequencing-grade layer a chip
can't reach. Confidence is graded by ClinVar review stars (★); only ≥1★ P/LP
findings are shown. {cv['n_scanned']:,} VCF variants scanned.
</p>
{headline}
<div style="background:#fff6e5;border-left:4px solid #d29922;border-radius:0 8px 8px 0;
     padding:12px 16px;margin:10px 0;line-height:1.55;font-size:.92em">
  ⚠️ <strong>{_esc(cv.get('negative_disclaimer',''))}</strong></div>
{groups_html}
{vus_html}
<div class="anc-caveat" style="margin-top:14px">{_esc(cv.get('disclaimer',''))}</div>
</section>
"""


def build_novel_variants_html(nv: dict | None) -> str:
    """Phase-3 computational novel/rare-variant screen. Every finding is a MODEL
    PREDICTION, not a clinical call — the framing, per-finding label, and mandatory
    disclaimers all say so. Never presented with ClinVar-style certainty."""
    if not nv:
        return ""
    if not nv.get("available"):
        reason = _esc(nv.get("reason", "not available"))
        return f"""
<section class="clinvar-section" id="novel-variants">
<h2>Novel &amp; Rare Variants (predicted) <span class="pro-pill">Phase 3</span></h2>
<div class="anc-caveat">Computational predictor screen not run: {reason}</div>
<div class="anc-caveat" style="margin-top:8px">{_esc(nv.get('negative_disclaimer',''))}</div>
</section>"""

    conf_color = {"higher": "#b3261e", "moderate": "#d29922", "low": "#8a94a3"}
    cat_meta = [
        ("predicted_splice_disrupting", "🧬 Predicted splice-disrupting (SpliceAI)", "#b3261e"),
        ("predicted_pathogenic_rare", "🚨 Predicted pathogenic · rare / absent in gnomAD", "#b3261e"),
        ("predicted_pathogenic_uncommon", "Predicted pathogenic · uncommon", "#d29922"),
        ("predicted_pathogenic_common", "Predicted pathogenic · common (down-weighted)", "#5b6673"),
        ("ambiguous", "Ambiguous (predictor uncertain)", "#5b6673"),
    ]
    buckets = nv.get("buckets", {})
    groups_html = ""
    for key, label, color in cat_meta:
        items = buckets.get(key) or []
        if not items:
            continue
        rows = ""
        for f in items:
            gene = _esc(f.get("gene") or "—")
            cc = conf_color.get(f.get("confidence"), "#8a94a3")
            rows += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {color};border-radius:6px;
     padding:10px 13px;margin:7px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{gene}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(f.get('rarity','unknown'))}</span></span>
    <span style="font-size:.78em;color:{cc};font-weight:700" title="predictor consensus confidence">
      {_esc(f.get('confidence','').upper())} confidence</span>
  </div>
  <div style="color:#4a5560;line-height:1.5;font-size:.9em;margin-top:3px">{_esc(f.get('interpretation',''))}</div>
  <div style="font-size:.74em;color:#9aa4b0;margin-top:3px">
    {_esc(f['chrom'])}:{f['pos']} {_esc(f['ref'])}&gt;{_esc(f['alt'])} · {_esc(f.get('zygosity','?'))}<br>
    {_esc(f.get('evidence',''))}</div>
</div>"""
        groups_html += f'<h3 style="margin:14px 0 4px;color:{color}">{label} ({len(items)})</h3>{rows}'

    n_pred = nv.get("n_predicted_pathogenic", 0)
    none_found = n_pred == 0
    headline = (
        '<div style="background:linear-gradient(135deg,#f2f9f4,#eef4fb);border:1px solid #cfe3d6;'
        'border-radius:10px;padding:14px 18px;margin:10px 0">'
        '<strong>No carried variant was predicted damaging by the available models.</strong> '
        'Reassuring but limited — this is a computational screen; read the note below.</div>'
        if none_found else
        f'<div style="background:linear-gradient(135deg,#fdf3f3,#f7f4fc);border:1px solid #f0cfcf;'
        f'border-radius:10px;padding:14px 18px;margin:10px 0">'
        f'<strong>{n_pred}</strong> carried variant(s) predicted damaging '
        f'({nv.get("n_rare_damaging", 0)} rare/absent). <strong>These are model '
        f'predictions, not diagnoses</strong> — most predicted-damaging variants are '
        f'benign in reality.</div>')

    used = nv.get("predictors_used", [])
    pred_line = ", ".join(f"{p['name']} ({p['license']})" for p in used) or "none"
    safe_note = (" <span class='pro-pill'>commercial-safe</span>"
                 if nv.get("commercial_safe") else "")

    return f"""
<section class="clinvar-section" id="novel-variants">
<h2>Novel &amp; Rare Variants (predicted) <span class="pro-pill">Phase 3</span>{safe_note}</h2>
<p class="anc-intro">
Beyond ClinVar: every carried missense/splice SNV that is <em>not</em> a known
ClinVar hit is scored by offline computational predictors, and the
predicted-damaging, rare ones are surfaced. {nv.get('n_queried', 0):,} carried
variants queried. <strong>Predictors used:</strong> {_esc(pred_line)}.
</p>
{headline}
<div style="background:#fff6e5;border-left:4px solid #d29922;border-radius:0 8px 8px 0;
     padding:12px 16px;margin:10px 0;line-height:1.55;font-size:.92em">
  ⚠️ <strong>{_esc(nv.get('negative_disclaimer',''))}</strong></div>
{groups_html}
<div class="anc-caveat" style="margin-top:14px">{_esc(nv.get('disclaimer',''))}</div>
</section>
"""


def _render_wgs_decision(w: dict | None) -> str:
    """Whether sequencing is worth buying, and why the other figure says $0.

    The report shows a chip-to-sequencing marginal value that is structurally
    zero for anyone with array data. Left unexplained it reads as "sequencing
    would add nothing", which is not what it means and is not true.
    """
    if not w or not w.get("available"):
        return ""

    def money(v):
        try:
            v = round(float(v))
        except (TypeError, ValueError):
            return "&mdash;"
        return f"-${abs(v):,}" if v < 0 else f"${v:,}"

    accent = "#177a54" if w.get("worth_it") else "#b06a00"
    zero_note = ""
    if w.get("why_retrospective_is_zero"):
        zero_note = f"""
  <div style="border:1px solid #f0dcc0;background:#fffaf2;border-radius:8px;
              padding:10px 12px;margin:10px 0;font-size:.87em;color:#6b5330">
    <strong>Why the other sequencing figure on this page reads
    {money(w.get('retrospective_value'))}.</strong>
    {_esc(w['why_retrospective_is_zero'])}</div>"""

    return f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {accent};
            border-radius:10px;padding:14px 16px;margin:14px 0;background:#fbfcfe">
  <div style="display:flex;justify-content:space-between;align-items:baseline;
              gap:12px;flex-wrap:wrap">
    <div style="font-weight:700;color:{accent}">
      Is whole-genome sequencing worth buying?</div>
    <div style="font-size:.72em;color:#8a94a3;border:1px solid #dfe4ea;
                border-radius:20px;padding:2px 9px;white-space:nowrap">
      prospective &middot; before testing</div>
  </div>
  <div style="font-size:.95em;color:#2b3440;margin-top:8px">
    {_esc(w.get('plain',''))}</div>
  {zero_note}
  <table style="width:100%;border-collapse:collapse;font-size:.87em;margin-top:8px">
    <tr style="border-bottom:1px solid #f0f2f5">
      <td style="padding:4px">Chance of a serious actionable finding your array missed</td>
      <td style="text-align:right;padding:4px"><strong>1 in
        {w.get('number_needed_to_sequence') or '&mdash;'}</strong></td></tr>
    <tr style="border-bottom:1px solid #f0f2f5">
      <td style="padding:4px">Value if you are that person</td>
      <td style="text-align:right;padding:4px">{money(w.get('value_per_finding'))}</td></tr>
    <tr style="border-bottom:1px solid #f0f2f5">
      <td style="padding:4px">Expected value, averaged over everyone</td>
      <td style="text-align:right;padding:4px">{money(w.get('gross_expected_value'))}</td></tr>
    <tr style="border-bottom:1px solid #f0f2f5">
      <td style="padding:4px">&nbsp;&nbsp;of which pharmacogenomics
        <span style="color:#8a94a3">(mostly already on your array)</span></td>
      <td style="text-align:right;padding:4px;color:#8a94a3">
        {money(w.get('pgx_incremental_value'))}</td></tr>
    <tr style="border-bottom:1px solid #f0f2f5">
      <td style="padding:4px">Extra cost over an array</td>
      <td style="text-align:right;padding:4px">{money(w.get('incremental_cost'))}</td></tr>
    <tr><td style="padding:4px"><strong>Net expected value</strong></td>
      <td style="text-align:right;padding:4px;color:{accent}">
        <strong>{money(w.get('net_expected_value'))}</strong></td></tr>
  </table>
  <div style="font-size:.85em;color:#48545f;margin-top:9px">
    <strong>Where the value comes from.</strong> {_esc(w.get('biggest_driver',''))}</div>
  <div style="font-size:.84em;color:#8a94a3;margin-top:7px;font-style:italic">
    {_esc(w.get('caveat',''))}</div>
</div>"""


def _render_plain_summary(pl: dict | None) -> str:
    """The plain-English answer, placed before any of the technical output.

    Everything else in this section is correct and unreadable without training.
    This card carries the same conclusions in the order a person actually asks
    them — is it worth it, what do I get, what does it cost, how sure are we,
    what should I do, and what would change the answer — with the hedging
    attached to each sentence rather than exiled to a footnote.
    """
    if not pl or not pl.get("available"):
        return ""

    v = pl.get("verdict") or {}
    accent = {"positive": "#177a54", "negative": "#b03a2e"}.get(
        v.get("tone"), "#5b6673")
    conf = pl.get("confidence") or {}
    money = pl.get("money") or {}
    payback = pl.get("payback") or {}
    time_gain = pl.get("healthy_time") or {}

    actions = "".join(f"""
  <li style="margin-bottom:7px">
    <strong>{_esc(a.get('what',''))}</strong>
    <span style="color:#8a94a3">&mdash; {_esc(a.get('value',''))} of modelled value</span>
    <div style="color:#5b6673">{_esc(a.get('why',''))}</div>
    <div style="color:#8a94a3;font-size:.94em">{_esc(a.get('scale',''))}</div>
  </li>""" for a in (pl.get("actions") or []))

    nns = "".join(f"""
  <li style="margin-bottom:4px;color:#5b6673">{_esc(r.get('plain',''))}</li>"""
        for r in (pl.get("number_needed_to_screen") or [])[:5])

    return f"""
<div style="border:2px solid {accent};border-radius:12px;padding:16px 18px;
            margin:16px 0;background:#fcfdfe">
  <div style="display:flex;justify-content:space-between;align-items:baseline;
              gap:12px;flex-wrap:wrap">
    <div style="font-weight:700;color:{accent};font-size:1.1em">
      {_esc(v.get('headline',''))}</div>
    <div style="font-size:.72em;color:#8a94a3;border:1px solid #dfe4ea;
                border-radius:20px;padding:2px 9px;white-space:nowrap">
      in plain English</div>
  </div>

  <div style="font-size:.95em;color:#2b3440;margin-top:8px">
    {_esc(v.get('plain',''))}</div>

  <div style="display:flex;gap:12px;flex-wrap:wrap;margin:14px 0">
    <div style="flex:1;min-width:190px;background:#f7f9fb;border:1px solid #e3e7ec;
                border-radius:10px;padding:11px 13px">
      <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                  letter-spacing:.04em">What you get</div>
      <div style="font-size:.95em;color:#2b3440;margin-top:3px">
        {_esc(time_gain.get('plain',''))}</div>
    </div>
    <div style="flex:1;min-width:190px;background:#f7f9fb;border:1px solid #e3e7ec;
                border-radius:10px;padding:11px 13px">
      <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                  letter-spacing:.04em">What it costs</div>
      <div style="font-size:.95em;color:#2b3440;margin-top:3px">
        {_esc(money.get('plain',''))}</div>
    </div>
    <div style="flex:1;min-width:190px;background:#f7f9fb;border:1px solid #e3e7ec;
                border-radius:10px;padding:11px 13px">
      <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                  letter-spacing:.04em">When it pays for itself</div>
      <div style="font-size:.95em;color:#2b3440;margin-top:3px">
        {_esc(payback.get('plain',''))}</div>
    </div>
  </div>

  {f'''<div style="margin-top:4px">
    <div style="font-weight:600;color:#5b6673">How sure is this?</div>
    <div style="font-size:.92em;color:#48545f;margin-top:3px">
      {_esc(conf.get('plain',''))} {_esc(conf.get('range_plain',''))}</div>
  </div>''' if conf.get("available") else ""}

  {f'''<div style="margin-top:12px">
    <div style="font-weight:600;color:#5b6673">What to actually do, most valuable first</div>
    <ul style="font-size:.92em;margin:6px 0 0 18px;padding:0">{actions}</ul>
  </div>''' if actions else ""}

  {f'''<div style="margin-top:12px">
    <div style="font-weight:600;color:#5b6673">How often it makes a difference</div>
    <ul style="font-size:.9em;margin:6px 0 0 18px;padding:0">{nns}</ul>
    <div style="font-size:.85em;color:#8a94a3;margin-top:5px">
      {_esc(time_gain.get('caveat') or '')}</div>
  </div>''' if nns else ""}

  <div style="border:1px solid #f0dcc0;background:#fffaf2;border-radius:8px;
              padding:10px 12px;margin-top:13px;font-size:.88em;color:#6b5330">
    <strong>What would change this answer</strong><br>
    {_esc(pl.get('what_would_change_it',''))}</div>

  <div style="font-size:.82em;color:#9aa4b0;margin-top:10px;font-style:italic">
    {_esc(pl.get('disclaimer',''))}</div>
</div>"""


def _render_decision_layer(d: dict | None) -> str:
    """Render value of information, breakeven, frontier, budget and equity.

    These answer the questions that follow a cost-effectiveness result rather
    than restating it: what is it worth to resolve the guesswork, how far can
    an assumption move before the recommendation changes, which option wins
    once dominance is applied properly, what it does to a budget, and whether
    the benefit lands evenly.
    """
    if not d or not d.get("available"):
        return ""

    def money(v):
        try:
            v = round(float(v))
        except (TypeError, ValueError):
            return "—"
        return f"-${abs(v):,}" if v < 0 else f"${v:,}"

    blocks = ""

    # ── Value of information ──
    v = d.get("evpi") or {}
    pop = d.get("population_evpi") or {}
    if v.get("available"):
        rows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8">
  <td style="padding:3px 4px">{_esc(r['parameter'])}
    <span style="color:{'#b06a00' if r.get('tier') == 'assumption' else '#9aa4b0'};
                 font-size:.85em"> [{_esc(r.get('tier',''))}]</span></td>
  <td style="text-align:right;padding:3px 4px">{money(r['evppi_per_person'])}</td>
</tr>""" for r in (d.get("evppi") or []))
        near_zero = (v.get("evpi_per_person") or 0) < 50
        blocks += f"""
<div style="margin-top:10px"><div style="font-weight:600;color:#5b6673">
  What is it worth to stop guessing?</div>
<div style="font-size:.85em;color:#5b6673;margin:5px 0">
  Perfect information about every parameter would be worth
  <strong>{money(v.get('evpi_per_person'))}</strong> per person
  ({money(pop.get('population_evpi'))} across
  {pop.get('population', 0):,} people over {pop.get('years', 0)} years). The
  current recommendation is wrong in
  <strong>{v.get('p_current_choice_wrong', 0):.1%}</strong> of simulations.
  {"A figure this small does not mean the model is precise — it means the "
   "uncertainty does not straddle the decision. You would act the same way "
   "either way, so resolving it changes nothing you would do."
   if near_zero else ""}</div>
{f'''<table style="width:100%;border-collapse:collapse;font-size:.85em">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th style="padding:4px">Value of resolving this parameter alone</th>
    <th style="text-align:right">Per person</th></tr></thead>
  <tbody>{rows}</tbody></table>''' if rows else ""}</div>"""

    # ── Breakeven ──
    be = [b for b in (d.get("breakeven") or []) if b.get("available")]
    if be:
        items = "".join(f"""
<li style="margin-bottom:4px"><strong>{_esc(b['parameter'])}</strong>
  <span style="color:{'#b06a00' if b.get('tier') == 'assumption' else '#9aa4b0'};
               font-size:.85em">[{_esc(b.get('tier',''))}]</span> &mdash;
  <span style="color:{'#b03a2e' if b.get('crosses_within_range') else '#177a54'}">
  {_esc(b.get('interpretation',''))}</span></li>""" for b in be)
        blocks += f"""
<div style="margin-top:12px"><div style="font-weight:600;color:#5b6673">
  How wrong could each assumption be before the answer changes?</div>
<ul style="font-size:.85em;color:#5b6673;margin:5px 0 0 18px">{items}</ul></div>"""

    # ── Efficiency frontier ──
    fr = d.get("frontier") or {}
    if fr.get("available"):
        srows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8;
           color:{'#9aa4b0' if s['status'] != 'on frontier' else '#2b3440'}">
  <td style="padding:3px 4px">{_esc(s['name'])}</td>
  <td style="text-align:right;padding:3px 4px">{money(s['cost'])}</td>
  <td style="text-align:right;padding:3px 4px">{s['qaly']:.3f}</td>
  <td style="text-align:right;padding:3px 4px">{
      money(s['icer']) + '/QALY' if s.get('icer') is not None else '&mdash;'}</td>
  <td style="padding:3px 4px;font-size:.9em">{_esc(s['status'])}</td>
</tr>""" for s in fr.get("strategies", []))
        blocks += f"""
<div style="margin-top:12px"><div style="font-weight:600;color:#5b6673">
  Which strategy, not just whether to act</div>
<table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:4px">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th style="padding:4px">Strategy</th><th style="text-align:right">Cost</th>
    <th style="text-align:right">QALYs</th><th style="text-align:right">ICER</th>
    <th>Status</th></tr></thead><tbody>{srows}</tbody></table>
<div style="font-size:.8em;color:#8a94a3;margin-top:4px">
  Recommended at this threshold: <strong>{_esc(fr.get('recommended',''))}</strong>.
  Costs are net of averted cost, on the same basis as the incremental-cost
  figure above. {_esc(fr.get('note',''))}
  {f"<div style='margin-top:4px;color:#b06a00'>{_esc(fr['wgs_not_estimable'])}</div>"
   if fr.get('wgs_not_estimable') else ""}</div></div>"""

    # ── Budget impact ──
    # Deliberately NOT rendered here. The payer budget-impact analysis with
    # its per-member-per-month figures is already reported in the HEOR
    # deliverables block below, computed by markov_model under ISPOR
    # conventions. Showing a second one here gave the reader two answers to
    # one affordability question.
    bi = d.get("budget_impact") or {}
    if bi.get("available"):
        blocks += f"""
<div style="margin-top:12px;font-size:.85em;color:#5b6673">
  <strong>Affordability</strong> is reported separately, under the payer
  budget-impact analysis further down &mdash; peak
  <strong>${bi.get('peak_pmpm', 0):.4f}</strong> per member per month in year
  {bi.get('peak_year', '&mdash;')}. Cost-effectiveness and affordability are
  different questions and the conventions differ, so they are not merged.</div>"""

    # ── Subgroups + equity ──
    sub = d.get("subgroups") or {}
    dist = d.get("distributional") or {}
    if sub.get("available"):
        best, worst = sub.get("best") or {}, sub.get("worst") or {}
        eq = ""
        if dist.get("available"):
            eq = (f" Equity: the Gini coefficient of remaining health moves "
                  f"from {dist.get('gini_before')} to {dist.get('gini_after')}, "
                  f"so the programme "
                  f"<strong>{'narrows' if dist.get('reduces_inequality') else 'slightly widens'}"
                  f"</strong> the gap.")
        blocks += f"""
<div style="margin-top:12px"><div style="font-weight:600;color:#5b6673">
  Does it work the same for everyone?</div>
<div style="font-size:.85em;color:#5b6673;margin:5px 0">
  Net benefit for {_esc(sub.get('condition',''))} ranges from
  <strong>{money(best.get('inmb'))}</strong> ({_esc(str(best.get('sex','')))},
  age {best.get('age','')}) down to <strong>{money(worst.get('inmb'))}</strong>
  ({_esc(str(worst.get('sex','')))}, age {worst.get('age','')}) &mdash; a spread of
  {money(sub.get('spread'))}. {_esc(sub.get('note',''))}{eq}</div></div>"""

    if not blocks:
        return ""
    return f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Decision analysis &mdash; what to resolve, what to choose, what it costs</summary>
{blocks}</details>"""


def _assumption_dominance_note(tornado: list[dict]) -> str:
    """Say plainly when judgement calls, not evidence, drive the conclusion.

    A reader can work this out from the tornado table, but only by noticing
    the tier tags and doing the arithmetic. If most of the swing comes from
    declared assumptions, that is the single most important thing to know
    about the result and it should not be left as an exercise.
    """
    if not tornado:
        return ""
    total = sum(r.get("swing", 0) for r in tornado) or 1
    assumed = sum(r.get("swing", 0) for r in tornado
                  if r.get("tier") == "assumption")
    share = 100.0 * assumed / total
    if share < 25.0:
        return ""
    top = [r["parameter"] for r in tornado[:3] if r.get("tier") == "assumption"]
    named = ", ".join(f"<code>{_esc(t)}</code>" for t in top)
    return f"""
<div style="border:1px solid #f0dcc0;background:#fffaf2;border-radius:8px;
            padding:10px 12px;margin-top:8px;font-size:.85em;color:#6b5330">
  <strong>Judgement, not evidence, is driving this result.</strong>
  {share:.0f}% of the total swing in net benefit comes from parameters with no
  published anchor{f" &mdash; chiefly {named}" if named else ""}. The
  conclusion is therefore about as reliable as those judgements are, and a
  reader who disagrees with them should expect a materially different number.
  They are stated in full under parameter provenance below.</div>"""


def _render_pooled_economics(p: dict | None) -> str:
    """Render the pooled cost-effectiveness result.

    Four things this has to communicate that the old output did not: that
    findings bearing on one condition were combined rather than summed, that
    costs and quality-adjusted life-years are different objects reported
    separately, that the perspective is a choice with itemised consequences,
    and that a stated share of the model rests on judgement rather than
    literature.
    """
    if not p or not p.get("available"):
        return ""

    def money(v):
        try:
            v = round(float(v))
        except (TypeError, ValueError):
            return "—"
        return f"-${abs(v):,}" if v < 0 else f"${v:,}"

    cea = p.get("cea") or {}
    dc = p.get("double_counting") or {}
    adh = p.get("adherence") or {}
    prov = p.get("provenance") or {}

    # ── Per-condition table, with the correction made visible ──
    rows = ""
    for c in sorted(p.get("conditions") or [],
                    key=lambda r: r.get("inmb", 0), reverse=True):
        # Three tiers, in the order the corrections are applied: what naive
        # addition would have claimed, what the pooled trial evidence supports,
        # and what this cohort would actually realise. Showing only the last
        # would hide two judgement calls behind one number.
        naive = c.get("naive_additive_rrr", 0)
        eff = c.get("pooled_efficacy_rrr", c.get("combined_rrr", 0))
        pooled = c.get("combined_rrr", 0)
        stacked = c.get("n_findings", 0) > 1
        steps = ([f"{naive:.0%}"] if stacked else []) + [f"{eff:.0%}"]
        corr = (f'<span style="color:#b06a00">{" &rarr; ".join(steps)} &rarr; '
                f'<strong>{pooled:.0%}</strong></span>'
                if stacked or eff > pooled else f"{pooled:.0%}")
        srcs = ", ".join(c.get("sources") or []) or "&mdash;"
        rows += f"""
<tr style="border-bottom:1px solid #f0f2f5">
  <td style="padding:5px 4px"><strong>{_esc(c.get('condition',''))}</strong>
    <div style="font-size:.78em;color:#8a94a3">{_esc(srcs)}</div></td>
  <td style="text-align:right;padding:5px 4px">{c.get('n_findings',0)}</td>
  <td style="text-align:right;padding:5px 4px">{c.get('baseline_risk',0):.1%}</td>
  <td style="text-align:right;padding:5px 4px">{corr}</td>
  <td style="text-align:right;padding:5px 4px">{c.get('adherence', 1):.0%}
    <div style="font-size:.72em;color:#8a94a3">{_esc(
      (c.get('adherence_archetype') or '').replace('adherence_', ''))}</div></td>
  <td style="text-align:right;padding:5px 4px">{money(c.get('cost_averted'))}</td>
  <td style="text-align:right;padding:5px 4px">{c.get('qaly_gained',0):.3f}</td>
  <td style="text-align:right;padding:5px 4px">{money(c.get('inmb'))}</td>
</tr>"""

    # ── The correction banner — only when it actually bit ──
    corr_html = ""
    if dc.get("inflation_removed", 0) > 0:
        corr_html = f"""
<div style="border:1px solid #f0dcc0;background:#fffaf2;border-radius:10px;
            padding:12px 14px;margin:12px 0">
  <div style="font-weight:700;color:#8a5a00">Double-counting correction applied</div>
  <div style="font-size:.87em;color:#6b5330;margin-top:5px">
    Adding every finding's risk reduction separately would have claimed
    <strong>{money(dc.get('naive_cost_averted'))}</strong> of avoided cost. Several
    findings speak to the same condition &mdash; a polygenic score, a biomarker and a
    causal estimate for one trait are re-measurements of one liability, not
    independent chances to prevent three separate events. Combining them on the risk
    scale and charging each condition's cost once gives
    <strong>{money(dc.get('pooled_cost_averted'))}</strong>, removing
    <strong>{money(dc.get('inflation_removed'))}</strong>
    ({dc.get('pct_removed', 0)}%).</div>
</div>"""

    # ── Efficacy vs. effectiveness ──
    if adh.get("qaly_lost_to_non_adherence", 0) > 0:
        # Only claim the cost per QALY moves when a ratio is actually being
        # reported. In the dominance quadrants the ICER is suppressed, and a
        # card explaining why a withheld number got worse is nonsense.
        _icer_clause = (
            " That is why the cost per QALY worsens even though the "
            "intervention's own value for money barely moves."
            if (p.get("cea") or {}).get("icer") is not None else
            " Because it is spread over fewer realised QALYs, it is the part "
            "of this analysis most sensitive to how many people follow through.")
        corr_html += f"""
<div style="border:1px solid #d8e2ee;background:#f7fafd;border-radius:10px;
            padding:12px 14px;margin:12px 0">
  <div style="font-weight:700;color:#2f5d8a">Adherence: trial efficacy vs. real-world
    effectiveness</div>
  <div style="font-size:.87em;color:#425468;margin-top:5px">
    Every effect size above is measured in people who followed a protocol. Roughly half
    of people stop long-term preventive medication, and sustained behaviour change fares
    worse than that. Charged against the benefit, that gap costs
    <strong>{adh.get('qaly_lost_to_non_adherence', 0):.3f} QALYs</strong> and
    <strong>{money(adh.get('value_lost_to_non_adherence'))}</strong> of avoided cost
    &mdash; <strong>{adh.get('pct_of_benefit_lost', 0)}%</strong> of the health benefit
    the trials promise. The ongoing intervention cost is discounted by the same factor,
    because someone who stops taking a statin stops paying for it; the one-off
    {money(adh.get('fixed_test_cost'))} test cost is not.{_icer_clause}</div>
  <div style="font-size:.78em;color:#7b8794;margin-top:6px">
    {_esc(adh.get('src', ''))}. Screening uptake and behavioural maintenance are
    declared assumptions, varied in the sensitivity analysis.</div>
</div>"""

    # ── Disaggregated headline: cost and health kept apart ──
    icer = cea.get("icer")
    icer_txt = (money(icer) + "/QALY" if icer is not None
                else _esc(cea.get("icer_note") or "not applicable"))
    cards = f"""
<div style="display:flex;gap:12px;flex-wrap:wrap;margin:12px 0">
  <div style="flex:1;min-width:160px;background:#f7f9fb;border:1px solid #e3e7ec;
              border-radius:10px;padding:11px 13px">
    <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                letter-spacing:.04em">Incremental cost</div>
    <div style="font-size:1.3em;font-weight:700;color:#2b3440">
      {money(cea.get('incremental_cost'))}</div>
    <div style="font-size:.76em;color:#8a94a3">healthcare sector, negative = saves money</div>
  </div>
  <div style="flex:1;min-width:160px;background:#f7f9fb;border:1px solid #e3e7ec;
              border-radius:10px;padding:11px 13px">
    <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                letter-spacing:.04em">Incremental QALYs</div>
    <div style="font-size:1.3em;font-weight:700;color:#2b3440">
      {cea.get('incremental_qaly', 0):.3f}</div>
    <div style="font-size:.76em;color:#8a94a3">health gain, not monetised here</div>
  </div>
  <div style="flex:1;min-width:160px;background:#f7f9fb;border:1px solid #e3e7ec;
              border-radius:10px;padding:11px 13px">
    <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                letter-spacing:.04em">ICER</div>
    <div style="font-size:1.3em;font-weight:700;color:#2b3440">{icer_txt}</div>
    <div style="font-size:.76em;color:#8a94a3">cost per QALY gained</div>
  </div>
  <div style="flex:1;min-width:160px;background:#f7f9fb;border:1px solid #e3e7ec;
              border-radius:10px;padding:11px 13px">
    <div style="font-size:.76em;color:#8a94a3;text-transform:uppercase;
                letter-spacing:.04em">Net monetary benefit</div>
    <div style="font-size:1.3em;font-weight:700;color:#177a54">
      {money(cea.get('inmb'))}</div>
    <div style="font-size:.76em;color:#8a94a3">at {money(cea.get('wtp'))}/QALY</div>
  </div>
</div>
<div style="font-size:.87em;color:#48545f;margin:-2px 0 8px">
  <strong>Verdict:</strong> {_esc(cea.get('verdict',''))}</div>"""

    # ── Structural (Markov) re-estimate ──
    struct = [s for s in (p.get("structural") or []) if s.get("available")]
    struct_html = ""
    if struct:
        srows = "".join(f"""
<tr style="border-bottom:1px solid #f0f2f5">
  <td style="padding:4px">{_esc(s.get('condition',''))}</td>
  <td style="text-align:right;padding:4px">{money(s.get('incremental_cost'))}</td>
  <td style="text-align:right;padding:4px">{s.get('incremental_qaly',0):.3f}</td>
  <td style="text-align:right;padding:4px">{s.get('life_years_gained',0):.3f}</td>
  <td style="text-align:right;padding:4px">{
      money(s['icer']) + '/QALY' if s.get('icer') is not None else '&mdash;'}</td>
</tr>""" for s in struct)
        struct_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Structural re-estimate &mdash; cohort state-transition model</summary>
<div style="font-size:.85em;color:#5b6673;margin:7px 0">
  The figures above discount the whole horizon at its midpoint, which cannot
  express the fact that a person has to survive long enough to collect a
  prevented event. These rerun the top conditions through a three-state
  Well/Diseased/Dead cohort model against US life-table mortality, with
  Simpson's 1/3 within-cycle correction. Where the two disagree, the
  structural figure is the more careful one.</div>
<table style="width:100%;border-collapse:collapse;font-size:.85em">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th style="padding:4px">Condition</th><th style="text-align:right">Δ cost</th>
    <th style="text-align:right">Δ QALY</th><th style="text-align:right">Δ life-years</th>
    <th style="text-align:right">ICER</th></tr></thead>
  <tbody>{srows}</tbody></table></details>"""

    # ── Uncertainty: PSA interval, CEAC, tornado ──
    psa = p.get("psa") or {}
    unc_html = ""
    if psa.get("available"):
        curve = p.get("ceac") or []
        curve_html = ""
        if curve:
            pts = "".join(
                f"<tr><td style='padding:3px 4px'>{money(c['wtp'])}/QALY</td>"
                f"<td style='text-align:right;padding:3px 4px'>"
                f"{c['p_cost_effective']:.0%}</td></tr>" for c in curve)
            curve_html = f"""
<div style="margin-top:8px"><div style="font-weight:600;color:#5b6673">
  Probability cost-effective, by threshold</div>
<table style="border-collapse:collapse;font-size:.85em;margin-top:4px">{pts}</table></div>"""

        trows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8">
  <td style="padding:3px 4px">{_esc(t['parameter'])}
    <span style="color:{'#b06a00' if t['tier'] == 'assumption' else '#9aa4b0'};
                 font-size:.85em"> [{_esc(t['tier'])}]</span></td>
  <td style="text-align:right;padding:3px 4px;color:#8a94a3">
    {t['low_value']:g} &ndash; {t['high_value']:g}</td>
  <td style="text-align:right;padding:3px 4px"><strong>{money(t['swing'])}</strong></td>
</tr>""" for t in (p.get("tornado") or []))

        unc_html = f"""
<details style="margin-top:10px" open><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Uncertainty &mdash; {psa.get('n_iterations', 0):,} simulations over
  {psa.get('n_parameters_varied', 0)} parameters</summary>
<div style="font-size:.85em;color:#5b6673;margin:7px 0">
  Each iteration draws every parameter with a documented spread from its own
  distribution and re-runs the model. Net monetary benefit
  <strong>{money(psa.get('mean_inmb'))}</strong>
  (95% interval {money(psa.get('inmb_ci_low'))} to
  {money(psa.get('inmb_ci_high'))}); cost-effective in
  <strong>{psa.get('p_cost_effective', 0):.0%}</strong> of simulations at
  {money(psa.get('wtp'))}/QALY, cost-<em>saving</em> in
  <strong>{psa.get('p_cost_saving', 0):.0%}</strong>.</div>
<div style="font-size:.8em;color:#8a94a3;font-style:italic">
  {_esc(psa.get('note', ''))}</div>
{curve_html}
<div style="margin-top:10px"><div style="font-weight:600;color:#5b6673">
  What drives the answer &mdash; swing in net benefit across each parameter's range</div>
<table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:4px">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th style="padding:4px">Parameter</th><th style="text-align:right">Range</th>
    <th style="text-align:right">Swing</th></tr></thead>
  <tbody>{trows}</tbody></table>
<div style="font-size:.8em;color:#8a94a3;margin-top:5px">
  A parameter high on this list is one the conclusion depends on. Anything
  tagged <span style="color:#b06a00">[assumption]</span> is a judgement call
  carrying real weight &mdash; worth disagreeing with first.</div>
{_assumption_dominance_note(p.get("tornado") or [])}</div></details>"""

    dec_html = _render_decision_layer(p.get("decision") or {})
    plain_html = _render_plain_summary(p.get("plain") or {})
    wgs_html = _render_wgs_decision(p.get("wgs_decision") or {})

    # ── Dual perspective ──
    dp = p.get("dual_perspective") or {}
    dp_html = ""
    if dp:
        adds = "".join(f"""
<tr><td style="padding:3px 4px">{_esc(a.get('item',''))}</td>
  <td style="text-align:right;padding:3px 4px">{money(a.get('value'))}</td>
  <td style="padding:3px 4px;font-size:.9em;color:#8a94a3">{_esc(a.get('basis',''))}</td></tr>"""
            for a in (dp.get("societal_additions") or []))
        dp_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Perspective &mdash; healthcare sector vs societal</summary>
<div style="font-size:.85em;color:#5b6673;margin:7px 0">
  The reference case counts healthcare-sector costs only
  (<strong>{money((dp.get('healthcare_sector') or {}).get('cost_averted'))}</strong>
  averted). The societal perspective adds
  <strong>{money(dp.get('delta'))}</strong> more, itemised below so that a reader
  who rejects any one of these valuations can subtract it rather than having to
  discard the whole figure.</div>
<table style="width:100%;border-collapse:collapse;font-size:.85em">{adds}</table></details>"""

    # ── Impact inventory ──
    inv = p.get("impact_inventory") or []
    inv_html = ""
    if inv:
        irows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8">
  <td style="padding:3px 4px">{_esc(r.get('sector',''))}</td>
  <td style="padding:3px 4px">{_esc(r.get('item',''))}</td>
  <td style="padding:3px 4px;color:{'#177a54' if r.get('healthcare')=='included' else '#9aa4b0'}">
    {_esc(r.get('healthcare',''))}</td>
  <td style="padding:3px 4px;color:{'#177a54' if r.get('societal')=='included' else '#9aa4b0'}">
    {_esc(r.get('societal',''))}</td>
  <td style="padding:3px 4px;font-size:.9em;color:#8a94a3">{_esc(r.get('note',''))}</td>
</tr>""" for r in inv)
        inv_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Impact inventory &mdash; what is counted and what is not</summary>
<table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:6px">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th style="padding:4px">Sector</th><th>Item</th><th>Healthcare</th>
    <th>Societal</th><th>Note</th></tr></thead>
  <tbody>{irows}</tbody></table></details>"""

    # ── Validation ──
    checks = p.get("validation") or []
    val_html = ""
    if checks:
        n_pass = sum(1 for c in checks if c.get("pass"))
        crows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8">
  <td style="padding:4px;color:{'#177a54' if c.get('pass') else '#b03a2e'};
             font-weight:700;white-space:nowrap">{'PASS' if c.get('pass') else 'FAIL'}</td>
  <td style="padding:4px">{_esc(c.get('check',''))}
    <div style="font-size:.85em;color:#8a94a3">{_esc(c.get('detail',''))}</div></td>
</tr>""" for c in checks)
        val_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Model validation &mdash; {n_pass}/{len(checks)} checks pass</summary>
<div style="font-size:.85em;color:#5b6673;margin:7px 0">
  Internal consistency checks plus one external cross-validation against
  published statin cost-effectiveness analyses. A model that has never been
  asked whether it obeys its own constraints is not validated by having a lot
  of methods in it.</div>
<table style="width:100%;border-collapse:collapse;font-size:.87em">{crows}</table></details>"""

    # ── Provenance + declared assumptions + references ──
    assums = "".join(f"""
<li style="margin-bottom:5px"><strong>{_esc(a.get('key',''))}</strong> =
  {_esc(str(a.get('value','')))} {_esc(a.get('units',''))}<br>
  <span style="color:#8a94a3">{_esc(a.get('note',''))}</span></li>"""
        for a in (p.get("declared_assumptions") or []))
    refs = "".join(f"""
<li style="margin-bottom:4px">{_esc(r.get('source',''))}
  <span style="color:#8a94a3">{_esc(r.get('citation',''))}</span>
  <span style="color:#b0b7c0;font-size:.9em"> &mdash; {_esc(r.get('params',''))}</span></li>"""
        for r in (p.get("references") or []))
    prov_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  Parameter provenance &mdash; {prov.get('n_parameters', 0)} registered parameters,
  {prov.get('pct_sourced', 0)}% cited</summary>
<div style="font-size:.85em;color:#5b6673;margin:7px 0">
  Of the {prov.get('n_parameters', 0)} <em>registered</em> parameters &mdash; method
  conventions, cost-of-illness anchors, effect sizes and utilities &mdash;
  <strong>{prov.get('n_published', 0)}</strong> are read directly from a cited
  source, <strong>{prov.get('n_derived', 0)}</strong> are derived from one by a
  stated arithmetic step, and <strong>{prov.get('n_assumption', 0)}</strong> are
  judgement calls with no published anchor. The last group is listed in full
  rather than blended in with the others &mdash; they are the parts of this model
  most worth arguing with.</div>
<div style="border:1px solid #dfe7f0;background:#f7fafd;border-radius:8px;
            padding:10px 12px;margin:8px 0;font-size:.85em;color:#3f5266">
  <strong>Whole-model provenance.</strong> The percentage above describes the
  registry only. Counting the
  <strong>{prov.get('n_unregistered', 0)}</strong> per-finding figures in the
  curated module tables as well &mdash; cost, outcome value, prevalence and QALY
  fields &mdash; gives {prov.get('n_total_known', 0)} parameters in total, of which
  <strong>{prov.get('model_pct_resolvable', 0)}%</strong> carry a citation you can
  follow directly (PMID or DOI) and
  <strong>{prov.get('model_pct_attributed_or_better', 0)}%</strong> carry at least a
  named literature attribution. Only
  <strong>{prov.get('model_pct_unsourced', 0)}%</strong> rest on judgement alone,
  and those are the declared assumptions listed below.
  <div style="margin-top:6px">The remaining gap is
  <strong>{prov.get('n_curated_attributed', 0)}</strong> figures whose source is
  named but whose identifier has not been verified, across
  {prov.get('n_unresolved_sources', 0)} distinct sources. They are real citations,
  not invented numbers &mdash; they just cannot be resolved automatically yet.</div></div>
<details style="margin-top:6px"><summary style="cursor:pointer;color:#8a94a3;
    font-size:.84em">Unverified sources, largest first &mdash; the work queue</summary>
<ul style="font-size:.82em;color:#6a7683;margin:5px 0 0 18px">{
  "".join(f"<li>{_esc(u.get('source',''))} "
          f"<span style='color:#b0b7c0'>&mdash; {u.get('n_params',0)} figures</span></li>"
          for u in (prov.get('unresolved_sources') or []))
}</ul></details>
<div style="font-weight:600;color:#5b6673;margin-top:8px">Declared assumptions</div>
<ul style="font-size:.85em;color:#5b6673;margin:5px 0 0 18px">{assums}</ul>
<div style="font-weight:600;color:#5b6673;margin-top:10px">References</div>
<ul style="font-size:.82em;color:#5b6673;margin:5px 0 0 18px">{refs}</ul></details>"""

    # ── CHEERS 2022 ──
    cheers = p.get("cheers") or []
    cheers_html = ""
    if cheers:
        crows = "".join(f"""
<tr style="border-bottom:1px solid #f5f6f8">
  <td style="padding:3px 4px;white-space:nowrap;color:#5b6673">{_esc(c.get('item',''))}</td>
  <td style="padding:3px 4px">{_esc(c.get('response',''))}</td></tr>"""
            for c in cheers)
        cheers_html = f"""
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">
  CHEERS 2022 reporting checklist</summary>
<table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:6px">
  {crows}</table></details>"""

    return f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #177a54;border-radius:10px;
            padding:14px 16px;margin:16px 0;background:#fbfcfe">
  <div style="display:flex;justify-content:space-between;align-items:baseline;
              gap:12px;flex-wrap:wrap">
    <div style="font-weight:700;color:#177a54">Pooled cost-effectiveness analysis</div>
    <div style="font-size:.72em;color:#8a94a3;border:1px solid #dfe4ea;border-radius:20px;
                padding:2px 9px;white-space:nowrap">costs and QALYs reported separately</div>
  </div>
  {plain_html}
  <div style="font-size:.87em;color:#48545f;margin-top:6px">
    Findings are grouped by the condition they bear on, combined on the risk
    scale rather than added, and each condition is charged its cost of illness
    once.</div>
  <div style="font-size:.83em;color:#6a7683;margin-top:6px;padding:8px 10px;
              background:#f4f7fa;border-radius:8px">
    <strong>Which number is the answer?</strong> This one. The report shows three
    totals because they answer three questions: this pooled analysis is the
    <em>reference case</em> &mdash; what acting on the findings costs and gains.
    The <em>value-of-information</em> figure above it prices something different:
    what it was worth knowing any of this, which is why it is smaller and why
    resolving one uncertain parameter can move it. The
    <em>personal economic-impact</em> page itemises the same findings
    line-by-line for readability; its total runs higher because per-finding
    figures come from the curated module tables, which are not yet on the
    provenance registry.</div>
  {corr_html}
  {cards}
  <table style="width:100%;border-collapse:collapse;font-size:.88em;margin-top:4px">
    <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
      <th style="padding:5px 4px">Condition</th>
      <th style="text-align:right">Findings</th>
      <th style="text-align:right">Baseline risk</th>
      <th style="text-align:right">Risk reduction</th>
      <th style="text-align:right">Adherence</th>
      <th style="text-align:right">Cost averted</th>
      <th style="text-align:right">QALYs</th>
      <th style="text-align:right">NMB</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  {wgs_html}{unc_html}{dec_html}{struct_html}{dp_html}{inv_html}{val_html}{prov_html}{cheers_html}
</div>"""


def _build_personalized_voi_html(p: dict | None) -> str:
    """Individually-relevant HEOR panels for this genome: a personal cost-
    effectiveness frontier + CEAC, the genome as an appreciating data asset, and the
    person's percentile in the population value distribution. Market-level analyses
    are deliberately excluded — they belong in the separate business report."""
    if not p or not p.get("available"):
        return ""

    def _dm(v):
        try:
            return f"${int(round(float(v))):,}"
        except Exception:
            return "—"

    fr = p.get("frontier") or {}
    ceac = p.get("ceac") or {}
    ltv = p.get("ltv") or {}
    frows = "".join(
        f"<tr><td>{_esc(o['name'])}</td>"
        f"<td style='text-align:right'>${o['cost']:.0f}</td>"
        f"<td style='text-align:right'>{o['qaly']:.3f}</td>"
        f"<td style='text-align:right'>{(_dm(o['icer_vs_previous']) + '/QALY') if o.get('icer_vs_previous') else '—'}</td>"
        f"<td>{_esc(o.get('status',''))}</td></tr>"
        for o in (fr.get("frontier") or []))
    ruled = ", ".join(fr.get("ruled_out") or []) or "none"
    winner = fr.get("recommended_strategy", "—")
    # Report the CONFIDENCE IN THE RECOMMENDED strategy (from the CEAC at $100k), so
    # the probabilistic panel is consistent with the point-estimate recommendation.
    _at100 = next((r for r in (ceac.get("ceac") or []) if r.get("wtp") == 100_000), None)
    prob = (_at100 or {}).get("p_optimal", {}).get(winner)
    pct = p.get("population_percentile")
    return f"""
<div style="margin-top:14px;border:1px solid #cfe0d6;background:#f4faf6;border-radius:10px;padding:14px 16px">
  <div style="font-weight:700;color:#22683f;font-size:1.02em">For your genome specifically</div>
  <div style="font-size:.86em;color:#48545f;margin:4px 0 10px">{_esc(p.get('plain_english',''))}</div>
  <div style="display:flex;gap:14px;flex-wrap:wrap">
    <div style="flex:2;min-width:300px">
      <div style="font-weight:600;color:#5b6673;font-size:.9em;margin-bottom:4px">
        Your efficient testing choice (cost-effectiveness frontier)</div>
      <table style="width:100%;border-collapse:collapse;font-size:.84em">
        <thead><tr style="text-align:left;color:#8a94a3">
          <th>Option</th><th style="text-align:right">Cost</th>
          <th style="text-align:right">QALYs</th><th style="text-align:right">ICER</th>
          <th>Status</th></tr></thead>
        <tbody>{frows}</tbody></table>
      <div style="font-size:.8em;color:#9aa4b0;margin-top:4px">
        Ruled out as inefficient for you: {_esc(ruled)}.</div>
    </div>
    <div style="flex:1;min-width:180px">
      <div style="background:#fff;border:1px solid #e3e7ec;border-radius:8px;padding:10px 12px;margin-bottom:8px">
        <div style="font-size:1.3em;font-weight:700;color:#22683f">{_esc(winner)}</div>
        <div style="font-size:.78em;color:#8a94a3">efficient choice{f" · optimal in {prob:.0%} of simulations" if prob is not None else ""}</div>
      </div>
      <div style="background:#fff;border:1px solid #e3e7ec;border-radius:8px;padding:10px 12px;margin-bottom:8px">
        <div style="font-size:1.3em;font-weight:700;color:#22683f">{pct}<span style="font-size:.5em">th</span></div>
        <div style="font-size:.78em;color:#8a94a3">percentile of population value</div>
      </div>
      <div style="background:#fff;border:1px solid #e3e7ec;border-radius:8px;padding:10px 12px">
        <div style="font-size:1.3em;font-weight:700;color:#22683f">+{(ltv.get('appreciation_premium') or 0):.0%}</div>
        <div style="font-size:.78em;color:#8a94a3">appreciation over {ltv.get('years', 10)}y (free re-analysis)</div>
      </div>
    </div>
  </div>
  <div style="font-size:.76em;color:#b0868a;margin-top:8px">{_esc(p.get('caveat',''))}</div>
</div>"""


def build_voi_html(voi: dict | None) -> str:
    """Value-of-Information health-economics section: the discipline-grade headline
    answering 'what is knowing your genome worth?'. Every number is illustrative +
    probabilistic; market PRICE is shown separately from health-economic VALUE."""
    if not voi:
        return ""
    if not voi.get("available"):
        return f"""
<section class="clinvar-section" id="value-of-information">
<h2>Value of Information (health economics) <span class="pro-pill">VOI</span></h2>
<div class="anc-caveat">{_esc(voi.get('reason', 'not available'))}</div>
</section>"""

    def _m(v):
        return f"${v:,.0f}"

    mean = voi.get("voi_expost_mean", voi.get("voi_expost_point", 0))
    lo, hi = voi.get("voi_ci_low", 0), voi.get("voi_ci_high", 0)
    marginal = voi.get("marginal_chip_to_wgs", 0)
    p_ce = voi.get("prob_cost_effective", 0)
    var_95, cvar_95 = voi.get("var_95"), voi.get("cvar_95")

    # NMB table (top 8 by NMB).
    rows = ""
    for r in (voi.get("nmb_rows") or [])[:8]:
        tag = ("<span style='color:#8a5cf6;font-size:.8em'> · WGS-only</span>"
               if r.get("wgs_only") else "")
        color = "#2a9d8f" if r["nmb"] >= 0 else "#b3261e"
        rows += (f"<tr><td>{_esc(r['label'])}{tag}</td>"
                 f"<td style='text-align:right'>{r['dqaly']:.2f}</td>"
                 f"<td style='text-align:right'>{_m(r['dcost_averted'])}</td>"
                 f"<td style='text-align:right;color:{color};font-weight:700'>{_m(r['nmb'])}</td>"
                 f"<td style='color:#8a94a3;font-size:.85em'>{_esc(r['confidence'])}</td></tr>")

    # CEAC inline SVG.
    ceac = voi.get("ceac") or []
    ceac_svg = ""
    if ceac:
        W, H, PL, PB = 520, 170, 44, 26
        pts = []
        for c in ceac:
            x = PL + (c["lam"] / 200_000) * (W - PL - 10)
            y = (H - PB) - c["prob"] * (H - PB - 10)
            pts.append(f"{x:.1f},{y:.1f}")
        base_x = PL + (voi.get("wtp_base", 100_000) / 200_000) * (W - PL - 10)
        ceac_svg = f"""
<svg viewBox="0 0 {W} {H}" style="width:100%;max-width:540px;margin:8px 0">
  <line x1="{PL}" y1="{H-PB}" x2="{W-6}" y2="{H-PB}" stroke="#c8d0d8"/>
  <line x1="{PL}" y1="10" x2="{PL}" y2="{H-PB}" stroke="#c8d0d8"/>
  <line x1="{base_x:.0f}" y1="10" x2="{base_x:.0f}" y2="{H-PB}" stroke="#d29922" stroke-dasharray="4 3"/>
  <polyline points="{' '.join(pts)}" fill="none" stroke="#2a6df4" stroke-width="2.5"/>
  <text x="{PL}" y="{H-8}" font-size="10" fill="#8a94a3">$0</text>
  <text x="{W-60}" y="{H-8}" font-size="10" fill="#8a94a3">$200k/QALY</text>
  <text x="{PL-6}" y="16" font-size="10" fill="#8a94a3" text-anchor="end">1.0</text>
  <text x="{base_x:.0f}" y="8" font-size="9" fill="#d29922" text-anchor="middle">λ=$100k</text>
</svg>
<div style="font-size:.8em;color:#8a94a3">Cost-effectiveness acceptability curve — P(net benefit &gt; 0) vs willingness-to-pay.</div>"""

    # ROH -> carrier-panel prior. Rendered as a recommendation card, deliberately
    # without a dollar figure and visually separated from the valued findings, so
    # a qualitative judgement cannot be misread as a modelled amount.
    cp = voi.get("carrier_panel_prior") or {}
    panel_html = ""
    if cp.get("available"):
        _accent = {"none": "#5b6673", "founder": "#2a6df4",
                   "recent": "#8a5cf6"}.get(cp.get("tier"), "#5b6673")
        panel_html = f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {_accent};border-radius:10px;
            padding:14px 16px;margin:14px 0;background:#fbfcfe">
  <div style="display:flex;justify-content:space-between;align-items:baseline;gap:12px;flex-wrap:wrap">
    <div style="font-weight:700;color:{_accent}">Carrier-panel prior (from runs of homozygosity)</div>
    <div style="font-size:.72em;color:#8a94a3;border:1px solid #dfe4ea;border-radius:20px;
                padding:2px 9px;white-space:nowrap">not monetised &middot; no dollar value</div>
  </div>
  <div style="font-size:1.05em;font-weight:600;margin:7px 0 4px">{_esc(cp["recommendation"])}</div>
  <div style="font-size:.87em;color:#48545f">{_esc(cp["rationale"])}</div>
  <div style="font-size:.8em;color:#6a7683;margin-top:8px">
    <strong>Decision this informs:</strong> {_esc(cp["decision"])}</div>
  <div style="font-size:.78em;color:#8a94a3;margin-top:6px">
    F<sub>ROH</sub> {cp["f_roh"]:.4f} &middot; long-ROH {cp["f_roh_long"]:.4f}
    ({cp["n_long_runs"]} long run{"" if cp["n_long_runs"] == 1 else "s"})</div>
  <div style="font-size:.75em;color:#9aa4b0;margin-top:8px;font-style:italic">
    {_esc(cp["why_not_monetised"])}</div>
  <div style="font-size:.72em;color:#9aa4b0;margin-top:5px">{_esc(cp["src"])}</div>
</div>"""

    pooled_html = _render_pooled_economics(voi.get("pooled_economics") or {})

    price = voi.get("price", {})
    methods = "".join(f"<li>{_esc(m)}</li>" for m in (voi.get("methods") or []))

    # ── Extended economic framings: health capital, real options, risk-adjusted ──
    hc = voi.get("health_capital") or {}
    ro = voi.get("real_option") or {}
    ra = voi.get("risk_adjusted") or {}
    extended = ""
    if hc.get("available") or ro.get("available") or ra.get("available"):
        blocks = ""
        if hc.get("available"):
            defer = hc.get("years_floor_deferred")
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Health capital (Grossman 1972)</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Health is a stock that depreciates faster with age — δ rises from
    <strong>{hc.get('delta_at_age')}</strong> to <strong>{hc.get('delta_at_end')}</strong>/yr over the horizon.
    Information doesn't add health directly; it makes health <em>investment</em> more efficient, so its
    value compounds over the years you have left.</div>
  <div style="font-size:.86em;color:#177a54;margin-top:6px">
    Modelled gain: <strong>{hc.get('pv_health_capital_gain')}</strong> discounted health-capital-years
    {f"· morbidity threshold deferred ≈ <strong>{defer} yr</strong>" if defer else ""}</div>
</div>"""
        if ro.get("available"):
            ov = ro.get("option_value_of_waiting", 0)
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Optimal timing (real options)</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Sequencing is irreversible, prices fall (~{ro.get('assumed_cost_decline',0):.0%}/yr) and interpretation
    improves (~{ro.get('assumed_knowledge_growth',0):.0%}/yr) — but the data is a permanent asset you can
    re-analyse free, so waiting only forfeits protection.</div>
  <div style="font-size:.86em;color:#177a54;margin-top:6px">
    Option value of waiting: <strong>{_m(ov)}</strong> → {_esc(ro.get('recommendation',''))}</div>
</div>"""
        if ra.get("available"):
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Risk-adjusted view</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Return multiple <strong>{ra.get('roi_multiple')}×</strong> the test cost ·
    reward-to-variability <strong>{ra.get('reward_to_variability')}</strong> ·
    certainty equivalent <strong>{_m(ra.get('certainty_equivalent', 0))}</strong>
    (γ = {ra.get('risk_aversion_gamma')}).</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(ra.get('note',''))}</div>
</div>"""
        ev = voi.get("evpi") or {}
        if ev.get("available"):
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Expected Value of Perfect Information</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    The ceiling on what <em>any</em> further research could be worth:
    <strong>{_m(ev.get('evpi', 0))}</strong>. Current information already captures
    <strong>{(ev.get('share_of_information_captured') or 0) * 100:.1f}%</strong> of the
    achievable benefit.</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(ev.get('interpretation',''))}</div>
</div>"""
        ut = voi.get("utility") or {}
        if ut.get("available") and ut.get("by_gamma"):
            gam_rows = " · ".join(
                f"γ={g['gamma']:.0f}: {_m(g['certainty_equivalent'])}"
                for g in ut["by_gamma"])
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Risk preferences (Arrow–Pratt)</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    A risk-averse person values variance reduction above the expected value. Certainty
    equivalents — {gam_rows}.</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">Health information is a
    variance-reducing asset; that is why insurance markets exist (Arrow 1963).</div>
</div>"""
        pc = voi.get("penetrance_correction") or {}
        if pc.get("available"):
            blocks += f"""
<div style="flex:1;min-width:250px;background:#fdf9f3;border:1px solid #f0e2cf;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Ascertainment-corrected penetrance</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Literature penetrance <strong>{pc.get('prior_literature_penetrance')}</strong> →
    population-corrected <strong>{pc.get('population_corrected')}</strong> →
    posterior <strong>{pc.get('posterior_penetrance')}</strong>
    (×{pc.get('shrinkage_factor')}).</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(pc.get('note',''))}</div>
</div>"""
        ep = voi.get("evppi") or {}
        if ep.get("available"):
            top_rows = " · ".join(f"{r['parameter'].replace('_',' ')}: {_m(r['evppi'])}"
                                  for r in (ep.get("by_parameter") or [])[:3])
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Which uncertainty actually matters (EVPPI)</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">{_esc(top_rows)}</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(ep.get('interpretation',''))}</div>
</div>"""
        bh = voi.get("behavioural") or {}
        if bh.get("available"):
            blocks += f"""
<div style="flex:1;min-width:250px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Why people still don't test (behavioural)</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Normative value (3% exponential): <strong>{_m(bh.get('pv_exponential',0))}</strong> ·
    with present bias (β={bh.get('present_bias_beta')}): <strong>{_m(bh.get('pv_hyperbolic',0))}</strong> ·
    adoption gap <strong>{_m(bh.get('adoption_gap',0))}</strong>.</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(bh.get('interpretation',''))}</div>
</div>"""
        lon = voi.get("longevity") or {}
        if lon.get("available") and lon.get("scenarios"):
            lrows = "".join(
                f"<tr><td>{_esc(s['scenario'])}</td>"
                f"<td style='text-align:right'>{s['life_expectancy'] or '—'}</td>"
                f"<td style='text-align:right'>{s['blended']:.1%}</td>"
                f"<td style='text-align:right'>{(s['relative_uplift'] or 0):+.1%}</td></tr>"
                for s in lon["scenarios"])
            blocks += f"""
<div style="flex:1 1 100%;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">If people live longer, this is worth more</div>
  <table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:6px">
    <thead><tr style="text-align:left;color:#8a94a3"><th>Scenario</th>
      <th style="text-align:right">Life expectancy</th>
      <th style="text-align:right">Realised risk</th>
      <th style="text-align:right">vs today</th></tr></thead>
    <tbody>{lrows}</tbody></table>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(lon.get('direction',''))}</div>
</div>"""
        mk_ = voi.get("markov") or {}
        if mk_.get("available"):
            sc, gg = mk_.get("standard_care", {}), mk_.get("genomic_guided", {})
            icer_txt = (f"{_m(mk_['icer'])}/QALY" if mk_.get("icer") is not None
                        else "n/a")
            val = (mk_.get("validation") or {}).get("all_passed")
            blocks += f"""
<div style="flex:1 1 100%;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Markov cohort model (state-transition CEA)</div>
  <table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:6px">
    <thead><tr style="text-align:left;color:#8a94a3"><th>Strategy</th>
      <th style="text-align:right">Cost</th><th style="text-align:right">QALYs</th>
      <th style="text-align:right">Life-years</th></tr></thead>
    <tbody>
      <tr><td>Standard care</td><td style="text-align:right">{_m(sc.get('total_cost',0))}</td>
        <td style="text-align:right">{sc.get('total_qaly',0):.3f}</td>
        <td style="text-align:right">{sc.get('total_life_years',0):.3f}</td></tr>
      <tr><td>Genotype-guided</td><td style="text-align:right">{_m(gg.get('total_cost',0))}</td>
        <td style="text-align:right">{gg.get('total_qaly',0):.3f}</td>
        <td style="text-align:right">{gg.get('total_life_years',0):.3f}</td></tr>
      <tr style="font-weight:700;border-top:1px solid #e3e7ec">
        <td>Incremental</td><td style="text-align:right">{_m(mk_.get('incremental_cost',0))}</td>
        <td style="text-align:right">{mk_.get('incremental_qaly',0):+.3f}</td>
        <td style="text-align:right">{mk_.get('incremental_life_years',0):+.3f}</td></tr>
    </tbody></table>
  <div style="font-size:.86em;color:#177a54;margin-top:6px">
    ICER <strong>{icer_txt}</strong> · NMB {_m(mk_.get('nmb_at_wtp',0))} ·
    <strong>{_esc(mk_.get('verdict',''))}</strong></div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:4px">
    Annual cycles · half-cycle correction · p = 1−e<sup>−rΔt</sup> · age-dependent competing
    mortality · both costs and QALYs discounted{' · structural validation passed' if val else ''}.</div>
</div>"""
        bi = voi.get("budget_impact") or {}
        if bi.get("available"):
            birows = "".join(
                f"<tr><td>Year {r['year']}</td>"
                f"<td style='text-align:right'>{r['uptake']:.0%}</td>"
                f"<td style='text-align:right'>{r['tested_cumulative']:,}</td>"
                f"<td style='text-align:right'>{_m(r['offsets'])}</td>"
                f"<td style='text-align:right'>{_m(r['net_budget_impact'])}</td>"
                f"<td style='text-align:right'>${r['pmpm']:.4f}</td></tr>"
                for r in bi.get("rows", []))
            blocks += f"""
<div style="flex:1 1 100%;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Budget impact — the payer's question</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Plan of {bi.get('plan_members',0):,} members, {bi.get('eligible_population',0):,} eligible.
    Undiscounted, uptake phased in — <strong>PMPM</strong> is the decision metric.</div>
  <table style="width:100%;border-collapse:collapse;font-size:.85em;margin-top:6px">
    <thead><tr style="text-align:left;color:#8a94a3"><th>Year</th>
      <th style="text-align:right">Uptake</th><th style="text-align:right">Tested</th>
      <th style="text-align:right">Offsets</th><th style="text-align:right">Net impact</th>
      <th style="text-align:right">PMPM</th></tr></thead>
    <tbody>{birows}</tbody></table>
  <div style="font-size:.86em;color:#177a54;margin-top:6px">
    Peak PMPM <strong>${bi.get('peak_pmpm',0):.4f}</strong> (year {bi.get('peak_year')}) ·
    5-year cumulative {_m(bi.get('cumulative_net',0))}
    {'· becomes cost-saving within the horizon' if bi.get('becomes_cost_saving') else ''}</div>
</div>"""
        wf = voi.get("welfare") or {}
        if wf.get("available"):
            verdict = ("local analysis preferred" if wf.get("local_preferred")
                       else "centralised would win at this breach risk")
            blocks += f"""
<div style="flex:1 1 100%;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
  <div style="font-weight:700;color:#5b6673">Welfare: local vs centralised analysis</div>
  <div style="font-size:.86em;color:#5b6673;margin-top:4px">
    Conceding centralised platforms a <strong>+{wf.get('capability_premium_assumed',0):.0%}</strong>
    capability advantage (gap {_m(wf.get('capability_gap',0))}), local is preferred only if the
    expected privacy cost exceeds it. Expected privacy cost
    <strong>{_m(wf.get('expected_privacy_cost',0))}</strong>
    ({wf.get('prob_exposure_over_horizon',0):.0%} chance of exposure over the horizon — a genome
    cannot be revoked) → <strong>S<sub>local</sub> {_m(wf.get('surplus_local',0))}</strong> vs
    S<sub>central</sub> {_m(wf.get('surplus_central',0))}: <strong>{_esc(verdict)}</strong>.</div>
  <div style="font-size:.86em;color:#177a54;margin-top:6px">
    Break-even annual breach probability:
    <strong>{(wf.get('breakeven_annual_breach_prob') or 0):.2%}</strong> — below this, centralised
    analysis wins; above it, local does.</div>
  <div style="font-size:.78em;color:#9aa4b0;margin-top:6px">{_esc(wf.get('interpretation',''))}</div>
  <div style="font-size:.78em;color:#b0868a;margin-top:4px">{_esc(wf.get('caveat',''))}</div>
</div>"""
        ie = voi.get("information_economics") or {}
        ie_block = ""
        if ie.get("available"):
            ie_block = f"""
<div style="background:#f7f9fb;border-left:4px solid #6b3fd1;border-radius:0 8px 8px 0;
     padding:12px 16px;margin-top:12px;font-size:.86em;color:#5b6673;line-height:1.5">
  <strong>Information economics.</strong> {_esc(ie.get('adverse_selection_note',''))}
  {_esc(ie.get('discrimination_note',''))}
  <em>{_esc(ie.get('privacy_as_economic_design',''))}</em>
</div>"""
        extended = f"""
<h3 style="margin:16px 0 6px;color:#5b6673">Extended economic framings</h3>
<div style="display:flex;gap:14px;flex-wrap:wrap">{blocks}</div>{ie_block}"""

    return f"""
<section class="clinvar-section" id="value-of-information">
<h2>Value of Information — Health Economics <span class="pro-pill">VOI</span></h2>
<p class="anc-intro">
The return-on-health question, in the author's discipline: <strong>what is knowing
your genome worth, in expectation?</strong> A decision-analytic model over every
actionable finding — discounted net monetary benefit, cost-of-illness averted, and
pharmacogenomic averted-adverse-reactions — with a Monte-Carlo sensitivity analysis.
</p>
<div style="display:flex;gap:14px;flex-wrap:wrap;margin:10px 0">
  <div style="flex:1;min-width:220px;background:linear-gradient(135deg,#eef7f1,#eef4fb);
       border:1px solid #cfe3d6;border-radius:12px;padding:16px">
    <div style="color:#5b6673;font-size:.85em">Expected value of your genome (net of test cost)</div>
    <div style="font-size:2em;font-weight:800;color:#177a54">{_m(mean)}</div>
    <div style="color:#8a94a3;font-size:.82em">95% CI {_m(lo)} – {_m(hi)} · P(cost-effective @ $100k/QALY) = {p_ce:.0%}</div>
  </div>
  <div style="flex:1;min-width:220px;background:linear-gradient(135deg,#f3effe,#eef4fb);
       border:1px solid #ddd0f5;border-radius:12px;padding:16px">
    <div style="color:#5b6673;font-size:.85em">Marginal value of upgrading chip → whole genome</div>
    <div style="font-size:2em;font-weight:800;color:#6b3fd1">{_m(marginal)}</div>
    <div style="color:#8a94a3;font-size:.82em">added expected value the WGS-only findings unlock</div>
  </div>
{f'''  <div style="flex:1;min-width:220px;background:linear-gradient(135deg,#fdf3f3,#fdf7ee);
       border:1px solid #f0d9cf;border-radius:12px;padding:16px">
    <div style="color:#5b6673;font-size:.85em">Downside risk — a bad 1-in-20 outcome</div>
    <div style="font-size:2em;font-weight:800;color:#b3541e">{_m(cvar_95)}</div>
    <div style="color:#8a94a3;font-size:.82em">CVaR&#8329;&#8325; (expected shortfall) · VaR&#8329;&#8325; = {_m(var_95)} &mdash; the
      average and threshold of the worst 5% of simulated outcomes</div>
  </div>''' if var_95 is not None else ''}
</div>
<table style="width:100%;border-collapse:collapse;font-size:.9em;margin-top:6px">
  <thead><tr style="text-align:left;color:#5b6673;border-bottom:1px solid #e3e7ec">
    <th>Finding</th><th style="text-align:right">ΔQALY</th><th style="text-align:right">Cost averted</th>
    <th style="text-align:right">Net monetary benefit</th><th>Conf.</th></tr></thead>
  <tbody>{rows}</tbody>
</table>
{ceac_svg}
{pooled_html}
{panel_html}
{_build_personalized_voi_html(voi.get("personalized"))}
{extended}
<div style="display:flex;gap:14px;flex-wrap:wrap;margin-top:10px">
  <div style="flex:1;min-width:240px;background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;padding:12px 14px">
    <div style="font-weight:700;color:#5b6673">Market <em>price</em> (what these tests cost to buy)</div>
    <div style="font-size:.86em;color:#8a94a3;margin:4px 0">
      à-la-carte ≈ <strong>{_m(price.get('a_la_carte_total',0))}</strong> vs one consolidated run ≈
      <strong>{_m(price.get('consolidated',0))}</strong>.</div>
    <div style="font-size:.8em;color:#b0868a">{_esc(price.get('note',''))}</div>
  </div>
</div>
<details style="margin-top:10px"><summary style="cursor:pointer;color:#5b6673;font-weight:600">Methods &amp; assumptions</summary>
<ul style="font-size:.85em;color:#5b6673;margin:6px 0 0 18px">{methods}</ul></details>
<div class="anc-caveat" style="margin-top:12px">{_esc(voi.get('disclaimer',''))}</div>
</section>
"""


def build_family_planning_html(fp: dict | None) -> str:
    """In-report reproductive-genetics section — carrier compound risk with a
    random partner, dominant transmission (kept separate from penetrance),
    sex-gated mtDNA, and hereditary-cancer partner-screening guidance."""
    if not fp or not fp.get("available"):
        return ""

    def _pct(x):
        if x is None:
            return "—"
        if x >= 0.1:
            return f"{x*100:.0f}%"
        if x >= 0.01:
            return f"{x*100:.1f}%"
        return f"{x*100:.2f}%"

    # ── Recessive compound-risk cards ────
    rec_html = ""
    for it in fp.get("recessive_items", []):
        cr = it.get("child_two_copy_risk")
        clin = it.get("child_clinical_risk")
        risk_line = ""
        if cr is not None:
            clin_txt = (f" → ~{_pct(clin[0])}–{_pct(clin[1])} after penetrance"
                        if clin else "")
            risk_line = (
                f'<div style="margin:6px 0;padding:8px 12px;background:#eef4fb;'
                f'border-radius:6px;font-size:.92em">'
                f'<strong>Child two-copy risk with a random {_esc(it["ancestry"])} '
                f'partner:</strong> ~{_pct(cr)}{clin_txt}'
                f'<div style="color:#8a94a3;font-size:.85em;margin-top:2px">'
                f'= partner carrier freq ({_pct(it.get("partner_carrier_freq"))}) '
                f'× {"25% (you carrier)" if it["status"]=="carrier" else "50% (you affected)"}</div>'
                f'</div>')
        badge = "semi-dominant" if it.get("semidominant") else "recessive"
        rec_html += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #2a9d8f;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{_esc(it['disease'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(it['gene'])} {_esc(it['variant'])}</span></span>
    <span style="font-size:.78em;color:#2a9d8f;font-weight:600">{badge} · you: {_esc(it['status'])}</span>
  </div>
  {risk_line}
  <div style="color:#4a5560;line-height:1.55;font-size:.9em;margin-top:4px">{_esc(it['note'])}</div>
  <div style="font-size:.85em;color:#2b5f8e;margin-top:4px"><strong>Partner testing:</strong> {_esc(it.get('partner_test') or '—')}</div>
</div>"""

    # ── Dominant transmission cards (transmission ≠ disease) ────
    dom_html = ""
    for it in fp.get("dominant_items", []):
        dom_html += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #d29922;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{_esc(it['disease'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(it['gene'])} {_esc(it['variant'])}</span></span>
    <span style="font-size:.78em;color:#d29922;font-weight:600">dominant · you: {_esc(it['status'])}</span>
  </div>
  <div style="display:flex;gap:16px;margin:6px 0;flex-wrap:wrap">
    <div style="padding:6px 12px;background:#fff6e5;border-radius:6px">
      <div style="font-size:.72em;color:#8a6d3b;text-transform:uppercase">Transmission</div>
      <div style="font-weight:700;color:#8a6d3b">50% per child</div></div>
    <div style="padding:6px 12px;background:#eef2f7;border-radius:6px">
      <div style="font-size:.72em;color:#5b6673;text-transform:uppercase">Penetrance</div>
      <div style="font-weight:700;color:#5b6673">{_esc(it['penetrance_text'])}</div></div>
  </div>
  <div style="color:#8a94a3;font-size:.85em;font-style:italic">Inheriting the variant is not the same as developing the condition — most low-penetrance carriers never do.</div>
  <div style="color:#4a5560;line-height:1.55;font-size:.9em;margin-top:4px">{_esc(it['note'])}</div>
</div>"""

    # ── mtDNA (sex-gated) ────
    mt = fp.get("mtdna") or {}
    if mt.get("pathogenic_variants"):
        mt_vars = "".join(
            f'<li>{_esc(v["label"])} ({_esc(v["rsid"])}, {_esc(v["genotype"])}) — {_esc(v["disease"])}</li>'
            for v in mt["pathogenic_variants"])
        mt_vars = f'<ul style="margin:6px 0 0 18px">{mt_vars}</ul>'
    else:
        mt_vars = ('<div style="color:#8a94a3;font-size:.88em;margin-top:4px">'
                   'No pathogenic mtDNA variants detected among the few testable '
                   'on this chip.</div>')
    mt_html = f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #7c5cbf;border-radius:8px;
     padding:12px 14px;margin:8px 0;background:#fff">
  <div style="font-weight:700">Mitochondrial (maternal-line) transmission
    {f'· haplogroup {_esc(mt.get("haplogroup"))}' if mt.get("haplogroup") else ''}</div>
  <div style="color:#33404d;line-height:1.55;margin:6px 0">{_esc(mt['transmission_note'])}</div>
  {mt_vars}
  <div style="color:#8a94a3;font-size:.82em;margin-top:6px">{_esc(mt['screening_note'])}</div>
</div>"""

    # ── Hereditary cancer partner-screening ────
    hc_rows = ""
    for c in fp.get("hereditary_cancer", []):
        rel_color = {"HIGH": "#b3261e", "low": "#8a94a3"}.get(c["partner_relevance"], "#8a94a3")
        hc_rows += f"""
<tr>
  <td><strong>{_esc(c['syndrome'])}</strong><div style="color:#8a94a3;font-size:.85em">{_esc(c['genes'])}</div></td>
  <td style="white-space:nowrap">{_esc(c['inheritance'])}</td>
  <td style="color:{rel_color};font-weight:700;white-space:nowrap">{_esc(c['partner_relevance'])}</td>
  <td style="font-size:.9em">{_esc(c['note'])}</td>
</tr>"""
    hc_html = f"""
<h3 style="margin:18px 0 6px">Hereditary-cancer syndromes &amp; partner screening</h3>
<p style="color:#5b6673;font-size:.9em">Most hereditary-cancer syndromes are
<em>dominant</em> — a carrier parent transmits to each child with 50% probability
regardless of partner, so partner screening isn't the lever (cascade testing of
relatives + enhanced screening is). The exception is <strong>recessive</strong>
MUTYH, where partner screening genuinely matters.</p>
<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>
  <th>Syndrome / genes</th><th>Inheritance</th><th>Partner screening</th><th>Notes</th>
</tr></thead><tbody>{hc_rows}</tbody></table></div>"""

    return f"""
<section class="family-planning-section" id="family-planning">
<h2>Family Planning <span class="pro-pill">V6.15</span></h2>
<p class="anc-intro">
Your genome read as a reproductive-planning document. Risk figures assume an
unrelated <strong>{_esc(fp.get('ancestry_assumption','European'))}</strong>
partner (Hardy-Weinberg) unless noted. Two disciplines throughout: (1)
<em>transmission probability is not disease probability</em> — kept as separate
numbers; (2) <em>mtDNA transmission is sex-gated</em> — only mothers pass it on.
</p>
<div style="background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1px solid #dbe3ec;
     border-radius:10px;padding:14px 18px;margin:10px 0;line-height:1.6">{_esc(fp.get('summary',''))}</div>
{f'<h3 style="margin:16px 0 6px">Recessive / semi-dominant — child risk depends on partner</h3>{rec_html}' if rec_html else ''}
{f'<h3 style="margin:16px 0 6px">Dominant susceptibility — 50% transmission (≠ disease)</h3>{dom_html}' if dom_html else ''}
<h3 style="margin:16px 0 6px">Mitochondrial line</h3>
{mt_html}
{hc_html}
<div class="anc-caveat" style="margin-top:14px">
Chip-based and far from a complete carrier screen. Before or during early
pregnancy, a board-certified genetic counsellor + a clinical carrier-screening
panel (250+ recessive conditions via sequencing) is the standard of care.
Directory: findageneticcounselor.com.
</div>
</section>
"""


def build_addiction_genetics_html(ag: dict | None) -> str:
    """Addiction genetics — alcohol / opioid / nicotine / cannabis susceptibility."""
    if not ag or not ag.get("available"):
        return ""
    c = ag["composite"]

    tier_color = {
        "Strongly protected": "#1a7f37", "Modestly protected": "#2a9d8f",
        "Baseline (typical European)": "#8b949e", "Baseline": "#8b949e",
        "Mixed": "#d29922",
        "Modestly susceptible": "#d29922", "Elevated susceptibility": "#b3261e",
        "Low overall susceptibility": "#1a7f37",
    }
    alc_color = tier_color.get(c["alcohol_tier"], "#8b949e")
    overall_color = tier_color.get(c["overall_tier"], "#8b949e")

    hero = f"""
<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:12px;
     margin:12px 0">
  <div style="background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1.5px solid {alc_color};
       border-radius:12px;padding:14px 18px">
    <div style="font-size:.75em;color:#5b6673;text-transform:uppercase;letter-spacing:.05em">Alcohol tier</div>
    <div style="font-size:1.3em;font-weight:800;color:{alc_color};margin:4px 0">{_esc(c['alcohol_tier'])}</div>
    <div style="color:#4a5560;font-size:.9em;line-height:1.5">{_esc(c['alcohol_narrative'])}</div>
  </div>
  <div style="background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1.5px solid {overall_color};
       border-radius:12px;padding:14px 18px">
    <div style="font-size:.75em;color:#5b6673;text-transform:uppercase;letter-spacing:.05em">Overall susceptibility</div>
    <div style="font-size:1.3em;font-weight:800;color:{overall_color};margin:4px 0">{_esc(c['overall_tier'])}</div>
    <div style="color:#4a5560;font-size:.9em;line-height:1.5">{_esc(c['overall_narrative'])}</div>
    <div style="color:#8a94a3;font-size:.78em;margin-top:6px">
      {c['n_protective']} protective · {c['n_susceptible']} susceptible · {c['n_conditional']} conditional</div>
  </div>
</div>"""

    flags_html = ""
    if c["clinical_flags"]:
        rows = "".join(
            f'<div style="border-left:4px solid #12467a;background:#eef4fb;'
            f'padding:10px 14px;border-radius:0 8px 8px 0;margin:6px 0">'
            f'<div style="font-weight:700;color:#12467a">{_esc(f["title"])}</div>'
            f'<div style="color:#4a5560;line-height:1.5;margin-top:3px">{_esc(f["text"])}</div>'
            f'</div>' for f in c["clinical_flags"])
        flags_html = (f'<h3 style="margin:18px 0 6px">💡 Clinically useful flags</h3>'
                      f'{rows}')

    impact_color = {"protective": "#3fb950", "susceptible": "#f85149",
                    "clinically-relevant": "#12467a", "informational": "#8b949e",
                    "conditional": "#d29922", "neutral": "#8b949e"}
    impact_emoji = {"protective": "✅", "susceptible": "⚠", "clinically-relevant": "🩹",
                    "informational": "ℹ", "conditional": "◐", "neutral": "·"}

    domain_html = ""
    for cat in ag.get("categories", []):
        rows = ""
        for f in ag["by_category"].get(cat, []):
            border = impact_color.get(f["impact"], "#8b949e")
            emoji = impact_emoji.get(f["impact"], "·")
            rows += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid {border};
     border-radius:6px;padding:11px 13px;margin:8px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{emoji} {_esc(f['name'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(f['gene'])}</span></span>
    <span style="font-size:.8em;color:{border};font-weight:600">{_esc(f['impact'])}</span>
  </div>
  <div style="font-weight:600;margin:4px 0;color:#33404d">{_esc(f['verdict'])}</div>
  <div style="color:#4a5560;line-height:1.55;font-size:.9em">{_esc(f['mechanism'])}</div>
  <div style="font-size:.86em;color:#2b5f8e;margin-top:4px"><strong>Action:</strong> {_esc(f['action'])}</div>
  <div style="font-size:.75em;color:#9aa4b0;margin-top:3px">
    genotype {_esc(f['genotype'])} · {_esc(f['rsid'])} · {_esc(f['confidence'])} confidence · 📖 {_esc(f['citation'])}</div>
</div>"""
        domain_html += f'<h3 style="margin:16px 0 4px">{_esc(cat)}</h3>{rows}'

    return f"""
<section class="addiction-section" id="addiction-genetics">
<h2>Addiction Genetics <span class="pro-pill">V6.14</span></h2>
<p class="anc-intro">
Alcohol, opioid, nicotine, and cannabis susceptibility from your genome, plus
stress × substance-use interaction loci. Composite tier + clinically-useful
flags (naltrexone response, opioid dosing, never-smoke warning, trauma-informed-
care alerts). Grounded in the primary literature (Edenberg 2004; Higuchi 1994;
Anton 2008; Thorgeirsson 2008; Binder 2008).
</p>
{hero}
{flags_html}
{domain_html}
<div class="anc-caveat" style="margin-top:14px">
Genetics contributes ~50% of the variance in substance-use-disorder risk in
population studies, but that's variance across the population — not
deterministic within an individual. Behaviour, environment, life circumstances,
and trauma history dominate individual outcomes. This section is educational.
</div>
</section>
"""


def build_neurochemistry_html(nc: dict | None) -> str:
    """Neurochemistry — COMT axis + composite phenotype recommendations."""
    if not nc or not nc.get("available"):
        return ""
    c = nc["composite"]

    # ── Hero: composite phenotype summary ────
    def _pill(label, value, color):
        return (f'<div style="text-align:center;padding:10px 14px;background:#fff;'
                f'border:1px solid #e3e7ec;border-radius:8px">'
                f'<div style="font-size:.72em;color:#8a94a3;text-transform:uppercase;'
                f'letter-spacing:.05em">{label}</div>'
                f'<div style="font-size:1.05em;font-weight:700;color:{color}">{value}</div></div>')
    comt_col = {"warrior": "#b3261e", "worrier": "#2b5f8e",
                "middle": "#1a7f37"}.get(c["comt_class"], "#33404d")
    maoa_col = {"MAOA-H": "#1a7f37", "MAOA-L": "#d29922",
                "heterozygous": "#33404d"}.get(c["maoa_class"], "#33404d")
    bdnf_col = "#1a7f37" if "Val/Val" in c["bdnf_class"] else "#d29922"
    hero_pills = (
        _pill("COMT axis", c["comt_class"].title(), comt_col) +
        _pill("MAOA", c["maoa_class"], maoa_col) +
        _pill("BDNF plasticity", c["bdnf_class"].split(" ")[0], bdnf_col))

    subs = "".join(f"<li style='margin:4px 0'>{_esc(s)}</li>" for s in c.get("substance_flags", []))
    subs_html = (f"<h4 style='margin:10px 0 4px'>Substance / addiction flags</h4>"
                 f"<ul style='margin:0 0 6px 18px'>{subs}</ul>") if subs else ""

    hero = f"""
<div style="background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1px solid #dbe3ec;
     border-radius:12px;padding:16px 20px;margin:12px 0">
  <div style="display:flex;gap:10px;flex-wrap:wrap;margin-bottom:10px">{hero_pills}</div>
  <div style="line-height:1.6;color:#33404d">
    <div style="margin-bottom:4px"><strong>Stress-response profile:</strong> {_esc(c['stress_response_profile'])}</div>
    <div style="margin-bottom:4px"><strong>Plasticity tier:</strong> {_esc(c['plasticity_tier'])}</div>
  </div>
</div>
<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:10px;margin:12px 0">
  <div style="border:1px solid #dbe3ec;border-radius:8px;padding:11px 14px;background:#fff">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">Stimulant response</div>
    <div style="line-height:1.55">{_esc(c['stimulant_response'])}</div></div>
  <div style="border:1px solid #dbe3ec;border-radius:8px;padding:11px 14px;background:#fff">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">SSRI response</div>
    <div style="line-height:1.55">{_esc(c['ssri_response'])}</div></div>
  <div style="border:1px solid #dbe3ec;border-radius:8px;padding:11px 14px;background:#fff">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">Caffeine protocol</div>
    <div style="line-height:1.55">{_esc(c['caffeine_protocol'])}</div></div>
  <div style="border:1px solid #dbe3ec;border-radius:8px;padding:11px 14px;background:#fff">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">Meditation-style fit</div>
    <div style="line-height:1.55">{_esc(c['meditation_fit'])}</div></div>
  <div style="border:1px solid #dbe3ec;border-radius:8px;padding:11px 14px;background:#fff;grid-column:1/-1">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">Career neurotype</div>
    <div style="line-height:1.55">{_esc(c['career_neurotype'])}</div>{subs_html}</div>
</div>"""

    # ── Domain findings ────
    domain_html = ""
    for cat in nc.get("categories", []):
        rows = ""
        for f in nc["by_category"].get(cat, []):
            rows += f"""
<div style="border:1px solid #e3e7ec;border-left:4px solid #12467a;border-radius:6px;
     padding:11px 13px;margin:8px 0;background:#fff">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{_esc(f['name'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(f['gene'])}</span></span>
    <span style="font-size:.8em;color:#8a94a3">{_esc(f['confidence'])} confidence</span>
  </div>
  <div style="font-weight:600;color:#33404d;margin:4px 0">{_esc(f['phenotype'])}</div>
  <div style="line-height:1.55;color:#4a5560;font-size:.9em">{_esc(f['mechanism'])}</div>
  <div style="font-size:.86em;color:#2b5f8e;margin-top:4px"><strong>Action:</strong> {_esc(f['action'])}</div>
  <div style="font-size:.75em;color:#9aa4b0;margin-top:3px">
    genotype {_esc(f['genotype'])} · {_esc(f['rsid'])} · 📖 {_esc(f['citation'])}</div>
</div>"""
        domain_html += f'<h3 style="margin:16px 0 4px">{_esc(cat)}</h3>{rows}'

    return f"""
<section class="neurochemistry-section" id="neurochemistry">
<h2>Neurochemistry <span class="pro-pill">V6.12</span></h2>
<p class="anc-intro">
Your dopamine, serotonin, plasticity, stress, and reward-signalling profile — packaged
into a composite phenotype with concrete stimulant, SSRI, caffeine, meditation, and
career recommendations. Grounded in the primary literature (COMT: Egan 2001 PNAS,
Diamond 2007; BDNF: Egan 2003 Cell, Chen 2006 Science; MAOA: Caspi 2002 Science).
</p>
{hero}
{domain_html}
<div class="anc-caveat" style="margin-top:14px">
Behavioural genetics of common variants operates at the level of <em>population
tendencies</em>, not individual determinism. Effect sizes for most of these loci
are modest (OR 1.1-1.5 per allele). Use as a decision-support layer for practices
you'd already consider — not a verdict on who you are.
</div>
</section>
"""


def build_immunogenetics_html(ig: dict | None) -> str:
    """Immunogenetics — viral/bacterial/parasitic resistance + Historical Selection."""
    if not ig or not ig.get("available"):
        return ""
    impact_color = {"protective": "#3fb950", "susceptible": "#f85149",
                    "intermediate": "#d29922", "neutral": "#8b949e",
                    "informational": "#8b949e"}
    impact_emoji = {"protective": "🛡", "susceptible": "⚠",
                    "intermediate": "◐", "neutral": "·", "informational": "ℹ"}

    # Headline resistances — hero cards
    headlines_html = ""
    if ig.get("headlines"):
        cards = ""
        for h in ig["headlines"]:
            cards += f"""
<div style="border:1.5px solid #3fb950;border-radius:10px;padding:12px 14px;
     background:linear-gradient(135deg,#f2f9f4,#eef4fb)">
  <div style="font-size:1.4em">🛡</div>
  <div style="font-weight:700;color:#12467a">{_esc(h['name'])}</div>
  <div style="color:#5b6673;font-size:.9em;margin:4px 0">{_esc(h['verdict'])}</div>
  <div style="color:#8a94a3;font-size:.78em">{_esc(h['gene'])} · {_esc(h['rsid'])} · genotype {_esc(h['genotype'])}</div>
</div>"""
        headlines_html = (f'<h3 style="margin:14px 0 8px;color:#1a7f37">🥇 Headline Resistances</h3>'
                          f'<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:10px;margin-bottom:14px">{cards}</div>')

    # All findings by category
    domain_html = ""
    for cat in ig.get("categories", []):
        items = ig["by_category"].get(cat, [])
        if not items:
            continue
        rows = ""
        for f in items:
            border = impact_color.get(f["impact"], "#8b949e")
            emoji = impact_emoji.get(f["impact"], "·")
            hist = (f'<div style="margin-top:6px;font-size:.82em;color:#12467a;'
                    f'background:#eef4fb;border-radius:4px;padding:4px 8px">'
                    f'📜 <strong>Historical:</strong> {_esc(f["historical"])}</div>') if f.get("historical") else ""
            rows += f"""
<div style="border:1px solid var(--bdr,#e3e7ec);border-left:4px solid {border};
     border-radius:6px;padding:11px 13px;margin:8px 0;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <span style="font-weight:700">{emoji} {_esc(f['name'])}
      <span style="font-weight:400;color:#8a94a3;font-size:.85em"> · {_esc(f['gene'])}</span></span>
    <span style="font-size:.8em;color:{border};font-weight:600">{_esc(f['impact'])}</span>
  </div>
  <div style="font-weight:600;margin:4px 0;color:#33404d">{_esc(f['verdict'])}</div>
  <div style="color:#4a5560;line-height:1.5;font-size:.9em">{_esc(f['mechanism'])}</div>
  <div style="font-size:.86em;color:#2b5f8e;margin-top:4px"><strong>Action:</strong> {_esc(f['action'])}</div>
  <div style="font-size:.75em;color:#9aa4b0;margin-top:3px">
    genotype {_esc(f['genotype'])} · {_esc(f['rsid'])} · {_esc(f['confidence'])} confidence · 📖 {_esc(f['citation'])}</div>
  {hist}
</div>"""
        domain_html += f'<h3 style="margin:16px 0 4px">{_esc(cat)}</h3>{rows}'

    # Historical Selection Timeline
    timeline_html = ""
    tl = ig.get("historical_timeline") or []
    if tl:
        rows = ""
        for ev in tl:
            rows += f"""
<tr>
  <td style="white-space:nowrap;color:#8a94a3">{_esc(ev['epoch'])}</td>
  <td><strong>{_esc(ev['driver'])}</strong></td>
  <td>{_esc(ev['finding'])}</td>
  <td style="color:#4a5560">{_esc(ev['verdict'])}</td>
</tr>"""
        timeline_html = f"""
<h3 style="margin:20px 0 6px;color:#8a4900">📜 Historical Selection Timeline</h3>
<p style="color:#5b6673;font-size:.9em">Each protective variant in your genome maps to a
historical pressure that selected for it. Reading down this list is reading a summary of
what your ancestors survived.</p>
<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>
  <th>Epoch</th><th>Selection pressure</th><th>Your variant</th><th>Effect</th></tr></thead>
<tbody>{rows}</tbody></table></div>"""

    return f"""
<section class="immunogenetics-section" id="immunogenetics">
<h2>Immunogenetics <span class="pro-pill">V6.11</span></h2>
<p class="anc-intro">
Your genetic resistance and susceptibility to human pathogens — viral, bacterial,
parasitic — plus a Historical Selection Timeline that connects your protective
variants to the pandemics that shaped them. {ig['n_findings']} findings:
<span style="color:#1a7f37"><strong>{ig['n_protective']} protective</strong></span> ·
<span style="color:#b3261e"><strong>{ig['n_susceptible']} susceptible</strong></span> ·
{ig['n_intermediate']} intermediate.
</p>
{headlines_html}
{domain_html}
{timeline_html}
<div class="anc-caveat" style="margin-top:16px">
Consumer-chip-based resistance / susceptibility estimates. Population-average
effects; individual outcomes depend on exposure, vaccination, HLA, and behaviour.
Vaccination and standard infection prevention still apply regardless of genotype.
</div>
</section>
"""


def build_ancestral_story_html(story: dict | None) -> str:
    """Long-form Ancestral Story — deterministic chapters + optional AI narrative."""
    if not story or not story.get("available"):
        return ""

    ai_html = ""
    if story.get("ai_used") and story.get("ai_text"):
        # Convert markdown paragraphs → HTML paragraphs
        paragraphs = story["ai_text"].split("\n\n")
        body = ""
        for p in paragraphs:
            p = p.strip()
            if not p:
                continue
            if p.startswith("### "):
                body += f'<h3 style="color:#12467a;margin:12px 0 4px">{_esc(p[4:].strip())}</h3>'
            elif p.startswith("## "):
                body += f'<h2 style="color:#12467a;margin:14px 0 6px;border:none">{_esc(p[3:].strip())}</h2>'
            elif p.startswith("**Chapter"):
                body += f'<h3 style="color:#12467a;margin:14px 0 4px">{_esc(p.strip("*"))}</h3>'
            else:
                # Preserve inline bold
                import re as _re
                p_html = _esc(p)
                p_html = _re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', p_html)
                body += f'<p style="line-height:1.65;margin:8px 0">{p_html}</p>'
        ai_html = f"""
<div style="border-left:4px solid #12467a;padding-left:16px;margin:14px 0;
     background:linear-gradient(180deg,#fbfcfd,#f5f8fb);padding:16px 20px;border-radius:0 8px 8px 0">
  <div style="color:#12467a;font-weight:700;margin-bottom:6px">✨ AI-Enhanced Narrative</div>
  {body}
</div>"""

    # Deterministic chapters
    chapters_html = ""
    for ch in (story.get("template") or {}).get("chapters", []):
        # Markdown-lite → HTML
        body = _esc(ch["body"])
        import re as _re
        body = _re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', body)
        body = _re.sub(r'### (.+)', r'<h4 style="color:#12467a;margin:10px 0 4px">\1</h4>', body)
        body = body.replace("\n\n", "</p><p>")
        chapters_html += f"""
<article style="border:1px solid #e3e7ec;border-radius:10px;padding:16px 20px;margin:12px 0;background:#fff">
  <h3 style="color:#12467a;margin:0 0 8px">{_esc(ch['title'])}</h3>
  <div style="line-height:1.65"><p>{body}</p></div>
</article>"""

    return f"""
<section class="ancestral-story-section" id="ancestral-story">
<h2>The Ancestral Story <span class="pro-pill">V6.11</span></h2>
<p class="anc-intro">
A long-form narrative weaving together your Y-DNA and mtDNA haplogroups, your
autosomal ancestry components (Neanderthal, Yamnaya, EEF, WHG), and the
historical-selection events visible in your immune genome. This is the story
of what your ancestors likely lived through, what they believed, what they ate
and drank, and what pandemics shaped them. Educational and historical — a
plausibility-weighted narrative from anonymous DNA-derived data, not a
documented lineage.
</p>
{ai_html}
{chapters_html}
</section>
"""


def build_blood_type_html(bt: dict | None) -> str:
    """Blood type inference — ABO + RhD + FUT2 secretor status."""
    if not bt or not bt.get("available"):
        return ""
    combined = bt.get("combined") or "—"
    a = bt["abo"]
    d = bt["rhd"]
    s = bt["secretor"]

    # Pick a colour by rarity
    rarity = {"O+": "#3fb950", "A+": "#3fb950", "B+": "#d29922", "AB+": "#d29922",
              "O-": "#f85149", "A-": "#f85149", "B-": "#b3261e", "AB-": "#b3261e"}
    hero_color = rarity.get(combined, "#8b949e")

    abo_evidence = "".join(
        f'<tr><td class="rsid-cell"><a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(e["rsid"])}" '
        f'target="_blank" rel="noopener">{_esc(e["rsid"])}</a></td>'
        f'<td class="gt-cell">{_esc(e["gt"])}</td>'
        f'<td>{_esc(e["interpretation"])}</td></tr>'
        for e in a.get("evidence", []))

    hidden_o = ("<div style='margin-top:6px;color:#12467a'><strong>Carries a hidden O allele</strong> "
                "— you can pass either A or O to a child.</div>"
                if a.get("carries_hidden_O") else "")

    secretor_html = ""
    if s.get("available"):
        secretor_html = f"""
<h3 style="margin:18px 0 6px">FUT2 Secretor Status</h3>
<div style="border:1px solid var(--bdr);border-left:4px solid #12467a;background:var(--bg2);
     border-radius:6px;padding:11px 13px">
  <div style="font-weight:700">{_esc(s['secretor_status'])}</div>
  {"".join(f'<div style="color:var(--muted);font-size:.9em;margin-top:3px">{_esc(e["interpretation"])}</div>' for e in s['evidence'])}
</div>"""

    return f"""
<section class="bloodtype-section" id="blood-type">
<h2>Blood Type <span class="pro-pill">V6.10</span></h2>
<div style="display:flex;gap:20px;align-items:center;flex-wrap:wrap;margin:12px 0;
     background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1.5px solid {hero_color};
     border-radius:14px;padding:20px 24px">
  <div style="text-align:center;min-width:120px">
    <div style="font-size:.75em;color:var(--muted);text-transform:uppercase;letter-spacing:.05em">Predicted</div>
    <div style="font-size:3.4em;font-weight:800;color:{hero_color};line-height:1">{_esc(combined)}</div>
  </div>
  <div style="flex:1;min-width:220px">
    <div style="font-size:1.1em"><strong>ABO:</strong> {_esc(a['phenotype'])}
      <span style="color:var(--muted);font-size:.85em"> · genotype {_esc(a['genotype'])} · {_esc(a['confidence'])} confidence</span></div>
    <div style="font-size:1em;margin-top:2px"><strong>Rh(D):</strong> {_esc(d['status'])}
      <span style="color:var(--muted);font-size:.85em"> · {_esc(d['confidence'])} confidence</span></div>
    {hidden_o}
    <div style="color:var(--muted);font-size:.85em;margin-top:8px">{_esc(bt.get('population_context', ''))}</div>
  </div>
</div>
<h3 style="margin:14px 0 4px">ABO evidence</h3>
<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>
  <th>rsID</th><th>Genotype</th><th>Interpretation</th></tr></thead>
<tbody>{abo_evidence}</tbody></table></div>
<h3 style="margin:14px 0 4px">Rh(D) inference</h3>
<div style="color:var(--muted);font-size:.9em">{_esc(d['method'])}</div>
{secretor_html}
<div class="anc-caveat" style="margin-top:12px">{_esc(bt.get('disclaimer',''))}</div>
</section>
"""


def build_deep_ancestry_html(da: dict | None) -> str:
    """State-of-the-art deep ancestry — Neanderthal, ancient populations,
    N-S European axis, and haplogroup migration timelines."""
    if not da or not da.get("available"):
        return ""

    # ── Neanderthal panel ────────────────────────────────────────────────
    n = da.get("neanderthal") or {}
    neanderthal_html = ""
    if n.get("available"):
        rows = ""
        for v in n["variants"]:
            col = "#b3261e" if v["n_alleles"] >= 2 else ("#d29922" if v["n_alleles"] == 1 else "#8b949e")
            rows += (
                f'<tr><td class="rsid-cell"><a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(v["rsid"])}" '
                f'target="_blank" rel="noopener">{_esc(v["rsid"])}</a></td>'
                f'<td><strong>{_esc(v["gene"])}</strong><div style="color:#8a94a3;font-size:.82em">{_esc(v["trait"])}</div></td>'
                f'<td class="gt-cell">{_esc(v["genotype"])}</td>'
                f'<td>{_esc(v["neanderthal_allele"])}</td>'
                f'<td style="color:{col};font-weight:700">{v["n_alleles"]}</td></tr>')
        tier_color = {"Below average": "#3fb950", "Average non-African": "#8b949e",
                      "Above average": "#d29922", "High": "#b3261e"}.get(n["tier"], "#8b949e")
        neanderthal_html = f"""
<h3 style="margin:18px 0 6px">🦴 Neanderthal Introgression</h3>
<div style="display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin-bottom:8px;
     background:linear-gradient(135deg,#f4f8fc,#eef2f7);border:1px solid #dbe3ec;
     border-radius:10px;padding:14px 18px">
  <div style="text-align:center;min-width:110px">
    <div style="font-size:.75em;color:#8a94a3;text-transform:uppercase">Affinity</div>
    <div style="font-size:2.3em;font-weight:800;color:{tier_color}">~{n['approx_pct']}%</div>
    <div style="font-size:.8em;color:#8a94a3">{n['n_neanderthal_alleles']}/{n['max_possible']} tagged alleles</div>
  </div>
  <div style="flex:1;min-width:220px">
    <div style="font-weight:700;color:{tier_color}">{_esc(n['tier'])}</div>
    <div style="color:#4a5560;line-height:1.55">{_esc(n['tier_note'])}</div>
    <div style="color:#8a94a3;font-size:.78em;margin-top:4px">
      Non-Africans average ~2% genome-wide Neanderthal ancestry. This is a curated-panel
      <em>affinity</em> — not a genome-wide percentage — comparable across users of this tool.
      Consumer chips tag only a subset of the ~6,000 Neanderthal-introgressed SNPs.</div>
  </div>
</div>
<details><summary>{n['n_typed']} Neanderthal-tagged markers checked</summary>
<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>
  <th>rsID</th><th>Locus / trait</th><th>Genotype</th><th>N-allele</th><th>Dose</th></tr></thead>
<tbody>{rows}</tbody></table></div></details>
<div style="font-size:.78em;color:#9aa4b0;margin:6px 0 4px">📖 {_esc(n['citation'])}</div>"""

    # ── Ancient populations ─────────────────────────────────────────────
    ap = da.get("ancient_populations") or {}
    ancient_html = ""
    if ap.get("available"):
        cards = ""
        top_short = ap.get("top")
        for p in ap["populations"]:
            highlight = p["short"] == top_short
            border = "#12467a" if highlight else "#dbe3ec"
            cards += f"""
<div style="border:1.5px solid {border};border-radius:10px;padding:12px 14px;background:#fff;break-inside:avoid">
  <div style="display:flex;justify-content:space-between;gap:8px;align-items:baseline;flex-wrap:wrap">
    <span style="font-weight:700;color:#12467a">{_esc(p['name'])}</span>
    <span style="font-size:.85em;color:#8a94a3">{p['n_carried']}/{p['n_max']} alleles</span>
  </div>
  <div style="font-size:1.9em;font-weight:800;color:#12467a;margin:3px 0">{p['affinity']*100:.0f}%</div>
  {_bar(p['affinity'], 1.0, '#12467a')}
  <div style="font-size:.86em;color:#4a5560;line-height:1.5;margin-top:6px">{_esc(p['narrative'])}</div>
</div>"""
        ancient_html = f"""
<h3 style="margin:18px 0 6px">🏹 Ancient-Population Affinity</h3>
<p style="color:#556;margin:4px 0 8px">Your affinity to three major ancient European
gene pools, scored on the trait alleles ancient-DNA studies have traced to each:</p>
<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:10px">{cards}</div>
<div style="font-size:.78em;color:#9aa4b0;margin-top:4px">📖 {_esc(ap['citation'])}</div>"""

    # ── N-S European axis ───────────────────────────────────────────────
    ea = da.get("european_axis") or {}
    axis_html = ""
    if ea.get("available"):
        idx = ea["index"] * 100
        # Position a marker on a horizontal S↔N axis
        marker = f"""
<div style="position:relative;height:36px;margin:8px 0">
  <div style="position:absolute;top:14px;left:0;right:0;height:8px;border-radius:4px;
       background:linear-gradient(90deg,#e76f51,#e9c46a,#2a9d8f)"></div>
  <div style="position:absolute;left:{idx:.1f}%;top:6px;transform:translateX(-50%);
       width:20px;height:20px;border-radius:50%;background:#12467a;border:3px solid #fff;
       box-shadow:0 0 0 1px rgba(0,0,0,.15)"></div>
  <div style="position:absolute;left:0;top:32px;font-size:.75em;color:#8a94a3">Southern</div>
  <div style="position:absolute;right:0;top:32px;font-size:.75em;color:#8a94a3">Northern</div>
</div>"""
        axis_html = f"""
<h3 style="margin:18px 0 6px">🧭 Sub-Continental European Axis</h3>
<div style="border:1px solid #dbe3ec;border-radius:10px;padding:14px 16px;background:#fbfcfd">
  <div style="font-weight:700;color:#12467a">{_esc(ea['lean'])}
    <span style="font-weight:400;color:#8a94a3;font-size:.82em"> · index {ea['index']}</span></div>
  {marker}
  <div style="color:#4a5560;line-height:1.5;margin-top:4px">{_esc(ea['note'])}</div>
  <div style="font-size:.8em;color:#8a94a3;margin-top:8px">
    Uses LCT (lactase persistence — Yamnaya-derived, Northern-enriched), HERC2 (blue eyes),
    TYR and MC1R. A soft axis; most modern individuals are admixed across it.</div>
</div>"""

    # ── Haplogroup migration timeline ───────────────────────────────────
    tl = da.get("haplogroup_timeline") or {}
    tl_html = ""

    def _tl_card(entry, kind):
        if not entry:
            return ""
        story = entry.get("story") or ""
        return f"""
<div style="border-left:4px solid #12467a;background:#fbfcfd;padding:10px 14px;border-radius:0 8px 8px 0;margin:6px 0">
  <div style="font-weight:700;color:#12467a">{_esc(kind)} {_esc(entry['haplogroup'])}
    <span style="font-weight:400;color:#8a94a3;font-size:.85em">
      · TMRCA ~{entry['tmrca_kya']:g} kya · {_esc(entry['origin'])}</span></div>
  {f'<div style="color:#4a5560;line-height:1.55;margin-top:4px">{_esc(story)}</div>' if story else ''}
</div>"""

    if tl.get("y") or tl.get("mt"):
        tl_html = f"""
<h3 style="margin:18px 0 6px">🗺️ Haplogroup Migration Timeline</h3>
{_tl_card(tl.get('y'), 'Y-DNA (paternal)')}
{_tl_card(tl.get('mt'), 'mtDNA (maternal)')}
<div style="font-size:.78em;color:#9aa4b0;margin-top:2px">
  TMRCAs = time to most recent common ancestor, from ISOGG / YFull.
</div>"""

    return f"""
<section class="deep-ancestry-section" id="deep-ancestry" style="margin:8px 0">
<h2 style="font-size:1.35em;border-bottom:2px solid #e3e8ee;padding-bottom:4px;color:#12467a">
  Deep Ancestry <span style="font-size:.6em;color:#8a94a3;font-weight:400">
  · archaic + ancient-population + migration</span></h2>
<p style="color:#667;margin:6px 0 10px">Beyond the continental estimate: your Neanderthal
introgression, affinity to the three major ancient European gene pools (Yamnaya-Steppe,
Anatolian Neolithic, Western Hunter-Gatherer), a Northern-vs-Southern European axis, and
the migration timeline of your Y-DNA and mtDNA haplogroups.</p>
{neanderthal_html}
{ancient_html}
{axis_html}
{tl_html}
</section>
"""


def build_ancestry_html(anc: dict | None) -> str:
    if not anc or not anc.get("available"):
        return ""
    proportions = anc.get("proportions", {})
    sorted_props = anc.get("sorted_proportions", sorted(proportions.items(), key=lambda x: -x[1]))
    bars = ""
    for sp, p in sorted_props:
        pct = p * 100
        long_name = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
                     "SAS": "South Asian", "AMR": "Admixed American"}.get(sp, sp)
        bars += f"""
<div class="anc-row">
  <div class="anc-pop">{_esc(sp)} — {_esc(long_name)}</div>
  <div class="anc-bar-wrap">
    <div class="anc-bar anc-bar-{sp.lower()}" style="width:{pct:.1f}%"></div>
  </div>
  <div class="anc-pct">{pct:.1f}%</div>
</div>
"""
    plot_html = ""
    if anc.get("plot_png_b64"):
        plot_html = (
            f'<div class="anc-plot">'
            f'<img alt="Ancestry PCA scatter plot" '
            f'src="data:image/png;base64,{anc["plot_png_b64"]}" />'
            f'</div>'
        )
    primary = anc.get("primary_population")
    primary_long = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
                    "SAS": "South Asian", "AMR": "Admixed American"}.get(primary, primary)
    n_indep = anc.get("n_aims_independent", anc.get("n_aims_used", 0))
    n_exp = anc.get("n_aims_expected", "—")
    conf_badge = _confidence_badge(anc.get("confidence"), anc.get("confidence_note"))

    # DEMOTED heuristic mode: show a prominent banner and skip the "best match"
    # ancestry claim entirely. A selection-confounded ~10-marker panel cannot call
    # ancestry (it mis-classified a European sample as East Asian); the markers are
    # presented as individual pigmentation/trait signals instead.
    suppressed_banner = ""
    if anc.get("ancestry_call_suppressed"):
        suppressed_banner = (
            '<div class="anc-ambiguous"><strong>No ancestry call.</strong> '
            + _esc(anc.get("suppression_reason", "")) +
            ' The markers below are individual pigmentation and trait variants, '
            'shown for interest — not a population assignment.</div>'
        )

    ambiguous_banner = ""
    if anc.get("ambiguous") and not anc.get("ancestry_call_suppressed"):
        runner = anc.get("runner_up_population")
        runner_long = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
                       "SAS": "South Asian", "AMR": "Admixed American"}.get(runner, runner)
        ambiguous_banner = (
            f'<div class="anc-ambiguous">⚠️ <strong>Ambiguous call.</strong> '
            f'{_esc(primary_long)} and {_esc(runner_long)} fit your markers almost '
            'equally well; this small panel cannot reliably tell them apart. Treat '
            'the result as "best guess," not a determination.</div>'
        )

    # Evidence margin + the actual markers behind the call. The module labels
    # evidence_margin_nats "the honest measure of confidence"; used_aims are the
    # AIMs typed (some flagged LD-redundant and not double-counted). All computed
    # by ancestry_pca — previously only leaked textually in the ambiguous branch.
    margin = anc.get("evidence_margin_nats")
    margin_html = ""
    if margin is not None:
        margin_html = (
            f'<div class="anc-margin"><em>Evidence margin:</em> {margin} nats — '
            f"the log-likelihood gap between the best and runner-up population "
            f"(larger = more confident; &lt;2.3 nats ≈ within 10×, treated as "
            f"ambiguous).</div>"
        )
    aims = anc.get("used_aims") or []
    aims_html = ""
    if aims:
        def _aim_role(a: dict) -> str:
            if a.get("palindromic"):
                return "palindrome — shown only, cannot be strand-oriented"
            if not a.get("counted", True):
                return "LD-redundant (not double-counted)"
            if a.get("selection"):
                return "counted (selection-influenced marker)"
            return "counted"
        arows = "".join(
            f'<tr class="{"" if a.get("counted", True) else "anc-aim-redundant"}">'
            f'<td class="rsid-cell">'
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(a.get("rsid",""))}" '
            f'target="_blank" rel="noopener">{_esc(a.get("rsid",""))}</a></td>'
            f'<td><strong>{_esc(a.get("gene",""))}</strong></td>'
            f'<td class="gt-cell">{_esc(a.get("genotype",""))}</td>'
            f'<td>{_esc(a.get("effect_allele",""))}</td>'
            f'<td>{a.get("dosage","")}</td>'
            f'<td>{_esc(_aim_role(a))}</td>'
            f"</tr>"
            for a in aims
        )
        aims_html = (
            f'<details class="anc-aims"><summary>Ancestry-informative markers '
            f"used ({len(aims)})</summary>"
            f'<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>'
            f"<th>rsID</th><th>Gene</th><th>Genotype</th><th>Effect allele</th>"
            f"<th>Copies</th><th>Role</th></tr></thead><tbody>{arows}</tbody>"
            f"</table></div></details>"
        )

    crosscheck_html = _build_crosscheck_html(anc.get("haplogroup_crosscheck"))
    n_palin = anc.get("n_aims_palindromic", 0)
    palin_note = ""
    if n_palin:
        palin_note = (
            f'<div class="anc-margin"><em>Strand safety:</em> {n_palin} '
            f"palindromic marker(s) (A/T or C/G, e.g. SLC45A2) were displayed but "
            f"excluded from scoring — their strand can't be recovered from "
            f"genotype alone, and mis-orienting one was a known cause of spurious "
            f"non-European calls. Remaining markers are read strand-aware.</div>"
        )

    return f"""
<section class="anc-section" id="ancestry">
<h2>Ancestry Estimation <span class="pro-pill">V3</span></h2>
<p class="anc-intro">
{_esc(anc.get('method', 'Genotype-based ancestry estimate'))}.
Based on <strong>{n_indep}</strong> independent ancestry-informative markers
(of {n_exp} on this panel).
</p>
{suppressed_banner}
{conf_badge}
{margin_html}
{palin_note}
{crosscheck_html}
{ambiguous_banner}
<div class="anc-caveat">
<strong>Informational use only — this is not a genealogical or clinical
ancestry test.</strong> It is a rough estimate next to commercial products
(which use tens of thousands of markers). The bars below are <em>relative
affinities</em> — the probability that your genotype best matches each
<em>single</em> 1000 Genomes population — <strong>not admixture
proportions</strong>; they should not be read as "% of your DNA." It is most
meaningful for recent single-continent ancestry. Note that "Admixed American
(AMR)" is itself a recently-admixed reference group, so it can appear as a
spurious match for intermediate genotypes.
</div>
{plot_html}
<div class="anc-bars-label">Relative affinity to each reference population (not admixture %):</div>
<div class="anc-bars">{bars}</div>
{aims_html}
</section>
"""


# ── Detoxification & metal / oxidative sections ───────────────────────────────

_DETOX_IMPACT_STYLE = {
    "higher-load":      ("#f85149", "activates faster"),
    "reduced-clearance":("#f85149", "clears slower"),
    "reduced":          ("#d29922", "reduced"),
    "intermediate":     ("#8b949e", "intermediate"),
    "typical":          ("#3fb950", "typical"),
    "protective":       ("#3fb950", "protective"),
}


def build_urologic_html(ur: dict | None) -> str:
    """Urologic & Genitourinary panel — OAB, BPH, prostate cancer, kidney stones,
    testicular / reproductive, DHT metabolism."""
    if not ur or not ur.get("available"):
        return ""
    impact_color = {
        "higher-load": "#f85149", "reduced": "#d29922",
        "reduced-clearance": "#f85149", "intermediate": "#8b949e",
        "typical": "#3fb950", "protective": "#3fb950",
    }
    conf_col = {"high": "#3fb950", "moderate": "#d29922", "low": "#8b949e"}

    domains_html = ""
    for cat in ur.get("categories", []):
        items = ur["by_category"].get(cat, [])
        rows = ""
        for f in items:
            border = impact_color.get(f["impact"], "#8b949e")
            rows += f"""
<div style="border:1px solid var(--bdr);border-left:4px solid {border};border-radius:6px;
     padding:12px;margin:8px 0;background:var(--bg2)">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <div style="font-weight:700">{_esc(f['trait'])}
      <span style="font-weight:400;color:var(--muted);font-size:.85em"> · {_esc(f['gene'])}</span></div>
    <div><span style="color:{border};font-weight:600;font-size:.85em">{_esc(f['impact'])}</span>
      <span style="color:{conf_col.get(f['confidence'], '#8b949e')};font-size:.8em"> · {_esc(f['confidence'])} confidence</span></div>
  </div>
  <div style="margin:6px 0;line-height:1.5">{_esc(f['result'])}</div>
  <div style="font-size:.9em;color:var(--acc2)"><strong>Action:</strong> {_esc(f['action'])}</div>
  <div style="font-size:.78em;color:var(--dim);margin-top:4px;font-family:var(--mono)">
    {_esc(f['rsid'])} · genotype {_esc(f['genotype'])} · {_esc(f['evidence'])}</div>
</div>"""
        domains_html += f'<h3 style="margin-top:18px">{_esc(cat)}</h3>{rows}'

    return f"""
<section class="urologic-section" id="urologic">
<h2>Urologic &amp; Genitourinary Panel <span class="pro-pill">V6.8</span></h2>
<p class="anc-intro">
Genotype-based screen for common urologic conditions — <strong>bladder function
(OAB)</strong>, <strong>BPH &amp; 5α-reductase</strong>,
<strong>prostate cancer</strong> (including HOXB13, the top hereditary marker),
<strong>kidney stones</strong>, <strong>testicular germ-cell cancer</strong>,
and <strong>androgen bioavailability</strong>. {ur['n_findings']} findings across
{len(ur.get('categories', []))} sub-panels.
</p>
<div class="urologic-domains">{domains_html}</div>
<div class="anc-caveat" style="margin-top:14px">
<strong>Consumer-chip screen, not a clinical diagnostic.</strong> A HOXB13 G84E
positive on this chip is worth confirming with clinical-grade sequencing; family
history and PSA/urine studies remain the primary clinical tools. Testicular
self-exam remains the single highest-yield screening action here.
</div>
</section>
"""


def build_detox_html(dx: dict | None) -> str:
    """Detoxification & Environmental Resilience — smoke / PAH / heavy metals."""
    if not dx or not dx.get("available"):
        return ""

    sr = dx.get("smoke_resilience", {})
    sr_color = sr.get("color", "#8b949e")
    mismatch_badge = ""
    if sr.get("activate_but_dont_clear"):
        mismatch_badge = (
            '<div style="margin-top:8px;padding:8px 12px;border-radius:6px;'
            'background:rgba(248,81,73,.12);border:1px solid #f8514955;font-size:.9em">'
            '⚠ <strong>Activate-but-don’t-clear pattern detected:</strong> your '
            'Phase I enzymes bioactivate combustion toxicants faster than your '
            'Phase II / antioxidant system clears them — the genotype combination '
            'most worth protecting against smoke exposure.</div>'
        )

    score_card = f"""
<div style="border:1.5px solid {sr_color};border-radius:10px;padding:16px;margin:14px 0;
     background:linear-gradient(180deg,var(--bg3),var(--bg2))">
  <div style="display:flex;align-items:baseline;gap:12px;flex-wrap:wrap">
    <div style="font-size:1.4em;font-weight:800;color:{sr_color}">{_esc(sr.get('tier','—'))}</div>
    <div style="color:var(--muted);font-size:.9em">Wildfire-smoke resilience index</div>
  </div>
  <p style="margin:8px 0 0;line-height:1.55">{_esc(sr.get('headline',''))}</p>
  <div style="display:flex;gap:16px;margin-top:10px;flex-wrap:wrap;font-size:.86em;color:var(--muted)">
    <span>Phase I activation hits: <strong>{sr.get('activation_hits',0)}</strong></span>
    <span>Phase II clearance deficits: <strong>{sr.get('clearance_deficit_hits',0)}</strong></span>
    <span>Antioxidant deficits: <strong>{sr.get('antioxidant_deficit_hits',0)}</strong></span>
  </div>
  {mismatch_badge}
</div>
"""

    # Findings grouped by domain
    domains_html = ""
    by_cat = dx.get("by_category", {})
    for cat in dx.get("category_order", list(by_cat.keys())):
        items = by_cat.get(cat, [])
        if not items:
            continue
        rows = ""
        for f in items:
            color, label = _DETOX_IMPACT_STYLE.get(f.get("impact"), ("#8b949e", f.get("impact", "")))
            rows += f"""
<div style="border:1px solid var(--bdr);border-left:4px solid {color};border-radius:6px;
     padding:12px;margin:8px 0;background:var(--bg2)">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <div style="font-weight:700">{_esc(f.get('trait',''))}
      <span style="font-weight:400;color:var(--muted);font-size:.85em"> · {_esc(f.get('gene',''))}</span></div>
    <div><span style="color:{color};font-weight:600;font-size:.85em">{_esc(label)}</span>
      <span style="color:var(--dim);font-size:.8em"> · {_esc(f.get('confidence',''))} confidence</span></div>
  </div>
  <div style="margin:6px 0;line-height:1.5">{_esc(f.get('result',''))}</div>
  <div style="font-size:.9em;color:var(--acc2)"><strong>Action:</strong> {_esc(f.get('action',''))}</div>
  <div style="font-size:.78em;color:var(--dim);margin-top:4px;font-family:var(--mono)">
    {_esc(f.get('rsid',''))} · genotype {_esc(f.get('genotype',''))} · {_esc(f.get('evidence',''))}</div>
</div>
"""
        domains_html += (
            f'<h3 style="margin-top:18px">{_esc(cat)}</h3>{rows}'
        )

    # Personalised protocol
    proto = dx.get("protocol", {})

    def _proto_block(title, items, key_item="item", key_detail="detail"):
        if not items:
            return ""
        lis = ""
        for it in items:
            emph = it.get("emphasis")
            star = (' <span style="background:var(--red);color:#fff;border-radius:4px;'
                    'padding:0 6px;font-size:.7em;vertical-align:middle">PRIORITY</span>'
                    if emph == "high" else "")
            lis += (
                f'<li style="margin:8px 0"><strong>{_esc(it.get(key_item,""))}</strong>{star}'
                f'<div style="color:var(--muted);font-size:.92em;line-height:1.5;margin-top:2px">'
                f'{_esc(it.get(key_detail,""))}</div></li>'
            )
        return (f'<h3 style="margin-top:18px">{_esc(title)}</h3>'
                f'<ul style="list-style:none;padding-left:0;margin:0">{lis}</ul>')

    protocol_html = (
        _proto_block("Your personalised nutrition levers", proto.get("nutrition"))
        + _proto_block("Behavioural protection during smoke events", proto.get("behavioural"))
        + _proto_block("Heavy-metal exposure & measurement", proto.get("metal"))
    )

    michigan = dx.get("michigan_context", "")
    michigan_html = (
        f'<div style="border-left:4px solid var(--acc);background:var(--bg2);'
        f'padding:12px 14px;border-radius:6px;margin:14px 0;line-height:1.55">'
        f'📍 {_esc(michigan)}</div>' if michigan else ""
    )

    return f"""
<section class="detox-section" id="detoxification">
<h2>Detoxification &amp; Environmental Resilience <span class="pro-pill">V10</span></h2>
<p class="anc-intro">
How your genome handles the toxicant load from <strong>wildfire / wood / tobacco
smoke</strong> and <strong>heavy metals</strong> — the balance between Phase I
activation, Phase II conjugation, and the antioxidant response, plus a
genotype-tuned action plan.
</p>
{michigan_html}
{score_card}
<div class="detox-domains">{domains_html}</div>
<div class="detox-protocol">
  <h2 style="font-size:1.25em;margin-top:24px">Your Personalised Detox Protocol</h2>
  {protocol_html}
</div>
<div class="anc-caveat" style="margin-top:16px">
<strong>Informational, not diagnostic.</strong> These are consumer-chip
genotype tendencies, not measurements of your toxicant burden or detox capacity.
A true GSTM1/GSTT1-null status needs a PCR/CNV assay; blood-lead, urine
heavy-metals and lung function are lab/clinical tests. The behavioural
protections during smoke events benefit everyone regardless of genotype — the
genetics only personalise <em>emphasis</em>. Discuss supplements (especially NAC
and selenium) with a clinician before starting.
</div>
</section>
"""


def build_metal_oxidative_html(mx: dict | None) -> str:
    """Metal-handling, oxidative-defense & neurodegeneration panel (wires in the
    previously-unrendered metal_oxidative module)."""
    if not mx or not mx.get("predictions"):
        return ""
    conf_color = {"high": "#3fb950", "moderate": "#d29922", "low": "#8b949e"}
    domains_html = ""
    for cat in mx.get("categories", list(mx.get("by_category", {}).keys())):
        items = mx.get("by_category", {}).get(cat, [])
        if not items:
            continue
        rows = ""
        for f in items:
            c = conf_color.get(f.get("confidence"), "#8b949e")
            clinical = ""
            if f.get("clinical_variant"):
                cv = f["clinical_variant"]
                clinical = (
                    f'<div style="margin-top:6px;padding:6px 10px;border-radius:6px;'
                    f'background:rgba(248,81,73,.10);border:1px solid #f8514955;font-size:.85em">'
                    f'🧬 <strong>Clinically-reportable:</strong> {_esc(cv.get("gene",""))} '
                    f'{_esc(cv.get("variant",""))} — {_esc(cv.get("clinical_significance",""))} '
                    f'({_esc(cv.get("inheritance",""))})</div>'
                )
            rows += f"""
<div style="border:1px solid var(--bdr);border-left:4px solid {c};border-radius:6px;
     padding:12px;margin:8px 0;background:var(--bg2)">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
    <div style="font-weight:700">{_esc(f.get('trait',''))}</div>
    <div style="color:var(--dim);font-size:.8em">{_esc(f.get('confidence',''))} confidence</div>
  </div>
  <div style="margin:6px 0;line-height:1.5">{_esc(f.get('result',''))}</div>
  <div style="font-size:.9em;color:var(--acc2)"><strong>Action:</strong> {_esc(f.get('action',''))}</div>
  <div style="font-size:.78em;color:var(--dim);margin-top:4px;font-family:var(--mono)">{_esc(f.get('evidence',''))}</div>
  {clinical}
</div>
"""
        domains_html += f'<h3 style="margin-top:18px">{_esc(cat)}</h3>{rows}'

    return f"""
<section class="metalox-section" id="metal-oxidative">
<h2>Metal Handling, Oxidative Defense &amp; Neurodegeneration <span class="pro-pill">V9</span></h2>
<p class="anc-intro">
A functional look at copper/iron/zinc handling, red-cell and cellular
antioxidant defenses, and the LRRK2/GBA Parkinson's-risk loci — a complement to
the Detoxification section above and the Wellness oxidative-stress trait.
</p>
<div class="metalox-domains">{domains_html}</div>
<div class="anc-caveat" style="margin-top:14px">
<strong>Mostly research-grade.</strong> Metallothionein, ZIP-transporter and
GST proxy signals are hints, not clinical metal/detox tests. HFE C282Y
homozygosity and G6PD deficiency are the higher-confidence, clinically
actionable exceptions and are flagged as such.
</div>
</section>
"""


# ── V5: premium section renderers ────────────────────────────────────────────

def build_hla_html(h: dict | None) -> str:
    if not h:
        return ""
    cards = ""
    for a in h.get("alleles", []):
        status = a["status"]
        if status == "no_tag_available":
            badge = '<span class="hla-status hla-na">No tag available</span>'
            cards += (
                f'<div class="hla-card">'
                f'<div class="hla-head"><span class="hla-allele">{_esc(a["allele"])}</span>{badge}</div>'
                f'<div class="hla-detail">No common tag SNP available for chip-based imputation. '
                f'Clinical HLA typing required for accurate detection.</div>'
                f'</div>'
            )
            continue
        if status == "untested":
            badge = '<span class="hla-status hla-na">Untested (tag SNPs not on chip)</span>'
            # Name the specific tag SNPs that would have been needed, so the gap
            # is concrete rather than an opaque "untested".
            untested = a.get("untested_tags") or []
            untested_html = ""
            if untested:
                chips = "".join(
                    f'<code class="hla-tag hla-tag-missing">{_esc(t.get("rsid", ""))}</code>'
                    for t in untested
                    if isinstance(t, dict) and t.get("rsid")
                )
                if chips:
                    untested_html = (
                        f'<div class="hla-detail">Tag SNP(s) not typed on this '
                        f"chip: {chips}</div>"
                    )
            cards += (
                f'<div class="hla-card">'
                f'<div class="hla-head"><span class="hla-allele">{_esc(a["allele"])}</span>{badge}</div>'
                f"{untested_html}"
                f'</div>'
            )
            continue
        if status == "negative":
            badge = '<span class="hla-status hla-neg">Negative</span>'
        elif status == "homozygous":
            badge = '<span class="hla-status hla-pos">Homozygous</span>'
        else:
            badge = '<span class="hla-status hla-pos">Carrier</span>'

        clinical_html = ""
        if status in ("carrier (heterozygous)", "homozygous"):
            for label, text in a.get("clinical", []):
                clinical_html += (
                    f'<div class="hla-clinical">'
                    f'<div class="hla-clinical-label">{_esc(label)}</div>'
                    f'<div class="hla-clinical-text">{_esc(text)}</div>'
                    f'</div>'
                )

        tag_chips = ""
        for t in a.get("called_tags", []):
            tag_chips += (
                f'<code class="hla-tag" title="{_esc(t["ld_quality"])}">'
                f'{_esc(t["rsid"])}={_esc(t.get("genotype","?"))}</code>'
            )

        caveat = ""
        if a.get("chip_caveat"):
            caveat = ('<div class="hla-caveat">'
                      '⚠ This allele cannot be reliably ruled in or out from common-SNP tags. '
                      'Direct HLA typing required.</div>')

        cards += f"""
<div class="hla-card">
  <div class="hla-head">
    <span class="hla-allele">{_esc(a["allele"])}</span>
    {badge}
    <span class="hla-conf">conf: {_esc(a.get("confidence","-"))}</span>
  </div>
  <div class="hla-freq">Population frequency — {_esc(a.get("frequency","—"))}</div>
  <div class="hla-tags">{tag_chips}</div>
  {clinical_html}
  {caveat}
</div>
"""

    transplant_html = (
        f'<details class="hla-transplant"><summary>Transplant compatibility context</summary>'
        f'<p>{_esc(h.get("transplant_context",""))}</p></details>'
    )

    infection_lines = "".join(
        f'<li><strong>{_esc(it["hla"])}</strong>: {_esc(it["association"])}</li>'
        for it in h.get("infection_context", [])
    )
    infection_html = (
        f'<details class="hla-infection"><summary>HLA & infection response (general context)</summary>'
        f'<ul>{infection_lines}</ul></details>'
    )

    return f"""
<section class="hla-section" id="hla-immune">
<h2>Immune Genotype — HLA Imputation <span class="pro-pill">V5</span></h2>
<p class="hla-intro">
HLA genes are the most polymorphic in the human genome and govern immune
response — including drug hypersensitivity, autoimmune disease susceptibility,
infection control, and transplant compatibility. This module imputes
clinically actionable HLA alleles from surrounding tag SNPs.
{h.get("n_alleles_called", 0)}/{h.get("n_alleles_tested", 0)} alleles
called; {h.get("n_carrier_alleles", 0)} carrier finding(s).
</p>
<div class="hla-caveat-block">
<strong>Caveats:</strong> Tag-SNP imputation is approximate — confirm any
positive call by direct HLA typing before clinical action (especially for
drug-hypersensitivity decisions). HLA-B*15:02 has no reliable common-SNP tag
outside Asian populations.
</div>
<div class="hla-grid">{cards}</div>
{transplant_html}
{infection_html}
</section>
"""


def build_roh_html(r: dict | None) -> str:
    if not r:
        return ""
    # Degrade gracefully if the ROH analysis failed or returned a partial result
    # (e.g. a degenerate/near-empty genome, or an upstream error). Missing core
    # stats must not crash the whole report — skip the section instead.
    if r.get("available") is False or not all(
            k in r for k in ("n_runs", "total_roh_mb", "f_roh",
                             "n_short", "n_medium", "n_long")):
        return ""
    runs = r.get("runs", [])
    ideogram_svg = render_ideogram_svg(runs) if (render_ideogram_svg and runs) else ""

    rows = ""
    for run in sorted(runs, key=lambda x: -x["length_mb"])[:25]:
        genes = ""
        if run.get("disease_genes"):
            genes = " · " + ", ".join(
                f'<strong>{_esc(g["gene"])}</strong> ({_esc(g["condition"])})'
                for g in run["disease_genes"]
            )
        cls = ("roh-long" if run["length_mb"] >= 8
               else ("roh-med" if run["length_mb"] >= 2 else "roh-short"))
        rows += (
            f'<tr class="{cls}">'
            f'<td>chr{_esc(run["chrom"])}</td>'
            f'<td>{run["start_bp"]/1e6:.1f}-{run["end_bp"]/1e6:.1f} Mb</td>'
            f'<td>{run["length_mb"]:.2f} Mb</td>'
            f'<td>{run["n_snps"]}</td>'
            f'<td class="sum-cell">{genes or "—"}</td>'
            f'</tr>'
        )

    return f"""
<section class="roh-section" id="runs-of-homozygosity">
<h2>Runs of Homozygosity (ROH) <span class="pro-pill">V5</span></h2>
<p class="roh-intro">
Contiguous stretches of identical alleles on both chromosomes reveal
population history, geographic isolation, and at the longest scales recent
parental relatedness. This is detected by scanning the full chip data for
extended homozygous segments.
</p>
<div class="roh-stats">
  <div class="roh-stat"><div class="roh-n">{r['n_runs']}</div><div class="roh-l">Total ROH</div></div>
  <div class="roh-stat"><div class="roh-n">{r['total_roh_mb']:.1f} Mb</div><div class="roh-l">Total ROH length</div></div>
  <div class="roh-stat"><div class="roh-n">{r['f_roh']:.4f}</div><div class="roh-l">F_ROH coefficient</div></div>
  <div class="roh-stat"><div class="roh-n">{r['n_short']} / {r['n_medium']} / {r['n_long']}</div><div class="roh-l">Short / Med / Long</div></div>
</div>
<div class="roh-context">{_esc(r.get("population_context",""))}</div>
{f'<div class="roh-ideogram">{ideogram_svg}</div>' if ideogram_svg else ''}
<details class="roh-table">
  <summary>Largest {min(25, len(runs))} ROH regions</summary>
  <div class="tbl-wrap">
  <table class="snp-tbl">
    <thead><tr><th>Chr</th><th>Position</th><th>Length</th><th>SNPs</th><th>Disease-relevant genes</th></tr></thead>
    <tbody>{rows or '<tr><td colspan="5">No ROH detected at default thresholds.</td></tr>'}</tbody>
  </table>
  </div>
</details>
</section>
"""


def build_local_ancestry_html(la: dict | None) -> str:
    if not la or not la.get("available"):
        return ""
    svg = ""
    if render_chromosome_painting_svg:
        svg = render_chromosome_painting_svg(la)
    by_call = la.get("by_call", {})
    by_call_html = ""
    total = sum(by_call.values()) or 1
    for sp, n in sorted(by_call.items(), key=lambda x: -x[1]):
        by_call_html += (
            f'<div class="la-call"><span class="la-sp la-sp-{sp.lower()}">{_esc(sp)}</span> '
            f'{n} windows ({100*n/total:.0f}%)</div>'
        )
    deviant_html = ""
    if la.get("deviant_segments"):
        items = ""
        for w in la["deviant_segments"][:10]:
            genes = ""
            if w.get("genes"):
                genes = " — " + ", ".join(_esc(g) for g in w["genes"])
            # Per-window quantitative evidence: the confidence gap and the
            # per-superpopulation log-likelihoods that produced the call.
            ll = w.get("log_liks") or {}
            ll_txt = ""
            if ll:
                ranked = sorted(ll.items(), key=lambda kv: -kv[1])
                ll_txt = " · ".join(f"{_esc(sp)} {v:.1f}" for sp, v in ranked)
            gap_txt = f", gap {w['gap']:.1f}" if w.get("gap") is not None else ""
            aims_txt = f", {w['n_aims']} AIMs" if w.get("n_aims") is not None else ""
            evidence = (
                f'<div class="la-ll">log-lik: {ll_txt}{gap_txt}{aims_txt}</div>'
                if ll_txt
                else ""
            )
            items += (
                f'<li>chr{_esc(w["chrom"])} {w["start_mb"]:.0f}-{w["end_mb"]:.0f} Mb: '
                f'<span class="la-sp la-sp-{w["call"].lower()}">{_esc(w["call"])}</span> '
                f'({w.get("confidence","?")}){genes}{evidence}</li>'
            )
        deviant_html = (
            f'<details class="la-deviant"><summary>Ancestry-discordant segments ({la["n_deviant"]})</summary>'
            f'<ul>{items}</ul></details>'
        )
    return f"""
<section class="la-section" id="local-ancestry">
<h2>Local Ancestry — Chromosome Painting <span class="pro-pill">V5</span></h2>
<p class="la-intro">
Per-window estimation of which ancestral superpopulation each chromosomal
segment most likely came from. {la['n_windows_called']}/{la['n_windows_total']}
windows could be confidently called given the AIM panel available.
</p>
<div class="la-caveat-block">{_esc(la.get("limitations",""))}</div>
<div class="la-summary">{by_call_html}</div>
{f'<div class="la-paint">{svg}</div>' if svg else ''}
{deviant_html}
</section>
"""


def build_phewas_html(p: dict | None) -> str:
    if not p:
        return ""
    by_cat = p.get("by_category", {})
    cats_html = ""
    for cat in sorted(by_cat.keys()):
        trait_names = by_cat[cat]
        rows = ""
        for tname in trait_names:
            t = p["traits"][tname]
            res = t["result"]
            if res["status"] != "ok":
                continue
            unit = t.get("unit", "")
            tier_cls = res["tier_class"]
            pct = res["percentile"]
            val = res.get("predicted_value")

            # Evidence behind the estimate: GWAS reference, the standardised
            # Z-score against the reference distribution, and the specific
            # variants that drove the score (rsID, effect allele, your dose,
            # published per-allele effect). All produced by phewas._score_trait.
            bits = []
            ref = t.get("reference")
            if ref:
                bits.append(f'<div class="ph-ref"><em>GWAS source:</em> {_esc(ref)}</div>')
            z = res.get("z_score")
            if z is not None:
                bits.append(
                    f'<div class="ph-z"><em>Z-score:</em> {z:+.2f} '
                    f"(reference mean {t.get('mean','—')} {_esc(unit)}, "
                    f"SD {t.get('sd','—')})</div>"
                )
            uv = res.get("used_variants") or []
            if uv:
                vrows = "".join(
                    f'<tr><td class="rsid-cell">'
                    f'<a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(u.get("rsid",""))}" '
                    f'target="_blank" rel="noopener">{_esc(u.get("rsid",""))}</a></td>'
                    f'<td>{_esc(u.get("effect_allele",""))}</td>'
                    f'<td>{u.get("dose",0)}</td>'
                    f'<td>{u.get("beta",0.0):+.3f}</td></tr>'
                    for u in uv
                )
                bits.append(
                    f'<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>'
                    f"<th>rsID</th><th>Effect allele</th><th>Your copies</th>"
                    f"<th>Per-allele β</th></tr></thead><tbody>{vrows}</tbody>"
                    f"</table></div>"
                )
            miss = res.get("missing_variants") or []
            if miss:
                mtxt = ", ".join(_esc(m) for m in miss[:20])
                extra = f" +{len(miss) - 20} more" if len(miss) > 20 else ""
                bits.append(
                    f'<div class="ph-missing">{len(miss)} panel SNP(s) not typed '
                    f"(excluded from the score): {mtxt}{extra}.</div>"
                )
            detail = ""
            if bits:
                detail = (
                    f'<details class="ph-detail"><summary>Evidence &amp; '
                    f'contributing variants</summary>{"".join(bits)}</details>'
                )

            # Variance the panel explains of the ACTUAL trait, and the honest
            # trait-prediction percentile (regresses to the mean as the panel weakens).
            r2 = res.get("variance_explained_pct", 0.0)
            signal = res.get("signal_strength", "negligible")
            trait_pct = res.get("trait_percentile", 50.0)
            r2_cls = ("sig-modest" if signal == "modest"
                      else "sig-weak" if signal == "weak" else "sig-negligible")
            explains_cell = (f'<span class="{r2_cls}">{r2:.2f}%</span>'
                             + ('' if signal != "negligible"
                                else '<br><span style="font-size:.8em;color:#a06800">'
                                     'marker only — not a trait measurement</span>'))
            rows += (
                f'<tr>'
                f'<td>{_esc(tname)}{detail}</td>'
                f'<td><span class="prs-tier {tier_cls}">{_esc(res["tier"])}</span></td>'
                f'<td>{pct:.0f}th</td>'
                f'<td>{explains_cell}</td>'
                f'<td class="gt-cell">{val} {_esc(unit)} '
                f'<span style="color:#8a94a3">(~{trait_pct:.0f}th)</span></td>'
                f'<td>{res["n_used"]}/{res["n_total"]}</td>'
                f'</tr>'
            )
        if not rows:
            continue
        cats_html += f"""
<details class="phewas-cat" open>
  <summary><span class="phewas-cat-name">{_esc(cat)}</span>
    <span class="phewas-cat-count">{len([t for t in trait_names if p['traits'][t]['result']['status']=='ok'])} traits</span></summary>
  <div class="tbl-wrap">
  <table class="snp-tbl">
    <thead><tr><th>Trait</th><th>Marker-score tier</th><th>Marker-score %ile</th>
      <th>Explains (of real trait)</th><th>Predicted value (trait %ile)</th>
      <th>Coverage</th></tr></thead>
    <tbody>{rows}</tbody>
  </table></div>
</details>
"""

    headline_html = ""
    if p.get("headline"):
        def _hl_item(h):
            # Pull R² for this trait so an extreme MARKER score can't read as an
            # extreme TRAIT value (e.g. "99th percentile brain volume" from a single
            # SNP explaining 0.1% of it).
            res = (p.get("traits", {}).get(h["trait"], {}) or {}).get("result", {})
            r2 = res.get("variance_explained_pct")
            r2_txt = (f' · explains only {r2:.2f}% of the actual trait'
                      if r2 is not None and r2 < 5 else '')
            return (
                f'<li><strong>{_esc(h["trait"])}</strong> — '
                f'high genetic-marker score ({h["percentile"]:.0f}th percentile on the '
                f'marker panel{r2_txt})</li>')
        items = "".join(_hl_item(h) for h in p["headline"][:15])
        headline_html = (
            f'<div class="phewas-headline"><strong>Notable marker-score extremes '
            f'(genotype ranks, not trait measurements):</strong>'
            f'<ul>{items}</ul></div>')

    return f"""
<section class="phewas-section" id="phewas">
<h2>Predicted Biomarker & Trait Profile <span class="pro-pill">V5</span></h2>
<p class="phewas-intro">
Phenome-wide scoring across {p['n_traits']} traits using curated GWAS effect
sizes — biomarkers (lipids, HbA1c, CRP, vitamin D, hormones), anthropometric,
hematological, behavioural. {p['n_scored']}/{p['n_traits']} traits had sufficient
SNP coverage to score.
</p>
<div class="phewas-intro" style="background:#fff8e6;border:1px solid #f0e0a8;border-radius:8px;padding:10px 14px;font-size:.88em">
<strong>How to read this — two different numbers.</strong> The <em>marker-score
percentile</em> is where your genotype ranks on each small SNP panel — a real fact
about which alleles you carry. The <em>"Explains"</em> column is how much of the
<em>actual</em> trait that panel accounts for; here it ranges from ~6% down to under
0.1%. So a high marker-score percentile does <strong>not</strong> mean a high value of
the real trait — e.g. a 99th-percentile brain-volume marker score comes from a single
variant that explains ~0.1% of actual brain size. The predicted-value column shows the
honest trait estimate (which regresses toward average precisely because these panels
are weak). These are genetic-propensity signals against a European reference, not
measurements, and not medical advice.
</div>
{headline_html}
{cats_html}
</section>
"""


def build_mr_html(m: dict | None) -> str:
    if not m or not m.get("findings"):
        return ""
    items = ""
    sig = sorted(
        [f for f in m["findings"] if f["status"] == "ok"],
        key=lambda f: abs(f.get("outcome_shift_log_or") or 0),
        reverse=True,
    )
    for f in sig:
        rr = f["outcome_relative_risk"]
        rr_class = "tier-high" if rr > 1.20 else ("tier-low" if rr < 0.85 else "tier-average")
        pct = f["exposure_percentile"]
        items += f"""
<div class="mr-finding">
  <div class="mr-head">
    <span class="mr-exposure">{_esc(f["exposure"])}</span>
    <span class="mr-arrow">→</span>
    <span class="mr-outcome">{_esc(f["outcome"])}</span>
  </div>
  <div class="mr-stats">
    <span>Your exposure PRS: <strong>{f["exposure_z"]:+.2f}</strong> SD ({pct:.0f}th pct)</span>
    <span>Causal projection: <span class="mr-rr {rr_class}">RR ≈ {rr}</span></span>
    <span class="mr-coverage">{f["n_used"]} SNPs scored</span>
  </div>
  <div class="mr-explain">{_esc(f["explanation"])}</div>
  <div class="mr-ref"><em>Source: {_esc(f["reference"])} · MR β = {f["mr_beta"]:+.2f} {_esc(f["mr_unit"])}</em></div>
</div>
"""
    return f"""
<section class="mr-section" id="mendelian-randomization">
<h2>Mendelian Randomization — Causal Projections <span class="pro-pill">V5</span></h2>
<p class="mr-intro">
For key exposure → outcome pairs with published MR causal estimates,
we project the outcome shift implied by your personal exposure PRS.
This is the "so what does this actually mean" layer — translating
single-trait findings into projected health implications.
</p>
<div class="mr-caveat-block">
<strong>Limitations:</strong> Population-level causal estimates may not map
one-to-one to individuals. MR assumes the instruments affect the outcome
only through the exposure. Educational only — not for clinical decisions.
</div>
<div class="mr-list">{items}</div>
</section>
"""


def build_genetic_age_html(g: dict | None) -> str:
    if not g or not g.get("available"):
        return ""
    long_pct = g["longevity"]["percentile"]
    years = g.get("longevity_years_offset", 0)

    def _panel_meta(panel: dict) -> str:
        # Z-score + variant coverage as inner text. These panels are small (the
        # telomere proxy can be ≤2 variants), so coverage is essential honesty —
        # a percentile built on 1–2 SNPs is far weaker than the full panel.
        nu, nt = panel.get("n_used"), panel.get("n_total")
        z = panel.get("z")
        parts = []
        if z is not None:
            parts.append(f"Z {z:+.2f}")
        if nu is not None and nt:
            parts.append(f"{nu}/{nt} variants")
        return " · ".join(parts)

    def _meta_div(panel: dict) -> str:
        m = _panel_meta(panel)
        return f'<div class="ga-sub-meta">{m}</div>' if m else ""

    sub_blocks = ""
    if g.get("telomere"):
        sub_blocks += (
            f'<div class="ga-sub">'
            f'<div class="ga-sub-name">Telomere Length (genetic proxy)</div>'
            f'<div class="ga-sub-val">{g["telomere"]["percentile"]:.0f}th percentile</div>'
            f'{_meta_div(g["telomere"])}'
            f'</div>'
        )
    if g.get("skin_aging"):
        sub_blocks += (
            f'<div class="ga-sub">'
            f'<div class="ga-sub-name">Skin Aging</div>'
            f'<div class="ga-sub-val">{g["skin_aging"]["percentile"]:.0f}th percentile</div>'
            f'{_meta_div(g["skin_aging"])}'
            f'</div>'
        )
    return f"""
<section class="ga-section" id="genetic-age">
<h2>Genetic Longevity Profile <span class="pro-pill">V5</span></h2>
<div class="ga-headline">
  <div class="ga-headline-pct">{long_pct:.0f}<span class="ga-pct-th">th</span></div>
  <div class="ga-headline-meta">
    <div>percentile genetic longevity</div>
    <div class="ga-years">{years:+.1f} years vs European mean</div>
    {_meta_div(g["longevity"])}
  </div>
</div>
<div class="ga-narr">{_esc(g.get("narrative",""))}</div>
<div class="ga-subs">{sub_blocks}</div>
<div class="ga-disc">{_esc(g.get("disclaimer",""))}</div>
</section>
"""


def build_pgx_sim_html(s: dict | None) -> str:
    if not s or not s.get("available") or not s.get("drugs"):
        return ""
    items = ""
    for d in s["drugs"]:
        ae = d["combined_ae_rr"]
        df = d["combined_dose_factor"]
        ae_color = "var(--red)" if ae > 2 else ("var(--ora)" if ae > 1.3 else "var(--grn)")
        df_color = "var(--red)" if df == 0 else ("var(--ora)" if df < 0.8 else "var(--grn)")
        def _gene_metrics(g: dict) -> str:
            # Per-gene contribution to the combined estimate. For multi-gene
            # drugs (e.g. warfarin VKORC1×CYP2C9) the combined metrics are a
            # product of these; showing each gene's share makes the result
            # auditable. Values are computed by pgx_simulation, not shown until now.
            bits = []
            if g.get("clearance") is not None:
                bits.append(f'clearance {g["clearance"]}%')
            if g.get("dose_factor") is not None:
                bits.append(f'dose ×{g["dose_factor"]}')
            if g.get("ae_rr") is not None:
                bits.append(f'AE risk ×{g["ae_rr"]}')
            return (
                f' <span class="sim-gene-metrics">({" · ".join(bits)})</span>'
                if bits
                else ""
            )

        gene_lines = "".join(
            f'<div class="sim-gene-line"><strong>{_esc(g["gene"])}</strong> '
            f'<span class="pgx-pheno {("pheno-pm" if g["phenotype_code"]=="PM" else "pheno-um" if g["phenotype_code"]=="UM" else "pheno-im" if g["phenotype_code"]=="IM" else "pheno-rm" if g["phenotype_code"]=="RM" else "pheno-nm")}">{_esc(g["phenotype_code"])}</span>'
            f'{_gene_metrics(g)} '
            f'{_esc(g.get("note",""))}</div>'
            for g in d["gene_findings"]
        )
        ddi_lines = "".join(
            f'<li>{_esc(label)}: {_esc(note)}</li>'
            for label, note in d.get("ddi_notes", [])
        )
        items += f"""
<div class="sim-drug">
  <div class="sim-head">
    <span class="sim-name">{_esc(d["drug"])}</span>
    <span class="sim-gene-pill">primary: {_esc(d["primary_gene"])}</span>
  </div>
  <div class="sim-metrics">
    <div class="sim-metric">
      <div class="sim-metric-name">Relative clearance</div>
      <div class="sim-metric-val">{d["combined_clearance_pct"]}%</div>
    </div>
    <div class="sim-metric">
      <div class="sim-metric-name">Suggested dose factor</div>
      <div class="sim-metric-val" style="color:{df_color}">×{df}</div>
    </div>
    <div class="sim-metric">
      <div class="sim-metric-name">Adverse-event relative risk</div>
      <div class="sim-metric-val" style="color:{ae_color}">×{ae}</div>
    </div>
  </div>
  <div class="sim-genes">{gene_lines}</div>
  {f'<details class="sim-ddi"><summary>Drug-drug-gene interactions</summary><ul>{ddi_lines}</ul></details>' if ddi_lines else ''}
</div>
"""
    return f"""
<section class="sim-section" id="pgx-simulation">
<h2>Quantitative Drug-Response Simulation <span class="pro-pill">V5</span></h2>
<p class="sim-intro">
For each major drug, your specific gene-phenotype profile is translated into
estimated clearance rate, dose-adjustment factor, and side-effect relative
risk vs population average. <strong>These are illustrative simulations, not
prescriptions — actual dosing requires clinical assessment.</strong>
</p>
<div class="sim-list">{items}</div>
</section>
"""


def build_reproductive_html(r: dict | None) -> str:
    if not r or not r.get("scenarios"):
        return ""
    items = ""
    for s in r["scenarios"]:
        pop_rows = ""
        for pp in s.get("by_population", []):
            pop_rows += (
                f'<tr><td>{_esc(pp["population"])}</td>'
                f'<td>{pp["partner_carrier_freq"]*100:.2f}%</td>'
                f'<td><strong>{_esc(pp["p_affected_str"])}</strong></td>'
                f'<td>{_esc(pp.get("p_carrier_str",""))}</td></tr>'
            )
        dom_note = ""
        if s.get("dominant_note"):
            dom_note = f'<div class="rep-dom">{_esc(s["dominant_note"])}</div>'
        xnote = ""
        if s.get("xlinked_note"):
            xnote = f'<div class="rep-x">{_esc(s["xlinked_note"])}</div>'
        items += f"""
<div class="rep-scenario">
  <div class="rep-head">
    <span class="rep-gene">{_esc(s["gene"])}</span>
    <span class="rep-variant">{_esc(s["variant"])}</span>
    <span class="rep-inh">{_esc(s["inheritance"])}</span>
  </div>
  <div class="rep-disease">{_esc(s["disease"])}</div>
  {f'<div class="rep-pop-tbl"><table class="snp-tbl"><thead><tr><th>Partner population</th><th>Partner carrier freq</th><th>Affected child</th><th>Carrier child</th></tr></thead><tbody>{pop_rows}</tbody></table></div>' if pop_rows else ''}
  {dom_note}
  {xnote}
  <div class="rep-advice"><strong>Partner testing:</strong> {_esc(s.get("partner_testing_advice",""))}</div>
  <div class="rep-outlook"><strong>Treatment / outlook:</strong> {_esc(s.get("treatment_outlook",""))}</div>
</div>
"""
    roh_context = ""
    if r.get("roh_context"):
        roh_context = f'<div class="rep-roh-context">{_esc(r["roh_context"])}</div>'

    return f"""
<section class="rep-section" id="reproductive">
<h2>Reproductive Genetics Simulator <span class="pro-pill">V5</span></h2>
<p class="rep-intro">
For each carrier finding, modelled probability of an affected child given
different partner-ancestry carrier frequencies. Use this to prioritise
which conditions a prospective partner should be tested for.
</p>
{roh_context}
<div class="rep-list">{items}</div>
</section>
"""


def build_wellness_html(w: dict | None) -> str:
    if not w or not w.get("predictions"):
        return ""
    cat_order = ["Nutrition", "Sleep & Circadian", "Fitness & Recovery",
                 "Stress & Mood", "Aging & Longevity"]
    cards_by_cat = w.get("by_category", {})
    sections = ""
    for cat in cat_order:
        if cat not in cards_by_cat:
            continue
        items = ""
        for p in cards_by_cat[cat]:
            conf = p.get("confidence", "low")
            items += f"""
<div class="wellness-card">
  <div class="wellness-trait">{_esc(p.get("trait",""))}</div>
  <div class="wellness-result">{_esc(p.get("result",""))}</div>
  <div class="wellness-action"><strong>Action:</strong> {_esc(p.get("action",""))}</div>
  <div class="wellness-evidence">{_esc(p.get("evidence",""))} · <span class="conf-{conf}">conf: {conf}</span></div>
</div>
"""
        sections += f"""
<details class="wellness-cat" open>
  <summary><span class="wellness-cat-name">{_esc(cat)}</span>
    <span class="wellness-cat-count">{len(cards_by_cat[cat])} predictions</span></summary>
  <div class="wellness-grid">{items}</div>
</details>
"""
    return f"""
<section class="wellness-section" id="wellness">
<h2>Wellness Insights <span class="pro-pill">V4</span></h2>
<p class="wellness-intro">
Genotype-driven wellness signals across nutrition, sleep, fitness, stress, and
aging — focused on <strong>actionable lifestyle</strong> implications rather than
disease. {w['n_predictions']} predictions across {len(w['categories'])} categories.
</p>
{sections}
</section>
"""


def build_medications_html(med: dict | None) -> str:
    if not med or not med.get("reviews"):
        return ""
    cards = ""
    for r in med["reviews"]:
        if r["status"] == "unknown_drug":
            cards += f"""
<div class="med-card med-unknown">
  <div class="med-head">
    <span class="med-drug">{_esc(r["input"])}</span>
    <span class="med-note">Not in PGx catalogue</span>
  </div>
  <div class="med-msg">{_esc(r["message"])}</div>
</div>
"""
            continue
        findings_html = ""
        for f in r["findings"]:
            phen = _esc(f.get("phenotype") or "—")
            phen_cls = f.get("phenotype_class", "pheno-na")
            rec = _esc(f.get("recommendation") or "No drug-specific recommendation; see gene's general profile.")
            findings_html += f"""
<div class="med-finding">
  <div class="med-finding-line">
    <span class="med-gene">{_esc(f["gene"])}</span>
    <span class="med-pheno {phen_cls}">{phen}</span>
    {f'<span class="med-as">activity {_esc(f["activity_score"])}</span>' if f.get("activity_score") is not None else ''}
  </div>
  <div class="med-pathway"><em>{_esc(f.get("pathway_note",""))}</em></div>
  <div class="med-rec">{rec}</div>
  {f'<div class="med-guide">{_esc(f.get("guideline",""))}</div>' if f.get("guideline") else ''}
</div>
"""
        cards += f"""
<div class="med-card">
  <div class="med-head">
    <span class="med-drug">{_esc(r["input"])}</span>
    {('<span class="med-generic">(' + _esc(r["generic"]) + ')</span>') if r["generic"] != r["input"].lower() else ''}
  </div>
  {findings_html}
</div>
"""
    return f"""
<section class="med-section" id="medication-review">
<h2>Medication Review <span class="pro-pill">V3</span></h2>
<p class="med-intro">
Cross-referenced {med['n_input']} medication{'s' if med['n_input'] != 1 else ''}
against your pharmacogenomic phenotypes. <strong>Educational only — discuss
with your prescriber before any medication change.</strong>
</p>
<div class="med-list">{cards}</div>
</section>
"""


def build_economics_html(economics_result: dict | None) -> str:
    """Render the health-economics section: findings ranked by ROI, a clinic
    ROI dashboard and a payer-impact brief. Returns "" when there is nothing
    to show (no module / no findings), matching every sibling renderer."""
    if not economics_result:
        return ""
    findings = economics_result.get("findings_with_economics") or []
    if not findings:
        return ""

    def _money(v):
        try:
            return f"${float(v):,.0f}"
        except (TypeError, ValueError):
            return "—"

    # ── Main summary: high-confidence findings only (per the report's
    # confidence-tiering convention). The detailed table below keeps every
    # finding with its confidence badge.
    high_conf = economics_result.get("high_confidence") or [
        f for f in findings if f.get("confidence") == "high"
    ]
    headline_html = ""
    if high_conf:
        top = sorted(high_conf, key=lambda f: (f.get("roi") or 0), reverse=True)[:3]
        items = "".join(
            f'<li><strong>{_esc(f.get("finding",""))}</strong> — '
            f'ROI <span class="econ-roi">{f.get("roi")}:1</span>, '
            f'value {_money(f.get("outcome_value"))} for {_money(f.get("intervention_cost"))} '
            f'({_esc(f.get("cost_basis",""))})</li>'
            for f in top
        )
        headline_html = (
            '<div class="econ-headline"><h3>Top high-confidence interventions</h3>'
            f'<ul class="econ-headline-list">{items}</ul></div>'
        )

    # ── Findings table (already ROI-ranked by the module) ──
    rows = ""
    for f in findings:
        conf = f.get("confidence", "n/a")
        roi = f.get("roi")
        roi_txt = f"{roi}:1" if roi is not None else "—"
        payback = f.get("payback_months")
        payback_txt = f"{payback} mo" if payback is not None else "—"
        # QALY gain + cost-per-QALY are computed per finding but were dropped in
        # favour of the payer aggregate. Cost/QALY is the standard cost-
        # effectiveness metric (≈ <$50k/QALY is conventionally "cost-effective").
        qaly = f.get("qaly_gain")
        qaly_txt = f"{qaly:.3f}".rstrip("0").rstrip(".") if isinstance(qaly, (int, float)) else "—"
        cpq = f.get("cost_per_qaly")
        cpq_sub = (
            f'<div class="econ-sub">{_money(cpq)}/QALY</div>'
            if isinstance(cpq, (int, float))
            else ""
        )
        rows += (
            f'<tr class="conf-{_esc(conf)}">'
            f'<td><strong>{_esc(f.get("finding",""))}</strong>'
            f'<div class="econ-benefit">{_esc(f.get("clinical_benefit",""))}</div></td>'
            f'<td>{_money(f.get("intervention_cost"))}'
            f'<div class="econ-sub">{_esc(f.get("cost_basis",""))}</div></td>'
            f'<td>{_money(f.get("outcome_value"))}</td>'
            f'<td>{qaly_txt}{cpq_sub}</td>'
            f'<td class="econ-roi"><strong>{roi_txt}</strong></td>'
            f'<td>{payback_txt}</td>'
            f'<td>{_money(f.get("npv_3year"))}</td>'
            f'<td><span class="econ-conf econ-conf-{_esc(conf)}">{_esc(conf)}</span></td>'
            f'</tr>\n'
        )

    # ── Clinic dashboard ──
    clinic = economics_result.get("clinic_dashboard") or {}
    clinic_html = ""
    if clinic.get("n_findings"):
        clinic_html = f"""
<div class="econ-panel">
  <h3>Clinic ROI Dashboard <span class="econ-pop">{clinic.get('patient_count')} patients</span></h3>
  <div class="econ-cards">
    <div class="econ-card"><div class="econ-card-v">{_money(clinic.get('avg_cost_per_patient'))}</div><div class="econ-card-l">Avg cost / patient</div></div>
    <div class="econ-card"><div class="econ-card-v">{_money(clinic.get('avg_benefit_per_patient'))}</div><div class="econ-card-l">Avg benefit / patient</div></div>
    <div class="econ-card"><div class="econ-card-v">{clinic.get('avg_roi')}:1</div><div class="econ-card-l">Average ROI</div></div>
    <div class="econ-card"><div class="econ-card-v">{clinic.get('payback_period_months')} mo</div><div class="econ-card-l">Payback period</div></div>
  </div>
  <p class="econ-model">Subscription model: {_money(clinic.get('revenue_model_monthly'))}/patient/month
  at {int(float(clinic.get('gross_margin', 0)) * 100)}% gross margin.
  Across {clinic.get('patient_count')} patients — cost {_money(clinic.get('total_cost'))},
  modeled benefit {_money(clinic.get('total_benefit'))}.</p>
</div>"""

    # ── Payer impact ──
    payer = economics_result.get("payer_impact") or {}
    payer_html = ""
    if payer.get("affected_members"):
        cpq = payer.get("cost_per_qaly")
        cpq_txt = _money(cpq) if cpq is not None else "n/a"
        payer_html = f"""
<div class="econ-panel">
  <h3>Payer Impact <span class="econ-pop">{payer.get('member_population'):,} members</span></h3>
  <div class="econ-cards">
    <div class="econ-card"><div class="econ-card-v">{payer.get('affected_members'):,}</div><div class="econ-card-l">Members affected</div></div>
    <div class="econ-card"><div class="econ-card-v">{_money(payer.get('total_cost'))}</div><div class="econ-card-l">Total intervention cost</div></div>
    <div class="econ-card"><div class="econ-card-v">{_money(payer.get('total_benefit'))}</div><div class="econ-card-l">Modeled savings</div></div>
    <div class="econ-card"><div class="econ-card-v">{payer.get('roi')}:1</div><div class="econ-card-l">Aggregate ROI</div></div>
    <div class="econ-card"><div class="econ-card-v">{cpq_txt}</div><div class="econ-card-l">Cost per QALY</div></div>
    <div class="econ-card"><div class="econ-card-v">{_money(payer.get('net_savings'))}</div><div class="econ-card-l">Net savings</div></div>
  </div>
</div>"""

    disclaimer = _esc(economics_result.get("disclaimer", ""))

    return f"""
<section class="econ-section" id="health-economics">
<h2>Health Economics <span class="pro-pill">ROI</span></h2>
<p class="econ-intro">
Clinical and payer return-on-investment for acting on your genomic findings.
Each intervention's cost is weighed against the modeled value of the adverse
outcome it averts, with payback period and 3-year NPV (discounted at 3%).
</p>
{headline_html}
<div class="tbl-wrap">
<table class="econ-tbl"><thead><tr>
<th>Finding &amp; intervention</th><th>Cost</th><th>Outcome value</th>
<th>QALY gain</th><th>ROI</th><th>Payback</th><th>NPV (3y)</th><th>Confidence</th>
</tr></thead><tbody>{rows}</tbody></table>
</div>
{clinic_html}
{payer_html}
<p class="econ-disclaimer"><strong>Estimates only.</strong> {disclaimer}</p>
</section>
"""


def build_pharmgkb_html(pk: dict | None) -> str:
    """ClinPGx/PharmGKB clinical-annotation positions typed on this chip.

    High-evidence (Level 1A/1B/2A/2B) positions are shown open; the weak,
    unreplicated Level 3/4 long tail is behind a collapsed, explicitly-labelled
    disclosure. Every row is framed as an *annotated position* (with the user's
    genotype), NOT a direction-of-effect call — the source table has no risk
    allele. See module docstring for the accuracy contract.
    """
    if not pk or not pk.get("available"):
        return ""

    def _level_badge(lvl: str) -> str:
        cls = "pk-lvl-high" if lvl in ("1A", "1B", "2A", "2B") else "pk-lvl-low"
        return f'<span class="pk-lvl {cls}">{_esc(lvl)}</span>'

    def _rows(entries: list[dict]) -> str:
        out = ""
        for e in entries:
            pgx_note = (
                ' <span class="pk-pgx" title="Also analysed at allele level in '
                'the Pharmacogenomics section">↗ PGx</span>'
                if e.get("pgx_gene")
                else ""
            )
            # One row per annotation (drug group) so every association is listed.
            for a in e["annotations"]:
                drugs = ", ".join(_esc(d) for d in a.get("drugs", [])) or "—"
                phen = ", ".join(_esc(p) for p in a.get("phenotypes", [])) or "—"
                out += (
                    f"<tr>"
                    f'<td class="rsid-cell">'
                    f'<a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(e["rsid"])}" '
                    f'target="_blank" rel="noopener">{_esc(e["rsid"])}</a></td>'
                    f'<td><strong>{_esc(e["gene"])}</strong>{pgx_note}</td>'
                    f'<td class="gt-cell">{_esc(e["genotype"])}</td>'
                    f"<td>{_level_badge(a.get('level',''))}</td>"
                    f"<td>{drugs}</td>"
                    f'<td>{_esc(a.get("type",""))}</td>'
                    f"<td>{phen}</td>"
                    f"</tr>\n"
                )
        return out

    def _table(entries: list[dict]) -> str:
        return (
            f'<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>'
            f"<th>rsID</th><th>Gene</th><th>Your genotype</th><th>Evidence</th>"
            f"<th>Drug(s)</th><th>Association</th><th>Phenotype</th>"
            f"</tr></thead><tbody>{_rows(entries)}</tbody></table></div>"
        )

    high_html = _table(pk["high"]) if pk.get("high") else ""
    low_html = ""
    if pk.get("low"):
        low_html = (
            f'<details class="pk-lowtail"><summary>Weak / unreplicated '
            f"annotations — ClinPGx Level 3 & 4 ({pk['n_low']} positions). "
            f"Single-study, case-report, in-vitro, or non-significant; shown "
            f"for completeness, not for action.</summary>{_table(pk['low'])}</details>"
        )

    return f"""
<section class="pk-section" id="pharmgkb-clinical">
<h2>Pharmacogenomic Variant Annotations <span class="pro-pill">ClinPGx</span></h2>
<p class="pk-intro">
Every variant on this chip cross-referenced against the
<strong>ClinPGx / PharmGKB clinical-variant annotation</strong> database
({pk.get('n_dataset_rsid_rows',0):,} annotated rsIDs). You carry a genotype at
<strong>{pk['n_typed_variants']}</strong> annotated positions
({pk['n_high']} high-evidence, {pk['n_low']} weak) spanning
{pk['n_drugs']} drugs, across {pk['n_annotation_rows']} annotation records.
</p>
<div class="pk-caveat">
<strong>Read this as "positions worth knowing about," not a phenotype call.</strong>
Evidence levels are the <a href="https://www.clinpgx.org/page/clinAnnLevels"
target="_blank" rel="noopener">ClinPGx/PharmGKB clinical-annotation levels</a>
(1A/1B strong & replicated → 3/4 weak/unreplicated) — <em>not</em> CPIC
guideline strength. The source table carries no risk allele or direction of
effect, so this panel flags that you have a genotype at an annotated position;
it does not tell you which way it pushes. For the clinically-actionable genes
(marked <span class="pk-pgx">↗ PGx</span>), the dedicated <a
href="#pharmacogenomics">Pharmacogenomics</a> section — with star-allele calls
and CPIC dosing — is authoritative.
</div>
<h3 class="pk-h3">High-evidence annotations (ClinPGx Level 1A / 1B / 2A / 2B)</h3>
{high_html or '<p class="pk-none">No high-evidence annotated positions typed on this chip.</p>'}
{low_html}
</section>
"""


def build_top_drugs_html(td: dict | None) -> str:
    """Top-prescribed-medications pharmacogenomic screen.

    Four tiers: drugs your genotype may affect (open); PGx drugs where your
    relevant gene is typed and normal; PGx drugs where the gene is unresolved
    on this chip; and drugs with no known PGx interaction — the last three
    collapsed. Gene phenotypes come from pgx.py; drug↔gene links and evidence
    levels from CPIC + ClinPGx/PharmGKB. Nothing here is a dosing instruction —
    the Pharmacogenomics section is authoritative for direction.
    """
    if not td or not td.get("available"):
        return ""

    def _pheno_cell(e: dict) -> str:
        gp = e.get("gene_phenotypes") or []
        if not gp:
            return "—"
        return "<br>".join(
            f'{_esc(p["gene"])}: <strong>{_esc(p.get("phenotype") or p.get("code") or "")}</strong>'
            for p in gp
        )

    def _lvl(e: dict) -> str:
        cp = e.get("cpic_level")
        cl = e.get("clin_level")
        bits = []
        if cp:
            bits.append(f'<span class="td-lvl td-cpic">CPIC {_esc(cp)}</span>')
        if cl:
            bits.append(f'<span class="td-lvl td-clin">PharmGKB {_esc(cl)}</span>')
        return " ".join(bits) or "—"

    def _table(rows: list[dict], with_pheno: bool) -> str:
        body = ""
        for e in rows:
            genes = ", ".join(_esc(g) for g in e.get("genes", [])) or "—"
            pheno = f"<td>{_pheno_cell(e)}</td>" if with_pheno else ""
            body += (
                f"<tr><td><strong>{_esc(e['generic'])}</strong>"
                f'<div class="td-brand">{_esc(e.get("brand",""))}</div></td>'
                f"<td>{_esc(e.get('class',''))}</td>"
                f"<td>{genes}</td>"
                f"{pheno}"
                f"<td>{_lvl(e)}</td></tr>\n"
            )
        pheno_th = "<th>Your metabolizer status</th>" if with_pheno else ""
        return (
            f'<div class="tbl-wrap"><table class="snp-tbl"><thead><tr>'
            f"<th>Drug</th><th>Class</th><th>PGx gene(s)</th>{pheno_th}"
            f"<th>Evidence</th></tr></thead><tbody>{body}</tbody></table></div>"
        )

    actionable_html = (
        _table(td["actionable"], True)
        if td.get("actionable")
        else '<p class="td-none">No top-prescribed drug maps to a gene where your '
             "genotype is a non-normal metabolizer (given the genes typed on this chip).</p>"
    )
    tn = (
        f'<details class="td-fold"><summary>Standard dosing likely — relevant gene '
        f'typed & normal ({td["n_typed_normal"]})</summary>{_table(td["typed_normal"], True)}</details>'
        if td.get("typed_normal") else ""
    )
    pr = (
        f'<details class="td-fold"><summary>Pharmacogenomically relevant, but the '
        f'gene was not resolved on this chip ({td["n_pgx_relevant"]}) — imputation '
        f'(--impute) would resolve many</summary>{_table(td["pgx_relevant"], True)}</details>'
        if td.get("pgx_relevant") else ""
    )
    npx = (
        f'<details class="td-fold"><summary>No known pharmacogenomic interaction '
        f'in the bundled CPIC/PharmGKB databases ({td["n_no_pgx"]})</summary>'
        f'{_table(td["no_pgx"], False)}</details>'
        if td.get("no_pgx") else ""
    )

    return f"""
<section class="td-section" id="top-drugs">
<h2>Top Prescribed Drugs — Pharmacogenomic Screen <span class="pro-pill">PGx</span></h2>
<p class="td-intro">
Each of the <strong>{td['n_screened']}</strong> most commonly prescribed
medications, screened against your genome using the bundled CPIC/DPWG drug
database and the ClinPGx/PharmGKB drug table.
<strong>{td['n_with_pgx']}</strong> have documented pharmacogenomic data;
<strong>{td['n_actionable']}</strong> map to a gene where your genotype is a
non-normal metabolizer and may warrant a dosing discussion.
</p>
<div class="td-caveat">
<strong>Not a prescription or dosing instruction.</strong> Drug↔gene links and
evidence levels are from CPIC and ClinPGx/PharmGKB; your metabolizer status is
from the star-allele calls in the <a href="#pharmacogenomics">Pharmacogenomics</a>
section, which is authoritative for direction and dosing. Drugs shown as "no
known interaction" simply have no entry in the bundled databases — absence of a
gene-drug pair here is not proof of safety. Always consult your prescriber or
pharmacist.
</div>
<h3 class="td-h3">Review with your prescriber — your genotype may change dosing ({td['n_actionable']})</h3>
{actionable_html}
{tn}
{pr}
{npx}
</section>
"""


# ── HTML report builder ───────────────────────────────────────────────────────

def build_html_report(
    tier1_results: list[dict],
    apoe_genotype: str | None,
    ai_results: dict[str, str],
    exec_summary: str | None,
    dna_filepath: str,
    no_ai: bool,
    model: str,
    y_result: dict | None = None,
    mt_result: dict | None = None,
    cross_cat_synthesis: str | None = None,
    qc_result: dict | None = None,
    prs_result: dict | None = None,
    pgx_result: dict | None = None,
    interactions_result: dict | None = None,
    carrier_result: dict | None = None,
    traits_result: dict | None = None,
    counseling_result: dict | None = None,
    references_used: list[dict] | None = None,
    imputation_info: dict | None = None,
    expanded_pgs_result: dict | None = None,
    ancestry_result: dict | None = None,
    medications_result: dict | None = None,
    wellness_result: dict | None = None,
    hla_result: dict | None = None,
    roh_result: dict | None = None,
    local_ancestry_result: dict | None = None,
    phewas_result: dict | None = None,
    mr_result: dict | None = None,
    genetic_age_result: dict | None = None,
    pgx_sim_result: dict | None = None,
    reproductive_result: dict | None = None,
    economics_result: dict | None = None,
    pharmgkb_result: dict | None = None,
    top_drugs_result: dict | None = None,
    metal_oxidative_result: dict | None = None,
    detox_result: dict | None = None,
    urologic_result: dict | None = None,
    deep_ancestry_result: dict | None = None,
    blood_type_result: dict | None = None,
    immunogenetics_result: dict | None = None,
    ancestral_story_result: dict | None = None,
    neurochemistry_result: dict | None = None,
    holistic_synthesis_result: dict | None = None,
    addiction_genetics_result: dict | None = None,
    family_planning_result: dict | None = None,
    polygenic_traits_result: dict | None = None,
    tnrc18_result: dict | None = None,
    environmental_optimization_result: dict | None = None,
    life_stage_playbook_result: dict | None = None,
    clinical_variants_result: dict | None = None,
    novel_variants_result: dict | None = None,
    voi_result: dict | None = None,
    module_ai: dict | None = None,
) -> str:
    # build_category_map expands cross-referenced SNPs into each relevant
    # category so they render in multiple sections with the appropriate
    # implication text.
    categories_map = build_category_map(tier1_results)

    sorted_cats = [
        (c, categories_map[c]) for c in CATEGORY_ORDER if c in categories_map
    ]
    for c, s in categories_map.items():
        if c not in CATEGORY_ORDER:
            sorted_cats.append((c, s))

    report_date = datetime.datetime.now().strftime("%B %d, %Y at %H:%M")
    dna_filename = Path(dna_filepath).name

    apoe_info = (
        APOE_INFO.get(
            apoe_genotype,
            {
                "risk_label": "Unknown",
                "color_class": "apoe-gray",
                "description": "Could not determine APOE genotype from available data.",
            },
        )
        if apoe_genotype
        else None
    )

    total_matched = len(tier1_results)
    high_risk_list = sorted(
        [r for r in tier1_results if r["risk_copies"] > 0 and r["significance"] == "high"],
        key=lambda x: x["category"],
    )
    moderate_risk_list = sorted(
        [r for r in tier1_results if r["risk_copies"] > 0 and r["significance"] == "moderate"],
        key=lambda x: x["category"],
    )
    risk_variants = len([r for r in tier1_results if r["risk_copies"] > 0])
    high_count = len(high_risk_list)

    # ── Category sections ──
    cat_sections_html = ""
    for cat, snps in sorted_cats:
        cat_id = _cat_id(cat)

        rows = ""
        for s in sorted(snps, key=lambda x: (-x["risk_copies"], x["gene"])):
            xref_badge = ""
            row_classes = []
            if s["risk_copies"] > 0:
                row_classes.append("risk-row")
            if s.get("is_cross_ref"):
                row_classes.append("xref-row")
                xref_badge = (
                    f' <span class="xref-badge" title="Cross-referenced from '
                    f'{s.get("primary_category", "")}">↗ {s.get("primary_category", "")}</span>'
                )
            row_cls = " ".join(row_classes)

            # Surface the verified per-variant context that the curated
            # database already carries but the summary column alone hid:
            # actionable recommendation, chip-coverage caveats, and the
            # cross-system implications. All author-curated text — nothing
            # is synthesised here.
            detail_bits = []
            rec = s.get("recommendation")
            if rec:
                detail_bits.append(
                    f'<div class="rec-line"><span class="rec-label">'
                    f'What to do</span> {_esc(rec)}</div>'
                )
            cov = s.get("chip_coverage_note")
            if cov:
                detail_bits.append(
                    f'<div class="rec-line rec-cov"><span class="rec-label">'
                    f'Chip coverage</span> {_esc(cov)}</div>'
                )
            xrefs = s.get("cross_references") or []
            if xrefs and not s.get("is_cross_ref"):
                xref_items = "".join(
                    f'<li><strong>{_esc(x.get("category", ""))}:</strong> '
                    f'{_esc(x.get("implication", ""))}</li>'
                    for x in xrefs
                    if isinstance(x, dict) and x.get("implication")
                )
                if xref_items:
                    detail_bits.append(
                        f'<div class="rec-line"><span class="rec-label">'
                        f'Related systems</span><ul class="rec-xref">'
                        f"{xref_items}</ul></div>"
                    )
            detail_html = ""
            if detail_bits:
                detail_html = (
                    '<details class="var-detail"><summary>'
                    "Recommendation &amp; context</summary>"
                    f'<div class="var-detail-body">{"".join(detail_bits)}</div>'
                    "</details>"
                )

            # Imputed calls are statistical, not measured — flag them (with R²
            # when available) so an expanded --impute run stays honest.
            imp_badge = ""
            if s.get("source") == "imputed":
                r2 = s.get("r2")
                r2txt = f" r²={r2}" if r2 is not None else ""
                imp_badge = (
                    f' <span class="imp-badge" title="Statistically imputed '
                    f'(1000G), not directly typed on the chip">imputed{r2txt}</span>'
                )
            rows += (
                f'<tr class="{row_cls}">'
                f'<td class="rsid-cell">'
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{s["rsid"]}" '
                f'target="_blank" rel="noopener">{s["rsid"]}</a></td>'
                f'<td><strong>{s["gene"]}</strong></td>'
                f'<td>{s["variant_name"]}{xref_badge}</td>'
                f'<td class="gt-cell">{s["my_genotype"]}{imp_badge}</td>'
                f'<td>{risk_indicator(s["risk_copies"])}</td>'
                f'<td>{significance_badge(s["significance"])}</td>'
                f'<td class="sum-cell">{s["summary"]}{detail_html}</td>'
                f"</tr>\n"
            )

        ai_block = ""
        if not no_ai and cat in ai_results:
            ai_html = md_to_html(ai_results[cat])
            ai_block = (
                f'<div class="ai-section">'
                f'<div class="ai-title">AI Interpretation</div>'
                f'<div class="ai-content">{ai_html}</div>'
                f"</div>"
            )

        # V4: count risk-carrying variants for the section header
        n_total = len(snps)
        n_risk = sum(1 for s in snps if s["risk_copies"] > 0)
        n_xref = sum(1 for s in snps if s.get("is_cross_ref"))
        cat_sections_html += (
            f'<details class="cat-section" id="{cat_id}" open>'
            f'<summary class="cat-summary">'
            f'<span class="cat-h2-inline">{cat}</span> '
            f'<span class="cat-count">{n_total} variant{"s" if n_total != 1 else ""}'
            f'{f" · {n_risk} with risk allele" if n_risk else ""}'
            f'{f" · {n_xref} cross-ref" if n_xref else ""}</span>'
            f'</summary>'
            f'<div class="tbl-wrap">'
            f'<table class="snp-tbl"><thead><tr>'
            f"<th>rsID</th><th>Gene</th><th>Variant</th>"
            f"<th>Genotype</th><th>Risk Alleles</th>"
            f"<th>Significance</th><th>Summary</th>"
            f"</tr></thead><tbody>{rows}</tbody></table>"
            f"</div>{ai_block}</details>\n"
        )

    # ── APOE section ──
    if apoe_genotype and apoe_info:
        apoe_html = (
            f'<section class="apoe-section {apoe_info["color_class"]}" id="apoe">'
            f'<h2>APOE Genotype: <span class="apoe-gt">{apoe_genotype}</span> '
            f'<span class="apoe-rl">{apoe_info["risk_label"]}</span></h2>'
            f'<p class="apoe-desc">{apoe_info["description"]}</p>'
            f'<p class="apoe-note">Determined by combining rs429358 (E4 marker) '
            f"and rs7412 (E2 marker). APOE genotype is the strongest common "
            f"genetic predictor of late-onset Alzheimer's disease and also "
            f"influences cardiovascular and cholesterol metabolism.</p>"
            f"</section>"
        )
    else:
        apoe_html = (
            '<section class="apoe-section apoe-gray" id="apoe">'
            "<h2>APOE Genotype: <span class=\"apoe-gt\">Not Determined</span></h2>"
            "<p>rs429358 and/or rs7412 were not called in your DNA file, "
            "so the APOE genotype cannot be determined. Your DNA provider may "
            "not include these SNPs, or they may have had insufficient coverage.</p>"
            "</section>"
        )

    # ── Executive summary ──
    exec_html = ""
    if not no_ai and exec_summary:
        exec_html = (
            '<section class="exec-section" id="exec-summary">'
            '<h2>Executive Summary <span class="ai-pill">AI</span></h2>'
            f'<div class="exec-content">{md_to_html(exec_summary)}</div>'
            "</section>"
        )

    # ── Recommendations section ──
    rec_items = ""
    for r in high_risk_list[:12]:
        rec_items += (
            f"<li><strong>{r['gene']} ({r['variant_name']})</strong>: "
            f"{r['recommendation']}</li>\n"
        )
    for r in moderate_risk_list[:10]:
        rec_items += (
            f"<li>{r['gene']} ({r['variant_name']}): "
            f"{r['recommendation']}</li>\n"
        )

    rec_section = (
        '<section class="rec-section" id="recommendations">'
        "<h2>Actionable Recommendations Summary</h2>"
        "<p>Based on Tier 1 deterministic lookup. High-significance findings first.</p>"
        f'<ul class="rec-list">{rec_items}</ul>'
        "</section>"
    )

    # ── TOC nav ──
    nav_links = ""
    for cat, snps in sorted_cats:
        cat_id = (
            cat.lower()
            .replace(" ", "-")
            .replace("/", "-")
            .replace("'", "")
            .replace("&", "")
            .replace(".", "")
        )
        # V4: variant counts inline with nav links
        nav_links += (
            f'<a href="#{cat_id}" class="nl">'
            f'{cat} <span class="nl-count">{len(snps)}</span>'
            f'</a>\n'
        )

    ai_tier_note = (
        f"Tier 1 + Tier 2 AI ({model})"
        if not no_ai
        else "Tier 1 only (--no-ai)"
    )
    qc_chip = ""
    if qc_result:
        qc_chip = (
            f' &nbsp;&middot;&nbsp; '
            f'<span class="hdr-qc {qc_result["grade_class"]}">QC: {qc_result["grade"]}</span>'
        )

    # ── Y-DNA section ──
    ydna_section_html = build_ydna_html(y_result) if y_result else ""

    # ── mtDNA section ──
    mtdna_section_html = build_mtdna_html(mt_result) if mt_result else ""

    # ── Cross-Category Interactions section ──
    cross_cat_html = ""
    if not no_ai and cross_cat_synthesis:
        cross_cat_html = (
            '<section class="cross-cat-section" id="cross-category">'
            '<h2>Cross-Category Interactions <span class="ai-pill">AI</span></h2>'
            '<p class="cross-cat-intro">How findings from different categories '
            'compound or interact in ways that wouldn\'t be obvious from any '
            'single-category view.</p>'
            f'<div class="cross-cat-content">{md_to_html(cross_cat_synthesis)}</div>'
            "</section>"
        )

    # ── Professional-grade sections ──
    qc_html = build_qc_html(qc_result)
    prs_html = build_prs_html(prs_result)
    pgx_html = build_pgx_html(pgx_result)
    interactions_html = build_interactions_html(interactions_result)
    carrier_html = build_carrier_html(carrier_result)
    traits_html = build_traits_html(traits_result)
    counseling_html = build_counseling_html(counseling_result)
    references_html = build_references_html(references_used)
    # ── V3 sections ──
    imputation_html = build_imputation_html(imputation_info)
    expanded_pgs_html = build_expanded_pgs_html(expanded_pgs_result)
    ancestry_html = build_ancestry_html(ancestry_result)
    deep_ancestry_html = build_deep_ancestry_html(deep_ancestry_result)
    blood_type_html = build_blood_type_html(blood_type_result)
    immunogenetics_html = build_immunogenetics_html(immunogenetics_result)
    neurochemistry_html = build_neurochemistry_html(neurochemistry_result)
    holistic_synthesis_html = build_holistic_synthesis_html(holistic_synthesis_result)
    addiction_genetics_html = build_addiction_genetics_html(addiction_genetics_result)
    family_planning_html = build_family_planning_html(family_planning_result)
    clinical_variants_html = build_clinical_variants_html(clinical_variants_result)
    novel_variants_html = build_novel_variants_html(novel_variants_result)
    voi_html = build_voi_html(voi_result)
    polygenic_traits_html = build_polygenic_traits_html(polygenic_traits_result)
    try:
        from tnrc18_marker import build_tnrc18_html as _build_tnrc18
        tnrc18_html = _build_tnrc18(tnrc18_result)
    except Exception:
        tnrc18_html = ""
    environmental_optimization_html = build_environmental_optimization_html(environmental_optimization_result)
    life_stage_playbook_html = build_life_stage_playbook_html(life_stage_playbook_result)
    ancestral_story_html = build_ancestral_story_html(ancestral_story_result)

    # AI on all tiers — attach each module's local-AI interpretation under its heading.
    _mai = module_ai or {}
    deep_ancestry_html = _prepend_ai_interpretation(deep_ancestry_html, _mai.get("deep-ancestry"))
    blood_type_html = _prepend_ai_interpretation(blood_type_html, _mai.get("blood-type"))
    immunogenetics_html = _prepend_ai_interpretation(immunogenetics_html, _mai.get("immunogenetics"))
    neurochemistry_html = _prepend_ai_interpretation(neurochemistry_html, _mai.get("neurochemistry"))
    holistic_synthesis_html = _prepend_ai_interpretation(holistic_synthesis_html, _mai.get("holistic-synthesis"))
    clinical_variants_html = _prepend_ai_interpretation(clinical_variants_html, _mai.get("clinical-variants"))
    novel_variants_html = _prepend_ai_interpretation(novel_variants_html, _mai.get("novel-variants"))
    voi_html = _prepend_ai_interpretation(voi_html, _mai.get("value-of-information"))
    addiction_genetics_html = _prepend_ai_interpretation(addiction_genetics_html, _mai.get("addiction-genetics"))
    family_planning_html = _prepend_ai_interpretation(family_planning_html, _mai.get("family-planning"))
    polygenic_traits_html = _prepend_ai_interpretation(polygenic_traits_html, _mai.get("polygenic-traits"))
    environmental_optimization_html = _prepend_ai_interpretation(environmental_optimization_html, _mai.get("environmental-optimization"))
    medications_html = build_medications_html(medications_result)
    # ── V4 sections ──
    wellness_html = build_wellness_html(wellness_result)
    # ── V5 sections ──
    hla_html = build_hla_html(hla_result)
    roh_html = build_roh_html(roh_result)
    local_ancestry_html = build_local_ancestry_html(local_ancestry_result)
    phewas_html = build_phewas_html(phewas_result)
    mr_html = build_mr_html(mr_result)
    genetic_age_html = build_genetic_age_html(genetic_age_result)
    pgx_sim_html = build_pgx_sim_html(pgx_sim_result)
    reproductive_html = build_reproductive_html(reproductive_result)
    economics_html = build_economics_html(economics_result)
    pharmgkb_html = build_pharmgkb_html(pharmgkb_result)
    top_drugs_html = build_top_drugs_html(top_drugs_result)
    # ── V9/V10 sections ──
    detox_html = build_detox_html(detox_result)
    urologic_html = build_urologic_html(urologic_result)
    metal_oxidative_html = build_metal_oxidative_html(metal_oxidative_result)

    # ── Full HTML ──
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>DNA Analysis Report — {report_date}</title>
<style>
/* ==== RESET & ROOT ==== */
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
:root{{
  --bg:#0d1117;--bg2:#161b22;--bg3:#21262d;--bg4:#30363d;
  --bdr:#30363d;--txt:#e6edf3;--muted:#8b949e;--dim:#656d76;
  --acc:#58a6ff;--acc2:#79c0ff;--grn:#3fb950;--ora:#d29922;
  --red:#f85149;--yel:#e3b341;--pur:#bc8cff;--blu:#58a6ff;
  --font:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  --mono:"SF Mono","Fira Code",Consolas,monospace;
  --r:8px;
}}
@media(prefers-color-scheme:light){{
  :root{{
    --bg:#ffffff;--bg2:#f6f8fa;--bg3:#eaeef2;--bg4:#d0d7de;
    --bdr:#d0d7de;--txt:#1f2328;--muted:#57606a;--dim:#6e7781;
    --acc:#0969da;--acc2:#0550ae;--grn:#1a7f37;--ora:#9a6700;
    --red:#cf222e;--yel:#9a6700;--pur:#8250df;--blu:#0550ae;
  }}
}}
html{{scroll-behavior:smooth}}
body{{font-family:var(--font);background:var(--bg);color:var(--txt);
  line-height:1.65;font-size:15px}}
a{{color:var(--acc);text-decoration:none}}
a:hover{{color:var(--acc2);text-decoration:underline}}
code{{font-family:var(--mono);font-size:0.85em;background:var(--bg3);
  padding:1px 5px;border-radius:4px}}

/* ==== LAYOUT ==== */
.wrap{{max-width:1120px;margin:0 auto;padding:0 24px 80px}}

/* ==== STICKY HEADER ==== */
.hdr{{background:var(--bg2);border-bottom:1px solid var(--bdr);
  padding:18px 24px;position:sticky;top:0;z-index:100;
  box-shadow:0 2px 10px rgba(0,0,0,.35)}}
.hdr-inner{{max-width:1120px;margin:0 auto;
  display:flex;justify-content:space-between;align-items:center;
  gap:12px;flex-wrap:wrap}}
.hdr-title{{font-size:19px;font-weight:700;letter-spacing:-.3px}}
.hdr-meta{{font-size:12px;color:var(--muted)}}

/* ==== DISCLAIMER ==== */
.disc{{background:rgba(88,166,255,.07);border:1px solid var(--bdr);
  border-left:4px solid var(--acc);border-radius:var(--r);
  padding:14px 18px;margin:28px 0;font-size:13px;color:var(--muted)}}
.disc strong{{color:var(--txt)}}

/* ==== STATS BAR ==== */
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px;margin-bottom:28px}}
.stat{{background:var(--bg2);border:1px solid var(--bdr);border-radius:var(--r);
  padding:16px;text-align:center}}
.stat-n{{font-size:34px;font-weight:800;color:var(--acc)}}
.stat-l{{font-size:11px;color:var(--muted);text-transform:uppercase;
  letter-spacing:.6px;margin-top:4px}}

/* ==== TOC ==== */
.toc{{background:var(--bg2);border:1px solid var(--bdr);border-radius:var(--r);
  padding:18px;margin-bottom:32px}}
.toc h3{{font-size:11px;text-transform:uppercase;letter-spacing:.8px;
  color:var(--muted);margin-bottom:10px}}
.toc-links{{display:flex;flex-wrap:wrap;gap:7px}}
.nl{{font-size:12px;background:var(--bg3);border:1px solid var(--bdr);
  border-radius:20px;padding:3px 11px;color:var(--txt);transition:background .18s}}
.nl:hover{{background:var(--bg4);text-decoration:none}}

/* ==== SECTIONS ==== */
section{{margin-bottom:48px}}
h2{{font-size:21px;font-weight:700;margin-bottom:14px;padding-bottom:8px;
  border-bottom:1px solid var(--bdr);letter-spacing:-.3px}}

/* ==== APOE ==== */
.apoe-section{{border-radius:var(--r);padding:22px 26px;margin-bottom:36px;
  border:1px solid var(--bdr)}}
.apoe-section h2{{border-bottom-color:transparent;display:flex;
  align-items:center;gap:10px;flex-wrap:wrap}}
.apoe-gt{{font-size:27px;font-weight:800}}
.apoe-rl{{font-size:13px;font-weight:600;padding:3px 11px;border-radius:20px;
  background:rgba(255,255,255,.1)}}
.apoe-desc{{font-size:15px;line-height:1.75;margin-bottom:10px}}
.apoe-note{{font-size:12px;color:var(--muted);font-style:italic}}
.apoe-green{{background:linear-gradient(135deg,rgba(63,185,80,.12),rgba(63,185,80,.04));
  border-color:rgba(63,185,80,.35)}}
.apoe-blue{{background:linear-gradient(135deg,rgba(88,166,255,.12),rgba(88,166,255,.04));
  border-color:rgba(88,166,255,.35)}}
.apoe-yellow{{background:linear-gradient(135deg,rgba(227,179,65,.12),rgba(227,179,65,.04));
  border-color:rgba(227,179,65,.35)}}
.apoe-orange{{background:linear-gradient(135deg,rgba(248,153,0,.13),rgba(248,153,0,.04));
  border-color:rgba(248,153,0,.4)}}
.apoe-red{{background:linear-gradient(135deg,rgba(248,81,73,.16),rgba(248,81,73,.05));
  border-color:rgba(248,81,73,.45)}}
.apoe-gray{{background:var(--bg2);border-color:var(--bdr)}}

/* ==== EXEC SUMMARY ==== */
.exec-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:22px 26px;margin-bottom:36px}}
.ai-pill{{font-size:10px;background:rgba(188,140,255,.2);color:var(--pur);
  padding:2px 7px;border-radius:10px;font-weight:600;vertical-align:middle;
  margin-left:6px}}
.exec-content p,.ai-content p{{margin-bottom:10px;line-height:1.75;font-size:14px}}
.exec-content ul,.ai-content ul{{padding-left:20px;margin-bottom:10px}}
.exec-content li,.ai-content li{{margin-bottom:5px;font-size:14px;line-height:1.6}}
.exec-content h4,.ai-content h4{{margin:14px 0 6px;font-size:14px;
  color:var(--acc);border:none;padding:0}}
.exec-content strong,.ai-content strong{{color:var(--txt)}}

/* ==== CATEGORY ==== */
.cat-h2{{font-size:19px}}

/* ==== TABLE ==== */
.tbl-wrap{{overflow-x:auto;border-radius:var(--r);border:1px solid var(--bdr);
  margin-bottom:18px}}
.snp-tbl{{width:100%;border-collapse:collapse;font-size:13px}}
.snp-tbl th{{background:var(--bg3);color:var(--muted);font-weight:600;
  font-size:10.5px;text-transform:uppercase;letter-spacing:.5px;
  padding:9px 13px;text-align:left;white-space:nowrap;
  border-bottom:1px solid var(--bdr)}}
.snp-tbl td{{padding:9px 13px;border-bottom:1px solid var(--bdr);vertical-align:top}}
.snp-tbl tr:last-child td{{border-bottom:none}}
.snp-tbl tr:hover{{background:var(--bg2)}}
.risk-row{{background:rgba(248,81,73,.04)}}
.risk-row:hover{{background:rgba(248,81,73,.09)!important}}
.rsid-cell{{font-family:var(--mono);font-size:12px;white-space:nowrap}}
.gt-cell{{font-family:var(--mono);font-weight:700;font-size:14px}}
.sum-cell{{font-size:12px;color:var(--muted);max-width:280px;line-height:1.5}}
.var-detail{{margin-top:6px}}
.var-detail>summary{{cursor:pointer;font-size:11px;font-weight:600;color:var(--accent,#2563eb);list-style:none;user-select:none}}
.var-detail>summary::-webkit-details-marker{{display:none}}
.var-detail>summary::before{{content:"\\25B8 ";font-size:10px}}
.var-detail[open]>summary::before{{content:"\\25BE "}}
.var-detail-body{{margin-top:6px;padding:8px 10px;background:var(--card2,rgba(127,127,127,.06));border-left:2px solid var(--accent,#2563eb);border-radius:4px}}
.rec-line{{font-size:12px;line-height:1.55;color:var(--txt);margin:4px 0}}
.rec-line .rec-label{{display:inline-block;font-weight:700;font-size:10px;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);margin-right:6px}}
.rec-cov{{color:var(--warn,#b45309)}}
.rec-xref{{margin:4px 0 0 0;padding-left:16px}}
.rec-xref li{{margin:3px 0}}

/* ==== BADGES ==== */
.badge{{display:inline-block;padding:2px 8px;border-radius:10px;
  font-size:11px;font-weight:600}}
.badge-high{{background:rgba(248,81,73,.2);color:var(--red)}}
.badge-moderate{{background:rgba(210,153,34,.2);color:var(--ora)}}
.badge-low{{background:rgba(63,185,80,.2);color:var(--grn)}}

/* ==== RISK INDICATORS ==== */
.risk-none{{color:var(--grn);font-size:12px;white-space:nowrap}}
.risk-one{{color:var(--ora);font-size:12px;white-space:nowrap;font-weight:600}}
.risk-two{{color:var(--red);font-size:12px;white-space:nowrap;font-weight:700}}

/* ==== AI BLOCK ==== */
.ai-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-left:3px solid var(--pur);border-radius:var(--r);
  padding:18px 22px;margin-top:2px}}
.ai-title{{font-size:12px;font-weight:700;color:var(--pur);
  text-transform:uppercase;letter-spacing:.6px;margin-bottom:12px}}

/* ==== RECOMMENDATIONS ==== */
.rec-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:22px 26px}}
.rec-list{{padding-left:20px}}
.rec-list li{{margin-bottom:9px;line-height:1.65;font-size:14px}}
.rec-list li strong{{color:var(--acc)}}

/* ==== METHOD ==== */
.meth-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:22px 26px}}
.meth-section p{{margin-bottom:10px;font-size:14px;color:var(--muted);
  line-height:1.7}}
.meth-section strong{{color:var(--txt)}}

/* ==== FOOTER ==== */
.ftr{{text-align:center;padding:28px;color:var(--dim);font-size:12px;
  border-top:1px solid var(--bdr);margin-top:40px}}
.ftr-cite{{margin-top:12px;color:var(--muted);font-size:12px;line-height:1.5}}
.ftr-copy{{margin-top:4px;color:var(--dim);font-size:11px}}

/* ==== PRINT ==== */
@media print{{
  .hdr{{position:static;box-shadow:none}}
  .toc,.nl{{display:none}}
  body{{background:#fff;color:#000;font-size:11pt}}
  .apoe-section,.exec-section,.ai-section,.rec-section,.meth-section{{
    border:1px solid #ccc!important;background:#f9f9f9!important}}
  h2{{page-break-after:avoid}}
  .cat-section{{page-break-inside:avoid}}
  a{{color:#000;text-decoration:underline}}
  .snp-tbl{{font-size:9pt}}
}}

/* ==== RESPONSIVE ==== */
@media(max-width:768px){{
  .wrap{{padding:0 14px 60px}}
  .hdr-inner{{flex-direction:column;align-items:flex-start}}
  .stat-n{{font-size:26px}}
  .apoe-gt{{font-size:22px}}
}}

/* ==== Y-DNA HAPLOGROUP ==== */
.ydna-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:22px 26px;margin-bottom:36px}}
.ydna-badge{{font-size:12px;font-weight:600;padding:3px 10px;
  border-radius:12px;vertical-align:middle;margin-left:8px}}
.ydna-badge-green{{background:rgba(63,185,80,.2);color:var(--grn)}}
.ydna-badge-amber{{background:rgba(210,153,34,.2);color:var(--ora)}}
.ydna-badge-blue{{background:rgba(88,166,255,.2);color:var(--acc)}}
.ydna-badge-gray{{background:var(--bg3);color:var(--muted)}}

/* breadcrumb path */
.ydna-path-wrap{{background:var(--bg);border:1px solid var(--bdr);
  border-radius:var(--r);padding:16px 18px;margin-bottom:16px}}
.ydna-path-label{{font-size:10px;text-transform:uppercase;letter-spacing:.7px;
  color:var(--muted);margin-bottom:10px}}
.ydna-crumbs{{display:flex;flex-wrap:wrap;align-items:center;gap:0;
  font-family:var(--mono);font-size:13px;margin-bottom:10px}}
.crumb-arrow{{color:var(--muted);margin:0 6px;font-size:16px;line-height:1}}
.crumb-ok{{display:inline-flex;flex-direction:column;align-items:center;
  background:rgba(63,185,80,.12);border:1px solid rgba(63,185,80,.35);
  color:var(--grn);font-weight:700;border-radius:6px;padding:4px 10px;
  min-width:52px;text-align:center;cursor:default}}
.crumb-gap{{display:inline-flex;flex-direction:column;align-items:center;
  background:rgba(210,153,34,.1);border:1px dashed rgba(210,153,34,.45);
  color:var(--ora);font-weight:600;border-radius:6px;padding:4px 10px;
  min-width:52px;text-align:center;cursor:default}}
.crumb-snp{{font-size:9px;font-weight:400;color:var(--muted);margin-top:2px}}
.ydna-path-note{{font-size:11px;color:var(--muted);margin-top:6px}}

/* info blocks */
.ydna-gap{{background:rgba(210,153,34,.08);border:1px solid rgba(210,153,34,.3);
  border-radius:var(--r);padding:12px 16px;margin:12px 0;font-size:13px}}
.ydna-branches{{background:var(--bg3);border:1px solid var(--bdr);
  border-radius:var(--r);padding:14px 16px;margin:12px 0;font-size:13px}}
.ydna-branches ul{{padding-left:18px;margin-top:6px}}
.ydna-branches li{{margin-bottom:6px;line-height:1.55}}
.ydna-migration{{background:linear-gradient(135deg,rgba(88,166,255,.06),
  rgba(188,140,255,.06));border:1px solid rgba(88,166,255,.25);
  border-radius:var(--r);padding:16px 20px;margin:14px 0}}
.ydna-migration-title{{font-size:11px;text-transform:uppercase;letter-spacing:.7px;
  color:var(--acc);font-weight:700;margin-bottom:8px}}
.ydna-migration p{{font-size:14px;line-height:1.75;color:var(--txt)}}
.ydna-further{{background:rgba(188,140,255,.07);border:1px solid rgba(188,140,255,.25);
  border-radius:var(--r);padding:14px 18px;margin:12px 0;font-size:13px;
  color:var(--txt)}}
.ydna-further p{{margin-top:6px;line-height:1.6;white-space:pre-line;font-size:13px}}
.ydna-ancient{{background:linear-gradient(135deg,rgba(63,185,80,.06),rgba(88,166,255,.05));
  border:1px solid rgba(63,185,80,.25);border-radius:var(--r);
  padding:14px 18px;margin:14px 0}}
.ydna-ancient-title{{font-size:11px;text-transform:uppercase;letter-spacing:.7px;
  color:var(--grn);font-weight:700;margin-bottom:8px}}
.ydna-ancient p{{font-size:13px;line-height:1.7;color:var(--txt)}}

/* ==== mtDNA (mirrors Y-DNA styling) ==== */
.mtdna-section{{background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:22px 26px;margin-bottom:36px}}
.mtdna-badge{{font-size:12px;font-weight:600;padding:3px 10px;
  border-radius:12px;vertical-align:middle;margin-left:8px}}
.mtdna-badge-pink{{background:rgba(248,81,249,.18);color:#e879f9}}
.mtdna-badge-gray{{background:var(--bg3);color:var(--muted)}}
.mtdna-result{{background:var(--bg);border:1px solid var(--bdr);
  border-radius:var(--r);padding:16px 18px;margin-bottom:16px}}
.mtdna-call{{font-family:var(--mono);font-size:22px;font-weight:800;
  color:#e879f9;margin-bottom:6px}}
.mtdna-evidence{{font-size:12px;color:var(--muted);margin-top:6px}}
.mtdna-migration{{background:linear-gradient(135deg,rgba(248,81,249,.07),
  rgba(188,140,255,.06));border:1px solid rgba(248,81,249,.25);
  border-radius:var(--r);padding:16px 20px;margin:14px 0}}
.mtdna-migration-title{{font-size:11px;text-transform:uppercase;letter-spacing:.7px;
  color:#e879f9;font-weight:700;margin-bottom:8px}}
.mtdna-migration p{{font-size:14px;line-height:1.75;color:var(--txt);margin-bottom:8px}}

/* ==== Cross-Category Interactions ==== */
.cross-cat-section{{background:linear-gradient(135deg,rgba(188,140,255,.07),
  rgba(88,166,255,.05));border:1px solid rgba(188,140,255,.3);
  border-radius:var(--r);padding:22px 26px;margin-bottom:36px}}
.cross-cat-intro{{font-size:13px;color:var(--muted);margin-bottom:14px;
  font-style:italic}}
.cross-cat-content p{{margin-bottom:10px;line-height:1.75;font-size:14px}}
.cross-cat-content h4{{margin:18px 0 8px;font-size:14px;color:var(--pur);
  border:none;padding:0;font-weight:700}}
.cross-cat-content strong{{color:var(--txt)}}
.cross-cat-content ul{{padding-left:20px;margin-bottom:10px}}
.cross-cat-content li{{margin-bottom:5px;font-size:14px;line-height:1.6}}

/* ==== Cross-reference badge on SNP rows ==== */
.xref-badge{{display:inline-block;font-size:10px;background:rgba(188,140,255,.18);
  color:var(--pur);padding:1px 7px;border-radius:10px;margin-left:6px;
  font-weight:600;cursor:default;vertical-align:middle}}
.xref-row{{background:rgba(188,140,255,.04)}}
.xref-row:hover{{background:rgba(188,140,255,.09)!important}}

/* ==== Pro pill / version chip ==== */
.pro-pill{{font-size:10px;background:linear-gradient(135deg,#3b82f6,#8b5cf6);
  color:#fff;padding:2px 8px;border-radius:10px;font-weight:700;
  vertical-align:middle;margin-left:8px;letter-spacing:.5px}}
.hdr-version{{font-size:10px;color:var(--muted);font-weight:500;margin-left:6px}}
.hdr-qc{{font-size:11px;padding:2px 8px;border-radius:10px;font-weight:600}}
.nl-pro{{background:linear-gradient(135deg,rgba(59,130,246,.18),rgba(139,92,246,.18))!important;
  border-color:rgba(139,92,246,.4)!important;color:var(--txt)!important}}

/* ==== QC SECTION ==== */
.qc-section{{margin-bottom:36px}}
.qc-grade-card{{display:flex;gap:24px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:20px 24px;margin-bottom:20px;align-items:flex-start}}
.qc-grade{{font-size:32px;font-weight:800;padding:14px 22px;border-radius:var(--r);
  white-space:nowrap;letter-spacing:-.5px}}
.qc-grade-excellent{{background:rgba(63,185,80,.16);color:var(--grn);border:1px solid rgba(63,185,80,.4)}}
.qc-grade-good{{background:rgba(88,166,255,.16);color:var(--acc);border:1px solid rgba(88,166,255,.4)}}
.qc-grade-fair{{background:rgba(210,153,34,.16);color:var(--ora);border:1px solid rgba(210,153,34,.4)}}
.qc-grade-limited{{background:rgba(248,81,73,.16);color:var(--red);border:1px solid rgba(248,81,73,.4)}}
.qc-grade-text{{flex:1}}
.qc-grade-note{{font-size:14px;margin-bottom:10px;line-height:1.55}}
.qc-grade-stats{{font-size:12px;color:var(--muted);margin-bottom:4px;line-height:1.55}}
.qc-grade-stats strong{{color:var(--txt)}}
.qc-h3{{font-size:14px;font-weight:700;color:var(--muted);text-transform:uppercase;
  letter-spacing:.6px;margin:18px 0 10px}}
.qc-tbl{{width:100%;border-collapse:collapse;font-size:12px}}
.qc-tbl th{{background:var(--bg3);font-size:10.5px;text-transform:uppercase;
  letter-spacing:.5px;padding:8px 12px;text-align:left;color:var(--muted)}}
.qc-tbl td{{padding:8px 12px;border-bottom:1px solid var(--bdr);vertical-align:middle}}
.qc-domain-name{{font-weight:600;white-space:nowrap}}
.qc-bar-wrap{{display:inline-block;width:140px;height:8px;background:var(--bg3);
  border-radius:4px;overflow:hidden;vertical-align:middle;margin-right:8px}}
.qc-bar{{height:100%;border-radius:4px}}
.qc-bar-full{{background:var(--grn)}}
.qc-bar-good{{background:var(--acc)}}
.qc-bar-fair{{background:var(--ora)}}
.qc-bar-poor{{background:var(--red)}}
.qc-bar-label{{font-size:11px;color:var(--muted)}}
.qc-chips{{display:flex;flex-wrap:wrap;gap:4px}}
.qc-var-on,.qc-var-off{{display:inline-block;font-size:10px;padding:2px 7px;
  border-radius:8px;font-family:var(--mono);font-weight:600}}
.qc-var-on{{background:rgba(63,185,80,.16);color:var(--grn)}}
.qc-var-off{{background:rgba(248,81,73,.12);color:var(--red);text-decoration:line-through;
  text-decoration-color:rgba(248,81,73,.5)}}

/* ==== PRS SECTION ==== */
.prs-section{{margin-bottom:42px}}
.prs-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:12px}}
.prs-caveat{{background:rgba(210,153,34,.07);border-left:3px solid var(--ora);
  border-radius:6px;padding:12px 16px;font-size:12.5px;color:var(--muted);
  line-height:1.65;margin-bottom:18px}}
.prs-caveat strong{{color:var(--ora)}}
.prs-headline{{background:rgba(248,81,73,.08);border:1px solid rgba(248,81,73,.3);
  border-radius:var(--r);padding:14px 18px;margin-bottom:18px;font-size:14px}}
.prs-headline strong{{color:var(--red)}}
.prs-headline ul{{padding-left:22px;margin-top:6px}}
.prs-headline li{{margin-bottom:4px;font-size:13px}}
.prs-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:14px}}
.prs-card{{background:var(--bg2);border:1px solid var(--bdr);border-radius:var(--r);
  padding:18px 20px;display:flex;flex-direction:column;gap:10px}}
.prs-card-na{{opacity:.6}}
.prs-card-head{{display:flex;justify-content:space-between;align-items:flex-start;
  gap:12px;flex-wrap:wrap}}
.prs-name{{font-size:15px;font-weight:700;line-height:1.3}}
.prs-short{{font-size:10px;background:var(--bg3);padding:2px 7px;border-radius:8px;
  color:var(--muted);font-weight:600;margin-left:6px;vertical-align:middle}}
.prs-tier{{font-size:11px;font-weight:700;padding:4px 10px;border-radius:10px;
  white-space:nowrap;text-transform:uppercase;letter-spacing:.5px}}
.tier-high{{background:rgba(248,81,73,.2);color:var(--red)}}
.tier-elevated{{background:rgba(210,153,34,.2);color:var(--ora)}}
.tier-average{{background:rgba(88,166,255,.18);color:var(--acc)}}
.tier-below{{background:rgba(63,185,80,.18);color:var(--grn)}}
.tier-low{{background:rgba(63,185,80,.22);color:var(--grn)}}
.prs-desc{{font-size:12.5px;color:var(--muted);line-height:1.6}}
.prs-na{{font-size:12px;color:var(--muted);font-style:italic;margin-top:8px}}
.prs-track{{position:relative;height:10px;background:linear-gradient(90deg,
  rgba(63,185,80,.4) 0%,rgba(88,166,255,.4) 20%,rgba(88,166,255,.4) 80%,
  rgba(210,153,34,.4) 95%,rgba(248,81,73,.4) 100%);border-radius:5px;margin:6px 0}}
.prs-track-pointer{{position:absolute;top:-3px;width:3px;height:16px;
  background:var(--txt);transform:translateX(-50%);border-radius:1px}}
.prs-details{{display:flex;gap:14px;font-size:12px;flex-wrap:wrap}}
.prs-lab{{color:var(--muted)}}
.prs-context{{font-size:12px;color:var(--muted);line-height:1.6;background:var(--bg);
  border:1px solid var(--bdr);border-radius:6px;padding:10px 12px}}
.prs-context strong{{color:var(--txt)}}
.prs-ref{{font-size:11px;color:var(--dim);margin-top:4px}}
.prs-variants{{margin:8px 0}}
.prs-variants>summary{{cursor:pointer;font-size:11px;font-weight:600;color:var(--accent,#2563eb)}}
.prs-variants table{{margin-top:6px;font-size:11px}}
.prs-missing{{font-size:11px;color:var(--warn,#b45309);margin-top:6px;line-height:1.5}}

/* ==== PGX SECTION ==== */
.pgx-section{{margin-bottom:42px}}
.pgx-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:12px}}
.pgx-caveat{{background:rgba(210,153,34,.07);border-left:3px solid var(--ora);
  border-radius:6px;padding:12px 16px;font-size:12.5px;color:var(--muted);
  line-height:1.65;margin-bottom:18px}}
.pgx-actionable{{background:rgba(248,81,73,.06);border:1px solid rgba(248,81,73,.25);
  border-radius:var(--r);padding:16px 20px;margin-bottom:22px}}
.pgx-actionable h3{{font-size:13px;text-transform:uppercase;letter-spacing:.6px;
  color:var(--red);margin-bottom:10px;font-weight:700}}
.actionable-list{{list-style:none;padding:0}}
.actionable-item{{padding:10px 0;border-bottom:1px solid var(--bdr)}}
.actionable-item:last-child{{border-bottom:none}}
.actionable-head{{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:5px}}
.actionable-gene{{font-weight:700}}
.actionable-drug{{color:var(--acc);font-weight:600}}
.actionable-rec{{font-size:13px;line-height:1.55;color:var(--txt)}}
.pgx-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:14px}}
.pgx-card{{background:var(--bg2);border:1px solid var(--bdr);border-radius:var(--r);
  padding:18px 22px}}
.pgx-card-head{{margin-bottom:12px}}
.pgx-gene-name{{font-size:18px;font-weight:800;letter-spacing:-.3px}}
.pgx-gene-long{{font-size:11px;color:var(--muted);text-transform:uppercase;
  letter-spacing:.5px;margin-top:2px}}
.pgx-score{{display:flex;align-items:center;gap:10px;margin-bottom:8px;flex-wrap:wrap}}
.pgx-activity-label{{font-size:11px;color:var(--muted);text-transform:uppercase;
  letter-spacing:.4px}}
.pgx-activity{{font-size:18px;font-weight:800;font-family:var(--mono);color:var(--acc)}}
.pgx-activity-base{{font-size:11px;color:var(--dim)}}
.pgx-pheno{{display:inline-block;font-size:12px;font-weight:700;padding:4px 11px;
  border-radius:14px;letter-spacing:.3px}}
.pheno-pm{{background:rgba(248,81,73,.18);color:var(--red)}}
.pheno-im{{background:rgba(210,153,34,.18);color:var(--ora)}}
.pheno-nm{{background:rgba(63,185,80,.18);color:var(--grn)}}
.pheno-rm{{background:rgba(88,166,255,.18);color:var(--acc)}}
.pheno-um{{background:rgba(188,140,255,.18);color:var(--pur)}}
.pheno-indeterminate{{background:var(--bg3);color:var(--muted);border:1px dashed var(--bdr)}}
.pgx-callability{{font-size:11px;color:var(--muted);margin-bottom:10px}}
/* Shared confidence badge used by every score section */
.conf-row{{display:flex;align-items:baseline;gap:8px;flex-wrap:wrap;margin:8px 0}}
.conf-badge{{font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:.4px;
  padding:2px 8px;border-radius:10px;white-space:nowrap}}
.conf-high{{background:rgba(63,185,80,.18);color:var(--grn)}}
.conf-moderate{{background:rgba(210,153,34,.18);color:var(--ora)}}
.conf-low{{background:rgba(248,81,73,.16);color:var(--red)}}
.conf-none{{background:var(--bg3);color:var(--muted);border:1px dashed var(--bdr)}}
.conf-na{{background:var(--bg3);color:var(--muted)}}
.conf-note{{font-size:11px;color:var(--muted);line-height:1.45}}
.prs-lowconf,.anc-ambiguous{{font-size:11.5px;color:var(--ora);background:rgba(210,153,34,.08);
  border-left:3px solid var(--ora);padding:6px 10px;border-radius:4px;margin:6px 0;line-height:1.45}}
.anc-bars-label{{font-size:11px;color:var(--muted);margin:10px 0 4px}}
.ga-sub-meta{{font-size:10.5px;color:var(--muted);margin-top:2px;font-variant-numeric:tabular-nums}}
.hla-tag-missing{{opacity:.7;text-decoration:line-through}}
.ref-clinvar{{font-size:11px;color:var(--muted);font-weight:600}}
.la-ll{{font-size:10.5px;color:var(--muted);margin-top:2px;font-variant-numeric:tabular-nums}}
.pk-intro{{font-size:12.5px;color:var(--txt);line-height:1.6;margin:6px 0}}
.pk-caveat{{font-size:11.5px;color:var(--muted);line-height:1.55;background:var(--bg2);border:1px solid var(--bdr);border-left:3px solid var(--accent,#2563eb);padding:8px 12px;border-radius:5px;margin:8px 0}}
.pk-caveat strong{{color:var(--txt)}}
.pk-h3{{font-size:13px;text-transform:uppercase;letter-spacing:.5px;margin:16px 0 8px;font-weight:700}}
.pk-lvl{{display:inline-block;padding:1px 7px;border-radius:9px;font-size:10.5px;font-weight:700}}
.pk-lvl-high{{background:rgba(34,153,84,.15);color:var(--grn,#2e8b57)}}
.pk-lvl-low{{background:rgba(127,127,127,.13);color:var(--muted)}}
.pk-pgx{{font-size:10px;color:var(--accent,#2563eb);font-weight:600}}
.pk-lowtail{{margin-top:12px}}
.pk-lowtail>summary{{cursor:pointer;font-size:11.5px;font-weight:600;color:var(--warn,#b45309);line-height:1.5}}
.pk-none{{font-size:12px;color:var(--muted)}}
.imp-badge{{display:inline-block;font-size:9.5px;font-weight:700;padding:1px 5px;border-radius:8px;background:rgba(147,112,219,.18);color:#7c5cbf;letter-spacing:.02em;vertical-align:middle}}
.td-intro{{font-size:12.5px;color:var(--txt);line-height:1.6;margin:6px 0}}
.td-caveat{{font-size:11.5px;color:var(--muted);line-height:1.55;background:var(--bg2);border:1px solid var(--bdr);border-left:3px solid var(--accent,#2563eb);padding:8px 12px;border-radius:5px;margin:8px 0}}
.td-caveat strong{{color:var(--txt)}}
.td-h3{{font-size:13px;text-transform:uppercase;letter-spacing:.5px;margin:16px 0 8px;font-weight:700;color:var(--red,#c0392b)}}
.td-brand{{font-size:10px;color:var(--muted)}}
.td-lvl{{display:inline-block;padding:1px 6px;border-radius:8px;font-size:10px;font-weight:700;margin-right:3px}}
.td-cpic{{background:rgba(34,153,84,.16);color:var(--grn,#2e8b57)}}
.td-clin{{background:rgba(127,127,127,.13);color:var(--muted)}}
.td-fold{{margin-top:12px}}
.td-fold>summary{{cursor:pointer;font-size:12px;font-weight:600;color:var(--accent,#2563eb);line-height:1.5}}
.td-none{{font-size:12px;color:var(--muted)}}
.anc-margin{{font-size:11px;color:var(--muted);margin:6px 0;line-height:1.45}}
.anc-margin em{{color:var(--txt);font-style:normal;font-weight:600}}
.anc-aims{{margin-top:8px}}
.anc-aims>summary{{cursor:pointer;font-size:11px;font-weight:600;color:var(--accent,#2563eb)}}
.anc-aims table{{margin-top:6px;font-size:11px}}
.anc-aim-redundant{{opacity:.6}}
.ydna-disclaimer{{font-size:11px;color:var(--muted);font-style:italic;margin-top:14px;
  padding-top:10px;border-top:1px solid var(--bdr);line-height:1.5}}
.pgx-details summary{{cursor:pointer;font-size:12px;color:var(--muted);padding:4px 0;
  user-select:none}}
.pgx-details summary:hover{{color:var(--txt)}}
.pgx-tbl{{width:100%;border-collapse:collapse;font-size:11.5px}}
.pgx-tbl th{{background:var(--bg3);font-size:10px;text-transform:uppercase;
  letter-spacing:.4px;padding:6px 10px;text-align:left;color:var(--muted)}}
.pgx-tbl td{{padding:6px 10px;border-bottom:1px solid var(--bdr)}}
.pgx-uncalled{{opacity:.5}}
.pgx-impact{{font-family:var(--mono);font-weight:700}}
.pgx-na{{color:var(--dim);font-style:italic}}
.pgx-drugs{{margin-top:14px}}
.pgx-drugs-title{{font-size:11px;text-transform:uppercase;letter-spacing:.5px;
  color:var(--muted);font-weight:700;margin-bottom:8px}}
.pgx-drug-row{{background:var(--bg);border:1px solid var(--bdr);border-radius:6px;
  padding:10px 12px;margin-bottom:8px}}
.pgx-drug-name{{font-weight:700;font-size:13px;margin-bottom:6px;color:var(--acc)}}
.pgx-rec-line{{font-size:12px;margin-bottom:3px;line-height:1.55}}
.pgx-rec-pheno{{display:inline-block;font-size:10px;font-weight:700;padding:1px 6px;
  border-radius:6px;background:var(--bg3);color:var(--muted);margin-right:6px;
  font-family:var(--mono);min-width:32px;text-align:center}}
.pgx-rec-text{{color:var(--muted)}}
.pgx-rec-active{{background:rgba(88,166,255,.08);padding:6px 8px;border-radius:5px;
  border-left:2px solid var(--acc)}}
.pgx-rec-active .pgx-rec-pheno{{background:var(--acc);color:#fff}}
.pgx-rec-active .pgx-rec-text{{color:var(--txt);font-weight:500}}
.pgx-caveat{{font-size:11px;color:var(--muted);font-style:italic;margin-top:10px;
  padding:8px 10px;background:var(--bg3);border-radius:5px;line-height:1.6}}
.pgx-guideline{{font-size:11px;color:var(--dim);margin-top:10px}}

/* ==== INTERACTIONS ==== */
.inter-section{{margin-bottom:42px}}
.inter-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:8px}}
.inter-stats{{font-size:12px;color:var(--muted);margin-bottom:14px}}
.inter-list{{display:flex;flex-direction:column;gap:12px}}
.inter-finding{{background:var(--bg2);border:1px solid var(--bdr);border-left:4px solid var(--ora);
  border-radius:var(--r);padding:14px 18px}}
.inter-head{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:8px}}
.inter-sev{{font-size:10px;font-weight:800;letter-spacing:.6px;padding:3px 9px;
  border-radius:8px}}
.sev-high{{background:var(--red);color:#fff}}
.sev-mod{{background:var(--ora);color:#fff}}
.sev-low{{background:var(--acc);color:#fff}}
.inter-title{{font-size:15px;font-weight:700}}
.inter-vars{{margin-bottom:10px;display:flex;gap:6px;flex-wrap:wrap}}
.inter-var{{font-size:11px;background:var(--bg3);padding:2px 8px;border-radius:8px}}
.inter-body{{font-size:13px;line-height:1.65}}
.inter-interp{{margin-bottom:8px;color:var(--txt)}}
.inter-action{{color:var(--txt)}}
.inter-interp strong,.inter-action strong{{color:var(--acc)}}

/* ==== CARRIER ==== */
.carr-section{{margin-bottom:42px}}
.carr-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.carr-h3{{font-size:13px;text-transform:uppercase;letter-spacing:.6px;
  margin:18px 0 10px;font-weight:700}}
.carr-aff{{color:var(--red)}}
.carr-car{{color:var(--ora)}}
.carr-nc,.carr-ut{{font-size:12px;color:var(--muted);line-height:1.6;margin-top:14px;
  padding:10px 14px;background:var(--bg2);border:1px solid var(--bdr);border-radius:6px}}
.carr-nc strong,.carr-ut strong{{color:var(--txt)}}
.carr-freq{{font-size:11px;color:var(--muted);margin-top:4px}}
.carr-freq em{{color:var(--txt);font-style:normal;font-weight:600}}
.carr-caveat{{font-size:11px;color:var(--warn,#b45309);margin-top:4px;line-height:1.5}}

/* ==== TRAITS ==== */
.traits-section{{margin-bottom:42px}}
.traits-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.traits-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:10px}}
.trait-card{{background:var(--bg2);border:1px solid var(--bdr);border-radius:var(--r);
  padding:14px 16px;display:flex;flex-direction:column;gap:5px}}
.trait-na{{opacity:.55}}
.trait-name{{font-size:11px;text-transform:uppercase;letter-spacing:.6px;
  color:var(--muted);font-weight:700}}
.trait-result{{font-size:14px;font-weight:700;line-height:1.4}}
.trait-evidence{{font-size:11px;color:var(--dim);font-family:var(--mono)}}
.trait-detail{{font-size:12px;color:var(--muted);line-height:1.6;
  background:var(--bg);border-left:2px solid var(--acc);padding:6px 10px;
  border-radius:4px;margin-top:4px}}
.trait-conf{{font-size:10px;text-transform:uppercase;letter-spacing:.6px;
  margin-top:4px}}
.conf-high{{color:var(--grn)}}
.conf-moderate{{color:var(--acc)}}
.conf-low{{color:var(--ora)}}
.conf-n\\/a{{color:var(--dim)}}
.conf-error{{color:var(--red)}}

/* ==== COUNSELING ==== */
.couns-section{{margin-bottom:42px;background:linear-gradient(135deg,
  rgba(248,81,73,.05),rgba(210,153,34,.04));border:1px solid rgba(248,81,73,.2);
  border-radius:var(--r);padding:22px 26px}}
.couns-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.couns-list{{display:flex;flex-direction:column;gap:12px}}
.couns-finding{{background:var(--bg2);border:1px solid var(--bdr);border-left:3px solid var(--red);
  border-radius:6px;padding:12px 16px}}
.couns-trigger{{font-size:14px;font-weight:700;margin-bottom:6px;color:var(--red)}}
.couns-urgency,.couns-specialist{{font-size:12px;color:var(--muted);margin-bottom:3px}}
.couns-urgency strong,.couns-specialist strong{{color:var(--txt)}}
.couns-reason{{font-size:12.5px;color:var(--muted);line-height:1.6;margin-top:6px}}

/* ==== REFERENCES ==== */
.ref-section{{margin-bottom:42px}}
.ref-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.ref-tbl{{width:100%;border-collapse:collapse;font-size:12px}}
.ref-tbl th{{background:var(--bg3);font-size:10px;text-transform:uppercase;
  letter-spacing:.4px;padding:8px 12px;text-align:left;color:var(--muted)}}
.ref-tbl td{{padding:8px 12px;border-bottom:1px solid var(--bdr);vertical-align:top}}
.ref-cat-row td{{background:var(--bg2);font-size:11px;text-transform:uppercase;
  letter-spacing:.5px;color:var(--muted);padding:8px 12px!important}}
.ref-level{{display:inline-block;font-size:10px;font-weight:700;padding:2px 7px;
  border-radius:8px;white-space:nowrap}}
.lvl-cpic-a{{background:rgba(63,185,80,.2);color:var(--grn)}}
.lvl-cpic-b{{background:rgba(63,185,80,.14);color:var(--grn)}}
.lvl-cpic-c{{background:rgba(63,185,80,.1);color:var(--grn)}}
.lvl-clin{{background:rgba(88,166,255,.2);color:var(--acc)}}
.lvl-clin-mod{{background:rgba(88,166,255,.14);color:var(--acc)}}
.lvl-gwas-a{{background:rgba(188,140,255,.18);color:var(--pur)}}
.lvl-pop{{background:rgba(210,153,34,.15);color:var(--ora)}}
.lvl-other{{background:var(--bg3);color:var(--muted)}}
.ref-pmids{{font-size:11px;line-height:1.7}}
.ref-pmids a{{margin-right:6px;font-family:var(--mono)}}
.ref-no-pmid{{color:var(--dim);font-style:italic}}

/* ==== V3 sections ==== */
.nl-v3{{background:linear-gradient(135deg,rgba(16,185,129,.18),rgba(59,130,246,.18))!important;
  border-color:rgba(16,185,129,.4)!important;color:var(--txt)!important}}

/* ==== IMPUTATION ==== */
.impute-section{{margin-bottom:36px;background:var(--bg2);border:1px solid var(--bdr);
  border-left:4px solid #10b981;border-radius:var(--r);padding:20px 24px}}
.impute-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.impute-stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));
  gap:12px;margin:14px 0}}
.impute-stat{{background:var(--bg);border:1px solid var(--bdr);border-radius:6px;
  padding:12px;text-align:center}}
.impute-stat-imp{{border-color:rgba(16,185,129,.4)}}
.impute-n{{font-size:24px;font-weight:800;color:#10b981}}
.impute-l{{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.5px;margin-top:4px}}
.impute-warn{{background:rgba(210,153,34,.07);border-left:3px solid var(--ora);
  border-radius:5px;padding:10px 14px;margin-top:10px;font-size:12.5px;color:var(--muted)}}

/* ==== EXPANDED PGS (reuses .prs-* classes) ==== */
.epgs-section{{margin-bottom:42px}}
.epgs-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.epgs-warn{{background:rgba(210,153,34,.07);border-left:3px solid var(--ora);
  border-radius:5px;padding:10px 14px;font-size:12.5px;color:var(--muted)}}

/* ==== ANCESTRY ==== */
.anc-section{{margin-bottom:42px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:20px 24px}}
.anc-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:10px}}
.anc-caveat{{background:rgba(210,153,34,.05);border-left:3px solid var(--ora);
  border-radius:5px;padding:10px 14px;font-size:12px;color:var(--muted);
  margin-bottom:16px;line-height:1.65}}
.anc-plot{{text-align:center;margin:16px 0}}
.anc-plot img{{max-width:100%;border-radius:6px;border:1px solid var(--bdr);
  background:#fff}}
.anc-bars{{display:flex;flex-direction:column;gap:6px;margin-top:10px}}
.anc-row{{display:grid;grid-template-columns:160px 1fr 60px;gap:10px;align-items:center;
  font-size:13px}}
.anc-pop{{font-weight:600}}
.anc-bar-wrap{{height:14px;background:var(--bg3);border-radius:7px;overflow:hidden}}
.anc-bar{{height:100%;border-radius:7px;transition:width .4s}}
.anc-bar-eur{{background:#3b82f6}}
.anc-bar-afr{{background:#10b981}}
.anc-bar-eas{{background:#f59e0b}}
.anc-bar-sas{{background:#ef4444}}
.anc-bar-amr{{background:#a855f7}}
.anc-pct{{text-align:right;font-family:var(--mono);font-weight:600}}
.anc-foot{{font-size:11px;color:var(--dim);margin-top:8px;text-align:right}}

/* ==== MEDICATIONS ==== */
.med-section{{margin-bottom:42px;background:linear-gradient(135deg,
  rgba(16,185,129,.04),rgba(59,130,246,.03));border:1px solid rgba(16,185,129,.25);
  border-radius:var(--r);padding:22px 26px}}
.med-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.med-list{{display:flex;flex-direction:column;gap:12px}}
.med-card{{background:var(--bg2);border:1px solid var(--bdr);border-radius:8px;
  padding:14px 18px}}
.med-unknown{{opacity:.7;border-style:dashed}}
.med-head{{display:flex;gap:8px;align-items:baseline;flex-wrap:wrap;margin-bottom:10px;
  padding-bottom:8px;border-bottom:1px solid var(--bdr)}}
.med-drug{{font-size:16px;font-weight:800;color:#10b981}}
.med-generic,.med-note{{font-size:12px;color:var(--muted);font-style:italic}}
.med-msg{{font-size:13px;color:var(--muted);line-height:1.6}}
.med-finding{{padding:8px 0;border-top:1px solid var(--bdr);font-size:13px}}
.med-finding:first-of-type{{border-top:none}}
.med-finding-line{{display:flex;gap:8px;flex-wrap:wrap;align-items:baseline;margin-bottom:4px}}
.med-gene{{font-weight:700;font-family:var(--mono);font-size:13px}}
.med-pheno{{display:inline-block;font-size:11px;font-weight:700;padding:2px 9px;
  border-radius:10px}}
.med-as{{font-size:11px;color:var(--muted);font-family:var(--mono)}}
.med-pathway{{font-size:11.5px;color:var(--muted);margin-bottom:4px}}
.med-rec{{font-size:13px;line-height:1.6;color:var(--txt)}}
.med-guide{{font-size:11px;color:var(--dim);margin-top:4px}}
.pheno-na{{background:var(--bg3);color:var(--muted)}}

/* ==== V4: collapsible category sections ==== */
details.cat-section{{margin-bottom:18px;border-radius:var(--r);
  border:1px solid var(--bdr);background:var(--bg2);overflow:hidden}}
details.cat-section[open]{{padding-bottom:18px}}
.cat-summary{{padding:14px 22px;cursor:pointer;list-style:none;
  display:flex;align-items:center;justify-content:space-between;gap:12px;
  font-size:18px;font-weight:700;border-bottom:1px solid transparent;
  user-select:none;transition:background .15s}}
.cat-summary::-webkit-details-marker{{display:none}}
.cat-summary:hover{{background:var(--bg3)}}
details.cat-section[open] .cat-summary{{border-bottom-color:var(--bdr)}}
.cat-summary::after{{content:"▾";font-size:14px;color:var(--muted);
  transition:transform .15s;margin-left:auto}}
details.cat-section[open] .cat-summary::after{{transform:rotate(180deg)}}
.cat-h2-inline{{font-size:18px;font-weight:700;letter-spacing:-.3px}}
.cat-count{{font-size:11px;font-weight:600;color:var(--muted);
  text-transform:uppercase;letter-spacing:.4px}}
details.cat-section .tbl-wrap{{margin:18px 22px 0}}
details.cat-section .ai-section{{margin:18px 22px 0}}

/* ==== V4: nav-link variant counts ==== */
.nl-count{{display:inline-block;background:var(--bg);color:var(--muted);
  font-size:9px;padding:1px 5px;border-radius:8px;margin-left:4px;
  font-weight:700;vertical-align:middle;font-family:var(--mono)}}
.nl-v4{{background:linear-gradient(135deg,rgba(244,114,182,.18),rgba(168,85,247,.18))!important;
  border-color:rgba(168,85,247,.4)!important;color:var(--txt)!important}}

/* ==== V4: Wellness section ==== */
.wellness-section{{margin-bottom:42px}}
.wellness-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:18px}}
details.wellness-cat{{margin-bottom:14px;background:var(--bg2);
  border:1px solid var(--bdr);border-radius:var(--r);overflow:hidden}}
details.wellness-cat summary{{cursor:pointer;list-style:none;
  padding:12px 18px;display:flex;justify-content:space-between;
  align-items:center;font-size:15px;font-weight:700;
  background:linear-gradient(135deg,rgba(244,114,182,.04),rgba(168,85,247,.03));
  user-select:none}}
details.wellness-cat summary::-webkit-details-marker{{display:none}}
details.wellness-cat summary:hover{{filter:brightness(1.1)}}
.wellness-cat-name{{color:#a855f7}}
.wellness-cat-count{{font-size:11px;color:var(--muted);font-weight:600;
  text-transform:uppercase;letter-spacing:.5px}}
.wellness-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));
  gap:10px;padding:14px 18px 18px}}
.wellness-card{{background:var(--bg);border:1px solid var(--bdr);
  border-left:3px solid #a855f7;border-radius:6px;padding:12px 14px}}
.wellness-trait{{font-size:11px;font-weight:700;color:var(--muted);
  text-transform:uppercase;letter-spacing:.5px;margin-bottom:5px}}
.wellness-result{{font-size:14px;font-weight:700;line-height:1.4;margin-bottom:6px}}
.wellness-action{{font-size:12.5px;color:var(--txt);line-height:1.55;
  background:var(--bg2);border-left:2px solid #a855f7;padding:6px 10px;
  border-radius:4px;margin-bottom:6px}}
.wellness-action strong{{color:#a855f7}}
.wellness-evidence{{font-size:10.5px;color:var(--dim);font-family:var(--mono)}}

/* ==== V5 PREMIUM — common ==== */
.nl-v5{{background:linear-gradient(135deg,rgba(245,158,11,.18),rgba(239,68,68,.18))!important;
  border-color:rgba(239,68,68,.4)!important;color:var(--txt)!important}}

/* ==== V5: HLA ==== */
.hla-section{{margin-bottom:42px}}
.hla-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:12px}}
.hla-caveat-block{{background:rgba(210,153,34,.07);border-left:3px solid var(--ora);
  padding:10px 14px;border-radius:5px;font-size:12.5px;color:var(--muted);
  margin-bottom:14px;line-height:1.65}}
.hla-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:12px}}
.hla-card{{background:var(--bg2);border:1px solid var(--bdr);border-radius:8px;padding:14px 18px}}
.hla-head{{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:6px}}
.hla-allele{{font-weight:800;font-size:15px;color:#f59e0b}}
.hla-status{{font-size:10px;font-weight:700;padding:2px 8px;border-radius:10px;
  text-transform:uppercase;letter-spacing:.5px}}
.hla-pos{{background:rgba(239,68,68,.18);color:var(--red)}}
.hla-neg{{background:rgba(63,185,80,.18);color:var(--grn)}}
.hla-na{{background:var(--bg3);color:var(--muted)}}
.hla-conf{{font-size:10px;color:var(--muted);font-family:var(--mono)}}
.hla-freq{{font-size:11px;color:var(--muted);margin-bottom:6px}}
.hla-tags{{display:flex;gap:5px;flex-wrap:wrap;margin-bottom:8px}}
.hla-tag{{font-size:10px;background:var(--bg3);padding:2px 7px;border-radius:6px}}
.hla-clinical{{background:var(--bg);border-left:3px solid #f59e0b;padding:8px 12px;
  border-radius:4px;margin-bottom:6px}}
.hla-clinical-label{{font-size:10px;font-weight:700;color:#f59e0b;text-transform:uppercase;
  letter-spacing:.4px;margin-bottom:3px}}
.hla-clinical-text{{font-size:12.5px;line-height:1.55;color:var(--txt)}}
.hla-caveat{{font-size:11px;color:var(--ora);margin-top:6px;font-style:italic}}
.hla-transplant,.hla-infection{{margin-top:14px;font-size:13px;color:var(--muted);
  background:var(--bg2);border:1px solid var(--bdr);border-radius:6px;padding:10px 14px}}
.hla-transplant summary,.hla-infection summary{{cursor:pointer;font-weight:600}}
.hla-transplant p,.hla-infection p,.hla-infection ul{{padding:8px 0 0;line-height:1.7}}
.hla-infection li{{margin-bottom:4px}}

/* ==== V5: ROH ==== */
.roh-section{{margin-bottom:42px}}
.roh-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.roh-stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:10px;margin-bottom:14px}}
.roh-stat{{background:var(--bg2);border:1px solid var(--bdr);border-radius:6px;padding:12px;
  text-align:center}}
.roh-n{{font-size:22px;font-weight:800;color:#f59e0b}}
.roh-l{{font-size:10.5px;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;
  margin-top:4px}}
.roh-context{{background:rgba(245,158,11,.06);border-left:3px solid #f59e0b;
  padding:12px 16px;border-radius:5px;font-size:13.5px;line-height:1.7;margin-bottom:14px}}
.roh-ideogram{{background:#fafafa;border:1px solid var(--bdr);border-radius:8px;
  padding:14px;margin-bottom:14px;overflow:auto}}
@media(prefers-color-scheme:dark){{.roh-ideogram{{background:#0d1117}}}}
.roh-table summary{{cursor:pointer;font-weight:600;padding:10px;background:var(--bg2);
  border:1px solid var(--bdr);border-radius:6px}}
tr.roh-long td{{background:rgba(248,81,73,.08)}}
tr.roh-med td{{background:rgba(210,153,34,.05)}}

/* ==== V5: Local Ancestry ==== */
.la-section{{margin-bottom:42px}}
.la-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:10px}}
.la-caveat-block{{font-size:12px;color:var(--muted);font-style:italic;line-height:1.6;
  margin-bottom:12px;background:var(--bg3);padding:10px 14px;border-radius:5px}}
.la-summary{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:14px;font-size:13px}}
.la-call{{padding:6px 10px;background:var(--bg2);border:1px solid var(--bdr);border-radius:6px}}
.la-sp{{display:inline-block;font-weight:700;padding:1px 7px;border-radius:8px;font-size:11px;color:#fff}}
.la-sp-eur{{background:#3b82f6}}
.la-sp-afr{{background:#f59e0b}}
.la-sp-eas{{background:#10b981}}
.la-sp-sas{{background:#a855f7}}
.la-sp-amr{{background:#ef4444}}
.la-sp-unk{{background:#94a3b8}}
.la-paint{{background:#fafafa;border:1px solid var(--bdr);border-radius:8px;padding:14px;
  margin-bottom:12px;overflow:auto}}
@media(prefers-color-scheme:dark){{.la-paint{{background:#0d1117}}}}
.la-deviant{{margin-top:8px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:6px;padding:10px 14px}}
.la-deviant summary{{cursor:pointer;font-weight:600;font-size:13px}}
.la-deviant ul{{padding:8px 0 0 18px;font-size:12.5px;line-height:1.6}}
.la-deviant li{{margin-bottom:4px}}

/* ==== V5: PheWAS ==== */
.phewas-section{{margin-bottom:42px}}
.phewas-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:12px}}
.phewas-headline{{background:rgba(245,158,11,.06);border-left:3px solid #f59e0b;
  padding:12px 16px;border-radius:5px;font-size:13px;margin-bottom:14px}}
.phewas-headline ul{{padding-left:22px;margin-top:6px}}
.phewas-headline li{{margin-bottom:4px}}
details.phewas-cat{{margin-bottom:10px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:6px}}
.ph-detail{{margin-top:6px}}
.ph-detail>summary{{cursor:pointer;font-size:11px;font-weight:600;color:var(--accent,#2563eb)}}
.ph-ref,.ph-z{{font-size:11px;color:var(--muted);margin:4px 0}}
.ph-ref em,.ph-z em{{color:var(--txt);font-style:normal;font-weight:600}}
.ph-detail table{{margin-top:6px;font-size:11px}}
.ph-missing{{font-size:11px;color:var(--warn,#b45309);margin-top:6px;line-height:1.5}}
details.phewas-cat summary{{cursor:pointer;list-style:none;padding:11px 16px;
  display:flex;justify-content:space-between;font-weight:700}}
details.phewas-cat summary::-webkit-details-marker{{display:none}}
.phewas-cat-name{{color:#f59e0b}}
.phewas-cat-count{{font-size:11px;color:var(--muted);font-weight:600}}

/* ==== V5: MR ==== */
.mr-section{{margin-bottom:42px}}
.mr-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:8px}}
.mr-caveat-block{{font-size:12px;color:var(--muted);font-style:italic;line-height:1.6;
  margin-bottom:14px;background:rgba(210,153,34,.06);padding:10px 14px;
  border-left:3px solid var(--ora);border-radius:5px}}
.mr-list{{display:flex;flex-direction:column;gap:12px}}
.mr-finding{{background:var(--bg2);border:1px solid var(--bdr);border-left:4px solid #ef4444;
  border-radius:6px;padding:12px 16px}}
.mr-head{{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:6px}}
.mr-exposure{{font-weight:700;font-size:14px}}
.mr-arrow{{color:var(--muted)}}
.mr-outcome{{font-weight:700;font-size:14px;color:#ef4444}}
.mr-stats{{display:flex;gap:14px;font-size:12px;color:var(--muted);
  margin-bottom:6px;flex-wrap:wrap}}
.mr-stats strong{{color:var(--txt)}}
.mr-rr{{display:inline-block;padding:1px 8px;border-radius:8px;font-weight:700;font-size:11px}}
.mr-coverage{{font-family:var(--mono);font-size:11px}}
.mr-explain{{font-size:12.5px;line-height:1.65;color:var(--txt);margin-top:4px}}
.mr-ref{{font-size:10.5px;color:var(--dim);margin-top:6px}}

/* ==== V5: Genetic Age ==== */
.ga-section{{margin-bottom:42px;background:linear-gradient(135deg,
  rgba(245,158,11,.08),rgba(239,68,68,.05));border:1px solid rgba(245,158,11,.35);
  border-radius:8px;padding:22px 26px}}
.ga-headline{{display:flex;align-items:center;gap:18px;margin-bottom:14px}}
.ga-headline-pct{{font-size:60px;font-weight:800;color:#f59e0b;line-height:1;
  letter-spacing:-1.5px;font-family:var(--mono)}}
.ga-pct-th{{font-size:24px;color:var(--muted);font-weight:600}}
.ga-headline-meta{{font-size:13px;line-height:1.55}}
.ga-years{{font-size:18px;font-weight:700;color:var(--acc);margin-top:4px}}
.ga-narr{{font-size:14px;line-height:1.75;margin-bottom:14px}}
.ga-subs{{display:flex;gap:10px;margin-bottom:12px;flex-wrap:wrap}}
.ga-sub{{flex:1;min-width:160px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:5px;padding:10px}}
.ga-sub-name{{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.4px}}
.ga-sub-val{{font-size:18px;font-weight:700;color:#f59e0b;margin-top:3px}}
.ga-disc{{font-size:11px;color:var(--muted);font-style:italic;line-height:1.65}}

/* ==== V5: PGx Simulation ==== */
.sim-section{{margin-bottom:42px}}
.sim-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:14px}}
.sim-list{{display:flex;flex-direction:column;gap:12px}}
.sim-drug{{background:var(--bg2);border:1px solid var(--bdr);border-radius:8px;
  padding:14px 18px}}
.sim-head{{display:flex;justify-content:space-between;align-items:center;
  margin-bottom:10px;flex-wrap:wrap;gap:8px}}
.sim-name{{font-size:16px;font-weight:800;color:#ef4444}}
.sim-gene-pill{{font-size:10px;background:var(--bg3);padding:3px 9px;border-radius:10px;
  color:var(--muted);font-family:var(--mono)}}
.sim-metrics{{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin-bottom:10px}}
.sim-metric{{background:var(--bg);border:1px solid var(--bdr);border-radius:5px;
  padding:10px;text-align:center}}
.sim-metric-name{{font-size:10px;color:var(--muted);text-transform:uppercase;
  letter-spacing:.4px;margin-bottom:4px}}
.sim-metric-val{{font-size:20px;font-weight:800;font-family:var(--mono)}}
.sim-genes{{margin-bottom:10px;font-size:12.5px;line-height:1.7}}
.sim-gene-line{{margin-bottom:4px}}
.sim-gene-metrics{{font-size:11px;color:var(--muted);font-variant-numeric:tabular-nums}}
.sim-ddi summary{{cursor:pointer;font-size:12px;color:var(--muted);font-weight:600}}
.sim-ddi ul{{padding:8px 0 0 20px;font-size:12px;color:var(--muted)}}
.sim-ddi li{{margin-bottom:4px;line-height:1.55}}

/* ==== V5: Reproductive ==== */
.rep-section{{margin-bottom:42px}}
.rep-intro{{font-size:13px;color:var(--muted);line-height:1.7;margin-bottom:10px}}
.rep-roh-context{{background:rgba(245,158,11,.08);border-left:3px solid #f59e0b;
  padding:10px 14px;border-radius:5px;font-size:13px;margin-bottom:14px;line-height:1.65}}
.rep-list{{display:flex;flex-direction:column;gap:12px}}
.rep-scenario{{background:var(--bg2);border:1px solid var(--bdr);
  border-left:4px solid #f472b6;border-radius:8px;padding:12px 16px}}
.rep-head{{display:flex;gap:8px;align-items:baseline;flex-wrap:wrap;margin-bottom:4px}}
.rep-gene{{font-weight:800;font-size:14px;color:#f472b6}}
.rep-variant{{font-size:12px;color:var(--muted)}}
.rep-inh{{font-size:10px;text-transform:uppercase;letter-spacing:.4px;color:var(--muted);
  background:var(--bg3);padding:2px 7px;border-radius:8px}}
.rep-disease{{font-size:13.5px;font-weight:600;margin-bottom:8px}}
.rep-pop-tbl{{margin-bottom:8px}}
.rep-dom,.rep-x{{font-size:12.5px;color:var(--muted);font-style:italic;
  background:var(--bg3);padding:8px 12px;border-radius:4px;margin-bottom:6px}}
.rep-advice,.rep-outlook{{font-size:12.5px;color:var(--txt);line-height:1.6;margin-top:4px}}
.rep-advice strong,.rep-outlook strong{{color:#f472b6}}
/* ==== Health Economics ==== */
.econ-section{{margin:30px 0}}
.econ-intro{{font-size:13px;color:var(--muted);line-height:1.6;margin-bottom:14px;max-width:760px}}
.econ-tbl{{width:100%;border-collapse:collapse;font-size:12.5px}}
.econ-tbl th{{text-align:left;padding:8px 10px;background:var(--bg3);color:var(--muted);
  font-size:11px;text-transform:uppercase;letter-spacing:.5px;border-bottom:1px solid var(--bdr)}}
.econ-tbl td{{padding:9px 10px;border-bottom:1px solid var(--bdr);vertical-align:top}}
.econ-benefit{{font-size:11.5px;color:var(--muted);margin-top:3px;line-height:1.5}}
.econ-sub{{font-size:10.5px;color:var(--dim);text-transform:uppercase;letter-spacing:.4px}}
.econ-roi{{color:var(--grn);font-size:14px;white-space:nowrap}}
.econ-headline{{background:var(--bg2);border:1px solid var(--bdr);border-left:3px solid var(--grn);
  border-radius:var(--r);padding:14px 18px;margin-bottom:16px}}
.econ-headline h3{{font-size:13px;text-transform:uppercase;letter-spacing:.6px;
  color:var(--muted);margin-bottom:8px}}
.econ-headline-list{{margin:0;padding-left:18px;font-size:13px;line-height:1.7}}
.econ-conf{{font-size:10px;text-transform:uppercase;letter-spacing:.4px;padding:2px 7px;
  border-radius:8px;background:var(--bg3);color:var(--muted)}}
.econ-conf-high{{background:rgba(63,185,80,.18);color:var(--grn)}}
.econ-conf-moderate{{background:rgba(210,153,34,.18);color:var(--ora)}}
.econ-conf-low{{background:var(--bg3);color:var(--dim)}}
.econ-panel{{margin-top:20px;background:var(--bg2);border:1px solid var(--bdr);
  border-radius:var(--r);padding:16px 18px}}
.econ-panel h3{{font-size:14px;margin-bottom:12px;display:flex;align-items:center;gap:10px}}
.econ-pop{{font-size:11px;font-weight:500;color:var(--muted);background:var(--bg3);
  padding:2px 9px;border-radius:8px}}
.econ-cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:10px}}
.econ-card{{background:var(--bg3);border-radius:6px;padding:12px;text-align:center}}
.econ-card-v{{font-size:18px;font-weight:700;color:var(--acc2)}}
.econ-card-l{{font-size:11px;color:var(--muted);margin-top:4px}}
.econ-model{{font-size:12px;color:var(--muted);line-height:1.6;margin-top:12px}}
.econ-disclaimer{{font-size:11.5px;color:var(--dim);line-height:1.6;margin-top:16px;
  font-style:italic;border-top:1px solid var(--bdr);padding-top:12px}}
</style>
</head>
<body>

<header class="hdr">
<div class="hdr-inner">
  <div>
    <div class="hdr-title">DNA Analysis Report <span class="hdr-version">v{REPORT_VERSION}</span></div>
    <div class="hdr-meta">Generated: {report_date}&nbsp;&middot;&nbsp;
    Source: {dna_filename}&nbsp;&middot;&nbsp;{ai_tier_note}{qc_chip}</div>
  </div>
</div>
</header>

<div class="wrap">

<div class="disc">
<strong>&#9888;&#65039; Informational and educational use only &mdash; not a
clinical diagnostic.</strong>
This report integrates a curated variant database (400 SNPs), CPIC-style
pharmacogenomic phenotyping, curated polygenic risk scores, compound
heterozygosity detection, carrier-status analysis, trait predictions, and
optional local-AI interpretation. <strong>It is not a CLIA/CAP-certified
clinical laboratory result.</strong> Findings — especially pharmacogenomic
phenotypes and pathogenic variants — should be confirmed by a clinical
laboratory before any medical decision is made. Polygenic risk scores are
curated-variant approximations, not full clinical-grade PGS. Always consult
a qualified physician, clinical pharmacist, or board-certified genetic
counselor before acting on any finding in this report.
<br><br>
<strong>Reading the results:</strong> every score (polygenic, pharmacogenomic,
ancestry, and haplogroup) is labelled with an explicit <em>confidence</em> level
and the number of SNPs actually typed versus expected. Results marked
<em>low confidence</em> or <em>indeterminate</em> were computed on too few
markers to be reliable and should not be interpreted; they are shown only for
transparency.
</div>

<div class="stats">
  <div class="stat"><div class="stat-n">{total_matched}</div>
    <div class="stat-l">Variants Analyzed</div></div>
  <div class="stat"><div class="stat-n">{risk_variants}</div>
    <div class="stat-l">Risk Alleles Found</div></div>
  <div class="stat"><div class="stat-n">{high_count}</div>
    <div class="stat-l">High Significance</div></div>
  <div class="stat"><div class="stat-n">{len(sorted_cats)}</div>
    <div class="stat-l">Categories</div></div>
</div>

<nav class="toc">
<h3>Quick Navigation</h3>
<div class="toc-links">
  {"<a href='#quality-control' class='nl nl-pro'>Quality Control</a>" if qc_html else ""}
  {"<a href='#imputation' class='nl nl-v3'>Imputation</a>" if imputation_html else ""}
  <a href="#apoe" class="nl">APOE Genotype</a>
  {"<a href='#ancestry' class='nl nl-v3'>Ancestry</a>" if ancestry_html else ""}
  {"<a href='#deep-ancestry' class='nl nl-v5'>Deep Ancestry</a>" if deep_ancestry_html else ""}
  {"<a href='#blood-type' class='nl nl-v5'>Blood Type</a>" if blood_type_html else ""}
  {"<a href='#immunogenetics' class='nl nl-v5'>Immunogenetics</a>" if immunogenetics_html else ""}
  {"<a href='#neurochemistry' class='nl nl-v5'>Neurochemistry</a>" if neurochemistry_html else ""}
  {"<a href='#addiction-genetics' class='nl nl-v5'>Addiction Genetics</a>" if addiction_genetics_html else ""}
  {"<a href='#clinical-variants' class='nl nl-v5'>Clinical Variants</a>" if clinical_variants_html else ""}
  {"<a href='#novel-variants' class='nl nl-v5'>Novel Variants</a>" if novel_variants_html else ""}
  {"<a href='#value-of-information' class='nl nl-v5'>Value of Information</a>" if voi_html else ""}
  {"<a href='#family-planning' class='nl nl-v5'>Family Planning</a>" if family_planning_html else ""}
  {"<a href='#polygenic-traits' class='nl nl-v5'>Trait Genetics</a>" if polygenic_traits_html else ""}
  {"<a href='#environmental-optimization' class='nl nl-v5'>Environmental Optimization</a>" if environmental_optimization_html else ""}
  {"<a href='#life-stage-playbook' class='nl nl-v5'>Life-Stage Playbook</a>" if life_stage_playbook_html else ""}
  {"<a href='#ancestral-story' class='nl nl-v5'>Ancestral Story</a>" if ancestral_story_html else ""}
  <a href="#y-haplogroup" class="nl">Y-DNA Haplogroup</a>
  {"<a href='#mt-haplogroup' class='nl'>mtDNA Haplogroup</a>" if mt_result else ""}
  {"<a href='#counseling-triggers' class='nl nl-pro'>Consultation Triggers</a>" if counseling_html else ""}
  {"<a href='#exec-summary' class='nl'>Executive Summary</a>" if not no_ai and exec_summary else ""}
  {"<a href='#holistic-synthesis' class='nl nl-v5'>Holistic Synthesis</a>" if holistic_synthesis_html else ""}
  {"<a href='#polygenic-risk-scores' class='nl nl-pro'>Polygenic Risk Scores</a>" if prs_html else ""}
  {"<a href='#expanded-pgs' class='nl nl-v3'>Expanded PGS</a>" if expanded_pgs_html else ""}
  {"<a href='#pharmacogenomics' class='nl nl-pro'>Pharmacogenomics</a>" if pgx_html else ""}
  {"<a href='#pharmgkb-clinical' class='nl nl-pro'>PGx Variant Annotations</a>" if pharmgkb_html else ""}
  {"<a href='#top-drugs' class='nl nl-pro'>Top Drugs Screen</a>" if top_drugs_html else ""}
  {"<a href='#medication-review' class='nl nl-v3'>Medication Review</a>" if medications_html else ""}
  {"<a href='#health-economics' class='nl nl-pro'>Health Economics</a>" if economics_html else ""}
  {"<a href='#variant-interactions' class='nl nl-pro'>Variant Interactions</a>" if interactions_html else ""}
  {"<a href='#carrier-status' class='nl nl-pro'>Carrier Status</a>" if carrier_html else ""}
  {"<a href='#trait-predictions' class='nl nl-pro'>Trait Predictions</a>" if traits_html else ""}
  {"<a href='#wellness' class='nl nl-v4'>Wellness</a>" if wellness_html else ""}
  {"<a href='#detoxification' class='nl nl-v5'>Detoxification</a>" if detox_html else ""}
  {"<a href='#urologic' class='nl nl-v5'>Urologic Panel</a>" if urologic_html else ""}
  {"<a href='#metal-oxidative' class='nl nl-v5'>Metal &amp; Oxidative</a>" if metal_oxidative_html else ""}
  {"<a href='#genetic-age' class='nl nl-v5'>Genetic Age</a>" if genetic_age_html else ""}
  {"<a href='#hla-immune' class='nl nl-v5'>HLA / Immune</a>" if hla_html else ""}
  {"<a href='#pgx-simulation' class='nl nl-v5'>Drug Simulation</a>" if pgx_sim_html else ""}
  {"<a href='#phewas' class='nl nl-v5'>PheWAS</a>" if phewas_html else ""}
  {"<a href='#mendelian-randomization' class='nl nl-v5'>MR Projections</a>" if mr_html else ""}
  {"<a href='#reproductive' class='nl nl-v5'>Reproductive</a>" if reproductive_html else ""}
  {"<a href='#runs-of-homozygosity' class='nl nl-v5'>ROH</a>" if roh_html else ""}
  {"<a href='#local-ancestry' class='nl nl-v5'>Local Ancestry</a>" if local_ancestry_html else ""}
  {nav_links}
  {"<a href='#cross-category' class='nl'>Cross-Category Interactions</a>" if cross_cat_html else ""}
  <a href="#recommendations" class="nl">Recommendations</a>
  {"<a href='#references' class='nl nl-pro'>References &amp; Evidence</a>" if references_html else ""}
  <a href="#methodology" class="nl">Methodology</a>
</div>
</nav>

{qc_html}

{imputation_html}

{apoe_html}

{ancestry_html}

{deep_ancestry_html}

{blood_type_html}

{immunogenetics_html}

{neurochemistry_html}

{addiction_genetics_html}

{clinical_variants_html}
{novel_variants_html}
{voi_html}

{family_planning_html}

{polygenic_traits_html}
{tnrc18_html}

{environmental_optimization_html}

{life_stage_playbook_html}

{ancestral_story_html}

{ydna_section_html}

{mtdna_section_html}

{counseling_html}

{exec_html}

{holistic_synthesis_html}

{prs_html}

{expanded_pgs_html}

{pgx_html}

{pharmgkb_html}

{top_drugs_html}

{medications_html}

{economics_html}

{interactions_html}

{carrier_html}

{traits_html}

{wellness_html}

{detox_html}

{urologic_html}

{metal_oxidative_html}

{genetic_age_html}

{hla_html}

{pgx_sim_html}

{phewas_html}

{mr_html}

{reproductive_html}

{roh_html}

{local_ancestry_html}

{cat_sections_html}

{cross_cat_html}

{rec_section}

{references_html}

<section class="meth-section" id="methodology">
<h2>Methodology</h2>
<p><strong>Tier 1 &mdash; Deterministic SNP Lookup:</strong> Your raw DNA file
was parsed using the <code>snps</code> Python library, which auto-detects
formats from 23andMe, AncestryDNA, TellmeGen, MyHeritage, and others. Each
SNP in your file was matched against a curated database of clinically and
scientifically studied variants, with risk allele counts and significance
levels based on published research.</p>
<p><strong>Tier 2 &mdash; Local AI Interpretation:</strong>
{"Results were sent in category batches to the <strong>" + model + "</strong> model "
"running locally via Ollama. No genetic data was transmitted to any external server. "
"The AI provides educational interpretation of variant interactions, prioritization, "
"and recommendations &mdash; not diagnosis." if not no_ai else
"Tier 2 AI analysis was skipped (<code>--no-ai</code> flag). "
"Run without <code>--no-ai</code> with Ollama running (<code>ollama serve</code>) "
"to get AI interpretation."}</p>
<p><strong>APOE Determination:</strong> APOE genotype is derived by combining
rs429358 (C allele = E4-defining) and rs7412 (T allele = E2-defining), covering
all common haplotypes (E2/E2 through E4/E4).</p>
<p><strong>Y-DNA Haplogroup:</strong> Y-chromosome SNPs are filtered from the
parsed data and matched against a decision tree of haplogroup-defining markers
focused on Macro-haplogroup K and all downstream subclades (K1, K2a, K2b &rarr;
N, O, Q, P &rarr; R, R1, R1a, R1b and its European branches: U106, P312, L21,
DF27, U152). Markers are looked up by rsID first and by GRCh37 chromosome Y
position as fallback. When a defining SNP is absent from the chip, the tree
continues into children and marks the gap node as &ldquo;inferred.&rdquo; Consumer
chip Y-coverage varies significantly by provider and array version; many
deep haplogroup markers require FTDNA Big Y-700 or equivalent long-read
sequencing for resolution.</p>
<p><strong>Privacy:</strong> All analysis runs entirely on your local machine.
No genetic data is transmitted to any external service. The only network call
(if Tier 2 is enabled) is to <code>localhost:11434</code> (Ollama).</p>
<p><strong>Polygenic Risk Scores:</strong> Each panel is a curated-variant
weighted log-OR sum normalised to an expected Hardy-Weinberg distribution
(European reference allele frequencies), translated into a Z-score and
percentile. Variants are drawn from CARDIoGRAMplusC4D, DIAGRAM, IGAP,
BCAC, PRACTICAL, AFGen, GIANT and related GWAS consortia. These are
<em>curated</em>, not full clinical PGS; effect magnitudes likely
under-estimate full polygenic burden.</p>
<p><strong>Pharmacogenomic Phenotyping:</strong> CPIC activity-score model.
Per-gene, each contributing variant's function impact is summed
(LOF = −1.0, reduced = −0.5, increased = +0.5) from baseline 2.0 to yield
an activity score, classified into Poor / Intermediate / Normal / Rapid /
Ultra-rapid metabolizer phenotypes. Drug recommendations follow published
CPIC guidelines. CYP2D6 ultra-rapid metabolizer status requires CNV
detection that SNP arrays cannot perform. HLA-B*57:01 calls from proxy
SNPs MUST be confirmed by direct HLA typing before any clinical use.</p>
<p><strong>Compound Heterozygosity / Variant Interactions:</strong>
Rule-based detection of clinically meaningful multi-variant patterns —
MTHFR C677T+A1298C compound, HFE C282Y/H63D, F5 Leiden + Prothrombin,
multi-locus thrombophilia, multi-9p21 + Lp(a), Crohn's polygenic
patterns, ALDH2 + ADH1B, and other compound effects.</p>
<p><strong>Carrier Status:</strong> Chip-detectable carrier and affected
status for selected recessive/dominant pathogenic variants. Not equivalent
to clinical carrier-screening panels (which use sequencing to cover ~250+
recessive conditions).</p>
<p><strong>Trait Predictions:</strong> Concrete genotype-based phenotype
calls for lactose persistence, alcohol flush, caffeine metabolism, bitter
taste, earwax type, eye/hair color, chronotype, short-sleeper allele,
muscle fiber composition, vitamin D synthesis efficiency, caffeine-anxiety
susceptibility, and nicotine-dependence susceptibility.</p>
<p><strong>Privacy:</strong> All analysis runs entirely on your local
machine. No genetic data is transmitted to any external service. The only
network call (if Tier 2 AI is enabled) is to <code>localhost:11434</code>
(Ollama).</p>
<p><strong>Limitations:</strong> Findings depend on chip coverage — see the
Quality Control section for callability per clinical domain. Many genetic
influences require deeper testing: rare pathogenic variants (clinical
sequencing); CYP2D6 CNVs (specialised PGx assays); full mitochondrial
sequencing (FTDNA mtFull); deep Y-DNA subclades (FTDNA Big Y-700);
clinical-grade PGS (large biobank-derived scores); imprinting, epigenetics,
and gene-environment interactions. Significance ratings reflect current
evidence, which continues to evolve.</p>
</section>

</div>

<footer class="ftr">
DNA Analysis Report &nbsp;&middot;&nbsp; Generated locally &amp; privately &nbsp;&middot;&nbsp;
Educational use only &mdash; not medical advice
<div class="ftr-cite">
  <strong>Suggested citation:</strong> Aque, C.&nbsp;R. ({datetime.datetime.now():%Y}).
  <em>DNA Analysis Tool</em> (Version {REPORT_VERSION}) [Computer software].
</div>
<div class="ftr-copy">&copy; {datetime.datetime.now():%Y} Conall&nbsp;R.&nbsp;Aque. All rights reserved.</div>
</footer>
</body>
</html>"""

    return html

