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
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
) -> Tuple[str, bool]:
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
Y_ANCIENT_DNA_REFS: Dict[str, str] = {
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


def build_ydna_ancient_block(terminal_hg: str, path: List[Dict]) -> str:
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


def build_ydna_html(y_result: Dict) -> str:
    """Render the Y-DNA haplogroup section for the HTML report."""
    status = y_result.get("status", "no_y_data")
    path = y_result.get("path", [])
    terminal = y_result.get("terminal_haplogroup", "Unknown")
    haplogroup_path = y_result.get("haplogroup_path", "Unknown")
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

    return f"""
<section class="ydna-section" id="y-haplogroup">
<h2>Y-DNA Haplogroup {header_badge}</h2>

<div class="ydna-path-wrap">
  <div class="ydna-path-label">Haplogroup path</div>
  <div class="ydna-crumbs">{crumbs if crumbs else "<em>Could not determine path</em>"}</div>
  <div class="ydna-path-note">
    Terminal haplogroup: <strong>{terminal}</strong>{y_count_note}
    &nbsp;&middot;&nbsp;
    <span class="crumb-ok" style="padding:2px 6px">&#9679; confirmed SNP</span>
    <span class="crumb-gap" style="padding:2px 6px">&#9675; inferred (not on chip)</span>
  </div>
</div>

{marker_table}
{gap_html}
{branch_html}
{migration_html}
{ancient_html}
{further_html}
</section>
"""


# ── mtDNA HTML builder ────────────────────────────────────────────────────────

def build_mtdna_html(mt_result: Dict) -> str:
    """Render the mtDNA haplogroup section for the HTML report."""
    status = mt_result.get("status", "no_data")
    haplogroup = mt_result.get("haplogroup", "Unknown")
    confidence = mt_result.get("confidence", "low")
    matched = mt_result.get("matched_markers", [])
    migration = mt_result.get("migration", "")
    ancient_dna = mt_result.get("ancient_dna", "")
    further = mt_result.get("further_testing", "")
    n_mt_snps = mt_result.get("mt_snp_count", 0)

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

    return f"""
<section class="mtdna-section" id="mt-haplogroup">
<h2>mtDNA Haplogroup <span class="mtdna-badge {conf_badge_cls}">{conf_label}</span></h2>

<div class="mtdna-result">
  <div class="mtdna-call">{haplogroup}</div>
  <div class="mtdna-evidence">
    Based on {len(matched)} mtDNA marker(s) on this chip{mt_count_note}.
    Note: mtDNA haplogroup calls from autosomal-chip data are approximate — full
    mitochondrial sequencing remains the gold standard.
  </div>
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


def build_qc_html(qc: Optional[Dict]) -> str:
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
      Chip: <strong>{_esc(qc['file_format'])}</strong> &nbsp;·&nbsp;
      Inferred sex: <strong>{_esc(qc['inferred_sex'])}</strong>
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


def build_prs_html(prs: Optional[Dict]) -> str:
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
                f'<div class="prs-na">{_esc(result["reason"])}</div>'
                f'</div>'
            )
            continue
        tier = result["tier"]
        tier_class = result["tier_class"]
        pct = result["percentile"]
        z = result["z_score"]
        callability = result["callability"]
        details = (
            f'<div class="prs-details">'
            f'<div><span class="prs-lab">Z-score:</span> <strong>{z:+.2f}</strong></div>'
            f'<div><span class="prs-lab">Percentile:</span> <strong>{pct:.0f}th</strong></div>'
            f'<div><span class="prs-lab">Callable variants:</span> '
            f'{len(result["used"])}/{len(result["used"])+len(result["missing"])} '
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
        panel_cards += f"""
<div class="prs-card">
  <div class="prs-card-head">
    <div class="prs-name">{_esc(name)} <span class="prs-short">{_esc(panel.get("trait_short",""))}</span></div>
    <div class="prs-tier {tier_class}">{_esc(tier)}</div>
  </div>
  <div class="prs-desc">{_esc(panel["description"])}</div>
  {track}
  {details}
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


def build_pgx_html(pgx: Optional[Dict]) -> str:
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
        if result.get("is_binary"):
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
  <div class="pgx-callability">
    Callability: {result['callable_variants']}/{result['total_variants']}
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
</section>
"""


def build_interactions_html(inter: Optional[Dict]) -> str:
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


def build_carrier_html(carr: Optional[Dict]) -> str:
    if not carr:
        return ""

    def render_group(items, title, css_class):
        if not items:
            return ""
        rows = ""
        for c in items:
            rows += f"""
<tr>
<td><strong>{_esc(c["gene"])}</strong></td>
<td>{_esc(c["variant"])}</td>
<td>{_esc(c["disease"])}</td>
<td class="gt-cell">{_esc(c.get("genotype",""))}</td>
<td>{_esc(c["inheritance"])}</td>
<td class="sum-cell">{_esc(c.get("carrier_implication") if items is carr["carriers"] else (c.get("affected_implication") if items is carr["affected"] else "—"))}</td>
</tr>
"""
        return f"""
<h3 class="carr-h3 {css_class}">{_esc(title)} ({len(items)})</h3>
<div class="tbl-wrap">
<table class="snp-tbl">
<thead><tr><th>Gene</th><th>Variant</th><th>Condition</th><th>Genotype</th>
<th>Inheritance</th><th>Implication</th></tr></thead>
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


def build_traits_html(tr: Optional[Dict]) -> str:
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


def build_counseling_html(c: Optional[Dict]) -> str:
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


def build_references_html(refs: Optional[List[Dict]]) -> str:
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
        rows += f"""
<tr>
<td class="rsid-cell"><a href="https://www.ncbi.nlm.nih.gov/snp/{_esc(r["rsid"])}" target="_blank" rel="noopener">{_esc(r["rsid"])}</a></td>
<td><strong>{_esc(r["gene"])}</strong> {_esc(r["variant_name"])}</td>
<td><span class="ref-level {level_class(r["evidence_level"])}">{_esc(r["evidence_level"])}</span></td>
<td class="sum-cell">{_esc(r["evidence_summary"])}{guidelines}</td>
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

def build_imputation_html(imp: Optional[Dict]) -> str:
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


def build_expanded_pgs_html(epgs: Optional[Dict]) -> str:
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
  <div class="prs-details">
    <div><span class="prs-lab">Percentile:</span> <strong>{pct:.0f}th</strong></div>
    <div><span class="prs-lab">Z-score:</span> <strong>{result.get('z_score', 0):+.2f}</strong></div>
    <div><span class="prs-lab">Coverage:</span>
      <strong>{cov.get('pct_callable',0)}%</strong>
      ({cov.get('chip',0)} chip + {cov.get('imputed',0)} imputed / {cov.get('total',0)})</div>
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
                f'<div class="prs-na">{_esc(reason)}</div></div>'
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
<div class="prs-grid">{cards}</div>
</section>
"""


def build_ancestry_html(anc: Optional[Dict]) -> str:
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
    return f"""
<section class="anc-section" id="ancestry">
<h2>Ancestry Estimation <span class="pro-pill">V3</span></h2>
<p class="anc-intro">
{_esc(anc.get('method', 'Genotype-based ancestry estimate'))}.
Based on {anc.get('n_aims_used', 0)} ancestry-informative markers.
Confidence: <strong>{_esc(anc.get('confidence','—'))}</strong>.
</p>
<div class="anc-caveat">
<strong>Limitations:</strong> This is a rough estimate compared to commercial
ancestry products (which use proprietary panels of tens of thousands of
markers and reference populations). It is most accurate for recent
single-continent ancestry; admixed individuals (and especially those
without recent European ancestry) may see higher uncertainty.
</div>
{plot_html}
<div class="anc-bars">{bars}</div>
<div class="anc-foot"><em>{_esc(anc.get('confidence_note',''))}</em></div>
</section>
"""


# ── V5: premium section renderers ────────────────────────────────────────────

def build_hla_html(h: Optional[Dict]) -> str:
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
            cards += (
                f'<div class="hla-card">'
                f'<div class="hla-head"><span class="hla-allele">{_esc(a["allele"])}</span>{badge}</div>'
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


def build_roh_html(r: Optional[Dict]) -> str:
    if not r:
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


def build_local_ancestry_html(la: Optional[Dict]) -> str:
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
            items += (
                f'<li>chr{_esc(w["chrom"])} {w["start_mb"]:.0f}-{w["end_mb"]:.0f} Mb: '
                f'<span class="la-sp la-sp-{w["call"].lower()}">{_esc(w["call"])}</span> '
                f'({w.get("confidence","?")}){genes}</li>'
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


def build_phewas_html(p: Optional[Dict]) -> str:
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
            rows += (
                f'<tr>'
                f'<td>{_esc(tname)}</td>'
                f'<td><span class="prs-tier {tier_cls}">{_esc(res["tier"])}</span></td>'
                f'<td>{pct:.0f}th</td>'
                f'<td class="gt-cell">{val} {_esc(unit)}</td>'
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
    <thead><tr><th>Trait</th><th>Tier</th><th>Percentile</th><th>Predicted value</th><th>Coverage</th></tr></thead>
    <tbody>{rows}</tbody>
  </table></div>
</details>
"""

    headline_html = ""
    if p.get("headline"):
        items = "".join(
            f'<li><strong>{_esc(h["trait"])}</strong> — '
            f'<span class="prs-tier tier-{"high" if "high" in h["tier"].lower() else "low"}">{_esc(h["tier"])}</span> '
            f'({h["percentile"]:.0f}th percentile)</li>'
            for h in p["headline"][:15]
        )
        headline_html = (
            f'<div class="phewas-headline"><strong>Notable extremes:</strong>'
            f'<ul>{items}</ul></div>'
        )

    return f"""
<section class="phewas-section" id="phewas">
<h2>Predicted Biomarker & Trait Profile <span class="pro-pill">V5</span></h2>
<p class="phewas-intro">
Phenome-wide scoring across {p['n_traits']} traits using curated GWAS effect
sizes — biomarkers (lipids, HbA1c, CRP, vitamin D, hormones), anthropometric,
hematological, behavioural. <strong>Predicted values are genetic-propensity
estimates against a European reference distribution, not measurements.</strong>
{p['n_scored']}/{p['n_traits']} traits had sufficient SNP coverage to score.
</p>
{headline_html}
{cats_html}
</section>
"""


def build_mr_html(m: Optional[Dict]) -> str:
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


def build_genetic_age_html(g: Optional[Dict]) -> str:
    if not g or not g.get("available"):
        return ""
    long_pct = g["longevity"]["percentile"]
    years = g.get("longevity_years_offset", 0)
    direction = g.get("longevity_direction", "")
    sub_blocks = ""
    if g.get("telomere"):
        sub_blocks += (
            f'<div class="ga-sub">'
            f'<div class="ga-sub-name">Telomere Length (genetic proxy)</div>'
            f'<div class="ga-sub-val">{g["telomere"]["percentile"]:.0f}th percentile</div>'
            f'</div>'
        )
    if g.get("skin_aging"):
        sub_blocks += (
            f'<div class="ga-sub">'
            f'<div class="ga-sub-name">Skin Aging</div>'
            f'<div class="ga-sub-val">{g["skin_aging"]["percentile"]:.0f}th percentile</div>'
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
  </div>
</div>
<div class="ga-narr">{_esc(g.get("narrative",""))}</div>
<div class="ga-subs">{sub_blocks}</div>
<div class="ga-disc">{_esc(g.get("disclaimer",""))}</div>
</section>
"""


def build_pgx_sim_html(s: Optional[Dict]) -> str:
    if not s or not s.get("available") or not s.get("drugs"):
        return ""
    items = ""
    for d in s["drugs"]:
        ae = d["combined_ae_rr"]
        df = d["combined_dose_factor"]
        ae_color = "var(--red)" if ae > 2 else ("var(--ora)" if ae > 1.3 else "var(--grn)")
        df_color = "var(--red)" if df == 0 else ("var(--ora)" if df < 0.8 else "var(--grn)")
        gene_lines = "".join(
            f'<div class="sim-gene-line"><strong>{_esc(g["gene"])}</strong> '
            f'<span class="pgx-pheno {("pheno-pm" if g["phenotype_code"]=="PM" else "pheno-um" if g["phenotype_code"]=="UM" else "pheno-im" if g["phenotype_code"]=="IM" else "pheno-rm" if g["phenotype_code"]=="RM" else "pheno-nm")}">{_esc(g["phenotype_code"])}</span> '
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


def build_reproductive_html(r: Optional[Dict]) -> str:
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


def build_wellness_html(w: Optional[Dict]) -> str:
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


def build_medications_html(med: Optional[Dict]) -> str:
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


# ── HTML report builder ───────────────────────────────────────────────────────

def build_html_report(
    tier1_results: List[Dict],
    apoe_genotype: Optional[str],
    ai_results: Dict[str, str],
    exec_summary: Optional[str],
    dna_filepath: str,
    no_ai: bool,
    model: str,
    y_result: Optional[Dict] = None,
    mt_result: Optional[Dict] = None,
    cross_cat_synthesis: Optional[str] = None,
    qc_result: Optional[Dict] = None,
    prs_result: Optional[Dict] = None,
    pgx_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    traits_result: Optional[Dict] = None,
    counseling_result: Optional[Dict] = None,
    references_used: Optional[List[Dict]] = None,
    imputation_info: Optional[Dict] = None,
    expanded_pgs_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
    medications_result: Optional[Dict] = None,
    wellness_result: Optional[Dict] = None,
    hla_result: Optional[Dict] = None,
    roh_result: Optional[Dict] = None,
    local_ancestry_result: Optional[Dict] = None,
    phewas_result: Optional[Dict] = None,
    mr_result: Optional[Dict] = None,
    genetic_age_result: Optional[Dict] = None,
    pgx_sim_result: Optional[Dict] = None,
    reproductive_result: Optional[Dict] = None,
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
            rows += (
                f'<tr class="{row_cls}">'
                f'<td class="rsid-cell">'
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{s["rsid"]}" '
                f'target="_blank" rel="noopener">{s["rsid"]}</a></td>'
                f'<td><strong>{s["gene"]}</strong></td>'
                f'<td>{s["variant_name"]}{xref_badge}</td>'
                f'<td class="gt-cell">{s["my_genotype"]}</td>'
                f'<td>{risk_indicator(s["risk_copies"])}</td>'
                f'<td>{significance_badge(s["significance"])}</td>'
                f'<td class="sum-cell">{s["summary"]}</td>'
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
.pgx-callability{{font-size:11px;color:var(--muted);margin-bottom:10px}}
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
<strong>&#9888;&#65039; Research/Educational Use &mdash; Not a Clinical Diagnostic.</strong>
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
  <a href="#y-haplogroup" class="nl">Y-DNA Haplogroup</a>
  {"<a href='#mt-haplogroup' class='nl'>mtDNA Haplogroup</a>" if mt_result else ""}
  {"<a href='#counseling-triggers' class='nl nl-pro'>Consultation Triggers</a>" if counseling_html else ""}
  {"<a href='#exec-summary' class='nl'>Executive Summary</a>" if not no_ai and exec_summary else ""}
  {"<a href='#polygenic-risk-scores' class='nl nl-pro'>Polygenic Risk Scores</a>" if prs_html else ""}
  {"<a href='#expanded-pgs' class='nl nl-v3'>Expanded PGS</a>" if expanded_pgs_html else ""}
  {"<a href='#pharmacogenomics' class='nl nl-pro'>Pharmacogenomics</a>" if pgx_html else ""}
  {"<a href='#medication-review' class='nl nl-v3'>Medication Review</a>" if medications_html else ""}
  {"<a href='#variant-interactions' class='nl nl-pro'>Variant Interactions</a>" if interactions_html else ""}
  {"<a href='#carrier-status' class='nl nl-pro'>Carrier Status</a>" if carrier_html else ""}
  {"<a href='#trait-predictions' class='nl nl-pro'>Trait Predictions</a>" if traits_html else ""}
  {"<a href='#wellness' class='nl nl-v4'>Wellness</a>" if wellness_html else ""}
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

{ydna_section_html}

{mtdna_section_html}

{counseling_html}

{exec_html}

{prs_html}

{expanded_pgs_html}

{pgx_html}

{medications_html}

{interactions_html}

{carrier_html}

{traits_html}

{wellness_html}

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
</footer>
</body>
</html>"""

    return html

