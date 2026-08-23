"""
Natural Language Full Report
============================

After all analyses complete, this module compiles a structured summary of
the most important findings and sends them to Ollama for a single
comprehensive narrative — written as if a genetic counselor were
explaining results to a patient.

Output: narrative_report.html (standalone, self-contained).
"""

from __future__ import annotations

import datetime as _dt
import re
from pathlib import Path


def _esc(s):
    if s is None:
        return ""
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


SYSTEM_PROMPT = (
    "You are a board-certified genetic counsellor writing a comprehensive but "
    "READABLE personal genetic report for an individual patient. You are NOT "
    "writing a technical journal article — you are writing as a counsellor "
    "would explain results in person.\n\n"
    "Tone: warm, professional, empowering — not alarming. Always explain what "
    "a genotype means in plain language. Be honest about uncertainty. Avoid "
    "absolute claims; genetics is one input among many.\n\n"
    "Structure your response with clear markdown headings (## Section Name):\n\n"
    "1. ## Top Priorities — Things to Act On\n"
    "   The 3-5 most important actionable findings that could affect health "
    "decisions today.\n\n"
    "2. ## Pharmacogenomic Profile\n"
    "   What this person's DNA says about how they metabolise major medications. "
    "Highlight any drug to avoid or use cautiously.\n\n"
    "3. ## Carrier Status & Family Planning\n"
    "   Recessive conditions they may transmit, with practical partner-testing "
    "advice if relevant.\n\n"
    "4. ## Disease Susceptibility\n"
    "   Major polygenic and single-variant risk findings (cardiovascular, "
    "diabetes, cancer, neurodegeneration). Frame each as 'genetic context, "
    "not destiny'.\n\n"
    "5. ## Personal Traits & Wellness\n"
    "   The genuinely interesting day-to-day findings — caffeine metabolism, "
    "alcohol tolerance, sleep patterns, fitness type, taste preferences.\n\n"
    "6. ## Ancestry & Heritage\n"
    "   What the DNA says about the user's heritage at a high level.\n\n"
    "7. ## Your Top 10 Action Items\n"
    "   A prioritised, specific, numbered list. Concrete actions — not vague "
    "advice. e.g. 'Get one-time Lp(a) test', 'Methylfolate supplement instead "
    "of folic acid'. Aim for 10 items, ranked by life impact.\n\n"
    "Always include: 'This is educational information, not medical advice. "
    "Discuss with your healthcare provider before any medication change or "
    "major health decision.'\n\n"
    "Reference specific genotypes when discussing findings, e.g. 'Your "
    "APOE-ε3/ε4 genotype …' — but always immediately translate the "
    "implication into plain English."
)


def _build_context(
    tier1_results: list[dict],
    apoe_genotype: str | None,
    y_result: dict | None,
    mt_result: dict | None,
    pgx_result: dict | None,
    prs_result: dict | None,
    interactions_result: dict | None,
    carrier_result: dict | None,
    traits_result: dict | None,
    wellness_result: dict | None,
    ancestry_result: dict | None,
    hla_result: dict | None,
    roh_result: dict | None,
    phewas_result: dict | None,
    mr_result: dict | None,
    genetic_age_result: dict | None,
    counseling_result: dict | None,
) -> str:
    """Build the long context string we send to the LLM."""
    out: list[str] = []
    out.append("=== PATIENT GENETIC SUMMARY ===\n")

    # Overall
    if apoe_genotype:
        out.append(f"APOE genotype: {apoe_genotype}")
    if y_result:
        out.append(f"Y-DNA haplogroup: {y_result.get('haplogroup_path', 'Unknown')}")
    if mt_result:
        out.append(f"mtDNA haplogroup: {mt_result.get('haplogroup', 'Unknown')}")

    # Risk-carrying variants
    risk_variants = [r for r in tier1_results if r.get("risk_copies", 0) > 0]
    high = [r for r in risk_variants if r.get("significance") == "high"]
    mod = [r for r in risk_variants if r.get("significance") == "moderate"]
    if high:
        out.append("\nHIGH-significance variants present:")
        for r in high[:20]:
            out.append(f"  - {r['gene']} {r['variant_name']} (cat {r['category']}, "
                       f"{r['risk_copies']} risk allele(s))")
    if mod:
        out.append(f"\nModerate-significance variants ({len(mod)} total, top 15):")
        for r in mod[:15]:
            out.append(f"  - {r['gene']} {r['variant_name']} (cat {r['category']}, "
                       f"{r['risk_copies']} risk allele(s))")

    # PGx
    if pgx_result:
        out.append("\nPharmacogenomic phenotypes:")
        for gene, r in pgx_result.get("per_gene", {}).items():
            phen = r.get("phenotype", "")
            if phen:
                out.append(f"  - {gene}: {phen}")
        actionable = pgx_result.get("actionable_findings", [])
        if actionable:
            out.append(f"  ({len(actionable)} actionable drug-gene findings)")

    # PRS tiers
    if prs_result:
        out.append("\nPolygenic risk score tiers:")
        for name, p in prs_result.get("panels", {}).items():
            tier = p.get("result", {}).get("tier")
            pct = p.get("result", {}).get("percentile")
            if tier:
                out.append(f"  - {name}: {tier} ({pct}th percentile)")

    # Carrier
    if carrier_result:
        n_aff = carrier_result.get("n_affected", 0)
        n_car = carrier_result.get("n_carriers", 0)
        if n_aff or n_car:
            out.append(f"\nCarrier/affected findings: {n_aff} affected, {n_car} carrier(s)")
            for c in carrier_result.get("affected", []):
                out.append(f"  AFFECTED: {c['gene']} {c['variant']} - {c['disease']}")
            for c in carrier_result.get("carriers", []):
                out.append(f"  Carrier: {c['gene']} {c['variant']} - {c['disease']}")

    # Compound interactions
    if interactions_result and interactions_result.get("findings"):
        out.append("\nCompound / variant interactions:")
        for f in interactions_result["findings"][:8]:
            out.append(f"  [{f['severity'].upper()}] {f['title']}")

    # HLA
    if hla_result and hla_result.get("carrier_alleles"):
        out.append(f"\nHLA carrier alleles: {', '.join(hla_result['carrier_alleles'])}")

    # Traits / wellness highlights
    if traits_result and traits_result.get("predictions"):
        out.append("\nTrait predictions (selected high-confidence):")
        for t in traits_result["predictions"]:
            if t.get("confidence") in ("high", "moderate") and t.get("result") != "Not tested":
                out.append(f"  - {t.get('trait')}: {t.get('result')}")

    if wellness_result and wellness_result.get("predictions"):
        out.append("\nWellness signals:")
        for p in wellness_result["predictions"]:
            out.append(f"  - [{p.get('category')}] {p.get('trait')}: {p.get('result')}")

    # PheWAS headline
    if phewas_result and phewas_result.get("headline"):
        out.append("\nNotable predicted biomarkers (top extremes):")
        for h in phewas_result["headline"][:10]:
            out.append(f"  - {h['trait']}: {h['tier']} ({h['percentile']}th percentile)")

    # MR highlights
    if mr_result and mr_result.get("findings"):
        out.append("\nMR causal projections (top):")
        sig = sorted(
            [f for f in mr_result["findings"] if f["status"] == "ok"],
            key=lambda f: abs(f.get("outcome_shift_log_or") or 0), reverse=True,
        )[:5]
        for f in sig:
            out.append(f"  - {f['exposure']} → {f['outcome']}: "
                       f"RR ≈ {f['outcome_relative_risk']}")

    # Genetic age
    if genetic_age_result and genetic_age_result.get("available"):
        out.append(f"\nGenetic longevity profile: {genetic_age_result['longevity']['percentile']}th "
                   f"percentile, {genetic_age_result['longevity_years_offset']:+.1f}y "
                   f"vs European mean.")

    # ROH
    if roh_result and roh_result.get("f_roh") is not None:
        out.append(f"\nROH burden: F_ROH = {roh_result['f_roh']:.4f} "
                   f"({roh_result.get('context_tier','')})")

    # Ancestry — these are relative single-population affinities, NOT admixture
    # proportions, and the heuristic is low/moderate confidence. Frame it so the
    # AI does not narrate it as "X% of your DNA".
    if ancestry_result and ancestry_result.get("primary_population"):
        conf = ancestry_result.get("confidence", "low")
        amb = " (ambiguous — cannot be reliably distinguished from the runner-up)" \
            if ancestry_result.get("ambiguous") else ""
        out.append(
            f"\nAncestry best single-population match: "
            f"{ancestry_result['primary_population']} "
            f"({conf} confidence{amb}). Note: rough affinity from a small marker "
            "panel, not admixture percentages."
        )

    # Counseling triggers
    if counseling_result and counseling_result.get("triggers"):
        out.append("\nProfessional consultation triggers:")
        for t in counseling_result["triggers"][:5]:
            out.append(f"  - {t.get('trigger')}")

    return "\n".join(out)


def generate_narrative_report(
    output_path: Path,
    file_label: str,
    file_hash: str,
    version: str,
    model: str,
    tier1_results: list[dict],
    apoe_genotype: str | None = None,
    y_result: dict | None = None,
    mt_result: dict | None = None,
    pgx_result: dict | None = None,
    prs_result: dict | None = None,
    interactions_result: dict | None = None,
    carrier_result: dict | None = None,
    traits_result: dict | None = None,
    wellness_result: dict | None = None,
    ancestry_result: dict | None = None,
    hla_result: dict | None = None,
    roh_result: dict | None = None,
    phewas_result: dict | None = None,
    mr_result: dict | None = None,
    genetic_age_result: dict | None = None,
    counseling_result: dict | None = None,
) -> str:
    """Send a comprehensive summary to Ollama and write narrative_report.html."""
    import requests

    context = _build_context(
        tier1_results=tier1_results,
        apoe_genotype=apoe_genotype, y_result=y_result, mt_result=mt_result,
        pgx_result=pgx_result, prs_result=prs_result,
        interactions_result=interactions_result, carrier_result=carrier_result,
        traits_result=traits_result, wellness_result=wellness_result,
        ancestry_result=ancestry_result, hla_result=hla_result,
        roh_result=roh_result, phewas_result=phewas_result,
        mr_result=mr_result, genetic_age_result=genetic_age_result,
        counseling_result=counseling_result,
    )

    prompt = SYSTEM_PROMPT + "\n\n" + context + "\n\nPlease produce the narrative report now."

    try:
        resp = requests.post(
            "http://localhost:11434/api/generate",
            json={
                "model": model,
                "messages": [{"role": "user", "content": prompt}],
                "stream": False,
                "options": {"temperature": 0.45, "num_ctx": 12288},
            },
            timeout=900,
        )
        resp.raise_for_status()
        narrative = resp.json()["message"]["content"]
        # Strip <think>...</think> blocks
        narrative = re.sub(r"<think>.*?</think>", "", narrative, flags=re.DOTALL).strip()
    except Exception as e:
        narrative = f"*Narrative generation failed: {e}*"

    html = _render_narrative_html(narrative, file_label, file_hash, version, model)
    output_path.write_text(html, encoding="utf-8")
    return narrative


def _render_narrative_html(narrative_md: str, file_label: str, file_hash: str,
                           version: str, model: str) -> str:
    body = _md_to_html(narrative_md)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Personal Genetic Narrative</title>
<style>
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:Georgia,"Times New Roman",serif;background:#fafaf6;color:#1f2328;
  line-height:1.75;font-size:16.5px}}
@media(prefers-color-scheme:dark){{body{{background:#0d1117;color:#e6edf3}}}}
.narr-wrap{{max-width:760px;margin:0 auto;padding:48px 28px 80px}}
.narr-header{{border-bottom:2px solid #c2a76e;padding-bottom:24px;margin-bottom:30px}}
.narr-pill{{display:inline-block;background:#c2a76e;color:#fff;padding:3px 12px;
  border-radius:20px;font-size:11px;font-weight:700;letter-spacing:1px;text-transform:uppercase;
  margin-bottom:10px}}
.narr-title{{font-size:38px;font-weight:700;letter-spacing:-.6px;margin-bottom:8px;font-family:Georgia,serif}}
.narr-sub{{font-size:14px;color:#888;font-style:italic}}
.narr-disc{{background:rgba(194,167,110,.10);border-left:4px solid #c2a76e;
  padding:16px 22px;margin-bottom:32px;font-size:14.5px;line-height:1.7;
  border-radius:4px;font-style:italic}}
.narr-body h2{{font-size:22px;font-weight:700;margin:36px 0 16px;
  padding-bottom:8px;border-bottom:1px solid #d0d7de;color:#5d4a1f}}
@media(prefers-color-scheme:dark){{.narr-body h2{{color:#c2a76e;border-color:#30363d}}}}
.narr-body h3{{font-size:17px;font-weight:700;margin:24px 0 10px}}
.narr-body p{{margin-bottom:14px;line-height:1.85}}
.narr-body ul,.narr-body ol{{padding-left:28px;margin-bottom:14px}}
.narr-body li{{margin-bottom:8px;line-height:1.7}}
.narr-body strong{{color:#5d4a1f;font-weight:700}}
@media(prefers-color-scheme:dark){{.narr-body strong{{color:#c2a76e}}}}
.narr-body code{{font-family:"SF Mono",Consolas,monospace;font-size:.85em;
  background:rgba(194,167,110,.15);padding:1px 6px;border-radius:3px}}
.narr-footer{{margin-top:40px;padding-top:18px;border-top:1px solid #d0d7de;
  font-size:11px;color:#888;text-align:center}}
.narr-footer code{{font-family:"SF Mono",Consolas,monospace;font-size:10px}}
</style>
</head>
<body>
<div class="narr-wrap">

<header class="narr-header">
  <div class="narr-pill">Personal Genetic Narrative</div>
  <h1 class="narr-title">Your Genetic Story</h1>
  <div class="narr-sub">Written as a genetic counsellor would explain in conversation.</div>
  <div class="narr-sub" style="margin-top:6px">Generated {_dt.datetime.now().strftime("%B %d, %Y")} · {_esc(version)} · model: {_esc(model)}</div>
</header>

<div class="narr-disc">
This narrative interprets your raw DNA data alongside curated variant
databases, polygenic scores, pharmacogenomic phenotypes, ancestry, and
trait predictions. It is educational only — not medical advice. Discuss
significant findings with a physician or board-certified genetic counsellor
before acting on them.
</div>

<div class="narr-body">
{body}
</div>

<footer class="narr-footer">
  <div>Source file: <code>{_esc(file_label)}</code> &middot; SHA-256: <code>{_esc(file_hash)}</code></div>
</footer>

</div>
</body>
</html>
"""


def _md_to_html(text: str) -> str:
    """Minimal markdown → HTML for the narrative body."""
    # Bold + italics
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*\n]+?)\*(?!\*)", r"<em>\1</em>", text)
    # Code spans
    text = re.sub(r"`([^`\n]+)`", r"<code>\1</code>", text)
    # Headings
    text = re.sub(r"^###\s+(.+)$", r"<h3>\1</h3>", text, flags=re.MULTILINE)
    text = re.sub(r"^##\s+(.+)$", r"<h2>\1</h2>", text, flags=re.MULTILINE)
    # Convert lists and paragraphs in a simple way
    lines = text.split("\n")
    html_lines = []
    in_ol = False
    in_ul = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if in_ol:
                html_lines.append("</ol>"); in_ol = False
            if in_ul:
                html_lines.append("</ul>"); in_ul = False
            continue
        m_ol = re.match(r"^\d+\.\s+(.+)$", stripped)
        m_ul = re.match(r"^[-*•]\s+(.+)$", stripped)
        if m_ol:
            if in_ul:
                html_lines.append("</ul>"); in_ul = False
            if not in_ol:
                html_lines.append("<ol>"); in_ol = True
            html_lines.append(f"<li>{m_ol.group(1)}</li>")
        elif m_ul:
            if in_ol:
                html_lines.append("</ol>"); in_ol = False
            if not in_ul:
                html_lines.append("<ul>"); in_ul = True
            html_lines.append(f"<li>{m_ul.group(1)}</li>")
        else:
            if in_ol:
                html_lines.append("</ol>"); in_ol = False
            if in_ul:
                html_lines.append("</ul>"); in_ul = False
            # Heading already converted
            if stripped.startswith("<h"):
                html_lines.append(stripped)
            else:
                html_lines.append(f"<p>{stripped}</p>")
    if in_ol:
        html_lines.append("</ol>")
    if in_ul:
        html_lines.append("</ul>")
    return "\n".join(html_lines)
