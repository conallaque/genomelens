"""
Family Planning / Carrier Report
================================

Generates a standalone "Carrier Status Report" focused on findings relevant
to family planning. Wraps carrier.py output with:

  * Severity-ordered organisation (serious childhood conditions first)
  * Carrier-frequency context per condition
  * Partner-testing guidance (which conditions to ask a prospective parent
    to test for)
  * Sensitivity-aware framing ('carrier status is information, not a
    diagnosis; consider genetic counselling for personalised advice')

Carrier frequency data drawn from ClinGen / OMIM / population genetics
literature for common European, Ashkenazi Jewish, African American, and
East Asian populations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List


# Per-condition severity + carrier frequency context. "severity" ranks:
#   "critical_childhood"  — serious, often fatal in childhood without intervention
#   "significant_adult"   — significant adult-onset disease
#   "moderate"            — clinically meaningful but variable penetrance
#   "informational"       — useful context but lower urgency

CONDITION_CONTEXT: Dict[str, Dict] = {
    "Hereditary Hemochromatosis (HH, type 1)": {
        "severity": "significant_adult",
        "carrier_freq": {
            "European": "~1 in 8 (very common; ~12% carrier rate)",
            "African": "Rare",
            "East Asian": "Very rare",
        },
        "child_risk": (
            "If partner is also a C282Y carrier, each child has a 25% chance "
            "of being homozygous (clinical hemochromatosis at-risk genotype). "
            "Penetrance is incomplete — 25-60% of homozygotes develop iron "
            "overload. Treatable with phlebotomy."
        ),
        "partner_testing": "Ask partner to test HFE C282Y (rs1800562) and H63D (rs1799945).",
        "treatment_outlook": "Excellent if detected early — therapeutic phlebotomy is curative for iron overload.",
    },
    "Hereditary Hemochromatosis (mild variant)": {
        "severity": "informational",
        "carrier_freq": {
            "European": "~1 in 4 (very common)",
            "African": "~1 in 20",
            "East Asian": "Uncommon",
        },
        "child_risk": (
            "H63D alone has very low penetrance. Most relevant for offspring "
            "if partner carries C282Y — compound heterozygosity has moderate "
            "iron-loading risk."
        ),
        "partner_testing": "Combine with HFE C282Y testing in partner.",
        "treatment_outlook": "Phlebotomy if iron overload develops.",
    },
    "Venous Thromboembolism (VTE) Susceptibility": {
        "severity": "significant_adult",
        "carrier_freq": {
            "European": "Factor V Leiden ~1 in 20; Prothrombin G20210A ~1 in 50",
            "African": "Rare",
            "East Asian": "Very rare",
        },
        "child_risk": (
            "These are autosomal dominant susceptibility variants — each child "
            "has 50% chance of inheriting if you're heterozygous. Homozygosity "
            "requires both parents to be carriers."
        ),
        "partner_testing": "Factor V Leiden (rs6025) + Prothrombin G20210A (rs1799963) testing in partner if family planning.",
        "treatment_outlook": (
            "Variants affect clot risk in surgery, pregnancy, hormonal contraception. "
            "Manageable with awareness — most carriers never have a VTE event."
        ),
    },
    "Late-Onset Alzheimer's Disease (susceptibility)": {
        "severity": "significant_adult",
        "carrier_freq": {
            "European": "TREM2 R47H ~1 in 200; APOE ε4 ~1 in 7",
            "African": "Different frequencies; less studied",
            "East Asian": "Lower APOE ε4 frequency",
        },
        "child_risk": (
            "Susceptibility variants — children inherit them but disease is "
            "polygenic, late-onset, and modifiable. NOT a deterministic risk."
        ),
        "partner_testing": "Optional. Most clinicians do not recommend testing children for adult-onset susceptibility variants.",
        "treatment_outlook": "Prevention via exercise, Mediterranean diet, sleep, BP/lipid control.",
    },
    "Atopic Dermatitis / Ichthyosis Vulgaris (filaggrin deficiency)": {
        "severity": "moderate",
        "carrier_freq": {
            "European": "FLG R501X ~1 in 25 (~4% carrier rate)",
            "African": "Different variants more common",
            "East Asian": "Different variants more common",
        },
        "child_risk": (
            "Heterozygous children have elevated eczema/atopic march risk. "
            "Homozygous (both parents carriers, 25% chance) typically causes "
            "severe ichthyosis vulgaris — manageable but lifelong."
        ),
        "partner_testing": "If strong atopic family history, FLG common-variant testing in partner.",
        "treatment_outlook": (
            "Treatable. Early aggressive moisturisation reduces atopic march "
            "and food allergy risk in children of carriers."
        ),
    },
    "Hereditary Breast/Colon/Prostate/Kidney Cancer Susceptibility": {
        "severity": "significant_adult",
        "carrier_freq": {
            "European": "CHEK2 I157T ~1 in 200 in Slavic populations; rarer elsewhere",
            "African": "Rare",
            "East Asian": "Rare",
        },
        "child_risk": (
            "Moderate-penetrance dominant susceptibility — 50% inheritance "
            "per child. Penetrance ~1.5–2× breast cancer risk."
        ),
        "partner_testing": "Not standardly recommended.",
        "treatment_outlook": (
            "Enhanced screening (earlier/more frequent mammography, "
            "colonoscopy) can detect cancer early. Full BRCA panel testing "
            "is the gold standard for hereditary breast/ovarian cancer."
        ),
    },
    "Ankylosing Spondylitis / Acute Anterior Uveitis": {
        "severity": "moderate",
        "carrier_freq": {
            "European": "HLA-B*27 ~6-8%",
            "African": "Very rare",
            "East Asian": "~2-6%",
        },
        "child_risk": (
            "Dominant susceptibility — 50% per child. Only ~6-8% of B27 carriers "
            "develop AS."
        ),
        "partner_testing": "Not standardly recommended.",
        "treatment_outlook": "Manageable with exercise, anti-inflammatory therapy, biologics if needed.",
    },
    "Abacavir Hypersensitivity": {
        "severity": "informational",
        "carrier_freq": {
            "European": "~6-8%",
            "African": "~3-4%",
            "East Asian": "<1%",
        },
        "child_risk": (
            "Only relevant if child ever needs HIV antiretroviral therapy. "
            "50% inheritance per child."
        ),
        "partner_testing": "Not relevant unless HIV-positive partner.",
        "treatment_outlook": "Easy to avoid — many alternative HIV drugs.",
    },
    "Broad Autoimmunity Susceptibility (T1D, RA, lupus, Graves')": {
        "severity": "moderate",
        "carrier_freq": {
            "European": "~10-15% carrier rate",
            "African": "Lower",
            "East Asian": "Lower",
        },
        "child_risk": (
            "Dominant susceptibility with low penetrance. Most carriers never "
            "develop autoimmune disease."
        ),
        "partner_testing": "Not standardly recommended.",
        "treatment_outlook": "Awareness; anti-inflammatory lifestyle; vitamin D sufficiency.",
    },
    "Celiac Disease Susceptibility": {
        "severity": "informational",
        "carrier_freq": {
            "European": "DQ2 ~30%, DQ8 ~10%",
            "African": "Lower",
            "East Asian": "Lower",
        },
        "child_risk": (
            "Necessary but not sufficient — 95%+ of celiac patients carry DQ2 "
            "or DQ8, but only ~1% of carriers develop celiac. 50% inheritance "
            "per child."
        ),
        "partner_testing": "Optional. Negative for both DQ2 + DQ8 essentially rules out celiac for offspring.",
        "treatment_outlook": "Gluten-free diet is fully effective treatment.",
    },
    "Celiac Disease Susceptibility (second-most-common allele)": {
        "severity": "informational",
        "carrier_freq": {
            "European": "~10%",
            "African": "Lower",
            "East Asian": "Lower",
        },
        "child_risk": "See DQ2.",
        "partner_testing": "See DQ2.",
        "treatment_outlook": "Same as DQ2.",
    },
}


SEVERITY_ORDER = ["critical_childhood", "significant_adult", "moderate", "informational"]
SEVERITY_LABEL = {
    "critical_childhood": "Critical / Childhood-Onset Conditions",
    "significant_adult":  "Significant Adult-Onset Conditions",
    "moderate":           "Moderate Susceptibility Variants",
    "informational":      "Informational Findings",
}


def build_carrier_report(carrier_result: Dict) -> Dict:
    """Reorganise carrier_result by severity + attach context."""
    sections: Dict[str, List[Dict]] = {sev: [] for sev in SEVERITY_ORDER}

    all_findings: List[Dict] = []
    for entry in carrier_result.get("affected", []):
        all_findings.append({**entry, "status_label": "AFFECTED (homozygous)"})
    for entry in carrier_result.get("carriers", []):
        all_findings.append({**entry, "status_label": "Carrier (heterozygous)"})

    for f in all_findings:
        ctx = CONDITION_CONTEXT.get(f["disease"], {})
        sev = ctx.get("severity", "informational")
        sections[sev].append({**f, "context": ctx})

    summary_text = _summary_text(sections, carrier_result)
    return {
        "sections": sections,
        "section_labels": SEVERITY_LABEL,
        "section_order": SEVERITY_ORDER,
        "summary": summary_text,
        "total_findings": len(all_findings),
        "n_untested": carrier_result.get("n_untested", 0),
    }


def _summary_text(sections: Dict, carrier_result: Dict) -> str:
    parts = []
    n_aff = carrier_result.get("n_affected", 0)
    n_car = carrier_result.get("n_carriers", 0)
    n_untested = carrier_result.get("n_untested", 0)
    if n_aff:
        parts.append(
            f"{n_aff} affected (homozygous) genotype{'s' if n_aff != 1 else ''} detected — "
            "see Personal Health Implications section."
        )
    if n_car:
        parts.append(
            f"{n_car} carrier (heterozygous) finding{'s' if n_car != 1 else ''} detected. "
            "Carriers are typically unaffected but can transmit variants to children."
        )
    if not n_aff and not n_car:
        parts.append("No carrier or affected findings detected among tested conditions.")
    parts.append(
        f"{n_untested} additional condition{'s' if n_untested != 1 else ''} could not "
        "be tested on this chip. Comprehensive clinical carrier screening (typically "
        "via sequencing) covers 250+ recessive conditions."
    )
    return " ".join(parts)


# ── Standalone HTML renderer ──────────────────────────────────────────────────

def render_carrier_html(carrier_report: Dict, file_label: str = "",
                        version: str = "v3.0.0") -> str:
    """Render the family-planning carrier report as a standalone HTML page."""
    import datetime as _dt

    sections_html = ""
    for sev in carrier_report["section_order"]:
        entries = carrier_report["sections"][sev]
        if not entries:
            continue
        label = carrier_report["section_labels"][sev]
        sev_class = {
            "critical_childhood": "cr-sev-critical",
            "significant_adult":  "cr-sev-significant",
            "moderate":           "cr-sev-moderate",
            "informational":      "cr-sev-info",
        }.get(sev, "cr-sev-info")
        items = ""
        for e in entries:
            ctx = e.get("context", {})
            freq_lines = ""
            for pop, freq in (ctx.get("carrier_freq", {}) or {}).items():
                freq_lines += f"<li><strong>{pop}:</strong> {_safe(freq)}</li>"
            child_risk = _safe(ctx.get("child_risk", "—"))
            partner_test = _safe(ctx.get("partner_testing", "—"))
            outlook = _safe(ctx.get("treatment_outlook", "—"))
            items += f"""
<div class="cr-finding">
  <div class="cr-finding-head">
    <span class="cr-finding-status">{_safe(e.get("status_label",""))}</span>
    <span class="cr-finding-gene">{_safe(e["gene"])}</span>
    <span class="cr-finding-variant">({_safe(e["variant"])})</span>
  </div>
  <div class="cr-finding-disease">{_safe(e["disease"])}</div>
  <div class="cr-finding-block">
    <strong>What this means:</strong> {_safe(e.get("affected_implication") if "AFFECTED" in (e.get("status_label","")) else e.get("carrier_implication",""))}
  </div>
  <div class="cr-finding-block">
    <strong>For future children:</strong> {child_risk}
  </div>
  <div class="cr-finding-block">
    <strong>Partner testing:</strong> {partner_test}
  </div>
  <div class="cr-finding-block">
    <strong>Treatment / outlook:</strong> {outlook}
  </div>
  <div class="cr-finding-block">
    <strong>Carrier frequencies by population:</strong>
    <ul class="cr-freq-list">{freq_lines or "<li>Not catalogued</li>"}</ul>
  </div>
</div>
"""
        sections_html += f"""
<section class="cr-section {sev_class}">
  <h2>{_safe(label)}</h2>
  <div class="cr-list">{items}</div>
</section>
"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Family Planning / Carrier Status Report</title>
<style>
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  background:#0d1117;color:#e6edf3;line-height:1.6;font-size:14.5px}}
@media(prefers-color-scheme:light){{
  body{{background:#fff;color:#1f2328}}
}}
.cr-wrap{{max-width:920px;margin:0 auto;padding:32px 24px 80px}}
.cr-hdr{{border-bottom:2px solid #30363d;padding-bottom:20px;margin-bottom:30px}}
.cr-title{{font-size:26px;font-weight:800;letter-spacing:-.4px;margin-bottom:6px}}
.cr-subtitle{{font-size:13px;color:#8b949e;font-style:italic}}
.cr-meta{{font-size:11.5px;color:#8b949e;margin-top:10px}}
.cr-summary{{background:#161b22;border:1px solid #30363d;border-radius:8px;
  padding:18px 22px;margin-bottom:28px;font-size:14px;line-height:1.7}}
@media(prefers-color-scheme:light){{
  .cr-summary{{background:#f6f8fa;border-color:#d0d7de}}
  .cr-hdr{{border-color:#d0d7de}}
}}
.cr-disc{{background:rgba(88,166,255,.08);border-left:4px solid #58a6ff;
  border-radius:6px;padding:14px 18px;margin-bottom:30px;font-size:13px;line-height:1.7}}
.cr-section{{margin-bottom:36px}}
.cr-section h2{{font-size:18px;padding-bottom:8px;border-bottom:1px solid #30363d;
  margin-bottom:14px}}
@media(prefers-color-scheme:light){{
  .cr-section h2{{border-color:#d0d7de}}
}}
.cr-sev-critical h2{{color:#f85149}}
.cr-sev-significant h2{{color:#d29922}}
.cr-sev-moderate h2{{color:#58a6ff}}
.cr-sev-info h2{{color:#8b949e}}
.cr-list{{display:flex;flex-direction:column;gap:14px}}
.cr-finding{{background:#161b22;border:1px solid #30363d;
  border-left:4px solid #58a6ff;border-radius:8px;padding:16px 20px}}
.cr-sev-critical .cr-finding{{border-left-color:#f85149}}
.cr-sev-significant .cr-finding{{border-left-color:#d29922}}
.cr-sev-moderate .cr-finding{{border-left-color:#58a6ff}}
.cr-sev-info .cr-finding{{border-left-color:#8b949e}}
@media(prefers-color-scheme:light){{
  .cr-finding{{background:#fff;border-color:#d0d7de}}
}}
.cr-finding-head{{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin-bottom:6px}}
.cr-finding-status{{font-size:10px;font-weight:800;letter-spacing:.6px;
  text-transform:uppercase;background:rgba(248,153,0,.18);color:#d29922;
  padding:3px 9px;border-radius:8px}}
.cr-finding-gene{{font-weight:700;font-size:15px}}
.cr-finding-variant{{font-size:13px;color:#8b949e}}
.cr-finding-disease{{font-size:14px;font-weight:600;margin-bottom:8px}}
.cr-finding-block{{font-size:13px;line-height:1.65;color:#c9d1d9;margin-bottom:6px}}
.cr-finding-block strong{{color:#e6edf3}}
@media(prefers-color-scheme:light){{
  .cr-finding-block{{color:#57606a}}
  .cr-finding-block strong{{color:#1f2328}}
}}
.cr-freq-list{{padding-left:20px;margin-top:4px;font-size:12.5px}}
.cr-counseling{{background:rgba(63,185,80,.07);border:1px solid rgba(63,185,80,.25);
  border-radius:8px;padding:16px 20px;margin-top:30px;font-size:13.5px;line-height:1.7}}
.cr-counseling strong{{color:#3fb950}}
</style>
</head>
<body>
<div class="cr-wrap">

<header class="cr-hdr">
  <div class="cr-title">Family Planning / Carrier Status Report</div>
  <div class="cr-subtitle">A focused look at variants that matter for reproductive decisions</div>
  <div class="cr-meta">Generated {_dt.datetime.now().strftime("%B %d, %Y")} from
  <code>{_safe(file_label)}</code> · {_safe(version)}</div>
</header>

<div class="cr-disc">
<strong>Reading this report:</strong> Being a <em>carrier</em> usually means you
do not have the disease but can pass the variant to children. Whether a child
is affected depends on what variants the other biological parent carries. This
report is informational — it is not a diagnosis, and it is far from a complete
carrier screen. A clinical carrier-screening panel via sequencing covers
~250+ recessive conditions; this chip-based analysis covers only the subset
detectable from common variants on the array.
</div>

<div class="cr-summary">{_safe(carrier_report["summary"])}</div>

{sections_html}

<div class="cr-counseling">
<strong>Considering having children?</strong> For any finding above marked as
significant or critical, talk with a board-certified genetic counsellor before
or during early pregnancy. They can: (a) recommend whether your partner should
be tested for the same conditions, (b) interpret any combined findings, and
(c) discuss reproductive options if both parents carry the same recessive
variant. The National Society of Genetic Counselors (NSGC) maintains a
directory at <code>findageneticcounselor.com</code>. Carrier screening is
typically covered by insurance for those planning pregnancy.
</div>

</div>
</body>
</html>
"""


def _safe(s: object) -> str:
    if s is None:
        return ""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;"))
