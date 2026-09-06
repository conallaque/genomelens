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
    diagnosis; consider genetic counselling for personalized advice')

Carrier frequency data drawn from ClinGen / OMIM / population genetics
literature for common European, Ashkenazi Jewish, African American, and
East Asian populations.
"""

from __future__ import annotations

from typing import Any

# Per-condition severity + carrier frequency context. "severity" ranks:
#   "critical_childhood"  — serious, often fatal in childhood without intervention
#   "significant_adult"   — significant adult-onset disease
#   "moderate"            — clinically meaningful but variable penetrance
#   "informational"       — useful context but lower urgency

CONDITION_CONTEXT: dict[str, dict] = {
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


def build_carrier_report(carrier_result: dict) -> dict:
    """Reorganise carrier_result by severity + attach context."""
    sections: dict[str, list[dict]] = {sev: [] for sev in SEVERITY_ORDER}

    all_findings: list[dict] = []
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


def _summary_text(sections: dict, carrier_result: dict) -> str:
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

def render_carrier_html(carrier_report: dict, file_label: str = "",
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
    <ul class="cr-freq-list">{freq_lines or "<li>Not cataloged</li>"}</ul>
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


# ══════════════════════════════════════════════════════════════════════════
# V35 — "Your genome as a family-planning document"
# ══════════════════════════════════════════════════════════════════════════
#
# A richer, in-report reproductive-genetics section (distinct from the
# standalone carrier page above). It reasons quantitatively about what a
# finding means for *future children*, separating three inheritance modes,
# and is careful about two things clinicians care about:
#
#   1. **Transmission probability ≠ disease probability.** A dominant variant
#      you carry is transmitted to a child with 50% probability, but whether
#      that child develops disease also depends on penetrance. These are shown
#      as separate numbers, never collapsed.
#   2. **mtDNA transmission is sex-gated.** Mitochondrial DNA is inherited
#      almost exclusively from the mother. A male passes essentially none of
#      his mtDNA to his children — so a male's mtDNA findings are relevant to
#      his own health and his maternal relatives, NOT to his offspring. A
#      female transmits her mtDNA to *all* of her children.
#
# The recessive child-affected risk uses the standard Hardy-Weinberg
# random-partner calculation:
#
#       P(child affected) ≈ P(partner is a carrier) × 0.25
#
# where P(partner carrier) is the population carrier frequency for the stated
# ancestry (European assumed unless noted). This assumes an unrelated partner
# of the same background and no consanguinity.

# Numeric European heterozygous-carrier frequencies for the recessive
# conditions the chip-based carrier panel can surface, plus clinical
# penetrance ranges. Sources: HFE — Allen 2008 (HealthIron), European
# C282Y carrier ~1 in 8-10; FLG — Palmer 2006, Sandilands 2007.
_RECESSIVE_CARRIER_FREQ = {
    # gene -> {ancestry -> het carrier frequency}, penetrance (low, high),
    #         partner_test, note
    "HFE": {
        "freq": {"European": 0.11, "African American": 0.03, "East Asian": 0.001},
        "penetrance": (0.25, 0.60),
        "partner_test": "HFE C282Y (rs1800562) and H63D (rs1799945)",
        "note": ("Hereditary hemochromatosis is recessive with incomplete "
                 "penetrance. A child who inherits two risk copies is at risk "
                 "of iron overload, but only 25-60% of homozygotes develop it, "
                 "and it is treatable with phlebotomy."),
    },
    "FLG": {
        "freq": {"European": 0.08, "African American": 0.01, "East Asian": 0.03},
        "penetrance": (0.50, 0.90),
        "partner_test": "FLG common loss-of-function panel (R501X, 2282del4)",
        "note": ("Filaggrin loss is semi-dominant: even one copy raises a "
                 "child's eczema / atopic-march risk; two copies typically "
                 "cause ichthyosis vulgaris. Aggressive early emollient use "
                 "in at-risk infants reduces the atopic march."),
    },
    # ── Remaining recessive conditions detectable by carrier.py ──────────────
    # Carrier frequencies are the standard published figures used in clinical
    # carrier screening (ACMG/ACOG panethnic screening guidance and GeneReviews
    # disease chapters). Ancestry keys match the labels the caller passes; the
    # European entry is the fallback when a person's ancestry is unknown.
    # Before these were added, only HFE and FLG had numeric frequencies, so a
    # CFTR/SMN1/HBB carrier — the conditions carrier screening exists for —
    # received no offspring-risk figure at all.
    "CFTR": {
        "freq": {"European": 0.040, "Ashkenazi Jewish": 0.042,
                 "Hispanic": 0.022, "African American": 0.016,
                 "East Asian": 0.011},
        "penetrance": (0.95, 1.0),
        "partner_test": "CFTR expanded panel (ACMG-23 minimum) or full-gene sequencing",
        "note": ("Cystic fibrosis is highly penetrant: a child with two "
                 "pathogenic copies will have CF, though severity varies by "
                 "variant combination. Partner screening is standard of care "
                 "before or early in pregnancy."),
    },
    "HBB": {
        "freq": {"European": 0.008, "Mediterranean": 0.033,
                 "African American": 0.083, "South Asian": 0.040,
                 "East Asian": 0.033},
        "penetrance": (0.95, 1.0),
        "partner_test": ("Hemoglobinopathy evaluation — CBC with MCV, "
                         "hemoglobin electrophoresis, HBB sequencing"),
        "note": ("Two pathogenic HBB copies cause sickle-cell disease or "
                 "beta-thalassaemia depending on the variants. Carrier "
                 "(trait) status is largely benign. Newborn screening "
                 "detects affected infants, and early penicillin "
                 "prophylaxis plus hydroxyurea substantially improve outcomes."),
    },
    "SMN1": {
        "freq": {"European": 0.021, "Ashkenazi Jewish": 0.020,
                 "African American": 0.014, "East Asian": 0.017,
                 "Hispanic": 0.017},
        "penetrance": (0.95, 1.0),
        "partner_test": "SMN1 copy-number (dosage) analysis — not a sequencing test",
        "note": ("Spinal muscular atrophy is highly penetrant, though copy "
                 "number of the SMN2 modifier gene shifts severity. Note "
                 "carrier testing is dosage-based: consumer chips do not "
                 "detect it reliably. Disease-modifying therapy now exists "
                 "and works best when started pre-symptomatically."),
    },
    "HEXA": {
        "freq": {"European": 0.003, "Ashkenazi Jewish": 0.037,
                 "French Canadian": 0.028, "African American": 0.003,
                 "East Asian": 0.003},
        "penetrance": (0.95, 1.0),
        "partner_test": ("HEXA enzyme (hexosaminidase A) assay plus targeted "
                         "variant panel — enzyme testing outperforms genotyping"),
        "note": ("Infantile Tay-Sachs is uniformly fatal in early childhood "
                 "and has no disease-modifying treatment, which is why it is "
                 "on every expanded carrier panel. Enzyme-based partner "
                 "screening is more sensitive than a variant panel."),
    },
    "GBA": {
        "freq": {"European": 0.010, "Ashkenazi Jewish": 0.067,
                 "African American": 0.005, "East Asian": 0.005},
        "penetrance": (0.30, 0.90),
        "partner_test": "GBA targeted panel (N370S, L444P, 84GG, IVS2+1) or sequencing",
        "note": ("Gaucher type 1 penetrance and severity vary widely — some "
                 "homozygotes stay asymptomatic for life. Enzyme replacement "
                 "therapy is effective. Separately, GBA heterozygosity is a "
                 "known Parkinson's risk factor for the carrier themself."),
    },
    "PAH": {
        "freq": {"European": 0.020, "Mediterranean": 0.023,
                 "African American": 0.008, "East Asian": 0.012},
        "penetrance": (0.95, 1.0),
        "partner_test": "PAH full-gene sequencing",
        "note": ("Phenylketonuria is highly penetrant but is the classic "
                 "treatable inborn error: newborn screening catches it "
                 "universally in the US, and dietary phenylalanine "
                 "restriction started early prevents intellectual "
                 "disability almost entirely."),
    },
    "ASPA": {
        "freq": {"European": 0.003, "Ashkenazi Jewish": 0.018,
                 "African American": 0.003, "East Asian": 0.003},
        "penetrance": (0.95, 1.0),
        "partner_test": "ASPA targeted panel (E285A, Y231X) or sequencing",
        "note": ("Canavan disease is a severe, untreatable "
                 "leukodystrophy of infancy. It is part of the standard "
                 "Ashkenazi Jewish carrier panel."),
    },
    "SMPD1": {
        "freq": {"European": 0.004, "Ashkenazi Jewish": 0.011,
                 "African American": 0.004, "East Asian": 0.004},
        "penetrance": (0.90, 1.0),
        "partner_test": "SMPD1 sequencing",
        "note": ("Niemann-Pick type A is severe and fatal in early "
                 "childhood; type B is far milder and compatible with adult "
                 "life. Which one results depends on the specific variant "
                 "pair, so partner results need specialist interpretation."),
    },
    "ACADM": {
        "freq": {"European": 0.015, "African American": 0.005,
                 "East Asian": 0.003},
        "penetrance": (0.80, 1.0),
        "partner_test": "ACADM sequencing (K304E is the common European variant)",
        "note": ("MCAD deficiency is detected by newborn screening and "
                 "managed simply — avoid prolonged fasting, treat illness "
                 "aggressively. Outcomes are excellent when known in "
                 "advance and potentially fatal when it is not."),
    },
    "GALT": {
        "freq": {"European": 0.014, "African American": 0.010,
                 "East Asian": 0.005},
        "penetrance": (0.95, 1.0),
        "partner_test": "GALT enzyme assay plus sequencing",
        "note": ("Classic galactosemia is caught by newborn screening and "
                 "treated with a galactose-restricted diet, which prevents "
                 "the acute neonatal crisis. Some long-term effects can "
                 "persist despite good dietary control."),
    },
    "BTD": {
        "freq": {"European": 0.008, "African American": 0.008,
                 "East Asian": 0.008},
        "penetrance": (0.80, 1.0),
        "partner_test": "Biotinidase enzyme activity assay plus BTD sequencing",
        "note": ("Profound biotinidase deficiency is one of the most "
                 "treatable inborn errors — lifelong oral biotin prevents "
                 "essentially all sequelae. Newborn screening covers it."),
    },
    "ATP7B": {
        "freq": {"European": 0.011, "East Asian": 0.017,
                 "African American": 0.008},
        "penetrance": (0.90, 1.0),
        "partner_test": "ATP7B full-gene sequencing",
        "note": ("Wilson disease presents anywhere from childhood to "
                 "middle age. It is treatable with chelation or zinc, and "
                 "pre-symptomatic detection substantially changes the "
                 "hepatic and neurologic course."),
    },
    "ALDOB": {
        "freq": {"European": 0.014, "African American": 0.005,
                 "East Asian": 0.005},
        "penetrance": (0.90, 1.0),
        "partner_test": "ALDOB targeted panel (A149P, A174D, N334K)",
        "note": ("Hereditary fructose intolerance is managed entirely by "
                 "dietary fructose/sucrose avoidance. Affected children "
                 "typically self-restrict, but unrecognised exposure "
                 "(including some IV fluids and medicines) is dangerous."),
    },
    "BCKDHA": {
        "freq": {"European": 0.007, "Ashkenazi Jewish": 0.008,
                 "African American": 0.005, "East Asian": 0.005},
        "penetrance": (0.95, 1.0),
        "partner_test": "BCKDHA / BCKDHB / DBT sequencing (MSUD is genetically heterogeneous)",
        "note": ("Maple syrup urine disease needs lifelong branched-chain "
                 "amino-acid restriction and carries real risk of metabolic "
                 "crisis during illness. Newborn screening detects it. Note "
                 "MSUD involves three genes — partner screening should cover "
                 "all of them, not BCKDHA alone."),
    },
    "MEFV": {
        "freq": {"European": 0.010, "Mediterranean": 0.100,
                 "Middle Eastern": 0.140, "Ashkenazi Jewish": 0.130,
                 "African American": 0.005, "East Asian": 0.003},
        "penetrance": (0.50, 0.80),
        "partner_test": "MEFV targeted panel (M694V, V726A, M680I, E148Q)",
        "note": ("Familial Mediterranean fever has notably incomplete "
                 "penetrance and is treatable — daily colchicine controls "
                 "attacks and prevents the amyloidosis that drives its "
                 "long-term morbidity. E148Q in particular is of "
                 "questionable clinical significance on its own."),
    },
    "SERPINA1": {
        "freq": {"European": 0.040, "Mediterranean": 0.030,
                 "African American": 0.010, "East Asian": 0.002},
        "penetrance": (0.30, 0.70),
        "partner_test": "SERPINA1 genotyping (S and Z alleles) plus A1AT serum level",
        "note": ("Alpha-1 antitrypsin deficiency is codominant rather than "
                 "cleanly recessive: ZZ individuals face emphysema and liver "
                 "disease risk, but penetrance depends heavily on smoking "
                 "status. Never smoking is by far the largest modifier, and "
                 "augmentation therapy exists for established disease."),
    },
}


# Hereditary cancer syndromes worth partner / cascade discussion, with the
# crucial dominant-vs-recessive distinction that determines whether partner
# screening is even relevant.
_HEREDITARY_CANCER = [
    {"syndrome": "Hereditary Breast & Ovarian Cancer (HBOC)",
     "genes": "BRCA1, BRCA2, PALB2",
     "inheritance": "autosomal dominant",
     "partner_relevance": "low",
     "note": ("Dominant — if you carry a pathogenic variant, each child has "
              "50% risk of inheriting it regardless of partner. Partner "
              "screening is not the point; cascade testing of blood relatives "
              "and your own enhanced screening are. Consumer chips do NOT "
              "reliably detect BRCA frameshift founder variants — clinical "
              "panel sequencing (Invitae/GeneDx/Ambry) is required.")},
    {"syndrome": "Lynch Syndrome (hereditary colorectal/endometrial)",
     "genes": "MLH1, MSH2, MSH6, PMS2, EPCAM",
     "inheritance": "autosomal dominant",
     "partner_relevance": "low",
     "note": ("Dominant, high-penetrance. Cascade testing + early/frequent "
              "colonoscopy is the management. Not chip-detectable — needs "
              "clinical sequencing if family history suggests it.")},
    {"syndrome": "MUTYH-Associated Polyposis (MAP)",
     "genes": "MUTYH",
     "inheritance": "autosomal recessive",
     "partner_relevance": "HIGH",
     "note": ("**The hereditary-cancer syndrome where partner screening "
              "genuinely matters** — MAP is recessive. If you are a MUTYH "
              "carrier, a child is affected only if the partner is also a "
              "carrier (~1-2% of Europeans). Worth partner MUTYH screening if "
              "you're a known carrier.")},
    {"syndrome": "Li-Fraumeni Syndrome",
     "genes": "TP53",
     "inheritance": "autosomal dominant",
     "partner_relevance": "low",
     "note": ("Dominant, very high penetrance, broad cancer spectrum. Rare. "
              "Clinical sequencing only.")},
    {"syndrome": "Hereditary Diffuse Gastric Cancer",
     "genes": "CDH1",
     "inheritance": "autosomal dominant",
     "partner_relevance": "low",
     "note": "Dominant. Clinical sequencing if family history of diffuse gastric / lobular breast cancer."},
]


# A small set of pathogenic mtDNA variants that occasionally appear on
# consumer arrays. mt_haplogroup classifies lineage, not disease, so this is
# a best-effort screen — real mtDNA-disease screening needs mtDNA sequencing.
_MT_DISEASE_VARIANTS = {
    "rs199476104": ("MT-ND4 m.11778G>A", "Leber Hereditary Optic Neuropathy (LHON)"),
    "rs199476112": ("MT-ND1 m.3460G>A", "Leber Hereditary Optic Neuropathy (LHON)"),
    "rs199476118": ("MT-ND6 m.14484T>C", "Leber Hereditary Optic Neuropathy (LHON)"),
}


def _pct(x: float) -> str:
    if x >= 0.1:
        return f"{x*100:.0f}%"
    if x >= 0.01:
        return f"{x*100:.1f}%"
    return f"{x*100:.2f}%"


def _dosage_of(entry: dict) -> int:
    d = entry.get("dosage")
    return d if isinstance(d, int) else (1 if "AFFECTED" not in entry.get("status_label", "") else 2)


def analyze_family_planning(carrier_result: dict | None = None,
                            # Accepted and ignored. The pipeline passes the
                            # tier-1 result LIST here against a Dict
                            # annotation; nothing in this function ever read
                            # it, so the mismatch was silent. Kept in the
                            # signature because callers pass it by keyword,
                            # but annotated honestly rather than removed in a
                            # cleanup commit.
                            tier1_results: Any | None = None,
                            mt_result: dict | None = None,
                            snps_df=None,
                            inferred_sex: str | None = None,
                            ancestry: str = "European") -> dict:
    """Reason quantitatively about what the person's genome means for future
    children. Returns recessive compound-risk items, dominant-transmission
    items (transmission vs penetrance kept separate), sex-gated mtDNA notes,
    and hereditary-cancer partner-screening guidance."""
    carrier_result = carrier_result or {}

    recessive_items: list[dict] = []
    dominant_items: list[dict] = []

    def _iter(entries, status):
        for e in entries or []:
            yield {**e, "_status": status}

    for e in list(_iter(carrier_result.get("carriers", []), "carrier")) + \
             list(_iter(carrier_result.get("affected", []), "affected")):
        inh = str(e.get("inheritance", "")).lower()
        gene = e.get("gene", "")
        disease = e.get("disease", "")
        variant = e.get("variant", "")
        ctx = CONDITION_CONTEXT.get(disease, {})

        is_recessive = "recessive" in inh
        is_semidominant = "semi-dominant" in inh
        is_dominant = "dominant" in inh and not is_semidominant

        if is_recessive or is_semidominant:
            spec = _RECESSIVE_CARRIER_FREQ.get(gene)
            partner_freq = None
            child_affected = None
            penetrance = None
            clinical_low = clinical_high = None
            if spec:
                partner_freq = spec["freq"].get(ancestry, spec["freq"].get("European"))
                penetrance = spec["penetrance"]
                if e["_status"] == "carrier":
                    # child two-copies risk = P(partner carrier) * 0.25
                    child_affected = partner_freq * 0.25
                else:  # affected (homozygous) parent
                    # child two-copies risk = P(partner carrier) * 0.5
                    child_affected = partner_freq * 0.5
                clinical_low = child_affected * penetrance[0]
                clinical_high = child_affected * penetrance[1]
            recessive_items.append({
                "gene": gene, "variant": variant, "disease": disease,
                "status": e["_status"], "inheritance": e.get("inheritance", ""),
                "ancestry": ancestry,
                "partner_carrier_freq": partner_freq,
                "child_two_copy_risk": child_affected,
                "penetrance": penetrance,
                "child_clinical_risk": (clinical_low, clinical_high)
                    if clinical_low is not None else None,
                "partner_test": (spec or {}).get("partner_test")
                    or ctx.get("partner_testing"),
                "note": (spec or {}).get("note") or ctx.get("child_risk", ""),
                "semidominant": is_semidominant,
            })
        elif is_dominant:
            # Extract a penetrance descriptor from the inheritance string
            if "incomplete" in inh:
                pen_txt = "incomplete / low penetrance"
            elif "moderate" in inh:
                pen_txt = "moderate penetrance"
            elif "susceptibility" in inh:
                pen_txt = "susceptibility only (low penetrance)"
            else:
                pen_txt = "penetrance varies"
            dominant_items.append({
                "gene": gene, "variant": variant, "disease": disease,
                "status": e["_status"], "inheritance": e.get("inheritance", ""),
                "transmission": 0.50,     # per child, if you are heterozygous
                "penetrance_text": pen_txt,
                "note": ctx.get("child_risk", ""),
                "partner_test": ctx.get("partner_testing", "Not standardly recommended."),
                "outlook": ctx.get("treatment_outlook", ""),
            })

    # ── mtDNA — sex-gated ────────────────────────────────────────────────
    sex = (inferred_sex or "").upper()[:1]
    mt_variants_found = []
    if snps_df is not None:
        for rsid, (label, disease) in _MT_DISEASE_VARIANTS.items():
            try:
                if rsid in snps_df.index:
                    gt = snps_df.loc[rsid]
                    if hasattr(gt, "iloc"):
                        gt = gt.iloc[0] if getattr(gt, "ndim", 1) > 1 else gt
                    g = str(gt.get("genotype") if hasattr(gt, "get") else gt).upper()
                    mt_variants_found.append({"rsid": rsid, "label": label,
                                              "disease": disease, "genotype": g})
            except Exception:
                continue
    if sex == "M":
        mt_transmission = ("You are male: your mitochondrial DNA is **not** "
                           "passed to your children (mtDNA is inherited "
                           "maternally). Any mtDNA finding is relevant to your "
                           "own health and to your mother's-line relatives "
                           "(sisters, maternal cousins), not to your offspring.")
    elif sex == "F":
        mt_transmission = ("You are female: **all** of your children inherit "
                           "your mitochondrial DNA. mtDNA disease variants, if "
                           "present, transmit to every child (with "
                           "heteroplasmy-dependent severity).")
    else:
        mt_transmission = ("Mitochondrial DNA is inherited maternally: only a "
                           "female transmits it to her children; a male does "
                           "not. (Sex was not determined for this sample.)")
    mt_block = {
        "sex": sex or "unknown",
        "haplogroup": (mt_result or {}).get("haplogroup"),
        "transmission_note": mt_transmission,
        "pathogenic_variants": mt_variants_found,
        "screening_note": ("Consumer chips test very few of the ~90 known "
                           "pathogenic mtDNA point mutations. A clear result "
                           "here does not rule out mtDNA disease — that needs "
                           "dedicated mtDNA sequencing."),
    }

    n_actionable = len([x for x in recessive_items
                        if x.get("child_two_copy_risk")]) + len(dominant_items)

    return {
        "available": bool(recessive_items or dominant_items or mt_variants_found
                          or carrier_result),
        "ancestry_assumption": ancestry,
        "recessive_items": recessive_items,
        "dominant_items": dominant_items,
        "mtdna": mt_block,
        "hereditary_cancer": _HEREDITARY_CANCER,
        "n_recessive": len(recessive_items),
        "n_dominant": len(dominant_items),
        "n_actionable": n_actionable,
        "summary": _fp_summary(recessive_items, dominant_items, mt_block, sex),
    }


def _fp_summary(recessive, dominant, mt_block, sex: str) -> str:
    bits = []
    if recessive:
        with_risk = [r for r in recessive if r.get("child_two_copy_risk")]
        if with_risk:
            top = max(with_risk, key=lambda r: r["child_two_copy_risk"])
            bits.append(
                f"You carry {len(recessive)} recessive/semi-dominant finding(s). "
                f"With a random {top['ancestry']} partner, the highest child "
                f"two-copy risk is ~{_pct(top['child_two_copy_risk'])} "
                f"({top['gene']} — {top['disease']}), before penetrance.")
        else:
            bits.append(f"You carry {len(recessive)} recessive/semi-dominant finding(s).")
    if dominant:
        bits.append(
            f"{len(dominant)} dominant susceptibility variant(s) transmit to each "
            "child with 50% probability — but transmission is not the same as "
            "disease; penetrance for these is low-to-moderate.")
    if sex == "M":
        bits.append("Your mtDNA does not pass to your children (you are male).")
    elif sex == "F":
        bits.append("Your mtDNA passes to all your children (you are female).")
    if not bits:
        bits.append("No reproductive-relevant carrier findings detected on this chip.")
    bits.append("This is chip-based and far from a complete carrier screen; a "
                "clinical panel covers 250+ recessive conditions.")
    return " ".join(bits)
