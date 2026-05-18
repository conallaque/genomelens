"""
Emergency Medical Genetics Card
===============================

A standalone one-page HTML document with ONLY the clinically actionable
findings an emergency clinician needs to know. Designed to:
  * Print on a wallet card
  * Be screenshotted on a phone
  * Display in landscape on a tablet
Print-friendly CSS optimised for laser/inkjet output.

Items shown (in priority order):
  1. Severe drug contraindications (HLA-B*57:01 → abacavir; B*15:02 → carbamazepine;
     B*58:01 → allopurinol; A*31:01 → carbamazepine; DPYD → 5-FU; BCHE → succinylcholine)
  2. Bleeding / clotting disorders (Factor V Leiden, Prothrombin G20210A)
  3. Malignant hyperthermia / anesthesia risks
  4. G6PD deficiency (hemolytic triggers)
  5. CYP extremes (PM/UM for codeine, tramadol, clopidogrel, warfarin)
  6. Hemochromatosis status
  7. APOE ε4/ε4 (Alzheimer's awareness — informational for clinicians)
  8. High-significance carrier findings

Items NOT shown: complex disease PRS, traits, ancestry — these aren't
emergency-relevant.
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Dict, List, Optional


def _esc(s):
    if s is None:
        return ""
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _filter_emergency_findings(
    carrier_result: Optional[Dict],
    pgx_result: Optional[Dict],
    pgx_sim_result: Optional[Dict],
    hla_result: Optional[Dict],
    interactions_result: Optional[Dict],
    apoe_genotype: Optional[str],
) -> List[Dict]:
    """Pull actionable, emergency-relevant findings."""
    items: List[Dict] = []

    # ── HLA drug hypersensitivities ──
    if hla_result:
        for a in hla_result.get("alleles", []):
            if a["status"] not in ("carrier (heterozygous)", "homozygous"):
                continue
            for label, text in a.get("clinical", []):
                if label.startswith("Drug"):
                    items.append({
                        "priority": 1,
                        "category": "Drug Contraindication",
                        "title": f"{a['allele']} — {label}",
                        "detail": text,
                    })

    # ── PGx extremes ──
    if pgx_result:
        for gene, r in pgx_result.get("per_gene", {}).items():
            code = r.get("phenotype_code")
            if code not in ("PM", "UM", "POS"):
                continue
            items.append({
                "priority": 2,
                "category": "PGx Extreme",
                "title": f"{gene} {r.get('phenotype', '')}",
                "detail": ", ".join(d["drug"] for d in r.get("drug_recs", []) if d.get(code)),
            })

    # ── Quantitative PGx simulation flags ──
    if pgx_sim_result and pgx_sim_result.get("drugs"):
        for ds in pgx_sim_result["drugs"]:
            if ds["combined_ae_rr"] >= 2.5 or ds["combined_dose_factor"] == 0:
                items.append({
                    "priority": 2,
                    "category": "Drug Response Alert",
                    "title": ds["drug"],
                    "detail": (
                        f"AE risk {ds['combined_ae_rr']}×, recommended dose factor "
                        f"{ds['combined_dose_factor']}×. " +
                        "; ".join(g["note"] for g in ds["gene_findings"])
                    ),
                })

    # ── Thrombophilia (from carrier or compound interactions) ──
    if carrier_result:
        for c in carrier_result.get("carriers", []) + carrier_result.get("affected", []):
            disease = c.get("disease", "").lower()
            if "thromboembolism" in disease or "thrombophilia" in disease:
                items.append({
                    "priority": 1,
                    "category": "Clotting Disorder",
                    "title": f"{c['gene']} {c['variant']} carrier",
                    "detail": c.get("affected_implication") if c.get("dosage") == 2 else c.get("carrier_implication"),
                })
            elif "g6pd" in disease:
                items.append({
                    "priority": 1,
                    "category": "Hemolytic Risk (G6PD)",
                    "title": f"G6PD deficiency — {c['variant']}",
                    "detail": "AVOID: sulfa drugs (sulfamethoxazole), primaquine, dapsone, nitrofurantoin, fava beans, naphthalene mothballs.",
                })
            elif "hemochromatosis" in disease and c.get("dosage") == 2:
                items.append({
                    "priority": 3,
                    "category": "Iron Overload Risk",
                    "title": f"HFE {c['variant']} homozygous",
                    "detail": "May have iron-overload phenotype. Avoid iron supplements and high-dose vitamin C with iron-rich meals.",
                })

    if interactions_result:
        for f in interactions_result.get("findings", []):
            if f["severity"] == "high" and "thrombophilia" in f["title"].lower():
                items.append({
                    "priority": 1,
                    "category": "Clotting Disorder",
                    "title": f["title"],
                    "detail": f["action"],
                })

    # APOE ε4/ε4
    if apoe_genotype == "E4/E4":
        items.append({
            "priority": 4,
            "category": "Informational — Genetic AD Risk",
            "title": "APOE ε4/ε4",
            "detail": "Highest-risk APOE genotype for Alzheimer's. Inform any new neurologist; consider risk-reduction discussion."
        })

    # Sort by priority
    items.sort(key=lambda x: x["priority"])
    return items


def build_emergency_card(
    file_label: str,
    file_hash: str,
    version: str,
    apoe_genotype: Optional[str] = None,
    carrier_result: Optional[Dict] = None,
    pgx_result: Optional[Dict] = None,
    pgx_sim_result: Optional[Dict] = None,
    hla_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    report_link: str = "report.html",
) -> str:
    items = _filter_emergency_findings(
        carrier_result, pgx_result, pgx_sim_result,
        hla_result, interactions_result, apoe_genotype,
    )

    if not items:
        items_html = (
            '<div class="ec-none">No emergency-relevant genetic findings detected. '
            'Routine clinical care applies.</div>'
        )
    else:
        cards = []
        for it in items:
            cls = {
                1: "ec-crit",
                2: "ec-warn",
                3: "ec-info",
                4: "ec-info",
            }.get(it["priority"], "ec-info")
            cards.append(f"""
<div class="ec-finding {cls}">
  <div class="ec-cat">{_esc(it["category"])}</div>
  <div class="ec-title">{_esc(it["title"])}</div>
  <div class="ec-detail">{_esc(it["detail"])}</div>
</div>
""")
        items_html = "".join(cards)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Emergency Medical Genetics Card</title>
<style>
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  background:#fff;color:#1f2328;line-height:1.45;font-size:13.5px;
  padding:14mm 12mm}}
.ec-header{{border-bottom:3px solid #cf222e;padding-bottom:10px;margin-bottom:14px;
  display:flex;justify-content:space-between;align-items:flex-end;gap:12px}}
.ec-banner{{font-size:11px;font-weight:800;letter-spacing:.8px;text-transform:uppercase;
  color:#cf222e;background:rgba(207,34,46,.08);padding:6px 12px;border-radius:4px;
  border:1.5px solid #cf222e;display:inline-block;margin-bottom:6px}}
.ec-title-main{{font-size:22px;font-weight:800;letter-spacing:-.3px;color:#1f2328}}
.ec-sub{{font-size:11px;color:#57606a;margin-top:3px}}
.ec-disclaimer{{background:#fff7e6;border-left:3px solid #d29922;padding:8px 12px;
  margin-bottom:14px;font-size:11px;line-height:1.55;color:#57606a;border-radius:3px}}
.ec-list{{display:flex;flex-direction:column;gap:8px}}
.ec-finding{{border-radius:5px;padding:8px 12px;border-left:4px solid #999;
  background:#f6f8fa;page-break-inside:avoid}}
.ec-crit{{border-left-color:#cf222e;background:#fff5f5}}
.ec-warn{{border-left-color:#d29922;background:#fff8e6}}
.ec-info{{border-left-color:#0969da;background:#f0f7ff}}
.ec-cat{{font-size:9.5px;font-weight:800;color:#57606a;text-transform:uppercase;
  letter-spacing:.4px;margin-bottom:2px}}
.ec-crit .ec-cat{{color:#cf222e}}
.ec-warn .ec-cat{{color:#d29922}}
.ec-info .ec-cat{{color:#0969da}}
.ec-title{{font-size:14px;font-weight:700;margin-bottom:3px;color:#1f2328}}
.ec-detail{{font-size:12px;color:#57606a;line-height:1.55}}
.ec-none{{padding:14px;background:#f6f8fa;border-radius:5px;color:#57606a;font-size:12px}}
.ec-footer{{margin-top:18px;padding-top:10px;border-top:1px solid #d0d7de;
  font-size:10px;color:#888;display:flex;justify-content:space-between;gap:12px;
  flex-wrap:wrap}}
.ec-footer code{{font-family:"SF Mono",Consolas,monospace;font-size:9px}}
@media print {{
  body {{padding:8mm 8mm}}
  .ec-finding {{box-shadow:none;border:1px solid #999;border-left-width:4px}}
}}
</style>
</head>
<body>

<div class="ec-header">
  <div>
    <div class="ec-banner">⚠ Show to Emergency Medical Personnel</div>
    <div class="ec-title-main">Emergency Medical Genetics Card</div>
    <div class="ec-sub">Generated {_dt.datetime.now().strftime("%B %d, %Y")} · {_esc(version)}</div>
  </div>
  <div style="text-align:right;font-size:10px;color:#888">
    See full report: <code>{_esc(report_link)}</code>
  </div>
</div>

<div class="ec-disclaimer">
<strong>Clinically actionable genetic findings only.</strong> This card is not a
substitute for clinical judgment. Pharmacogenomic and HLA proxies from a
consumer SNP array should be CONFIRMED by clinical testing before use in
prescribing — particularly for HLA-mediated drug hypersensitivity decisions.
</div>

<div class="ec-list">
{items_html}
</div>

<div class="ec-footer">
  <span>Source file: <code>{_esc(file_label)}</code></span>
  <span>File hash (SHA256, first 16): <code>{_esc(file_hash)}</code></span>
</div>

</body>
</html>"""
