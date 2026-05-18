"""
Genetic Counseling Triggers
---------------------------

Automated detection of findings that, in clinical genetics practice,
warrant referral to a board-certified genetic counselor or
medical-genetics specialist.

Triggers reflect ACMG / NSGC / ClinGen / NCCN guidance for when
chip-based findings cross from 'informative' into 'professional review
warranted'.
"""

from typing import Dict, List


def evaluate_counseling_triggers(
    tier1_results: List[Dict],
    apoe_genotype: str,
    pgx_results: Dict,
    prs_results: Dict,
    carrier_results: Dict,
    interactions_results: Dict,
) -> Dict:
    """Aggregate counseling triggers from all analyses."""
    triggers: List[Dict] = []

    # ── APOE ε4/ε4 ───────────────────────────────────────────────────────────
    if apoe_genotype == "E4/E4":
        triggers.append({
            "trigger": "APOE ε4/ε4 — highest-risk Alzheimer's genotype",
            "urgency": "Discuss within 6 months",
            "specialist": "Neurologist or AD-prevention clinic; consider genetic counselor",
            "reason": (
                "Lifetime AD risk substantially elevated. Pre-symptomatic counseling "
                "can address risk-disclosure preferences, insurance implications, "
                "and prevention planning."
            ),
        })

    # ── Hemochromatosis homozygous ───────────────────────────────────────────
    for entry in carrier_results.get("affected", []):
        if "C282Y" in entry.get("variant", ""):
            triggers.append({
                "trigger": "HFE C282Y/C282Y — classical hemochromatosis genotype",
                "urgency": "Within 1–3 months",
                "specialist": "Hepatologist or hematologist",
                "reason": (
                    "Test serum ferritin and transferrin saturation. If elevated, "
                    "evaluation for therapeutic phlebotomy. Periodic monitoring "
                    "if currently normal."
                ),
            })

    # ── Multiple thrombophilia variants ─────────────────────────────────────
    for f in interactions_results.get("findings", []):
        if "Thrombophilia" in f["title"] and f["severity"] == "high":
            triggers.append({
                "trigger": f["title"],
                "urgency": "Before any surgery, pregnancy, or hormone therapy",
                "specialist": "Hematologist",
                "reason": (
                    "Multi-variant thrombophilia with multiplicative risk. "
                    "Decisions about prophylactic anticoagulation in high-risk "
                    "situations require hematology input."
                ),
            })

    # ── CHEK2 / ATM / BRCA-region findings + family history hook ────────────
    for r in tier1_results:
        if r["rsid"] in ("rs17879961", "rs1801516", "rs1799966") and r["risk_copies"] > 0:
            triggers.append({
                "trigger": f"{r['gene']} {r['variant_name']} carrier",
                "urgency": "Within 6 months",
                "specialist": "Cancer genetic counselor (board-certified, NSGC)",
                "reason": (
                    "Moderate-penetrance cancer-susceptibility variant detected. "
                    "Full panel sequencing (BRCA1/2/ATM/PALB2/CHEK2/RAD51 etc.) "
                    "with a board-certified genetic counselor is the standard of "
                    "care to clarify family-history-adjusted risk and screening "
                    "intensity. Consider for self AND blood relatives."
                ),
            })
            break  # one trigger covers the whole panel

    # ── TREM2 R47H ────────────────────────────────────────────────────────────
    for r in tier1_results:
        if r["rsid"] == "rs75932628" and r["risk_copies"] > 0:
            triggers.append({
                "trigger": "TREM2 R47H carrier",
                "urgency": "Within 6 months",
                "specialist": "Neurologist or AD-prevention clinic",
                "reason": (
                    "Rare moderately-penetrant Alzheimer's risk variant. Comparable "
                    "in effect size to one APOE-ε4 allele. Pre-symptomatic counseling "
                    "and prevention planning are appropriate."
                ),
            })
            break

    # ── Pharmacogenomic high-impact phenotypes ──────────────────────────────
    pgx_actionable = pgx_results.get("actionable_findings", [])
    high_impact_drugs = {"Clopidogrel (Plavix)", "Codeine, Tramadol",
                        "Azathioprine, 6-Mercaptopurine, Thioguanine", "Warfarin",
                        "Tamoxifen (endocrine breast cancer therapy)", "Abacavir (HIV)"}
    high_impact_pgx = [a for a in pgx_actionable if any(d in a["drug"] for d in high_impact_drugs)]
    if high_impact_pgx:
        unique_genes = sorted({a["gene"] for a in high_impact_pgx})
        triggers.append({
            "trigger": f"High-impact pharmacogenomic phenotypes: {', '.join(unique_genes)}",
            "urgency": "Before initiating relevant medications",
            "specialist": "Clinical pharmacist, prescriber, or pharmacogenomics-trained clinician",
            "reason": (
                "Detected PGx phenotypes substantially alter drug choice or dose for "
                "high-impact medications. SHARE this report with prescribing clinicians. "
                "For HLA-B*57:01 specifically, confirm with direct HLA typing before "
                "any abacavir consideration."
            ),
        })

    # ── High-PRS findings ───────────────────────────────────────────────────
    high_prs = []
    panels = prs_results.get("panels", {})
    for name, p in panels.items():
        if p.get("result", {}).get("tier") == "High":
            high_prs.append(name)
    if high_prs:
        triggers.append({
            "trigger": f"High polygenic risk: {'; '.join(high_prs)}",
            "urgency": "Discuss with primary care; consider preventive-specialist referral",
            "specialist": "Primary care + condition-specific specialist (cardiology, oncology)",
            "reason": (
                "Top-decile polygenic risk for one or more conditions justifies "
                "discussion of enhanced screening or earlier preventive intervention. "
                "Polygenic risk should be integrated with family history and "
                "established clinical risk factors."
            ),
        })

    # ── Multi-variant FH possibility ────────────────────────────────────────
    cad_panel = panels.get("Coronary Artery Disease", {}).get("result", {})
    if cad_panel.get("tier") in ("Elevated", "High"):
        # Check for personal-action triggers around lipids
        has_lpa = any(r["rsid"] in ("rs10455872", "rs3798220") and r["risk_copies"] > 0
                      for r in tier1_results)
        if has_lpa:
            triggers.append({
                "trigger": "Lp(a) elevation likely + elevated CAD PRS",
                "urgency": "Within 3–6 months",
                "specialist": "Preventive cardiologist or lipidologist",
                "reason": (
                    "Lp(a) is an independent CV risk factor not lowered by statins. "
                    "Measure Lp(a) once (mass or molar). If elevated, aggressive "
                    "non-Lp(a) risk-factor control + monitoring; emerging "
                    "Lp(a)-targeted therapies (siRNA, antisense) are in late "
                    "trials and may be eligible options in the next 1–3 years."
                ),
            })

    return {
        "triggers": triggers,
        "n_triggers": len(triggers),
        "urgent_count": sum(1 for t in triggers if "Within 1" in t.get("urgency", "") or "Before any" in t.get("urgency", "")),
    }
