"""
Quantitative Pharmacogenomic Simulation
=======================================

Layered ON TOP of pgx.py's phenotype calls. Adds:
  * Relative clearance vs population (%)
  * Approximate effective-dose adjustment factor
  * Side-effect probability vs population (×)
  * Drug-drug-gene interaction notes
  * "Virtual prescriber" summary AI prompt context

This is illustrative — the numeric estimates are derived from CPIC
guideline ranges and published PK studies, applied at the activity-score
level. Real prescribing requires clinical assessment.
"""

from __future__ import annotations

from typing import Dict, List, Optional


# Per-drug PK / response model. Each entry:
#   drug: name
#   gene_pheno_map: gene → phenotype → (rel_clearance, dose_factor, ae_rr, note)
#   ddi: list of (other_drug_class, modifier_description)

DRUG_MODELS: Dict[str, Dict] = {
    "Vyvanse (lisdexamfetamine)": {
        "primary_gene": "CYP2D6",
        "gene_pheno": {
            "CYP2D6": {
                "PM": {"clearance": 60, "dose_factor": 0.7,
                       "ae_rr": 1.4, "note": "Slower amphetamine clearance; reduce starting dose."},
                "IM": {"clearance": 80, "dose_factor": 0.85,
                       "ae_rr": 1.2, "note": "Modestly slower clearance."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
                "UM": {"clearance": 130, "dose_factor": 1.2,
                       "ae_rr": 0.9, "note": "Faster clearance; higher dose may be needed."},
            },
        },
        "ddi": [
            ("Strong CYP2D6 inhibitors (paroxetine, fluoxetine, bupropion)",
             "Effectively converts intermediate → poor metabolism; reduce dose 30-50% or avoid."),
            ("MAOIs", "AVOID — serotonin syndrome / hypertensive crisis risk regardless of genotype."),
        ],
    },
    "Warfarin": {
        "primary_gene": "VKORC1",
        "gene_pheno": {
            "VKORC1": {
                "PM": {"clearance": 100, "dose_factor": 0.50,
                       "ae_rr": 3.0, "note": "VKORC1 A/A — high sensitivity; reduce dose ~50%."},
                "IM": {"clearance": 100, "dose_factor": 0.70,
                       "ae_rr": 1.8, "note": "VKORC1 G/A — moderate sensitivity."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard sensitivity."},
            },
            "CYP2C9": {
                "PM": {"clearance": 50, "dose_factor": 0.50,
                       "ae_rr": 2.5, "note": "Slow metabolism; reduce dose."},
                "IM": {"clearance": 75, "dose_factor": 0.75,
                       "ae_rr": 1.6, "note": "Intermediate."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard."},
            },
        },
        "ddi": [
            ("Amiodarone, fluconazole, metronidazole",
             "Inhibit CYP2C9 — compound with genetic slow metabolism. Dose reduction needed."),
            ("Antibiotics generally",
             "Disrupt gut vitamin K — INR fluctuations during therapy."),
        ],
    },
    "Clopidogrel": {
        "primary_gene": "CYP2C19",
        "gene_pheno": {
            "CYP2C19": {
                "PM": {"clearance": 30, "dose_factor": 0.0,
                       "ae_rr": 3.0, "note": "AVOID — minimal antiplatelet effect, "
                                              "high stent thrombosis risk. Use prasugrel or ticagrelor."},
                "IM": {"clearance": 60, "dose_factor": 0.0,
                       "ae_rr": 2.0, "note": "Reduced effect; consider alternatives."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
                "RM": {"clearance": 130, "dose_factor": 1.0,
                       "ae_rr": 1.2, "note": "Slightly enhanced antiplatelet effect; bleed risk slightly up."},
                "UM": {"clearance": 150, "dose_factor": 1.0,
                       "ae_rr": 1.3, "note": "Enhanced antiplatelet effect."},
            },
        },
        "ddi": [
            ("Omeprazole, esomeprazole",
             "Inhibit CYP2C19 — compound with PM/IM status. Use pantoprazole instead "
             "if a PPI is needed."),
        ],
    },
    "Sertraline / Citalopram / Escitalopram (SSRIs)": {
        "primary_gene": "CYP2C19",
        "gene_pheno": {
            "CYP2C19": {
                "PM": {"clearance": 50, "dose_factor": 0.5,
                       "ae_rr": 1.6, "note": "Higher plasma levels; reduce starting dose."},
                "IM": {"clearance": 75, "dose_factor": 0.75,
                       "ae_rr": 1.3, "note": "Modestly elevated levels."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
                "RM": {"clearance": 130, "dose_factor": 1.2,
                       "ae_rr": 0.9, "note": "Lower levels — may need higher dose for efficacy."},
                "UM": {"clearance": 160, "dose_factor": 1.3,
                       "ae_rr": 0.8, "note": "May underrespond — consider non-2C19 alternative (sertraline → fluoxetine)."},
            },
        },
        "ddi": [
            ("Strong CYP2C19 inhibitors", "Omeprazole, fluvoxamine — compound."),
            ("Tramadol/codeine + SSRIs", "Serotonin syndrome risk; reduce or avoid combination."),
        ],
    },
    "Simvastatin / Atorvastatin (statins)": {
        "primary_gene": "SLCO1B1",
        "gene_pheno": {
            "SLCO1B1": {
                "PM": {"clearance": 50, "dose_factor": 0.5,
                       "ae_rr": 4.5, "note": "Substantial myopathy risk at 80mg; cap simvastatin at 20mg or "
                                              "switch to rosuvastatin/pravastatin."},
                "IM": {"clearance": 75, "dose_factor": 0.75,
                       "ae_rr": 2.5, "note": "Cap simvastatin at 40mg; monitor for muscle symptoms."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
            },
        },
        "ddi": [
            ("Macrolides, azole antifungals, cyclosporine",
             "CYP3A4 / OATP inhibitors — markedly elevate exposure. Hold statin or switch."),
            ("Grapefruit juice",
             "CYP3A4 inhibition — avoid with simvastatin/atorvastatin."),
        ],
    },
    "Codeine / Tramadol (opioids)": {
        "primary_gene": "CYP2D6",
        "gene_pheno": {
            "CYP2D6": {
                "PM": {"clearance": 10, "dose_factor": 0.0,
                       "ae_rr": 1.0, "note": "AVOID — pro-drug not activated. No analgesia. Use morphine or non-opioid."},
                "IM": {"clearance": 50, "dose_factor": 0.7,
                       "ae_rr": 1.0, "note": "Reduced analgesia."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
                "UM": {"clearance": 200, "dose_factor": 0.5,
                       "ae_rr": 4.0, "note": "AVOID OR REDUCE — supratherapeutic morphine levels; "
                                              "respiratory depression risk, lethal in children."},
            },
        },
        "ddi": [
            ("Strong CYP2D6 inhibitors", "Effectively turn UM → NM, NM → IM, IM → PM. Major efficacy shifts."),
        ],
    },
    "Tamoxifen (endocrine therapy)": {
        "primary_gene": "CYP2D6",
        "gene_pheno": {
            "CYP2D6": {
                "PM": {"clearance": 40, "dose_factor": 0.0,
                       "ae_rr": 1.0, "note": "Major concern — tamoxifen activates to endoxifen via CYP2D6. "
                                              "PMs may have reduced efficacy. Discuss aromatase inhibitor alternative with oncologist."},
                "IM": {"clearance": 70, "dose_factor": 0.0,
                       "ae_rr": 1.0, "note": "Possible reduced efficacy. Discuss with oncology."},
                "NM": {"clearance": 100, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing."},
                "UM": {"clearance": 130, "dose_factor": 1.0,
                       "ae_rr": 1.0, "note": "Standard dosing; expected efficacy."},
            },
        },
        "ddi": [
            ("Strong CYP2D6 inhibitors (paroxetine, fluoxetine)",
             "AVOID — converts NM → PM efficacy. Major implication for breast cancer treatment."),
        ],
    },
}


def simulate_pgx(pgx_result: Optional[Dict]) -> Dict:
    """Layer quantitative simulation over PGx phenotype calls."""
    if not pgx_result or not pgx_result.get("per_gene"):
        return {"available": False, "drugs": []}

    pheno_by_gene = {g: r.get("phenotype_code") for g, r in pgx_result["per_gene"].items()}
    drug_sims: List[Dict] = []
    for drug, model in DRUG_MODELS.items():
        relevant_phenos: List[Dict] = []
        for gene, pheno_map in model["gene_pheno"].items():
            pheno_code = pheno_by_gene.get(gene)
            if not pheno_code:
                continue
            pheno_data = pheno_map.get(pheno_code)
            if not pheno_data:
                continue
            relevant_phenos.append({"gene": gene, "phenotype_code": pheno_code,
                                    **pheno_data})
        if relevant_phenos:
            # Combined effect (multiplicative for multiple-gene drugs like warfarin)
            combined_clearance = 100
            combined_dose_factor = 1.0
            combined_ae = 1.0
            for p in relevant_phenos:
                combined_clearance = int(combined_clearance * (p["clearance"] / 100))
                combined_dose_factor *= p["dose_factor"] if p["dose_factor"] > 0 else 0
                combined_ae *= p["ae_rr"]
            drug_sims.append({
                "drug": drug,
                "primary_gene": model["primary_gene"],
                "gene_findings": relevant_phenos,
                "combined_clearance_pct": combined_clearance,
                "combined_dose_factor": round(combined_dose_factor, 2),
                "combined_ae_rr": round(combined_ae, 2),
                "ddi_notes": model["ddi"],
            })
    return {"available": True, "drugs": drug_sims}
