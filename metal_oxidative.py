"""
Metal Handling, Oxidative Defense & Neurodegeneration Panel
===========================================================

Covers a cluster of genes the other modules did not surface as a section:

  * Neurodegeneration (Parkinson's) — LRRK2, GBA
  * Metal handling                  — MT1A, MT2A, ATP7B
  * Oxidative defense               — CAT, G6PD

Two of these SNPs (GBA rs76763715, G6PD rs1050828) are already registered
and used by carrier.py for *recessive disease* status. This module reuses
them for a different functional angle — Parkinson's risk (GBA) and red-cell
oxidative fragility (G6PD) — rather than re-declaring them.

The metallothionein (MT1A/MT2A) and ATP7B-K832R signals are research-grade
on consumer chips and are flagged ``confidence: low``. They are not a
substitute for clinical copper studies or heavy-metal panels.

Output shape matches wellness.py so the renderer can group by category:
  {category, trait, result, action, evidence, confidence}

This complements — does not replace — the wellness module's existing
"Oxidative Stress Defense" trait (SOD2/GPX1); the catalase/G6PD signals
here are an additional axis, not a contradiction.
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path
from typing import Dict, List, Optional
import pandas as pd

import snp_registry  # cross-check via audit_against_registry below


def _gt(snps_df: pd.DataFrame, rsid: str) -> Optional[str]:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> Optional[int]:
    gt = _gt(snps_df, rsid)
    if gt is None or len(gt) != 2:
        return None
    return gt.count(allele.upper())


CAT_PD = "Neurodegeneration (Parkinson's)"
CAT_METAL = "Metal Handling"
CAT_OXID = "Oxidative Defense"


# ─── Neurodegeneration (Parkinson's) ──────────────────────────────────────

def _lrrk2_g2019s(snps):
    dose = _dose(snps, "rs34637584", "A")
    if dose is None:
        return None
    finding = {"category": CAT_PD, "trait": "LRRK2 G2019S (Parkinson's)",
               "evidence": f"rs34637584 A-allele dose: {dose}",
               "confidence": "moderate"}
    if dose >= 1:
        finding["result"] = ("Carries the LRRK2 G2019S risk allele — the most "
                  "common dominant Parkinson's variant. Penetrance is "
                  "incomplete (~28-74% by age 80); many carriers never "
                  "develop PD.")
        finding["action"] = ("Confirm with clinical-grade sequencing before "
                  "acting. Discuss with a neurologist/genetic counselor; "
                  "consider exercise and avoiding pesticide exposure "
                  "(modifiable risk).")
        # Clinically-reportable dominant variant — flagged for the FHIR
        # export. Only attached when the risk allele is actually present, so
        # negative calls and wellness-tier findings are never exported.
        finding["clinical_variant"] = {
            "gene": "LRRK2",
            "rsid": "rs34637584",
            "variant": "G2019S (p.Gly2019Ser)",
            "dose": dose,
            "clinical_significance": "Pathogenic",
            "inheritance": "autosomal dominant (incomplete penetrance)",
            "note": ("LRRK2 G2019S is a ClinVar-pathogenic, autosomal-dominant "
                     "Parkinson's risk variant with incomplete penetrance. "
                     "Consumer-chip call — confirm with clinical-grade "
                     "sequencing before any clinical action."),
        }
    else:
        finding["result"] = "No LRRK2 G2019S risk allele detected at this position."
        finding["action"] = "No LRRK2-specific action indicated."
    return finding


def _gba_parkinsons(snps):
    # Reuses the registered GBA N370S SNP (carrier.py uses it for Gaucher
    # status). C = derived risk allele; here interpreted for PD risk.
    dose = _dose(snps, "rs76763715", "C")
    if dose is None:
        return None
    if dose >= 1:
        result = ("Carries GBA N370S — beyond Gaucher carrier status, GBA "
                  "variants raise Parkinson's risk roughly 5-fold (still a "
                  "low absolute lifetime risk). Often associated with earlier "
                  "onset if PD does occur.")
        action = ("Awareness only; no validated preventive therapy. Regular "
                  "exercise is the best-evidenced modifiable factor. See the "
                  "Carrier Status section for the Gaucher-disease angle.")
        conf = "moderate"
    else:
        result = "No GBA N370S allele detected (this single SNP only)."
        action = "No GBA-specific action indicated from this variant."
        conf = "moderate"
    return {"category": CAT_PD, "trait": "GBA N370S (Parkinson's risk)",
            "result": result, "action": action,
            "evidence": f"rs76763715 C-allele dose: {dose}", "confidence": conf}


# ─── Metal handling ───────────────────────────────────────────────────────

def _atp7b_copper(snps):
    dose = _dose(snps, "rs1061472", "G")
    if dose is None:
        return None
    result = ("ATP7B K832R is a common copper-transporter polymorphism. It "
              "is NOT the Wilson's disease diagnostic variant (H1069Q) — "
              "this SNP only flags minor variation in copper handling.")
    action = ("No action from this SNP alone. If Wilson's disease is "
              "suspected clinically, request serum ceruloplasmin, 24h urinary "
              "copper, and ATP7B sequencing.")
    return {"category": CAT_METAL, "trait": "ATP7B Copper Transport (K832R)",
            "result": result, "action": action,
            "evidence": f"rs1061472 G-allele dose: {dose}", "confidence": "low"}


def _mt_heavy_metal(snps):
    mt1a = _dose(snps, "rs8052394", "G")
    mt2a = _dose(snps, "rs28366003", "G")
    if mt1a is None and mt2a is None:
        return None
    score = (mt1a or 0) + (mt2a or 0)
    if score >= 2:
        result = ("Metallothionein variants (MT1A/MT2A) that some cohorts "
                  "associate with higher cadmium/lead retention and lower "
                  "heavy-metal buffering capacity. Research-grade signal.")
    else:
        result = ("Metallothionein genotype not associated with reduced "
                  "heavy-metal buffering in the variants tested.")
    action = ("Research-grade only — do not treat as a heavy-metal toxicity "
              "test. General prudence: limit known cadmium/lead exposure "
              "(smoking, contaminated water); ensure adequate dietary zinc, "
              "which induces metallothionein.")
    parts = []
    if mt1a is not None:
        parts.append(f"MT1A rs8052394 G-dose {mt1a}")
    if mt2a is not None:
        parts.append(f"MT2A rs28366003 G-dose {mt2a}")
    return {"category": CAT_METAL, "trait": "Metallothionein Heavy-Metal Binding",
            "result": result, "action": action,
            "evidence": "; ".join(parts), "confidence": "low"}


# ─── Oxidative defense ────────────────────────────────────────────────────

def _cat_catalase(snps):
    dose = _dose(snps, "rs1001179", "T")
    if dose is None:
        return None
    if dose >= 2:
        result = ("CAT -262 TT — lower catalase expression, reduced "
                  "hydrogen-peroxide clearance and higher oxidative-stress "
                  "susceptibility.")
        conf = "moderate"
    elif dose == 1:
        result = ("CAT -262 CT — intermediate catalase expression.")
        conf = "moderate"
    else:
        result = ("CAT -262 CC — typical catalase expression and H2O2 "
                  "clearance.")
        conf = "moderate"
    action = ("Supports the antioxidant-rich-diet rationale (cruciferous veg "
              "→ Nrf2, adequate selenium for glutathione peroxidase). See the "
              "Wellness 'Oxidative Stress Defense' trait for the SOD2/GPX1 axis.")
    return {"category": CAT_OXID, "trait": "Catalase (CAT) Antioxidant Capacity",
            "result": result, "action": action,
            "evidence": f"rs1001179 T-allele dose: {dose}", "confidence": conf}


def _g6pd_oxidative(snps):
    # Reuses the registered G6PD V68M SNP (carrier.py uses it for the
    # recessive/X-linked deficiency call). Here: red-cell oxidative angle.
    dose = _dose(snps, "rs1050828", "T")
    if dose is None:
        return None
    if dose >= 1:
        result = ("G6PD A- deficiency allele present — red blood cells have "
                  "reduced protection against oxidative damage, raising "
                  "hemolysis risk from fava beans and oxidant drugs.")
        action = ("Avoid known oxidant triggers (fava beans, primaquine, "
                  "rasburicase, high-dose vitamin C, naphthalene). Flag G6PD "
                  "status to prescribers — see the Emergency Card.")
        conf = "high"
    else:
        result = ("No G6PD A- deficiency allele at this SNP — normal red-cell "
                  "oxidative protection from this variant.")
        action = "No G6PD-specific precautions indicated from this variant."
        conf = "high"
    return {"category": CAT_OXID, "trait": "G6PD Red-Cell Oxidative Protection",
            "result": result, "action": action,
            "evidence": f"rs1050828 T-allele dose: {dose}", "confidence": conf}


# ─── Master analyzer ──────────────────────────────────────────────────────

def analyze_metal_oxidative(snps_df: pd.DataFrame) -> Dict:
    analyzers = [
        # Neurodegeneration
        _lrrk2_g2019s, _gba_parkinsons,
        # Metal handling
        _atp7b_copper, _mt_heavy_metal,
        # Oxidative defense
        _cat_catalase, _g6pd_oxidative,
    ]
    predictions: List[Dict] = []
    for a in analyzers:
        try:
            r = a(snps_df)
            if r is not None:
                predictions.append(r)
        except Exception:
            continue

    by_category: Dict[str, List[Dict]] = {}
    for p in predictions:
        by_category.setdefault(p["category"], []).append(p)

    return {
        "predictions": predictions,
        "by_category": by_category,
        "n_predictions": len(predictions),
        "categories": list(by_category.keys()),
    }


# ── Cross-check against the unified SNP registry ──────────────────────────

def _scan_rsids_referenced() -> List[str]:
    """Every rsID literal referenced in this module must be registered."""
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> Dict[str, List[str]]:
    """Presence-only audit: every rsID referenced here must be in the
    registry. Returns ``{"registered": [...], "missing": [...]}``."""
    registered, missing = [], []
    for rsid in _scan_rsids_referenced():
        (registered if snp_registry.get(rsid) is not None else missing).append(rsid)
    return {"registered": registered, "missing": missing}
