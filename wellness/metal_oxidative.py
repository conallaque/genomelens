"""
Metal Handling, Oxidative Defense & Neurodegeneration Panel
===========================================================

Covers a cluster of genes the other modules did not surface as a section:

  * Neurodegeneration (Parkinson's) — LRRK2, GBA
  * Metal handling                  — MT1A, MT2A, ATP7B, HFE, SLC39A8,
                                       SLC39A14, ABCG1
  * Oxidative defense               — CAT, G6PD, GSTM1, GSTT1

Two of these SNPs (GBA rs76763715, G6PD rs1050828) are already registered
and used by carrier.py for *recessive disease* status. This module reuses
them for a different functional angle — Parkinson's risk (GBA) and red-cell
oxidative fragility (G6PD) — rather than re-declaring them.

The metallothionein (MT1A/MT2A) and ATP7B-K832R signals are research-grade
on consumer chips and are flagged ``confidence: low``. They are not a
substitute for clinical copper studies or heavy-metal panels. The same
``low`` caveat applies to the ZIP transporters (SLC39A8/SLC39A14) and ABCG1.

GSTM1 and GSTT1 are glutathione-S-transferase genes whose functional
variant is a whole-gene *deletion* (the "null" genotype), not a point SNP.
Consumer arrays do not genotype the deletion directly; the proxy SNPs used
here only hint at it, and a true null call requires a PCR/CNV assay — so
these are flagged ``confidence: low`` and read conservatively.

HFE (C282Y/H63D) is the exception: it is a clinically-validated
iron-overload locus reused from the carrier records, and a C282Y homozygous
call is surfaced at higher confidence.

Output shape matches wellness.py so the renderer can group by category:
  {category, trait, result, action, evidence, confidence}

This complements — does not replace — the wellness module's existing
"Oxidative Stress Defense" trait (SOD2/GPX1); the catalase/G6PD signals
here are an additional axis, not a contradiction.
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path

import pandas as pd

from core import snp_registry  # cross-check via audit_against_registry below


def _gt(snps_df: pd.DataFrame, rsid: str) -> str | None:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> int | None:
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


def _hfe_iron(snps):
    # Reuses the registered HFE C282Y / H63D SNPs (carrier.py uses them for
    # hereditary-hemochromatosis carrier status). Here: the iron-overload
    # angle for the metal-handling section.
    c282y = _dose(snps, "rs1800562", "A")   # derived A = C282Y
    h63d = _dose(snps, "rs1799945", "G")     # derived G = H63D
    if c282y is None and h63d is None:
        return None
    parts = []
    if c282y is not None:
        parts.append(f"C282Y rs1800562 A-dose {c282y}")
    if h63d is not None:
        parts.append(f"H63D rs1799945 G-dose {h63d}")
    finding = {"category": CAT_METAL,
               "trait": "HFE Iron Overload (Hemochromatosis)",
               "evidence": "; ".join(parts), "confidence": "low"}
    if (c282y or 0) >= 2:
        finding["result"] = ("Homozygous for HFE C282Y — the classical "
                  "hereditary-hemochromatosis genotype, the strongest common "
                  "genetic cause of iron overload. Penetrance is incomplete, "
                  "especially in women; many homozygotes never accumulate "
                  "clinically significant iron.")
        finding["action"] = ("Confirm with serum ferritin and transferrin "
                  "saturation; discuss with a clinician. Periodic phlebotomy "
                  "is the established treatment if iron studies are elevated.")
        finding["confidence"] = "high"
        finding["clinical_variant"] = {
            "gene": "HFE",
            "rsid": "rs1800562",
            "variant": "C282Y (p.Cys282Tyr) homozygous",
            "dose": c282y,
            "clinical_significance": "Pathogenic (when homozygous)",
            "inheritance": "autosomal recessive (incomplete penetrance)",
            "note": ("HFE C282Y homozygosity is the classical hereditary "
                     "hemochromatosis genotype. Consumer-chip call — confirm "
                     "with iron studies before any clinical action."),
        }
    elif (c282y or 0) == 1 and (h63d or 0) >= 1:
        finding["result"] = ("Compound heterozygous C282Y/H63D — a minority "
                  "develop mild iron overload. Lower risk than C282Y "
                  "homozygosity.")
        finding["action"] = ("Consider a one-time ferritin / transferrin-"
                  "saturation check; routine monitoring is usually "
                  "unnecessary unless iron studies are elevated.")
        finding["confidence"] = "moderate"
    elif (c282y or 0) >= 1 or (h63d or 0) >= 1:
        finding["result"] = ("Carries a single HFE variant (heterozygous "
                  "C282Y or H63D). Generally not associated with clinically "
                  "significant iron overload on its own.")
        finding["action"] = ("No HFE-specific action indicated. See the "
                  "Carrier Status section for the hemochromatosis angle.")
    else:
        finding["result"] = ("No HFE C282Y or H63D variant detected — typical "
                  "iron-handling genotype at these loci.")
        finding["action"] = "No HFE-specific action indicated."
    return finding


def _zip_transporters(snps):
    # SLC39A8 (ZIP8) and SLC39A14 (ZIP14) — divalent-metal importers handling
    # zinc/manganese (and iron for ZIP14). Collapsed into one finding.
    a8 = _dose(snps, "rs13107325", "T")   # derived T = A391T, reduced transport
    a14 = _dose(snps, "rs896378", "C")    # common ZIP14 coding allele
    if a8 is None and a14 is None:
        return None
    parts = []
    if a8 is not None:
        parts.append(f"SLC39A8 rs13107325 T-dose {a8}")
    if a14 is not None:
        parts.append(f"SLC39A14 rs896378 C-dose {a14}")
    if (a8 or 0) >= 1:
        result = ("Carries the SLC39A8 (ZIP8) A391T allele, which lowers "
                  "zinc/manganese transport and is one of the most pleiotropic "
                  "common variants (linked to blood pressure, lipids, "
                  "neuropsychiatric traits and IBD in large GWAS). Effect per "
                  "trait is small. SLC39A14 (ZIP14) variation co-reported.")
    else:
        result = ("No reduced-transport SLC39A8 A391T allele detected; "
                  "ZIP8/ZIP14 zinc/manganese transport genotype is typical "
                  "for the variants tested.")
    action = ("Research-grade only — not a clinical zinc/manganese test. "
              "Ensure adequate but not excessive dietary zinc and manganese; "
              "do not megadose. Discuss any specific metabolic concern with a "
              "clinician.")
    return {"category": CAT_METAL,
            "trait": "ZIP Zinc/Manganese Transport (SLC39A8/SLC39A14)",
            "result": result, "action": action,
            "evidence": "; ".join(parts), "confidence": "low"}


def _abcg1_efflux(snps):
    dose = _dose(snps, "rs1893590", "C")  # derived C = -204A>C, lower HDL
    if dose is None:
        return None
    if dose >= 1:
        result = ("Carries the ABCG1 −204A>C allele, associated in candidate-"
                  "gene studies with reduced cholesterol/sterol efflux and "
                  "lower HDL-C. ABCG1 moves cholesterol and oxysterols out of "
                  "cells, complementing the lipid modules. Research-grade.")
    else:
        result = ("No ABCG1 −204A>C allele detected — typical cholesterol-"
                  "efflux genotype at this locus.")
    action = ("No action from this SNP alone. See the lipid/cardiometabolic "
              "sections for HDL and cholesterol guidance; standard heart-"
              "healthy measures apply regardless of this variant.")
    return {"category": CAT_METAL, "trait": "ABCG1 Cholesterol Efflux",
            "result": result, "action": action,
            "evidence": f"rs1893590 C-allele dose: {dose}", "confidence": "low"}


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


def _gst_detox(snps):
    # GSTM1 / GSTT1 glutathione-S-transferases — phase-II conjugation enzymes
    # central to xenobiotic detox and oxidative defense. The functional null
    # is a whole-gene DELETION, which consumer arrays don't call directly;
    # these proxy SNPs only hint at it (a no-call can mean the deletion OR
    # simply missing data). Strictly research-grade.
    m1 = _dose(snps, "rs4147565", "A")   # within-GSTM1 proxy marker
    t1 = _dose(snps, "rs4630", "A")       # within-GSTT1 proxy marker
    if m1 is None and t1 is None:
        return None
    parts = []
    if m1 is not None:
        parts.append(f"GSTM1 rs4147565 A-dose {m1}")
    if t1 is not None:
        parts.append(f"GSTT1 rs4630 A-dose {t1}")
    result = ("Glutathione-S-transferase (GSTM1/GSTT1) proxy genotypes. These "
              "enzymes conjugate glutathione to detoxify carcinogens, drugs "
              "and oxidative by-products. The clinically relevant 'null' "
              "(absent-enzyme) state is a whole-gene deletion that these "
              "consumer-chip SNPs cannot call directly — treat as a weak hint "
              "only, never as a confirmed null.")
    action = ("Research-grade only — do NOT treat as a detox-capacity or "
              "cancer-risk test. A true GSTM1/GSTT1-null determination needs a "
              "PCR/CNV assay. General prudence: avoid tobacco smoke and other "
              "carcinogen exposures, and eat cruciferous vegetables (induce "
              "phase-II detox enzymes).")
    return {"category": CAT_OXID,
            "trait": "Glutathione-S-Transferase Detox (GSTM1/GSTT1)",
            "result": result, "action": action,
            "evidence": "; ".join(parts), "confidence": "low"}


# ─── Master analyzer ──────────────────────────────────────────────────────

def analyze_metal_oxidative(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        # Neurodegeneration
        _lrrk2_g2019s, _gba_parkinsons,
        # Metal handling
        _atp7b_copper, _mt_heavy_metal, _hfe_iron,
        _zip_transporters, _abcg1_efflux,
        # Oxidative defense
        _cat_catalase, _g6pd_oxidative, _gst_detox,
    ]
    predictions: list[dict] = []
    for a in analyzers:
        try:
            r = a(snps_df)
            if r is not None:
                predictions.append(r)
        except Exception:
            continue

    by_category: dict[str, list[dict]] = {}
    for p in predictions:
        by_category.setdefault(p["category"], []).append(p)

    return {
        "predictions": predictions,
        "by_category": by_category,
        "n_predictions": len(predictions),
        "categories": list(by_category.keys()),
    }


# ── Cross-check against the unified SNP registry ──────────────────────────

def _scan_rsids_referenced() -> list[str]:
    """Every rsID literal referenced in this module must be registered."""
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> dict[str, list[str]]:
    """Presence-only audit: every rsID referenced here must be in the
    registry. Returns ``{"registered": [...], "missing": [...]}``."""
    registered, missing = [], []
    for rsid in _scan_rsids_referenced():
        (registered if snp_registry.get(rsid) is not None else missing).append(rsid)
    return {"registered": registered, "missing": missing}
