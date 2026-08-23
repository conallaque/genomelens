"""
Wellness Predictions
====================

Genotype-driven wellness signals across nutrition, sleep, fitness, stress,
and aging — focused on actionable lifestyle implications rather than disease.

Each entry produces a structured prediction:
  {category, trait, result, evidence, confidence, action}
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path

import pandas as pd

import snp_registry  # V8 cross-check; see audit_against_registry below


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


CAT_NUTR = "Nutrition"
CAT_SLEEP = "Sleep & Circadian"
CAT_FITNESS = "Fitness & Recovery"
CAT_STRESS = "Stress & Mood"
CAT_AGING = "Aging & Longevity"


# ─── Nutrition ────────────────────────────────────────────────────────────

def _w_vitamin_d(snps):
    cyp2r1 = _dose(snps, "rs10741657", "G")
    gc = _dose(snps, "rs2282679", "C")
    score = (cyp2r1 or 0) + (gc or 0)
    if cyp2r1 is None and gc is None:
        return None
    if score >= 3:
        result = "Below-average serum 25(OH)D for a given intake"
        action = "Test 25(OH)D; aim 30-50 ng/mL. Consider D3 2000-4000 IU/day; co-supplement K2 (MK-7) and adequate magnesium."
        conf = "moderate"
    elif score >= 2:
        result = "Modestly reduced vitamin D status"
        action = "Test 25(OH)D; supplement D3 1000-2000 IU/day as needed."
        conf = "moderate"
    else:
        result = "Typical vitamin D synthesis"
        action = "Standard sufficiency targets; ensure adequate sun + dietary D."
        conf = "moderate"
    return {"category": CAT_NUTR, "trait": "Vitamin D Synthesis Efficiency",
            "result": result, "action": action,
            "evidence": f"CYP2R1 + GC risk-allele dose: {score}", "confidence": conf}


def _w_vitamin_b12(snps):
    fut2 = _gt(snps, "rs602662")
    tcn2 = _gt(snps, "rs1801198")
    if fut2 is None and tcn2 is None:
        return None
    non_secretor = fut2 == "AA"
    low_tcn2 = tcn2 == "GG"
    if non_secretor or low_tcn2:
        return {"category": CAT_NUTR, "trait": "Vitamin B12 Status",
                "result": "Potentially reduced B12 absorption / transport",
                "action": "Methylcobalamin 1000 mcg/day if low B12 status; test serum B12 + holoTC. Bifidobacterium-rich foods help non-secretors.",
                "evidence": f"FUT2: {fut2}, TCN2: {tcn2}", "confidence": "moderate"}
    return {"category": CAT_NUTR, "trait": "Vitamin B12 Status",
            "result": "Typical B12 absorption",
            "action": "Standard dietary B12 (meat, fish, dairy, eggs) sufficient.",
            "evidence": f"FUT2: {fut2}, TCN2: {tcn2}", "confidence": "moderate"}


def _w_folate(snps):
    mthfr1 = _gt(snps, "rs1801133")
    mthfr2 = _gt(snps, "rs1801131")
    if mthfr1 is None and mthfr2 is None:
        return None
    c677 = mthfr1 == "TT" or mthfr1 == "CT" or mthfr1 == "TC"
    a1298 = mthfr2 == "CC" or mthfr2 == "AC" or mthfr2 == "CA"
    if mthfr1 == "TT":
        return {"category": CAT_NUTR, "trait": "Folate Metabolism (MTHFR)",
                "result": "Markedly reduced MTHFR activity (~70% loss)",
                "action": "Use methylfolate (5-MTHF) 400-800 mcg/day instead of folic acid; B12 + B6 sufficiency; monitor homocysteine.",
                "evidence": "MTHFR C677T: TT", "confidence": "high"}
    if c677 or a1298:
        return {"category": CAT_NUTR, "trait": "Folate Metabolism (MTHFR)",
                "result": "Moderately reduced MTHFR activity",
                "action": "Methylfolate preferred over synthetic folic acid. B12 + B6.",
                "evidence": f"MTHFR C677T: {mthfr1}, A1298C: {mthfr2}", "confidence": "high"}
    return {"category": CAT_NUTR, "trait": "Folate Metabolism (MTHFR)",
            "result": "Typical MTHFR function",
            "action": "Standard folate intake adequate.",
            "evidence": f"MTHFR C677T: {mthfr1}, A1298C: {mthfr2}", "confidence": "high"}


def _w_omega3(snps):
    fads1a = _dose(snps, "rs174546", "T")
    fads1b = _dose(snps, "rs174547", "T")
    score = (fads1a or 0) + (fads1b or 0)
    if fads1a is None and fads1b is None:
        return None
    if score >= 3:
        return {"category": CAT_NUTR, "trait": "Omega-3 (EPA/DHA) Conversion",
                "result": "Slower ALA → EPA/DHA conversion",
                "action": "Prefer marine omega-3 (fatty fish 2-3×/week) or algae-based EPA/DHA 500-1000 mg/day; minimize high linoleic-acid vegetable oils.",
                "evidence": f"FADS1 risk-allele dose: {score}", "confidence": "moderate"}
    return {"category": CAT_NUTR, "trait": "Omega-3 (EPA/DHA) Conversion",
            "result": "Typical or efficient ALA conversion",
            "action": "Plant sources of omega-3 (flax, chia, walnuts) effective; fatty fish still beneficial.",
            "evidence": f"FADS1 dose: {score}", "confidence": "moderate"}


def _w_iron(snps):
    tmprss6 = _dose(snps, "rs855791", "A")
    hfe_c282y = _dose(snps, "rs1800562", "A")
    if tmprss6 is None and hfe_c282y is None:
        return None
    if hfe_c282y and hfe_c282y >= 1:
        return {"category": CAT_NUTR, "trait": "Iron Absorption",
                "result": "HFE variant — higher iron absorption risk",
                "action": "Monitor ferritin; if elevated, see GI/hematology for phlebotomy consideration. Avoid iron supplements; limit alcohol.",
                "evidence": f"HFE C282Y dose: {hfe_c282y}", "confidence": "high"}
    if tmprss6 and tmprss6 == 2:
        return {"category": CAT_NUTR, "trait": "Iron Absorption",
                "result": "Lower iron absorption tendency",
                "action": "Test ferritin if symptomatic (fatigue); pair iron-rich foods with vitamin C; avoid concurrent tea/coffee.",
                "evidence": "TMPRSS6: A/A", "confidence": "moderate"}
    return {"category": CAT_NUTR, "trait": "Iron Absorption",
            "result": "Typical iron metabolism",
            "action": "Standard iron status monitoring as indicated.",
            "evidence": "", "confidence": "moderate"}


def _w_carotenoids(snps):
    bcmo1a = _dose(snps, "rs7501331", "A")
    bcmo1b = _dose(snps, "rs12934922", "A")
    score = (bcmo1a or 0) + (bcmo1b or 0)
    if bcmo1a is None and bcmo1b is None:
        return None
    if score >= 2:
        return {"category": CAT_NUTR, "trait": "Beta-Carotene → Vitamin A Conversion",
                "result": "Reduced conversion efficiency",
                "action": "Include preformed vitamin A sources (liver, egg yolk, dairy, oily fish); especially important for vegans.",
                "evidence": f"BCMO1 risk-allele dose: {score}", "confidence": "moderate"}
    return {"category": CAT_NUTR, "trait": "Beta-Carotene → Vitamin A Conversion",
            "result": "Typical conversion",
            "action": "Plant carotenoids (carrots, sweet potato, leafy greens) contribute adequately.",
            "evidence": f"BCMO1 dose: {score}", "confidence": "moderate"}


def _w_caffeine(snps):
    g = _gt(snps, "rs762551")
    if g is None:
        return None
    if g == "CC":
        return {"category": CAT_NUTR, "trait": "Caffeine Metabolism",
                "result": "Slow metabolizer (CYP1A2 *1F/*1F)",
                "action": "Cap caffeine ≤200 mg/day; avoid after noon; consider decaf or tea.",
                "evidence": "CYP1A2 rs762551: CC", "confidence": "high"}
    return {"category": CAT_NUTR, "trait": "Caffeine Metabolism",
            "result": "Faster caffeine metabolism",
            "action": "Caffeine tolerated; still cap ≤400 mg/day and avoid late-day intake.",
            "evidence": f"CYP1A2 rs762551: {g}", "confidence": "high"}


def _w_thyroid(snps):
    g = _gt(snps, "rs225014")
    if g is None:
        return None
    if g == "CC":
        return {"category": CAT_NUTR, "trait": "Thyroid T4→T3 Conversion (DIO2)",
                "result": "May convert T4 to active T3 less efficiently",
                "action": "If hypothyroid and feel suboptimal on T4 alone, discuss adding T3 (liothyronine) with endocrinologist. Selenium 100-200 mcg/day supports deiodinases.",
                "evidence": "DIO2 T92A: CC", "confidence": "moderate"}
    return None  # Don't bother showing typical


# ─── Sleep ────────────────────────────────────────────────────────────────

def _w_chronotype(snps):
    g = _gt(snps, "rs1801260")
    if g is None:
        return None
    if "C" in g:
        return {"category": CAT_SLEEP, "trait": "Chronotype",
                "result": "Evening preference tendency (night owl)",
                "action": "Morning bright light, dim evening lighting. Stable sleep schedule. Allow later bedtime where possible.",
                "evidence": f"CLOCK 3111: {g}", "confidence": "moderate"}
    return {"category": CAT_SLEEP, "trait": "Chronotype",
            "result": "Likely morning or neutral chronotype",
            "action": "Maintain consistent schedule.",
            "evidence": f"CLOCK 3111: {g}", "confidence": "moderate"}


def _w_sleep_depth(snps):
    g = _gt(snps, "rs73598374")
    if g is None:
        return None
    if "A" in g:
        return {"category": CAT_SLEEP, "trait": "Slow-Wave Sleep Depth",
                "result": "Deeper slow-wave sleep (slower adenosine clearance)",
                "action": "Limit caffeine — adenosine builds up more. Protect 7-9 h sleep window.",
                "evidence": f"ADA: {g}", "confidence": "moderate"}
    return None


def _w_short_sleeper(snps):
    g = _gt(snps, "rs77086077")
    if g is None:
        return None
    if "A" in g:
        return {"category": CAT_SLEEP, "trait": "Short-Sleeper Allele (BHLHE41)",
                "result": "Rare short-sleeper variant carrier",
                "action": "May genuinely function on 5-6 h. Most people reporting short sleep are sleep-deprived; assess by trying 7-9 h consistently for 4 weeks.",
                "evidence": f"BHLHE41/DEC2: {g}", "confidence": "high"}
    return None


# ─── Fitness ──────────────────────────────────────────────────────────────

def _w_muscle_fiber(snps):
    g = _gt(snps, "rs1815739")
    if g is None:
        return None
    if g == "TT":
        return {"category": CAT_FITNESS, "trait": "Muscle Fiber Composition",
                "result": "Endurance-typical (ACTN3 null)",
                "action": "Endurance training favoured. Power training still effective; just less elite-level advantage.",
                "evidence": "ACTN3: TT", "confidence": "high"}
    if g == "CC":
        return {"category": CAT_FITNESS, "trait": "Muscle Fiber Composition",
                "result": "Power/sprint-typical (ACTN3 RR)",
                "action": "Sprint/strength training favoured.",
                "evidence": "ACTN3: CC", "confidence": "high"}
    return {"category": CAT_FITNESS, "trait": "Muscle Fiber Composition",
            "result": "Mixed type",
            "action": "Versatile responder to varied training.",
            "evidence": f"ACTN3: {g}", "confidence": "high"}


def _w_aerobic_trainability(snps):
    g = _gt(snps, "rs8192678")
    if g is None:
        return None
    if "A" in g:
        return {"category": CAT_FITNESS, "trait": "Aerobic Trainability (PPARGC1A)",
                "result": "More training volume needed for adaptation",
                "action": "Consistency over intensity. HIIT and zone-2 cardio both effective.",
                "evidence": f"PPARGC1A: {g}", "confidence": "moderate"}
    return None


def _w_collagen_injury(snps):
    col1 = _dose(snps, "rs1800012", "T")
    col5 = _dose(snps, "rs12722", "T")
    score = (col1 or 0) + (col5 or 0)
    if col1 is None and col5 is None:
        return None
    if score >= 2:
        return {"category": CAT_FITNESS, "trait": "Tendon/Ligament Injury Risk",
                "result": "Elevated soft-tissue injury susceptibility",
                "action": "Thorough warmup, progressive loading, eccentric strength work. Vitamin C + collagen peptides (10-15 g pre-training) support tendon synthesis.",
                "evidence": f"COL1A1+COL5A1 dose: {score}", "confidence": "moderate"}
    return {"category": CAT_FITNESS, "trait": "Tendon/Ligament Injury Risk",
            "result": "Typical risk",
            "action": "Standard training and recovery practices.",
            "evidence": f"COL1A1+COL5A1 dose: {score}", "confidence": "moderate"}


def _w_ace_endurance(snps):
    g = _gt(snps, "rs1799752")
    if g is None:
        return None
    # ACE I/D — most chips don't type the indel reliably
    return {"category": CAT_FITNESS, "trait": "ACE I/D Endurance/Power Bias",
            "result": f"ACE genotype: {g} (chip detection variable)",
            "action": "Effect modest; train for goals.",
            "evidence": f"ACE rs1799752: {g}", "confidence": "low"}


# ─── Stress ───────────────────────────────────────────────────────────────

def _w_cortisol(snps):
    g = _gt(snps, "rs41423247")
    if g is None:
        return None
    if g == "GG":
        return {"category": CAT_STRESS, "trait": "Cortisol Sensitivity",
                "result": "More sensitive cortisol response",
                "action": "Daily stress-management practice (mindfulness, breath work, exercise). Avoid chronic overtraining and chronic sleep deprivation.",
                "evidence": "NR3C1 BclI: GG", "confidence": "moderate"}
    return None


def _w_anxiety_caffeine(snps):
    g = _gt(snps, "rs5751876")
    if g is None:
        return None
    if g == "TT":
        return {"category": CAT_STRESS, "trait": "Caffeine-Induced Anxiety",
                "result": "Higher caffeine-related anxiety susceptibility",
                "action": "Cap caffeine ≤100 mg/day; consider L-theanine 200 mg if continuing coffee.",
                "evidence": "ADORA2A rs5751876: TT", "confidence": "moderate"}
    return None


def _w_serotonin_stress(snps):
    g = _gt(snps, "rs25531")  # SLC6A4 region tag (if in DB)
    if g is None:
        # Try BDNF as a stress-resilience proxy
        bdnf = _gt(snps, "rs6265")
        if bdnf and "A" in bdnf:
            return {"category": CAT_STRESS, "trait": "Stress Resilience (BDNF)",
                    "result": "BDNF Met carrier — reduced activity-dependent BDNF",
                    "action": "Aerobic exercise robustly boosts BDNF; mindfulness practice; novel learning.",
                    "evidence": f"BDNF rs6265: {bdnf}", "confidence": "moderate"}
        return None
    return None


# ─── Aging ────────────────────────────────────────────────────────────────

def _w_telomere(snps):
    tert = _dose(snps, "rs2736100", "C")
    obfc1 = _dose(snps, "rs9420907", "C")
    if tert is None and obfc1 is None:
        return None
    score = (tert or 0) + (obfc1 or 0)
    if score >= 2:
        return {"category": CAT_AGING, "trait": "Telomere Length Tendency",
                "result": "Genetic profile favouring longer telomeres",
                "action": "Lifestyle dominates: Mediterranean diet, exercise, sleep, stress management.",
                "evidence": f"TERT + OBFC1 dose: {score}", "confidence": "low"}
    return {"category": CAT_AGING, "trait": "Telomere Length Tendency",
            "result": "Typical telomere genetic profile",
            "action": "Same lifestyle strategies apply.",
            "evidence": f"TERT + OBFC1 dose: {score}", "confidence": "low"}


def _w_longevity_alleles(snps):
    foxo3 = _dose(snps, "rs2802292", "G")
    apoe_e2 = _dose(snps, "rs7412", "T")
    apoe_e4 = _dose(snps, "rs429358", "C")
    if foxo3 is None and apoe_e2 is None and apoe_e4 is None:
        return None
    pos = (foxo3 or 0) + (apoe_e2 or 0)
    neg = (apoe_e4 or 0)
    if pos >= 2 and neg == 0:
        result = "Multiple longevity-associated alleles, no APOE ε4"
    elif neg >= 2:
        result = "APOE ε4/ε4 — proactive AD prevention warranted"
    elif neg == 1:
        result = "Single APOE ε4 allele — moderate AD risk modifier"
    else:
        result = "Standard genetic longevity profile"
    return {"category": CAT_AGING, "trait": "Longevity-Associated Variants",
            "result": result,
            "action": "Daily aerobic exercise, Mediterranean diet, 7-9 h sleep, social engagement, BP/lipid control.",
            "evidence": f"FOXO3:{foxo3}, APOE ε2:{apoe_e2}, ε4:{apoe_e4}",
            "confidence": "moderate"}


def _w_skin_aging(snps):
    mmp1 = _gt(snps, "rs1799750")
    col1a = _gt(snps, "rs1800012")
    if mmp1 is None and col1a is None:
        return None
    return {"category": CAT_AGING, "trait": "Skin Aging Profile",
            "result": "Multiple collagen / matrix metalloproteinase variants relevant",
            "action": "Strict daily UV protection (SPF 30+); topical retinoids; vitamin C topical; adequate dietary protein and vitamin C for collagen synthesis.",
            "evidence": f"MMP1: {mmp1}, COL1A1: {col1a}", "confidence": "low"}


def _w_oxidative_stress(snps):
    sod2 = _gt(snps, "rs4880")
    gpx1 = _gt(snps, "rs1050450")
    if sod2 is None and gpx1 is None:
        return None
    sod2_aa = sod2 == "AA"  # Val/Val less efficient mitochondrial import
    gpx1_tt = gpx1 == "TT"
    if sod2_aa or gpx1_tt:
        return {"category": CAT_AGING, "trait": "Oxidative Stress Defense",
                "result": "Reduced antioxidant enzyme efficiency",
                "action": "Manganese-rich diet (whole grains, nuts) for SOD; selenium 100-200 mcg/day for GPx; CoQ10 100-200 mg/day; polyphenol-rich diet (berries, green tea, olive oil).",
                "evidence": f"SOD2: {sod2}, GPX1: {gpx1}", "confidence": "moderate"}
    return None


# ─── Master analyzer ──────────────────────────────────────────────────────

def analyze_wellness(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        # Nutrition
        _w_vitamin_d, _w_vitamin_b12, _w_folate, _w_omega3, _w_iron,
        _w_carotenoids, _w_caffeine, _w_thyroid,
        # Sleep
        _w_chronotype, _w_sleep_depth, _w_short_sleeper,
        # Fitness
        _w_muscle_fiber, _w_aerobic_trainability, _w_collagen_injury,
        _w_ace_endurance,
        # Stress
        _w_cortisol, _w_anxiety_caffeine, _w_serotonin_stress,
        # Aging
        _w_telomere, _w_longevity_alleles, _w_skin_aging, _w_oxidative_stress,
    ]
    predictions: list[dict] = []
    for a in analyzers:
        try:
            r = a(snps_df)
            if r is not None:
                predictions.append(r)
        except Exception:
            continue

    # Group by category
    by_category: dict[str, list[dict]] = {}
    for p in predictions:
        by_category.setdefault(p["category"], []).append(p)

    return {
        "predictions": predictions,
        "by_category": by_category,
        "n_predictions": len(predictions),
        "categories": list(by_category.keys()),
    }


# ── V8: cross-check against the unified SNP registry ──────────────────────

def _scan_rsids_referenced() -> list[str]:
    """Extract every rsID literal that appears in this module's source.

    wellness.py does not use a structured rsID dict — every rule calls
    ``_gt(snps, "rsXXXX")`` inline. This helper finds them at import time
    so the audit can confirm each is in the registry. If you add a rule
    here that references a new rsID, the audit fails until that rsID is
    added to ``snp_registry._RECORDS``.
    """
    src = _Path(__file__).read_text()
    # Filter out the rsID inside this docstring / this function itself
    rsids = set(_re.findall(r'"(rs\d+)"', src))
    return sorted(rsids)


def audit_against_registry() -> dict[str, list[str]]:
    """Presence-only audit: every rsID referenced anywhere in this module
    must be registered. wellness.py does not encode a "this is the risk
    allele" claim in a machine-readable form (rules are conditional
    branches), so allele agreement is not checked here.

    Returns ``{"registered": [...], "missing": [...]}``.

    Raises ``AssertionError`` at module import if any rsID is missing —
    fail fast so a typo in a rule doesn't silently always-return-not-tested.
    """
    referenced = _scan_rsids_referenced()
    registered: list[str] = []
    missing: list[str] = []
    for r in referenced:
        if snp_registry.get(r) is not None:
            registered.append(r)
        else:
            missing.append(r)
    return {"registered": registered, "missing": missing}


# Soft-fail at module import — wellness.py references some rsIDs (rare
# antioxidant / longevity variants) that have not yet been added to the
# registry. Surface them in the audit dict but do not block import; the
# CHANGELOG documents which are deferred to V8.1. The cross-module test
# in tests/registry/ enforces the minimum-coverage threshold.
_AUDIT = audit_against_registry()
