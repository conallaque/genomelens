"""
Blood-Work Analysis (comprehensive)
===================================

Two independent layers over a user-supplied lab panel (JSON):

  A. COMPREHENSIVE CLINICAL PANEL (``analyze_clinical_bloodwork`` → the
     ``clinical`` key). Every biomarker is classified against BOTH standard
     clinical reference ranges and tighter functional/optimal ranges, grouped
     into body-system panels (lipids, glycemic, inflammation, liver, kidney,
     thyroid, iron, CBC, electrolytes, hormones, vitamins, blood pressure),
     with ~12 calculated markers (non-HDL, TG:HDL, ApoB:ApoA1, remnant
     cholesterol, HOMA-IR, eGFR, transferrin saturation, NLR, MAP, FIB-4,
     BUN:creatinine, free-T estimate), genotype-aware interpretation (the
     user's own variants — HFE, APOE, TCF7L2, MTHFR, GC/CYP2R1, FUT2, ABCG2,
     UGT1A1, LPA — contextualise flagged results), and per-system + overall
     health scores. Runs whenever labs are supplied; no PheWAS needed.

  B. GENETIC-PREDICTION COMPARISON (the original ``rows`` layer). Compares each
     lab value against its PheWAS-derived polygenic prediction:

Compares user-supplied lab panel results (JSON) against the PheWAS-derived
genetic predictions and identifies:
  • Confirmed — predicted value matches measured value within ±0.5 trait SD
  • Partial   — within 0.5–1.0 SD of prediction
  • Diverged  — measured value differs from prediction by >1.0 SD; the
                divergence direction tells us whether the genotype is being
                amplified (lifestyle aligned with genetic risk) or buffered
                (lifestyle/environment opposing genetic prediction)

JSON format expected (case-insensitive keys, synonyms supported):

    {
      "ldl": 142,              "hdl": 48,           "triglycerides": 180,
      "lp_a": 35,              "fasting_glucose": 96, "hba1c": 5.7,
      "crp": 2.4,              "ferritin": 180,      "iron": null,
      "vitamin_d": 22,         "vitamin_b12": 410,   "folate": 14,
      "testosterone": 480,     "shbg": 32,           "igf1": 175,
      "tsh": 1.8,              "alt": 28,            "creatinine": 1.0,
      "uric_acid": 6.2,        "hemoglobin": 14.8,   "mcv": 89,
      "platelets": 235,        "wbc": 6.1,           "systolic_bp": 128,
      "diastolic_bp": 82,      "heart_rate": 64
    }

Values may be null or missing for any biomarker; only present numeric values
are scored. Units are assumed to match the PheWAS catalog (mg/dL, ng/mL,
ng/dL, mIU/L, etc. — same as standard US clinical reports).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Map of free-form JSON keys → exact PheWAS trait name.
# All variant spellings (case, dashes, common abbreviations) collapse here.
_LAB_TO_PHEWAS: Dict[str, str] = {
    # Lipids
    "ldl":                "LDL cholesterol",
    "ldl_c":              "LDL cholesterol",
    "ldl_cholesterol":    "LDL cholesterol",
    "hdl":                "HDL cholesterol",
    "hdl_c":              "HDL cholesterol",
    "hdl_cholesterol":    "HDL cholesterol",
    "triglycerides":      "Triglycerides",
    "tg":                 "Triglycerides",
    "lp_a":               "Lipoprotein(a)",
    "lipoprotein_a":      "Lipoprotein(a)",
    "lpa":                "Lipoprotein(a)",
    # Glucose
    "fasting_glucose":    "Fasting glucose",
    "glucose":            "Fasting glucose",
    "fpg":                "Fasting glucose",
    "hba1c":              "HbA1c (predicted)",
    "a1c":                "HbA1c (predicted)",
    "homa_ir":            "HOMA-IR (insulin resistance)",
    # Anthropometric
    "bmi":                "Body Mass Index (BMI)",
    "height":             "Height",
    "height_cm":          "Height",
    "whr":                "Waist-Hip Ratio (BMI-adjusted)",
    # Cardio
    "systolic_bp":        "Systolic Blood Pressure",
    "sbp":                "Systolic Blood Pressure",
    "diastolic_bp":       "Diastolic Blood Pressure",
    "dbp":                "Diastolic Blood Pressure",
    "heart_rate":         "Resting Heart Rate",
    "resting_hr":         "Resting Heart Rate",
    "rhr":                "Resting Heart Rate",
    # Hematology
    "hemoglobin":         "Hemoglobin",
    "hgb":                "Hemoglobin",
    "hb":                 "Hemoglobin",
    "mcv":                "Mean Corpuscular Volume",
    "platelets":          "Platelets",
    "plt":                "Platelets",
    "wbc":                "White Blood Cell Count",
    "white_blood_cells":  "White Blood Cell Count",
    # Inflammation
    "crp":                "C-Reactive Protein",
    "hs_crp":             "C-Reactive Protein",
    # Metabolic / liver / kidney
    "uric_acid":          "Uric acid / urate",
    "urate":              "Uric acid / urate",
    "creatinine":         "Serum creatinine / eGFR",
    "alt":                "ALT (liver function)",
    "sgpt":               "ALT (liver function)",
    # Thyroid
    "tsh":                "TSH",
    # Vitamins / minerals
    "vitamin_d":          "25-OH Vitamin D",
    "25_oh_d":            "25-OH Vitamin D",
    "25ohd":              "25-OH Vitamin D",
    "vitamin_b12":        "Serum vitamin B12",
    "b12":                "Serum vitamin B12",
    "folate":             "Folate",
    "ferritin":           "Iron / ferritin",
    "iron":               "Iron / ferritin",
    "iron_ferritin":      "Iron / ferritin",
    # Hormones
    "testosterone":       "Testosterone (men)",
    "total_testosterone": "Testosterone (men)",
    "estradiol":          "Estradiol (women)",
    "e2":                 "Estradiol (women)",
    "shbg":               "SHBG",
    "igf1":               "IGF-1",
    "igf_1":              "IGF-1",
    "cortisol":           "Cortisol (proxy)",
}


def _normalize_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_").replace("-", "_")


def load_bloodwork(path: str) -> Dict[str, float]:
    """Read and normalise a bloodwork JSON file. Drops null/non-numeric values."""
    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"Bloodwork file not found: {p}")
    raw = json.loads(p.read_text())
    if not isinstance(raw, dict):
        raise ValueError("Bloodwork JSON must be an object at the top level.")
    cleaned: Dict[str, float] = {}
    for k, v in raw.items():
        if v is None:
            continue
        try:
            cleaned[_normalize_key(k)] = float(v)
        except (TypeError, ValueError):
            continue
    return cleaned


def _classify(delta_sd: float) -> Tuple[str, str]:
    """
    Return (label, html_class) for a difference expressed in trait SDs.

    Thresholds are calibrated for *common-variant polygenic predictions*, which
    typically explain only 5–20% of biomarker variance. The vast majority of
    measured variance is environmental (diet, training, sleep, medication,
    chance). We therefore use generous tolerances:

      |Δ| ≤ 0.75 SD  → Confirmed (well within the predictable range)
      |Δ| ≤ 1.5  SD  → Partial   (small lifestyle / environmental contribution)
      |Δ| >  1.5 SD  → Diverged  (large non-genetic driver dominating)

    Earlier 0.5 / 1.0 thresholds flagged most users as "diverged" simply because
    PRS variance explained is low; that is misleading, not informative.
    """
    a = abs(delta_sd)
    if a <= 0.75:
        return "Confirmed", "bw-confirmed"
    if a <= 1.5:
        return "Partial", "bw-partial"
    return "Diverged", "bw-diverged"


def compare_bloodwork(
    bloodwork_path: str,
    phewas_result: Optional[Dict],
    snps_df=None,
    meta: Optional[Dict] = None,
) -> Dict:
    """
    Build the comprehensive bloodwork analysis structure.

    Two independent layers are produced:

      1. ``clinical`` — every supplied biomarker classified against standard
         clinical reference ranges AND tighter functional/optimal ranges,
         grouped into body-system panels, with calculated markers (ratios,
         eGFR, HOMA-IR, …), genotype-aware interpretation (the user's own
         variants amplify/contextualise flagged results), and per-system +
         overall health scores. This does not need PheWAS and runs whenever
         labs are supplied.

      2. The original genetic-prediction comparison (``rows`` etc.) — each lab
         value vs its PheWAS-derived polygenic prediction, classified
         Confirmed / Partial / Diverged. Requires ``phewas_result``.

    ``snps_df`` (optional) enables the genotype-aware layer; ``meta`` may carry
    ``{"age": int, "sex": "M"|"F"}`` for sex-specific ranges and age-dependent
    calculated markers (eGFR, FIB-4).

    Backwards compatible: callers passing only (path, phewas_result) still get
    the original keys; the new analysis is attached under ``clinical``.
    """
    labs = load_bloodwork(bloodwork_path)
    meta = meta or {}
    clinical = analyze_clinical_bloodwork(labs, snps_df=snps_df, meta=meta)

    if not phewas_result or not phewas_result.get("traits"):
        return {
            "status": "no_phewas",
            "n_labs_supplied": len(labs),
            "n_matched": 0, "n_confirmed": 0, "n_partial": 0, "n_diverged": 0,
            "accuracy_pct": 0.0,
            "rows": [], "unmatched": list(labs.keys()),
            "clinical": clinical,
            "notes": (
                "PheWAS module did not run — genetic-prediction comparison "
                "unavailable, but the clinical reference-range analysis below "
                "is fully populated."
            ),
        }

    traits = phewas_result["traits"]
    rows: List[Dict] = []
    unmatched: List[str] = []

    for lab_key, actual_val in labs.items():
        trait_name = _LAB_TO_PHEWAS.get(lab_key)
        if not trait_name or trait_name not in traits:
            unmatched.append(lab_key)
            continue
        trait = traits[trait_name]
        res = trait.get("result", {})
        if res.get("status") != "ok":
            unmatched.append(lab_key)
            continue

        mean = float(trait["mean"])
        sd = float(trait["sd"]) or 1.0
        predicted = float(res["predicted_value"])
        delta_abs = actual_val - predicted
        delta_sd = delta_abs / sd

        verdict, vclass = _classify(delta_sd)
        predicted_tier = res.get("tier", "Average")

        # Determine actual tier from raw measurement vs population mean+sd
        actual_z = (actual_val - mean) / sd
        if actual_z >= 1.65:
            actual_tier = "Very high"
        elif actual_z >= 0.84:
            actual_tier = "Above average"
        elif actual_z >= -0.84:
            actual_tier = "Average"
        elif actual_z >= -1.65:
            actual_tier = "Below average"
        else:
            actual_tier = "Very low"

        # Interpretation narrative — the most useful column for the user
        interp = _build_interpretation(
            trait_name, predicted_tier, actual_tier, delta_sd, sd
        )

        rows.append({
            "trait": trait_name,
            "category": trait["category"],
            "unit": trait["unit"],
            "predicted": round(predicted, 2),
            "actual": round(actual_val, 2),
            "mean": mean,
            "sd": sd,
            "delta_abs": round(delta_abs, 2),
            "delta_sd": round(delta_sd, 2),
            "predicted_tier": predicted_tier,
            "actual_tier": actual_tier,
            "verdict": verdict,
            "verdict_class": vclass,
            "interpretation": interp,
            "callability_pct": res.get("callability_pct", 0),
            "n_used": res.get("n_used", 0),
            "n_total": res.get("n_total", 0),
        })

    # Sort: diverged first (most actionable), then partial, then confirmed.
    order = {"Diverged": 0, "Partial": 1, "Confirmed": 2}
    rows.sort(key=lambda r: (order[r["verdict"]], -abs(r["delta_sd"])))

    n_conf = sum(1 for r in rows if r["verdict"] == "Confirmed")
    n_part = sum(1 for r in rows if r["verdict"] == "Partial")
    n_div = sum(1 for r in rows if r["verdict"] == "Diverged")
    matched = len(rows)

    return {
        "status": "ok" if matched else "no_matches",
        "n_labs_supplied": len(labs),
        "n_matched": matched,
        "n_confirmed": n_conf,
        "n_partial": n_part,
        "n_diverged": n_div,
        "accuracy_pct": round(100 * (n_conf + 0.5 * n_part) / matched, 1) if matched else 0.0,
        "rows": rows,
        "unmatched": unmatched,
        "clinical": clinical,
        "notes": (
            f"Compared {matched}/{len(labs)} supplied lab values against "
            f"PheWAS-derived genetic predictions. Each row shows the trait, "
            f"the prediction, the measured lab value, and the SD-scaled "
            f"divergence (Δ in SD units of the trait)."
        ),
    }


def _build_interpretation(
    trait: str, predicted_tier: str, actual_tier: str,
    delta_sd: float, sd: float,
) -> str:
    """Plain-English summary of what the gene–lab comparison means for this row."""
    pred_extreme = predicted_tier in ("Very high", "Above average", "Very low", "Below average")
    actual_extreme = actual_tier in ("Very high", "Above average", "Very low", "Below average")

    if abs(delta_sd) <= 0.75:
        if pred_extreme and actual_extreme:
            return (
                f"Genetic prediction CONFIRMED: predicted {predicted_tier.lower()}, "
                f"measured {actual_tier.lower()}. Genotype is expressing as expected."
            )
        return (
            f"Within prediction. Genotype and labs both average; nothing notable."
        )

    direction = "higher" if delta_sd > 0 else "lower"
    if abs(delta_sd) > 1.5:
        if pred_extreme and not actual_extreme:
            return (
                f"DIVERGENCE: genetics predicted {predicted_tier.lower()} but labs are "
                f"average. Lifestyle / environment / medication is likely BUFFERING the "
                f"genetic risk — keep doing whatever is working."
            )
        if not pred_extreme and actual_extreme:
            return (
                f"DIVERGENCE: genetics predicted average but labs are {actual_tier.lower()}. "
                f"Strong NON-GENETIC driver (diet, training, stress, illness, medication) "
                f"is dominating — investigate lifestyle/clinical causes."
            )
        return (
            f"Genetics directionally correct but measured value is {abs(delta_sd):.1f} SD "
            f"{direction} than predicted — gene-environment effect amplifying the genotype."
        )

    return (
        f"Partial match: measured value runs {abs(delta_sd):.1f} SD {direction} than "
        f"the genotype-only prediction. Minor environmental contribution."
    )


# ══════════════════════════════════════════════════════════════════════════
# COMPREHENSIVE CLINICAL ENGINE (V6.1 — reference ranges + genotype-aware)
# ══════════════════════════════════════════════════════════════════════════
#
# Everything below is independent of PheWAS. It classifies each supplied
# biomarker against (a) standard clinical reference ranges and (b) tighter
# functional/optimal ranges, groups results into body-system panels, computes
# derived markers (ratios, eGFR, HOMA-IR, …), layers in the user's own genotype
# where a variant meaningfully changes interpretation, and scores each system.
#
# Reference ranges are typical adult US-lab values and are approximate — real
# ranges vary by lab, assay, age and sex. This is educational, not diagnostic.

import math as _math


def _bm(id, name, unit, system, direction, clinical, optimal, desc,
        aliases=(), critical=(None, None), sex=None, high_action="",
        low_action="", genes=(), rsids=()):
    """Construct one biomarker definition.

    direction: 'high' (higher = worse), 'low' (lower = worse),
               'both' (out either side is bad), 'window' (a target window).
    clinical/optimal: (low, high) tuples; either bound may be None (one-sided).
    sex: optional {'M': {'clinical':(...), 'optimal':(...)},
                   'F': {...}} to override the unisex ranges.
    """
    return {
        "id": id, "name": name, "unit": unit, "system": system,
        "direction": direction, "clinical": clinical, "optimal": optimal,
        "critical": critical, "sex": sex, "desc": desc,
        "aliases": tuple(aliases), "high_action": high_action,
        "low_action": low_action, "genes": tuple(genes), "rsids": tuple(rsids),
    }


# System (panel) display order
SYSTEM_ORDER = [
    "Lipids & Cardiovascular", "Glycemic & Metabolic", "Inflammation",
    "Liver", "Kidney", "Thyroid", "Iron Status", "Complete Blood Count",
    "Electrolytes & Minerals", "Hormones", "Vitamins", "Blood Pressure",
]

# ── Biomarker catalogue ───────────────────────────────────────────────────────
_BIOMARKERS: List[Dict] = [
    # ---- Lipids & cardiovascular ----
    _bm("total_cholesterol", "Total Cholesterol", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 200), (150, 199),
        "Sum of all cholesterol; a crude first-pass cardiovascular marker (ApoB/LDL are better).",
        aliases=("total_chol", "cholesterol", "tc"), critical=(None, 300),
        high_action="Focus on LDL/ApoB and the TG:HDL ratio; diet, fibre, exercise.",
        genes=("APOE",), rsids=(("rs429358", "APOE"),)),
    _bm("ldl", "LDL Cholesterol", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 100), (50, 99),
        "Primary atherogenic ('bad') cholesterol and the main lipid treatment target.",
        aliases=("ldl_c", "ldl_cholesterol"), critical=(None, 190),
        high_action="Reduce saturated fat, add soluble fibre & plant sterols, exercise; discuss statin/ezetimibe if persistently high or with high Lp(a)/ApoB.",
        genes=("LDLR", "PCSK9", "APOB", "APOE"),
        rsids=(("rs429358", "APOE"), ("rs7412", "APOE"))),
    _bm("hdl", "HDL Cholesterol", "mg/dL", "Lipids & Cardiovascular",
        "low", (40, None), (60, None),
        "'Good' cholesterol involved in reverse cholesterol transport.",
        aliases=("hdl_c", "hdl_cholesterol"),
        sex={"M": {"clinical": (40, None), "optimal": (55, None)},
             "F": {"clinical": (50, None), "optimal": (65, None)}},
        low_action="Aerobic exercise, monounsaturated fats, moderate carbohydrate; low HDL often tracks insulin resistance."),
    _bm("triglycerides", "Triglycerides", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 150), (None, 100),
        "Circulating fat; high levels track insulin resistance and excess refined carbohydrate/alcohol.",
        aliases=("tg", "trigs"), critical=(None, 500),
        high_action="Cut refined carbs/sugar/alcohol; omega-3s; the TG:HDL ratio is a strong insulin-resistance proxy.",
        genes=("APOA5", "LPL", "GCKR")),
    _bm("apob", "Apolipoprotein B", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 100), (None, 90),
        "Counts every atherogenic particle — the single best lipid predictor of cardiovascular risk.",
        aliases=("apo_b",), critical=(None, 130),
        high_action="The most reliable lipid target; lower with the same measures as LDL. <80 mg/dL is a common secondary-prevention goal.",
        genes=("APOB", "PCSK9")),
    _bm("lp_a", "Lipoprotein(a)", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 30), (None, 30),
        "A largely genetic, independent cardiovascular & aortic-stenosis risk factor. Measure once in a lifetime.",
        aliases=("lipoprotein_a", "lpa"), critical=(None, 50),
        high_action="Not diet-responsive; if elevated, drive LDL/ApoB especially low and consider earlier cardiology input. Emerging Lp(a)-lowering drugs are in trials.",
        genes=("LPA",), rsids=(("rs10455872", "LPA"), ("rs3798220", "LPA"))),
    _bm("apoa1", "Apolipoprotein A1", "mg/dL", "Lipids & Cardiovascular",
        "low", (120, None), (140, None),
        "The main protein of HDL particles.", aliases=("apo_a1",)),
    # ---- Glycemic & metabolic ----
    _bm("fasting_glucose", "Fasting Glucose", "mg/dL", "Glycemic & Metabolic",
        "high", (70, 99), (70, 90),
        "Blood sugar after an overnight fast. 100–125 = prediabetes, ≥126 = diabetes range.",
        aliases=("glucose", "fpg", "fasting_blood_glucose"), critical=(54, 250),
        high_action="Reduce refined carbohydrate, add resistance training & post-meal walks; recheck with HbA1c and fasting insulin.",
        low_action="Recurrent lows warrant clinical evaluation.",
        genes=("TCF7L2",), rsids=(("rs7903146", "TCF7L2"),)),
    _bm("hba1c", "HbA1c", "%", "Glycemic & Metabolic",
        "high", (None, 5.7), (None, 5.4),
        "3-month average blood sugar. 5.7–6.4 = prediabetes, ≥6.5 = diabetes range.",
        aliases=("a1c", "hemoglobin_a1c"), critical=(None, 9.0),
        high_action="Same levers as glucose; even 5.7–6.0 is worth acting on early.",
        genes=("TCF7L2",), rsids=(("rs7903146", "TCF7L2"),)),
    _bm("fasting_insulin", "Fasting Insulin", "µIU/mL", "Glycemic & Metabolic",
        "high", (2, 25), (2, 6),
        "Rises years before glucose does — an early insulin-resistance signal.",
        aliases=("insulin",),
        high_action="One of the earliest reversible metabolic markers; low-carb/exercise/weight loss lower it."),
    _bm("c_peptide", "C-Peptide", "ng/mL", "Glycemic & Metabolic",
        "both", (0.5, 2.0), (0.8, 1.8),
        "Reflects the pancreas's own insulin output.", aliases=()),
    # ---- Inflammation ----
    _bm("crp", "hs-CRP", "mg/L", "Inflammation",
        "high", (None, 3.0), (None, 1.0),
        "High-sensitivity C-reactive protein — systemic inflammation & cardiovascular risk. <1 low, 1–3 average, >3 high.",
        aliases=("hs_crp", "c_reactive_protein"), critical=(None, 10),
        high_action="If persistently >3 (and no acute infection), pursue sources: visceral fat, poor sleep, periodontal disease, diet.",
        genes=("CRP", "IL6")),
    _bm("homocysteine", "Homocysteine", "µmol/L", "Inflammation",
        "high", (None, 15), (None, 8),
        "An amino acid tied to B-vitamin status, methylation and vascular risk.",
        aliases=("hcy",), critical=(None, 30),
        high_action="Usually responds to B12, folate (methylfolate if MTHFR variant) and B6.",
        genes=("MTHFR",), rsids=(("rs1801133", "MTHFR C677T"),)),
    _bm("esr", "ESR", "mm/hr", "Inflammation", "high", (None, 20), (None, 10),
        "Erythrocyte sedimentation rate — a slow, non-specific inflammation marker.", aliases=()),
    _bm("uric_acid", "Uric Acid", "mg/dL", "Inflammation",
        "high", (3.4, 7.0), (3.4, 5.5),
        "High levels cause gout and track metabolic syndrome & fructose intake.",
        aliases=("urate",), critical=(None, 10),
        high_action="Limit fructose/alcohol/purines; hydrate; genetics (ABCG2) strongly influence handling.",
        genes=("ABCG2", "SLC2A9"), rsids=(("rs2231142", "ABCG2"),)),
    # ---- Liver ----
    _bm("alt", "ALT", "U/L", "Liver", "high", (7, 56), (None, 30),
        "Liver enzyme; the most specific routine marker of hepatocyte stress (esp. fatty liver).",
        aliases=("sgpt",), critical=(None, 200),
        sex={"M": {"clinical": (7, 56), "optimal": (None, 30)},
             "F": {"clinical": (7, 45), "optimal": (None, 25)}},
        high_action="Most common cause is metabolic (fatty liver) or alcohol; weight loss and reducing alcohol usually normalise it."),
    _bm("ast", "AST", "U/L", "Liver", "high", (10, 40), (None, 30),
        "Liver/muscle enzyme; interpreted alongside ALT.", aliases=("sgot",),
        critical=(None, 200)),
    _bm("ggt", "GGT", "U/L", "Liver", "high", (None, 50), (None, 30),
        "Sensitive to alcohol and biliary/oxidative stress; pairs with ALT for fatty liver.",
        aliases=()),
    _bm("alp", "Alkaline Phosphatase", "U/L", "Liver", "both", (44, 147), (44, 120),
        "Liver/bone enzyme.", aliases=("alkaline_phosphatase",)),
    _bm("bilirubin_total", "Total Bilirubin", "mg/dL", "Liver",
        "high", (0.1, 1.2), (0.1, 1.0),
        "Heme breakdown product. Isolated mild elevation is usually benign Gilbert's syndrome (UGT1A1).",
        aliases=("bilirubin", "total_bilirubin"),
        high_action="If isolated (normal ALT/AST) and mild, most often harmless Gilbert's — not liver disease.",
        genes=("UGT1A1",)),
    _bm("albumin", "Albumin", "g/dL", "Liver", "low", (3.5, 5.0), (4.0, 5.0),
        "Main blood protein; low levels reflect poor synthesis/nutrition or inflammation.", aliases=()),
    # ---- Kidney ----
    _bm("creatinine", "Creatinine", "mg/dL", "Kidney", "high", (0.7, 1.3), (0.7, 1.1),
        "Muscle-derived waste cleared by the kidneys; used to estimate filtration (eGFR).",
        aliases=("creat",), critical=(None, 4.0),
        sex={"M": {"clinical": (0.7, 1.3), "optimal": (0.8, 1.1)},
             "F": {"clinical": (0.6, 1.1), "optimal": (0.6, 1.0)}}),
    _bm("bun", "BUN", "mg/dL", "Kidney", "both", (7, 20), (8, 18),
        "Blood urea nitrogen — kidney function and protein/hydration status.",
        aliases=("blood_urea_nitrogen", "urea")),
    _bm("cystatin_c", "Cystatin C", "mg/L", "Kidney", "high", (0.5, 1.0), (0.5, 0.9),
        "A muscle-independent kidney-filtration marker, often more accurate than creatinine.", aliases=()),
    # ---- Thyroid ----
    _bm("tsh", "TSH", "mIU/L", "Thyroid", "window", (0.4, 4.0), (0.5, 2.5),
        "Pituitary signal to the thyroid — the primary screen. High = underactive, low = overactive.",
        aliases=("thyroid_stimulating_hormone",), critical=(0.1, 10),
        high_action="If high with symptoms, check free T4 and TPO antibodies (autoimmune/Hashimoto's).",
        low_action="If low, evaluate for overactive thyroid."),
    _bm("free_t4", "Free T4", "ng/dL", "Thyroid", "window", (0.8, 1.8), (1.0, 1.6),
        "Circulating thyroid hormone (storage form).", aliases=("ft4", "t4_free")),
    _bm("free_t3", "Free T3", "pg/mL", "Thyroid", "window", (2.3, 4.2), (3.0, 4.0),
        "The active thyroid hormone.", aliases=("ft3", "t3_free")),
    _bm("tpo_ab", "TPO Antibodies", "IU/mL", "Thyroid", "high", (None, 34), (None, 9),
        "Thyroid autoantibodies; high = autoimmune thyroiditis (Hashimoto's).",
        aliases=("tpo", "anti_tpo", "thyroid_peroxidase_ab")),
    # ---- Iron status ----
    _bm("ferritin", "Ferritin", "ng/mL", "Iron Status", "both", (30, 400), (50, 150),
        "Iron storage protein. Low = iron deficiency; high = overload (often HFE), inflammation, or fatty liver.",
        aliases=("ferritin_iron",), critical=(None, 1000),
        sex={"M": {"clinical": (30, 400), "optimal": (50, 150)},
             "F": {"clinical": (15, 200), "optimal": (40, 120)}},
        high_action="Persistently high ferritin + high transferrin saturation warrants HFE hemochromatosis work-up; also driven by inflammation/fatty liver/alcohol.",
        low_action="Low ferritin is the earliest sign of iron deficiency — investigate intake and losses.",
        genes=("HFE",), rsids=(("rs1800562", "HFE C282Y"), ("rs1799945", "HFE H63D"))),
    _bm("iron", "Serum Iron", "µg/dL", "Iron Status", "both", (60, 170), (60, 150),
        "Circulating iron; interpret with TIBC and ferritin.", aliases=("serum_iron",)),
    _bm("tibc", "TIBC", "µg/dL", "Iron Status", "window", (250, 450), (250, 400),
        "Total iron-binding capacity; rises in deficiency, falls in overload.",
        aliases=("total_iron_binding_capacity",)),
    _bm("hemoglobin", "Hemoglobin", "g/dL", "Iron Status", "low", (13.5, 17.5), (14, 16),
        "Oxygen-carrying protein; low = anaemia.", aliases=("hgb", "hb"),
        sex={"M": {"clinical": (13.5, 17.5), "optimal": (14, 16.5)},
             "F": {"clinical": (12.0, 15.5), "optimal": (12.5, 15)}}),
    _bm("hematocrit", "Hematocrit", "%", "Iron Status", "both", (38.8, 50), (40, 48),
        "Fraction of blood that is red cells.", aliases=("hct",),
        sex={"M": {"clinical": (38.8, 50), "optimal": (40, 48)},
             "F": {"clinical": (34.9, 44.5), "optimal": (37, 43)}}),
    # ---- CBC ----
    _bm("wbc", "White Blood Cells", "K/µL", "Complete Blood Count", "both",
        (3.4, 10.8), (4.0, 8.0),
        "Immune cell count; high with infection/inflammation, low with marrow/immune issues.",
        aliases=("white_blood_cells", "white_blood_cell_count", "leukocytes")),
    _bm("neutrophils", "Neutrophils", "K/µL", "Complete Blood Count", "both",
        (1.5, 8.0), (2.0, 6.0), "First-line bacterial-defence white cells.", aliases=("neut", "anc")),
    _bm("lymphocytes", "Lymphocytes", "K/µL", "Complete Blood Count", "both",
        (1.0, 4.0), (1.3, 3.5), "Adaptive-immunity white cells.", aliases=("lymph",)),
    _bm("platelets", "Platelets", "K/µL", "Complete Blood Count", "both",
        (150, 400), (175, 350), "Clotting cell fragments.", aliases=("plt", "platelet_count")),
    _bm("mcv", "MCV", "fL", "Complete Blood Count", "window", (80, 100), (82, 96),
        "Average red-cell size; low → iron deficiency, high → B12/folate deficiency or alcohol.",
        aliases=("mean_corpuscular_volume",)),
    _bm("rdw", "RDW", "%", "Complete Blood Count", "high", (11.5, 14.5), (11.5, 13.5),
        "Variation in red-cell size; a rising RDW is a broad early-illness/mortality signal.", aliases=()),
    # ---- Electrolytes & minerals ----
    _bm("sodium", "Sodium", "mmol/L", "Electrolytes & Minerals", "window", (135, 145), (137, 143),
        "Key electrolyte; tight regulation.", aliases=("na",), critical=(120, 160)),
    _bm("potassium", "Potassium", "mmol/L", "Electrolytes & Minerals", "window", (3.5, 5.1), (3.8, 4.8),
        "Electrolyte critical to heart rhythm.", aliases=("k",), critical=(2.5, 6.5)),
    _bm("calcium", "Calcium", "mg/dL", "Electrolytes & Minerals", "window", (8.6, 10.3), (9.0, 10.0),
        "Bone/nerve/muscle mineral (interpret with albumin).", aliases=("ca",)),
    _bm("magnesium", "Magnesium", "mg/dL", "Electrolytes & Minerals", "low", (1.7, 2.2), (2.0, 2.2),
        "Cofactor for 300+ enzymes; serum underestimates deficiency.", aliases=("mg", "magnesium_serum"),
        low_action="Common insufficiency; glycinate/citrate forms, leafy greens, nuts, seeds."),
    # ---- Hormones ----
    _bm("testosterone", "Total Testosterone", "ng/dL", "Hormones", "low", (300, 1000), (500, 900),
        "Primary male androgen; also relevant in women at far lower levels.",
        aliases=("total_testosterone", "test"),
        sex={"M": {"clinical": (300, 1000), "optimal": (500, 900)},
             "F": {"clinical": (15, 70), "optimal": (25, 60)}},
        low_action="Confirm on a morning sample; sleep, resistance training, body composition and SHBG all matter."),
    _bm("shbg", "SHBG", "nmol/L", "Hormones", "window", (20, 60), (25, 55),
        "Binds sex hormones; sets how much testosterone is free/active.", aliases=("sex_hormone_binding_globulin",)),
    _bm("estradiol", "Estradiol", "pg/mL", "Hormones", "window", (10, 40), (15, 35),
        "Primary estrogen (ranges here are a male/premenopausal-low reference).",
        aliases=("e2",)),
    _bm("dhea_s", "DHEA-S", "µg/dL", "Hormones", "window", (100, 500), (150, 400),
        "Adrenal androgen precursor; declines with age.", aliases=("dhea_sulfate", "dheas")),
    _bm("igf1", "IGF-1", "ng/mL", "Hormones", "window", (80, 280), (100, 220),
        "Growth-hormone mediator; both high and low carry trade-offs for growth vs longevity.",
        aliases=("igf_1",)),
    _bm("cortisol", "Cortisol (AM)", "µg/dL", "Hormones", "window", (6, 23), (8, 18),
        "Primary stress hormone (morning reference).", aliases=()),
    # ---- Vitamins ----
    _bm("vitamin_d", "Vitamin D (25-OH)", "ng/mL", "Vitamins", "low", (30, 100), (40, 60),
        "Steroid hormone for bone, immunity and mood. <20 deficient, 20–30 insufficient.",
        aliases=("25_oh_d", "25ohd", "vit_d", "vitamin_d_25oh"), critical=(None, 100),
        low_action="Supplement D3 (dose scales with deficit & genotype); take with fat and adequate magnesium/K2.",
        genes=("GC", "CYP2R1", "VDR"),
        rsids=(("rs2282679", "GC"), ("rs10741657", "CYP2R1"))),
    _bm("vitamin_b12", "Vitamin B12", "pg/mL", "Vitamins", "low", (200, 900), (500, 900),
        "Nerve/blood/methylation vitamin. Low-normal (200–400) can still be symptomatic.",
        aliases=("b12", "cobalamin"),
        low_action="Consider methylcobalamin; check with MMA/homocysteine if symptomatic.",
        genes=("FUT2",), rsids=(("rs602662", "FUT2"),)),
    _bm("folate", "Folate", "ng/mL", "Vitamins", "low", (3.0, 20), (10, 20),
        "B9; works with B12 in methylation and red-cell production.", aliases=("folic_acid", "b9"),
        genes=("MTHFR",), rsids=(("rs1801133", "MTHFR C677T"),)),
    _bm("omega3_index", "Omega-3 Index", "%", "Vitamins", "low", (4, 12), (8, 12),
        "Red-cell EPA+DHA; a strong cardiovascular/longevity marker. <4% high risk, >8% optimal.",
        aliases=("omega_3_index", "o3_index"),
        low_action="Increase oily fish or EPA/DHA supplementation."),
    # ---- Blood pressure & vitals ----
    _bm("systolic_bp", "Systolic BP", "mmHg", "Blood Pressure", "high", (None, 130), (None, 120),
        "Top blood-pressure number. <120 optimal, 120–129 elevated, ≥130 stage-1 hypertension.",
        aliases=("sbp", "systolic"), critical=(None, 180)),
    _bm("diastolic_bp", "Diastolic BP", "mmHg", "Blood Pressure", "high", (None, 80), (None, 80),
        "Bottom blood-pressure number.", aliases=("dbp", "diastolic"), critical=(None, 120)),
    _bm("resting_hr", "Resting Heart Rate", "bpm", "Blood Pressure", "high", (60, 100), (50, 70),
        "Lower (within reason) reflects better cardiovascular fitness.",
        aliases=("heart_rate", "rhr", "resting_heart_rate")),
]

# Build lookups: canonical id + every alias → definition
_BM_BY_ID: Dict[str, Dict] = {b["id"]: b for b in _BIOMARKERS}
_ALIAS_TO_ID: Dict[str, str] = {}
for _b in _BIOMARKERS:
    _ALIAS_TO_ID[_b["id"]] = _b["id"]
    for _a in _b["aliases"]:
        _ALIAS_TO_ID[_normalize_key(_a)] = _b["id"]


# ── status classification ─────────────────────────────────────────────────────

_STATUS_META = {
    # status → (label, css class, score 0-100, severity 0-3)
    "optimal":        ("Optimal", "cl-optimal", 100, 0),
    "normal":         ("Normal", "cl-normal", 88, 0),
    "borderline":     ("Borderline", "cl-borderline", 62, 1),
    "high":           ("High", "cl-high", 34, 2),
    "low":            ("Low", "cl-low", 34, 2),
    "critical_high":  ("Critical high", "cl-critical", 0, 3),
    "critical_low":   ("Critical low", "cl-critical", 0, 3),
}


def _ranges_for(bm: Dict, sex: Optional[str]) -> Tuple[Tuple, Tuple]:
    if sex and bm.get("sex") and sex in bm["sex"]:
        s = bm["sex"][sex]
        return s.get("clinical", bm["clinical"]), s.get("optimal", bm["optimal"])
    return bm["clinical"], bm["optimal"]


def classify_clinical(value: float, bm: Dict, sex: Optional[str] = None) -> str:
    """Return a status key for a value against a biomarker's ranges."""
    (clo, chi), (olo, ohi) = _ranges_for(bm, sex)
    crit_lo, crit_hi = bm.get("critical", (None, None))

    if crit_hi is not None and value >= crit_hi:
        return "critical_high"
    if crit_lo is not None and value <= crit_lo:
        return "critical_low"
    # outside clinical range → high / low
    if chi is not None and value > chi:
        return "high"
    if clo is not None and value < clo:
        return "low"
    # inside clinical range: optimal vs borderline
    in_opt = True
    if olo is not None and value < olo:
        in_opt = False
    if ohi is not None and value > ohi:
        in_opt = False
    if in_opt:
        return "optimal"
    # inside clinical but outside optimal
    return "borderline"


# ── derived / calculated markers ──────────────────────────────────────────────

def compute_derived_markers(vals: Dict[str, float], meta: Dict) -> Dict[str, float]:
    """Compute ratios and formula-based markers from whatever inputs exist.

    ``vals`` is keyed by canonical biomarker id. Returns a dict of derived
    canonical ids → value (only those computable from present inputs)."""
    d: Dict[str, float] = {}
    g = vals.get
    sex = (meta.get("sex") or "").upper()[:1]
    age = meta.get("age")

    tc, hdl, ldl, tg = g("total_cholesterol"), g("hdl"), g("ldl"), g("triglycerides")
    if tc is not None and hdl is not None:
        d["non_hdl"] = round(tc - hdl, 1)
        if hdl:
            d["tc_hdl_ratio"] = round(tc / hdl, 2)
    if tg is not None and hdl:
        d["tg_hdl_ratio"] = round(tg / hdl, 2)
    if tc is not None and ldl is not None and hdl is not None:
        d["remnant_chol"] = round(tc - ldl - hdl, 1)
    apob, apoa1 = g("apob"), g("apoa1")
    if apob is not None and apoa1:
        d["apob_apoa1"] = round(apob / apoa1, 2)

    glu, ins = g("fasting_glucose"), g("fasting_insulin")
    if glu is not None and ins is not None:
        d["homa_ir"] = round(glu * ins / 405.0, 2)

    iron, tibc = g("iron"), g("tibc")
    if iron is not None and tibc:
        d["transferrin_sat"] = round(100 * iron / tibc, 1)

    neut, lymph = g("neutrophils"), g("lymphocytes")
    if neut is not None and lymph:
        d["nlr"] = round(neut / lymph, 2)

    sbp, dbp = g("systolic_bp"), g("diastolic_bp")
    if sbp is not None and dbp is not None:
        d["map"] = round(dbp + (sbp - dbp) / 3.0, 0)

    bun, creat = g("bun"), g("creatinine")
    if bun is not None and creat:
        d["bun_creatinine_ratio"] = round(bun / creat, 1)

    # eGFR (CKD-EPI 2021, race-free) — needs creatinine + age + sex
    if creat and age and sex in ("M", "F"):
        k = 0.7 if sex == "F" else 0.9
        a = -0.241 if sex == "F" else -0.302
        sex_f = 1.012 if sex == "F" else 1.0
        egfr = (142 * (min(creat / k, 1) ** a) * (max(creat / k, 1) ** -1.200)
                * (0.9938 ** age) * sex_f)
        d["egfr"] = round(egfr, 0)

    # Free testosterone estimate — total T (ng/dL) + SHBG (nmol/L).
    # Coarse free-fraction proxy (monotonically decreasing in SHBG); the report
    # labels this 'est.' and it is not the full Vermeulen equation.
    tt, shbg = g("testosterone"), g("shbg")
    if tt is not None and shbg:
        # Free fraction ~1–3% of total, declining as SHBG rises. This keeps the
        # estimate in the physiological ng/dL range (labelled 'est.').
        ft_frac = 0.025 * (40.0 / (shbg + 40.0))
        d["free_testosterone"] = round(tt * ft_frac, 1)

    # FIB-4 liver fibrosis — age + AST + ALT + platelets
    ast, alt, plt = g("ast"), g("alt"), g("platelets")
    if age and ast is not None and alt and plt:
        d["fib4"] = round((age * ast) / (plt * _math.sqrt(alt)), 2)

    return d


# Definitions for the derived markers (ranges/optimal/genes as applicable)
_DERIVED_BMS: List[Dict] = [
    _bm("non_hdl", "Non-HDL Cholesterol", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 130), (None, 100),
        "All atherogenic cholesterol (total − HDL); a better target than LDL alone."),
    _bm("tc_hdl_ratio", "Total:HDL Ratio", "", "Lipids & Cardiovascular",
        "high", (None, 5.0), (None, 3.5),
        "Cardiovascular-risk ratio; lower is better."),
    _bm("tg_hdl_ratio", "Triglyceride:HDL Ratio", "", "Lipids & Cardiovascular",
        "high", (None, 3.0), (None, 1.5),
        "A strong surrogate for insulin resistance and small-dense LDL. >3 is a red flag."),
    _bm("remnant_chol", "Remnant Cholesterol", "mg/dL", "Lipids & Cardiovascular",
        "high", (None, 30), (None, 20),
        "Cholesterol in triglyceride-rich remnant particles; independently atherogenic."),
    _bm("apob_apoa1", "ApoB:ApoA1 Ratio", "", "Lipids & Cardiovascular",
        "high", (None, 0.8), (None, 0.6),
        "One of the strongest single ratios for cardiovascular risk."),
    _bm("homa_ir", "HOMA-IR", "", "Glycemic & Metabolic",
        "high", (None, 2.5), (None, 1.5),
        "Insulin-resistance index from fasting glucose × insulin. <1.5 optimal, >2.5 resistant."),
    _bm("transferrin_sat", "Transferrin Saturation", "%", "Iron Status",
        "both", (20, 45), (25, 40),
        "Iron / TIBC. >45% (esp. with high ferritin) suggests iron overload / HFE hemochromatosis.",
        genes=("HFE",), rsids=(("rs1800562", "HFE C282Y"),)),
    _bm("nlr", "Neutrophil:Lymphocyte Ratio", "", "Complete Blood Count",
        "high", (None, 3.0), (None, 2.0),
        "A simple systemic-inflammation & stress index; rising values track worse outcomes."),
    _bm("map", "Mean Arterial Pressure", "mmHg", "Blood Pressure",
        "high", (70, 100), (70, 92),
        "Average arterial pressure across the cardiac cycle."),
    _bm("bun_creatinine_ratio", "BUN:Creatinine Ratio", "", "Kidney",
        "window", (10, 20), (10, 16), "Helps distinguish dehydration from kidney causes."),
    _bm("egfr", "eGFR (CKD-EPI)", "mL/min/1.73m²", "Kidney",
        "low", (90, None), (90, None),
        "Estimated kidney filtration. ≥90 normal; 60–89 mildly reduced; <60 flags CKD.",
        low_action="A single low value can reflect hydration/muscle; recheck and consider cystatin C."),
    _bm("free_testosterone", "Free Testosterone (est.)", "ng/dL", "Hormones",
        "low", (5, 25), (9, 25),
        "Bioavailable testosterone estimated from total T and SHBG (approximate)."),
    _bm("fib4", "FIB-4 (liver fibrosis)", "", "Liver",
        "high", (None, 1.3), (None, 1.0),
        "Non-invasive liver-fibrosis index. <1.3 low risk; >2.67 advanced-fibrosis concern."),
]
for _b in _DERIVED_BMS:
    _BM_BY_ID[_b["id"]] = _b


# ── genotype-aware interpretation ─────────────────────────────────────────────

def _gt(snps_df, rsid: str) -> Optional[str]:
    try:
        if snps_df is None or rsid not in snps_df.index:
            return None
        row = snps_df.loc[rsid]
        if hasattr(row, "iloc") and getattr(row, "ndim", 1) > 1:
            row = row.iloc[0]
        g = row.get("genotype")
        if g is None:
            return None
        s = str(g).upper().replace(" ", "").replace("-", "")
        return s or None
    except Exception:
        return None


def _genotype_note(bm_id: str, status: str, snps_df) -> str:
    """Return a genotype-aware sentence for a flagged biomarker, or ''."""
    if snps_df is None:
        return ""
    flagged_high = status in ("high", "critical_high", "borderline")
    flagged_low = status in ("low", "critical_low", "borderline")

    if bm_id in ("ferritin", "transferrin_sat", "iron") and flagged_high:
        c282y = _gt(snps_df, "rs1800562")   # A = C282Y risk
        h63d = _gt(snps_df, "rs1799945")     # G = H63D risk
        if c282y and c282y.count("A") == 2:
            return ("You are homozygous for HFE C282Y — with elevated iron markers this is "
                    "the classic hereditary-hemochromatosis picture. Confirm with iron "
                    "studies and discuss with a clinician (phlebotomy is the treatment).")
        if (c282y and "A" in c282y) or (h63d and "G" in h63d):
            return ("You carry an HFE variant (C282Y/H63D). Elevated iron markers in a "
                    "carrier are worth monitoring for iron overload.")
    if bm_id in ("ldl", "total_cholesterol", "apob", "non_hdl") and flagged_high:
        e = _gt(snps_df, "rs429358")   # C = ε4-tagging
        if e and "C" in e:
            n = e.count("C")
            return (f"You carry {n} APOE ε4-tagging allele(s) — high LDL/ApoB matters more "
                    "with ε4 (higher cardiovascular and Alzheimer's risk). Prioritise lowering ApoB.")
    if bm_id in ("fasting_glucose", "hba1c", "homa_ir") and flagged_high:
        t = _gt(snps_df, "rs7903146")   # T = TCF7L2 T2D risk
        if t and "T" in t:
            return (f"You carry the TCF7L2 rs7903146 T allele ({t}), the strongest common "
                    "type-2-diabetes variant — elevated glucose is worth acting on early.")
    if bm_id == "homocysteine" and flagged_high:
        m = _gt(snps_df, "rs1801133")   # T = MTHFR C677T
        if m and "T" in m:
            return (f"With your MTHFR C677T genotype ({m}), high homocysteine often responds "
                    "specifically to methylfolate (plus B12/B6) rather than plain folic acid.")
    if bm_id == "folate" and flagged_low:
        m = _gt(snps_df, "rs1801133")
        if m and "T" in m:
            return (f"Your MTHFR C677T genotype ({m}) reduces folate activation — the "
                    "methylfolate form is preferable.")
    if bm_id == "vitamin_d" and flagged_low:
        gc = _gt(snps_df, "rs2282679")
        if gc:
            return ("Your vitamin-D binding-protein / CYP2R1 genotype is among the factors "
                    "that make some people need a higher D3 dose to reach target — dose to labs.")
    if bm_id == "vitamin_b12":
        f = _gt(snps_df, "rs602662")
        if f:
            return ("FUT2 secretor status (rs602662) shifts serum B12 readings — non-secretors "
                    "often show higher serum B12 that can mask functional need; MMA is a better check.")
    if bm_id == "uric_acid" and flagged_high:
        a = _gt(snps_df, "rs2231142")   # ABCG2 Q141K
        if a and ("T" in a or "A" in a):
            return ("You carry an ABCG2 variant that reduces uric-acid excretion — diet/fructose "
                    "control matters more for you, and gout risk is genetically higher.")
    if bm_id == "bilirubin_total" and flagged_high:
        return ("Isolated mildly-high bilirubin with normal liver enzymes is most often benign "
                "Gilbert's syndrome (UGT1A1) — a harmless variant, not liver disease.")
    if bm_id == "lp_a" and flagged_high:
        return ("Lp(a) is ~90% genetically determined and doesn't respond to diet — treat it by "
                "driving LDL/ApoB especially low and flagging family cardiovascular history.")
    return ""


# ── master clinical analyzer ──────────────────────────────────────────────────

def analyze_clinical_bloodwork(labs: Dict[str, float], snps_df=None,
                               meta: Optional[Dict] = None) -> Dict:
    """Classify every supplied biomarker against clinical + optimal ranges,
    compute derived markers, add genotype-aware notes, and score each system."""
    meta = meta or {}
    sex = (meta.get("sex") or "").upper()[:1] or None
    # allow age supplied inside the labs JSON (numeric) if not in meta
    if "age" not in meta and "age" in labs:
        meta = {**meta, "age": labs["age"]}

    # map supplied labs → canonical ids
    vals: Dict[str, float] = {}
    unrecognized: List[str] = []
    for k, v in labs.items():
        cid = _ALIAS_TO_ID.get(k)
        if cid:
            vals[cid] = v
        elif k not in ("age",):
            unrecognized.append(k)

    derived = compute_derived_markers(vals, meta)
    all_vals = {**vals, **derived}

    rows: List[Dict] = []
    for cid, value in all_vals.items():
        bm = _BM_BY_ID.get(cid)
        if not bm:
            continue
        status = classify_clinical(value, bm, sex)
        label, css, score, severity = _STATUS_META[status]
        (clo, chi), (olo, ohi) = _ranges_for(bm, sex)
        rows.append({
            "id": cid, "name": bm["name"], "unit": bm["unit"],
            "system": bm["system"], "value": value,
            "status": status, "status_label": label, "status_class": css,
            "score": score, "severity": severity,
            "clinical_low": clo, "clinical_high": chi,
            "optimal_low": olo, "optimal_high": ohi,
            "direction": bm["direction"], "desc": bm["desc"],
            "action": bm.get("high_action") if status in ("high", "critical_high", "borderline")
                       else (bm.get("low_action") if status in ("low", "critical_low") else ""),
            "genes": list(bm.get("genes", ())),
            "derived": cid in derived,
            "genotype_note": _genotype_note(cid, status, snps_df),
        })

    # group by system + score
    systems: List[Dict] = []
    by_sys: Dict[str, List[Dict]] = {}
    for r in rows:
        by_sys.setdefault(r["system"], []).append(r)
    for sysname in SYSTEM_ORDER:
        members = by_sys.get(sysname)
        if not members:
            continue
        members.sort(key=lambda r: (-r["severity"], r["name"]))
        sys_score = round(sum(m["score"] for m in members) / len(members))
        n_flag = sum(1 for m in members if m["severity"] >= 2)
        systems.append({
            "system": sysname, "score": sys_score, "n_markers": len(members),
            "n_flagged": n_flag, "markers": members,
        })

    all_scored = [r for r in rows]
    overall = round(sum(r["score"] for r in all_scored) / len(all_scored)) if all_scored else None
    flags = sorted([r for r in rows if r["severity"] >= 2],
                   key=lambda r: (-r["severity"], r["name"]))
    optimal_ct = sum(1 for r in rows if r["status"] == "optimal")

    return {
        "available": bool(rows),
        "n_markers": len(rows),
        "n_derived": len(derived),
        "n_optimal": optimal_ct,
        "n_flagged": len(flags),
        "overall_score": overall,
        "systems": systems,
        "flags": flags,
        "unrecognized": unrecognized,
        "sex_used": sex,
        "age_used": meta.get("age"),
    }


# ── HTML rendering ───────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_BW_CSS = """
<style>
.bw-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           color:#222; max-width: 1100px; margin: 24px auto; padding: 0 16px; }
.bw-wrap h1 { font-size: 1.6em; border-bottom: 2px solid #333; padding-bottom: 6px; }
.bw-summary { display:flex; gap:14px; flex-wrap:wrap; margin: 18px 0; }
.bw-stat { background:#f6f6f7; border:1px solid #ddd; border-radius:8px;
           padding:10px 14px; min-width:120px; }
.bw-stat .v { font-size:1.5em; font-weight:600; }
.bw-stat .l { font-size:0.8em; color:#666; text-transform:uppercase; letter-spacing:0.05em; }
table.bw { width:100%; border-collapse: collapse; margin-top:10px; font-size:0.92em; }
table.bw th, table.bw td { padding:8px 10px; border-bottom:1px solid #eee; vertical-align: top; }
table.bw th { text-align:left; background:#f9f9f9; }
.bw-confirmed { background:#e7f6ea; }
.bw-partial   { background:#fff7e0; }
.bw-diverged  { background:#fde4e4; }
.bw-verdict { font-weight:600; }
.bw-verdict.bw-confirmed { color:#2c7a30; }
.bw-verdict.bw-partial   { color:#a06800; }
.bw-verdict.bw-diverged  { color:#a32a2a; }
.bw-num { font-variant-numeric: tabular-nums; }
.bw-unmatched { color:#666; font-size:0.9em; margin-top:14px; }

/* ---- comprehensive clinical layer ---- */
.cl-hero { display:flex; gap:20px; align-items:center; flex-wrap:wrap;
           background:linear-gradient(135deg,#f8fafc,#eef2f7); border:1px solid #dde3ea;
           border-radius:14px; padding:18px 22px; margin:16px 0; }
.cl-score-ring { width:96px; height:96px; border-radius:50%; display:flex;
           align-items:center; justify-content:center; font-size:1.9em; font-weight:800;
           color:#fff; flex:0 0 auto; }
.cl-hero .cl-summary { flex:1; min-width:220px; }
.cl-hero h2 { border:none; margin:0 0 4px; font-size:1.35em; }
.cl-chips { display:flex; gap:8px; flex-wrap:wrap; margin-top:8px; }
.cl-chip { font-size:.8em; padding:3px 10px; border-radius:20px; background:#fff; border:1px solid #dde3ea; }
.cl-sys { border:1px solid #e3e7ec; border-radius:12px; margin:14px 0; overflow:hidden; }
.cl-sys-head { display:flex; align-items:center; gap:12px; padding:11px 16px;
           background:#f7f9fb; border-bottom:1px solid #e3e7ec; cursor:pointer; }
.cl-sys-head h3 { margin:0; font-size:1.05em; flex:1; }
.cl-sys-bar { width:120px; height:8px; border-radius:5px; background:#e6e9ee; overflow:hidden; }
.cl-sys-bar > div { height:100%; }
.cl-sys-score { font-weight:700; font-variant-numeric:tabular-nums; min-width:34px; text-align:right; }
table.cl { width:100%; border-collapse:collapse; font-size:0.9em; }
table.cl th, table.cl td { padding:8px 12px; border-bottom:1px solid #eef1f4; vertical-align:top; text-align:left; }
table.cl th { background:#fbfcfd; font-size:.82em; text-transform:uppercase; letter-spacing:.04em; color:#667; }
.cl-name { font-weight:600; } .cl-name .cl-derived { font-size:.7em; color:#8a94a3; font-weight:500;
           background:#eef2f7; border-radius:4px; padding:1px 5px; margin-left:6px; vertical-align:middle; }
.cl-badge { display:inline-block; font-size:.78em; font-weight:700; padding:2px 9px; border-radius:20px; white-space:nowrap; }
.cl-optimal   { background:#e7f6ea; color:#1a7f37; }
.cl-normal    { background:#eef2f7; color:#41505f; }
.cl-borderline{ background:#fff4d6; color:#8a6100; }
.cl-high, .cl-low { background:#fde4e4; color:#b3261e; }
.cl-critical  { background:#b3261e; color:#fff; }
.cl-ref { color:#78828e; font-size:.85em; }
.cl-note { margin-top:5px; font-size:.86em; color:#2b5f8e; background:#eef5fb;
           border-left:3px solid #7ab3e0; padding:5px 9px; border-radius:0 6px 6px 0; }
.cl-action { margin-top:4px; font-size:.86em; color:#444; }
.cl-track { position:relative; height:8px; border-radius:5px; background:#e9edf1; margin:6px 0 2px; }
.cl-opt-band { position:absolute; height:100%; background:#bfe6c8; border-radius:5px; }
.cl-dot { position:absolute; top:-3px; width:14px; height:14px; border-radius:50%;
          border:2px solid #fff; box-shadow:0 0 0 1px rgba(0,0,0,.15); transform:translateX(-50%); }
.cl-flags { background:#fff8f8; border:1px solid #f3d6d6; border-radius:12px; padding:6px 16px 12px; margin:14px 0; }
.cl-flags h3 { color:#b3261e; }
.cl-disclaimer { color:#7a828c; font-size:.82em; margin-top:8px; }
</style>
"""


def _score_color(score: Optional[int]) -> str:
    if score is None:
        return "#8a94a3"
    if score >= 85:
        return "#1a7f37"
    if score >= 65:
        return "#8a6100"
    return "#b3261e"


def _range_text(r: Dict) -> str:
    def fmt(lo, hi):
        if lo is not None and hi is not None:
            return f"{lo:g}–{hi:g}"
        if hi is not None:
            return f"&lt;{hi:g}"
        if lo is not None:
            return f"&gt;{lo:g}"
        return "—"
    clin = fmt(r["clinical_low"], r["clinical_high"])
    opt = fmt(r["optimal_low"], r["optimal_high"])
    return (f"<span class='cl-ref'>ref {clin} {_esc(r['unit'])}"
            f" · optimal {opt}</span>")


def _range_bar(r: Dict) -> str:
    """A small track showing the optimal band + a dot at the measured value."""
    v = r["value"]
    lo = r["optimal_low"] if r["optimal_low"] is not None else r["clinical_low"]
    hi = r["optimal_high"] if r["optimal_high"] is not None else r["clinical_high"]
    # establish a display scale
    lo_ref = r["clinical_low"] if r["clinical_low"] is not None else (lo if lo is not None else 0)
    hi_ref = r["clinical_high"] if r["clinical_high"] is not None else (hi if hi is not None else v * 1.6 or 1)
    span_pts = [p for p in (lo_ref, hi_ref, lo, hi, v) if p is not None]
    tmin, tmax = min(span_pts), max(span_pts)
    pad = (tmax - tmin) * 0.18 or (abs(tmax) * 0.2 or 1)
    tmin -= pad
    tmax += pad
    if tmax <= tmin:
        return ""

    def pct(x):
        return max(0.0, min(100.0, 100 * (x - tmin) / (tmax - tmin)))

    band = ""
    if r["optimal_low"] is not None or r["optimal_high"] is not None:
        bl = pct(r["optimal_low"]) if r["optimal_low"] is not None else 0.0
        bh = pct(r["optimal_high"]) if r["optimal_high"] is not None else 100.0
        band = f"<div class='cl-opt-band' style='left:{bl:.1f}%;width:{max(0,bh-bl):.1f}%'></div>"
    dot_color = _STATUS_META[r["status"]][1]
    dot_bg = {"cl-optimal": "#1a7f37", "cl-normal": "#41505f", "cl-borderline": "#d29922",
              "cl-high": "#b3261e", "cl-low": "#b3261e", "cl-critical": "#7a0f0f"}.get(dot_color, "#41505f")
    return (f"<div class='cl-track'>{band}"
            f"<div class='cl-dot' style='left:{pct(v):.1f}%;background:{dot_bg}'></div></div>")


def _render_clinical_section(clinical: Dict) -> str:
    if not clinical or not clinical.get("available"):
        return ""
    overall = clinical.get("overall_score")
    ring = (f"<div class='cl-score-ring' style='background:{_score_color(overall)}'>"
            f"{overall if overall is not None else '—'}</div>")
    hero = f"""
    <div class="cl-hero">
      {ring}
      <div class="cl-summary">
        <h2>Comprehensive Blood Panel</h2>
        <div style="color:#556">Each biomarker scored against standard clinical ranges and tighter
        functional/optimal targets, grouped by body system, with calculated markers and your own
        genotype folded in.</div>
        <div class="cl-chips">
          <span class="cl-chip"><strong>{clinical['n_markers']}</strong> markers</span>
          <span class="cl-chip"><strong>{clinical['n_derived']}</strong> calculated</span>
          <span class="cl-chip" style="color:#1a7f37"><strong>{clinical['n_optimal']}</strong> optimal</span>
          <span class="cl-chip" style="color:#b3261e"><strong>{clinical['n_flagged']}</strong> flagged</span>
        </div>
      </div>
    </div>
    """

    # priority flags
    flags_html = ""
    if clinical["flags"]:
        items = ""
        for f in clinical["flags"][:12]:
            note = f"<div class='cl-note'>🧬 {_esc(f['genotype_note'])}</div>" if f.get("genotype_note") else ""
            action = f"<div class='cl-action'>{_esc(f['action'])}</div>" if f.get("action") else ""
            items += (f"<div style='margin:9px 0'><span class='cl-badge {f['status_class']}'>"
                      f"{_esc(f['status_label'])}</span> <strong>{_esc(f['name'])}</strong> "
                      f"<span class='bw-num'>{_esc(f['value'])} {_esc(f['unit'])}</span> "
                      f"{_range_text(f)}{action}{note}</div>")
        flags_html = f"<div class='cl-flags'><h3>⚑ Priority — out-of-range markers</h3>{items}</div>"

    # per-system panels
    sys_html = ""
    for s in clinical["systems"]:
        rows = ""
        for m in s["markers"]:
            note = f"<div class='cl-note'>🧬 {_esc(m['genotype_note'])}</div>" if m.get("genotype_note") else ""
            action = f"<div class='cl-action'>{_esc(m['action'])}</div>" if m.get("action") else ""
            derived_tag = "<span class='cl-derived'>calculated</span>" if m["derived"] else ""
            rows += f"""
            <tr>
              <td style="width:34%"><span class="cl-name">{_esc(m['name'])}{derived_tag}</span>
                  <div style="color:#8a94a3;font-size:.82em">{_esc(m['desc'])}</div></td>
              <td style="width:16%" class="bw-num"><strong>{_esc(m['value'])}</strong> {_esc(m['unit'])}<br>
                  {_range_text(m)}</td>
              <td style="width:18%"><span class="cl-badge {m['status_class']}">{_esc(m['status_label'])}</span>
                  {_range_bar(m)}</td>
              <td>{action}{note}{'' if (action or note) else '<span style=color:#aeb6c0>—</span>'}</td>
            </tr>"""
        col = _score_color(s["score"])
        sys_html += f"""
        <details class="cl-sys" open>
          <summary class="cl-sys-head">
            <h3>{_esc(s['system'])}</h3>
            <span style="color:#8a94a3;font-size:.85em">{s['n_flagged']} flagged / {s['n_markers']}</span>
            <div class="cl-sys-bar"><div style="width:{s['score']}%;background:{col}"></div></div>
            <span class="cl-sys-score" style="color:{col}">{s['score']}</span>
          </summary>
          <table class="cl">
            <tr><th>Biomarker</th><th>Value / range</th><th>Status</th><th>What it means / action</th></tr>
            {rows}
          </table>
        </details>"""

    unrec = ""
    if clinical.get("unrecognized"):
        unrec = (f"<p class='bw-unmatched'>Supplied values not recognised as biomarkers: "
                 f"{', '.join(_esc(u) for u in clinical['unrecognized'])}.</p>")

    return (hero + flags_html + sys_html + unrec +
            "<p class='cl-disclaimer'>Reference and optimal ranges are typical adult values and "
            "vary by lab, assay, age and sex. This is educational, not a diagnosis — discuss "
            "abnormal or borderline results with your clinician.</p>")


def render_bloodwork_html(result: Dict, file_label: str = "") -> str:
    """Render the bloodwork analysis as a standalone HTML document: the
    comprehensive clinical panel first, then the genetic-prediction comparison."""
    clinical_html = _render_clinical_section(result.get("clinical"))

    gene_header = (
        "<h1 style='margin-top:34px'>Genetic-Prediction Comparison</h1>"
        "<p style='color:#666'>How your measured labs line up with your "
        "<em>polygenic</em> (PheWAS) predictions — a different lens from the "
        "clinical ranges above.</p>"
    ) if clinical_html else ""

    if result.get("status") == "no_phewas":
        body = gene_header + (
            f"<p>PheWAS predictions weren't available, so the genetic-comparison "
            f"layer is skipped. The comprehensive clinical panel above covers all "
            f"{result['n_labs_supplied']} supplied values.</p>"
        )
    elif not result.get("rows"):
        body = gene_header + (
            f"<p>None of the {result['n_labs_supplied']} supplied lab values matched "
            f"a scored PheWAS biomarker.</p>"
        )
    else:
        rows_html = []
        for r in result["rows"]:
            rows_html.append(f"""
              <tr class="{r['verdict_class']}">
                <td><strong>{_esc(r['trait'])}</strong><br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['category'])}</span></td>
                <td class="bw-num">{_esc(r['predicted'])} {_esc(r['unit'])}<br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['predicted_tier'])}</span>
                    {f'<br><span style="color:#999;font-size:0.78em">genetic coverage {r["callability_pct"]:.0f}% ({r.get("n_used",0)} SNPs)</span>' if r.get('callability_pct') else ''}</td>
                <td class="bw-num">{_esc(r['actual'])} {_esc(r['unit'])}<br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['actual_tier'])}</span></td>
                <td class="bw-num">{r['delta_abs']:+.2f}<br>
                    <span style="color:#777;font-size:0.85em">{r['delta_sd']:+.2f} SD</span></td>
                <td class="bw-verdict {r['verdict_class']}">{_esc(r['verdict'])}</td>
                <td>{_esc(r['interpretation'])}</td>
              </tr>
            """)
        body = gene_header + _GENE_HOWTO_HTML + f"""
          <div class="bw-summary">
            <div class="bw-stat"><div class="v">{result['n_matched']}</div>
                <div class="l">Compared</div></div>
            <div class="bw-stat"><div class="v" style="color:#2c7a30">{result['n_confirmed']}</div>
                <div class="l">Confirmed</div></div>
            <div class="bw-stat"><div class="v" style="color:#a06800">{result['n_partial']}</div>
                <div class="l">Partial</div></div>
            <div class="bw-stat"><div class="v" style="color:#a32a2a">{result['n_diverged']}</div>
                <div class="l">Diverged</div></div>
            <div class="bw-stat"><div class="v">{result['accuracy_pct']}%</div>
                <div class="l">Accuracy</div></div>
          </div>
          <table class="bw">
            <tr>
              <th>Biomarker</th><th>Predicted (genetics)</th><th>Actual (lab)</th>
              <th>Δ</th><th>Verdict</th><th>Interpretation</th>
            </tr>
            {''.join(rows_html)}
          </table>
        """
        if result.get("unmatched"):
            body += (
                f"<p class='bw-unmatched'>Supplied labs not matched to a scored "
                f"genetic prediction: {', '.join(_esc(u) for u in result['unmatched'])}.</p>"
            )
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Comprehensive Blood Work Analysis{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_BW_CSS}</head><body><div class="bw-wrap">
<h1>Comprehensive Blood Work Analysis</h1>
<p style="color:#666">{_esc(result.get('notes',''))}</p>
{clinical_html}
{body}
</div></body></html>"""


# Genetic-comparison explainer (attached to the genetic section, not the top).
_GENE_HOWTO_HTML = """
<details style="background:#f6f6f7;border:1px solid #ddd;border-radius:8px;
                padding:10px 14px;margin:12px 0;">
  <summary style="cursor:pointer;font-weight:600">
    How to read this table — what "Diverged" really means
  </summary>
  <p style="margin-top:10px;color:#444">
    Polygenic predictions built from common GWAS variants typically explain
    only <strong>5–20% of the variance</strong> in any single biomarker. The
    remaining 80–95% is driven by diet, exercise, sleep, medications, illness,
    age, sex, season, and chance. The thresholds used here are calibrated to
    that reality:
  </p>
  <ul style="color:#444">
    <li><strong>Confirmed</strong> (|Δ| &le; 0.75 SD) — measured value sits well
        within where genetics alone would predict.</li>
    <li><strong>Partial</strong> (|Δ| &le; 1.5 SD) — small lifestyle / environmental
        contribution on top of genetics.</li>
    <li><strong>Diverged</strong> (|Δ| &gt; 1.5 SD) — large non-genetic driver
        dominating; investigate diet, training, medication, stress, illness.
        This is the most actionable row in the table.</li>
  </ul>
  <p style="color:#444">
    "Diverged" never means the genetic prediction is <em>wrong</em> — it means
    a real-world driver is overriding it. That driver is what's modifiable.
  </p>
</details>
"""
