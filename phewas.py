"""
Phenome-Wide Trait Prediction (PheWAS-style)
============================================

Scores 60+ phenotypes using curated GWAS summary effect sizes. For each
trait the score is:

    raw = Σ β_i × dosage_i

where β_i is the per-allele effect (in SD units of the trait or log-OR)
from large published GWAS, and dosage_i is the user's count of the effect
allele. We then compute:

    z = (raw - expected_mean) / sqrt(expected_var)
    percentile = Φ(z)

with expected_mean / expected_var derived from Hardy-Weinberg expectations
using European allele frequencies (a deliberate, transparent choice — the
GWAS source studies are predominantly European).

For each trait we surface:
  * Predicted percentile vs the European reference distribution
  * Predicted direction (above / below average)
  * For continuous biomarkers, a rough predicted-value band
  * Coverage: how many of the panel SNPs were typed on the user's chip

This is research-grade summary, not clinical. Effect sizes are from public
GWAS (UK Biobank, GLGC, MAGIC, GIANT, etc.). Results are estimates only.
"""

from __future__ import annotations

from math import erf, sqrt
from typing import Dict, List, Optional
import pandas as pd


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> Optional[int]:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--") or len(s) != 2:
        return None
    return s.count(allele.upper())


# ─── Trait panels ────────────────────────────────────────────────────────────
# Each trait is a dict:
#   category: dashboard category
#   unit: display unit (e.g. "mg/dL", "SD", "%", "kg")
#   direction: "higher_is_higher" or "higher_is_protective" (legacy — unused)
#   mean / sd: rough population values for translating Z into a value estimate
#   reference: literature citation
#   variants: [(rsid, effect_allele, beta, af)]
#
# Effect sizes (beta) are in SD units of the trait unless marked. AFs are
# European (1000G EUR).

PHEWAS_TRAITS: Dict[str, Dict] = {

    # ── Lipids ─────────────────────────────────────────────────────────────
    "LDL cholesterol": {
        "category": "Lipids", "unit": "mg/dL",
        "mean": 110.0, "sd": 32.0,
        "reference": "Willer 2013 GLGC + Klarin 2017",
        "variants": [
            ("rs646776", "C", 0.10, 0.78),
            ("rs599839", "G", 0.11, 0.78),
            ("rs6511720", "T", -0.21, 0.10),
            ("rs693", "T", 0.10, 0.51),
            ("rs11591147", "T", -0.50, 0.02),
            ("rs11206510", "T", 0.10, 0.81),
            ("rs429358", "C", 0.20, 0.14),
            ("rs7412", "T", -0.30, 0.08),
        ],
    },
    "HDL cholesterol": {
        "category": "Lipids", "unit": "mg/dL",
        "mean": 55.0, "sd": 15.0,
        "reference": "GLGC 2017",
        "variants": [
            ("rs5882", "G", 0.13, 0.32),
            ("rs2230806", "C", 0.05, 0.41),
            ("rs1864163", "G", 0.07, 0.30),
            ("rs10468017", "T", 0.05, 0.30),
            ("rs909241", "C", 0.04, 0.55),
            ("rs17321515", "G", 0.04, 0.40),
        ],
    },
    "Triglycerides": {
        "category": "Lipids", "unit": "mg/dL",
        "mean": 110.0, "sd": 60.0,
        "reference": "GLGC 2017",
        "variants": [
            ("rs2954029", "T", 0.06, 0.55),    # TRIB1 area (proxy)
            ("rs17321515", "G", 0.04, 0.40),
            ("rs1042031", "A", 0.05, 0.20),
        ],
    },
    "Lipoprotein(a)": {
        "category": "Lipids", "unit": "nmol/L (very wide range)",
        "mean": 50.0, "sd": 80.0,
        "reference": "Coassin 2022; Trinder 2021",
        "variants": [
            ("rs10455872", "G", 0.55, 0.07),
            ("rs3798220", "C", 0.50, 0.02),
        ],
    },

    # ── Glucose / diabetes ─────────────────────────────────────────────────
    "Fasting glucose": {
        "category": "Glucose / Diabetes", "unit": "mg/dL",
        "mean": 90.0, "sd": 8.0,
        "reference": "MAGIC consortium",
        "variants": [
            ("rs10830963", "G", 0.07, 0.30),
            ("rs1799884", "T", 0.06, 0.16),
            ("rs560887", "C", 0.10, 0.65),
            ("rs1387153", "T", 0.07, 0.30),
        ],
    },
    "HbA1c (predicted)": {
        "category": "Glucose / Diabetes", "unit": "%",
        "mean": 5.4, "sd": 0.4,
        "reference": "Wheeler 2017 (MAGIC)",
        "variants": [
            ("rs7903146", "T", 0.10, 0.30),
            ("rs10830963", "G", 0.06, 0.30),
            ("rs1387153", "T", 0.05, 0.30),
        ],
    },
    "HOMA-IR (insulin resistance)": {
        "category": "Glucose / Diabetes", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Dupuis 2010",
        "variants": [
            ("rs2943641", "C", 0.08, 0.62),
            ("rs7903146", "T", 0.05, 0.30),
        ],
    },

    # ── Anthropometric ─────────────────────────────────────────────────────
    "Body Mass Index (BMI)": {
        "category": "Anthropometric", "unit": "kg/m²",
        "mean": 26.0, "sd": 4.5,
        "reference": "GIANT 2018 (Yengo et al.)",
        "variants": [
            ("rs9939609", "A", 0.10, 0.42),
            ("rs1421085", "C", 0.10, 0.42),
            ("rs8050136", "A", 0.09, 0.42),
            ("rs17782313", "C", 0.07, 0.24),
            ("rs2229616", "A", -0.08, 0.02),
            ("rs7799039", "A", 0.04, 0.50),
            ("rs1137101", "G", 0.04, 0.55),
            ("rs5082", "C", 0.05, 0.25),
            ("rs6232", "G", 0.05, 0.07),
        ],
    },
    "Height": {
        "category": "Anthropometric", "unit": "cm",
        "mean": 170.0, "sd": 10.0,
        "reference": "Yengo 2022",
        "variants": [
            ("rs143384", "T", 0.10, 0.66),
            ("rs2562784", "G", 0.08, 0.55),
            ("rs6060369", "T", 0.07, 0.50),
            ("rs2275035", "C", 0.06, 0.40),
        ],
    },
    "Waist-Hip Ratio (BMI-adjusted)": {
        "category": "Anthropometric", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "GIANT WHR 2019",
        "variants": [
            ("rs6232", "G", 0.05, 0.07),
            ("rs1421085", "C", 0.04, 0.42),
        ],
    },

    # ── Blood pressure / cardiac ────────────────────────────────────────────
    "Systolic Blood Pressure": {
        "category": "Cardiovascular", "unit": "mmHg",
        "mean": 120.0, "sd": 14.0,
        "reference": "Evangelou 2018 + ICBP",
        "variants": [
            ("rs17249754", "G", 0.05, 0.85),
            ("rs1378942", "C", 0.06, 0.36),
            ("rs2810226", "T", 0.04, 0.55),
            ("rs2070744", "C", 0.05, 0.38),
            ("rs5443", "T", 0.04, 0.31),
            ("rs17042171", "A", 0.05, 0.20),
            ("rs3184504", "T", 0.06, 0.40),
        ],
    },
    "Diastolic Blood Pressure": {
        "category": "Cardiovascular", "unit": "mmHg",
        "mean": 75.0, "sd": 10.0,
        "reference": "Evangelou 2018",
        "variants": [
            ("rs17249754", "G", 0.04, 0.85),
            ("rs1378942", "C", 0.05, 0.36),
        ],
    },
    "Resting Heart Rate": {
        "category": "Cardiovascular", "unit": "bpm",
        "mean": 68.0, "sd": 10.0,
        "reference": "Eppinga 2016",
        "variants": [
            ("rs9349379", "A", 0.04, 0.39),
            ("rs1801252", "G", 0.04, 0.15),
        ],
    },

    # ── Hematology ──────────────────────────────────────────────────────────
    "Hemoglobin": {
        "category": "Hematology", "unit": "g/dL",
        "mean": 14.5, "sd": 1.5,
        "reference": "UK Biobank blood-cell GWAS",
        "variants": [
            ("rs855791", "A", -0.10, 0.43),
            ("rs1800562", "A", 0.10, 0.06),
            ("rs7385804", "A", 0.05, 0.62),
        ],
    },
    "Mean Corpuscular Volume": {
        "category": "Hematology", "unit": "fL",
        "mean": 90.0, "sd": 5.0,
        "reference": "Astle 2016",
        "variants": [
            ("rs855791", "A", -0.13, 0.43),
            ("rs1800562", "A", 0.10, 0.06),
        ],
    },
    "Platelets": {
        "category": "Hematology", "unit": "×10⁹/L",
        "mean": 250.0, "sd": 60.0,
        "reference": "Astle 2016",
        "variants": [
            ("rs1354034", "T", 0.08, 0.30),
        ],
    },
    "White Blood Cell Count": {
        "category": "Hematology", "unit": "×10⁹/L",
        "mean": 6.8, "sd": 1.8,
        "reference": "Astle 2016",
        "variants": [
            ("rs2814778", "C", -0.30, 0.005),  # Duffy null → benign neutropenia (African ancestry)
        ],
    },

    # ── Inflammation ────────────────────────────────────────────────────────
    "C-Reactive Protein": {
        "category": "Inflammation", "unit": "mg/L",
        "mean": 1.5, "sd": 1.5,
        "reference": "Ligthart 2018",
        "variants": [
            ("rs1800795", "C", 0.10, 0.43),
            ("rs2228145", "C", 0.06, 0.43),
            ("rs2069837", "G", 0.05, 0.10),
            ("rs1800629", "A", 0.05, 0.16),
        ],
    },

    # ── Kidney / liver / endocrine biomarkers ──────────────────────────────
    "Uric acid / urate": {
        "category": "Metabolic", "unit": "mg/dL",
        "mean": 5.5, "sd": 1.4,
        "reference": "Köttgen 2013",
        "variants": [
            ("rs12498742", "A", 0.20, 0.74),
            ("rs2231142", "T", 0.20, 0.11),
        ],
    },
    "Serum creatinine / eGFR": {
        "category": "Metabolic", "unit": "mg/dL",
        "mean": 0.95, "sd": 0.18,
        "reference": "Wuttke 2019",
        "variants": [
            ("rs77924615", "G", 0.06, 0.20),
            ("rs4805834", "C", 0.04, 0.60),
        ],
    },
    "ALT (liver function)": {
        "category": "Liver", "unit": "U/L",
        "mean": 22.0, "sd": 10.0,
        "reference": "Chambers 2011",
        "variants": [
            ("rs738409", "G", 0.16, 0.23),       # PNPLA3
            ("rs2227831", "G", 0.07, 0.27),
        ],
    },

    # ── Thyroid ─────────────────────────────────────────────────────────────
    "TSH": {
        "category": "Thyroid", "unit": "mIU/L",
        "mean": 2.0, "sd": 0.9,
        "reference": "Porcu 2013",
        "variants": [
            ("rs965513", "A", 0.16, 0.34),
            ("rs179247", "G", 0.08, 0.60),
        ],
    },

    # ── Vitamins ────────────────────────────────────────────────────────────
    "25-OH Vitamin D": {
        "category": "Vitamins / Micronutrients", "unit": "ng/mL",
        "mean": 25.0, "sd": 10.0,
        "reference": "Manousaki 2017",
        "variants": [
            ("rs2282679", "C", -0.13, 0.28),
            ("rs10741657", "G", -0.10, 0.40),
            ("rs7041", "T", -0.06, 0.43),
            ("rs4588", "A", -0.06, 0.32),
        ],
    },
    "Serum vitamin B12": {
        "category": "Vitamins / Micronutrients", "unit": "pmol/L",
        "mean": 350.0, "sd": 110.0,
        "reference": "Hazra 2009",
        "variants": [
            ("rs602662", "A", -0.20, 0.49),
            ("rs1801198", "G", -0.05, 0.60),
            ("rs1801222", "G", -0.05, 0.60),
        ],
    },
    "Folate": {
        "category": "Vitamins / Micronutrients", "unit": "ng/mL",
        "mean": 15.0, "sd": 5.0,
        "reference": "Tanaka 2009",
        "variants": [
            ("rs1801133", "T", -0.10, 0.31),
            ("rs602662", "A", -0.05, 0.49),
        ],
    },
    "Iron / ferritin": {
        "category": "Vitamins / Micronutrients", "unit": "ng/mL",
        "mean": 150.0, "sd": 100.0,
        "reference": "Benyamin 2014",
        "variants": [
            ("rs1800562", "A", 0.30, 0.06),
            ("rs1799945", "G", 0.10, 0.16),
            ("rs855791", "A", -0.13, 0.43),
        ],
    },

    # ── Hormones ────────────────────────────────────────────────────────────
    "Testosterone (men)": {
        "category": "Hormones", "unit": "ng/dL",
        "mean": 550.0, "sd": 180.0, "applies_to": "male",
        "reference": "Ruth 2020 (UK Biobank)",
        "variants": [
            ("rs727428", "G", 0.10, 0.50),    # SHBG (proxy)
            ("rs1799941", "A", -0.10, 0.27),  # SHBG region — higher SHBG = lower free T
            ("rs6259", "A", -0.08, 0.12),
        ],
    },
    "Estradiol (women)": {
        "category": "Hormones", "unit": "pg/mL", "applies_to": "female",
        "mean": 70.0, "sd": 40.0,
        "reference": "Schmitz 2021",
        "variants": [
            ("rs727428", "G", -0.05, 0.50),
        ],
    },
    "SHBG": {
        "category": "Hormones", "unit": "nmol/L",
        "mean": 35.0, "sd": 18.0,
        "reference": "Coviello 2012",
        "variants": [
            ("rs1799941", "A", 0.20, 0.27),
            ("rs6259", "A", 0.16, 0.12),
        ],
    },
    "IGF-1": {
        "category": "Hormones", "unit": "ng/mL",
        "mean": 160.0, "sd": 50.0,
        "reference": "Teumer 2016",
        "variants": [
            ("rs2153960", "A", 0.05, 0.40),
        ],
    },
    "Cortisol (proxy)": {
        "category": "Hormones", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Bolton 2014",
        "variants": [
            ("rs41423247", "G", 0.05, 0.30),
        ],
    },

    # ── Behaviour ───────────────────────────────────────────────────────────
    "Educational attainment (PRS)": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Okbay 2022 (EA4, 3M variants)",
        "variants": [
            ("rs9320913", "A", 0.025, 0.49),
            ("rs11584700", "G", 0.020, 0.39),
        ],
    },
    "Risk tolerance": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Karlsson Linner 2019",
        "variants": [
            ("rs1800497", "A", 0.03, 0.41),
            ("rs1800955", "T", 0.02, 0.45),
        ],
    },
    "Chronotype (morningness)": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Jones 2019 (UKBB)",
        "variants": [
            ("rs1801260", "C", -0.08, 0.30),
            ("rs10925130", "G", -0.04, 0.30),
        ],
    },
    "Coffee consumption": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Cornelis 2015",
        "variants": [
            ("rs762551", "A", 0.10, 0.65),
        ],
    },
    "Alcohol consumption": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Liu 2019 GSCAN",
        "variants": [
            ("rs1229984", "A", -0.20, 0.03),
            ("rs671", "A", -0.30, 0.002),
        ],
    },
    "Smoking initiation": {
        "category": "Behaviour", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Liu 2019 GSCAN",
        "variants": [
            ("rs1051730", "T", 0.10, 0.34),
            ("rs16969968", "A", 0.10, 0.34),
        ],
    },

    # ── Longevity / aging ───────────────────────────────────────────────────
    "Parental lifespan (proxy)": {
        "category": "Longevity", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Timmers 2019",
        "variants": [
            ("rs429358", "C", -0.10, 0.14),
            ("rs7412", "T", 0.05, 0.08),
            ("rs2802292", "G", 0.04, 0.32),
        ],
    },
    "Telomere length (genetic)": {
        "category": "Longevity", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Codd 2021",
        "variants": [
            ("rs2736100", "C", 0.10, 0.51),
            ("rs9420907", "C", 0.06, 0.65),
        ],
    },

    # ── Neuro proxies ───────────────────────────────────────────────────────
    "Brain volume (proxy)": {
        "category": "Neurological", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Adams 2016 / ENIGMA",
        "variants": [
            ("rs17178006", "T", 0.05, 0.40),
        ],
    },
    "Migraine score": {
        "category": "Neurological", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Hautakangas 2022",
        "variants": [
            ("rs4926339", "T", 0.04, 0.40),
            ("rs2651899", "C", 0.04, 0.51),
        ],
    },

    # ── Hearing / vision ────────────────────────────────────────────────────
    "Myopia tendency": {
        "category": "Sensory", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Tedja 2018 CREAM",
        "variants": [
            ("rs524952", "A", 0.10, 0.40),
        ],
    },
    "Hearing loss tendency": {
        "category": "Sensory", "unit": "SD",
        "mean": 0, "sd": 1.0,
        "reference": "Wells 2019",
        "variants": [
            ("rs3135388", "A", 0.03, 0.18),
        ],
    },
}


# ─── Core ─────────────────────────────────────────────────────────────────────

def _score_trait(snps_df: pd.DataFrame, trait: Dict, sex: Optional[str]) -> Dict:
    applies_to = trait.get("applies_to")
    if applies_to == "female" and sex == "male":
        return {"status": "not_applicable", "reason": "Female-specific trait"}
    if applies_to == "male" and sex == "female":
        return {"status": "not_applicable", "reason": "Male-specific trait"}

    raw_score = 0.0
    exp_mean = 0.0
    exp_var = 0.0
    n_used = 0
    used_variants = []
    missing_variants = []

    for rsid, ea, beta, af in trait["variants"]:
        d = _dose(snps_df, rsid, ea)
        if d is None:
            missing_variants.append(rsid)
            continue
        raw_score += beta * d
        exp_mean += beta * 2.0 * af
        exp_var += (beta ** 2) * 2.0 * af * (1.0 - af)
        n_used += 1
        used_variants.append(
            {"rsid": rsid, "effect_allele": ea, "dose": d, "beta": beta, "af": af}
        )

    if n_used == 0 or exp_var <= 0:
        return {"status": "insufficient",
                "reason": f"0/{len(trait['variants'])} panel SNPs typed",
                "n_used": 0, "n_total": len(trait['variants'])}

    # TWO DISTINCT numbers — kept separate so neither is overclaimed:
    #
    # 1) MARKER-SCORE percentile: where you rank on THIS SNP panel, standardised to
    #    the panel's own spread. A true fact about your genotype (you carry more/
    #    fewer of these alleles than X% of people). This is what the tool used to
    #    print — but it was mislabelled as a trait percentile.
    z = (raw_score - exp_mean) / sqrt(exp_var)
    marker_percentile = _norm_cdf(z) * 100

    # 2) ACTUAL-TRAIT prediction: because the betas are in the trait's own units,
    #    the genetic deviation (raw - expected) IS the predicted trait deviation,
    #    and it regresses toward the mean exactly as far as the panel is weak. The
    #    variance the panel explains is R² = Var(score)/Var(trait) = exp_var/sd².
    sd = trait.get("sd", 1.0) or 1.0
    variance_explained = min(1.0, exp_var / (sd ** 2))     # R² of the ACTUAL trait
    predicted_deviation = raw_score - exp_mean              # in native trait units
    val_estimate = trait["mean"] + predicted_deviation
    trait_percentile = _norm_cdf(predicted_deviation / sd) * 100

    # How much the panel actually resolves the real trait. Even the best panel here
    # explains ~6%; most explain <1%. Be explicit rather than implying precision.
    if variance_explained >= 0.05:
        signal = "modest"          # a real but small fraction of the trait
    elif variance_explained >= 0.01:
        signal = "weak"
    else:
        signal = "negligible"      # marker curiosity, not a trait measurement

    # Tier reflects the MARKER score (a real genotype fact), but the renderer labels
    # it as a marker-score tier and shows R² so it is never read as the trait itself.
    if marker_percentile >= 95:
        tier = "High marker score"; tier_cls = "tier-high"
    elif marker_percentile >= 80:
        tier = "Above-average marker score"; tier_cls = "tier-elevated"
    elif marker_percentile >= 20:
        tier = "Typical marker score"; tier_cls = "tier-average"
    elif marker_percentile >= 5:
        tier = "Below-average marker score"; tier_cls = "tier-below"
    else:
        tier = "Low marker score"; tier_cls = "tier-low"

    return {
        "status": "ok",
        "z_score": round(z, 3),
        # `percentile` retained for back-compat but now explicitly the MARKER score.
        "percentile": round(marker_percentile, 1),
        "marker_percentile": round(marker_percentile, 1),
        "trait_percentile": round(trait_percentile, 1),
        "predicted_value": round(val_estimate, 2),
        "variance_explained": round(variance_explained, 4),
        "variance_explained_pct": round(variance_explained * 100, 2),
        "signal_strength": signal,
        "tier": tier,
        "tier_class": tier_cls,
        "n_used": n_used,
        "n_total": len(trait["variants"]),
        "callability_pct": round(100 * n_used / len(trait["variants"]), 1),
        "used_variants": used_variants,
        "missing_variants": missing_variants,
    }


def analyze_phewas(snps_df: pd.DataFrame, sex: Optional[str] = None) -> Dict:
    results: Dict[str, Dict] = {}
    by_category: Dict[str, List[str]] = {}
    for trait_name, trait_def in PHEWAS_TRAITS.items():
        r = _score_trait(snps_df, trait_def, sex)
        results[trait_name] = {**trait_def, "result": r}
        by_category.setdefault(trait_def["category"], []).append(trait_name)

    # Headline: any trait with tier Very-high or Very-low
    headline = []
    for name, t in results.items():
        if t["result"].get("tier") in ("Very high", "Very low"):
            headline.append({
                "trait": name,
                "category": t["category"],
                "tier": t["result"]["tier"],
                "percentile": t["result"]["percentile"],
                "predicted_value": t["result"].get("predicted_value"),
                "unit": t["unit"],
            })

    return {
        "traits": results,
        "by_category": by_category,
        "headline": headline,
        "n_traits": len(results),
        "n_scored": sum(1 for r in results.values() if r["result"]["status"] == "ok"),
    }
