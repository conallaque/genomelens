"""
Polygenic / Curated-Variant Risk Scores
---------------------------------------

Per-trait weighted risk scores computed from genotyped variants with
published effect sizes (log-OR from GWAS / replicated meta-analyses).

This is NOT a million-SNP clinical-grade PGS. It is a transparent,
literature-derived curated-variant risk score (CVRS): a weighted sum
of effect alleles at well-replicated GWAS hits, expressed as a Z-score
against a European reference distribution and translated into a
percentile and risk tier.

Math:
    raw_score = Σ log(OR_i) × dosage_i
    expected_mean = Σ log(OR_i) × 2 × AF_i        (Hardy-Weinberg)
    expected_var  = Σ log(OR_i)² × 2 × AF_i × (1-AF_i)
    z_score = (raw_score - expected_mean) / sqrt(expected_var)
    percentile ≈ Φ(z_score)   (standard normal CDF)

Missing variants are dropped from BOTH the score and the expectation,
so callability does not bias the result.

Limitations (called out in the report):
  * Curated subset, not full PGS — effect magnitude understates true PRS.
  * Reference distribution is European; non-European ancestry interpretation
    is less reliable.
  * Lifetime-risk estimates are population-derived approximations.
  * Not a clinical diagnostic.
"""

from math import sqrt, erf, log
from typing import Dict, List, Optional
import pandas as pd


# ─── Helper: standard-normal CDF (so we don't need scipy) ─────────────────────
def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


# ─── PRS panels ───────────────────────────────────────────────────────────────
#
# Each variant entry:
#   rsid               — rsID
#   effect_allele      — the risk-raising allele (or protective if log_or < 0)
#   log_or             — natural log of per-allele OR from large GWAS
#   af                 — European reference allele frequency for effect_allele
#   gene               — gene symbol for display
#
# Effect sizes derived from published GWAS meta-analyses / consortium reports
# (CARDIoGRAMplusC4D, DIAGRAM/MAGIC, BCAC, PRACTICAL, IGAP, AFGen, GIANT). Where
# the same variant appears in multiple panels (e.g. APOE in AD), entries are
# panel-specific.
#

PRS_PANELS: Dict[str, Dict] = {
    "Coronary Artery Disease": {
        "trait_short": "CAD",
        "description": (
            "Lifetime risk of clinical coronary artery disease (myocardial "
            "infarction, angina, coronary revascularisation). Combines "
            "well-replicated CAD-locus variants with classic Lp(a) and "
            "LDL-pathway effects."
        ),
        "reference": "CARDIoGRAMplusC4D Consortium 2015 + Inouye 2018 + Klarin 2017",
        "population_lifetime_risk": "~30–40% lifetime (men 50%+, women ~30%)",
        "high_tier_implication": (
            "Top 5% PRS confers ~3-fold relative CAD risk. Aggressive "
            "primary prevention is well-supported by clinical trials."
        ),
        "variants": [
            {"rsid": "rs10455872", "effect_allele": "G", "log_or": 0.378, "af": 0.07, "gene": "LPA"},      # OR 1.46
            {"rsid": "rs3798220",  "effect_allele": "C", "log_or": 0.531, "af": 0.02, "gene": "LPA"},      # OR 1.7
            {"rsid": "rs10757278", "effect_allele": "G", "log_or": 0.262, "af": 0.50, "gene": "9p21"},     # OR 1.3
            {"rsid": "rs1333049",  "effect_allele": "C", "log_or": 0.247, "af": 0.49, "gene": "9p21 CDKN2B-AS1"},
            {"rsid": "rs2383206",  "effect_allele": "G", "log_or": 0.231, "af": 0.49, "gene": "9p21"},
            {"rsid": "rs646776",   "effect_allele": "C", "log_or": 0.140, "af": 0.78, "gene": "SORT1"},
            {"rsid": "rs599839",   "effect_allele": "G", "log_or": 0.124, "af": 0.78, "gene": "PSRC1"},
            {"rsid": "rs12740374", "effect_allele": "T", "log_or": 0.148, "af": 0.78, "gene": "CELSR2"},
            {"rsid": "rs11591147", "effect_allele": "T", "log_or": -0.598, "af": 0.02, "gene": "PCSK9 R46L"},  # protective
            {"rsid": "rs11206510", "effect_allele": "T", "log_or": 0.140, "af": 0.81, "gene": "PCSK9"},
            {"rsid": "rs6511720",  "effect_allele": "T", "log_or": -0.163, "af": 0.10, "gene": "LDLR"},        # protective
            {"rsid": "rs693",      "effect_allele": "T", "log_or": 0.075, "af": 0.51, "gene": "APOB"},
            {"rsid": "rs17465637", "effect_allele": "C", "log_or": 0.110, "af": 0.74, "gene": "MIA3"},
            {"rsid": "rs7173743",  "effect_allele": "C", "log_or": 0.085, "af": 0.42, "gene": "MORF4L1"},
            {"rsid": "rs17228212", "effect_allele": "C", "log_or": 0.060, "af": 0.30, "gene": "SMAD3"},
            {"rsid": "rs5882",     "effect_allele": "G", "log_or": -0.080, "af": 0.32, "gene": "CETP I405V"}, # protective
        ],
    },

    "Type 2 Diabetes": {
        "trait_short": "T2D",
        "description": (
            "Lifetime risk of type 2 diabetes mellitus. Captures the "
            "strongest replicated common T2D variants from DIAGRAM/MAGIC "
            "meta-analyses."
        ),
        "reference": "DIAGRAM Consortium 2018 + Mahajan 2018",
        "population_lifetime_risk": "~30–40% lifetime in US/UK populations",
        "high_tier_implication": (
            "Top decile PRS ~2× T2D risk; with obesity, ~4×. Lifestyle "
            "prevention (DPP-style) reduces incidence ~58% at any genetic risk."
        ),
        "variants": [
            {"rsid": "rs7903146",  "effect_allele": "T", "log_or": 0.336, "af": 0.30, "gene": "TCF7L2"},   # OR 1.4
            {"rsid": "rs9939609",  "effect_allele": "A", "log_or": 0.166, "af": 0.42, "gene": "FTO"},      # OR 1.18
            {"rsid": "rs8050136",  "effect_allele": "A", "log_or": 0.149, "af": 0.42, "gene": "FTO"},
            {"rsid": "rs10830963", "effect_allele": "G", "log_or": 0.105, "af": 0.30, "gene": "MTNR1B"},
            {"rsid": "rs1387153",  "effect_allele": "T", "log_or": 0.092, "af": 0.30, "gene": "MTNR1B"},
            {"rsid": "rs13266634", "effect_allele": "C", "log_or": 0.140, "af": 0.72, "gene": "SLC30A8"},
            {"rsid": "rs5219",     "effect_allele": "T", "log_or": 0.122, "af": 0.40, "gene": "KCNJ11"},
            {"rsid": "rs7754840",  "effect_allele": "C", "log_or": 0.140, "af": 0.30, "gene": "CDKAL1"},
            {"rsid": "rs1111875",  "effect_allele": "C", "log_or": 0.140, "af": 0.55, "gene": "HHEX"},
            {"rsid": "rs10811661", "effect_allele": "T", "log_or": 0.182, "af": 0.83, "gene": "CDKN2A/B"},
            {"rsid": "rs2241766",  "effect_allele": "G", "log_or": 0.122, "af": 0.13, "gene": "ADIPOQ"},
            {"rsid": "rs1501299",  "effect_allele": "T", "log_or": 0.075, "af": 0.32, "gene": "ADIPOQ"},
            {"rsid": "rs864745",   "effect_allele": "T", "log_or": 0.092, "af": 0.50, "gene": "JAZF1"},
            {"rsid": "rs2943641",  "effect_allele": "C", "log_or": 0.110, "af": 0.62, "gene": "IRS1"},
            {"rsid": "rs2287019",  "effect_allele": "C", "log_or": 0.080, "af": 0.20, "gene": "GIPR"},
            {"rsid": "rs10010131", "effect_allele": "G", "log_or": 0.075, "af": 0.60, "gene": "WFS1"},
            {"rsid": "rs1801214",  "effect_allele": "C", "log_or": 0.060, "af": 0.40, "gene": "WFS1"},
            {"rsid": "rs10885122", "effect_allele": "G", "log_or": 0.058, "af": 0.86, "gene": "ADRA2A"},
        ],
    },

    "Breast Cancer": {
        "trait_short": "BC",
        "applies_to": "female",
        "description": (
            "Lifetime risk of female breast cancer from common-variant "
            "polygenic burden. Does NOT capture high-penetrance BRCA1/2/PALB2 "
            "mutations — those require separate panel testing."
        ),
        "reference": "BCAC + DRIVE 2017 + Mavaddat 2019 (313-SNP PRS basis)",
        "population_lifetime_risk": "~12% (1 in 8) in Western populations",
        "high_tier_implication": (
            "Top 5% PRS confers ~3-fold relative risk (~30% lifetime). "
            "Eligibility for enhanced screening (annual mammogram ± MRI, "
            "earlier start age) follows risk-stratification guidelines."
        ),
        "variants": [
            {"rsid": "rs2981582",  "effect_allele": "A", "log_or": 0.236, "af": 0.38, "gene": "FGFR2"},  # OR 1.27
            {"rsid": "rs2981579",  "effect_allele": "A", "log_or": 0.236, "af": 0.41, "gene": "FGFR2"},
            {"rsid": "rs1219648",  "effect_allele": "G", "log_or": 0.182, "af": 0.40, "gene": "FGFR2"},
            {"rsid": "rs3803662",  "effect_allele": "T", "log_or": 0.182, "af": 0.26, "gene": "TOX3"},
            {"rsid": "rs17879961", "effect_allele": "C", "log_or": 0.470, "af": 0.005, "gene": "CHEK2 I157T"}, # OR ~1.6
            {"rsid": "rs1799966",  "effect_allele": "G", "log_or": 0.040, "af": 0.40, "gene": "BRCA1 common"},
            {"rsid": "rs1801516",  "effect_allele": "A", "log_or": 0.090, "af": 0.16, "gene": "ATM D1853N"},
            {"rsid": "rs1042522",  "effect_allele": "G", "log_or": 0.050, "af": 0.45, "gene": "TP53 P72R"},
            {"rsid": "rs1801320",  "effect_allele": "C", "log_or": 0.080, "af": 0.05, "gene": "RAD51 G135C"},
            {"rsid": "rs861539",   "effect_allele": "T", "log_or": 0.060, "af": 0.35, "gene": "XRCC3 T241M"},
        ],
    },

    "Prostate Cancer": {
        "trait_short": "PCa",
        "applies_to": "male",
        "description": (
            "Lifetime risk of clinical prostate cancer from common-variant "
            "polygenic burden. Most additive risk in the 8q24 region."
        ),
        "reference": "PRACTICAL Consortium 2018 + Schumacher 2018",
        "population_lifetime_risk": "~12% lifetime, much of which is indolent",
        "high_tier_implication": (
            "Top decile ~2.5× risk for clinically significant disease. "
            "Discuss informed PSA screening starting age 40–45 with urologist."
        ),
        "variants": [
            {"rsid": "rs1447295",   "effect_allele": "A", "log_or": 0.262, "af": 0.10, "gene": "8q24"},   # OR 1.3
            {"rsid": "rs6983267",   "effect_allele": "G", "log_or": 0.182, "af": 0.50, "gene": "8q24"},
            {"rsid": "rs10505477",  "effect_allele": "A", "log_or": 0.140, "af": 0.50, "gene": "CASC8/8q24"},
            {"rsid": "rs10993994",  "effect_allele": "T", "log_or": 0.182, "af": 0.39, "gene": "MSMB"},
            {"rsid": "rs2735839",   "effect_allele": "G", "log_or": 0.095, "af": 0.85, "gene": "KLK3 (PSA)"},
            {"rsid": "rs401681",    "effect_allele": "C", "log_or": 0.070, "af": 0.55, "gene": "CLPTM1L/TERT"},
            {"rsid": "rs2736098",   "effect_allele": "A", "log_or": 0.075, "af": 0.30, "gene": "TERT"},
            {"rsid": "rs10086908",  "effect_allele": "T", "log_or": 0.080, "af": 0.40, "gene": "8q24"},   # placeholder unused
        ],
    },

    "Alzheimer's Disease (Late-onset)": {
        "trait_short": "LOAD",
        "description": (
            "Lifetime risk of late-onset Alzheimer's disease (age 65+). "
            "Dominated by APOE-ε4 dosage; modulated by the polygenic "
            "background from IGAP-replicated loci."
        ),
        "reference": "IGAP 2019 + Lambert 2013 + Kunkle 2019",
        "population_lifetime_risk": "~10% by age 80, ~30% by age 90",
        "high_tier_implication": (
            "APOE-ε4 homozygotes have ~8–12× AD risk. High polygenic-plus-"
            "APOE risk warrants aggressive prevention: aerobic exercise, "
            "MIND diet, sleep, BP control, hearing care, social engagement."
        ),
        "variants": [
            # APOE is contributed via rs429358 / rs7412 — they appear here as
            # separate weights; the actual ε4 dosage effect is captured by
            # rs429358's strong log-OR.
            {"rsid": "rs429358",    "effect_allele": "C", "log_or": 1.099, "af": 0.14, "gene": "APOE ε4"},  # OR ~3
            {"rsid": "rs7412",      "effect_allele": "T", "log_or": -0.511, "af": 0.08, "gene": "APOE ε2"}, # protective
            {"rsid": "rs75932628",  "effect_allele": "T", "log_or": 1.099, "af": 0.003, "gene": "TREM2 R47H"}, # OR ~3
            {"rsid": "rs744373",    "effect_allele": "C", "log_or": 0.140, "af": 0.27, "gene": "BIN1"},
            {"rsid": "rs3851179",   "effect_allele": "T", "log_or": 0.095, "af": 0.36, "gene": "PICALM"},
            {"rsid": "rs10792832",  "effect_allele": "A", "log_or": 0.095, "af": 0.36, "gene": "PICALM"},
            {"rsid": "rs6656401",   "effect_allele": "A", "log_or": 0.140, "af": 0.20, "gene": "CR1"},
            {"rsid": "rs9331896",   "effect_allele": "T", "log_or": 0.095, "af": 0.38, "gene": "CLU"},
            {"rsid": "rs11136000",  "effect_allele": "T", "log_or": 0.075, "af": 0.62, "gene": "CLU"},
            {"rsid": "rs3764650",   "effect_allele": "G", "log_or": 0.182, "af": 0.10, "gene": "ABCA7"},
            {"rsid": "rs610932",    "effect_allele": "A", "log_or": 0.060, "af": 0.40, "gene": "MS4A6A"},
            {"rsid": "rs11218343",  "effect_allele": "C", "log_or": -0.182, "af": 0.04, "gene": "SORL1"},  # protective
            {"rsid": "rs10498633",  "effect_allele": "T", "log_or": 0.080, "af": 0.20, "gene": "SLC24A4"},
            {"rsid": "rs1476679",   "effect_allele": "T", "log_or": 0.060, "af": 0.30, "gene": "ZCWPW1"},
            {"rsid": "rs28834970",  "effect_allele": "C", "log_or": 0.090, "af": 0.36, "gene": "PTK2B"},
            {"rsid": "rs7274581",   "effect_allele": "C", "log_or": 0.060, "af": 0.08, "gene": "CASS4"},
            {"rsid": "rs9275595",   "effect_allele": "C", "log_or": 0.060, "af": 0.20, "gene": "HLA-DRB1"},
            {"rsid": "rs6857",      "effect_allele": "T", "log_or": 0.500, "af": 0.18, "gene": "TOMM40"},  # tags APOE haplotype
        ],
    },

    "Atrial Fibrillation": {
        "trait_short": "AFib",
        "description": (
            "Lifetime risk of atrial fibrillation (any-cause). The PITX2 "
            "locus on 4q25 dominates common-variant heritability."
        ),
        "reference": "AFGen 2017 + Roselli 2018",
        "population_lifetime_risk": "~25% lifetime over age 40",
        "high_tier_implication": (
            "Top decile ~3× AF risk. Aggressive risk-factor management: "
            "BP <130/80, body weight, sleep apnea screening, moderate "
            "alcohol (<7 drinks/wk), avoid extreme endurance loads."
        ),
        "variants": [
            {"rsid": "rs2200733",  "effect_allele": "T", "log_or": 0.531, "af": 0.12, "gene": "PITX2 (4q25)"},  # OR 1.7
            {"rsid": "rs2106261",  "effect_allele": "T", "log_or": 0.182, "af": 0.18, "gene": "ZFHX3"},
            {"rsid": "rs10033464", "effect_allele": "T", "log_or": 0.140, "af": 0.13, "gene": "KCNN3"},
        ],
    },

    "BMI / Obesity Tendency": {
        "trait_short": "BMI",
        "description": (
            "Continuous additive contribution to BMI from common variants. "
            "Reported in PRS standard deviations rather than disease risk; "
            "informative for weight-management planning."
        ),
        "reference": "GIANT Consortium 2018 + Locke 2015",
        "population_lifetime_risk": "Mean BMI ~26 in US adults; obesity ~40%",
        "high_tier_implication": (
            "Top decile genetic BMI tendency carries ~2 BMI-units of additive "
            "risk vs. average. Lifestyle (≥250 min/wk activity, higher-"
            "protein meals) overcomes most genetic tendency."
        ),
        "variants": [
            {"rsid": "rs9939609",  "effect_allele": "A", "log_or": 0.166, "af": 0.42, "gene": "FTO"},
            {"rsid": "rs8050136",  "effect_allele": "A", "log_or": 0.149, "af": 0.42, "gene": "FTO"},
            {"rsid": "rs17782313", "effect_allele": "C", "log_or": 0.110, "af": 0.24, "gene": "MC4R"},
            {"rsid": "rs2229616",  "effect_allele": "A", "log_or": -0.140, "af": 0.02, "gene": "MC4R V103I"}, # protective
            {"rsid": "rs7799039",  "effect_allele": "A", "log_or": 0.060, "af": 0.50, "gene": "LEP"},
            {"rsid": "rs1137101",  "effect_allele": "G", "log_or": 0.060, "af": 0.55, "gene": "LEPR Q223R"},
            {"rsid": "rs5082",     "effect_allele": "C", "log_or": 0.085, "af": 0.25, "gene": "APOA2"},
            {"rsid": "rs6232",     "effect_allele": "G", "log_or": 0.080, "af": 0.07, "gene": "PCSK1"},
        ],
    },
}


# Coverage thresholds for these small curated panels. Because each panel has
# only ~15-20 variants and the expected variance shrinks as variants drop out,
# a score built on a small fraction of the panel is unreliable.
PRS_MIN_USED = 3           # absolute floor of typed variants for any score
PRS_LOW_COVERAGE_PCT = 50  # below this the score is downgraded to "low"
PRS_HIGH_COVERAGE_PCT = 80  # at/above this coverage the score is "high"


def _prs_confidence(n_used: int, n_total: int) -> tuple:
    """Map curated-panel coverage to an explicit confidence level + note."""
    pct = 100.0 * n_used / max(n_total, 1)
    base = f"{n_used} of {n_total} panel variants typed ({pct:.0f}%)"
    if n_used < PRS_MIN_USED or pct < PRS_LOW_COVERAGE_PCT:
        return "low", (
            f"{base}. Too few variants for a reliable percentile — the score is "
            "shown for transparency only and should not be interpreted."
        )
    if pct < PRS_HIGH_COVERAGE_PCT:
        return "moderate", f"{base}. Missing variants add uncertainty to the percentile."
    return "high", f"{base}."


# ─── Core scoring ─────────────────────────────────────────────────────────────
def _dosage(genotype: str, effect_allele: str) -> Optional[int]:
    if not genotype:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--"):
        return None
    if len(gt) != 2:
        return None
    return gt.count(effect_allele.upper())


def _score_panel(snps_df: pd.DataFrame, panel: Dict, sex: Optional[str] = None) -> Dict:
    """Compute z-score / percentile / tier for a single panel."""
    applies_to = panel.get("applies_to")
    if applies_to == "female" and sex == "male":
        return {"status": "not_applicable", "confidence": "n/a",
                "reason": "Female-specific risk score; not applicable."}
    if applies_to == "male" and sex == "female":
        return {"status": "not_applicable", "confidence": "n/a",
                "reason": "Male-specific risk score; not applicable."}

    raw_score = 0.0
    expected_mean = 0.0
    expected_var = 0.0
    used: List[Dict] = []
    missing: List[Dict] = []

    for v in panel["variants"]:
        rsid = v["rsid"]
        if rsid not in snps_df.index:
            missing.append(v)
            continue
        row = snps_df.loc[rsid]
        gt = row.get("genotype")
        dose = _dosage(gt, v["effect_allele"])
        if dose is None:
            missing.append(v)
            continue
        beta = v["log_or"]
        af = v["af"]
        raw_score += beta * dose
        expected_mean += beta * 2.0 * af
        expected_var += (beta ** 2) * 2.0 * af * (1.0 - af)
        used.append({**v, "dosage": dose, "genotype": str(gt).upper()})

    n_total = len(panel["variants"])
    if not used or expected_var <= 0 or len(used) < PRS_MIN_USED:
        return {
            "status": "insufficient_data",
            "reason": (
                f"Only {len(used)} of {n_total} panel variants "
                f"were typed on this chip (min {PRS_MIN_USED}) — not enough for a "
                "score."
            ),
            "confidence": "none",
            "n_used": len(used),
            "n_expected": n_total,
            "callability": round(100.0 * len(used) / n_total, 1),
            "used": used,
            "missing": missing,
        }

    confidence, confidence_note = _prs_confidence(len(used), n_total)

    z_score = (raw_score - expected_mean) / sqrt(expected_var)
    percentile = _norm_cdf(z_score) * 100.0

    if percentile >= 95:
        tier = "High"
        tier_class = "tier-high"
    elif percentile >= 80:
        tier = "Elevated"
        tier_class = "tier-elevated"
    elif percentile >= 20:
        tier = "Average"
        tier_class = "tier-average"
    elif percentile >= 5:
        tier = "Below Average"
        tier_class = "tier-below"
    else:
        tier = "Low"
        tier_class = "tier-low"

    return {
        "status": "computed",
        "raw_score": round(raw_score, 4),
        "expected_mean": round(expected_mean, 4),
        "z_score": round(z_score, 3),
        "percentile": round(percentile, 1),
        "tier": tier,
        "tier_class": tier_class,
        "confidence": confidence,
        "confidence_note": confidence_note,
        "n_used": len(used),
        "n_expected": n_total,
        "used": used,
        "missing": missing,
        "callability": round(100.0 * len(used) / n_total, 1),
    }


def analyze_polygenic_scores(snps_df: pd.DataFrame, sex: Optional[str] = None) -> Dict:
    """Run all PRS panels. `sex` may be 'male', 'female', or None (auto)."""
    if sex is None:
        # Heuristic: count Y-chrom SNPs with non-trivial genotypes
        y = snps_df[snps_df.get("chrom") == "Y"] if "chrom" in snps_df.columns else pd.DataFrame()
        sex = "male" if len(y) > 100 else "female"

    panels: Dict[str, Dict] = {}
    for name, panel in PRS_PANELS.items():
        result = _score_panel(snps_df, panel, sex=sex)
        panels[name] = {
            "trait_short": panel.get("trait_short", ""),
            "description": panel["description"],
            "reference": panel["reference"],
            "population_lifetime_risk": panel.get("population_lifetime_risk", "—"),
            "high_tier_implication": panel.get("high_tier_implication", ""),
            "result": result,
        }

    # Headline findings: panels with tier Elevated or High — but only when
    # coverage is good enough to trust the percentile. A "High" tier built on a
    # handful of typed variants is not a headline finding.
    headline = [
        (name, p) for name, p in panels.items()
        if p["result"].get("tier") in ("Elevated", "High")
        and p["result"].get("confidence") in ("moderate", "high")
    ]
    return {
        "inferred_sex": sex,
        "panels": panels,
        "headline_findings": headline,
    }
