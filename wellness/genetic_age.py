"""
Genetic Longevity / Biological Age Proxy
========================================

Computes a "genetic longevity score" from a curated set of longevity-,
telomere-, and skin-aging-associated common variants. Compared against an
assumed European reference distribution, expressed as a percentile and a
rough biological-age offset (years vs population average).

This is NOT a clinical biological-age measurement (DNAm clocks like
Horvath, GrimAge, PhenoAge are the gold standards there). It is a
*genetic propensity* score — your underlying tendency before lifestyle.
Lifestyle (exercise, Mediterranean diet, sleep, social engagement, BP/
lipid control, no smoking) has 5-10× the impact of these variants and
can fully override an unfavorable score.
"""

from __future__ import annotations

from math import erf, sqrt

import pandas as pd


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def _dose(snps_df, rsid, allele):
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--") or len(s) != 2:
        return None
    return s.count(allele.upper())


# ─── Longevity score panel ────────────────────────────────────────────────────
LONGEVITY_VARIANTS = [
    # (rsid, effect_allele, beta_sd_units, af, gene, direction)
    ("rs2802292", "G", 0.10, 0.32, "FOXO3", "longevity favorable"),
    ("rs7412",    "T", 0.05, 0.08, "APOE ε2", "longevity favorable"),
    ("rs429358",  "C", -0.10, 0.14, "APOE ε4", "longevity unfavorable"),
    ("rs5882",    "G", 0.04, 0.32, "CETP I405V", "longevity favorable"),
    ("rs1800795", "C", -0.03, 0.43, "IL6 -174G/C", "inflammaging modifier"),
    ("rs2736100", "C", 0.04, 0.51, "TERT", "telomere length"),
    ("rs9420907", "C", 0.03, 0.65, "OBFC1", "telomere length"),
    ("rs10757278","G", -0.05, 0.50, "9p21 / CDKN2B-AS1", "CV risk modifier"),
    ("rs1333049", "C", -0.04, 0.49, "9p21", "CV risk modifier"),
]

# ─── Skin aging panel ────────────────────────────────────────────────────────
SKIN_AGING_VARIANTS = [
    ("rs1799750", "INS", 0.06, 0.50, "MMP1 promoter"),
    ("rs1800012", "T", 0.04, 0.20, "COL1A1 Sp1"),
    ("rs1042522", "G", 0.03, 0.45, "TP53 P72R"),
]

# ─── Telomere panel (genetic telomere length proxy) ──────────────────────────
TELOMERE_VARIANTS = [
    ("rs2736100", "C", 0.10, 0.51, "TERT"),
    ("rs9420907", "C", 0.06, 0.65, "OBFC1"),
]


def _score_panel(snps_df, panel):
    raw = exp_mean = exp_var = 0.0
    n_used = 0
    used = []
    for rsid, ea, beta, af, *rest in panel:
        d = _dose(snps_df, rsid, ea)
        if d is None:
            continue
        raw += beta * d
        exp_mean += beta * 2.0 * af
        exp_var += (beta ** 2) * 2.0 * af * (1.0 - af)
        n_used += 1
        used.append({"rsid": rsid, "dose": d, "beta": beta,
                     "gene": rest[0] if rest else ""})
    if n_used == 0 or exp_var <= 0:
        return None
    z = (raw - exp_mean) / sqrt(exp_var)
    pct = _norm_cdf(z) * 100
    return {"z": round(z, 3), "percentile": round(pct, 1),
            "n_used": n_used, "n_total": len(panel),
            "used": used}


def analyze_genetic_age(snps_df: pd.DataFrame) -> dict:
    longevity = _score_panel(snps_df, LONGEVITY_VARIANTS)
    skin = _score_panel(snps_df, SKIN_AGING_VARIANTS)
    telomere = _score_panel(snps_df, TELOMERE_VARIANTS)

    if not longevity:
        return {"available": False,
                "reason": "Insufficient longevity-panel coverage."}

    # Translate the longevity Z into a rough "years vs average" offset.
    # Calibration: each Z of 1.0 ≈ 2-3 years of difference at population level
    # (this is a rough mapping derived from Timmers 2019 parental lifespan MR;
    # an honest first-order proxy, not a clinical clock).
    years_offset = round(longevity["z"] * 2.5, 1)
    direction = "longer-lived" if years_offset > 0 else "shorter-lived"
    if abs(years_offset) < 0.5:
        narrative = (
            "Your genetic longevity profile is approximately average for European "
            f"reference (percentile {longevity['percentile']:.0f})."
        )
    else:
        narrative = (
            f"Your genetic longevity profile is in the {longevity['percentile']:.0f}th "
            f"percentile — suggesting a roughly {abs(years_offset):.1f}-year "
            f"{direction} genetic tendency vs the European population mean. "
            "Lifestyle factors typically dominate this — exercise, Mediterranean "
            "diet, sleep, BP/lipid control, social engagement add far more years "
            "than the genetic spread."
        )

    return {
        "available": True,
        "longevity": longevity,
        "longevity_years_offset": years_offset,
        "longevity_direction": direction,
        "narrative": narrative,
        "skin_aging": skin,
        "telomere": telomere,
        "disclaimer": (
            "This is a genetic propensity score, not a clinical biological age. "
            "DNAm-based clocks (Horvath, GrimAge, PhenoAge) measure actual "
            "epigenetic aging — those require methylation data not present here."
        ),
    }
