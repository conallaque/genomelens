"""
Mendelian Randomization Inference
=================================

Uses pre-defined exposure → outcome MR effect estimates from large
published studies, personalized by the user's own polygenic score for
the exposure trait. The output is "if your genetics predict X is shifted
by ΔSD, the literature-derived causal effect on outcome Y is approximately Z".

Conceptual basis:
    MR uses genetic variants as natural-experiment instrumental variables.
    Because variant assignment at conception is independent of confounders
    (a Mendelian randomisation), it provides causal — not merely
    correlational — estimates of an exposure's effect on an outcome.

Our implementation is a simplification: we take a published MR causal
estimate (β_MR per SD increase in exposure), score the user's exposure
PRS from phewas-style weights, and compute the projected outcome shift.

Limitations:
    * MR assumes the genetic instruments affect the outcome only through
      the exposure (no pleiotropy). Real published MR papers test this;
      our simplified module trusts the published estimates.
    * Population-level causal estimates may not translate one-to-one to
      individual phenotypes — these are projections, not predictions.
    * Educational only — not a clinical decision tool.
"""

from __future__ import annotations

from math import erf, sqrt

import pandas as pd


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> int | None:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--") or len(s) != 2:
        return None
    return s.count(allele.upper())


def _score_exposure(snps_df, exposure_snps):
    """Score the user's exposure PRS — return (z, callability)."""
    raw = exp_mean = exp_var = 0.0
    n_used = 0
    for rsid, ea, beta, af in exposure_snps:
        d = _dose(snps_df, rsid, ea)
        if d is None:
            continue
        raw += beta * d
        exp_mean += beta * 2.0 * af
        exp_var += (beta ** 2) * 2.0 * af * (1.0 - af)
        n_used += 1
    if n_used == 0 or exp_var <= 0:
        return None, 0
    z = (raw - exp_mean) / sqrt(exp_var)
    return z, n_used


# ─── MR causal-effect catalog ──────────────────────────────────────────────
# Each entry:
#   exposure: name + SNPs for scoring the user's exposure PRS
#   outcome: name + literature MR causal estimate (β_MR per SD increase in exposure)
#   units: human description of effect

MR_LIBRARY: list[dict] = [
    {
        "exposure": "Genetically predicted LDL cholesterol",
        "outcome": "Coronary Artery Disease",
        "mr_beta": 0.50,
        "mr_unit": "log-OR per 1 SD higher LDL",
        "reference": "Ference 2017 (Nature Reviews Cardiology)",
        "explanation": (
            "Decades of MR studies confirm that lifelong genetically elevated LDL "
            "causally raises CAD risk approximately log-linearly with cumulative "
            "exposure. This is the strongest piece of MR evidence in cardiology — "
            "and the foundation of lifelong-LDL-lowering strategies."
        ),
        "exposure_snps": [
            ("rs6511720", "T", -0.21, 0.10),
            ("rs646776", "C", 0.10, 0.78),
            ("rs693", "T", 0.10, 0.51),
            ("rs11591147", "T", -0.50, 0.02),
        ],
    },
    {
        "exposure": "Genetically predicted Lp(a)",
        "outcome": "Coronary Artery Disease",
        "mr_beta": 0.30,
        "mr_unit": "log-OR per 1 SD higher Lp(a)",
        "reference": "Burgess 2018 (JAMA Cardiology)",
        "explanation": (
            "Lp(a) is a residual cardiovascular risk factor that is largely "
            "unchanged by statins. MR studies provide strong evidence for a "
            "causal effect on CAD, stroke, and aortic valve stenosis. Therapies "
            "(siRNA, antisense oligonucleotides) that lower Lp(a) are in "
            "late-stage trials."
        ),
        "exposure_snps": [
            ("rs10455872", "G", 0.55, 0.07),
            ("rs3798220", "C", 0.50, 0.02),
        ],
    },
    {
        "exposure": "Genetically predicted BMI",
        "outcome": "Type 2 Diabetes",
        "mr_beta": 0.85,
        "mr_unit": "log-OR per 1 SD higher BMI",
        "reference": "Holmes 2014 (PLoS Genetics)",
        "explanation": (
            "MR strongly supports BMI as a causal driver of T2D, beyond mere "
            "correlation. Each SD of higher BMI (~4.5 kg/m²) raises T2D risk "
            "approximately 2.3-fold."
        ),
        "exposure_snps": [
            ("rs9939609", "A", 0.10, 0.42),
            ("rs1421085", "C", 0.10, 0.42),
            ("rs17782313", "C", 0.07, 0.24),
        ],
    },
    {
        "exposure": "Genetically predicted 25-OH Vitamin D",
        "outcome": "Multiple Sclerosis",
        "mr_beta": -0.50,
        "mr_unit": "log-OR per 1 SD higher vitamin D",
        "reference": "Mokry 2015 (PLoS Medicine)",
        "explanation": (
            "MR provides strong evidence that genetically elevated vitamin D "
            "lowers MS risk. The effect on other autoimmune diseases is weaker "
            "but suggestive."
        ),
        "exposure_snps": [
            ("rs2282679", "C", -0.13, 0.28),
            ("rs10741657", "G", -0.10, 0.40),
        ],
    },
    {
        "exposure": "Genetically predicted 25-OH Vitamin D",
        "outcome": "Bone Mineral Density",
        "mr_beta": 0.20,
        "mr_unit": "SD BMD per 1 SD higher vitamin D",
        "reference": "Trajanoska 2018",
        "explanation": (
            "MR supports a causal effect of vitamin D status on BMD, though "
            "smaller than the observational association suggests."
        ),
        "exposure_snps": [
            ("rs2282679", "C", -0.13, 0.28),
            ("rs10741657", "G", -0.10, 0.40),
        ],
    },
    {
        "exposure": "Genetically predicted serum testosterone (men)",
        "outcome": "Prostate Cancer Risk",
        "mr_beta": 0.20,
        "mr_unit": "log-OR per 1 SD higher T",
        "reference": "Mohammadi-Shemirani 2019",
        "applies_to": "male",
        "explanation": (
            "MR supports a modest causal effect of higher testosterone on "
            "prostate cancer risk. The observational correlation is widely "
            "discussed; MR suggests it is largely real but with smaller "
            "magnitude than naive correlation suggests."
        ),
        "exposure_snps": [
            ("rs727428", "G", 0.10, 0.50),
            ("rs1799941", "A", -0.10, 0.27),
            ("rs6259", "A", -0.08, 0.12),
        ],
    },
    {
        "exposure": "Genetically predicted coffee consumption",
        "outcome": "Coronary Artery Disease",
        "mr_beta": -0.05,
        "mr_unit": "log-OR per 1 SD more coffee",
        "reference": "Kwok 2016",
        "explanation": (
            "MR suggests modest cardio-protective effect of habitual coffee "
            "consumption — opposite of what some observational studies once "
            "suggested. Caffeine-anxiety and sleep-disruption side effects are "
            "separate concerns."
        ),
        "exposure_snps": [
            ("rs762551", "A", 0.10, 0.65),
        ],
    },
    {
        "exposure": "Genetically predicted alcohol consumption (ADH1B-driven)",
        "outcome": "Esophageal & Upper-GI Cancers",
        "mr_beta": 0.80,
        "mr_unit": "log-OR per 1 SD more alcohol",
        "reference": "Lewis 2018 + Yokoyama 2003",
        "explanation": (
            "MR using ADH1B and ALDH2 instruments provides strong causal "
            "evidence that even moderate alcohol intake elevates upper-GI cancer "
            "risk, especially in carriers of slow ALDH2 (ALDH2*2)."
        ),
        "exposure_snps": [
            ("rs1229984", "A", -0.20, 0.03),
            ("rs671", "A", -0.30, 0.002),
        ],
    },
    {
        "exposure": "Genetically predicted iron / ferritin",
        "outcome": "Cardiovascular Mortality (in iron-overload range)",
        "mr_beta": 0.10,
        "mr_unit": "log-OR per 1 SD higher iron (high range)",
        "reference": "Gill 2017",
        "explanation": (
            "Iron overload causally elevates cardiovascular mortality risk. "
            "Most relevant in HFE C282Y homozygotes; in the typical range "
            "iron status is more important to optimize downward than upward."
        ),
        "exposure_snps": [
            ("rs1800562", "A", 0.30, 0.06),
            ("rs1799945", "G", 0.10, 0.16),
        ],
    },
    {
        "exposure": "Genetically predicted SHBG",
        "outcome": "Type 2 Diabetes",
        "mr_beta": -0.40,
        "mr_unit": "log-OR per 1 SD higher SHBG",
        "reference": "Ding 2009",
        "explanation": (
            "MR supports a causal protective effect of higher SHBG on T2D, "
            "particularly in women. Mechanism involves altered free-hormone "
            "balance and hepatic insulin signaling."
        ),
        "exposure_snps": [
            ("rs1799941", "A", 0.20, 0.27),
            ("rs6259", "A", 0.16, 0.12),
        ],
    },
    {
        "exposure": "Genetically predicted CRP (chronic inflammation)",
        "outcome": "Coronary Artery Disease",
        "mr_beta": 0.05,
        "mr_unit": "log-OR per 1 SD higher CRP",
        "reference": "CRP CHD Coalition 2011",
        "explanation": (
            "MR studies suggest that the observational CRP-CAD correlation is "
            "largely NOT causal — most of the association is driven by confounding "
            "with metabolic risk factors. CRP is a marker, not a direct driver."
        ),
        "exposure_snps": [
            ("rs1800795", "C", 0.10, 0.43),
            ("rs2228145", "C", 0.06, 0.43),
        ],
    },
    {
        "exposure": "Genetically predicted BMI",
        "outcome": "Coronary Artery Disease",
        "mr_beta": 0.40,
        "mr_unit": "log-OR per 1 SD higher BMI",
        "reference": "Dale 2017",
        "explanation": (
            "MR confirms higher BMI causally raises CAD risk independent of "
            "lipids/BP. Each SD of higher BMI raises CAD risk ~1.5×."
        ),
        "exposure_snps": [
            ("rs9939609", "A", 0.10, 0.42),
            ("rs1421085", "C", 0.10, 0.42),
        ],
    },
]


def analyze_mr(snps_df: pd.DataFrame, sex: str | None = None) -> dict:
    """Compute personalized MR projections for each exposure-outcome pair."""
    findings: list[dict] = []
    for entry in MR_LIBRARY:
        if entry.get("applies_to") == "male" and sex == "female":
            continue
        if entry.get("applies_to") == "female" and sex == "male":
            continue
        z, n_used = _score_exposure(snps_df, entry["exposure_snps"])
        if z is None:
            findings.append({
                **entry,
                "exposure_z": None,
                "exposure_percentile": None,
                "outcome_shift_log_or": None,
                "outcome_relative_risk": None,
                "n_used": 0,
                "status": "insufficient_data",
            })
            continue
        # Projected outcome shift = mr_beta × exposure_z
        outcome_log_or = entry["mr_beta"] * z
        rr = round(2.71828 ** outcome_log_or, 2)
        findings.append({
            **entry,
            "exposure_z": round(z, 2),
            "exposure_percentile": round(_norm_cdf(z) * 100, 1),
            "outcome_shift_log_or": round(outcome_log_or, 3),
            "outcome_relative_risk": rr,
            "n_used": n_used,
            "status": "ok",
        })

    return {
        "findings": findings,
        "n_total": len(findings),
        "n_computed": sum(1 for f in findings if f["status"] == "ok"),
    }
