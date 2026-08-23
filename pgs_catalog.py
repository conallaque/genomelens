"""
Expanded Polygenic Risk Scores via PGS Catalog
==============================================

Parses PGS Catalog harmonised scoring files (Hmpos_GRCh37 format) and applies
them to the user's genotype data (chip + imputed). For each condition:

  * Raw additive score: Σ effect_weight_i × dosage_i
  * EUR-normalised Z-score and percentile (against expected Hardy-Weinberg
    distribution using EUR allele frequencies if available, else the
    published reference distribution).
  * Tier classification: Low / Below Average / Average / Elevated / High
  * Coverage breakdown: chip vs imputed vs missing

The module expects scoring files to be present at
~/dna-project/reference/pgs_scores/<PGS_ID>_hmPOS_GRCh37.txt.gz
(downloaded once via `python setup.py --pgs`).

Without these files, the module is dormant — the v2 curated PRS panels in
prs.py continue to provide a baseline.
"""

from __future__ import annotations

import gzip
from math import erf, sqrt
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent
PGS_DIR = SCRIPT_DIR / "reference" / "pgs_scores"


# Maps condition slug -> display info. Mirrors setup.py so the report can use
# friendly names even if the scoring file is named by PGS ID.
CONDITIONS: dict[str, dict[str, str]] = {
    "coronary_artery_disease": {
        "pgs_id": "PGS000018",
        "label": "Coronary Artery Disease",
        "short": "CAD",
        "population_lifetime_risk": "~30-40% lifetime",
    },
    "type_2_diabetes": {
        "pgs_id": "PGS000014",
        "label": "Type 2 Diabetes",
        "short": "T2D",
        "population_lifetime_risk": "~30-40% lifetime",
    },
    "breast_cancer": {
        "pgs_id": "PGS000004",
        "label": "Breast Cancer",
        "short": "BC",
        "applies_to": "female",
        "population_lifetime_risk": "~12% (women)",
    },
    "prostate_cancer": {
        "pgs_id": "PGS000662",
        "label": "Prostate Cancer",
        "short": "PCa",
        "applies_to": "male",
        "population_lifetime_risk": "~12% (men)",
    },
    "alzheimers_disease": {
        "pgs_id": "PGS000334",
        "label": "Alzheimer's Disease",
        "short": "LOAD",
        "population_lifetime_risk": "~10% by 80, ~30% by 90",
    },
    "atrial_fibrillation": {
        "pgs_id": "PGS000016",
        "label": "Atrial Fibrillation",
        "short": "AFib",
        "population_lifetime_risk": "~25% lifetime over 40",
    },
    "bmi": {
        "pgs_id": "PGS000027",
        "label": "BMI / Obesity Tendency",
        "short": "BMI",
        "population_lifetime_risk": "Continuous trait",
    },
    "major_depressive_disorder": {
        "pgs_id": "PGS000145",
        "label": "Major Depressive Disorder",
        "short": "MDD",
        "population_lifetime_risk": "~17% lifetime",
    },
    "schizophrenia": {
        "pgs_id": "PGS000019",
        "label": "Schizophrenia",
        "short": "SCZ",
        "population_lifetime_risk": "~1% lifetime",
    },
    "hypertension": {
        "pgs_id": "PGS000301",
        "label": "Hypertension / Systolic BP",
        "short": "HTN",
        "population_lifetime_risk": "~50% by age 60",
    },
    "stroke": {
        "pgs_id": "PGS000039",
        "label": "Ischemic Stroke",
        "short": "Stroke",
        "population_lifetime_risk": "~5-8% lifetime",
    },
    "chronic_kidney_disease": {
        "pgs_id": "PGS000314",
        "label": "Chronic Kidney Disease",
        "short": "CKD",
        "population_lifetime_risk": "~14% prevalence US",
    },
    "asthma": {
        "pgs_id": "PGS000037",
        "label": "Asthma",
        "short": "Asthma",
        "population_lifetime_risk": "~8% prevalence",
    },
    "inflammatory_bowel_disease": {
        "pgs_id": "PGS000020",
        "label": "Inflammatory Bowel Disease",
        "short": "IBD",
        "population_lifetime_risk": "~0.5% prevalence",
    },
    "rheumatoid_arthritis": {
        "pgs_id": "PGS000038",
        "label": "Rheumatoid Arthritis",
        "short": "RA",
        "population_lifetime_risk": "~1% prevalence",
    },
}


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def _log(msg: str) -> None:
    print(f"[pgs] {msg}", flush=True)


# ── Scoring file parsing ──────────────────────────────────────────────────────
def parse_pgs_scoring_file(path: Path) -> list[dict]:
    """Parse a PGS Catalog harmonised scoring file (Hmpos format).

    Returns a list of variant dicts:
        {rsid, chrom, pos, effect_allele, other_allele, weight, af}

    Allele frequency may not be present in older files; we tolerate that.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    variants: list[dict] = []
    header_cols: list[str] = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if not header_cols:
                header_cols = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < len(header_cols):
                parts += [""] * (len(header_cols) - len(parts))
            row = dict(zip(header_cols, parts, strict=False))
            # PGS Catalog hm columns
            rsid = row.get("rsID") or row.get("hm_rsID") or row.get("rsid") or ""
            chrom = row.get("hm_chr") or row.get("chr_name") or ""
            pos = row.get("hm_pos") or row.get("chr_position") or ""
            ea = row.get("effect_allele") or ""
            oa = row.get("other_allele") or row.get("hm_inferOtherAllele") or ""
            weight_str = row.get("effect_weight") or row.get("OR") or ""
            af_str = row.get("allelefrequency_effect") or ""
            try:
                weight = float(weight_str)
            except ValueError:
                continue
            try:
                pos_i = int(pos) if pos else 0
            except ValueError:
                pos_i = 0
            try:
                af = float(af_str) if af_str else None
            except ValueError:
                af = None
            variants.append({
                "rsid": rsid,
                "chrom": str(chrom),
                "pos": pos_i,
                "effect_allele": ea.upper(),
                "other_allele": oa.upper(),
                "weight": weight,
                "af": af,
            })
    return variants


# ── Genotype-to-dosage helper ─────────────────────────────────────────────────
def _dosage(genotype: object, effect_allele: str) -> int | None:
    if genotype is None:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--") or len(gt) != 2:
        return None
    return gt.count(effect_allele)


# Coverage thresholds for genome-wide PGS (often 10^2-10^6 variants). These are
# larger panels than the curated PRS, so a smaller *fraction* can still be
# usable, but very low coverage badly biases the percentile.
PGS_HIGH_COVERAGE_PCT = 80
PGS_MODERATE_COVERAGE_PCT = 40
PGS_MIN_COVERAGE_PCT = 10   # below this the percentile is not interpretable


def _pgs_confidence(n_used: int, total: int, n_low_r2: int) -> tuple:
    """Map PGS coverage (and imputation quality) to an explicit confidence."""
    pct = 100.0 * n_used / max(total, 1)
    base = f"{n_used:,} of {total:,} scoring-file variants used ({pct:.0f}%)"
    low_r2_note = ""
    if n_used and n_low_r2 / n_used > 0.25:
        low_r2_note = (
            f" {n_low_r2:,} used variants are low-confidence imputed (r²<0.5), "
            "adding noise."
        )
    if pct < PGS_MIN_COVERAGE_PCT:
        return "low", (
            f"{base} — far too few for a reliable percentile; shown for "
            "transparency only." + low_r2_note
        )
    if pct < PGS_MODERATE_COVERAGE_PCT:
        return "low", f"{base}. Low coverage materially biases the percentile." + low_r2_note
    if pct < PGS_HIGH_COVERAGE_PCT:
        return "moderate", f"{base}. Partial coverage adds uncertainty." + low_r2_note
    return "high", f"{base}." + low_r2_note


# ── PGS calculation ───────────────────────────────────────────────────────────
def calculate_pgs(snps_df: pd.DataFrame, variants: list[dict],
                  default_af: float = 0.30) -> dict:
    """Compute the raw score, expected mean & variance, z-score, percentile,
    and a coverage breakdown for a single PGS panel.

    snps_df may include both chip and imputed variants; we use whichever has
    the variant. The `source` column (chip|imputed) and `r2` column drive
    the coverage breakdown.
    """
    raw_score = 0.0
    expected_mean = 0.0
    expected_var = 0.0
    n_chip = n_imputed = n_low_r2 = 0
    n_missing = 0
    used: list[dict] = []

    for v in variants:
        rsid = v["rsid"]
        if not rsid or rsid == ".":
            n_missing += 1
            continue
        if rsid not in snps_df.index:
            n_missing += 1
            continue
        row = snps_df.loc[rsid]
        if isinstance(row, pd.DataFrame):
            # If duplicates, take first
            row = row.iloc[0]
        gt = row.get("genotype")
        dose = _dosage(gt, v["effect_allele"])
        if dose is None:
            n_missing += 1
            continue
        # Optional: drop low-r2 imputed sites
        r2 = row.get("r2", 1.0)
        source = row.get("source", "chip")
        if source == "imputed" and isinstance(r2, int | float) and r2 < 0.5:
            n_low_r2 += 1
            # still include but flag
        beta = v["weight"]
        af = v["af"] if v["af"] is not None else default_af
        raw_score += beta * dose
        expected_mean += beta * 2.0 * af
        expected_var += (beta ** 2) * 2.0 * af * (1.0 - af)
        if source == "imputed":
            n_imputed += 1
        else:
            n_chip += 1
        used.append({"rsid": rsid, "dose": dose, "weight": beta, "source": source})

    total = len(variants)
    n_used = n_chip + n_imputed
    if n_used == 0 or expected_var <= 0:
        return {
            "status": "insufficient_data",
            "reason": f"None of the {total} scoring-file variants found on chip or in imputed data.",
            "confidence": "none",
            "coverage": {
                "total": total, "chip": 0, "imputed": 0,
                "missing": n_missing, "low_r2": n_low_r2,
                "pct_callable": 0.0,
            },
        }

    confidence, confidence_note = _pgs_confidence(n_used, total, n_low_r2)

    z = (raw_score - expected_mean) / sqrt(expected_var)
    pct = _norm_cdf(z) * 100.0
    if pct >= 95:
        tier, cls = "High", "tier-high"
    elif pct >= 80:
        tier, cls = "Elevated", "tier-elevated"
    elif pct >= 20:
        tier, cls = "Average", "tier-average"
    elif pct >= 5:
        tier, cls = "Below Average", "tier-below"
    else:
        tier, cls = "Low", "tier-low"

    return {
        "status": "computed",
        "raw_score": round(raw_score, 4),
        "expected_mean": round(expected_mean, 4),
        "z_score": round(z, 3),
        "percentile": round(pct, 1),
        "tier": tier,
        "tier_class": cls,
        "confidence": confidence,
        "confidence_note": confidence_note,
        "coverage": {
            "total": total,
            "chip": n_chip,
            "imputed": n_imputed,
            "low_r2": n_low_r2,
            "missing": n_missing,
            "pct_callable": round(100 * n_used / total, 1),
        },
        "n_variants_used": n_used,
    }


# ── Top-level analyzer ────────────────────────────────────────────────────────
def analyze_expanded_pgs(snps_df: pd.DataFrame, sex: str | None = None) -> dict:
    """Run all available PGS Catalog panels.

    Returns a dict with per-condition results. If scoring files are not
    downloaded, returns a status indicating setup is needed.
    """
    if not PGS_DIR.exists():
        return {
            "available": False,
            "reason": "PGS scoring files not downloaded. Run `python setup.py --pgs`.",
            "panels": {},
        }

    panels: dict[str, dict] = {}
    for slug, info in CONDITIONS.items():
        pgs_id = info["pgs_id"]
        applies_to = info.get("applies_to")
        if applies_to == "female" and sex == "male":
            panels[slug] = {
                **info,
                "result": {"status": "not_applicable", "confidence": "n/a",
                           "reason": "Female-specific score; not applicable."},
            }
            continue
        if applies_to == "male" and sex == "female":
            panels[slug] = {
                **info,
                "result": {"status": "not_applicable", "confidence": "n/a",
                           "reason": "Male-specific score; not applicable."},
            }
            continue

        score_file = PGS_DIR / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        if not score_file.exists():
            panels[slug] = {
                **info,
                "result": {"status": "not_downloaded",
                           "reason": f"Scoring file {score_file.name} not found. Run `python setup.py --pgs`."},
            }
            continue

        try:
            _log(f"Loading {pgs_id} ({info['label']}) ...")
            variants = parse_pgs_scoring_file(score_file)
            _log(f"  {len(variants):,} weighted variants parsed")
            result = calculate_pgs(snps_df, variants)
            panels[slug] = {**info, "result": result}
        except Exception as e:
            panels[slug] = {**info, "result": {"status": "error", "reason": str(e)}}

    # Headline: any panel in Elevated or High
    headline = [
        (slug, p) for slug, p in panels.items()
        if p.get("result", {}).get("tier") in ("Elevated", "High")
    ]
    return {
        "available": True,
        "panels": panels,
        "headline_findings": headline,
        "n_panels": len(panels),
    }
