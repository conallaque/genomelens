"""
Expanded Polygenic Risk Scores via PGS Catalog
==============================================

Parses PGS Catalog harmonised scoring files (Hmpos_GRCh37 format) and applies
them to the user's genotype data (chip + imputed). For each condition:

  * Raw additive score: Σ effect_weight_i × dosage_i
  * EUR-normalized Z-score and percentile (against expected Hardy-Weinberg
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

# .parent.parent: moved into the risk package, but reference/ stays at the
# repository root. Left at .parent these resolve inside risk/ and every
# consumer degrades quietly on a missing directory — so nothing raises
# and the screen silently reports no findings.
SCRIPT_DIR = Path(__file__).parent.parent
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
            # EUR-SPECIFIC FREQUENCY COLUMNS COUNT TOO. This read only the
            # plain `allelefrequency_effect`, so a file publishing
            # `allelefrequency_effect_European` — which PGS000662 does for all
            # 269 of its variants — parsed as "no frequency available" and every
            # variant fell back to the 0.30 default. Correcting only the
            # frequency source on that panel moves it from z +5.37 (100.0th) to
            # z -1.20 (11.6th). The module docstring already promises EUR
            # frequencies "if available"; they were available and unread.
            af_str = row.get("allelefrequency_effect") or ""
            if not af_str:
                for _k, _v in row.items():
                    if (_k.startswith("allelefrequency_effect")
                            and "european" in _k.lower() and _v):
                        af_str = _v
                        break
            if not af_str:
                for _k, _v in row.items():
                    if _k.startswith("allelefrequency_effect") and _v:
                        af_str = _v
                        break
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

# Ratio of carried dose to frequency-predicted dose above which the callable
# set is judged non-random. See the gate in calculate_pgs.
PGS_INFORMATIVE_MISSINGNESS_RATIO = 1.5


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
    # Dose actually carried vs dose the frequencies predict — see the
    # informative-missingness gate below.
    observed_dose_sum = 0.0
    expected_dose_sum = 0.0
    n_no_rsid = 0
    n_chip = n_imputed = n_low_r2 = 0
    n_missing = 0
    used: list[dict] = []

    for v in variants:
        rsid = v["rsid"]
        if not rsid or rsid == ".":
            n_missing += 1
            n_no_rsid += 1
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
        observed_dose_sum += dose
        expected_dose_sum += 2.0 * af
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
            # BLAME THE RIGHT THING. A scoring file with no rsID column — PGS000004
            # publishes coordinates and frequencies for all 313 of its variants and
            # no rsIDs — yields an empty rsid for every row, and lookup here is
            # rsID-only, so all 313 counted as "missing" and the message told the
            # reader their genome lacked the variants. It did not; this module
            # never looked them up. Coordinate matching would need a GRCh37->38
            # liftover the scoring files do not carry, so this states the real
            # cause rather than inventing a join.
            "reason": (
                f"None of the {total} scoring-file variants could be looked up."
                + ("  This scoring file publishes no rsIDs, and matching here is "
                   "rsID-only — the variants were not searched for, rather than "
                   "searched for and absent." if n_no_rsid == total and total
                   else "  Not found on chip or in imputed data.")
            ),
            "confidence": "none",
            "coverage": {
                "total": total, "chip": 0, "imputed": 0,
                "missing": n_missing, "low_r2": n_low_r2,
                "pct_callable": 0.0,
            },
        }

    confidence, confidence_note = _pgs_confidence(n_used, total, n_low_r2)

    # SUPPRESSION WAS GATED ON ZERO, NOT ON ENOUGH. The only path to
    # "insufficient_data" above is `n_used == 0` — literally none of the
    # scoring-file variants found. One variant out of 1.7 million fell through
    # and rendered a full percentile and a tier. PGS_MIN_COVERAGE_PCT existed
    # and said "below this the percentile is not interpretable", but nothing
    # ever read it except the confidence label, so the module printed a
    # percentile alongside its own text calling that percentile unreliable.
    # On a real genome that produced "High, 100th percentile" for coronary
    # artery disease off 28% coverage, beside a curated panel calling the same
    # person Average.
    #
    # The rule now matches the wording: if coverage is low enough that this
    # module would label the score low-confidence, it does not report a
    # percentile at all. JUDGMENT CALL — no enforced threshold existed in
    # either scoring module (the curated PRS gate is an absolute floor of
    # PRS_MIN_USED variants, not a percentage), so this reuses the declared
    # PGS_MODERATE_COVERAGE_PCT boundary rather than inventing a new number.
    # INFORMATIVE MISSINGNESS. Dropping absent variants is correct only when
    # absence is random, which holds for a genotyping chip. It does NOT hold for
    # a whole-genome gVCF, where rsIDs are attached only to non-reference sites:
    # absence there means homozygous reference, not "not measured". The callable
    # set is then conditioned on carrying an alt allele, every retained variant
    # scores high, and the score is compared against a null built for a random
    # sample. Measured on a real 30x genome: of 502,690 callable scoring
    # variants, ZERO were 0/0, and mean carried dose ran 1.00-1.14 against the
    # 0.60 the frequencies predict. That is the whole reason five unrelated
    # traits pinned at exactly the 100th percentile.
    #
    # This is NOT correctable from the callable set. The right value for an
    # absent variant is 2 if the effect allele is the reference allele and 0
    # otherwise; the gVCF's hom-ref records carry no rsID to look up, and the
    # scoring files are GRCh37 against a GRCh38 callset, so coordinates do not
    # join either. Scoring absent variants as 0 just produces the mirror
    # artifact — the same panels pinned at the 0th percentile instead. So the
    # percentile is withheld rather than corrected.
    #
    # JUDGMENT CALL, flagged like the coverage threshold: 1.5x. Random
    # missingness puts this ratio at ~1.0 (a chip export does not trip it);
    # the affected panels here sit at 1.7-2.3.
    dose_ratio = (observed_dose_sum / expected_dose_sum
                  if expected_dose_sum > 0 else 0.0)
    if dose_ratio > PGS_INFORMATIVE_MISSINGNESS_RATIO:
        return {
            "status": "insufficient_data",
            "reason": (
                f"Carried dose is {dose_ratio:.2f}x what the allele "
                "frequencies predict, so the variants that were callable are "
                "not a random sample of the panel — on a whole-genome callset "
                "an absent variant means homozygous reference, not unmeasured. "
                "A percentile computed here would be biased upward and cannot "
                "be corrected from the callable variants alone."
            ),
            "confidence": "none",
            "informative_missingness": True,
            "dose_ratio": round(dose_ratio, 3),
            "coverage": {
                "total": total, "chip": n_chip, "imputed": n_imputed,
                "missing": n_missing, "low_r2": n_low_r2,
                "pct_callable": round(100.0 * n_used / max(total, 1), 1),
            },
        }

    pct_callable = 100.0 * n_used / max(total, 1)
    if pct_callable < PGS_MODERATE_COVERAGE_PCT:
        return {
            "status": "insufficient_data",
            "reason": (
                f"Only {n_used:,} of {total:,} scoring-file variants were "
                f"callable ({pct_callable:.1f}%, minimum "
                f"{PGS_MODERATE_COVERAGE_PCT}%) — not enough for an "
                "interpretable percentile."
            ),
            "confidence": "none",
            "coverage": {
                "total": total, "chip": n_chip, "imputed": n_imputed,
                "missing": n_missing, "low_r2": n_low_r2,
                "pct_callable": round(pct_callable, 1),
            },
        }

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
        # Recorded on scores that PASS the gate as well, so a reader can see how
        # close a surviving score sat to it. Inflammatory bowel disease clears
        # the bar and still reports the 98th percentile; without this number
        # there is no way to judge that from the outside.
        "dose_ratio": round(dose_ratio, 3),
        "coverage": {
            "total": total,
            "chip": n_chip,
            "imputed": n_imputed,
            "low_r2": n_low_r2,
            "missing": n_missing,
            "pct_callable": round(100 * n_used / total, 1),
        },
        "n_variants_used": n_used,
        # Per-variant provenance. This list was built and appended to for every
        # variant and then never placed in either return dict — a dead
        # computation. Every PGS card runs on thin coverage (9-23% callable) and
        # prints a z-score whose two inputs were withheld, so the percentile
        # could not be sanity-checked. Capped: the point is auditability, not a
        # dump of a whole scoring file.
        "used": used[:200],
        "n_used_returned": min(len(used), 200),
    }


# ── Top-level analyzer ────────────────────────────────────────────────────────
def analyze_expanded_pgs(snps_df: pd.DataFrame, sex: str | None = None,
                         vcf_path: str | None = None,
                         build: str | None = None) -> dict:
    """Run all available PGS Catalog panels.

    When `vcf_path` and `build` are supplied and scoring files exist for that
    build, panels are scored by COORDINATE against the callset. That is the
    correct join for a whole genome: an absent rsID there means the site is
    homozygous reference, not unmeasured, and dropping those sites conditions
    the score on carrying an alt allele. Chip exports keep the rsID path, where
    missingness genuinely is close to random.
    """
    if not PGS_DIR.exists():
        return {
            "available": False,
            "reason": "PGS scoring files not downloaded. Run `python setup.py --pgs`.",
            "panels": {},
        }

    # ── coordinate pre-pass ──────────────────────────────────────────────
    # One read of the callset for every panel at once. Positions are collected
    # across all panels first so the file is streamed a single time rather than
    # fifteen; the per-panel arrays then index into the shared result.
    coord_ok = False
    needed = a1 = a2 = None
    panel_variants: dict[str, list] = {}
    if vcf_path and build:
        try:
            import numpy as np
            want = (build or "").lower()
            keys: set = set()
            for slug, info in CONDITIONS.items():
                f = PGS_DIR / f"{info['pgs_id']}_hmPOS_{'GRCh38' if want == 'grch38' else 'GRCh37'}.txt.gz"
                if not f.exists() or pgs_file_build(f) != want:
                    continue
                vs = parse_pgs_scoring_file(f)
                panel_variants[slug] = vs
                for v in vs:
                    k = _pack(v["chrom"], v["pos"]) if v.get("pos") else -1
                    if k > 0:
                        keys.add(k)
            if keys:
                _log(f"Coordinate join: {len(keys):,} positions across "
                     f"{len(panel_variants)} panel(s), build {want}")
                needed = np.array(sorted(keys), dtype=np.int64)
                a1, a2 = build_coordinate_genotypes(vcf_path, needed)
                seen = int((a1 > 0).sum())
                _log(f"  {seen:,} of {len(needed):,} positions present in the "
                     f"callset ({100.0 * seen / max(len(needed), 1):.1f}%)")
                coord_ok = True
        except Exception as e:  # fall back to the rsID path
            _log(f"  WARNING: coordinate join unavailable ({e}); using rsIDs")
            coord_ok = False

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
            if coord_ok and slug in panel_variants:
                result = calculate_pgs_by_coordinate(
                    panel_variants[slug], needed, a1, a2)
            else:
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


# ── cross-module reconciliation: curated PRS vs PGS Catalog ──────────────────
# Two independent modules score the same seven conditions from the same genome
# and their answers were rendered side by side with nothing saying they
# disagreed. On a real whole genome the curated panel put this person at the
# 96th percentile for Alzheimer's (High) while the PGS Catalog panel put them
# at the 13th (Below Average) — the same trait, the same DNA, opposite ends of
# the distribution. Neither is necessarily wrong: they are different variant
# sets, different training cohorts and different coverage. But a reader given
# both numbers and no flag will assume they corroborate each other.
#
# Same posture as the economics reconciliation gate: FLAG, do not suppress and
# do not silently pick a winner. Deciding which score is right is a modeling
# question this function is not entitled to answer.
_PRS_TO_PGS: dict[str, str] = {
    "Coronary Artery Disease": "coronary_artery_disease",
    "Type 2 Diabetes": "type_2_diabetes",
    "Breast Cancer": "breast_cancer",
    "Prostate Cancer": "prostate_cancer",
    "Alzheimer's Disease (Late-onset)": "alzheimers_disease",
    "Atrial Fibrillation": "atrial_fibrillation",
    "BMI / Obesity Tendency": "bmi",
}

# A percentile gap this wide means the two scores disagree about which side of
# the distribution the person sits on, not merely by how much.
PRS_PGS_DIVERGENCE_PCT = 40.0


def reconcile_prs_pgs(prs_result: dict, pgs_result: dict) -> list[dict]:
    """Compare curated-PRS and PGS-Catalog scores for conditions scored by both.

    Returns one record per condition where both modules produced a percentile,
    with `divergent` set when they disagree sharply. Conditions where either
    side was suppressed for coverage are skipped rather than compared: a
    withheld score is not a disagreement.
    """
    prs_panels = (prs_result or {}).get("panels") or {}
    pgs_panels = (pgs_result or {}).get("panels") or {}
    out: list[dict] = []
    for prs_name, pgs_slug in _PRS_TO_PGS.items():
        a = (prs_panels.get(prs_name) or {}).get("result") or {}
        b = (pgs_panels.get(pgs_slug) or {}).get("result") or {}
        if a.get("status") != "computed" or b.get("status") != "computed":
            continue
        pa, pb = a.get("percentile"), b.get("percentile")
        if pa is None or pb is None:
            continue
        gap = abs(float(pa) - float(pb))
        out.append({
            "condition": prs_name,
            "pgs_id": (pgs_panels.get(pgs_slug) or {}).get("pgs_id", ""),
            "curated_percentile": pa, "curated_tier": a.get("tier"),
            "curated_confidence": a.get("confidence"),
            "catalog_percentile": pb, "catalog_tier": b.get("tier"),
            "catalog_confidence": b.get("confidence"),
            "percentile_gap": round(gap, 1),
            "tiers_agree": a.get("tier") == b.get("tier"),
            "divergent": gap >= PRS_PGS_DIVERGENCE_PCT,
        })
    out.sort(key=lambda r: -r["percentile_gap"])
    return out


def format_prs_pgs_reconciliation(rows: list[dict]) -> str:
    """One line per compared condition; divergent ones marked."""
    if not rows:
        return "No condition is scored by both the curated PRS and PGS Catalog modules."
    lines = []
    for r in rows:
        mark = "  <-- DIVERGENT" if r["divergent"] else ""
        lines.append(
            f"{r['condition']:34} curated {r['curated_percentile']:>5}th "
            f"({r['curated_tier']}) vs catalog {r['catalog_percentile']:>5}th "
            f"({r['catalog_tier']})  gap {r['percentile_gap']:>5}{mark}"
        )
    n = sum(1 for r in rows if r["divergent"])
    lines.append(f"{n} of {len(rows)} compared conditions diverge by "
                 f"{PRS_PGS_DIVERGENCE_PCT:.0f} percentile points or more.")
    return "\n".join(lines)


# ── coordinate-joined scoring (whole-genome callsets) ────────────────────────
# WHY THIS EXISTS. Lookup here was rsID-only. On a gVCF that is not a coverage
# limitation, it is a selection effect: annotators attach rsIDs only to
# non-reference sites, so every hom-reference call — 65% of records in the
# genome this was built against, and 0 of 261,114 of them carrying an rsID —
# was invisible. The scored set was therefore conditioned on carrying an alt
# allele, which is what pinned five unrelated traits at the 100th percentile.
#
# The information was never missing. It was unjoinable: hom-ref records have no
# rsID, and the shipped scoring files are GRCh37 against a GRCh38 callset, so
# coordinates did not line up either. Joining on coordinates in the SAME build
# fixes both at once — absence stops meaning "unmeasured" and starts meaning
# "homozygous reference", which is a dose, not a gap.
_ALLELE_CODE = {"A": 1, "C": 2, "G": 3, "T": 4}
_CODE_ALLELE = {1: "A", 2: "C", 3: "G", 4: "T"}


def _chrom_num(c: str) -> int:
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.isdigit():
        return int(c)
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c.upper(), 0)


def _pack(chrom: str, pos: int) -> int:
    n = _chrom_num(chrom)
    return n * 1_000_000_000 + int(pos) if n else -1


def build_coordinate_genotypes(vcf_path, needed_packed):
    """Read a VCF once and return allele codes at exactly the positions wanted.

    `needed_packed` is a sorted numpy int64 array of packed chrom/pos keys.
    Returns (a1, a2) uint8 arrays parallel to it; 0 means the position was not
    seen. Hom-reference records are kept — they are the whole point.
    """
    import gzip as _gz

    import numpy as np

    n = len(needed_packed)
    a1 = np.zeros(n, dtype=np.uint8)
    a2 = np.zeros(n, dtype=np.uint8)
    opener = _gz.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if len(f) < 10:
                continue
            ref, alt = f[3], f[4]
            if len(ref) != 1 or ref not in _ALLELE_CODE:
                continue
            key = _pack(f[0], f[1]) if f[1].isdigit() else -1
            if key < 0:
                continue
            i = np.searchsorted(needed_packed, key)
            if i >= n or needed_packed[i] != key:
                continue
            gt = f[9].split(":", 1)[0].replace("|", "/")
            parts = gt.split("/")
            if len(parts) == 1:
                parts = [parts[0], parts[0]]
            if len(parts) != 2 or "." in parts:
                continue
            alts = [x for x in alt.split(",")] if alt not in (".", "") else []
            codes = []
            ok = True
            for p in parts:
                try:
                    idx = int(p)
                except ValueError:
                    ok = False
                    break
                if idx == 0:
                    seq = ref
                elif 1 <= idx <= len(alts):
                    seq = alts[idx - 1]
                else:
                    ok = False
                    break
                if len(seq) != 1 or seq.upper() not in _ALLELE_CODE:
                    ok = False
                    break
                codes.append(_ALLELE_CODE[seq.upper()])
            if ok and len(codes) == 2:
                a1[i], a2[i] = codes[0], codes[1]
    return a1, a2


def pgs_file_build(path) -> str | None:
    """Which reference build a scoring filename declares, or None."""
    low = str(path).lower()
    if "grch38" in low or "hg38" in low:
        return "grch38"
    if "grch37" in low or "hg19" in low:
        return "grch37"
    return None


def calculate_pgs_by_coordinate(variants: list[dict], needed_packed,
                                a1, a2, default_af: float = 0.30) -> dict:
    """Score one panel from coordinate-joined genotypes.

    Absence now means "position not present in the callset" rather than
    "person carries no alt allele", so the informative-missingness gate that
    guards the rsID path is not needed here — hom-reference sites contribute
    their real dose (2 when the effect allele IS the reference allele, 0
    otherwise) instead of being dropped.
    """
    from math import sqrt

    import numpy as np

    raw = mean = var = 0.0
    n_used = n_missing = n_homref = 0
    obs = exp = 0.0
    for v in variants:
        key = _pack(v["chrom"], v["pos"]) if v.get("pos") else -1
        if key < 0:
            n_missing += 1
            continue
        i = np.searchsorted(needed_packed, key)
        if i >= len(needed_packed) or needed_packed[i] != key or not a1[i]:
            n_missing += 1
            continue
        g = _CODE_ALLELE.get(int(a1[i]), "") + _CODE_ALLELE.get(int(a2[i]), "")
        dose = _dosage(g, v["effect_allele"])
        if dose is None:
            n_missing += 1
            continue
        af = v["af"] if v["af"] is not None else default_af
        beta = v["weight"]
        raw += beta * dose
        mean += beta * 2.0 * af
        var += (beta ** 2) * 2.0 * af * (1.0 - af)
        obs += dose
        exp += 2.0 * af
        if dose == 0:
            n_homref += 1
        n_used += 1

    total = len(variants)
    if n_used == 0 or var <= 0:
        return {"status": "insufficient_data", "confidence": "none",
                "reason": f"None of the {total} scoring-file variants were "
                          "present in the callset at matching coordinates.",
                "coverage": {"total": total, "chip": 0, "imputed": 0,
                             "missing": n_missing, "low_r2": 0,
                             "pct_callable": 0.0}}

    z = (raw - mean) / sqrt(var)
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
    pct_callable = 100.0 * n_used / max(total, 1)
    confidence, note = _pgs_confidence(n_used, total, 0)
    return {
        "status": "computed", "raw_score": round(raw, 4),
        "expected_mean": round(mean, 4), "z_score": round(z, 3),
        "percentile": round(pct, 1), "tier": tier, "tier_class": cls,
        "confidence": confidence, "confidence_note": note,
        "match_basis": "coordinate",
        # Proof the selection effect is gone: on the rsID path this was
        # structurally zero, because a hom-reference call carries no rsID.
        "n_hom_reference": n_homref,
        "dose_ratio": round(obs / exp, 3) if exp > 0 else None,
        "coverage": {"total": total, "chip": n_used, "imputed": 0,
                     "missing": n_missing, "low_r2": 0,
                     "pct_callable": round(pct_callable, 1)},
        "n_variants_used": n_used,
    }
