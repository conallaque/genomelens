"""
Unit tests for `supplements.py`.

The critical correctness invariants here:

  1. `_risk_dose` returns the right count regardless of which strand a chip
     reports on. The previous `count("T") + count("A")` antipattern would
     double-count and falsely flag a hypothetical "AT" genotype as homozygous.
  2. Chip gaps are surfaced explicitly — a missing rsID produces a `_chip_gap`
     placeholder, not a silent skip.
  3. Bloodwork-driven refinements bump matching supplements up a tier when
     measured labs confirm the genetic prediction.
"""

from __future__ import annotations

import pandas as pd
import pytest

from wellness import supplements as sup

# ── _risk_dose strand-awareness ──────────────────────────────────────────────

@pytest.mark.parametrize("genotype,expected", [
    # MTHFR C677T:  + strand uses C/T, − strand uses G/A
    ("TT", 2),    # homozygous risk on + strand
    ("AA", 2),    # homozygous risk on − strand
    ("CT", 1),    # heterozygous on + strand
    ("GA", 1),    # heterozygous on − strand
    ("CC", 0),    # homozygous ancestral on + strand
    ("GG", 0),    # homozygous ancestral on − strand
])
def test_risk_dose_strand_invariant(genotype: str, expected: int) -> None:
    """Dose is correct regardless of which strand the chip reports."""
    df = pd.DataFrame({"genotype": [genotype]}, index=["rs1801133"])
    assert sup._risk_dose(df, "rs1801133", risk_allele="T", ref_allele="C") == expected


def test_risk_dose_missing_rsid_returns_none() -> None:
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs999999"])
    assert sup._risk_dose(df, "rs1801133", risk_allele="T", ref_allele="C") is None


def test_risk_dose_no_call_returns_none() -> None:
    df = pd.DataFrame({"genotype": ["--"]}, index=["rs1801133"])
    assert sup._risk_dose(df, "rs1801133", risk_allele="T", ref_allele="C") is None


def test_risk_dose_ambiguous_falls_back_to_plus_strand() -> None:
    """If the genotype contains alleles from neither strand cleanly (rare noise),
    fall back to the + strand interpretation rather than crash."""
    df = pd.DataFrame({"genotype": ["TG"]}, index=["rs1801133"])
    # G is the − strand complement of C; T is + strand risk.
    # Ambiguous → counts T on + strand → 1
    assert sup._risk_dose(df, "rs1801133", risk_allele="T", ref_allele="C") == 1


# ── Strand-aware rules produce equivalent output across strands ──────────────

def test_mthfr_rule_strand_equivalence() -> None:
    """A user genotyped CC vs GG (same biological state, different strand) must
    produce the SAME supplement output."""
    df_plus = pd.DataFrame({"genotype": ["CC"]}, index=["rs1801133"])
    df_minus = pd.DataFrame({"genotype": ["GG"]}, index=["rs1801133"])
    assert sup._rule_mthfr(df_plus) == sup._rule_mthfr(df_minus)
    # Both should produce nothing (ancestral homozygous)
    assert sup._rule_mthfr(df_plus) == []


def test_mthfr_rule_homozygous_risk_is_essential() -> None:
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs1801133"])
    recs = sup._rule_mthfr(df)
    assert len(recs) == 1
    assert recs[0]["tier"] == sup.TIER_ESSENTIAL
    assert "800-1000 mcg" in recs[0]["dose"]


def test_mthfr_rule_heterozygous_is_recommended() -> None:
    df = pd.DataFrame({"genotype": ["CT"]}, index=["rs1801133"])
    recs = sup._rule_mthfr(df)
    assert len(recs) == 1
    assert recs[0]["tier"] == sup.TIER_RECOMMENDED


# ── Chip-gap surfacing ───────────────────────────────────────────────────────

def test_chip_gaps_surface_missing_rsids() -> None:
    """A near-empty chip should produce chip_gap entries for every key rule SNP."""
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs1801133"])
    result = sup.build_supplement_stack(snps_df=df)
    assert result["n_chip_gaps"] >= 14   # we have 16 tracked; this chip has only 1
    # And the gaps must contain the explicit "not typed" marker
    assert any("not typed" in g["snps"][0] for g in result["chip_gaps"])


def test_chip_gap_card_renders_with_grey_style() -> None:
    """Gaps surface in the HTML rather than being silently dropped."""
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs1801133"])
    result = sup.build_supplement_stack(snps_df=df)
    html = sup.render_supplements_html(result)
    assert "Chip gaps" in html
    assert "sp-chip_gap" in html
    assert "[Not tested]" in html


def test_no_chip_gaps_when_all_rules_typed(synthetic_snps_df) -> None:
    """The synthetic chip fixture covers every key rule SNP."""
    result = sup.build_supplement_stack(snps_df=synthetic_snps_df)
    assert result["n_chip_gaps"] == 0
    assert result["chip_gaps"] == []


# ── HFE / iron interaction ──────────────────────────────────────────────────

def test_hfe_carrier_disables_iron_recommendation() -> None:
    """A user with HFE H63D carrier status should see 'AVOID iron' in the output,
    not a gentle-iron recommendation, even if TMPRSS6 suggests low iron."""
    df = pd.DataFrame(
        {"genotype": ["AA"]},
        index=["rs855791"],  # TMPRSS6 low-iron allele
    )
    carrier = {
        "affected": [],
        "carriers": [{
            "gene": "HFE", "rsid": "rs1799945", "variant": "H63D",
            "pathogenic_allele": "G", "dosage": 1,
        }],
    }
    result = sup.build_supplement_stack(snps_df=df, carrier_result=carrier)
    iron_recs = [r for r in result["avoid"] if "iron" in r["name"].lower()]
    assert len(iron_recs) == 1
    assert iron_recs[0]["tier"] == "avoid"


# ── Bloodwork-driven escalation ─────────────────────────────────────────────

def test_high_crp_bumps_curcumin_to_essential(synthetic_snps_df) -> None:
    """If measured CRP is above-average, curcumin recommendation jumps to
    Essential tier even if the genotype-only call placed it in Recommended."""
    bw = {
        "status": "ok",
        "rows": [{
            "trait": "C-Reactive Protein", "predicted": 1.8, "actual": 4.6,
            "unit": "mg/L", "delta_sd": 1.87, "actual_tier": "Above average",
            "verdict": "Diverged", "interpretation": "test",
        }],
    }
    result = sup.build_supplement_stack(
        snps_df=synthetic_snps_df, bloodwork_result=bw,
    )
    curcumin = [r for tier in result["tiers"].values()
                for r in (tier if isinstance(tier, list) else [])
                if r["name"].startswith("Curcumin")]
    assert curcumin, "expected a Curcumin recommendation"
    assert curcumin[0]["tier"] == sup.TIER_ESSENTIAL
    assert "confirms" in curcumin[0]["reasoning"].lower()


# ── Aggregator behavior ─────────────────────────────────────────────────────

def test_status_no_data_when_no_input() -> None:
    result = sup.build_supplement_stack(snps_df=None)
    assert result["status"] == "no_data"
    assert result["n_supplements"] == 0


def test_monthly_cost_is_non_negative(synthetic_snps_df) -> None:
    result = sup.build_supplement_stack(snps_df=synthetic_snps_df)
    assert result["total_estimated_cost_usd_monthly"] >= 0
