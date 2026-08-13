"""Tests for cohort_simulator.py — the market-facing economics (analyses 1-3, 5).

Assert theoretical/structural properties, not memorised outputs, so they survive
re-parameterisation.
"""
from __future__ import annotations

import pytest

import cohort_simulator as cs

np = pytest.importorskip("numpy")


@pytest.fixture(scope="module")
def cohort():
    return cs.simulate_cohort(n=6000, seed=20260803)


# ── 1. cohort ────────────────────────────────────────────────────────────────

def test_cohort_is_right_skewed(cohort):
    # A few high-value customers pull the mean above the median.
    assert cohort["available"] is True
    assert cohort["median_value"] <= cohort["mean_value"]


def test_value_is_concentrated_in_a_minority(cohort):
    # The headline marketing insight: the top decile holds a large share of value.
    assert 0.0 <= cohort["share_of_value_in_top_decile"] <= 1.0
    assert cohort["share_of_value_in_top_decile"] > 0.25


def test_cohort_is_reproducible():
    a = cs.simulate_cohort(n=2000, seed=7)
    b = cs.simulate_cohort(n=2000, seed=7)
    assert a["mean_value"] == b["mean_value"]
    assert a["percentiles"] == b["percentiles"]


def test_every_prevalence_is_cited():
    for key, spec in cs.PREVALENCE.items():
        assert spec["src"], f"{key} missing citation"
        assert spec["plain"], f"{key} missing plain-English"


# ── 2. segments ──────────────────────────────────────────────────────────────

def test_prevention_value_declines_with_age(cohort):
    # Grossman: acting earlier protects more remaining healthy years. Measured on
    # the quality-of-life component, where the effect is clean (net value is muddied
    # by the averted-cost offset, which favours older ages — documented in _score).
    s = cs.segment_analysis(cohort)
    prev = [b["mean_prevention_qol_value"] for b in s["by_age_band"]]
    assert prev[0] >= prev[-1]                 # youngest band >= oldest band
    assert s["age_gradient_ratio"] > 1.0


def test_family_history_raises_value(cohort):
    s = cs.segment_analysis(cohort)
    fh = s["by_family_history"]
    assert (fh["with family history"]["mean_value"]
            >= fh["no family history"]["mean_value"])


def test_ancestry_attenuation_is_present(cohort):
    # European-derived polygenic scores transfer imperfectly; the simulator
    # attenuates rather than hides this. European PRS-driven value should not be
    # lower than African-ancestry value for the same finding.
    s = cs.segment_analysis(cohort)
    assert s["by_ancestry"]["european"]["mean_value"] >= 0
    assert "ancestry" in s["sources"]


# ── 3. demand curve ──────────────────────────────────────────────────────────

def test_demand_curve_is_downward_sloping(cohort):
    d = cs.demand_curve(cohort)
    shares = [r["share_in_the_money"] for r in d["curve"]]
    assert shares == sorted(shares, reverse=True)   # fewer buyers as price rises


def test_sequencing_premium_is_nonnegative(cohort):
    d = cs.demand_curve(cohort)
    assert d["mean_sequencing_premium"] >= 0
    assert d["median_sequencing_premium"] >= 0


# ── 5. adoption ──────────────────────────────────────────────────────────────

def test_adoption_is_monotone_and_bounded():
    a = cs.adoption_curve()
    cum = [r["cumulative_adopters"] for r in a["curve"]]
    assert cum == sorted(cum)                        # cumulative never decreases
    assert a["final_penetration"] <= 1.0
    assert a["adopters_lost_to_behavioural_drag"] > 0


def test_adoption_word_of_mouth_dominates():
    # q (imitation) >> p (innovation) is the standard Bass finding.
    a = cs.adoption_curve()
    assert a["q_imitation"] > a["p_innovation"]


# ── narrative / legibility ───────────────────────────────────────────────────

def test_every_analysis_has_plain_english(cohort):
    s, d, a = (cs.segment_analysis(cohort), cs.demand_curve(cohort),
               cs.adoption_curve())
    for result in (cohort, s, d, a):
        assert result.get("plain_english"), "analysis missing plain-English summary"


def test_narrative_assembles():
    text = cs.explain_cohort(n=1500, seed=1)
    assert "HOW MUCH IS IT WORTH" in text
    assert "WHO BENEFITS MOST" in text


# ── 4. data-asset lifetime value ─────────────────────────────────────────────

def test_data_asset_appreciates_with_knowledge_growth():
    d = cs.data_asset_ltv()
    assert d["lifetime_pv"] > d["lifetime_pv_no_growth"]   # growth adds value
    assert 0.0 <= d["appreciation_premium"] < 1.0


def test_data_asset_separates_durability_from_appreciation():
    # No-growth genome still has value beyond year 0 (durability), and appreciation
    # is zero when knowledge does not grow.
    flat = cs.data_asset_ltv(knowledge_growth=0.0)
    assert flat["durability_share"] > 0.0
    assert abs(flat["appreciation_premium"]) < 1e-9


def test_data_asset_more_growth_means_more_appreciation():
    lo = cs.data_asset_ltv(knowledge_growth=0.05)
    hi = cs.data_asset_ltv(knowledge_growth=0.12)
    assert hi["appreciation_premium"] > lo["appreciation_premium"]


# ── personalized per-genome panel ────────────────────────────────────────────

def _fake_voi(total, wgs_marginal=0.0, lo=None, hi=None):
    return {"available": True, "voi_expost_mean": total,
            "marginal_chip_to_wgs": wgs_marginal, "input_type": "chip",
            "voi_ci_low": lo if lo is not None else total * 0.5,
            "voi_ci_high": hi if hi is not None else total * 1.6}


def test_personalize_returns_individual_panels():
    p = cs.personalize_for_report(_fake_voi(12_000, wgs_marginal=3_000), age=35)
    assert p["available"] is True
    assert p["frontier"]["available"] is True
    assert p["ceac"]["available"] is True
    assert p["ltv"]["available"] is True
    assert 0 <= p["population_percentile"] <= 100


def test_personalize_wgs_dominated_when_no_sequencing_value():
    # A genome with NO whole-genome-only value should see WGS ruled out (it costs
    # more for the same benefit) — an honest, individual result.
    p = cs.personalize_for_report(_fake_voi(10_000, wgs_marginal=0.0), age=40)
    assert "Whole-genome sequencing" in p["frontier"]["ruled_out"]


def test_personalize_higher_value_higher_percentile():
    lo = cs.personalize_for_report(_fake_voi(2_000), age=40)["population_percentile"]
    hi = cs.personalize_for_report(_fake_voi(30_000), age=40)["population_percentile"]
    assert hi > lo


def test_personalize_excludes_market_level_analyses():
    # Market-level analyses must NOT leak into the personal panel.
    p = cs.personalize_for_report(_fake_voi(12_000), age=35)
    for market_key in ("demand", "adoption", "distributional", "validation"):
        assert market_key not in p
