"""Tests for heor_frontier.py — HEOR deliverables (analyses 6-9).

Theoretical-property assertions, plus a check that every output is legible
(carries a plain-English sentence and a citation).
"""
from __future__ import annotations

import pytest

from econ import frontier as hf

np = pytest.importorskip("numpy")


# ── 6. value-based price ─────────────────────────────────────────────────────

def test_vbp_rises_with_willingness_to_pay():
    v = hf.value_based_price(incremental_qaly=0.12, incremental_cost_offset=800)
    r = v["value_based_price_range"]
    assert r["at_50k"] < r["at_100k"] < r["at_150k"]


def test_vbp_includes_cost_offset():
    v = hf.value_based_price(incremental_qaly=0.10, incremental_cost_offset=1_000,
                             wtp=100_000)
    # 0.10 QALY * 100k + 1000 offset = 11,000
    assert v["value_based_price"] == 11_000


# ── 7. cost-effectiveness frontier ───────────────────────────────────────────

def test_frontier_icers_are_increasing():
    f = hf.cost_effectiveness_frontier()
    icers = [o["icer_vs_previous"] for o in f["frontier"]
             if o["icer_vs_previous"] is not None]
    assert icers == sorted(icers)             # defining property of the frontier


def test_frontier_rules_out_inefficient_strategies():
    # A strictly dominated option (costs more, does less) must be removed.
    strategies = [
        {"name": "None", "cost": 0.0, "qaly": 20.0},
        {"name": "Good", "cost": 300.0, "qaly": 20.10},
        {"name": "Bad", "cost": 500.0, "qaly": 20.05},   # dominated by Good
    ]
    f = hf.cost_effectiveness_frontier(strategies)
    assert "Bad" in f["ruled_out"]
    assert "Bad" not in [o["name"] for o in f["frontier"]]


def test_frontier_recommends_within_threshold():
    f = hf.cost_effectiveness_frontier(wtp=100_000)
    assert f["recommended_strategy"] in [o["name"] for o in f["frontier"]]


# ── 7B. probabilistic frontier / CEAC ────────────────────────────────────────

def test_ceac_probabilities_sum_to_one():
    p = hf.frontier_psa(n_mc=2000)
    for row in p["ceac"]:
        assert abs(sum(row["p_optimal"].values()) - 1.0) < 0.02


def test_ceac_no_testing_wins_at_zero_wtp():
    # At WTP = 0, only cost matters, so the cheapest (no testing) must win.
    p = hf.frontier_psa(n_mc=2000)
    at_zero = next(r for r in p["ceac"] if r["wtp"] == 0)
    assert at_zero["most_likely_optimal"] == "No testing"


def test_ceac_higher_wtp_favors_more_effective():
    # As willingness-to-pay rises, the most effective strategy should gain share.
    p = hf.frontier_psa(n_mc=3000)
    wgs = "Whole-genome sequencing"
    lo = next(r for r in p["ceac"] if r["wtp"] == 0)["p_optimal"][wgs]
    hi = next(r for r in p["ceac"] if r["wtp"] == 200_000)["p_optimal"][wgs]
    assert hi > lo


# ── 8. distributional CEA ────────────────────────────────────────────────────

def test_dcea_upweights_worse_off_groups():
    dc = hf.distributional_cea(inequality_aversion=2.0)
    # A group with below-average baseline health should get weight > 1.
    below = [g for g in dc["groups"] if g["equity_weight"] > 1.0]
    assert below, "expected at least one up-weighted worse-off group"


def test_dcea_surfaces_the_benefit_gap():
    dc = hf.distributional_cea()
    assert dc["benefit_gap_ratio"] > 1.0
    # The honest point: equity weighting cannot close a large portability gap.
    assert dc["equity_impact_ratio"] < dc["benefit_gap_ratio"]


def test_dcea_zero_benefit_group_does_not_crash():
    # REGRESSION: a group with ~zero transferable PRS benefit (an infinite gap — the
    # exact scenario this function exists to display) set benefit_gap_ratio=None,
    # which then crashed the plain_english f-string.
    dc = hf.distributional_cea([
        {"name": "A", "population_share": 0.5, "baseline_health": 68, "qaly_gain": 0.10},
        {"name": "B", "population_share": 0.5, "baseline_health": 64, "qaly_gain": 0.0}])
    assert dc["available"] is True
    assert dc["benefit_gap_ratio"] is None            # infinite gap, reported as None
    assert "essentially no benefit" in dc["plain_english"]


def test_dcea_negative_benefit_group_does_not_crash():
    dc = hf.distributional_cea([
        {"name": "A", "population_share": 0.5, "baseline_health": 68, "qaly_gain": 0.10},
        {"name": "B", "population_share": 0.5, "baseline_health": 64, "qaly_gain": -0.02}])
    assert dc["available"] is True


def test_dcea_more_aversion_weights_harder():
    low = hf.distributional_cea(inequality_aversion=1.0)
    high = hf.distributional_cea(inequality_aversion=4.0)
    worst_low = min(g["equity_weight"] for g in low["groups"])
    worst_high = min(g["equity_weight"] for g in high["groups"])
    # Higher aversion spreads the weights further from 1.
    assert (max(g["equity_weight"] for g in high["groups"]) - worst_high) >= \
           (max(g["equity_weight"] for g in low["groups"]) - worst_low)


# ── 9. validation ────────────────────────────────────────────────────────────

def test_validation_actually_computes_icers():
    # REGRESSION: our_icer must be COMPUTED by the engine, not hardcoded. Recomputing
    # from published inputs must reproduce the same value, and it must not equal a
    # round hand-picked number.
    val = hf.validate_against_published()
    cyp = next(c for c in val["cases"] if "CYP2C19" in c["scenario"])
    icer, _c, _q = hf._pgx_icer(test_cost=100.0, p_event=0.13, rrr=0.25,
                                event_cost=30_000.0, qaly_per_event=3.0,
                                benefiting_fraction=0.30, added_treatment_cost=900.0)
    assert cyp["our_icer"] == round(icer)          # computed, reproducible


def test_validation_agrees_on_direction_and_order_of_magnitude():
    val = hf.validate_against_published()
    # Direction correct in most cases; the PGx scenarios land within an order of
    # magnitude. (Lynch is a documented structural miss — see below.)
    assert val["n_direction_correct"] >= 2
    cyp = next(c for c in val["cases"] if "CYP2C19" in c["scenario"])
    assert cyp["within_order_of_magnitude"] is True
    dpyd = next(c for c in val["cases"] if "DPYD" in c["scenario"])
    assert dpyd["our_icer"] < 0                     # reproduced as cost-saving


def test_validation_reports_its_misses_honestly():
    # A validation where everything passes looks tuned. The Lynch structural
    # mismatch must be surfaced with a reason, not hidden.
    val = hf.validate_against_published()
    lynch = next(c for c in val["cases"] if "Lynch" in c["scenario"])
    if not lynch["within_order_of_magnitude"]:
        assert lynch.get("miss_reason"), "a miss must carry an explanation"


def test_validation_wrong_direction_not_counted_as_agreement():
    # REGRESSION: a published cost-saving result the engine calls COSTLY was labeled
    # "both dominant" and miscounted as direction-correct.
    val = hf.validate_against_published(
        [{"name": "X", "published_icer": -5000, "our_icer": 40000,
          "tolerance_pct": 0.5, "source": "..."}])
    row = val["cases"][0]
    assert row["ratio_to_published"] != "both dominant"
    assert row["direction_correct"] is False
    assert val["n_direction_correct"] == 0


def test_frontier_empty_list_no_crash():
    # REGRESSION: cost_effectiveness_frontier([]) indexed rows[0] and crashed.
    f = hf.cost_effectiveness_frontier([])
    assert f["available"] is False


# ── legibility ───────────────────────────────────────────────────────────────

def test_every_analysis_is_cited_and_legible():
    for result in (hf.value_based_price(0.1, 500),
                   hf.cost_effectiveness_frontier(),
                   hf.distributional_cea(),
                   hf.validate_against_published()):
        assert result.get("plain_english"), "missing plain-English summary"
        assert result.get("src"), "missing citation"
