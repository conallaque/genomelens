"""Unit tests for markov_model.py — HEOR cohort model + budget impact analysis.

Asserts the structural conventions an HTA reviewer checks, not memorised outputs.
"""
from __future__ import annotations

from econ import markov as mk

# ── rate → probability ───────────────────────────────────────────────────────

def test_rate_to_probability_never_exceeds_one():
    # The naive r*dt form breaks for large rates; p = 1-exp(-r*dt) cannot exceed 1.
    assert mk.rate_to_prob(1.5) < 1.0
    assert mk.rate_to_prob(10.0) < 1.0
    assert abs(mk.rate_to_prob(0.01) - 0.00995) < 1e-4      # ≈ r for small r
    assert mk.rate_to_prob(0.0) == 0.0


# ── transition matrix ────────────────────────────────────────────────────────

def test_transition_matrix_rows_sum_to_one():
    P = mk.build_transition_matrix(0.02, 0.03, 0.01)
    for row in P:
        assert abs(sum(row) - 1.0) < 1e-12


def test_death_state_is_absorbing():
    P = mk.build_transition_matrix(0.02, 0.03, 0.01)
    assert P[2] == [0.0, 0.0, 1.0]


def test_transition_matrix_handles_competing_exits_over_one():
    # Two large exit probabilities must be normalised, not left summing above 1.
    P = mk.build_transition_matrix(0.8, 0.1, 0.7)
    for row in P:
        assert abs(sum(row) - 1.0) < 1e-12
        assert all(0.0 <= x <= 1.0 for x in row)


# ── cohort model ─────────────────────────────────────────────────────────────

def test_cohort_is_conserved_and_structurally_valid():
    r = mk.markov_cost_effectiveness()
    v = mk.validate_markov(r)
    assert v["all_passed"] is True, v["checks"]


def test_intervention_reduces_disease_and_gains_qalys():
    r = mk.markov_cost_effectiveness()
    assert r["incremental_qaly"] > 0                     # prevention gains QALYs
    assert (r["genomic_guided"]["final_distribution"]["disease"]
            < r["standard_care"]["final_distribution"]["disease"])


def test_discounting_reduces_totals():
    undisc = mk.run_markov("standard_care", discount_rate=0.0)
    disc = mk.run_markov("standard_care", discount_rate=0.03)
    assert undisc["total_qaly"] > disc["total_qaly"]
    assert undisc["total_cost"] > disc["total_cost"]


def test_half_cycle_correction_changes_result():
    with_hcc = mk.run_markov("standard_care", half_cycle=True)["total_qaly"]
    without = mk.run_markov("standard_care", half_cycle=False)["total_qaly"]
    assert abs(with_hcc - without) > 1e-6


def test_icer_rises_with_intervention_cost():
    # The ratio is withheld in the dominance quadrants (a negative ICER reads
    # as a bargain whether the strategy is excellent or terrible), so compare
    # only the scenarios where a ratio is defined, and require the dominant
    # ones to come first — cheap interventions dominate, dear ones do not.
    results = [mk.markov_cost_effectiveness(cost_intervention_annual=c)
               for c in (250, 500, 1000, 2000)]
    defined = [r["icer"] for r in results if r["icer"] is not None]
    assert defined == sorted(defined), (
        f"costlier intervention should worsen the ICER: {defined}")
    dominant_flags = [r["dominant"] for r in results]
    assert dominant_flags == sorted(dominant_flags, reverse=True), (
        "dominance should be lost as the intervention gets more expensive, "
        "not regained")


def test_no_negative_icer_is_ever_reported():
    # A negative ratio is ambiguous by construction and HTA convention is to
    # state dominance instead. The verdict already said "dominant"; the ratio
    # was printed next to it anyway.
    for c in (0, 100, 250, 1000, 100_000):
        r = mk.markov_cost_effectiveness(cost_intervention_annual=c)
        assert r["icer"] is None or r["icer"] >= 0, (
            f"reported ICER {r['icer']} at intervention cost {c}")
        if r["icer"] is None:
            assert r["dominant"] or r["dominated"] or \
                abs(r["incremental_qaly"]) < 1e-9, (
                "ICER withheld without a dominance reason")


def test_verdict_flips_when_intervention_is_expensive_enough():
    cheap = mk.markov_cost_effectiveness(cost_intervention_annual=100)
    dear = mk.markov_cost_effectiveness(cost_intervention_annual=100_000)
    assert "dominant" in cheap["verdict"]
    assert "not cost-effective" in dear["verdict"]   # the model can say no


# ── budget impact ────────────────────────────────────────────────────────────

def test_budget_impact_scales_with_plan_size():
    small = mk.budget_impact(plan_members=100_000)
    big = mk.budget_impact(plan_members=1_000_000)
    # PMPM is per-member, so it should be ~invariant to plan size...
    assert abs(small["peak_pmpm"] - big["peak_pmpm"]) < 1e-6
    # ...while absolute dollars scale with membership.
    assert big["peak_net_budget_impact"] > small["peak_net_budget_impact"]


def test_budget_impact_front_loads_cost_then_offsets():
    b = mk.budget_impact()
    first, last = b["rows"][0], b["rows"][-1]
    assert first["offsets"] == 0                  # no offsets in year 1
    assert last["offsets"] > last["cost_testing"]  # offsets dominate by year 5
    assert last["net_budget_impact"] < first["net_budget_impact"]


def test_budget_impact_reports_pmpm():
    b = mk.budget_impact()
    for row in b["rows"]:
        expected = row["net_budget_impact"] / (b["plan_members"] * 12.0)
        assert abs(row["pmpm"] - expected) < 1e-3
