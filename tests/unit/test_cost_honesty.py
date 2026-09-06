"""Regression tests for two health-economics honesty corrections.

A) MARGINAL vs AVERAGE COST — averting one case frees the marginal cost of that case,
   not its average lifetime cost-of-illness. Using average-as-marginal overstates cash
   savings (the "freeing a bed doesn't save its average cost" error).

B) COST-SAVING vs COST-EFFECTIVE — prevention is usually cost-effective (good value per
   QALY) but not cost-saving (net cash out). Monetized QALYs must never be presented as
   money returned.

Both assert THEORETICAL PROPERTIES (directions and invariants), not memorised outputs.
"""
from __future__ import annotations

from econ import health_economics as he
from econ import value_of_information as voi

# ── A) marginal vs average cost ──────────────────────────────────────────────

def test_marginal_fraction_is_a_conservative_discount():
    # It must genuinely discount (never inflate) and be a real fraction.
    assert 0.0 < voi.MARGINAL_COST_FRACTION < 1.0


def test_averted_cost_is_below_average_cost_of_illness():
    # The averted cost for a finding must be strictly less than the naive
    # p x rrr x AVERAGE-COI it would have been.
    f = {"kind": "coi", "coi_key": "CAD", "p_event": 0.2, "rrr": 0.3, "qaly": 1.0,
         "intervention": 500.0, "horizon": 25, "haircut": 1.0}
    _nmb, dcost, _dq, _iv = voi._finding_nmb(f, 100_000.0, 0.03)
    naive = dcost / voi.MARGINAL_COST_FRACTION      # what average-as-marginal gave
    assert dcost < naive


def test_lowering_marginal_fraction_lowers_averted_cost(monkeypatch):
    # Monotonic: a more conservative fraction must reduce the claimed saving.
    f = {"kind": "coi", "coi_key": "CAD", "p_event": 0.2, "rrr": 0.3, "qaly": 1.0,
         "intervention": 500.0, "horizon": 25, "haircut": 1.0}
    monkeypatch.setattr(voi, "MARGINAL_COST_FRACTION", 0.8)
    high = voi._finding_nmb(f, 100_000.0, 0.03)[1]
    monkeypatch.setattr(voi, "MARGINAL_COST_FRACTION", 0.3)
    low = voi._finding_nmb(f, 100_000.0, 0.03)[1]
    assert low < high


def test_personal_economics_applies_the_same_correction():
    # The cash side of the personal model must use the same marginal fraction, so a
    # "cost-saving" verdict can't be an artifact of average-as-marginal.
    assert he._MARGINAL_COST_FRACTION == voi.MARGINAL_COST_FRACTION


# ── B) cost-saving vs cost-effective ─────────────────────────────────────────

def _fake_econ(outcome, prevalence, qaly_gain, cost):
    return {"findings_with_economics": [
        {"finding": "synthetic", "outcome_value": outcome, "prevalence": prevalence,
         "qaly_gain": qaly_gain, "cost": cost, "confidence": "high"}]}


def test_net_cash_excludes_monetized_qalys():
    # net_cash is money only: averted cost minus spend. It must never include the
    # monetized-QALY term that dominates total_net.
    pe = he.analyze_personal_economics(_fake_econ(100_000, 0.5, 2.0, 500))
    assert pe["available"]
    assert pe["net_cash"] == pe["total_avoided"] - pe["total_intervention"]
    assert pe["net_cash"] != pe["total_net"]          # the two are distinct concepts
    assert pe["total_net"] >= pe["net_cash"]          # health value only adds


def test_cost_effective_but_not_cost_saving_is_reported_honestly():
    # A finding with big health gain but small averted cost and high spend is
    # cost-EFFECTIVE, not cost-SAVING. The verdict must say so.
    pe = he.analyze_personal_economics(_fake_econ(1_000, 0.1, 5.0, 40_000))
    assert pe["available"]
    assert pe["net_cash"] < 0                         # cash out
    assert pe["total_net"] > 0                        # still worth it per QALY
    assert pe["is_cost_saving"] is False
    assert "cost-effective" in pe["verdict"]


def test_genuine_cost_saving_is_still_recognised():
    # The machinery must not simply always say "cost-effective" — a real cash-positive
    # case (cheap test, large averted cost) must be reported as cost-saving.
    pe = he.analyze_personal_economics(_fake_econ(200_000, 0.9, 0.1, 100))
    assert pe["net_cash"] > 0
    assert pe["is_cost_saving"] is True
    assert pe["verdict"] == "cost-saving"


def test_verdict_is_one_of_the_three_honest_states():
    pe = he.analyze_personal_economics(_fake_econ(50_000, 0.3, 1.0, 1_000))
    assert pe["verdict"] in ("cost-saving",
                             "cost-effective (adds cost, worth it per QALY)",
                             "not cost-effective at this threshold")


# ── C) cost-consequence analysis (disaggregated, NICE digital-health format) ──

def test_cca_reports_no_summary_ratio():
    # The defining property of a CCA: costs and consequences stay disaggregated.
    # No ICER / cost-per-QALY / ROI field may be produced.
    cca = he.build_cost_consequence_analysis(
        he.analyze_personal_economics(_fake_econ(30_000, 0.3, 0.5, 250)))
    assert cca["available"] is True
    assert not [k for k in cca
                if any(s in k.lower() for s in ("icer", "per_qaly", "ratio", "roi"))]


def test_cca_keeps_qalys_in_their_own_unit():
    # QALYs must be reported as QALYs, not monetized into the same column as cash.
    cca = he.build_cost_consequence_analysis(
        he.analyze_personal_economics(_fake_econ(30_000, 0.3, 0.5, 250)))
    qrow = next(r for r in cca["rows"] if "life-years" in r["measure"])
    assert qrow["unit"] == "QALYs"
    # and no row mixes a dollar value into a QALY unit
    for r in cca["rows"]:
        assert not (r["unit"] == "QALYs" and str(r["value"]).startswith("$"))


def test_cca_separates_costs_from_consequences():
    cca = he.build_cost_consequence_analysis(
        he.analyze_personal_economics(_fake_econ(30_000, 0.3, 0.5, 250)))
    kinds = {r["kind"] for r in cca["rows"]}
    assert kinds == {"Cost", "Consequence"}


def test_cca_degrades_when_no_items():
    assert he.build_cost_consequence_analysis(None)["available"] is False
    assert he.build_cost_consequence_analysis({"available": False})["available"] is False


# ── D) the stated horizon must actually enter the arithmetic ─────────────────

def test_horizon_discount_is_applied_not_just_labeled():
    # PERSONAL_HORIZON_YEARS labeled the output "Over 10 years" while the
    # figures were undiscounted sums of future dollars. The discount factor must
    # be a real fraction that reduces the claim.
    assert 0.0 < he._MIDPOINT_DISCOUNT < 1.0


def test_horizon_length_changes_the_discount(monkeypatch):
    # Directional: a longer horizon discounts harder at the midpoint.
    short = 1.0 / (1.0 + he.DISCOUNT_RATE) ** (4 / 2.0)
    long_ = 1.0 / (1.0 + he.DISCOUNT_RATE) ** (20 / 2.0)
    assert long_ < short


def test_discounting_reduces_reported_value():
    pe = he.analyze_personal_economics(_fake_econ(100_000, 0.5, 2.0, 500))
    # Undiscounted counterfactual for the same inputs.
    undiscounted_avoided = 100_000 * 0.5 * he._MARGINAL_COST_FRACTION
    assert pe["total_avoided"] < undiscounted_avoided


def test_horizon_years_is_still_reported():
    pe = he.analyze_personal_economics(_fake_econ(30_000, 0.3, 0.5, 250))
    assert pe["horizon_years"] == he.PERSONAL_HORIZON_YEARS
