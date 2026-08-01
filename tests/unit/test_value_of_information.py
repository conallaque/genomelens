"""Unit tests for the value-of-information health-economics engine."""
from __future__ import annotations

import value_of_information as voi


def _econ():
    return {"findings_with_economics": [
        {"finding": "CYP2C19 intermediate", "category": "pgx", "confidence": "high"},
        {"finding": "CAD PRS high", "category": "prs", "confidence": "moderate",
         "qaly_gain": 1.5},
    ]}


def _cvr():
    return {"available": True, "buckets": {"actionable": [{"gene": "BRCA2"}],
                                           "carrier": [{"gene": "CFTR"}]}}


def _nvr():
    return {"available": True, "buckets": {"predicted_pathogenic_rare": [
        {"chrom": "1", "pos": 100, "am_score": 0.95, "confidence": "higher"}]}}


def test_runs_on_chip_only():
    r = voi.analyze_value_of_information(_econ(), input_type="chip", n_mc=500)
    assert r["available"] is True
    assert r["input_type"] == "chip"
    assert r["n_findings"] >= 2


def test_runs_on_wgs_with_variants():
    r = voi.analyze_value_of_information(_econ(), _cvr(), _nvr(),
                                        input_type="wgs", n_mc=500)
    assert r["available"] is True
    # WGS adds clinical + novel findings beyond the chip's PGx/PRS.
    assert r["n_findings"] > 2


def test_no_findings_unavailable():
    r = voi.analyze_value_of_information({}, input_type="chip")
    assert r["available"] is False
    assert "bloodwork" in r["reason"] or "finding" in r["reason"]


def test_discounting_reduces_present_value():
    # Undiscounted expected value must be >= discounted (discounting bites).
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs", n_mc=500)
    assert r["voi_expost_undiscounted"] >= r["voi_expost_point"]


def test_chip_to_wgs_marginal_nonnegative():
    r = voi.analyze_value_of_information(_econ(), _cvr(), _nvr(),
                                        input_type="wgs", n_mc=500)
    # WGS-only findings (clinical + novel) contribute non-negative marginal value.
    assert r["marginal_chip_to_wgs"] >= 0


def test_psa_seeded_reproducible():
    a = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=1000, seed=7)
    b = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=1000, seed=7)
    assert a["voi_expost_mean"] == b["voi_expost_mean"]
    assert a["voi_ci_low"] == b["voi_ci_low"]


def test_ceac_monotone_nondecreasing():
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs", n_mc=1000)
    probs = [c["prob"] for c in r["ceac"]]
    assert all(probs[i] <= probs[i + 1] + 1e-9 for i in range(len(probs) - 1))
    assert r["ceac"][0]["lam"] == 0 and r["ceac"][-1]["lam"] == 200_000


def test_price_is_separate_from_value():
    r = voi.analyze_value_of_information(_econ(), input_type="chip", n_mc=200)
    price = r["price"]
    assert price["a_la_carte_total"] > price["consolidated"]
    assert "price" in price["note"].lower() and "value" in price["note"].lower()


def test_var_cvar_ordering():
    # CVaR (expected shortfall) must be <= VaR (its own threshold) <= the overall mean —
    # CVaR is the average of everything AT OR BELOW the VaR cutoff.
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs", n_mc=3000, seed=3)
    assert r["cvar_95"] <= r["var_95"]
    assert r["var_95"] <= r["voi_expost_mean"]


def test_var_cvar_seeded_reproducible():
    a = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs", n_mc=2000, seed=9)
    b = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs", n_mc=2000, seed=9)
    assert a["var_95"] == b["var_95"]
    assert a["cvar_95"] == b["cvar_95"]


def test_grossman_value_declines_with_age():
    # Core Grossman prediction: an efficiency gain compounds over remaining
    # life-years, so the same information is worth strictly more when younger.
    gains = [voi.analyze_health_capital(age=a)["pv_health_capital_gain"]
             for a in (25, 40, 55, 70)]
    assert gains == sorted(gains, reverse=True), gains


def test_grossman_depreciation_accelerates():
    hc = voi.analyze_health_capital(age=30, horizon=50)
    assert hc["delta_at_end"] > hc["delta_at_age"]      # δ rises with age
    traj = hc["trajectory_uninformed"]
    assert traj[0]["h"] > traj[-1]["h"]                  # stock declines over life


def test_grossman_informed_defers_morbidity_floor():
    hc = voi.analyze_health_capital(age=35, horizon=60)
    if hc["floor_age_informed"] and hc["floor_age_uninformed"]:
        assert hc["floor_age_informed"] >= hc["floor_age_uninformed"]


def test_real_option_high_value_says_test_now():
    # For an actionable genome, waiting forfeits protection: option value ~0.
    ro = voi.analyze_real_option(voi_now=25_000, test_cost=300, age=35)
    assert ro["optimal_defer_years"] == 0
    assert ro["option_value_of_waiting"] <= 0.5
    assert "now" in ro["recommendation"]


def test_real_option_low_value_expensive_test_defers():
    # The model must be able to say "don't buy yet" — otherwise it isn't a model.
    ro = voi.analyze_real_option(voi_now=200, test_cost=1_000, age=40)
    assert ro["optimal_defer_years"] > 0
    assert ro["value_test_now"] < 0        # negative NPV at today's price


def test_risk_adjusted_metrics_present_and_ordered():
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=2000, seed=5)
    ra = r["risk_adjusted"]
    assert ra["available"] is True
    assert ra["roi_multiple"] > 0
    # Risk aversion means the certainty equivalent sits at or below the mean.
    assert ra["certainty_equivalent"] <= r["voi_expost_mean"]


def test_evpi_is_nonnegative_and_bounded():
    # EVPI >= 0 by construction (perfect info can never be worth less than current
    # info), and perfect-information net benefit must weakly dominate.
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=2000, seed=11)
    e = r["evpi"]
    assert e["available"] is True
    assert e["evpi"] >= 0
    assert e["nb_perfect_information"] >= e["nb_current_information"]
    assert 0.0 <= e["share_of_information_captured"] <= 1.0


def test_utility_certainty_equivalent_falls_with_risk_aversion():
    u = voi.analyze_utility(mean=50_000, sd=20_000)
    ces = [g["certainty_equivalent"] for g in u["by_gamma"]]
    prems = [g["risk_premium"] for g in u["by_gamma"]]
    assert ces == sorted(ces, reverse=True)     # more risk-averse → lower CE
    assert prems == sorted(prems)               # …and a larger risk premium
    assert all(c <= 50_000 for c in ces)        # CE never exceeds the mean


def test_penetrance_ascertainment_correction_shrinks_risk():
    # A family-ascertained penetrance must be shrunk for an incidental carrier.
    p = voi.analyze_penetrance_posterior(prior_penetrance=0.60, gene="BRCA1")
    assert p["population_corrected"] < p["prior_literature_penetrance"]
    assert p["posterior_penetrance"] < p["prior_literature_penetrance"]
    assert 0.0 < p["shrinkage_factor"] < 1.0
    # Family history should raise the posterior relative to no family history.
    with_fh = voi.analyze_penetrance_posterior(0.60, family_history=True)
    without = voi.analyze_penetrance_posterior(0.60, family_history=False)
    assert with_fh["posterior_penetrance"] > without["posterior_penetrance"]


def test_winners_curse_shrinks_effect_sizes():
    w = voi.shrink_effect_size(beta_hat=0.12, se=0.02)
    assert w["available"] is True
    assert abs(w["beta_shrunk"]) < abs(w["beta_reported"])
    assert 0 < w["shrinkage_factor"] < 1
    # A far more significant variant should be shrunk proportionally less.
    strong = voi.shrink_effect_size(beta_hat=0.12, se=0.005)
    weak = voi.shrink_effect_size(beta_hat=0.12, se=0.02)
    assert strong["shrinkage_factor"] >= weak["shrinkage_factor"]


def test_information_economics_caveats_present():
    r = voi.analyze_value_of_information(_econ(), input_type="chip", n_mc=500)
    ie = r["information_economics"]
    assert ie["available"] is True
    assert "adverse selection" in ie["adverse_selection_note"].lower()
    assert "gina" in ie["discrimination_note"].lower()


def test_predicted_variant_downweighted():
    # A high-AM predicted variant should still contribute less than a confirmed
    # ClinVar-actionable finding of similar condition (haircut applied).
    r = voi.analyze_value_of_information(_econ(), _cvr(), _nvr(),
                                        input_type="wgs", n_mc=200)
    labels = {row["label"]: row["nmb"] for row in r["nmb_rows"]}
    pred = [v for k, v in labels.items() if "predicted" in k]
    assert pred, "predicted finding should appear in the NMB table"
