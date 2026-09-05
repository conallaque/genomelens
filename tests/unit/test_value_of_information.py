"""Unit tests for the value-of-information health-economics engine."""
from __future__ import annotations

from econ import value_of_information as voi


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


def test_real_health_economics_category_strings_are_recognised():
    # REGRESSION: health_economics.py labels findings with human-readable sources
    # ("Pharmacogenomics", "Polygenic Risk", "Genotype", ...). An earlier version
    # matched only short keys ("pgx"/"prs"/"apoe"), so every chip finding was
    # silently dropped and the engine reported "no findings" on chip input.
    real = {"Pharmacogenomics": "pgx", "Polygenic Risk": "coi",
            "Genotype": "coi", "Exercise / Lifestyle": "coi"}
    for cat, expected_kind in real.items():
        kind, _ = voi._classify_category(cat, "APOE e3/e4" if cat == "Genotype" else "")
        assert kind == expected_kind, f"{cat!r} → {kind!r}, expected {expected_kind!r}"
    # "Longevity" is deliberately NOT valued: the composite re-aggregates
    # variants already valued individually, so routing it to a cost-of-illness
    # anchor counted the same genotypes twice. It must resolve to no kind AND
    # carry a documented reason, so it reads as a decision, not an oversight.
    assert voi._classify_category("Longevity")[0] == ""
    assert voi._not_valued_reason("Longevity")
    # Short internal keys must keep working too (both conventions supported).
    assert voi._classify_category("pgx")[0] == "pgx"
    assert voi._classify_category("apoe")[1] == "Alzheimer"
    # Unknown categories are ignored rather than mis-bucketed.
    assert voi._classify_category("Something Unknown")[0] == ""


def test_chip_only_findings_produce_a_valuation():
    # The engine must work on a consumer chip, not only on a whole-genome VCF.
    chip = {"findings_with_economics": [
        {"finding": "CYP2C19 intermediate metabolizer",
         "category": "Pharmacogenomics", "confidence": "high"},
        {"finding": "CAD polygenic risk elevated", "category": "Polygenic Risk",
         "confidence": "moderate", "qaly_gain": 1.5},
        {"finding": "APOE e3/e4", "category": "Genotype",
         "confidence": "moderate", "qaly_gain": 1.5},
    ]}
    r = voi.analyze_value_of_information(chip, input_type="chip", n_mc=500)
    assert r["available"] is True, r.get("reason")
    assert r["n_findings"] == 3
    assert r["voi_expost_mean"] != 0


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


def test_evppi_bounded_by_evpi():
    # Theory: 0 <= EVPPI(any subset) <= EVPI. Violating this signals a broken estimator.
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=1500, seed=21)
    ev, evppi = r["evpi"], r["evppi"]
    assert evppi["available"] is True
    for row in evppi["by_parameter"]:
        assert row["evppi"] >= 0
        assert row["evppi"] <= max(ev["evpi"], 1) + 50    # tolerance for MC noise


def test_behavioural_present_bias_reduces_value():
    b = voi.analyze_behavioural(mean=40_000, sd=15_000, test_cost=300)
    assert b["available"] is True
    assert b["pv_hyperbolic"] <= b["pv_exponential"]     # present bias discounts more
    assert b["adoption_gap"] >= 0
    assert b["loss_aversion_lambda"] > 1                 # losses loom larger


def test_longevity_wired_into_voi_result():
    r = voi.analyze_value_of_information(_econ(), _cvr(), input_type="wgs",
                                        n_mc=800, seed=3)
    lon = r["longevity"]
    assert lon["available"] is True
    blended = [s["blended"] for s in lon["scenarios"]]
    assert blended == sorted(blended)     # longer life → more realised risk


def test_welfare_comparison_is_not_assumed_true():
    # The model must be able to conclude CENTRALISED wins — otherwise the local-wins
    # result is assumed by construction and carries no information.
    low_risk = voi.analyze_welfare_comparison(voi=25_000, test_cost=300,
                                              p_breach_annual=0.001)
    high_risk = voi.analyze_welfare_comparison(voi=25_000, test_cost=300,
                                               p_breach_annual=0.05)
    assert low_risk["local_preferred"] is False      # capability gap dominates
    assert high_risk["local_preferred"] is True      # privacy cost dominates


def test_welfare_breakeven_probability_is_the_pivot():
    w = voi.analyze_welfare_comparison(voi=25_000, test_cost=300)
    p_star = w["breakeven_annual_breach_prob"]
    assert 0.0 < p_star < 1.0
    just_below = voi.analyze_welfare_comparison(voi=25_000, test_cost=300,
                                                p_breach_annual=p_star * 0.5)
    just_above = voi.analyze_welfare_comparison(voi=25_000, test_cost=300,
                                                p_breach_annual=min(1.0, p_star * 2))
    assert just_below["local_preferred"] is False
    assert just_above["local_preferred"] is True


def test_welfare_privacy_cost_rises_with_exposure_horizon():
    # A genome cannot be revoked, so exposure hazard accumulates over the horizon.
    short = voi.analyze_welfare_comparison(voi=25_000, test_cost=300, horizon_years=5)
    long = voi.analyze_welfare_comparison(voi=25_000, test_cost=300, horizon_years=40)
    assert long["prob_exposure_over_horizon"] > short["prob_exposure_over_horizon"]


def test_welfare_access_channel_is_positive():
    # Higher participation under local analysis is itself a social-surplus gain.
    w = voi.analyze_welfare_comparison(voi=25_000, test_cost=300)
    assert w["access_channel_gain"] > 0
    assert w["social_surplus_local"] > w["social_surplus_central"]


def _nvr_gene(gene="BRCA1", **kw):
    """A predicted-pathogenic variant IN A NAMED GENE.

    `_nvr()` carries no gene, so before this fixture existed the gene-anchored
    branch of the predicted-variant path had no test at all — which is part of
    why it routed every prediction to one generic bucket unnoticed. The
    gene-less fixture is kept as its own case below, because the two exercise
    different paths.
    """
    f = {"chrom": "17", "pos": 43045711, "am_score": 0.95,
         "confidence": "higher", "gene": gene}
    f.update(kw)
    return {"available": True, "buckets": {"predicted_pathogenic_rare": [f]}}


def test_predicted_variant_in_anchored_gene_is_priced_and_downweighted():
    r = voi.analyze_value_of_information(_econ(), None, _nvr_gene(),
                                        input_type="wgs", n_mc=200)
    pred = [row for row in r["nmb_rows"] if "predicted" in row["label"]]
    assert pred, "an anchored predicted finding must reach the NMB table"


def test_predicted_variant_without_an_anchor_is_reported_not_priced():
    # PATH C. `_nvr()` names no gene, so there is no condition anchor. It must
    # be reported as unvaluable rather than priced off a generic bucket — the
    # behaviour that let nine unrelated findings share one figure.
    r = voi.analyze_value_of_information(_econ(), None, _nvr(),
                                        input_type="wgs", n_mc=200)
    assert not [row for row in r["nmb_rows"] if "predicted" in row["label"]]
    reasons = " ".join(str(u.get("reason", "")) for u in
                       (r.get("unvalued_findings") or []))
    assert "no registry-backed condition anchor" in reasons


def test_predicted_variant_valued_below_a_clinvar_assertion():
    """THE ANTI-INFLATION GUARD. A sequence model's call must never be worth as
    much as a curated clinical assertion for the same gene.

    Each finding is collected ALONE so it is the only member of its condition
    pool: ConditionPool.qaly_loss() returns max(overrides), so a Path A and a
    Path B finding sharing a pool would let Path A's larger override mask Path
    B entirely — the assertion would then pass because Path B was DISCARDED
    rather than down-weighted, which is not what this test is for.
    """
    wtp, rate = 100_000.0, 0.03
    asserted = [f for f in voi._collect(
        None, {"available": True, "buckets": {"actionable": [{"gene": "BRCA1"}]}},
        None) if "BRCA1" in f["label"]]
    predicted = [f for f in voi._collect(None, None, _nvr_gene("BRCA1"))
                 if "BRCA1" in f["label"]]
    assert len(asserted) == 1 and len(predicted) == 1

    a, b = asserted[0], predicted[0]
    nmb_a = voi._finding_nmb(a, wtp, rate)[0]
    nmb_b = voi._finding_nmb(b, wtp, rate)[0]
    assert nmb_b < nmb_a, (
        f"predicted BRCA1 ({nmb_b:,.0f}) must be worth strictly less than "
        f"ClinVar-asserted BRCA1 ({nmb_a:,.0f})")

    # Both discounts are present and neither is double-applied.
    ppv = float(voi._econ_params.value("predictor_ppv_no_clinvar"))
    assert b["haircut"] == ppv, "haircut must be the registered PPV"
    assert b["penetrance_corrected"] < b["penetrance_literature"], (
        "the ascertainment correction must apply to Path B too — an "
        "unconfirmed prediction is more incidentally identified, not less")
    # The anchor is NOT pre-multiplied: _finding_nmb applies the haircut.
    assert b["qaly"] == a["qaly"], (
        "the QALY anchor must be identical on both paths; PPV is applied by "
        "_finding_nmb via the haircut, and pre-multiplying would apply it twice")
    # Horizon deliberately shorter on the unconfirmed path.
    assert b["horizon"] == 25 and a["horizon"] == 30


def test_confirmation_cost_is_charged_on_predicted_variants():
    b = next(f for f in voi._collect(None, None, _nvr_gene("BRCA1"))
             if "BRCA1" in f["label"])
    expected = float(voi._econ_params.value(
        "intervention_cost_predicted_variant"))
    assert b["intervention"] == expected > 0
    # Charged in full: _finding_nmb haircuts the benefit, then subtracts the
    # intervention. A confirmation that cost nothing would make chasing a
    # prediction look free.
    assert voi._finding_nmb(b, 100_000.0, 0.03)[3] == expected
    assert "onfirm" in b["action"]
