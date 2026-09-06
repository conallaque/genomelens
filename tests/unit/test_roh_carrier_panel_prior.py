"""The one defensible economic link from ROH: a CATEGORICAL carrier-panel prior.

Runs of homozygosity previously contributed nothing to the economic model — the
result was computed, rendered, and never seen by either economics module. The
link wired here names a real decision (buy an expanded / ancestry-matched
carrier panel, or a minimal one) without attaching a dollar figure to it.

The central contract these tests defend: **this analysis must never produce a
monetary value.** A person's own F_ROH describes their parents' relatedness, not
their children's, and no published study maps F_ROH onto incremental panel
yield — so any number here would be fabricated rather than estimated.
"""
from __future__ import annotations

from econ import value_of_information as voi


def _roh(f_roh=0.0, f_long=0.0, n_long=0, tier="no_recent_relatedness"):
    return {"f_roh": f_roh, "f_roh_long": f_long, "n_long": n_long,
            "context_tier": tier}


# ── the non-monetization contract ────────────────────────────────────────────

def test_result_never_contains_a_monetary_field():
    for r in (_roh(), _roh(0.02), _roh(0.06, 0.05, 4, "first_cousin_or_closer")):
        out = voi.assess_carrier_panel_prior(r)
        bad = [k for k in out
               if any(s in k.lower() for s in
                      ("cost", "value", "usd", "dollar", "nmb", "price",
                       "saving", "benefit", "icer", "qaly"))]
        assert not bad, f"monetary-looking field(s) leaked: {bad}"


def test_no_numeric_field_could_be_read_as_currency():
    out = voi.assess_carrier_panel_prior(_roh(0.06, 0.05, 4))
    for k, v in out.items():
        if isinstance(v, int | float) and not isinstance(v, bool):
            assert v <= 100, f"{k}={v} is large enough to be mistaken for money"


def test_monetized_flag_is_false_and_explained():
    out = voi.assess_carrier_panel_prior(_roh(0.03))
    assert out["monetized"] is False
    assert len(out["why_not_monetized"]) > 80


def test_no_risk_score_is_emitted():
    # Consanguinity must never be scored — the framing is informational.
    out = voi.assess_carrier_panel_prior(_roh(0.06, 0.05, 4))
    assert not [k for k in out if "risk" in k.lower() and k != "why_not_monetized"]


# ── tiering behavior ────────────────────────────────────────────────────────

def test_outbred_profile_recommends_no_change():
    out = voi.assess_carrier_panel_prior(_roh(0.002))
    assert out["available"] and out["tier"] == "none"


def test_threshold_clears_the_observed_outbred_range():
    # CALIBRATION GUARD. Eleven real public genomes (PGP Harvard) measured
    # F_ROH 0.0113-0.0161 with zero long runs. An earlier 0.010 cutoff fired on
    # all eleven, making the founder tier meaningless. The threshold must stay
    # above that observed range or the recommendation becomes universal.
    assert voi.ROH_FOUNDER_F_THRESHOLD > 0.0161, (
        "threshold has dropped back into the ordinary outbred range")


def test_observed_outbred_values_do_not_trigger_founder_tier():
    for f_roh in (0.01134, 0.01243, 0.01388, 0.01485, 0.01614):   # real samples
        out = voi.assess_carrier_panel_prior(_roh(f_roh))
        assert out["tier"] == "none", f"F_ROH {f_roh} should not flag founder"


def test_founder_background_is_distinguished_from_recent_relatedness():
    # High total burden but NO long runs -> founder tier, not "recent".
    founder = voi.assess_carrier_panel_prior(_roh(0.030, 0.0, 0))
    assert founder["tier"] == "founder"
    assert "related" not in founder["recommendation"].lower()
    # Long runs present -> recent tier.
    recent = voi.assess_carrier_panel_prior(_roh(0.055, 0.050, 4))
    assert recent["tier"] == "recent"


def test_long_roh_dominates_even_at_modest_total_burden():
    out = voi.assess_carrier_panel_prior(_roh(0.009, 0.008, 1))
    assert out["tier"] == "recent"


def test_every_tier_names_the_decision_and_cites_a_source():
    for r in (_roh(0.001), _roh(0.03), _roh(0.06, 0.05, 3)):
        out = voi.assess_carrier_panel_prior(r)
        assert "panel" in out["decision"].lower()
        assert "Kirin" in out["src"]
        assert len(out["rationale"]) > 80


def test_founder_tier_does_not_imply_parental_relatedness():
    # The false positive the old blended-F_ROH tiers produced.
    out = voi.assess_carrier_panel_prior(_roh(0.030, 0.0, 0))
    assert "says nothing about" in out["rationale"].lower()


# ── degradation ──────────────────────────────────────────────────────────────

def test_missing_or_unavailable_roh_degrades_cleanly():
    for r in (None, {}, {"context_tier": "unavailable"}):
        out = voi.assess_carrier_panel_prior(r)
        assert out["available"] is False
        assert out["monetized"] is False


# ── integration with the economic model ──────────────────────────────────────

def _econ():
    return {"findings_with_economics": [
        {"finding": "CYP2C19 poor metaboliser", "category": "Pharmacogenomics",
         "confidence": "high", "outcome_value": 5000, "qaly_gain": 0.3}]}


def test_panel_prior_is_attached_to_the_voi_result():
    r = voi.analyze_value_of_information(_econ(), None, None,
                                         roh_result=_roh(0.03), n_mc=200, seed=1)
    assert r["carrier_panel_prior"]["available"] is True
    assert r["carrier_panel_prior"]["tier"] == "founder"


def test_panel_prior_does_not_change_any_dollar_total():
    kw = dict(economics_result=_econ(), clinical_variants_result=None,
              novel_variants_result=None, n_mc=600, seed=5)
    without = voi.analyze_value_of_information(**kw)
    with_roh = voi.analyze_value_of_information(
        roh_result=_roh(0.06, 0.05, 4), **kw)
    for key in ("voi_expost_point", "voi_expost_mean", "marginal_chip_to_wgs",
                "chip_value_point", "n_findings"):
        assert without.get(key) == with_roh.get(key), (
            f"{key} changed when ROH was supplied — the prior must be "
            f"qualitative and must not enter any total")


def test_panel_prior_stays_out_of_the_nmb_table():
    r = voi.analyze_value_of_information(_econ(), None, None,
                                         roh_result=_roh(0.06, 0.05, 4),
                                         n_mc=200, seed=1)
    labels = " ".join(row["label"].lower() for row in r["nmb_rows"])
    assert "homozygosity" not in labels and "roh" not in labels.split()


def test_missing_roh_does_not_break_the_model():
    r = voi.analyze_value_of_information(_econ(), None, None, n_mc=200, seed=1)
    assert r["available"] is True
    assert r["carrier_panel_prior"]["available"] is False
    assert r["fully_computed"] is True     # absent ROH is not a degradation
