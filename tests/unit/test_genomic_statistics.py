"""Unit tests for genomic_statistics.py — quantitative-genetics corrections.

These assert *theoretical properties* (calibration, monotonicity, competing-risk
bounds), not memorised outputs, so they stay meaningful if constants are retuned.
"""
from __future__ import annotations

import pytest

import genomic_statistics as gs

np = pytest.importorskip("numpy")


# ── liability-threshold model ────────────────────────────────────────────────

def test_liability_model_calibrates_to_prevalence_in_expectation():
    # Risk is convex in liability, so the MEAN over the PRS distribution equals the
    # prevalence while the MEDIAN individual sits below it (Jensen's inequality).
    K = 0.06
    qs = np.linspace(0.0005, 0.9995, 2000)
    mean_risk = float(np.mean([gs.liability_threshold_risk(q, K)["absolute_risk"]
                               for q in qs]))
    assert abs(mean_risk - K) < 0.003
    assert gs.liability_threshold_risk(0.5, K)["absolute_risk"] < K


def test_liability_risk_increases_with_percentile():
    K = 0.06
    risks = [gs.liability_threshold_risk(q, K)["absolute_risk"]
             for q in (0.10, 0.50, 0.90, 0.99)]
    assert risks == sorted(risks)


def test_liability_reports_absolute_not_just_relative():
    r = gs.liability_threshold_risk(0.99, 0.001)   # big RR, tiny absolute risk
    assert r["relative_risk"] > 3
    assert r["absolute_risk"] < 0.05               # still small in absolute terms
    assert r["absolute_risk_increase"] is not None


# ── age-dependent penetrance + competing risks ───────────────────────────────

def test_competing_risks_never_exceed_naive_km():
    for age in (20, 40, 60, 80):
        r = gs.age_dependent_penetrance(0.55, current_age=age)
        assert r["remaining_lifetime_risk"] <= r["naive_km_risk"] + 1e-9


def test_remaining_risk_declines_with_age():
    vals = [gs.age_dependent_penetrance(0.55, current_age=a)["remaining_lifetime_risk"]
            for a in (20, 35, 50, 65, 80)]
    assert vals == sorted(vals, reverse=True)


def test_longer_lifespan_raises_realised_risk():
    # Anchoring penetrance to a fixed reference age means extending the lifespan
    # assumption ADDS post-anchor risk rather than diluting the hazard.
    s = gs.longevity_sensitivity(0.55, current_age=35)
    blended = [row["blended"] for row in s["scenarios"]]
    assert blended == sorted(blended)
    assert s["scenarios"][0]["scenario"].startswith("2025")


# ── empirical Bayes / James–Stein ────────────────────────────────────────────

def test_empirical_bayes_shrinks_noisiest_estimate_most():
    eff = [0.20, 0.15, 0.05, 0.30, 0.02]
    ses = [0.02, 0.05, 0.04, 0.10, 0.03]
    eb = gs.empirical_bayes_shrinkage(eff, ses)
    assert eb["available"] is True
    moves = [abs(r - s) for r, s in zip(eb["effects_raw"], eb["effects_shrunk"])]
    assert moves.index(max(moves)) == ses.index(max(ses))   # noisiest moves most
    assert 0.0 <= eb["mean_shrinkage_weight"] <= 1.0


# ── PRS ancestry portability ─────────────────────────────────────────────────

def test_portability_pulls_toward_base_rate():
    eur = gs.prs_portability(0.95, 0.06, ancestry="european")
    afr = gs.prs_portability(0.95, 0.06, ancestry="african")
    assert afr["risk_ancestry_adjusted"] < eur["risk_ancestry_adjusted"]
    assert afr["attenuation"] > 0
    assert afr["r2_retained_fraction"] < eur["r2_retained_fraction"]


def test_unknown_ancestry_falls_back_safely():
    p = gs.prs_portability(0.9, 0.05, ancestry="not-a-real-population")
    assert p["available"] is True
    assert 0 < p["r2_retained_fraction"] <= 1.0


def test_correction_chain_reports_each_step():
    s = gs.summarise_genomic_corrections()
    assert s["available"] is True
    assert len(s["chain"]) == 3
    assert all("value" in step for step in s["chain"])
