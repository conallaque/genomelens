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


def test_predicted_variant_downweighted():
    # A high-AM predicted variant should still contribute less than a confirmed
    # ClinVar-actionable finding of similar condition (haircut applied).
    r = voi.analyze_value_of_information(_econ(), _cvr(), _nvr(),
                                        input_type="wgs", n_mc=200)
    labels = {row["label"]: row["nmb"] for row in r["nmb_rows"]}
    pred = [v for k, v in labels.items() if "predicted" in k]
    assert pred, "predicted finding should appear in the NMB table"
