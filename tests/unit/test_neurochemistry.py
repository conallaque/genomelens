"""Unit tests for neurochemistry.py."""

from __future__ import annotations

import pandas as pd

import neurochemistry as nc


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def test_comt_val_met_heterozygous_middle() -> None:
    r = nc.analyze_neurochemistry(_df({"rs4680": "AG"}))
    comt = next(f for f in r["findings"] if f["gene"] == "COMT")
    assert "middle" in comt["phenotype"].lower() or "adaptive" in comt["phenotype"].lower()
    assert r["composite"]["comt_class"] == "middle"


def test_comt_val_val_warrior() -> None:
    r = nc.analyze_neurochemistry(_df({"rs4680": "GG"}))
    assert r["composite"]["comt_class"] == "warrior"


def test_comt_met_met_worrier() -> None:
    r = nc.analyze_neurochemistry(_df({"rs4680": "AA"}))
    assert r["composite"]["comt_class"] == "worrier"


def test_maoa_high_activity_flagged() -> None:
    r = nc.analyze_neurochemistry(_df({"rs6323": "TT"}))
    assert r["composite"]["maoa_class"] == "MAOA-H"


def test_bdnf_val_val_full_plasticity() -> None:
    r = nc.analyze_neurochemistry(_df({"rs6265": "CC"}))
    assert "Val/Val" in r["composite"]["bdnf_class"]


def test_chrna5_flags_smoking_risk() -> None:
    r = nc.analyze_neurochemistry(_df({"rs16969968": "AG"}))
    subs = r["composite"]["substance_flags"]
    assert any("smoking" in s.lower() or "chrna5" in s.lower() for s in subs)


def test_composite_generates_full_recommendations() -> None:
    r = nc.analyze_neurochemistry(_df({
        "rs4680": "AG", "rs6323": "TT", "rs6265": "CC",
        "rs1799971": "AG", "rs16969968": "AG",
    }))
    c = r["composite"]
    for key in ("stress_response_profile", "plasticity_tier",
                "stimulant_response", "ssri_response", "caffeine_protocol",
                "meditation_fit", "career_neurotype"):
        assert c.get(key), f"missing composite field: {key}"


def test_empty_input() -> None:
    r = nc.analyze_neurochemistry(_df({}))
    assert r["available"] is False
