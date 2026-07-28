"""Unit tests for environmental_optimization.py."""

from __future__ import annotations
import pandas as pd
import environmental_optimization as eo


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def test_circadian_evening_type_gets_morning_light_protocol() -> None:
    r = eo.analyze_environmental_optimization(_df({"rs1801260": "CC"}))
    c = r["circadian"]
    assert "Evening" in c["lean"]
    assert any("light within" in p.lower() or "bright light" in p.lower()
               for p in c["protocol"])


def test_circadian_morning_type() -> None:
    r = eo.analyze_environmental_optimization(_df({"rs1801260": "TT"}))
    assert "Morning" in r["circadian"]["lean"]


def test_exercise_power_lean_from_actn3_rr_plus_ace_dd() -> None:
    r = eo.analyze_environmental_optimization(_df({"rs1815739": "CC", "rs4343": "GG"}))
    assert "Power" in r["exercise"]["lean"]
    assert r["exercise"]["power_score"] > r["exercise"]["endurance_score"]


def test_exercise_endurance_lean_from_actn3_xx() -> None:
    r = eo.analyze_environmental_optimization(_df({"rs1815739": "TT", "rs4343": "AA"}))
    assert "Endurance" in r["exercise"]["lean"]


def test_vitamin_d_latitude_drives_months() -> None:
    high = eo.analyze_environmental_optimization(_df({}), latitude=55.0)
    low = eo.analyze_environmental_optimization(_df({}), latitude=10.0)
    assert "October" in high["vitamin_d"]["protocol"][0] or high["vitamin_d"]["supplement_months"]
    # tropical latitude → little/no supplementation from latitude alone
    assert "not needed" in low["vitamin_d"]["supplement_months"].lower()


def test_vitamin_d_tendency_higher_with_low_alleles() -> None:
    r = eo.analyze_environmental_optimization(
        _df({"rs2282679": "GG", "rs10741657": "AA", "rs12785878": "GG", "rs2228570": "AA"}),
        latitude=43.0)
    assert "Higher" in r["vitamin_d"]["tendency"]


def test_available_true_even_with_no_genotypes_due_to_latitude_guidance() -> None:
    r = eo.analyze_environmental_optimization(_df({}), latitude=43.0)
    # vitamin-D latitude guidance always available
    assert r["available"] is True
    assert r["vitamin_d"] is not None
