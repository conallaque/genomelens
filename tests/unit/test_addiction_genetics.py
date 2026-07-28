"""Unit tests for addiction_genetics.py."""

from __future__ import annotations
import pandas as pd
import addiction_genetics as ag


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def test_adh1b_homozygous_protective_flag() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs1229984": "AA"}))
    adh = next(f for f in r["findings"] if f["gene"] == "ADH1B")
    assert adh["impact"] == "protective"
    assert "strong" in adh["verdict"].lower()


def test_aldh2_flush_variant_flags_esophageal_risk() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs671": "GA"}))
    aldh = next(f for f in r["findings"] if f["gene"] == "ALDH2")
    assert aldh["impact"] == "protective"
    flags = " ".join(f["title"] for f in r["composite"]["clinical_flags"])
    assert "esophageal" in flags.lower() or "cancer" in flags.lower()


def test_gabra2_risk_allele_flagged() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs279858": "CC"}))
    g = next(f for f in r["findings"] if f["gene"] == "GABRA2")
    assert g["impact"] == "susceptible"


def test_oprm1_g_carrier_flags_naltrexone_and_opioid() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs1799971": "AG"}))
    flag_titles = " ".join(f["title"] for f in r["composite"]["clinical_flags"])
    assert "naltrexone" in flag_titles.lower()
    assert "opioid dosing" in flag_titles.lower()


def test_chrna5_flags_never_smoke() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs16969968": "AG"}))
    flag_titles = " ".join(f["title"] for f in r["composite"]["clinical_flags"])
    assert "smok" in flag_titles.lower()


def test_alcohol_tier_strong_protection_with_adh1b_plus_aldh2() -> None:
    r = ag.analyze_addiction_genetics(_df({"rs1229984": "AA", "rs671": "AA"}))
    assert r["composite"]["alcohol_tier"] == "Strongly protected"


def test_empty_input() -> None:
    r = ag.analyze_addiction_genetics(_df({}))
    assert r["available"] is False
