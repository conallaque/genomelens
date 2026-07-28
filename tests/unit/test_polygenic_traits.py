"""Unit tests for polygenic_traits.py — honest genotype-level trait reporting."""

from __future__ import annotations
import json
import pandas as pd
import polygenic_traits as pt


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def test_tas2r38_strong_taster() -> None:
    r = pt.analyze_polygenic_traits(_df({"rs713598": "GG"}))
    t = next(f for f in r["findings"] if f["gene"] == "TAS2R38")
    assert "taster" in t["call"].lower()


def test_tas2r38_non_taster() -> None:
    r = pt.analyze_polygenic_traits(_df({"rs713598": "CC"}))
    t = next(f for f in r["findings"] if f["gene"] == "TAS2R38")
    assert "non-taster" in t["call"].lower()


def test_herc2_blue_eye_call() -> None:
    r = pt.analyze_polygenic_traits(_df({"rs12913832": "GG"}))
    e = next(f for f in r["findings"] if "HERC2" in f["gene"])
    assert "blue" in e["call"].lower()


def test_abcc11_dry_earwax() -> None:
    r = pt.analyze_polygenic_traits(_df({"rs17822931": "TT"}))
    e = next(f for f in r["findings"] if f["gene"] == "ABCC11")
    assert "dry" in e["call"].lower()


def test_no_percentile_fields_emitted() -> None:
    """The whole point of the module: it must never emit a numeric percentile
    or score field for any trait finding (only genotype-level calls)."""
    r = pt.analyze_polygenic_traits(_df({"rs713598": "GG", "rs12913832": "AG"}))
    for f in r["findings"]:
        assert "percentile" not in f
        assert "score" not in f
        assert "z_score" not in f


def test_polygenic_notes_have_no_scores() -> None:
    r = pt.analyze_polygenic_traits(_df({}))
    traits = {n["trait"] for n in r["polygenic_notes"]}
    assert "Height" in traits
    assert any("Cognitive" in t for t in traits)
    assert any("Personality" in t for t in traits)
    for n in r["polygenic_notes"]:
        # explanation only — no numeric fields
        assert "percentile" not in n
        assert "score" not in n
        assert n["why_no_number"] and n["honest_statement"]


def test_empty_input() -> None:
    r = pt.analyze_polygenic_traits(_df({}))
    assert r["available"] is False
    # even with no genotypes, the no-score polygenic notes are still present
    assert len(r["polygenic_notes"]) == 3
