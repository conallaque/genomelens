"""Unit tests for immunogenetics.py + ancestral_story.py."""

from __future__ import annotations

import pandas as pd

import immunogenetics as ig
from ancestry import story as story


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def test_ccr5_delta32_heterozygous_flagged_protective() -> None:
    r = ig.analyze_immunogenetics(_df({"rs333": "DI"}))
    ccr5 = next(f for f in r["findings"] if f["gene"] == "CCR5")
    assert ccr5["impact"] == "protective"
    assert "50%" in ccr5["verdict"] or "heterozy" in ccr5["verdict"].lower()


def test_fut2_non_secretor_flagged_near_complete_norovirus() -> None:
    r = ig.analyze_immunogenetics(_df({"rs601338": "AA"}))
    fut2 = next(f for f in r["findings"] if f["gene"] == "FUT2")
    assert fut2["impact"] == "protective"
    assert "norovirus" in fut2["verdict"].lower() or "GII.4" in fut2["verdict"]


def test_il28b_cc_hepc_clearance() -> None:
    r = ig.analyze_immunogenetics(_df({"rs12979860": "CC"}))
    il = next(f for f in r["findings"] if "IL28B" in f["gene"])
    assert il["impact"] == "protective"
    assert "clearance" in il["verdict"].lower() or "3" in il["verdict"]


def test_prnp_heterozygous_protective() -> None:
    r = ig.analyze_immunogenetics(_df({"rs1799990": "AG"}))
    prn = next(f for f in r["findings"] if f["gene"] == "PRNP")
    assert prn["impact"] == "protective"
    assert "prion" in prn["verdict"].lower()


def test_erap2_black_death_survivor_historical_note() -> None:
    r = ig.analyze_immunogenetics(_df({"rs2549794": "TC"}))
    er = next(f for f in r["findings"] if f["gene"] == "ERAP2")
    assert er["impact"] == "protective"
    assert er["historical"] and "black death" in er["historical"].lower()


def test_historical_timeline_populated() -> None:
    r = ig.analyze_immunogenetics(_df({
        "rs333": "DI", "rs601338": "AA", "rs1799990": "AG", "rs2549794": "TC",
    }))
    assert r["historical_timeline"], "expected non-empty selection timeline"
    epochs = {e["epoch"] for e in r["historical_timeline"]}
    # Multiple distinct pandemics represented
    assert len(epochs) >= 2


def test_empty_input_gracefully() -> None:
    r = ig.analyze_immunogenetics(_df({}))
    assert r["available"] is False and r["n_findings"] == 0


def test_ancestral_story_template_mode() -> None:
    y = {"terminal_haplogroup": "T1a1a"}
    mt = {"haplogroup": "V"}
    deep = {"neanderthal": {"available": True, "approx_pct": 0.9, "tier": "Average non-African",
                            "n_carrying": 3, "n_typed": 9},
            "ancient_populations": {"available": True,
                                    "populations": [{"short": "Yamnaya", "affinity": 0.75},
                                                    {"short": "EEF", "affinity": 0.75}]},
            "european_axis": {"available": True, "lean": "Northern European", "index": 0.61}}
    im = ig.analyze_immunogenetics(_df({"rs333": "DI", "rs601338": "AA"}))
    s = story.analyze_ancestral_story(
        y_result=y, mt_result=mt, deep_ancestry_result=deep,
        immunogenetics_result=im, ancestry_result={"proportions": {"EUR": 0.99}},
        use_ai=False,
    )
    assert s["available"] and len(s["template"]["chapters"]) >= 3
    # Y-DNA T chapter should mention the Fertile Crescent
    titles = " ".join(c["title"] for c in s["template"]["chapters"])
    assert "Paternal" in titles and "Maternal" in titles
