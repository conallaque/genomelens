"""Unit tests for deep_ancestry.py (Neanderthal / ancient-pop / N-S axis /
haplogroup timelines)."""

from __future__ import annotations

import pandas as pd

import deep_ancestry as da


def test_neanderthal_scores_present_alleles() -> None:
    # Two known Neanderthal loci: rs11024102 (C = N-derived), rs2549794 (C = N).
    df = pd.DataFrame({"genotype": {
        "rs11024102": "CC",
        "rs2549794": "TC",
        "rs35044562": "AA",   # non-N
    }})
    res = da.analyze_neanderthal(df)
    assert res["available"]
    # 2 (rs11024102 CC) + 1 (rs2549794 TC) + 0 = 3 N-alleles over 6 possible
    assert res["n_neanderthal_alleles"] == 3
    assert res["max_possible"] == 6
    assert abs(res["affinity"] - 0.5) < 0.001
    assert res["tier"] in ("Above average", "Average non-African", "High")


def test_neanderthal_empty_gracefully() -> None:
    res = da.analyze_neanderthal(pd.DataFrame({"genotype": {}}))
    assert res["available"] is False


def test_ancient_populations_ranks_yamnaya_when_carrying_lct_persistence() -> None:
    df = pd.DataFrame({"genotype": {
        "rs4988235": "AA",     # LCT persistence homozygous (Yamnaya signature)
        "rs16891982": "GG",    # SLC45A2 light skin
        "rs1426654": "AA",     # SLC24A5 EEF marker
        "rs12913832": "AA",    # HERC2 brown eyes — no WHG blue-eye copies
    }})
    res = da.analyze_ancient_populations(df)
    assert res["available"]
    # Yamnaya should score highest with LCT persistence + SLC45A2
    top = res["populations"][0]
    assert top["short"] == "Yamnaya"


def test_north_south_european_axis_leans_northern_with_lct() -> None:
    df = pd.DataFrame({"genotype": {
        "rs4988235": "AA", "rs12913832": "GG", "rs1042602": "AA", "rs1805007": "CC",
    }})
    res = da.analyze_north_south_europe(df)
    assert res["available"]
    # Full Northern signature except MC1R → high index
    assert res["index"] > 0.7
    assert "Northern" in res["lean"]


def test_haplogroup_timeline_longest_prefix_match() -> None:
    y = {"terminal_haplogroup": "T1a1a"}
    mt = {"haplogroup": "V"}
    tl = da.build_haplogroup_timeline(y, mt)
    # T1a1a should resolve via "T" (or a longer key if present)
    assert tl["y"] is not None and tl["y"]["tmrca_kya"] > 0
    assert tl["mt"] is not None and tl["mt"]["origin"]


def test_full_deep_ancestry_end_to_end() -> None:
    df = pd.DataFrame({"genotype": {
        "rs4988235": "AA", "rs12913832": "AG", "rs1426654": "AA",
        "rs16891982": "GG", "rs11024102": "CT", "rs2549794": "CC",
    }})
    y = {"terminal_haplogroup": "R1b"}
    mt = {"haplogroup": "H"}
    res = da.analyze_deep_ancestry(df, y_result=y, mt_result=mt)
    assert res["available"]
    assert res["neanderthal"]["available"]
    assert res["ancient_populations"]["available"]
    assert res["european_axis"]["available"]
    assert res["haplogroup_timeline"]["y"] is not None
