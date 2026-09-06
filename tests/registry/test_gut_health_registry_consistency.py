"""Cross-module consistency for gut_health.py vs snp_registry."""

from __future__ import annotations

import pandas as pd

from wellness import gut_health


def test_audit_returns_well_structured_dict() -> None:
    audit = gut_health.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}
    assert isinstance(audit["registered"], list)
    assert isinstance(audit["missing"], list)


def test_every_referenced_rsid_is_registered() -> None:
    """The panel references no rsID the registry doesn't know about — otherwise
    a typo would silently emit a chip-gap card forever."""
    audit = gut_health.audit_against_registry()
    assert audit["missing"] == [], (
        f"Unregistered rsIDs referenced by gut_health.py: {audit['missing']}. "
        "Add them to snp_registry._RECORDS."
    )


def test_panel_genes_all_present() -> None:
    """Lock the SNP set this section was created to cover. FUT2 and the celiac
    HLA-DQ tags are reused from existing records; LCT/AOC1/NOD2/IL23R were added
    alongside this module."""
    expected = {"rs4988235", "rs602662", "rs2187668", "rs7454108",
                "rs10156191", "rs2066844", "rs11209026"}
    assert expected.issubset(set(gut_health.audit_against_registry()["registered"]))


def test_all_categories_produced_with_genotypes() -> None:
    df = pd.DataFrame({"genotype": {
        "rs4988235": "GG", "rs602662": "AA", "rs2187668": "CT",
        "rs7454108": "TT", "rs10156191": "CT", "rs2066844": "CT",
        "rs11209026": "GA",
    }})
    res = gut_health.analyze_gut_health(df)
    assert res["n_predictions"] == 6
    assert set(res["categories"]) == {
        "Carbohydrate Digestion", "Microbiome Shaping", "Gluten & Celiac Risk",
        "Food Intolerance", "Inflammatory-Bowel Predisposition",
    }


def test_strand_flip_lactase_persistence() -> None:
    """rs4988235 is the classic strand-flip trap: chips often report it as C/T
    (- strand) while the registry is canonical + strand G/A. The registry's
    strand-aware dose must read a - strand TT as two persistence alleles."""
    minus_persistent = pd.DataFrame({"genotype": {"rs4988235": "TT"}})  # - strand
    plus_non_persist = pd.DataFrame({"genotype": {"rs4988235": "GG"}})  # + strand, no A
    p_yes = gut_health.analyze_gut_health(minus_persistent)["predictions"][0]
    p_no = gut_health.analyze_gut_health(plus_non_persist)["predictions"][0]
    assert "non-persistence" not in p_yes["result"]
    assert "non-persistence" in p_no["result"]


def test_missing_snp_surfaces_chip_gap_not_silence() -> None:
    """An absent SNP must produce an explicit chip-gap card, never a silent drop."""
    res = gut_health.analyze_gut_health(pd.DataFrame({"genotype": {}}))
    assert res["n_predictions"] == 6  # all six become chip-gap cards
    assert all(p["confidence"] == "none" for p in res["predictions"])
    assert all("Not evaluable" in p["result"] for p in res["predictions"])
