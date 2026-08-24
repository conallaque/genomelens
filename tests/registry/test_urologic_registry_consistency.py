"""Cross-module consistency for urologic.py vs snp_registry."""

from __future__ import annotations

import pandas as pd

from wellness import urologic


def test_audit_returns_well_structured_dict() -> None:
    audit = urologic.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}


def test_every_referenced_rsid_is_registered() -> None:
    """No urologic rsID may be missing from the registry — otherwise a typo
    would silently drop findings from the panel."""
    audit = urologic.audit_against_registry()
    assert audit["missing"] == [], (
        f"Unregistered rsIDs in urologic.py: {audit['missing']}. "
        "Add them to snp_registry._RECORDS."
    )


def test_panel_covers_expected_genes() -> None:
    """Lock the SNP set the panel was created to cover: OAB (ADRB3, CHRM3),
    BPH / 5α-reductase (SRD5A2), prostate cancer (HOXB13, 8q24, MSMB), kidney
    stones (CLDN14, SLC34A1, CASR), testicular (KITLG, SPRY4), SHBG."""
    expected = {
        "rs4994", "rs2229870",                                   # OAB
        "rs523349", "rs9282858",                                 # SRD5A2
        "rs138213197", "rs1447295", "rs6983267", "rs10993994",   # prostate cancer
        "rs219780", "rs4074995", "rs1042636", "rs2072499",       # kidney/PKD1
        "rs995030", "rs4324715",                                 # testicular
        "rs1799941",                                             # SHBG
    }
    assert expected.issubset(
        set(urologic.audit_against_registry()["registered"]))


def test_full_panel_runs_and_categorises() -> None:
    df = pd.DataFrame({"genotype": {
        # OAB / bladder
        "rs4994": "AG", "rs2229870": "CT",
        # NAT2 slow (cross-referenced from detox)
        "rs1801280": "CC", "rs1799930": "AA", "rs1799931": "GG",
        # Prostate — HOXB13 heterozygous positive (rare, high impact)
        "rs138213197": "CT",
        "rs523349": "CG", "rs9282858": "GG",
        "rs1447295": "AC", "rs6983267": "GG", "rs10993994": "TT",
        # Kidney stones
        "rs219780": "CC", "rs4074995": "AG", "rs1042636": "AG", "rs2072499": "AA",
        # Testicular
        "rs995030": "AG", "rs4324715": "TC",
        # Hormones
        "rs1799941": "AG",
    }})
    res = urologic.analyze_urologic(df)
    assert res["available"] and res["n_findings"] >= 12
    # All five sub-panels should be represented in the full-genotype scenario.
    cats = set(res["categories"])
    assert urologic.CAT_BLADDER in cats
    assert urologic.CAT_PROSTATE in cats
    assert urologic.CAT_STONES in cats
    assert urologic.CAT_TESTIS in cats
    assert urologic.CAT_HORMONES in cats
    # HOXB13 G84E heterozygous must be flagged with high confidence.
    hox = next(f for f in res["findings"] if f["rsid"] == "rs138213197")
    assert hox["impact"] == "higher-load" and hox["confidence"] == "high"


def test_empty_genotypes_produce_no_findings() -> None:
    res = urologic.analyze_urologic(pd.DataFrame({"genotype": {}}))
    assert res["available"] is False and res["n_findings"] == 0
