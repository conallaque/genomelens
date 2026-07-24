"""Cross-module consistency for detox.py vs snp_registry."""

from __future__ import annotations

import pandas as pd

import detox


def test_audit_returns_well_structured_dict() -> None:
    audit = detox.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}
    assert isinstance(audit["registered"], list)
    assert isinstance(audit["missing"], list)


def test_every_referenced_rsid_is_registered() -> None:
    """The detox panel references no rsID the registry doesn't know about —
    otherwise a typo would silently drop a finding."""
    audit = detox.audit_against_registry()
    assert audit["missing"] == [], (
        f"Unregistered rsIDs referenced by detox.py: {audit['missing']}. "
        "Add them to snp_registry._RECORDS."
    )


def test_panel_genes_all_present() -> None:
    """Lock the SNP set this section was created to cover: Phase I activation
    (CYP1A1/1B1/AHR), Phase II conjugation (EPHX1, NAT2, NQO1, GSTP1), the
    antioxidant axis (NRF2/SOD2/GPX1/CAT/HMOX1) and heavy-metal handling
    (ALAD, AS3MT, PON1, MT1A/MT2A)."""
    expected = {
        "rs1048943", "rs4646903", "rs1056836", "rs2066853",   # Phase I
        "rs1051740", "rs2234922", "rs1801280", "rs1799930",   # Phase II
        "rs1799931", "rs1800566", "rs1695",
        "rs6721961", "rs4880", "rs1050450", "rs1001179", "rs2071746",  # antioxidant
        "rs1800435", "rs11191439", "rs662", "rs8052394", "rs28366003",  # metals
    }
    assert expected.issubset(set(detox.audit_against_registry()["registered"]))


def test_full_panel_produces_all_domains_and_a_score() -> None:
    df = pd.DataFrame({"genotype": {
        # Phase I (fast activation)
        "rs1048943": "AG", "rs4646903": "TC", "rs1056836": "CG", "rs2066853": "GA",
        # Phase II (slow clearance)
        "rs1051740": "CC", "rs2234922": "AA", "rs1801280": "CC", "rs1799930": "AA",
        "rs1799931": "GG", "rs1800566": "TT", "rs1695": "GG",
        "rs4147565": "GA", "rs4630": "GA",
        # Antioxidant (reduced)
        "rs6721961": "TT", "rs4880": "TT", "rs1050450": "TT", "rs1001179": "TT",
        "rs2071746": "AT",
        # Metals
        "rs1800435": "CG", "rs11191439": "TC", "rs662": "AG",
        "rs8052394": "GG", "rs28366003": "GG",
    }})
    res = detox.analyze_detox(df)
    assert res["available"] is True
    assert res["n_findings"] >= 15
    assert set(res["by_category"].keys()) == {
        detox.CAT_ACTIVATION, detox.CAT_CONJUGATION,
        detox.CAT_ANTIOX, detox.CAT_METAL,
    }
    # Fast Phase I + slow Phase II/antioxidant ⇒ the high-susceptibility mismatch.
    sr = res["smoke_resilience"]
    assert sr["activate_but_dont_clear"] is True
    assert sr["tier"] == "Higher susceptibility"
    # Protocol should escalate the NRF2/glutathione levers to priority.
    assert any(n["emphasis"] == "high" for n in res["protocol"]["nutrition"])


def test_empty_genotypes_produce_no_findings() -> None:
    res = detox.analyze_detox(pd.DataFrame({"genotype": {}}))
    assert res["available"] is False
    assert res["n_findings"] == 0
