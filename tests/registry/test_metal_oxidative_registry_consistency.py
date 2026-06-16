"""Cross-module consistency for metal_oxidative.py vs snp_registry."""

from __future__ import annotations

import pandas as pd

import metal_oxidative


def test_audit_returns_well_structured_dict() -> None:
    audit = metal_oxidative.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}
    assert isinstance(audit["registered"], list)
    assert isinstance(audit["missing"], list)


def test_every_referenced_rsid_is_registered() -> None:
    """The panel references no rsID that the registry doesn't know about —
    otherwise the section would silently emit zero findings."""
    audit = metal_oxidative.audit_against_registry()
    assert audit["missing"] == [], (
        f"Unregistered rsIDs referenced by metal_oxidative.py: {audit['missing']}. "
        "Add them to snp_registry._RECORDS."
    )


def test_panel_genes_all_present() -> None:
    """Lock the gene set this section was created to cover. GBA and G6PD are
    reused from existing carrier records; the rest were added in V9."""
    expected = {"rs34637584", "rs76763715", "rs1061472",
                "rs8052394", "rs28366003", "rs1001179", "rs1050828"}
    assert expected.issubset(set(metal_oxidative.audit_against_registry()["registered"]))


def test_all_categories_produced_with_risk_genotypes() -> None:
    df = pd.DataFrame({"genotype": {
        "rs34637584": "GA", "rs76763715": "TC", "rs1061472": "AG",
        "rs8052394": "GG", "rs28366003": "AG", "rs1001179": "TT",
        "rs1050828": "CT",
    }})
    res = metal_oxidative.analyze_metal_oxidative(df)
    # MT1A + MT2A collapse into one metallothionein finding → 6, not 7.
    assert res["n_predictions"] == 6
    assert set(res["categories"]) == {
        "Neurodegeneration (Parkinson's)", "Metal Handling", "Oxidative Defense",
    }


def test_empty_genotypes_produce_no_findings() -> None:
    res = metal_oxidative.analyze_metal_oxidative(pd.DataFrame({"genotype": {}}))
    assert res["n_predictions"] == 0
