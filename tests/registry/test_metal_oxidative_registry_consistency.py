"""Cross-module consistency for metal_oxidative.py vs snp_registry."""

from __future__ import annotations

import pandas as pd

from wellness import metal_oxidative


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
    """Lock the gene set this section was created to cover. GBA, G6PD and HFE
    are reused from existing carrier records; the metallothionein/ATP7B/CAT
    SNPs were added in V9; the ZIP transporters (SLC39A8/A14), ABCG1 and the
    glutathione-S-transferases (GSTM1/GSTT1) were added later."""
    expected = {"rs34637584", "rs76763715", "rs1061472",
                "rs8052394", "rs28366003", "rs1001179", "rs1050828",
                "rs1800562", "rs1799945", "rs13107325", "rs896378",
                "rs1893590", "rs4147565", "rs4630"}
    assert expected.issubset(set(metal_oxidative.audit_against_registry()["registered"]))


def test_all_categories_produced_with_risk_genotypes() -> None:
    df = pd.DataFrame({"genotype": {
        "rs34637584": "GA", "rs76763715": "TC", "rs1061472": "AG",
        "rs8052394": "GG", "rs28366003": "AG", "rs1001179": "TT",
        "rs1050828": "CT",
        # HFE iron, ZIP transporters, ABCG1, and GST detox findings.
        "rs1800562": "GA", "rs1799945": "CG", "rs13107325": "CT",
        "rs896378": "TC", "rs1893590": "AC", "rs4147565": "GA",
        "rs4630": "GA",
    }})
    res = metal_oxidative.analyze_metal_oxidative(df)
    # Paired-gene analyzers each collapse into one finding: MT1A+MT2A,
    # SLC39A8+SLC39A14, and GSTM1+GSTT1. The 13 SNPs above therefore yield 10
    # findings (LRRK2, GBA, ATP7B, MT, HFE, ZIP, ABCG1, CAT, G6PD, GST).
    assert res["n_predictions"] == 10
    assert set(res["categories"]) == {
        "Neurodegeneration (Parkinson's)", "Metal Handling", "Oxidative Defense",
    }


def test_empty_genotypes_produce_no_findings() -> None:
    res = metal_oxidative.analyze_metal_oxidative(pd.DataFrame({"genotype": {}}))
    assert res["n_predictions"] == 0
