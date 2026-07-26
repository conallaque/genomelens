"""Unit tests for blood_type.py."""

from __future__ import annotations

import pandas as pd

import blood_type as bt


def _df(genotypes: dict, rhd_calls: int = 40, rhd_total: int = 42) -> pd.DataFrame:
    """Build a DataFrame with the requested ABO/FUT2 genotypes plus a
    controllable amount of RHD-locus coverage."""
    rows = [{"rsid": k, "chrom": "9", "pos": 136_131_000, "genotype": v}
            for k, v in genotypes.items()]
    for i in range(rhd_total):
        rows.append({"rsid": f"rhd_probe_{i}", "chrom": "1",
                     "pos": 25_600_000 + i * 100,
                     "genotype": "AA" if i < rhd_calls else "--"})
    df = pd.DataFrame(rows).set_index("rsid")
    return df


def test_abo_type_a_positive_from_A_over_O() -> None:
    df = _df({"rs8176719": "DI", "rs8176746": "GG", "rs8176747": "CC",
              "rs7853989": "GG", "rs601338": "AA"})
    r = bt.analyze_blood_type(df)
    assert r["abo"]["phenotype"] == "A"
    assert r["abo"]["genotype"] == "A/O"
    assert r["abo"]["carries_hidden_O"] is True
    assert r["rhd"]["status"] == "Rh-positive (Rh+)"
    assert r["combined"] == "A+"
    assert r["secretor"]["secretor_status"].startswith("Non-secretor")


def test_abo_type_o_from_homozygous_deletion() -> None:
    df = _df({"rs8176719": "DD", "rs8176746": "GG", "rs8176747": "CC"})
    r = bt.analyze_blood_type(df)
    assert r["abo"]["phenotype"] == "O"
    assert r["abo"]["genotype"] == "O/O"


def test_abo_type_ab() -> None:
    df = _df({"rs8176719": "II", "rs8176746": "GT", "rs8176747": "CG"})
    r = bt.analyze_blood_type(df)
    assert r["abo"]["phenotype"] == "AB"
    assert r["abo"]["genotype"] == "A/B"


def test_rhd_negative_from_low_coverage() -> None:
    # Simulate RHD-gene-deleted individual: mostly no-calls at RHD locus.
    df = _df({"rs8176719": "II", "rs8176746": "GG", "rs8176747": "CC"},
             rhd_calls=3, rhd_total=42)
    r = bt.analyze_blood_type(df)
    assert r["rhd"]["status"] == "Rh-negative (Rh-)"


def test_secretor_full_secretor() -> None:
    df = _df({"rs601338": "GG"})
    s = bt.analyze_secretor_bombay(df)
    assert s["available"] and s["secretor_status"] == "Secretor"


def test_empty_gracefully() -> None:
    df = pd.DataFrame(columns=["chrom", "pos", "genotype"])
    df.index.name = "rsid"
    r = bt.analyze_blood_type(df)
    # Should not crash; both sub-analyses return empty/indeterminate.
    assert "abo" in r and "rhd" in r
