"""ClinPGx/PharmGKB clinical-variant annotation module."""

import pandas as pd
import pharmgkb_clinical as pg


def test_no_data_returns_unavailable():
    assert pg.analyze_pharmgkb_clinical(None)["available"] is False


def test_only_clean_rsids_parsed():
    # star alleles / HGVS variants must be excluded (handled by pgx.py)
    import re
    for row in pg._load_table():
        assert re.fullmatch(r"rs\d+", row["rsid"]), row["rsid"]
        assert row["gene"], "empty-gene rows must be dropped"


def test_typed_variant_is_reported_and_tiered():
    tbl = pg._load_table()
    assert tbl, "clinicalVariants.tsv should be present"
    rsids = list(dict.fromkeys(t["rsid"] for t in tbl))[:200]
    df = pd.DataFrame({"genotype": ["AG"] * len(rsids)}, index=rsids)
    res = pg.analyze_pharmgkb_clinical(df)
    assert res["available"]
    assert res["n_typed_variants"] == res["n_high"] + res["n_low"]
    # every high entry's best level is a strong tier; low entries are 3/4
    for e in res["high"]:
        assert e["best_level"] in ("1A", "1B", "2A", "2B")
    for e in res["low"]:
        assert e["best_level"] in ("3", "4")
    # each entry carries a genotype and at least one annotation
    for e in res["high"] + res["low"]:
        assert e["genotype"] and e["annotations"]


def test_untyped_rsids_excluded():
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs_not_in_table_999999"])
    res = pg.analyze_pharmgkb_clinical(df)
    assert res["available"]
    assert res["n_typed_variants"] == 0
