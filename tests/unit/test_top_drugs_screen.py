"""Top-prescribed-drugs pharmacogenomic screen."""

import pandas as pd

import top_drugs_screen as t


def test_reference_list_present_and_sane():
    ref = t._load_reference()
    assert len(ref) >= 500, "expected the full most-prescribed reference menu"
    for d in ref[:20]:
        assert d["generic"] and d["class"]


def test_pgx_relevance_comes_from_real_data_not_the_menu():
    # A CYP2D6 poor-metabolizer must push CYP2D6 drugs (codeine, atomoxetine)
    # into the actionable tier; drugs with no PGx entry land in no_pgx.
    pgx = {"per_gene": {"CYP2D6": {"phenotype": "Poor Metabolizer",
                                   "phenotype_code": "PM"}}}
    df = pd.DataFrame({"genotype": ["AG"]}, index=["rs3892097"])
    res = t.analyze_top_drugs(df, pgx)
    assert res["available"]
    actionable_names = {e["generic"] for e in res["actionable"]}
    assert "codeine" in actionable_names
    # every actionable entry has a non-normal metabolizer phenotype recorded
    for e in res["actionable"]:
        codes = {p["code"] for p in e["gene_phenotypes"]}
        assert codes & {"PM", "IM", "UM", "RM", "POS"}
    # tiers partition the whole list
    assert (res["n_actionable"] + res["n_typed_normal"] + res["n_pgx_relevant"]
            + res["n_no_pgx"]) == res["n_screened"]


def test_no_pgx_drugs_are_labeled_not_fabricated():
    res = t.analyze_top_drugs(pd.DataFrame({"genotype": []}), {})
    # drugs absent from CPIC + PharmGKB must be reported as no-PGx, never as
    # a finding with an invented gene link
    for e in res["no_pgx"]:
        assert not e["genes"] and not e["cpic_level"] and not e["clin_level"]
