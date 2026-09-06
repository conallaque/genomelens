"""Catalog rendering: the per-variant "Recommendation & context" block.

Locks the depth pass that surfaces `recommendation`, `chip_coverage_note`,
and `cross_references` (all curated in data/snp_database.json) in the variant
tables. The synthetic test genome does not match any coverage-note or
cross-ref variant, so these paths need crafted records to be exercised.
"""

from report.renderers import build_html_report


def _record(**overrides):
    base = {
        "rsid": "rs0000001",
        "gene": "GENE1",
        "variant_name": "p.Test",
        "category": "Methylation & Folate",
        "my_genotype": "AG",
        "risk_allele": "A",
        "risk_copies": 1,
        "significance": "moderate",
        "summary": "A curated one-line summary.",
        "recommendation": "Do the recommended thing.",
        "cross_references": [],
        "chip_coverage_note": None,
    }
    base.update(overrides)
    return base


def _render(records):
    return build_html_report(
        tier1_results=records,
        apoe_genotype=None,
        ai_results={},
        exec_summary=None,
        dna_filepath="synthetic.txt",
        no_ai=True,
        model="none",
    )


def test_recommendation_block_renders():
    html = _render([_record(recommendation="Take 400mcg methylfolate daily.")])
    assert "var-detail" in html
    assert "What to do" in html
    assert "Take 400mcg methylfolate daily." in html
    # existing summary column is untouched
    assert "A curated one-line summary." in html


def test_chip_coverage_note_renders():
    note = "Consumer chips often omit this position; a negative cannot rule it out."
    html = _render([_record(chip_coverage_note=note)])
    assert "rec-cov" in html
    assert "Chip coverage" in html
    assert note in html


def test_cross_references_render():
    html = _render(
        [
            _record(
                cross_references=[
                    {
                        "category": "Male Fertility & Sperm",
                        "implication": "Impairs folate-dependent sperm DNA integrity.",
                    }
                ]
            )
        ]
    )
    assert "rec-xref" in html
    assert "Male Fertility &amp; Sperm" in html
    assert "Impairs folate-dependent sperm DNA integrity." in html


def test_no_context_block_when_all_fields_empty():
    html = _render([_record(recommendation="", cross_references=[], chip_coverage_note=None)])
    # No curated context to show -> no expandable block on that row
    assert "Recommendation &amp; context" not in html


def test_imputed_provenance_threads_and_renders():
    """--impute adds source/r2 columns to snps_df; tier1_lookup must record
    them and the catalog must flag imputed calls (never mistaking a
    statistical call for a measured one). Chip data is unaffected."""
    import pandas as pd

    import analyze
    db = {"rs1801133": {"gene": "MTHFR", "name": "C677T",
                        "category": "Methylation & Folate", "risk_allele": "T",
                        "normal_allele": "C", "significance": "high",
                        "summary": "s", "recommendation": "r"}}
    imp = pd.DataFrame({"genotype": ["CT"], "source": ["imputed"], "r2": [0.82]},
                       index=["rs1801133"])
    rec = analyze.tier1_lookup(imp, db)[0][0]
    assert rec["source"] == "imputed" and rec["r2"] == 0.82
    html = _render([rec])
    assert '<span class="imp-badge"' in html and "imputed r²=0.82" in html
    # chip data: no source column -> defaults to chip, no badge, no regression
    chip = pd.DataFrame({"genotype": ["CT"]}, index=["rs1801133"])
    rec2 = analyze.tier1_lookup(chip, db)[0][0]
    assert rec2["source"] == "chip" and rec2["r2"] is None
    assert '<span class="imp-badge"' not in _render([rec2])
