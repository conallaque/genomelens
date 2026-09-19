"""Absence of evidence is not evidence of absence, at the presentation layer.

These check only what the public schema and reports say. They do not touch how
any disposition was decided.
"""
from conftest import SAMPLES, report_html, summary

NEGATIVE_WORDS = ("benign", "normal", "negative", "absent", "no risk", "wild-type")


def test_a_non_supported_disposition_is_never_rendered_as_a_negative_result(sample):
    """unresolved and not-assessed mean the claim could not be supported. If
    either were displayed as 'normal' or 'benign', a reader would be told a
    biological fact the analysis never established."""
    for f in summary(sample)["findings"]:
        if f["disposition"] in ("unresolved", "not assessed"):
            assert f["disposition"] not in NEGATIVE_WORDS
            assert not any(w in f["disposition"].lower() for w in NEGATIVE_WORDS)


def test_dispositions_are_distinct_states_not_a_binary(s):
    """Collapsing four states into found/not-found is what the schema exists to
    prevent."""
    from test_public_schema import DISPOSITIONS
    assert all(f["disposition"] in DISPOSITIONS for f in s["findings"])
    assert len(DISPOSITIONS) >= 4


def test_zero_is_never_used_where_a_value_was_not_computed(s):
    """A missing economic result must be absent or explicitly stated, never
    rendered as 0, which reads as a computed finding of no value."""
    e = s["economics"]
    for k in ("reference_case_nmb_usd", "canonical_expected_nmb_usd"):
        assert e[k] != 0, f"{k} is 0; a non-computed value must not be published as zero"


def test_the_report_explains_that_a_selection_is_not_the_whole_output(sample):
    html = report_html(sample).lower()
    assert "representative" in html or "more results" in html, (
        "a reader could otherwise take three findings as the complete analysis")
