"""Two different quantities must stay two different quantities.

The central public claim of this project is that GenomeLens does not reduce a
genome to one number. The reference-case result and the canonical standardized
pathway sum answer different questions, and the schema has to keep them apart
rather than inviting a reader to add them.
"""
from conftest import SAMPLES, summary


def test_reference_case_and_canonical_are_separate_fields(s):
    e = s["economics"]
    assert "reference_case_nmb_usd" in e and "canonical_expected_nmb_usd" in e


def test_no_single_collapsed_genome_value_is_published(s):
    """There must be no field a reader could mistake for 'the value of this
    genome'. Constructing one requires adding quantities with different
    estimands, beneficiaries and evidence grades."""
    keys = {k.lower() for k in s["economics"]}
    for banned in ("total", "genome_value", "total_value", "sum", "combined",
                   "overall", "aggregate_value"):
        assert not any(banned in k for k in keys), (
            f"economics exposes a collapsed-total-looking field matching "
            f"{banned!r}: {sorted(keys)}")


def test_the_two_quantities_actually_differ(s):
    """If they were equal the separation would be decorative."""
    e = s["economics"]
    assert e["reference_case_nmb_usd"] != e["canonical_expected_nmb_usd"]


def test_the_public_report_states_they_are_not_summed(sample):
    from conftest import report_html
    html = report_html(sample).lower()
    assert "never summed" in html or "not summed" in html, (
        "the public report no longer tells the reader these are different "
        "quantities; the numbers alone invite addition")


def test_each_finding_keeps_its_own_disposition(s):
    """A finding's disposition is per finding, not a single verdict for the
    sample."""
    dispositions = [f["disposition"] for f in s["findings"]]
    assert dispositions, "no findings published"
    assert all(isinstance(d, str) and d for d in dispositions)
