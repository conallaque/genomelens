"""Published economic figures, and the identities they must satisfy.

These check the public numbers against each other and against the definitions
in docs/ECONOMICS.md. They do not reconstruct how any figure was produced.
"""
import pytest
from conftest import PUBLISHED, SAMPLES, summary


def test_published_headline_figures(s, sample):
    e = s["economics"]
    assert e["reference_case_nmb_usd"] == PUBLISHED[sample]["nmb"]
    assert e["canonical_expected_nmb_usd"] == PUBLISHED[sample]["canonical"]


def test_dominance_is_reported_as_a_status_never_as_a_ratio(s):
    """docs/ECONOMICS.md: a dominant strategy's ratio is finite but
    uninterpretable, so it is reported as a status. A number here would invite
    a reader to compare it against a threshold, which is meaningless."""
    e = s["economics"]
    dominant = e["incremental_cost_usd"] < 0 and e["incremental_qaly"] > 0
    if dominant:
        assert e["icer"] == "dominant", (
            "less costly and more effective, but a ratio was published")
    else:
        assert e["icer"] != "dominant"


def test_a_published_icer_is_consistent_with_its_own_components(s):
    """ICER = incremental cost / incremental QALYs. Published to the nearest
    dollar, so the tolerance is the rounding of its own inputs, not a fudge."""
    e = s["economics"]
    if not isinstance(e["icer"], (int, float)):
        pytest.skip("dominant or not applicable — no ratio to check")
    expected = e["incremental_cost_usd"] / e["incremental_qaly"]
    assert abs(e["icer"] - expected) <= abs(expected) * 0.01


def test_no_published_case_claims_dominance():
    """HG003 WAS the dominant case. It is not any more, and that is correct.

    Dominance means cheaper AND more effective. HG003's incremental cost was
    -$115 when whole-genome sequencing was priced at the vendor's promotional
    $300. At the $599 list price the same strategy costs $388.53 more than
    the comparator, so it buys health at a price rather than for free.

    The test is inverted rather than deleted. A dominance claim is the single
    most quotable thing this repository publishes and the one where a silent
    sign flip would be least visible, so the guard is kept pointing the other
    way: no sample may claim dominance while costing more.
    """
    e = summary("HG003")["economics"]
    assert e["incremental_cost_usd"] > 0, "HG003 no longer costs less"
    assert e["icer"] != "dominant"
    for x in SAMPLES:
        assert summary(x)["economics"]["icer"] != "dominant", (
            f"{x} claims dominance; no current published case is dominant")


def test_positive_incremental_cost_never_claims_dominance():
    for x in SAMPLES:
        e = summary(x)["economics"]
        if e["incremental_cost_usd"] > 0:
            assert e["icer"] != "dominant", f"{x} claims dominance while costing more"


def test_economic_figures_are_present_for_every_sample(s):
    e = s["economics"]
    for k in ("reference_case_nmb_usd", "canonical_expected_nmb_usd",
              "incremental_cost_usd", "incremental_qaly", "icer"):
        assert k in e, f"missing published economic field: {k}"
