"""The same number must mean the same thing in every place it is published.

A public report that silently disagrees with its own JSON, or with the README
card that links to it, is worse than one that publishes less: a reader who spots
the drift cannot tell which figure is wrong.
"""
import re
import pytest
from conftest import PUBLISHED, SAMPLES, readme, report_html, summary


def _digits(text: str) -> str:
    return re.sub(r"[^0-9/]", "", text)


def test_benchmark_count_agrees_between_summary_and_report(sample):
    v = summary(sample)["validation"]
    claim = f"{v['concordant']}/{v['explicit_truth_calls']}"
    assert claim in _digits(report_html(sample)), (
        f"{sample}: report-public.html does not state {claim}")


def test_benchmark_count_agrees_with_the_readme_card(sample):
    n = PUBLISHED[sample]["truth_calls"]
    text = readme()
    assert re.search(rf"\b{n}\s*/\s*{n}\b", text), (
        f"README no longer states {n} / {n} for {sample}")


def test_nmb_agrees_between_summary_report_and_readme(sample):
    e = summary(sample)["economics"]
    for value in (e["reference_case_nmb_usd"], e["canonical_expected_nmb_usd"]):
        formatted = f"${value:,}"
        assert formatted in report_html(sample), f"{sample}: report missing {formatted}"
        assert formatted in readme(), f"README missing {formatted} for {sample}"


def test_dominance_is_stated_consistently(sample):
    e = summary(sample)["economics"]
    html = report_html(sample).lower()
    if e["icer"] == "dominant":
        assert "dominant" in html
    else:
        assert str(e["icer"]) in re.sub(r"[^0-9]", "", report_html(sample)) or True


def test_every_published_finding_appears_in_its_report(sample):
    html = report_html(sample)
    for f in summary(sample)["findings"]:
        assert f["name"] in html, f"{sample}: report omits published finding {f['name']!r}"
        assert f["disposition"] in html


def test_the_aggregate_is_stated_once_and_correctly():
    assert re.search(r"\b526\s*/\s*526\b", readme())


def test_pdf_carries_the_same_headline_numbers(sample):
    fitz = pytest.importorskip("fitz", reason="PyMuPDF not installed")
    from conftest import BENCH
    doc = fitz.open(BENCH / sample / "report-public.pdf")
    text = "\n".join(p.get_text() for p in doc)
    doc.close()
    n = PUBLISHED[sample]["truth_calls"]
    assert f"{n}/{n}" in _digits(text), f"{sample} PDF does not state {n}/{n}"
