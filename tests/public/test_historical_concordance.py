"""The historical Genome in a Bottle concordance record, checked across the places it is published.

The claim this repository keeps is narrow on purpose: agreement with the GIAB truth set on positions where that truth set states a genotype,
inside defined confident regions, for an EARLIER engine build. It is not a v1.11 result and these tests hold the published numbers to that claim
and no wider. (The earlier build's per-genome economic reports are no longer published; the record that remains is this concordance.)
"""
import re

import pytest
from conftest import AGGREGATE, PUBLISHED, ROOT, SAMPLES, readme, validation


def _section() -> str:
    t = readme()
    i = t.index("## Earlier GIAB trio concordance")
    return t[i:t.index("## Public verification", i)]


@pytest.mark.parametrize("sample", SAMPLES)
def test_each_sample_is_fully_concordant_in_the_readme(sample):
    n = PUBLISHED[sample]["truth_calls"]
    row = re.search(rf"\|\s*{sample}\s*\|[^\n]*", _section())
    assert row, f"the README's concordance table has no row for {sample}"
    assert re.search(rf"\b{n}/{n}\b", row.group(0)) and "0 mismatches" in row.group(0), row.group(0)


def test_the_per_sample_counts_sum_to_the_published_aggregate():
    assert sum(v["truth_calls"] for v in PUBLISHED.values()) == AGGREGATE
    assert re.search(rf"\b{AGGREGATE}\s*/\s*{AGGREGATE}\b", _section()), f"README no longer states {AGGREGATE} / {AGGREGATE}"


def test_no_mismatches_or_normalization_failures_are_stated():
    sec = _section()
    assert re.search(r"\b0\s+mismatches", sec) and re.search(r"\b0\s+normalization failures", sec)


def test_benchmark_provenance_is_pinned_in_the_readme_and_in_the_validation_document():
    """An unversioned benchmark claim is ambiguous: GIAB has published releases after v4.2.1, so the release string is load-bearing."""
    for text in (_section(), validation()):
        for token in ("v4.2.1", "GRCh38", "chr1"):
            assert token in text, token


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_validation_document_states_each_per_sample_count(sample):
    v = validation()
    block = v[v.index(f"### {sample}"):][:1500]
    assert str(PUBLISHED[sample]["truth_calls"]) in block


def test_the_trio_is_not_three_independent_genomes():
    """HG002 is the son of HG003 and HG004, so the aggregate is never to be described as breadth."""
    assert "not statistically independent" in _section()
    assert "not 526 independent observations" in validation()


def test_the_record_is_labeled_historical_and_not_v111_validation():
    sec = _section().lower()
    assert "historical" in sec and "not representative of v1.11 validation" in sec


def test_no_per_genome_report_of_the_earlier_build_is_left_in_the_tree():
    assert not (ROOT / "artifacts").exists(), "the earlier trio's per-genome reports were removed from the tree"
    for md in sorted(ROOT.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        for gone in ("artifacts/public-benchmarks", "report-public.pdf", "report-public.html"):
            assert gone not in text, f"{md.relative_to(ROOT)} still points at a removed report ({gone})"
