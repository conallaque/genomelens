"""The external benchmark claim, checked against the published artifacts.

The claim this repository makes is narrow on purpose: agreement with the GIAB
truth set on positions where that truth set states a genotype, inside defined
confident regions. These tests hold the published numbers to that claim and no
wider.
"""
from conftest import AGGREGATE, PUBLISHED, SAMPLES, summary


def test_each_sample_is_fully_concordant(s, sample):
    v = s["validation"]
    assert v["concordant"] == v["explicit_truth_calls"] == PUBLISHED[sample]["truth_calls"]


def test_no_mismatches_or_normalization_failures(s):
    assert s["validation"]["mismatches"] == 0
    assert s["validation"]["normalization_failures"] == 0


def test_the_per_sample_counts_sum_to_the_published_aggregate():
    total = sum(summary(x)["validation"]["explicit_truth_calls"] for x in SAMPLES)
    assert total == AGGREGATE, (
        f"per-sample counts sum to {total}; README publishes {AGGREGATE}")


def test_benchmark_provenance_is_pinned(s):
    """An unversioned benchmark claim is ambiguous: GIAB has published releases
    after v4.2.1, so the release string is load-bearing, not decoration."""
    d = s["dataset"]
    assert d["source"] == "Genome in a Bottle"
    assert d["reference_build"] == "GRCh38"
    assert d["benchmark_release"] == "v4.2.1"
    assert d["benchmark_scope"] == "chr1-22"


def test_the_benchmark_scope_excludes_sex_chromosomes_and_mitochondria(s):
    """chr1-22 means chrX, chrY and chrM carry zero benchmark coverage. That is
    a different statement from the modality limits in LIMITATIONS.md, and
    conflating them would overstate what was tested."""
    assert s["dataset"]["benchmark_scope"] == "chr1-22"


def test_sample_identifiers_are_distinct_and_self_consistent():
    ids = [summary(x)["sample"] for x in SAMPLES]
    assert ids == list(SAMPLES), f"sample identifiers mis-propagated: {ids}"
    assert len(set(ids)) == 3


def test_the_trio_is_not_three_independent_genomes():
    """HG002 is the son of HG003 and HG004. Roughly half of HG002's alleles are
    present in a parent by descent, so 526 comparisons are not 526 independent
    observations. Asserted so the aggregate is never described as breadth."""
    from conftest import readme
    text = readme().lower()
    assert "trio" in text, (
        "README no longer identifies the samples as a trio; the pedigree "
        "relationship must stay visible next to the aggregate figure")
