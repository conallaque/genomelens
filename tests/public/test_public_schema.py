"""The public schema is closed: known-safe fields only.

DESIGN NOTE. The obvious way to write this is a list of prohibited internal
field names. That would be a mistake -- a public deny-list of production field
names is itself a map of the private schema, and publishing it hands a reader
the internal vocabulary for free. It also fails open: a field nobody thought to
ban passes.

So these tests assert the complement. Every key at every level must appear in an
explicit allowlist, and anything else fails regardless of what it is called. A
new engine field cannot reach the public artifact by being unanticipated; it
would have to be added here deliberately.
"""
import json
from conftest import BENCH, SAMPLES, summary

TOP = {"sample", "dataset", "validation", "findings", "economics", "limitations"}
DATASET = {"source", "reference_build", "benchmark_release", "benchmark_scope"}
VALIDATION = {"explicit_truth_calls", "concordant", "mismatches", "normalization_failures"}
FINDING = {"name", "category", "disposition", "evidence_basis"}
ECONOMICS = {"reference_case_nmb_usd", "canonical_expected_nmb_usd",
             "incremental_cost_usd", "incremental_qaly", "icer"}

CATEGORIES = {"carrier", "pharmacogenomics", "polygenic risk", "longevity",
              "wellness", "hereditary condition", "other"}
DISPOSITIONS = {"supported", "conditional", "unresolved", "not assessed"}


def test_top_level_keys_are_exactly_the_allowlist(s):
    assert set(s) == TOP, f"unexpected top-level keys: {set(s) ^ TOP}"


def test_nested_blocks_are_closed(s):
    assert set(s["dataset"]) == DATASET
    assert set(s["validation"]) == VALIDATION
    assert set(s["economics"]) == ECONOMICS


def test_every_finding_is_closed(s):
    for f in s["findings"]:
        assert set(f) == FINDING, f"unexpected finding keys: {set(f) ^ FINDING}"


def test_findings_use_public_vocabulary_only(s):
    for f in s["findings"]:
        assert f["category"] in CATEGORIES, f"non-public category {f['category']!r}"
        assert f["disposition"] in DISPOSITIONS, f"non-public disposition {f['disposition']!r}"


def test_no_nested_structure_smuggles_extra_data(s):
    """A finding must be a flat record of four strings. Nesting is where a
    payload fragment would hide."""
    for f in s["findings"]:
        for k, v in f.items():
            assert isinstance(v, str), f"{k} is {type(v).__name__}, expected str"


def test_the_summary_stays_small(sample):
    """Size is a proxy for scope creep. A public summary that grows into a
    machine-readable shadow of the engine has stopped being a summary."""
    n = (BENCH / sample / "summary.json").stat().st_size
    assert n < 4096, f"summary.json is {n} bytes; it is meant to be a presentation artifact"


def test_findings_are_a_representative_subset_not_a_dump(s):
    assert 1 <= len(s["findings"]) <= 5, (
        f"{len(s['findings'])} findings published; the summary is a selection, "
        "the full sanitized report is the complete output")


def test_limitations_are_published_with_every_sample(s):
    assert len(s["limitations"]) >= 3
    joined = " ".join(s["limitations"]).lower()
    assert "not clinical validation" in joined
    assert "not realized savings" in joined
