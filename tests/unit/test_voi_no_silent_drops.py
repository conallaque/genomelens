"""Regression: a finding must never leave the economic model without a trace.

``_classify_category`` returns ('', '') for any category string it does not
recognise, and ``_collect`` then skips that finding. Because health_economics.py
gains new ``source=`` labels regularly, an unmapped label silently removed
findings from the model while the report still claimed a complete computation.

The contract asserted here: every finding is either valued, or recorded — and a
category dropped by oversight (rather than by documented decision) makes
``fully_computed`` False.
"""
from __future__ import annotations

import re

from econ import health_economics as he
from econ import value_of_information as voi


def _he_source() -> str:
    """Read the module's own source, located from the module rather than from
    the working directory — a hardcoded ``health_economics.py`` broke the moment
    the file moved into the econ package."""
    return open(he.__file__, encoding="utf-8").read()


def _econ(*cats):
    return {"findings_with_economics": [
        {"finding": f"synthetic {c}", "category": c, "confidence": "moderate",
         "outcome_value": 3000, "qaly_gain": 0.2} for c in cats]}


def test_unmapped_category_is_recorded_not_dropped():
    r = voi.analyze_value_of_information(
        _econ("Pharmacogenomics", "Totally Made Up Source"),
        None, None, n_mc=200, seed=1)
    assert r["n_findings"] == 1
    assert r["n_unvalued"] == 1
    assert r["unvalued_findings"][0]["category"] == "Totally Made Up Source"


def test_oversight_makes_fully_computed_false():
    r = voi.analyze_value_of_information(
        _econ("Pharmacogenomics", "Totally Made Up Source"),
        None, None, n_mc=200, seed=1)
    assert r["fully_computed"] is False
    assert "Totally Made Up Source" in r["unmapped_categories"]


def test_documented_exclusion_does_not_count_as_incomplete():
    # A category listed in NOT_VALUED is a modelling decision, not a gap.
    r = voi.analyze_value_of_information(
        _econ("Pharmacogenomics", "Family Planning"),
        None, None, n_mc=200, seed=1)
    assert r["n_unvalued"] == 1
    assert r["n_unvalued_intentional"] == 1
    assert r["unmapped_categories"] == []
    assert r["fully_computed"] is True


def test_every_module_source_reaches_the_model():
    # THE POINT OF WIRING THE MODULES IN: each source label must either produce
    # a valued finding or be the one documented ethical exclusion.
    import re
    src = _he_source()
    cats = sorted(set(re.findall(r'source="([^"]+)"', src)))
    valued = [c for c in cats if voi._classify_category(c, "")[0] in ("pgx", "coi")]
    excluded = [c for c in cats if voi._not_valued_reason(c)]
    assert len(valued) + len(excluded) == len(cats)
    # Two categories are excluded on principle, each for a documented reason:
    # Family Planning because monetising a reproductive outcome prices a
    # prospective child, and Longevity because the composite re-aggregates
    # variants already valued individually (double counting) at a rate that
    # had no published source. Anything else appearing here is an oversight.
    assert excluded == ["Family Planning", "Longevity"], (
        "Family Planning and Longevity are the only categories excluded on "
        f"principle; got {excluded}")
    for cat in excluded:
        assert len(voi._not_valued_reason(cat)) > 80, (
            f"{cat} must carry a substantive documented reason, not a stub")


def test_weak_evidence_sources_are_discounted_not_dropped():
    # A hypothesis-generating panel must still enter the model, at a discount.
    weak = voi._evidence_haircut("PheWAS Biomarker")
    strong = voi._evidence_haircut("Clinical Variant (ClinVar)")
    assert 0.0 < weak < strong == 1.0


def test_haircut_ordering_tracks_evidence_strength():
    h = voi._evidence_haircut
    assert h("PheWAS Biomarker") < h("Polygenic Risk") < h("Pharmacogenomics")
    assert h("Wellness Prediction") < h("Addiction Genetics")


def test_weak_source_yields_less_value_than_strong_source():
    # End-to-end: same finding, different source label -> different NMB.
    strong = voi.analyze_value_of_information(
        _econ("Clinical Variant (ClinVar)"), None, None, n_mc=400, seed=2)
    weak = voi.analyze_value_of_information(
        _econ("PheWAS Biomarker"), None, None, n_mc=400, seed=2)
    assert weak["nmb_rows"][0]["nmb"] < strong["nmb_rows"][0]["nmb"]


def test_every_coi_anchor_used_by_the_router_exists():
    # A router pointing at a missing COI key would silently zero the finding.
    import re
    src = _he_source()
    for c in sorted(set(re.findall(r'source="([^"]+)"', src))):
        kind, key = voi._classify_category(c, "")
        if kind == "coi":
            assert key in voi.COI, f"{c!r} routes to unknown COI key {key!r}"


def test_new_coi_anchors_are_sourced():
    for key in ("SubstanceUse", "Depression", "Autoimmune", "Urologic",
                "IronOverload"):
        assert key in voi.COI
        assert voi.COI[key]["src"].strip(), f"{key} has no citation"
        assert voi.COI[key]["cost"] > 0


def test_every_intentional_exclusion_carries_a_reason():
    for cat, reason in voi.NOT_VALUED.items():
        assert reason.strip(), f"{cat} has no documented reason"
        assert len(reason) > 40, f"{cat} reason is too thin to be meaningful"


def test_every_emitted_category_is_mapped_or_documented():
    # THE INTEGRATION GUARD: health_economics.py and value_of_information.py
    # must not drift. Any new source= label has to be either given an economic
    # mapping or explicitly registered as not-valued.
    src = _he_source()
    for cat in sorted(set(re.findall(r'source="([^"]+)"', src))):
        kind, _ = voi._classify_category(cat, "")
        mapped = kind in ("pgx", "coi")
        documented = bool(voi._not_valued_reason(cat))
        assert mapped or documented, (
            f"category {cat!r} is emitted by health_economics.py but is neither "
            f"mapped in _classify_category nor registered in NOT_VALUED")


def test_clean_run_reports_fully_computed():
    r = voi.analyze_value_of_information(_econ("Pharmacogenomics"),
                                         None, None, n_mc=200, seed=1)
    assert r["n_unvalued"] == 0
    assert r["fully_computed"] is True
