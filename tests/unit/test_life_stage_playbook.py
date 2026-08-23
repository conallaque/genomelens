"""Unit tests for life_stage_playbook.py."""

from __future__ import annotations

import life_stage_playbook as lsp

HOLISTIC = {"genome_leverage": {"tier": "Favorable", "score": 74},
            "insights": [{"id": "glucose_hba1c_stress_discordance"}]}
ADDICTION = {"composite": {"clinical_flags": [
    {"title": "🚭 Never start smoking", "text": "chrna5"}]}}
NEURO = {"composite": {"bdnf_class": "Val/Val (full)", "comt_class": "middle"}}
IMMUNO = {"findings": [{"impact": "susceptible",
                        "name": "Severe influenza risk", "gene": "IFITM3"}]}
FAM = {"n_recessive": 2, "n_dominant": 1}


def _run(age):
    return lsp.analyze_life_stage_playbook(
        age=age, holistic_synthesis_result=HOLISTIC, immunogenetics_result=IMMUNO,
        addiction_genetics_result=ADDICTION, neurochemistry_result=NEURO,
        family_planning_result=FAM, tier1_results={"apoe_genotype": "e3/e3"})


def test_age_24_highlights_20s() -> None:
    r = _run(24)
    assert r["current_decade"] == "20s"
    twenties = next(d for d in r["decades"] if d["key"] == "20s")
    assert twenties["is_current"] is True


def test_age_47_highlights_40s() -> None:
    r = _run(47)
    assert r["current_decade"] == "40s"
    assert all(d["is_current"] == (d["key"] == "40s") for d in r["decades"])


def test_no_age_highlights_nothing() -> None:
    r = _run(None)
    assert r["current_decade"] is None
    assert r["age_known"] is False
    assert all(d["is_current"] is False for d in r["decades"])
    assert "not provided" in r["note"].lower()


def test_chrna5_flag_appears_in_20s() -> None:
    r = _run(24)
    twenties = next(d for d in r["decades"] if d["key"] == "20s")
    texts = " ".join(gi["text"].lower() for gi in twenties["genome_items"])
    assert "smoking" in texts


def test_leverage_framing_present_when_favorable() -> None:
    r = _run(24)
    twenties = next(d for d in r["decades"] if d["key"] == "20s")
    sources = {gi["source"] for gi in twenties["genome_items"]}
    assert "holistic_synthesis" in sources


def test_all_five_decades_present() -> None:
    r = _run(30)
    keys = [d["key"] for d in r["decades"]]
    assert keys == ["20s", "30s", "40s", "50s", "60s+"]


def test_resolve_age_prefers_explicit() -> None:
    assert lsp.resolve_age(33, {"clinical": {"age_used": 41.0}}) == 33


def test_resolve_age_falls_back_to_bloodwork_clinical_age_used() -> None:
    # Locks the true nested key path (top-level has no age); regressed twice.
    assert lsp.resolve_age(None, {"clinical": {"age_used": 41.0}}) == 41
    assert isinstance(lsp.resolve_age(None, {"clinical": {"age_used": 41.0}}), int)


def test_resolve_age_none_when_absent() -> None:
    assert lsp.resolve_age(None, None) is None
    assert lsp.resolve_age(None, {"clinical": {}}) is None
    assert lsp.resolve_age(None, {}) is None


def test_genome_items_are_source_tagged() -> None:
    r = _run(24)
    for d in r["decades"]:
        for gi in d["genome_items"]:
            assert gi.get("source"), "every genome-driven item must cite a source"
