"""Unit tests for holistic_synthesis.py."""

from __future__ import annotations

from wellness import holistic_synthesis as hs


def test_glucose_hba1c_discordance_detected() -> None:
    bw = {"clinical": {"systems": [{"markers": [
        {"name": "Fasting Glucose", "value": 105},
        {"name": "HbA1c", "value": 5.4},
    ]}], "flags": []}}
    r = hs.analyze_holistic_synthesis(bloodwork_result=bw)
    ids = {i["id"] for i in r["insights"]}
    assert "glucose_hba1c_stress_discordance" in ids


def test_fut2_crp_insight() -> None:
    immuno = {"findings": [
        {"gene": "FUT2", "phenotype": "non-secretor", "verdict": "resistance",
         "impact": "protective"}]}
    bw = {"clinical": {"systems": [{"markers": [
        {"name": "hs-CRP", "value": 2.0},
    ]}], "flags": []}}
    r = hs.analyze_holistic_synthesis(immunogenetics_result=immuno,
                                      bloodwork_result=bw)
    ids = {i["id"] for i in r["insights"]}
    assert "fut2_crp_baseline" in ids


def test_apoe_e4_lipid_amplification_triggers() -> None:
    tier1 = {"apoe_genotype": "e3/e4"}
    bw = {"clinical": {"systems": [{"markers": [
        {"name": "LDL Cholesterol", "value": 140},
    ]}], "flags": []}}
    r = hs.analyze_holistic_synthesis(tier1_results=tier1, bloodwork_result=bw)
    ids = {i["id"] for i in r["insights"]}
    assert "apoe_e4_lipid_amplification" in ids


def test_genome_leverage_favorable_when_lots_protective() -> None:
    r = hs.analyze_holistic_synthesis(
        tier1_results={"apoe_genotype": "e3/e3"},
        immunogenetics_result={"headlines": [{"name": "FUT2"}, {"name": "IL28B"},
                                              {"name": "PRNP"}]},
        neurochemistry_result={"composite": {
            "comt_class": "middle", "maoa_class": "MAOA-H",
            "bdnf_class": "Val/Val (full)"}},
    )
    assert r["genome_leverage"]["score"] >= 60
    assert r["genome_leverage"]["tier"] in ("Favorable", "Very favorable")


def test_chrna5_realised_prevention() -> None:
    r = hs.analyze_holistic_synthesis(
        neurochemistry_result={"findings": [
            {"gene": "CHRNA5", "phenotype": "A-carrier"}], "composite": {}},
        meta={"smoker": False},
    )
    ids = {i["id"] for i in r["insights"]}
    assert "chrna5_realised_prevention_success" in ids


def test_empty_input_gracefully() -> None:
    r = hs.analyze_holistic_synthesis()
    # Should still return a leverage score (neutral)
    assert r["genome_leverage"]["score"] == 50
