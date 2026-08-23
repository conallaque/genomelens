"""Tests for the AI-upgrade layer: shared module-context digest, per-module
AI specs, and the section AI-callout injector."""

from __future__ import annotations

import analyze
import renderers


def _sample_results():
    return dict(
        holistic_synthesis_result={"available": True,
            "genome_leverage": {"score": 74, "tier": "Favorable", "narrative": "env dominates"},
            "insights": [{"title": "Ancestral diet fit", "explanation": "Mediterranean+dairy"}]},
        immunogenetics_result={"available": True, "n_protective": 7, "n_susceptible": 3,
            "headlines": [{"name": "FUT2 non-secretor", "gene": "FUT2", "genotype": "AA",
                           "verdict": "norovirus resistance"}],
            "findings": [{"impact": "susceptible", "name": "IFITM3", "verdict": "severe flu"}]},
        neurochemistry_result={"composite": {"comt_class": "middle", "maoa_class": "MAOA-H",
            "bdnf_class": "Val/Val (full)", "stress_response_profile": "adaptive",
            "stimulant_response": "intermediate", "caffeine_protocol": "espresso am"}},
        addiction_genetics_result={"composite": {"alcohol_tier": "Modestly susceptible",
            "overall_tier": "Mixed", "clinical_flags": [{"title": "Never smoke", "text": "CHRNA5"}]}},
        blood_type_result={"available": True, "combined": "A+",
            "abo": {"phenotype": "A", "genotype": "A/O"}, "rhd": {"status": "Rh-positive"},
            "secretor": {"secretor_status": "Non-secretor"}},
        bloodwork_result={"clinical": {"advanced": {"biological_age": {"phenoage": 18.0, "accel": -6.5}},
            "flags": [{"name": "LDL", "value": 132, "note": "elevated"}]}},
    )


def test_module_context_includes_every_populated_module():
    ctx = analyze._build_module_context(**_sample_results())
    # Every populated module must keep representation (per-section budgeting).
    for tag in ["HOLISTIC SYNTHESIS", "IMMUNOGENETICS", "NEUROCHEMISTRY",
                "ADDICTION GENETICS", "BLOOD TYPE", "BLOODWORK"]:
        assert f"[{tag}" in ctx, f"missing section {tag}"  # prefix (titles may add detail)


def test_module_context_respects_total_budget():
    ctx = analyze._build_module_context(**_sample_results(), max_chars=400)
    assert len(ctx) <= 400 + 40  # allow the truncation marker


def test_module_context_empty_when_nothing_available():
    assert analyze._build_module_context() == ""


def test_module_ai_spec_keys_match_report_section_ids():
    # Every per-module AI key must be a real report section anchor id, else the
    # renderer can't attach the interpretation.
    valid_ids = {"holistic-synthesis", "immunogenetics", "neurochemistry",
                 "addiction-genetics", "deep-ancestry", "blood-type",
                 "family-planning", "polygenic-traits", "environmental-optimization",
                 "clinical-variants", "novel-variants", "value-of-information"}
    assert set(analyze._MODULE_AI_SPECS) <= valid_ids


def test_ai_grounding_guard_flags_fabricated_numbers():
    # Anti-hallucination guardrail: numbers absent from the source data are flagged.
    ctx = "CAD PRS high; value $57,000; WTP $100,000; 3% discount."
    clean = "Your CAD risk is elevated; the modelled value is about $57,000 at 3%."
    dirty = "Risk is 87% and cost $1,234,567 per a 2019 study; MTHFR raises it 42%."
    assert analyze._ground_ai_output(clean, ctx) == []
    flagged = analyze._ground_ai_output(dirty, ctx)
    assert set(flagged) >= {"87", "42", "1234567"}


def test_module_ai_spec_kwargs_are_valid_context_params():
    import inspect
    params = set(inspect.signature(analyze._build_module_context).parameters)
    for spec in analyze._MODULE_AI_SPECS.values():
        assert spec["kwarg"] in params, f"{spec['kwarg']} not a _build_module_context param"


def test_prepend_ai_injects_after_h2():
    section = '<section id="immunogenetics"><h2>Immunogenetics</h2><p>body</p></section>'
    out = renderers._prepend_ai_interpretation(section, "Some AI text.")
    assert "AI interpretation" in out
    # must appear after the heading, before the body
    assert out.index("AI interpretation") > out.index("</h2>")
    assert out.index("AI interpretation") < out.index("<p>body")


def test_prepend_ai_noop_without_text_or_section():
    section = "<section><h2>X</h2></section>"
    assert renderers._prepend_ai_interpretation(section, None) == section
    assert renderers._prepend_ai_interpretation("", "text") == ""
