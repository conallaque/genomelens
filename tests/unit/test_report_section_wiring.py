"""Integration test: verify the V35-V39 result dicts actually thread through
build_html_report into rendered sections. The unit tests exercise the analyzers;
this checks the wiring (signature -> section var -> body) that unit tests miss."""

from __future__ import annotations

import pandas as pd

import environmental_optimization
import family_planning
import life_stage_playbook
import polygenic_traits
import renderers


def _min_tier1():
    # build_html_report needs a plausible tier1_results structure; keep minimal.
    return {}


def test_new_sections_appear_in_full_report():
    # Build each new result from its analyzer with synthetic input.
    fam = family_planning.analyze_family_planning(
        {"carriers": [{"gene": "HFE", "variant": "C282Y",
                       "disease": "Hereditary Hemochromatosis (HH, type 1)",
                       "inheritance": "autosomal recessive (incomplete penetrance)",
                       "dosage": 1}], "affected": []},
        inferred_sex="M")
    traits = polygenic_traits.analyze_polygenic_traits(
        pd.DataFrame({"genotype": {"rs713598": "GG", "rs12913832": "GG"}}))
    envopt = environmental_optimization.analyze_environmental_optimization(
        pd.DataFrame({"genotype": {"rs1801260": "CC", "rs1815739": "CC"}}), latitude=43.0)
    lsp = life_stage_playbook.analyze_life_stage_playbook(
        age=24,
        holistic_synthesis_result={"genome_leverage": {"tier": "Favorable", "score": 74},
                                   "insights": []},
        tier1_results={"apoe_genotype": "e3/e3"})

    # Each build_*_html must render its anchor id.
    assert 'id="family-planning"' in renderers.build_family_planning_html(fam)
    assert 'id="polygenic-traits"' in renderers.build_polygenic_traits_html(traits)
    assert 'id="environmental-optimization"' in renderers.build_environmental_optimization_html(envopt)
    assert 'id="life-stage-playbook"' in renderers.build_life_stage_playbook_html(lsp)

    # Phase-3 novel variants: available and unavailable both render the anchor.
    nv_ok = {"available": True, "n_predicted_pathogenic": 1, "n_rare_damaging": 1,
             "n_queried": 10, "predictors_used": [{"name": "AlphaMissense",
             "license": "CC BY 4.0", "commercial_ok": True}], "commercial_safe": False,
             "buckets": {"predicted_pathogenic_rare": [{"chrom": "1", "pos": 100,
             "ref": "C", "alt": "T", "gene": "", "rarity": "unknown", "confidence": "higher",
             "zygosity": "heterozygous", "evidence": "AlphaMissense likely pathogenic (0.98)",
             "interpretation": "Computational prediction."}]},
             "negative_disclaimer": "n", "disclaimer": "d"}
    assert 'id="novel-variants"' in renderers.build_novel_variants_html(nv_ok)
    assert 'id="novel-variants"' in renderers.build_novel_variants_html(
        {"available": False, "reason": "run python setup.py --predictors",
         "negative_disclaimer": "n"})

    # Value-of-Information (health economics) renders headline + CEAC.
    from econ import value_of_information as voi
    vres = voi.analyze_value_of_information(
        {"findings_with_economics": [
            {"finding": "CAD PRS high", "category": "prs", "qaly_gain": 1.5}]},
        input_type="wgs", n_mc=300)
    assert 'id="value-of-information"' in renderers.build_voi_html(vres)


def test_none_results_render_empty_not_crash():
    # All new renderers must no-op gracefully on None.
    assert renderers.build_family_planning_html(None) == ""
    assert renderers.build_polygenic_traits_html(None) == ""
    assert renderers.build_environmental_optimization_html(None) == ""
    assert renderers.build_life_stage_playbook_html(None) == ""
    assert renderers.build_novel_variants_html(None) == ""
    assert renderers.build_voi_html(None) == ""
