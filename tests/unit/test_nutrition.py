"""Unit tests for nutrition.py."""

from __future__ import annotations

import pandas as pd

from wellness import nutrition as nu


def test_no_data_when_df_is_none() -> None:
    assert nu.analyze_nutrition(None) == {"status": "no_data"}


def test_macro_ratios_sum_to_100(synthetic_snps_df) -> None:
    result = nu.analyze_nutrition(synthetic_snps_df)
    m = result["macros"]
    assert m["pct_carbs"] + m["pct_fat"] + m["pct_protein"] == 100


def test_aldh2_homozygous_variant_avoids_alcohol_entirely() -> None:
    """rs671 AA — non-functional ALDH2, acetaldehyde accumulates dangerously."""
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs671"])
    result = nu.analyze_nutrition(df)
    assert result["alcohol"]["risk"] == "Avoid entirely"
    assert "abstinence" in result["alcohol"]["guidance"].lower()


def test_aldh2_wild_type_standard_moderation() -> None:
    df = pd.DataFrame({"genotype": ["GG"]}, index=["rs671"])
    result = nu.analyze_nutrition(df)
    assert result["alcohol"]["risk"] == "Standard"


def test_lct_persistence_means_dairy_ok() -> None:
    """rs4988235 T-allele = lactase persistence."""
    df = pd.DataFrame({"genotype": ["CT"]}, index=["rs4988235"])
    result = nu.analyze_nutrition(df)
    assert result["lactose"]["tolerance"] == "Persistent"


def test_lct_homozygous_ancestral_is_non_persistent() -> None:
    df = pd.DataFrame({"genotype": ["CC"]}, index=["rs4988235"])
    result = nu.analyze_nutrition(df)
    assert "non-persistent" in result["lactose"]["tolerance"].lower()


def test_cyp1a2_slow_metaboliser_caps_caffeine() -> None:
    df = pd.DataFrame({"genotype": ["CC"]}, index=["rs762551"])
    result = nu.analyze_nutrition(df)
    assert "Slow" in result["caffeine"]["metabolism"]
    assert result["caffeine"]["limit_mg"] <= 200


def test_daily_template_has_meals(synthetic_snps_df) -> None:
    result = nu.analyze_nutrition(synthetic_snps_df)
    meals = [m["meal"] for m in result["daily_template"]]
    assert "Breakfast" in meals
    assert "Lunch" in meals
    assert "Dinner" in meals


def test_render_html_includes_macro_bar(synthetic_snps_df) -> None:
    result = nu.analyze_nutrition(synthetic_snps_df)
    html = nu.render_nutrition_html(result)
    assert "Macronutrient Ratio" in html
    assert "nu-macro-bar" in html


def test_advanced_sections_not_gated_on_polygenic_scores() -> None:
    """Regression: _render_advanced_sections must not drop the whole advanced
    layer (dashboard/inflammation/detox/protocols) when polygenic_scores is
    absent. Each sub-section guards on its own data."""
    from wellness import nutrition as nu
    no_pgs = {
        "detoxification": {
            "phase1_typed": ["CYP1A2 fast"], "phase2_typed": ["GST null"],
            "cruciferous_target_servings_per_week": 5, "guidance": "Eat crucifers.",
        },
    }
    html = nu._render_advanced_sections(no_pgs)
    assert "Detoxification" in html, "detox section dropped when no pgs"
    assert "Polygenic Scores" not in html, "empty pgs table rendered with no pgs"
    # and pgs table still renders when present
    with_pgs = {"polygenic_scores": {"obesity": {"percentile": 70, "tier": "High", "coverage": 85}}}
    assert "Polygenic Scores" in nu._render_advanced_sections(with_pgs)
