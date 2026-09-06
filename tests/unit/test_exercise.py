"""Unit tests for exercise.py."""

from __future__ import annotations

import pandas as pd

from wellness import exercise as ex


def test_no_data_when_df_is_none() -> None:
    assert ex.analyze_exercise(None) == {"status": "no_data"}


def test_full_synthetic_chip_produces_all_pillars(synthetic_snps_df) -> None:
    result = ex.analyze_exercise(synthetic_snps_df)
    assert result["status"] == "ok"
    # Every pillar must populate
    for key in ("power_endurance", "injury_risk", "recovery",
                "chronotype", "cognitive", "weekly_template"):
        assert key in result, f"missing pillar: {key}"


def test_clock_evening_genotype_drives_evening_window() -> None:
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs1801260"])
    result = ex.analyze_exercise(df)
    assert "Evening" in result["chronotype"]["chronotype"]
    assert "16:00" in result["chronotype"]["optimal_window"]


def test_actn3_homozygous_x_is_endurance_biased() -> None:
    """ACTN3 TT (X/X) — α-actinin-3 deficiency, no Z-line αα → endurance bias."""
    df = pd.DataFrame({"genotype": ["TT"]}, index=["rs1815739"])
    result = ex.analyze_exercise(df)
    assert "Endurance" in result["power_endurance"]["bias"]


def test_power_endurance_ratio_sums_to_100(synthetic_snps_df) -> None:
    result = ex.analyze_exercise(synthetic_snps_df)
    pe = result["power_endurance"]
    assert pe["ratio_pct_power"] + pe["ratio_pct_endurance"] == 100


def test_bdnf_met_met_signals_extra_aerobic_importance() -> None:
    """rs6265 AA = Met/Met homozygous — reduced activity-dependent BDNF."""
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs6265"])
    result = ex.analyze_exercise(df)
    assert result["cognitive"]["genotype"] == "Met/Met"
    assert "especially important" in result["cognitive"]["summary"].lower()


def test_weekly_template_has_seven_days(synthetic_snps_df) -> None:
    result = ex.analyze_exercise(synthetic_snps_df)
    assert len(result["weekly_template"]) == 7
    days = [d["day"] for d in result["weekly_template"]]
    assert days == ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def test_render_html_smoke(synthetic_snps_df) -> None:
    result = ex.analyze_exercise(synthetic_snps_df)
    html = ex.render_exercise_html(result)
    assert "<html" in html.lower()
    assert "Personalized Exercise" in html
    assert result["power_endurance"]["bias"] in html


def test_sport_plan_weeks_render_all_schemas() -> None:
    """Regression: endurance plans (marathon/cycling/tri) use
    weekly_mileage_km/key_workouts, strength plans use schedule/accessory.
    The week renderer must surface every schema, never a blank body."""
    protocols = {
        "sport_specific_plans": [
            {"sport": "Marathon", "weeks": [
                {"phase": "Base (wk 1-6)",
                 "weekly_mileage_km": "40 to 60 km",
                 "key_workouts": "Long run +1km/wk; 1 tempo."}]},
            {"sport": "Powerlifting", "weeks": [
                {"phase": "Hypertrophy",
                 "schedule": "Squat 2x/wk 4x8.",
                 "accessory": "RDL, rows 3x10-12."}]},
        ]
    }
    html = ex._render_exercise_protocols({"protocols": protocols})
    for needle in ("40 to 60 km", "Long run +1km/wk", "Key workouts",
                   "Squat 2x/wk", "RDL, rows", "Accessory"):
        assert needle in html, f"missing: {needle}"
    # no phase header with an empty body
    assert "</strong> </li>" not in html
    assert "</strong></li>" not in html
