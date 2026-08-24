"""Unit tests for personalized_plan.py — the master dashboard synthesiser."""

from __future__ import annotations

from report import personalized_plan as pp

# ── Stub inputs (matching shapes from the V6 modules) ────────────────────────

def _supp_stub() -> dict:
    return {
        "status": "ok",
        "tiers": {
            "essential": [{
                "name": "Vitamin D3 + K2 (MK-7)", "dose": "4000 IU",
                "timing": "Morning with breakfast", "tier": "essential",
            }],
            "recommended": [{
                "name": "Curcumin", "dose": "500 mg",
                "timing": "With meals 1-2× daily", "tier": "recommended",
            }],
            "optional": [], "avoid": [], "chip_gap": [],
        },
    }


def _ex_stub() -> dict:
    return {
        "status": "ok",
        "power_endurance": {"bias": "Endurance-dominant",
                             "ratio_pct_power": 30, "ratio_pct_endurance": 70},
        "recovery": {"speed": "Slow"},
        "chronotype": {"chronotype": "Morning", "optimal_window": "07:00 – 11:00"},
    }


def _nu_stub() -> dict:
    return {
        "status": "ok",
        "macros": {"pct_carbs": 35, "pct_fat": 35, "pct_protein": 30},
        "caffeine": {"metabolism": "Slow", "limit_mg": 200, "cutoff_time": "12:00"},
        "alcohol": {"risk": "Standard"},
        "daily_template": [
            {"meal": "Breakfast", "example": "Eggs + spinach + coffee"},
        ],
    }


def _bw_stub() -> dict:
    return {
        "status": "ok", "n_matched": 4, "n_confirmed": 2,
        "n_partial": 0, "n_diverged": 2, "accuracy_pct": 50.0,
        "rows": [{
            "trait": "C-Reactive Protein", "verdict": "Diverged",
            "predicted": 1.8, "actual": 4.6, "unit": "mg/L",
            "delta_sd": 1.87,
            "interpretation": "Strong non-genetic inflammation driver.",
        }, {
            "trait": "25-OH Vitamin D", "verdict": "Diverged",
            "predicted": 19, "actual": 14, "unit": "ng/mL",
            "delta_sd": -1.55, "interpretation": "Below predicted.",
        }],
    }


def _phewas_stub() -> dict:
    return {"headline": [
        {"trait": "LDL cholesterol", "tier": "Above average",
         "percentile": 82, "predicted_value": 138, "unit": "mg/dL"},
    ]}


# ── Core synthesis ───────────────────────────────────────────────────────────

def test_morning_actions_pull_from_supplements_and_nutrition() -> None:
    plan = pp.build_personalized_plan(
        supplement_result=_supp_stub(),
        exercise_result=_ex_stub(),
        nutrition_result=_nu_stub(),
    )
    morning = plan["morning_actions"]
    # Morning-timed supplement
    assert any("Vitamin D3" in a for a in morning)
    # Caffeine cap (from nutrition.caffeine)
    assert any("200 mg" in a for a in morning)
    # Breakfast example
    assert any("Breakfast" in a for a in morning)
    # Morning training window
    assert any("07:00" in a for a in morning)


def test_evening_actions_pull_evening_timed_supplements() -> None:
    """Add an evening-timed supplement and verify it lands in evening_actions."""
    supp = _supp_stub()
    supp["tiers"]["recommended"].append({
        "name": "Magnesium Glycinate", "dose": "400 mg",
        "timing": "Evening (calming)", "tier": "recommended",
    })
    plan = pp.build_personalized_plan(supplement_result=supp,
                                       nutrition_result=_nu_stub())
    assert any("Magnesium" in a for a in plan["evening_actions"])


def test_evening_window_added_for_evening_chronotype() -> None:
    ex = _ex_stub()
    ex["chronotype"]["chronotype"] = "Evening"
    ex["chronotype"]["optimal_window"] = "16:00 – 20:00"
    plan = pp.build_personalized_plan(exercise_result=ex)
    assert any("16:00" in a for a in plan["evening_actions"])


def test_reconciliation_maps_diverged_biomarkers_to_supplements() -> None:
    """The headline V7 cross-module synthesis: for each Diverged lab row, find
    supplements in the user's stack that address it."""
    plan = pp.build_personalized_plan(
        supplement_result=_supp_stub(),
        bloodwork_result=_bw_stub(),
    )
    recon = plan["reconciliation"]
    assert len(recon) == 2

    crp_row = next(r for r in recon if r["biomarker"] == "C-Reactive Protein")
    crp_supp_names = [s["name"] for s in crp_row["supplements_in_play"]]
    assert any("Curcumin" in n for n in crp_supp_names)

    d_row = next(r for r in recon if r["biomarker"] == "25-OH Vitamin D")
    d_supp_names = [s["name"] for s in d_row["supplements_in_play"]]
    assert any("Vitamin D3" in n for n in d_supp_names)


def test_no_reconciliation_when_no_bloodwork() -> None:
    plan = pp.build_personalized_plan(supplement_result=_supp_stub())
    assert plan["reconciliation"] == []


def test_headlines_populate_for_each_pillar() -> None:
    plan = pp.build_personalized_plan(
        supplement_result=_supp_stub(), exercise_result=_ex_stub(),
        nutrition_result=_nu_stub(),  bloodwork_result=_bw_stub(),
        phewas_result=_phewas_stub(),
    )
    for k in ("supplements", "exercise", "nutrition", "bloodwork"):
        assert plan["headlines"][k] not in ("", "—")


def test_render_html_smoke() -> None:
    plan = pp.build_personalized_plan(
        supplement_result=_supp_stub(), exercise_result=_ex_stub(),
        nutrition_result=_nu_stub(), bloodwork_result=_bw_stub(),
        phewas_result=_phewas_stub(),
    )
    html = pp.render_plan_html(plan)
    assert "Personalised Plan" in html
    assert "Pillars at a glance" in html
    assert "supplements.html" in html
    # Reconciliation section appears only when there are diverged rows
    assert "highest-leverage" in html
