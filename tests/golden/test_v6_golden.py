"""
Golden-snapshot regression tests for the V6 modules.

These lock the exact structured output of each module against a committed JSON
snapshot. Any intentional behaviour change must be paired with a deliberate
snapshot update:

    pytest --snapshot-update tests/golden/test_v6_golden.py

…and the diff in `tests/snapshots/*.json` must be reviewed in the PR.

The fixtures normalise out volatile fields (timestamps, UUIDs) so the
comparisons are stable across runs.
"""

from __future__ import annotations

import pytest

import fhir_export as fx
import personalized_plan as pp
from wellness import bloodwork as bw
from wellness import exercise as ex
from wellness import nutrition as nu
from wellness import supplements as sup

pytestmark = pytest.mark.golden


# ── Single-module snapshots ──────────────────────────────────────────────────

def test_golden_supplements(synthetic_snps_df, assert_snapshot) -> None:
    result = sup.build_supplement_stack(snps_df=synthetic_snps_df)
    assert_snapshot("supplements_default", result)


def test_golden_exercise(synthetic_snps_df, assert_snapshot) -> None:
    result = ex.analyze_exercise(synthetic_snps_df)
    assert_snapshot("exercise_default", result)


def test_golden_nutrition(synthetic_snps_df, assert_snapshot) -> None:
    result = nu.analyze_nutrition(synthetic_snps_df)
    assert_snapshot("nutrition_default", result)


def test_golden_bloodwork(tmp_path, assert_snapshot) -> None:
    import json
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({
        "ldl": 145, "hdl": 52, "crp": 4.6, "vitamin_d": 14,
        "ferritin": 220, "tsh": 1.9,
    }))
    phewas = {"traits": {
        "LDL cholesterol": {
            "mean": 110, "sd": 32, "unit": "mg/dL", "category": "Lipids",
            "result": {"status": "ok", "predicted_value": 128,
                       "tier": "Above average", "callability_pct": 70,
                       "n_used": 6, "n_total": 8},
        },
        "HDL cholesterol": {
            "mean": 55, "sd": 15, "unit": "mg/dL", "category": "Lipids",
            "result": {"status": "ok", "predicted_value": 50,
                       "tier": "Average", "callability_pct": 80,
                       "n_used": 5, "n_total": 6},
        },
        "C-Reactive Protein": {
            "mean": 1.5, "sd": 1.5, "unit": "mg/L", "category": "Inflammation",
            "result": {"status": "ok", "predicted_value": 2.1,
                       "tier": "Average", "callability_pct": 75,
                       "n_used": 3, "n_total": 4},
        },
        "25-OH Vitamin D": {
            "mean": 25, "sd": 10, "unit": "ng/mL", "category": "Vitamins",
            "result": {"status": "ok", "predicted_value": 19,
                       "tier": "Below average", "callability_pct": 75,
                       "n_used": 3, "n_total": 4},
        },
        "Iron / ferritin": {
            "mean": 150, "sd": 100, "unit": "ng/mL", "category": "Vitamins",
            "result": {"status": "ok", "predicted_value": 170,
                       "tier": "Average", "callability_pct": 75,
                       "n_used": 3, "n_total": 3},
        },
        "TSH": {
            "mean": 2.0, "sd": 0.9, "unit": "mIU/L", "category": "Thyroid",
            "result": {"status": "ok", "predicted_value": 2.0,
                       "tier": "Average", "callability_pct": 100,
                       "n_used": 2, "n_total": 2},
        },
    }}
    result = bw.compare_bloodwork(str(p), phewas)
    assert_snapshot("bloodwork_default", result)


def test_golden_fhir_bundle(assert_snapshot) -> None:
    pgx = {
        "per_gene": {"CYP2D6": {
            "long_name": "CYP2D6", "phenotype": "Intermediate Metaboliser",
            "phenotype_code": "IM", "activity_score": 1.25,
            "callable_variants": 3, "total_variants": 5,
            "callability_pct": 60.0,
            "variant_calls": [
                {"rsid": "rs1065852", "called": True, "dosage": 1, "genotype": "GA"},
            ],
            "cpic_guideline": "CPIC v2024.1", "is_binary": False,
        }},
        "actionable_findings": [],
    }
    result = fx.build_fhir_bundle(
        pgx_result=pgx, apoe_genotype="ε3/ε4", file_label="synthetic.csv",
    )
    # Compare the bundle, not the timestamp-laden top-level wrapper
    assert_snapshot("fhir_bundle_default", result["bundle"])


# ── Cross-module integration snapshot ────────────────────────────────────────

def test_golden_personalized_plan_end_to_end(synthetic_snps_df, assert_snapshot,
                                              tmp_path) -> None:
    """
    Full V6 pipeline: supplements + exercise + nutrition feeding the master
    dashboard. Locks the cross-module synthesis output.
    """
    supp = sup.build_supplement_stack(snps_df=synthetic_snps_df)
    exercise = ex.analyze_exercise(synthetic_snps_df)
    nutrition = nu.analyze_nutrition(synthetic_snps_df)
    plan = pp.build_personalized_plan(
        supplement_result=supp,
        exercise_result=exercise,
        nutrition_result=nutrition,
    )
    assert_snapshot("plan_end_to_end", plan)
