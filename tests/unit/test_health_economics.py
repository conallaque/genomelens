"""
Unit tests for `health_economics.py`.

Critical invariants:
  1. ROI / payback / NPV math is exact on known cases. NPV uses a real 3%
     discount, so the warfarin case is ~42.1k, NOT the undiscounted 44.7k.
  2. Only *actionable* (atypical) PGx phenotypes and *elevated* PRS tiers
     become findings; Normal / Indeterminate / Average are suppressed.
  3. scale_to_clinic / scale_to_payer aggregate correctly and degrade
     gracefully to a `note` when there are no findings.
  4. Missing cost data and empty inputs never raise.
  5. End-to-end: a profile with findings produces an HTML section.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import health_economics as he
import renderers

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


# ── Core math ────────────────────────────────────────────────────────────────

def test_calculate_roi_known_case():
    # CYP2C9 warfarin: $15,000 averted / $300 test = 50:1
    assert he.calculate_roi(300, 15_000) == 50.0


def test_calculate_roi_zero_cost_is_none():
    assert he.calculate_roi(0, 15_000) is None
    assert he.calculate_roi(-5, 15_000) is None


def test_payback_months_known_case():
    # 300 / 15000 * 12 = 0.24 months
    assert he.calculate_payback_months(300, 15_000) == 0.24


def test_payback_zero_outcome_is_none():
    assert he.calculate_payback_months(300, 0) is None


def test_npv_one_time_cost_is_discounted_not_naive():
    # Benefit 15k/yr discounted at 3% over 3 yrs, minus 300 upfront.
    # Σ 15000/1.03^t (t=1..3) = 42429.17 ; − 300 = 42129.17.
    # Explicitly NOT the undiscounted 3*15000-300 = 44700.
    npv = he.calculate_npv(300, 15_000, recurring_cost=False)
    assert npv == pytest.approx(42129.17, abs=0.5)
    assert npv != 44_700


def test_npv_recurring_cost_discounts_both_streams():
    # Statin: 500/yr cost AND 250k/yr benefit, both discounted over 3 yrs.
    benefit = sum(250_000 / 1.03 ** t for t in range(1, 4))
    spend = sum(500 / 1.03 ** t for t in range(1, 4))
    assert he.calculate_npv(500, 250_000, recurring_cost=True) == pytest.approx(
        round(benefit - spend, 2), abs=0.5
    )


# ── Finding extraction ─────────────────────────────────────────────────────

def test_actionable_pgx_classification():
    assert he._is_actionable_pgx("Poor Metabolizer (PM)")
    assert he._is_actionable_pgx("Rapid Metabolizer (RM)")
    assert he._is_actionable_pgx("Non-expressor (PM)")
    # Not actionable:
    assert not he._is_actionable_pgx("Normal Metabolizer (NM)")
    assert not he._is_actionable_pgx("Normal Warfarin Sensitivity")
    assert not he._is_actionable_pgx("Indeterminate — no defining variants typed")
    assert not he._is_actionable_pgx("")


def test_pgx_normal_phenotype_suppressed():
    findings = {"pgx_summary": {"CYP2C9": {"phenotype": "Normal Metabolizer (NM)"}}}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 0
    assert res["status"] == "no_findings"


def test_pgx_actionable_phenotype_included_with_correct_roi():
    findings = {"pgx_summary": {"CYP2C9": {"phenotype": "Poor Metabolizer (PM)"}}}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 1
    f = res["findings_with_economics"][0]
    assert f["roi"] == 50.0
    assert f["confidence"] == "high"
    assert "warfarin" in f["finding"]


def test_unknown_gene_without_cost_data_is_suppressed():
    findings = {"pgx_summary": {"FAKE9": {"phenotype": "Poor Metabolizer (PM)"}}}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 0


def test_prs_only_elevated_tiers_fire():
    findings = {"prs_summary": {
        "Coronary Artery Disease": {"tier": "Elevated", "percentile": 89.0},
        "Type 2 Diabetes": {"tier": "Below Average", "percentile": 8.6},
    }}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 1
    assert "coronary" in res["findings_with_economics"][0]["finding"].lower()


def test_apoe_e4_genotype_fires():
    findings = {"apoe_genotype": "e3/e4"}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 1
    assert "APOE" in res["findings_with_economics"][0]["finding"]


def test_apoe_non_e4_genotype_suppressed():
    findings = {"apoe_genotype": "e3/e3"}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] == 0


# ── Scaling ─────────────────────────────────────────────────────────────────

def _five_finding_profile():
    return {
        "pgx_summary": {
            "CYP2C9": {"phenotype": "Poor Metabolizer (PM)"},
            "CYP2C19": {"phenotype": "Rapid Metabolizer (RM)"},
            "TPMT": {"phenotype": "Intermediate Metabolizer (IM)"},
        },
        "prs_summary": {
            "Coronary Artery Disease": {"tier": "Elevated", "percentile": 90},
            "Type 2 Diabetes": {"tier": "High", "percentile": 95},
        },
    }


def test_scale_to_clinic_100_patients_5_findings():
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    assert res["n_findings"] == 5
    clinic = he.scale_to_clinic(res, patient_count=100)
    assert clinic["patient_count"] == 100
    assert clinic["n_findings"] == 5
    # totals are per-patient average * count
    assert clinic["total_cost"] == pytest.approx(clinic["avg_cost_per_patient"] * 100)
    assert clinic["total_benefit"] == pytest.approx(clinic["avg_benefit_per_patient"] * 100)
    assert clinic["revenue_model_monthly"] == he.CLINIC_REVENUE_MONTHLY
    assert clinic["gross_margin"] == he.CLINIC_GROSS_MARGIN
    assert clinic["avg_roi"] > 0


def test_scale_to_payer_100k_members():
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    payer = he.scale_to_payer(res, member_population=100_000)
    assert payer["member_population"] == 100_000
    assert payer["affected_members"] > 0
    assert payer["total_benefit"] > payer["total_cost"]
    assert payer["roi"] > 0
    assert payer["cost_per_qaly"] is not None and payer["cost_per_qaly"] > 0
    # per-finding rows sum to the aggregate cost
    assert sum(p["total_cost"] for p in payer["per_finding"]) == pytest.approx(
        payer["total_cost"], abs=1
    )


# ── Edge cases ──────────────────────────────────────────────────────────────

def test_zero_findings_graceful():
    res = he.analyze_health_economics({}, pd.DataFrame())
    assert res["status"] == "no_findings"
    assert res["findings_with_economics"] == []
    assert "note" in res["clinic_dashboard"]
    assert "note" in res["payer_impact"]


def test_none_findings_does_not_raise():
    res = he.analyze_health_economics(None, pd.DataFrame())
    assert res["status"] == "no_findings"


def test_missing_cost_data_suppressed_not_error():
    # gene present in pgx summary but absent from PGX_ECONOMICS
    findings = {"pgx_summary": {"NOTACPIC": {"phenotype": "Poor Metabolizer (PM)"}}}
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["status"] == "no_findings"


def test_summary_markdown_renders_top3():
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    md = he.generate_economics_summary(res)
    assert "ROI" in md
    # disclaimer present
    assert "estimates" in md.lower()


def test_main_summary_gated_to_high_confidence():
    """Mixed profile: high-confidence PGx + moderate PRS/APOE. The main
    summary (markdown top-3 and HTML headline) must show only the
    high-confidence findings, while the detailed table keeps all."""
    profile = {
        "pgx_summary": {"CYP2C9": {"phenotype": "Poor Metabolizer (PM)"}},  # high
        "prs_summary": {"Coronary Artery Disease": {"tier": "Elevated", "percentile": 90}},  # moderate
        "apoe_genotype": "e3/e4",  # moderate
    }
    res = he.analyze_health_economics(profile, pd.DataFrame())
    assert res["n_findings"] == 3
    assert len(res["high_confidence"]) == 1  # only the PGx finding

    md = he.generate_economics_summary(res)
    assert "high-confidence" in md.lower()
    # CAD/APOE (moderate) must NOT appear in the gated top section
    headline_section = md.split("### Clinic")[0]
    assert "CYP2C9" in headline_section
    assert "coronary" not in headline_section.lower()
    assert "APOE" not in headline_section

    html = renderers.build_economics_html(res)
    # headline callout shows only the high-confidence finding...
    headline = html.split('class="tbl-wrap"')[0]
    assert "CYP2C9" in headline
    assert "APOE" not in headline
    # ...but the full table still lists all three findings.
    assert "APOE" in html
    assert "coronary" in html.lower()


# ── End-to-end: module → renderer → HTML section ─────────────────────────────

def test_renderer_empty_returns_blank():
    assert renderers.build_economics_html(None) == ""
    assert renderers.build_economics_html({"findings_with_economics": []}) == ""


def test_end_to_end_html_section_appears():
    """Run on the saved tier1_results.json (produced from test_genome.txt) and
    confirm the economics section is rendered into HTML."""
    tier1_path = PROJECT_ROOT / "tier1_results.json"
    if not tier1_path.exists():
        pytest.skip("tier1_results.json not present — run analyze.py first")
    d = json.loads(tier1_path.read_text())
    findings = {
        "pgx_summary": d.get("pgx_summary", {}),
        "prs_summary": d.get("prs_summary", {}),
        "apoe_genotype": d.get("apoe_genotype"),
    }
    res = he.analyze_health_economics(findings, pd.DataFrame())
    assert res["n_findings"] >= 1  # CAD PRS at 89th percentile guarantees ≥1
    html = renderers.build_economics_html(res)
    assert 'id="health-economics"' in html
    assert "Health Economics" in html
    assert "ROI" in html
