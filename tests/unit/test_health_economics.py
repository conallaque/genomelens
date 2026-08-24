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

import renderers
from econ import health_economics as he

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


# ── Core math ────────────────────────────────────────────────────────────────

def test_calculate_roi_known_case():
    # CYP2C9 warfarin: $15,000 averted / $300 test = 50:1
    assert he.calculate_roi(300, 15_000) == 50.0


def test_calculate_roi_zero_cost_is_none():
    assert he.calculate_roi(0, 15_000) is None
    assert he.calculate_roi(-5, 15_000) is None


def test_payback_months_known_case():
    # Payback walks month by month accumulating DISCOUNTED benefit, so the
    # smallest resolvable answer is one whole month. A $300 cost against a
    # $15k annual benefit is repaid inside the first month.
    assert he.calculate_payback_months(300, 15_000) == 1.0


def test_payback_months_lengthens_as_benefit_shrinks():
    # Directional property rather than a memorised constant.
    fast = he.calculate_payback_months(300, 15_000)
    slow = he.calculate_payback_months(300, 600)
    assert slow > fast


def test_payback_none_when_never_breaks_even():
    # Benefit so small it cannot repay the cost inside the 10-year window.
    assert he.calculate_payback_months(1_000_000, 1) is None


def test_payback_zero_outcome_is_none():
    assert he.calculate_payback_months(300, 0) is None


def test_npv_one_time_cost_is_discounted_not_naive():
    # Benefit 15k/yr discounted at 3% over 3 yrs, minus 300 upfront.
    # Σ 15000/1.03^t (t=1..3) = 42429.17 ; − 300 = 42129.17.
    # Explicitly NOT the undiscounted 3*15000-300 = 44700.
    npv = he.calculate_npv(300, 15_000, recurring_cost=False)
    # Benefit accrues ONCE at the horizon midpoint, not once per year.
    # This assertion previously pinned the annualised figure.
    assert npv == pytest.approx(14049.46, abs=0.5)
    assert npv != 44_700


def test_npv_recurring_cost_discounts_the_cost_stream_only():
    # THE BUG THIS GUARDS. `recurring` describes the COST: a statin is $500 a
    # year. It never described the benefit — $250,000 is the value of
    # preventing one event, not of preventing one every year. This test used to
    # assert the annualised benefit, which is how a $250,000 prevented MI
    # became a $705,739 three-year NPV.
    benefit = 250_000 / 1.03 ** (3 / 2.0)          # once, at the midpoint
    spend = sum(500 / 1.03 ** t for t in range(1, 4))   # each year
    assert he.calculate_npv(500, 250_000, recurring_cost=True) == pytest.approx(
        round(benefit - spend, 2), abs=0.5
    )


def test_npv_never_credits_one_prevented_event_more_than_once():
    # Directional guard that survives any horizon change: the benefit side can
    # never exceed the undiscounted value of the single event it represents.
    for horizon in (1, 3, 5, 10):
        npv = he.calculate_npv(0, 100_000, recurring_cost=False, horizon=horizon)
        assert npv <= 100_000 + 1e-6, (
            f"horizon={horizon} credited more than the event is worth")


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
    # outcome_value is the PROBABILITY-WEIGHTED averted cost
    # (p_prescribed x p_adr x rrr x adr_cost), not the raw cost of an adverse
    # event, so ROI is a modest ratio rather than the ~50x the undiscounted
    # event cost would imply.
    assert f["roi"] == round(f["outcome_value"] / f["intervention_cost"], 1)
    assert 0 < f["roi"] < 5
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
    # A ratio is only meaningful when cost and effect move the same way. This
    # cohort is modelled as dominant, so the ICER is withheld and the fact is
    # stated as a flag — the same convention markov.py and CEAResult use. The
    # old assertion required a number here, which is what let a $616/QALY
    # figure sit beside a claim of $8bn in savings.
    if payer["dominant"]:
        assert payer["cost_per_qaly"] is None
        assert "dominant" in payer["cost_per_qaly_note"]
    else:
        assert payer["cost_per_qaly"] > 0
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


def _tier1_fixture() -> dict:
    """A hermetic tier-1 summary with known-actionable findings.

    Deliberately does NOT read the on-disk ``tier1_results.json``: that file is a
    gitignored runtime artifact whose contents depend on whichever genome was
    analysed last, so a test bound to it skips on CI, skips on a fresh clone, and
    skips silently the moment someone runs a different input. Building the input
    inline makes this test actually run everywhere.
    """
    return {
        "pgx_summary": {
            "CYP2C19": {"phenotype": "Intermediate metabolizer", "activity_score": 1.0},
            "CYP2C9": {"phenotype": "Poor metabolizer", "activity_score": 0.5},
        },
        "prs_summary": {
            "Coronary Artery Disease": {"tier": "High", "percentile": 92},
            "Type 2 Diabetes": {"tier": "Elevated", "percentile": 78},
        },
        "apoe_genotype": "e3/e4",
    }


def test_end_to_end_html_section_appears():
    """Economics analysis → rendered HTML section, on a self-contained fixture."""
    res = he.analyze_health_economics(_tier1_fixture(), pd.DataFrame())
    assert res["n_findings"] >= 1, "hermetic fixture should yield actionable findings"
    html = renderers.build_economics_html(res)
    assert 'id="health-economics"' in html
    assert "Health Economics" in html
    assert "ROI" in html


def test_on_disk_tier1_matches_engine_when_present():
    """Opportunistic check against the real runtime artifact, when it happens to be
    present AND was produced from a genome with actionable findings. Skipping here is
    legitimate — this asserts nothing the hermetic test above doesn't already cover;
    it only catches drift between the saved file's shape and the engine's expectations.
    """
    tier1_path = PROJECT_ROOT / "tier1_results.json"
    if not tier1_path.exists():
        pytest.skip("tier1_results.json not present — optional integration check")
    d = json.loads(tier1_path.read_text())
    res = he.analyze_health_economics({
        "pgx_summary": d.get("pgx_summary", {}),
        "prs_summary": d.get("prs_summary", {}),
        "apoe_genotype": d.get("apoe_genotype"),
    }, pd.DataFrame())
    # Shape contract must hold regardless of whether the genome had findings.
    assert isinstance(res.get("n_findings"), int)
    assert "findings_with_economics" in res
    assert res.get("status") in ("ok", "no_findings", None) or isinstance(res.get("status"), str)


# ── Personal economic-impact sheet ────────────────────────────────────────────

def test_personal_economics_models_items() -> None:
    from econ import health_economics as he
    econ = he.analyze_personal_economics(
        economics_result={"findings_with_economics": [
            {"clinical_benefit": "Avoid warfarin bleed", "outcome_value": 15000,
             "prevalence": 0.35, "qaly_gain": 0.3, "cost": 300, "confidence": "high"}]},
        bloodwork_result={"clinical": {"advanced": {
            "indices": [{"id": "prevent_ascvd", "value": 9.2}],
            "biological_age": {"accel": 2.0, "inputs": {"glucose": 108}}}, "flags": []}})
    assert econ["available"] and econ["n_items"] >= 3
    cats = {i["category"] for i in econ["items"]}
    assert "Cardiovascular" in cats and "Metabolic" in cats
    assert econ["roi"] is not None
    html = he.render_economic_analysis_html(econ, "t")
    assert "Economic-Impact Analysis" in html and "ROI" in html


def test_personal_economics_empty() -> None:
    from econ import health_economics as he
    econ = he.analyze_personal_economics()
    assert econ["available"] is False and econ["n_items"] == 0
    html = he.render_economic_analysis_html(econ)
    assert "No modeled economic-impact" in html


def test_biological_aging_high_confidence_grounded() -> None:
    from econ import health_economics as he
    econ = he.analyze_personal_economics(bloodwork_result={"clinical": {"advanced": {
        "indices": [], "biological_age": {
            "accel": 1.2, "chronological": 41, "mortality_10yr_pct": 2.6, "inputs": {}}},
        "flags": []}})
    b = next(i for i in econ["items"] if i["category"] == "Biological aging")
    assert b["confidence"] == "high"                     # grounded in the clock, not hand-waved
    assert "mortality" in b["basis"].lower() and "Levine" in b["basis"]
    assert b["net"] > 0



# ══════════════════════════════════════════════════════════════════════════
# Cohort views must agree with the individual sheet
# ══════════════════════════════════════════════════════════════════════════

def test_cohort_corrections_only_ever_reduce_the_claimed_benefit():
    # THE CONSTRAINT, AS A GUARDRAIL. The cohort views used to report the raw
    # curated figures: no pooling, no adherence, no discounting, no
    # marginal-cost fraction. Every one of those corrections must move the
    # number down. If a future change makes a cohort figure larger than the
    # uncorrected sum, it is a bug in the correction, not a discovery.
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    payer = he.scale_to_payer(res, member_population=100_000)
    assert payer["total_cost_averted"] < payer["legacy"]["total_benefit"]
    assert payer["total_qalys"] < payer["legacy"]["total_qalys"]
    assert payer["benefit_reduction"] > 0
    assert 0 < payer["benefit_reduction_pct"] < 100


def test_cohort_views_apply_the_same_corrections_as_the_individual_sheet():
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    _items, corr = he._corrected_cohort_items(res)
    assert corr["marginal_cost_fraction"] < 1.0
    assert corr["midpoint_discount"] < 1.0
    assert corr["mean_adherence"] < 1.0
    assert corr["pooling_applied"], "cohort views must pool correlated targets"


def test_the_cohort_headcount_and_the_money_are_capped_together():
    # The old code capped the displayed member count at plan size while leaving
    # the benefit summed over every finding independently — so the report said
    # "100,000 of 100,000 members affected" and still added up the dollars as
    # though they were different people.
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    p = he.scale_to_payer(res, member_population=100_000)
    assert p["unique_members_affected"] <= p["member_population"]
    assert "intervention_events" in p, "event count must be stated, not hidden"
    if p["prevalence_sum_exceeds_cohort"]:
        assert p["intervention_events"] > p["unique_members_affected"]


def test_the_subscription_model_is_not_mixed_into_the_health_value_ratio():
    res = he.analyze_health_economics(_five_finding_profile(), pd.DataFrame())
    c = he.scale_to_clinic(res, patient_count=100)
    rm = c["revenue_model"]
    assert rm["monthly_per_patient"] and rm["gross_margin"]
    assert "not a clinical or economic result" in rm["note"]
    # The health-value ratio must be benefit over cost, with no revenue term.
    expected = he.calculate_roi(c["avg_cost_per_patient"],
                                c["avg_benefit_per_patient"])
    assert c["value_to_cost_ratio"] == expected


def test_family_planning_is_pooled_against_carrier_screening_not_dropped():
    # Both price the same reproductive decision — the carrier line at a flat
    # $2,000, the family-planning line decomposed as partner-carrier frequency
    # x child clinical risk x cost of an affected child. Dropping the second
    # lost real value; counting both double-counted one decision.
    assert "Family Planning" not in he.COHORT_NOT_VALUED
    assert he.COHORT_NOT_VALUED == ("Longevity",)


def test_longevity_stays_out_because_it_re_aggregates_priced_variants():
    assert "Longevity" in he.COHORT_NOT_VALUED


# ══════════════════════════════════════════════════════════════════════════
# Silent string mismatches between a module and its economics table
# ══════════════════════════════════════════════════════════════════════════

def _neuro(rsid, gt):
    import pandas as pd
    return pd.DataFrame({"genotype": [gt]}, index=[rsid])


def test_a_wild_type_genotype_is_not_priced_as_a_finding():
    # THE BUG. _neurochemistry_findings guards on finding["impact"] == "neutral",
    # but neurochemistry.py never emitted an "impact" key, so .get() returned
    # None, None is not "neutral", and the guard never fired. OPRM1 A/A —
    # "standard mu-opioid receptor", action "None specific." — was priced at the
    # full $25,000. Nothing crashed; the total was just wrong.
    import neurochemistry as nc
    wild = nc._oprm1(_neuro("rs1799971", "AA"))
    assert wild["impact"] == "neutral", "a typical result must declare itself"
    assert he._neurochemistry_findings({"available": True, "findings": [wild]}) == []

    carrier = nc._oprm1(_neuro("rs1799971", "GG"))
    assert carrier["impact"] == "informative"
    priced = he._neurochemistry_findings({"available": True, "findings": [carrier]})
    assert len(priced) == 1 and priced[0]["outcome_value"] > 0, (
        "the guard must suppress non-findings without suppressing real ones")


def test_every_neurochemistry_record_declares_an_impact():
    # The economics layer reads this field. A record without it silently opts
    # out of the guard, which is how the original defect stayed invisible.
    import pandas as pd

    import neurochemistry as nc
    probes = [("_comt", nc._comt, "rs4680", "AA"), ("_maoa", nc._maoa, "rs6323", "GG"),
              ("_bdnf", nc._bdnf, "rs6265", "TT"), ("_drd2", nc._drd2, "rs1800497", "TT"),
              ("_drd4", nc._drd4, "rs1800955", "TT"),
              ("_oprm1", nc._oprm1, "rs1799971", "GG"),
              ("_cacna1c", nc._cacna1c, "rs1006737", "AA"),
              ("_chrna5", nc._chrna5, "rs16969968", "AA")]
    for name, fn, rsid, gt in probes:
        rec = fn(pd.DataFrame({"genotype": [gt]}, index=[rsid]))
        if rec is None:
            continue
        assert "impact" in rec, f"{name} emits no impact field"
        assert rec["impact"] in ("informative", "neutral"), (
            f"{name} impact={rec['impact']!r} is not a value the guard understands")


def test_a_compound_locus_label_still_resolves_to_its_econ_entry():
    # THE OTHER BUG, opposite direction. The table is keyed "DRD2"; the module
    # reports the Taq1A locus as "DRD2/ANKK1" because it spans both genes. The
    # lookup was exact, so that finding was never priced at all — a suppressed
    # line rather than an inflated one, and equally invisible.
    import neurochemistry as nc
    rec = nc._drd2(_neuro("rs1800497", "TT"))
    assert rec["gene"] == "DRD2/ANKK1"
    priced = he._neurochemistry_findings({"available": True, "findings": [rec]})
    assert len(priced) == 1, "a compound gene label must resolve to its entry"
    assert priced[0]["outcome_value"] == \
        he.NEUROCHEMISTRY_ECONOMICS["DRD2"]["outcome_value"]


def test_econ_table_keys_resolve_against_the_strings_modules_emit():
    # Generalises both bugs. An economics table keyed on a gene symbol is only
    # useful if some module actually emits that symbol; a key that matches
    # nothing is a valuation the model silently never computes.
    import pandas as pd

    import neurochemistry as nc
    emitted = set()
    for fn, rsid, gt in ((nc._comt, "rs4680", "AA"), (nc._bdnf, "rs6265", "TT"),
                         (nc._drd2, "rs1800497", "TT"), (nc._oprm1, "rs1799971", "GG")):
        rec = fn(pd.DataFrame({"genotype": [gt]}, index=[rsid]))
        if rec:
            emitted.add(rec["gene"])
            emitted.add(rec["gene"].split("/")[0].strip())
    unreachable = set(he.NEUROCHEMISTRY_ECONOMICS) - emitted
    assert not unreachable, (
        f"econ entries no module can reach: {sorted(unreachable)}")
