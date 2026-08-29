"""The report payload and its consistency validator.

These tests are about *report drift*: the class of defect where every module
computes correctly and the assembled report still says two incompatible things.
They build payloads directly rather than running the pipeline, so a broken
invariant fails in a second rather than after a genome-scale run.

The validator's own first run against real fixture output produced a false
ERROR — it selected the budget-impact peak by absolute value, which reports a
year of maximum *savings* as the peak spend. That bug is pinned below, because a
validator that cries wolf gets switched off.
"""
from __future__ import annotations

import json

import pytest

from report.payload import (
    AnalysisBasis,
    ConditionResult,
    EconomicsReportPayload,
    FindingEconomics,
    PricingPath,
    ReferenceCase,
    ReportMetadata,
    TestingDecision,
    Uncertainty,
    build_report_payload,
    payload_to_json,
)
from report.validate import Severity, errors_in, format_report, validate_payload


def _checks(payload, name):
    return [f for f in validate_payload(payload) if f["check"] == name]


def _severities(payload, name):
    return {f["severity"] for f in _checks(payload, name)}


def _coherent() -> EconomicsReportPayload:
    """A payload that should raise no ERROR. λ·ΔQALY − ΔCost = 100000·0.032 + 4223."""
    return EconomicsReportPayload(
        metadata=ReportMetadata(is_synthetic=True, willingness_to_pay=100_000.0),
        reference_case=ReferenceCase(
            incremental_cost=-4223.0, incremental_qalys=0.032,
            nmb=7423.0, wtp=100_000.0, icer=None,
            icer_note="not defined — strategy dominates",
            dominance_status="dominant (more health, lower cost)"),
        uncertainty=Uncertainty(
            psa_available=True, psa_iterations=1500, psa_mean_nmb=7759.0,
            nmb_ci_low=1435.0, nmb_ci_high=20594.0,
            probability_cost_effective=0.9987, probability_cost_saving=0.956),
    )


# ── 1. schema construction ────────────────────────────────────────────────────

def test_empty_payload_constructs_and_validates():
    """Every section optional. A run that produced no economics must not crash
    the report; it must produce an empty payload that still serialises."""
    p = EconomicsReportPayload()
    assert p.findings == []
    assert p.reference_case.nmb == 0.0
    assert errors_in(validate_payload(p)) == []
    json.loads(payload_to_json(p))


def test_builder_tolerates_all_inputs_missing():
    p = build_report_payload(None, None, None)
    assert isinstance(p, EconomicsReportPayload)
    assert p.metadata.willingness_to_pay == 100_000.0


def test_payload_serialises_enums_as_strings():
    p = EconomicsReportPayload(
        findings=[FindingEconomics(display_name="x",
                                   pricing_path=PricingPath.VOI_PARAMETRIC)])
    blob = json.loads(payload_to_json(p))
    assert blob["findings"][0]["pricing_path"] == "voi_parametric"
    assert (blob["testing_decision"]["observed_basis"] == "observed_findings")


# ── 2. NMB identity ───────────────────────────────────────────────────────────

def test_nmb_identity_holds_on_coherent_payload():
    assert errors_in(validate_payload(_coherent())) == []


def test_nmb_identity_catches_a_broken_total():
    p = _coherent()
    p.reference_case.nmb = 9999.0
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "NMB identity" for e in errs), format_report(errs)


def test_nmb_identity_tolerates_engine_rounding():
    """The engine reports QALYs to four decimals; at a $100k threshold that is
    up to $5 of rounding on its own. A strict equality test would fire on
    correct arithmetic, so the tolerance has a floor."""
    p = _coherent()
    p.reference_case.nmb = 7423.0 - 4.0
    assert errors_in(validate_payload(p)) == []


# ── 3. dominance ──────────────────────────────────────────────────────────────

def test_dominance_requires_cheaper_and_better():
    p = _coherent()
    p.reference_case.incremental_cost = 500.0      # costs more
    p.reference_case.nmb = 100_000.0 * 0.032 - 500.0
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Dominance classification" for e in errs)


def test_dominance_rejects_dominated_label_on_a_winning_strategy():
    p = _coherent()
    p.reference_case.dominance_status = "dominated"
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Dominance classification" for e in errs)


def test_undefined_icer_without_a_note_warns():
    p = _coherent()
    p.reference_case.icer_note = ""
    assert Severity.WARNING in _severities(p, "ICER note")


# ── 4. net cash ───────────────────────────────────────────────────────────────

def test_net_cash_identity_per_finding():
    p = _coherent()
    p.findings = [FindingEconomics(
        display_name="APOE", pricing_path=PricingPath.CURATED_TABLE,
        medical_cost_averted=165.0, intervention_cost=100.0, net_cash=65.0)]
    assert errors_in(validate_payload(p)) == []

    p.findings[0].net_cash = 900.0
    errs = errors_in(validate_payload(p))
    assert any("Net-cash identity" in e["check"] for e in errs)


def test_health_economic_value_is_not_net_cash():
    """The two must stay distinct: one is money, the other includes health
    monetised at the threshold. Conflating them is the "your genome is worth $X"
    error the report is meant to avoid."""
    f = FindingEconomics(medical_cost_averted=84.0, intervention_cost=65.0,
                         net_cash=19.0, monetized_qaly_value=5607.0,
                         health_economic_value=5626.0)
    assert f.net_cash != f.health_economic_value
    assert f.health_economic_value == pytest.approx(
        f.net_cash + f.monetized_qaly_value)


# ── 5/6. pooling and the risk cap ─────────────────────────────────────────────

def test_pooling_must_reduce_never_inflate():
    p = _coherent()
    p.corrections.naive_cost_averted = 14692.0
    p.corrections.pooled_cost_averted = 5468.0
    assert errors_in(validate_payload(p)) == []

    p.corrections.pooled_cost_averted = 20000.0
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Pooling direction" for e in errs)


def test_adherence_cannot_exceed_trial_efficacy():
    p = _coherent()
    p.condition_results = [ConditionResult(
        condition="Pathogenic Carrier Screening", n_contributing_findings=3,
        naive_additive_rrr=0.94, pooled_efficacy_rrr=0.52,
        adherence_adjusted_rrr=0.36, adherence=0.65)]
    assert errors_in(validate_payload(p)) == []

    p.condition_results[0].adherence_adjusted_rrr = 0.80
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Adherence direction" for e in errs)


def test_pooled_efficacy_cannot_exceed_naive_sum():
    p = _coherent()
    p.condition_results = [ConditionResult(
        condition="X", n_contributing_findings=2, naive_additive_rrr=0.40,
        pooled_efficacy_rrr=0.55, adherence_adjusted_rrr=0.30, adherence=0.6)]
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Pooling direction" for e in errs)


def test_baseline_risk_and_adherence_must_be_probabilities():
    p = _coherent()
    p.condition_results = [ConditionResult(condition="X", baseline_risk=1.4,
                                           adherence=2.0)]
    errs = errors_in(validate_payload(p))
    assert sum(1 for e in errs if e["check"] == "Probability range") == 2


# ── 7. the WGS separation, which is the point of the exercise ─────────────────

def _wgs_payload() -> EconomicsReportPayload:
    p = _coherent()
    p.findings = [
        FindingEconomics(finding_id="a", display_name="predicted-pathogenic 17",
                         pricing_path=PricingPath.VOI_PARAMETRIC,
                         health_economic_value=2341.0, is_wgs_only=True),
        FindingEconomics(finding_id="b", display_name="predicted-pathogenic 12",
                         pricing_path=PricingPath.VOI_PARAMETRIC,
                         health_economic_value=1134.0, is_wgs_only=True),
    ]
    p.testing_decision = TestingDecision(
        incremental_chip_to_wgs_cost=200.0,
        prospective_expected_yield=0.018,
        prospective_number_needed_to_sequence=56,
        prospective_value_per_finding=25189.0,
        prospective_gross_expected_value=491.0,
        prospective_net_expected_value=291.0,
        observed_wgs_only_findings=2,
        observed_wgs_only_value=3475.0)
    return p


def test_observed_and_prospective_wgs_coexist_without_error():
    """The fixture genome is the case where both are non-zero. Neither
    suppresses the other, and the coexistence is recorded as INFO so the
    renderer is obliged to explain it."""
    p = _wgs_payload()
    assert errors_in(validate_payload(p)) == []
    assert Severity.INFO in _severities(p, "WGS basis coexistence")


def test_observed_wgs_value_must_be_a_sum_over_observed_findings():
    """The contamination guard. If a future change sources the observed figure
    from the population prior, the reconstruction stops matching and this
    fires."""
    p = _wgs_payload()
    p.testing_decision.observed_wgs_only_value = 491.0   # the prospective number
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "WGS observed value contamination" for e in errs)


def test_no_observed_findings_means_no_observed_value():
    p = _wgs_payload()
    p.findings = []
    p.testing_decision.observed_wgs_only_findings = 0
    p.testing_decision.observed_wgs_only_value = 3475.0
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "WGS observed value contamination" for e in errs)


def test_observed_count_must_match_the_findings_present():
    p = _wgs_payload()
    p.testing_decision.observed_wgs_only_findings = 5
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "WGS observed count" for e in errs)


def test_prospective_net_is_gross_minus_incremental_cost():
    p = _wgs_payload()
    p.testing_decision.prospective_net_expected_value = 450.0
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Prospective WGS identity" for e in errs)


def test_number_needed_to_sequence_matches_its_yield():
    p = _wgs_payload()
    p.testing_decision.prospective_number_needed_to_sequence = 12
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Number needed to sequence" for e in errs)


def test_wgs_basis_tags_cannot_be_swapped():
    p = _wgs_payload()
    p.testing_decision.prospective_basis = AnalysisBasis.OBSERVED_FINDINGS
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "WGS basis tagging" for e in errs)


# ── 8. budget impact ──────────────────────────────────────────────────────────

def _budget_rows():
    return [{"year": 1, "pmpm": 0.1167, "net_budget_impact": 1_400_000},
            {"year": 2, "pmpm": 0.1080, "net_budget_impact": 1_296_000},
            {"year": 3, "pmpm": -0.0133, "net_budget_impact": -160_000},
            {"year": 4, "pmpm": -0.0827, "net_budget_impact": -992_000},
            {"year": 5, "pmpm": -0.1633, "net_budget_impact": -1_960_000}]


def test_peak_budget_impact_is_the_highest_spend_not_the_largest_magnitude():
    """REGRESSION. The validator's first version selected the peak by absolute
    value and flagged correct engine output as an error: once a programme turns
    cost-saving, the final year has the largest magnitude and the smallest
    budget impact. Peak means worst year for the payer."""
    p = _coherent()
    p.advanced.budget_impact = {"available": True, "rows": _budget_rows(),
                                "peak_pmpm": 0.1167, "peak_year": 1}
    assert errors_in(validate_payload(p)) == []


def test_peak_pmpm_mismatch_is_an_error():
    p = _coherent()
    p.advanced.budget_impact = {"available": True, "rows": _budget_rows(),
                                "peak_pmpm": 0.05, "peak_year": 1}
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Budget-impact peak" for e in errs)


def test_peak_year_mismatch_is_an_error():
    p = _coherent()
    p.advanced.budget_impact = {"available": True, "rows": _budget_rows(),
                                "peak_pmpm": 0.1167, "peak_year": 4}
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Budget-impact peak year" for e in errs)


# ── the Markov cross-check must not become the headline ───────────────────────

def test_markov_flagged_as_reference_case_is_an_error():
    p = _coherent()
    p.structural_crosscheck.available = True
    p.structural_crosscheck.is_reference_case = True
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Structural cross-check misclassified" for e in errs)


def test_large_markov_divergence_warns_but_does_not_block():
    p = _coherent()
    p.structural_crosscheck.available = True
    p.structural_crosscheck.markov_qaly_gain = 0.436     # 13.6x the reference
    assert errors_in(validate_payload(p)) == []
    assert Severity.WARNING in _severities(p, "Structural cross-check divergence")


# ── differences that are correct ──────────────────────────────────────────────

def test_deterministic_and_psa_nmb_differ_without_being_an_error():
    p = _coherent()
    assert p.reference_case.nmb != p.uncertainty.psa_mean_nmb
    assert errors_in(validate_payload(p)) == []
    assert Severity.INFO in _severities(p, "Deterministic vs probabilistic NMB")


def test_per_finding_and_pooled_totals_differ_without_being_an_error():
    p = _coherent()
    p.findings = [FindingEconomics(display_name="a",
                                   pricing_path=PricingPath.VOI_PARAMETRIC,
                                   canonical_expected_nmb=5626.0)]
    p.condition_results = [ConditionResult(condition="c", nmb=1516.0)]
    assert errors_in(validate_payload(p)) == []
    assert Severity.INFO in _severities(p, "Per-finding vs pooled totals")


def test_retained_curated_figure_is_reported_as_audit_not_as_a_rival_nmb():
    """After canonicalisation a pathway carries ONE value. The curated figure
    survives only as an audit field, and a large divergence is INFO — worth
    explaining, not a defect, because the two quantities are different."""
    p = _coherent()
    p.findings = [FindingEconomics(
        display_name="PTPN22 R620W carrier",
        economic_pathway_id="pgx:ptpn22:screening:autoimmune",
        pathway_id_is_legacy=False,
        pricing_path=PricingPath.VOI_PARAMETRIC,
        canonical_expected_nmb=2509.0,
        economic_value_basis="parametric_expected_nmb",
        legacy_curated_value=5800.0,
        legacy_curated_value_basis="curated_prevalence_weighted_mixture")]
    assert errors_in(validate_payload(p)) == []
    assert Severity.INFO in _severities(p, "Legacy curated value diverges")


def test_canonical_records_on_legacy_identifiers_are_flagged():
    p = _coherent()
    p.findings = [FindingEconomics(
        display_name="APOE e4 carrier",
        economic_pathway_id="legacy:apoe:apoe_e4_carrier",
        pathway_id_is_legacy=True,
        pricing_path=PricingPath.VOI_PARAMETRIC,
        canonical_expected_nmb=3514.0)]
    assert errors_in(validate_payload(p)) == []
    assert Severity.WARNING in _severities(
        p, "Canonical records on legacy identifiers")


def test_unstandardised_findings_are_reported_not_hidden():
    p = _coherent()
    p.findings = [FindingEconomics(
        display_name="1 actionable wellness variant(s)", is_monetized=True,
        canonical_expected_nmb=None, legacy_curated_value=1767.0)]
    assert errors_in(validate_payload(p)) == []
    assert Severity.WARNING in _severities(
        p, "Findings without a standardised value")


# ── unmonetised, hypothetical, synthetic ──────────────────────────────────────

def test_unmonetized_findings_are_kept_not_dropped():
    p = EconomicsReportPayload(findings=[
        FindingEconomics(display_name="Reproductive genetics",
                         is_monetized=False,
                         reason_not_monetized="no proven intervention"),
        FindingEconomics(display_name="APOE", is_monetized=True,
                         health_economic_value=550.0)])
    assert len(p.unmonetized_findings) == 1
    assert p.unmonetized_findings[0].reason_not_monetized
    assert len(p.monetized_findings) == 1


def test_unmonetized_finding_contributes_no_value():
    f = FindingEconomics(display_name="x", is_monetized=False)
    assert f.health_economic_value == 0.0
    assert f.net_cash == 0.0


def test_hypothetical_findings_are_flagged_distinctly_from_unmonetized():
    """Awareness findings can carry a modelled value; not-monetised ones carry
    none. They are different states and the schema keeps them separate."""
    aware = FindingEconomics(display_name="a", is_monetized=True,
                             is_hypothetical_or_awareness=True,
                             health_economic_value=174.0)
    assert aware.is_monetized and aware.is_hypothetical_or_awareness


def test_synthetic_labelling_comes_from_metadata_not_prose():
    p = build_report_payload(None, None, None,
                             metadata={"is_synthetic": True,
                                       "input_label": "synthetic_wgs.vcf"})
    assert p.metadata.is_synthetic is True
    assert p.metadata.input_label == "synthetic_wgs.vcf"
    assert build_report_payload(None, None, None).metadata.is_synthetic is False


def test_provenance_floor_breach_is_an_error():
    p = _coherent()
    p.provenance.registry_n_parameters = 66
    p.provenance.registry_pct_sourced = 74.2
    errs = errors_in(validate_payload(p))
    assert any(e["check"] == "Registry provenance floor" for e in errs)


def test_two_provenance_denominators_are_reported_not_merged():
    p = _coherent()
    p.provenance.registry_pct_sourced = 75.4
    p.provenance.registry_n_parameters = 65
    p.provenance.model_pct_resolvable = 46.9
    p.provenance.model_n_total_known = 382
    info = _checks(p, "Two provenance denominators")
    assert info and "65" in info[0]["detail"] and "382" in info[0]["detail"]


# ── the builder against real engine output ────────────────────────────────────

def test_builder_maps_curated_keys_that_actually_exist():
    """REGRESSION. The first builder read `value` and `cost` from the personal
    sheet; the real keys are `avoided` and `intervention`, so every per-finding
    cash figure silently came through as zero and the net-cash identity check
    passed vacuously."""
    personal = {"items": [{"category": "Carrier Screening",
                           "finding": "PTPN22 R620W carrier",
                           "avoided": 84, "qaly": 0.07, "qaly_value": 5607,
                           "intervention": 65, "net": 5626,
                           "confidence": "moderate", "adherence": 0.65,
                           "pool_target": "topic:autoimmune"}]}
    p = build_report_payload(None, None, personal)
    f = p.findings[0]
    assert f.medical_cost_averted == 84.0
    assert f.intervention_cost == 65.0
    assert f.net_cash == 19.0
    assert f.monetized_qaly_value == 5607.0
    assert f.health_economic_value == 5626.0
    assert f.pool_target == "topic:autoimmune"
    assert errors_in(validate_payload(p)) == []


# ── canonicalisation: one dollar figure per pathway ───────────────────────────

def test_canonical_value_is_the_parametric_expected_nmb():
    """Registry-backed parametric NMB is canonical because it decomposes into
    explicit probabilities, effect sizes and costs, participates in the PSA, and
    responds to willingness-to-pay. The curated figure does none of those."""
    voi = {"nmb_rows": [{"label": "CYP2C19 IM (clopidogrel)",
                         "economic_pathway_id": "pgx:cyp2c19:clopidogrel:mace",
                         "confidence": "high", "nmb": 490,
                         "dcost_averted": 300, "dqaly": 0.002}]}
    personal = {"items": [{"finding": "Avoid clopidogrel non-response",
                           "economic_pathway_id": "pgx:cyp2c19:clopidogrel:mace",
                           "avoided": 23, "intervention": 100, "qaly": 0.05,
                           "qaly_value": 4529, "net": 4452,
                           "confidence": "high"}]}
    p = build_report_payload(None, voi, personal)
    assert len(p.findings) == 1, "one pathway must produce one record"
    f = p.findings[0]
    assert f.canonical_expected_nmb == 490.0
    assert f.economic_value_basis == "parametric_expected_nmb"
    assert f.legacy_curated_value == 4452.0
    assert f.legacy_curated_value_basis == "curated_prevalence_weighted_mixture"


def test_curated_value_never_carries_a_generic_economic_label():
    """It is neither a gross event value nor an expected NMB: outcome_value is
    already p_rx x p_adr x rrr x adr_cost, and the curated path then multiplies
    it — and the QALY — by genotype prevalence."""
    p = build_report_payload(
        None, {"nmb_rows": [{"label": "x", "economic_pathway_id": "pgx:a:b:c",
                             "nmb": 10}]},
        {"items": [{"finding": "y", "economic_pathway_id": "pgx:a:b:c",
                    "net": 999}]})
    basis = p.findings[0].legacy_curated_value_basis
    assert basis == "curated_prevalence_weighted_mixture"
    assert basis not in ("nmb", "net_value", "economic_value")


def test_unmapped_curated_finding_gets_null_not_an_invented_value():
    p = build_report_payload(None, None, {"items": [
        {"finding": "1 actionable wellness variant(s)", "net": 1767,
         "avoided": 362, "intervention": 105, "confidence": "low"}]})
    f = p.findings[0]
    assert f.canonical_expected_nmb is None
    assert f.economic_value_basis == "not_yet_standardised"
    assert f.legacy_curated_value == 1767.0
    assert f.display_name, "the genomic finding must still render"


def test_findings_sort_by_canonical_with_nulls_last():
    p = build_report_payload(None, {"nmb_rows": [
        {"label": "lo", "economic_pathway_id": "pgx:a:b:c", "nmb": 90,
         "confidence": "high"},
        {"label": "hi", "economic_pathway_id": "pgx:d:e:f", "nmb": 490,
         "confidence": "high"}]},
        {"items": [{"finding": "unstd", "net": 5000, "confidence": "high"}]})
    ordered = [f for _g, items in p.findings_page_groups() for f in items]
    vals = [f.canonical_expected_nmb for f in ordered]
    assert vals[0] == 490.0 and vals[1] == 90.0
    assert vals[-1] is None, "unstandardised findings sort last, never dropped"


def test_clinical_grouping_outranks_dollar_value():
    """A low-confidence finding must not outrank a high-confidence prescribing
    finding on the strength of a modelled dollar estimate."""
    p = build_report_payload(None, {"nmb_rows": [
        {"label": "CYP2C19 IM", "economic_pathway_id": "pgx:cyp2c19:c:mace",
         "nmb": 90, "confidence": "high"},
        {"label": "wellness thing", "economic_pathway_id": "well:x",
         "nmb": 9000, "confidence": "low"}]})
    groups = dict(p.findings_page_groups())
    assert "CYP2C19 IM" in [f.display_name
                            for f in groups["Medication & prescribing"]]
    assert "wellness thing" in [
        f.display_name for f in groups["Lower-confidence & exploratory"]]
    first_group = p.findings_page_groups()[0][0]
    assert first_group == "Medication & prescribing"


def test_mthfr_stays_negative_under_current_parametric_assumptions():
    p = build_report_payload(None, {"nmb_rows": [
        {"label": "Folate Metabolism (MTHFR)",
         "economic_pathway_id": "legacy:mthfr:folate", "nmb": -216,
         "dcost_averted": 219, "confidence": "low"}]}, None)
    assert p.findings[0].canonical_expected_nmb == -216.0


def test_pooled_reference_case_is_untouched_by_canonicalisation():
    voi = {"pooled_economics": {"cea": {"inmb": 7410, "incremental_cost": -4223,
                                        "incremental_qaly": 0.0319,
                                        "wtp": 100000}},
           "nmb_rows": [{"label": "a", "economic_pathway_id": "pgx:a:b:c",
                         "nmb": 490}]}
    p = build_report_payload(None, voi, {"items": [
        {"finding": "a", "economic_pathway_id": "pgx:a:b:c", "net": 4452}]})
    assert p.reference_case.nmb == 7410.0


def test_per_finding_canonical_values_do_not_sum_to_the_reference_case():
    """Standalone expected values are not additive; the reference case pools
    overlapping signals first. No pro-rata allocation is invented to force
    agreement."""
    voi = {"pooled_economics": {"cea": {"inmb": 7410, "wtp": 100000,
                                        "incremental_qaly": 0.0319,
                                        "incremental_cost": -4223}},
           "nmb_rows": [{"label": "a", "economic_pathway_id": "p:a", "nmb": 3514},
                        {"label": "b", "economic_pathway_id": "p:b", "nmb": 2509}]}
    p = build_report_payload(None, voi, None)
    assert sum(f.canonical_expected_nmb or 0 for f in p.findings) != p.reference_case.nmb


# ── identifiers ───────────────────────────────────────────────────────────────

def test_pathway_ids_are_semantic_not_slugified_display_text():
    from econ import identity as ident
    pid = ident.economic_pathway_id(kind="pgx", gene="CYP2C19",
                                    drug="clopidogrel", condition="mace")
    assert pid == "pgx:cyp2c19:clopidogrel:mace"
    assert ident.economic_pathway_id(
        kind="pgx", gene="SLCO1B1", drug="simvastatin",
        condition="myopathy") == "pgx:slco1b1:simvastatin:myopathy"


def test_display_wording_changes_do_not_move_a_semantic_id():
    """The whole point of not slugifying display text."""
    from econ import identity as ident
    a = ident.economic_pathway_id(kind="pgx", gene="CYP2C19",
                                  drug="clopidogrel", condition="mace")
    b = ident.economic_pathway_id(kind="pgx", gene="CYP2C19",
                                  drug="clopidogrel", condition="mace")
    assert a == b
    assert ident.finding_id(kind="pgx", gene="CYP2C19",
                            phenotype="Intermediate Metabolizer (IM)") == \
           ident.finding_id(kind="pgx", gene="CYP2C19",
                            phenotype="intermediate metabolizer")


def test_legacy_ids_are_flagged_as_such():
    from econ import identity as ident
    legacy = ident.legacy_pathway_id("Some finding text", "Category")
    assert ident.is_legacy(legacy)
    assert not ident.is_legacy("pgx:cyp2c19:clopidogrel:mace")


def test_pgx_extractor_emits_semantic_not_legacy_ids():
    from econ.health_economics import _pgx_findings
    rows = _pgx_findings({"CYP2C19": {"phenotype": "Intermediate Metabolizer (IM)"}})
    assert rows[0]["economic_pathway_id"] == "pgx:cyp2c19:clopidogrel:mace"
    assert rows[0]["pathway_id_is_legacy"] is False


# ── confidence vocabulary ─────────────────────────────────────────────────────

def test_unrecognised_confidence_is_normalised_for_display_and_flagged():
    """The novel-variants module emits "higher" for its most confident
    predictions; _CONFIDENCE_CONC has no such key, so the PSA silently gives it
    the moderate concentration."""
    p = build_report_payload(None, {"nmb_rows": [
        {"label": "predicted-pathogenic 17:7674220", "economic_pathway_id": "v:1",
         "nmb": 2341, "confidence": "higher"}]}, None)
    f = p.findings[0]
    assert f.evidence_confidence == "high"
    assert f.raw_evidence_confidence == "higher"
    assert Severity.WARNING in _severities(p, "Unrecognised confidence vocabulary")


# ── budget metric rename ──────────────────────────────────────────────────────

def test_budget_validator_reads_the_explicit_name():
    p = _coherent()
    p.advanced.budget_impact = {
        "available": True, "rows": _budget_rows(),
        "maximum_budget_burden_pmpm": 0.1167, "maximum_budget_burden_year": 1}
    assert errors_in(validate_payload(p)) == []


def test_budget_validator_still_accepts_the_legacy_name():
    p = _coherent()
    p.advanced.budget_impact = {"available": True, "rows": _budget_rows(),
                                "peak_pmpm": 0.1167, "peak_year": 1}
    assert errors_in(validate_payload(p)) == []
