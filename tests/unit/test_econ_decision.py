"""Correctness guarantees for the decision-analytic layer.

These outputs are easy to compute and easy to get subtly wrong in ways that
still look plausible — an EVPPI above its own EVPI, an efficiency frontier that
keeps an unchoosable option, a budget impact that credits decades of prevention
savings against a five-year window. Each of those happened while this module
was being written; the tests below are what caught them.
"""

import pytest

from econ import decision as ed
from econ import engine as ee
from econ import params as ep


def _build(conditions=(("CAD", 1.0), ("T2D", 1.0), ("Alzheimer", 1.0))):
    def build():
        return ee.pool_findings([
            ee.Finding(label=f"f{i}", coi_key=k,
                       p_event=ep.value("baseline_event_probability"),
                       rrr=ep.value("actionable_rrr"), haircut=h,
                       intervention_cost=ep.value("intervention_cost_standard"))
            for i, (k, h) in enumerate(conditions)])
    return build


def _marginal():
    """A finding set where the decision genuinely straddles zero."""
    def build():
        return ee.pool_findings([
            ee.Finding(label="f", coi_key="Urologic", p_event=0.05,
                       rrr=ep.value("actionable_rrr"), haircut=0.4,
                       intervention_cost=ep.value("intervention_cost_standard"))])
    return build


# ══════════════════════════════════════════════════════════════════════════
# Value of information
# ══════════════════════════════════════════════════════════════════════════

def test_evpi_is_non_negative_and_reports_its_own_iteration_count():
    r = ed.evpi(_build(), n=300, test_cost=100)
    assert r["evpi_per_person"] >= 0
    assert r["n_iterations"] == 300


def test_evpi_is_larger_when_the_decision_actually_straddles():
    # EVPI measures decision regret, not spread. A confident recommendation
    # should carry near-zero EVPI however wide its interval; a marginal one
    # should carry more.
    confident = ed.evpi(_build(), n=1200, test_cost=100)
    marginal = ed.evpi(_marginal(), n=1200, test_cost=300)
    assert marginal["p_current_choice_wrong"] > confident["p_current_choice_wrong"]
    assert marginal["evpi_per_person"] >= confident["evpi_per_person"]


def test_evppi_never_exceeds_evpi():
    # THE BUG THIS CATCHES. Estimating the unconditional mean from a separate
    # sample let Monte Carlo noise alone produce a positive EVPPI — an early
    # version reported $2,299 of partial information value against an EVPI of
    # $0. Partial information cannot be worth more than perfect information.
    build = _marginal()
    total = ed.evpi(build, n=1500, test_cost=300)["evpi_per_person"]
    rows = ed.evppi(build, parameters=["actionable_rrr",
                                       "baseline_event_probability",
                                       "coi_urologic"],
                    n_outer=40, n_inner=20, test_cost=300)
    assert rows
    for r in rows:
        assert r["evppi_per_person"] <= total + max(5, 0.25 * total), (
            f"EVPPI for {r['parameter']} is {r['evppi_per_person']} against "
            f"an EVPI of {total} — partial information cannot beat perfect")


def test_evppi_is_non_negative_and_sorted():
    rows = ed.evppi(_build(), parameters=["actionable_rrr", "coi_mace"],
                    n_outer=25, n_inner=12, test_cost=100)
    assert all(r["evppi_per_person"] >= 0 for r in rows)
    vals = [r["evppi_per_person"] for r in rows]
    assert vals == sorted(vals, reverse=True)


def test_evppi_ignores_unknown_parameters_rather_than_crashing():
    rows = ed.evppi(_build(), parameters=["not_a_parameter", "actionable_rrr"],
                    n_outer=10, n_inner=5, test_cost=100)
    assert [r["parameter"] for r in rows] == ["actionable_rrr"]


def test_population_evpi_scales_and_discounts():
    one = ed.population_evpi(100.0, population=100_000, years=10)
    ten = ed.population_evpi(1000.0, population=100_000, years=10)
    assert ten["population_evpi"] == pytest.approx(one["population_evpi"] * 10,
                                                   rel=1e-6)
    undisc = ed.population_evpi(100.0, population=100_000, years=10, rate=0.0)
    assert undisc["population_evpi"] > one["population_evpi"]


def test_zero_evpi_scales_to_zero_population_value():
    assert ed.population_evpi(0.0)["population_evpi"] == 0


# ══════════════════════════════════════════════════════════════════════════
# Breakeven
# ══════════════════════════════════════════════════════════════════════════

def test_breakeven_reports_no_crossing_when_the_sign_never_changes():
    r = ed.breakeven(_build(), parameter="actionable_rrr", test_cost=100)
    assert r["available"]
    if not r["crosses_within_range"]:
        assert r["breakeven_value"] is None
        assert "keeps its sign" in r["interpretation"]


def test_breakeven_finds_the_crossing_when_one_exists():
    # Push the test cost high enough that low effect sizes stop being worth it.
    r = ed.breakeven(_marginal(), parameter="actionable_rrr", test_cost=3_000)
    if r["crosses_within_range"]:
        lo, hi = r["range"]
        assert lo <= r["breakeven_value"] <= hi
        assert r["margin"] is not None


def test_breakeven_rejects_unknown_parameters():
    r = ed.breakeven(_build(), parameter="nope")
    assert r["available"] is False


# ══════════════════════════════════════════════════════════════════════════
# Efficiency frontier
# ══════════════════════════════════════════════════════════════════════════

def test_strictly_dominated_strategies_are_marked():
    f = ed.efficiency_frontier([
        {"name": "A", "cost": 0, "qaly": 10.0},
        {"name": "B", "cost": 1_000, "qaly": 10.1},
        {"name": "Worse", "cost": 2_000, "qaly": 10.05},
    ])
    by = {s["name"]: s for s in f["strategies"]}
    assert by["Worse"]["status"] == "dominated"


def test_extended_dominance_is_applied():
    # The step usually skipped. B is cheaper and less effective than C, and not
    # strictly dominated, but a mix of A and C beats it — so no decision-maker
    # should ever choose B, and leaving it on the frontier also corrupts C's ICER.
    f = ed.efficiency_frontier([
        {"name": "A", "cost": 0, "qaly": 10.00},
        {"name": "B", "cost": 1_400, "qaly": 10.13},
        {"name": "C", "cost": 1_800, "qaly": 10.20},
    ])
    by = {s["name"]: s for s in f["strategies"]}
    assert by["B"]["status"] == "extendedly dominated"
    assert by["C"]["status"] == "on frontier"


def test_frontier_icers_increase_along_the_frontier():
    f = ed.efficiency_frontier([
        {"name": "A", "cost": 0, "qaly": 10.0},
        {"name": "B", "cost": 500, "qaly": 10.10},
        {"name": "C", "cost": 2_000, "qaly": 10.20},
        {"name": "D", "cost": 6_000, "qaly": 10.25},
    ])
    icers = [s["icer"] for s in f["strategies"]
             if s["status"] == "on frontier" and s["icer"] is not None]
    assert icers == sorted(icers), (
        f"ICERs along the efficiency frontier must increase; got {icers}")


def test_frontier_recommends_the_highest_net_benefit_option():
    f = ed.efficiency_frontier([
        {"name": "Cheap", "cost": 0, "qaly": 10.0},
        {"name": "Best", "cost": 1_000, "qaly": 10.5},
    ], wtp=100_000)
    assert f["recommended"] == "Best"
    low = ed.efficiency_frontier([
        {"name": "Cheap", "cost": 0, "qaly": 10.0},
        {"name": "Best", "cost": 1_000, "qaly": 10.005},
    ], wtp=10_000)
    assert low["recommended"] == "Cheap"


# ══════════════════════════════════════════════════════════════════════════
# Budget impact
# ══════════════════════════════════════════════════════════════════════════

def test_budget_impact_delegates_rather_than_duplicating():
    # This module briefly carried its own budget-impact implementation
    # alongside the one in markov_model, so the report contained two
    # affordability analyses with different conventions and different answers.
    # One question, one calculation.
    from econ import markov as mk
    mine = ed.budget_impact(per_person_cost=300, population=1_000_000)
    theirs = mk.budget_impact(plan_members=1_000_000, test_cost=300)
    assert mine["rows"] == theirs["rows"]
    assert mine["peak_pmpm"] == theirs["peak_pmpm"]


def test_budget_impact_reports_per_member_per_month():
    # PMPM is the metric coverage decisions are actually argued on, and the
    # duplicate implementation did not produce it.
    bi = ed.budget_impact(per_person_cost=300, population=1_000_000)
    assert "peak_pmpm" in bi and bi["peak_pmpm"] > 0
    assert all("pmpm" in r for r in bi["rows"])


def test_budget_impact_is_undiscounted():
    # BIA answers affordability in nominal cash. Discounting it would silently
    # turn it into a small, wrong cost-effectiveness analysis.
    bi = ed.budget_impact(per_person_cost=300, population=1_000_000)
    for r in bi["rows"]:
        assert r["net_budget_impact"] == (
            r["cost_testing"] + r["cost_intervention"] - r["offsets"])


def test_budget_impact_scales_with_plan_size():
    small = ed.budget_impact(per_person_cost=300, population=100_000)
    big = ed.budget_impact(per_person_cost=300, population=1_000_000)
    assert big["cumulative_net"] == pytest.approx(small["cumulative_net"] * 10,
                                                  rel=1e-6)


def test_prevention_costs_a_payer_money_before_it_saves_any():
    # Offsets ramp as the intervention takes effect, so year one is pure spend.
    bi = ed.budget_impact(per_person_cost=300, population=1_000_000)
    assert bi["rows"][0]["offsets"] == 0
    assert bi["rows"][0]["net_budget_impact"] > 0


# ══════════════════════════════════════════════════════════════════════════
# Heterogeneity and equity
# ══════════════════════════════════════════════════════════════════════════

def test_benefit_falls_with_age_through_competing_mortality():
    s = ed.subgroup_analysis(ages=(40, 60, 80), sexes=("Female",))
    inmbs = [r["inmb"] for r in sorted(s["rows"], key=lambda r: r["age"])]
    assert inmbs == sorted(inmbs, reverse=True), (
        f"net benefit should fall with age as competing mortality rises: {inmbs}")


def test_subgroup_analysis_covers_both_sexes_and_reports_the_spread():
    s = ed.subgroup_analysis(ages=(40, 70), sexes=("Female", "Male"))
    assert {r["sex"] for r in s["rows"]} == {"Female", "Male"}
    assert s["spread"] > 0


def test_distributional_analysis_detects_a_pro_poor_allocation():
    # Health gains concentrated on the worst-off group must register as
    # inequality-reducing; the opposite allocation must not.
    base = {"worst off": 5.0, "middle": 10.0, "best off": 15.0}
    pro_poor = ed.distributional_cea(base, {"worst off": 2.0, "middle": 0.0,
                                            "best off": 0.0})
    pro_rich = ed.distributional_cea(base, {"worst off": 0.0, "middle": 0.0,
                                            "best off": 2.0})
    assert pro_poor["reduces_inequality"] is True
    assert pro_rich["reduces_inequality"] is False
    assert pro_poor["ede_gain"] > pro_rich["ede_gain"], (
        "with inequality aversion, a gain to the worst-off must be valued "
        "above an identical gain to the best-off")


def test_equally_distributed_equivalent_never_exceeds_the_mean():
    # Jensen's inequality: the EDE of an unequal distribution is below its mean.
    d = ed.distributional_cea({"a": 4.0, "b": 10.0, "c": 20.0},
                              {"a": 1.0, "b": 1.0, "c": 1.0})
    assert d["ede_before"] <= d["mean_qaly_before"] + 1e-9
    assert d["ede_after"] <= d["mean_qaly_after"] + 1e-9


def test_perfectly_equal_distribution_has_ede_equal_to_mean():
    d = ed.distributional_cea({"a": 10.0, "b": 10.0}, {"a": 0.0, "b": 0.0})
    assert d["ede_before"] == pytest.approx(d["mean_qaly_before"], rel=1e-6)
    assert d["gini_before"] == pytest.approx(0.0, abs=1e-6)


def test_distributional_analysis_degrades_on_empty_input():
    assert ed.distributional_cea({}, {})["available"] is False


# ══════════════════════════════════════════════════════════════════════════
# Orchestration
# ══════════════════════════════════════════════════════════════════════════

def test_decision_layer_assembles_every_section():
    d = ed.analyze_decision_layer(_build(), test_cost=100, fast=True)
    for key in ("evpi", "population_evpi", "evppi", "breakeven",
                "subgroups", "budget_impact", "distributional"):
        assert key in d, f"decision layer missing {key}"


def test_decision_layer_restores_the_registry():
    before = {k: p.value for k, p in ep.PARAMS.items()}
    ed.analyze_decision_layer(_build(), test_cost=100, fast=True)
    assert {k: p.value for k, p in ep.PARAMS.items()} == before


def test_decision_layer_prioritises_the_parameters_the_tornado_names():
    tor = [{"parameter": "coi_mace"}, {"parameter": "actionable_rrr"}]
    d = ed.analyze_decision_layer(_build(), tornado_rows=tor, test_cost=100,
                                  fast=True)
    assert {r["parameter"] for r in d["evppi"]} <= {"coi_mace", "actionable_rrr"}


# ══════════════════════════════════════════════════════════════════════════
# Integration: what value_of_information actually feeds the decision layer
# ══════════════════════════════════════════════════════════════════════════
# The unit tests above exercise these functions with synthetic inputs. They
# cannot see the wiring, which is exactly where an invented multiplier and a
# gross-vs-net cost mismatch survived several passes.

def _voi(input_type="chip", with_wgs_findings=False):
    from econ import value_of_information as voi
    econ = {"findings_with_economics": [
        {"finding": "CAD polygenic score", "category": "Polygenic Risk",
         "qaly_gain": 0.5, "confidence": "moderate"},
        {"finding": "CYP2C19 IM", "category": "Pharmacogenomics",
         "qaly_gain": 0.3, "confidence": "high"},
    ]}
    cvr = None
    if with_wgs_findings:
        cvr = {"available": True, "buckets": {
            "actionable": [{"gene": "BRCA2"}],
            "carrier": [{"gene": "CFTR"}]}}
    return voi.analyze_value_of_information(
        economics_result=econ, clinical_variants_result=cvr,
        input_type=input_type, n_mc=200, seed=5, log=lambda *a: None)


def test_frontier_costs_match_the_headline_incremental_cost():
    # THE BUG THIS CATCHES. The frontier was built on gross cost while the
    # CEA card beside it reported net, so the same strategy appeared as
    # +$2,100 in one table and -$6,440 in the other — same box, same
    # question, opposite signs.
    p = _voi()["pooled_economics"]
    fr = p["decision"]["frontier"]
    chip = next(s for s in fr["strategies"] if s["name"] == "Genotyping chip")
    assert chip["cost"] == pytest.approx(p["cea"]["incremental_cost"], abs=2), (
        f"frontier says {chip['cost']}, CEA card says "
        f"{p['cea']['incremental_cost']} — the two must share a cost basis")
    assert chip["qaly"] == pytest.approx(p["cea"]["incremental_qaly"], abs=1e-3)


def test_no_sequencing_arm_is_invented_from_chip_input():
    # The WGS arm was previously manufactured by multiplying the chip result
    # by 1.35 QALYs and 1.6 cost — two numbers typed into wiring code that
    # between them decided which strategy the report recommended.
    fr = _voi(input_type="chip")["pooled_economics"]["decision"]["frontier"]
    names = {s["name"] for s in fr["strategies"]}
    assert "Whole-genome sequencing" not in names, (
        "a sequencing arm cannot be estimated from chip input — only asserted")
    assert fr.get("wgs_not_estimable"), "the omission must be explained"


def test_sequencing_arm_appears_only_when_sequencing_findings_do():
    fr = _voi(with_wgs_findings=True)["pooled_economics"]["decision"]["frontier"]
    names = {s["name"] for s in fr["strategies"]}
    assert "Whole-genome sequencing" in names
    assert not fr.get("wgs_not_estimable")


def test_every_decision_layer_constant_is_registered():
    # A default argument that decides the sign of a reported result is a
    # parameter, and belongs in the registry with the rest of them.
    import inspect
    for fn in (ed.distributional_cea, ed.subgroup_analysis):
        for name, prm in inspect.signature(fn).parameters.items():
            if name in ("offset_realised_in_horizon", "inequality_aversion",
                        "annual_incidence"):
                assert prm.default is None, (
                    f"{fn.__name__}({name}=...) still carries a literal "
                    f"default; it should read from econ_params so it appears "
                    f"in the provenance count and the sensitivity analysis")


def test_registered_decision_parameters_drive_the_functions():
    with ep.overridden({"inequality_aversion": 0.0}):
        neutral = ed.distributional_cea({"a": 5.0, "b": 15.0},
                                        {"a": 1.0, "b": 0.0})
    with ep.overridden({"inequality_aversion": 20.0}):
        averse = ed.distributional_cea({"a": 5.0, "b": 15.0},
                                       {"a": 1.0, "b": 0.0})
    assert averse["ede_gain"] != neutral["ede_gain"], (
        "the registered inequality-aversion parameter must reach the "
        "calculation, or the provenance table is decorative")


def test_subgroup_table_follows_the_dominant_condition():
    d = ed.analyze_decision_layer(_build((("Alzheimer", 1.0),)),
                                  test_cost=100, fast=True)
    assert d["subgroups"]["condition"] == "Alzheimer", (
        "the subgroup table should describe the person's own dominant "
        "condition, not always CAD")


def test_subgroup_table_declares_itself_illustrative():
    s = ed.subgroup_analysis(ages=(40, 60), sexes=("Female",))
    assert s.get("illustrative") is True
    assert "not a personalised risk estimate" in s["note"]


# ══════════════════════════════════════════════════════════════════════════
# Plain-language layer
# ══════════════════════════════════════════════════════════════════════════
# A translation layer is where over-claiming creeps in: it is easy to write a
# confident English sentence on top of a hedged technical result. These tests
# check the translation stays as honest as the thing it translates.

from econ import plain as epl


def test_number_needed_to_screen_inverts_absolute_risk_reduction():
    r = epl.number_needed_to_screen([{"condition": "CAD",
                                      "cases_averted": 0.02}])
    assert r[0]["number_needed_to_screen"] == 50


def test_plain_language_uses_the_conditional_when_the_baseline_is_assumed():
    # The baseline risk behind these numbers is a registered assumption for
    # most findings. Stating "1 in 17 people avoid a serious inherited
    # condition" in the indicative turns that assumption into a claim — which
    # would make the plain-language layer less honest than the technical
    # output it summarises.
    r = epl.number_needed_to_screen([{"condition": "Pathogenic",
                                      "cases_averted": 0.06}])
    assert r[0]["baseline_is_assumed"] is True
    assert "if the model's assumptions hold" in r[0]["plain"].lower()


def test_plain_text_contains_no_markdown_or_jargon():
    pooled = {"available": True,
              "cea": {"incremental_cost": -6_440, "incremental_qaly": 0.069,
                      "icer": None, "inmb": 13_340, "wtp": 100_000,
                      "cost_averted": 8_540, "intervention_cost": 2_000,
                      "test_cost": 100, "horizon_years": 10},
              "conditions": [{"condition": "CAD", "cases_averted": 0.02,
                              "combined_rrr": 0.3, "n_findings": 2,
                              "inmb": 5_000}],
              "psa": {"available": True, "p_cost_effective": 0.98,
                      "n_iterations": 1500, "inmb_ci_low": 300,
                      "inmb_ci_high": 39_000, "note": "n"},
              "tornado": [{"parameter": "actionable_rrr", "tier": "assumption",
                           "swing": 100}]}
    s = epl.build_plain_summary(pooled)
    blob = " ".join([s["bottom_line"], s["verdict"]["plain"],
                     s["money"]["plain"], s["what_would_change_it"],
                     s["healthy_time"]["plain"]])
    assert "*" not in blob, "markdown emphasis does not render in HTML"
    for jargon in ("QALY", "INMB", "ICER", "incremental net monetary",
                   "coi_key", "EVPPI"):
        assert jargon not in blob, f"{jargon!r} leaked into plain-language text"


def test_verdict_distinguishes_cost_saving_from_merely_cost_effective():
    # Conflating these is the most common way a report like this misleads:
    # prevention usually improves health while COSTING money.
    saving = epl.plain_verdict({"incremental_cost": -5_000,
                                "incremental_qaly": 0.5, "icer": None})
    buying = epl.plain_verdict({"incremental_cost": 20_000,
                                "incremental_qaly": 0.5, "icer": 40_000})
    assert "saves money" in saving["headline"]
    assert "costs money" in buying["headline"]
    assert "not a saving" in buying["plain"]


def test_verdict_rejects_a_strategy_above_the_threshold():
    bad = epl.plain_verdict({"incremental_cost": 200_000,
                             "incremental_qaly": 0.5, "icer": 400_000},
                            wtp=100_000)
    assert bad["tone"] == "negative"


def test_healthy_time_scales_units_to_the_magnitude():
    assert "days" in epl.healthy_time_gained(0.05)["plain"]
    assert "months" in epl.healthy_time_gained(0.5)["plain"]
    assert "years" in epl.healthy_time_gained(4.0)["plain"]
    assert "rounding error" in epl.healthy_time_gained(0.001)["plain"]


def test_healthy_time_carries_the_averaging_caveat():
    # "25 extra days of healthy life" invites being read as a personal
    # promise. Most people get nothing from any one preventive action.
    assert "not a promise" in epl.healthy_time_gained(0.07)["caveat"]


def test_payback_is_honest_when_there_is_none():
    r = epl.payback_period(upfront_cost=2_000, annual_saving=0)
    assert r["pays_back"] is False
    assert "buying health, not savings" in r["plain"]
    slow = epl.payback_period(upfront_cost=2_000, annual_saving=10)
    assert slow["pays_back"] is False


def test_payback_reports_a_real_recovery_period():
    r = epl.payback_period(upfront_cost=2_100, annual_saving=854)
    assert r["pays_back"] is True
    assert r["years"] == pytest.approx(2.5, abs=0.1)


def test_confidence_is_expressed_as_a_count_out_of_a_hundred():
    c = epl.plain_confidence({"available": True, "p_cost_effective": 0.83,
                              "n_iterations": 1000, "inmb_ci_low": -100,
                              "inmb_ci_high": 5_000, "note": ""})
    assert c["n_in_100"] == 83
    assert "83 of every 100" in c["plain"]
    assert "not all" in c["plain"]


def test_confidence_says_so_when_the_result_is_a_coin_flip():
    c = epl.plain_confidence({"available": True, "p_cost_effective": 0.52,
                              "n_iterations": 1000, "inmb_ci_low": -5_000,
                              "inmb_ci_high": 5_000, "note": ""})
    assert "could go either way" in c["plain"]


def test_summary_warns_when_assumptions_drive_the_answer():
    pooled = {"available": True, "cea": {"incremental_qaly": 0.1, "wtp": 100_000},
              "conditions": [], "psa": {},
              "tornado": [{"parameter": "actionable_rrr",
                           "tier": "assumption", "swing": 900},
                          {"parameter": "coi_mace", "tier": "derived",
                           "swing": 100}]}
    s = epl.build_plain_summary(pooled)
    assert s["assumption_share"] >= 50
    assert "guesswork" in s["what_would_change_it"]


def test_summary_does_not_cry_wolf_when_evidence_dominates():
    pooled = {"available": True, "cea": {"incremental_qaly": 0.1, "wtp": 100_000},
              "conditions": [], "psa": {},
              "tornado": [{"parameter": "coi_mace", "tier": "derived",
                           "swing": 900},
                          {"parameter": "actionable_rrr", "tier": "assumption",
                           "swing": 50}]}
    s = epl.build_plain_summary(pooled)
    assert "published evidence" in s["what_would_change_it"]


def test_condition_codes_are_translated_to_english():
    for key in ("CAD", "T2D", "Alzheimer", "Pathogenic"):
        assert key not in epl.CONDITION_NAMES[key]
    assert epl._name("CAD") == "heart attack or stroke"


def test_summary_carries_a_disclaimer():
    pooled = {"available": True, "cea": {"incremental_qaly": 0.1},
              "conditions": [], "psa": {}, "tornado": []}
    s = epl.build_plain_summary(pooled)
    assert "not medical advice" in s["disclaimer"]


def test_summary_degrades_when_there_is_nothing_to_summarise():
    assert epl.build_plain_summary(None)["available"] is False
    assert epl.build_plain_summary({"available": False})["available"] is False


# ══════════════════════════════════════════════════════════════════════════
# Is sequencing worth buying?
# ══════════════════════════════════════════════════════════════════════════

def test_prospective_sequencing_value_is_not_structurally_zero():
    # THE BUG THIS ADDRESSES. `marginal_chip_to_wgs` sums findings tagged
    # wgs_only, which only exist in a VCF — so on array input it is $0 for
    # everybody, every time, and reads as "sequencing would add nothing".
    # A prospective estimate must not inherit that property.
    w = ed.wgs_marginal_value(n_wgs_only_findings=0)
    assert w["available"]
    assert w["gross_expected_value"] > 0, (
        "the prospective estimate must not be zero just because the input "
        "file was an array")
    assert w["number_needed_to_sequence"] and w["number_needed_to_sequence"] > 1


def test_the_zero_is_explained_when_there_are_no_sequencing_findings():
    w = ed.wgs_marginal_value(n_wgs_only_findings=0, wgs_only_findings_value=0)
    assert "statement about the input file" in w["why_retrospective_is_zero"]


def test_no_spurious_explanation_when_sequencing_findings_exist():
    w = ed.wgs_marginal_value(n_wgs_only_findings=3,
                              wgs_only_findings_value=12_000)
    assert w["why_retrospective_is_zero"] == ""
    assert w["retrospective_value"] == 12_000


def test_array_blindness_to_monogenic_disease_drives_the_value():
    # If arrays detected monogenic variants well, sequencing would add little.
    # This is the mechanism, and it should be visible in the arithmetic.
    with ep.overridden({"chip_detection_share_monogenic": 0.95}):
        good_chip = ed.wgs_marginal_value()
    with ep.overridden({"chip_detection_share_monogenic": 0.05}):
        blind_chip = ed.wgs_marginal_value()
    assert blind_chip["secondary_findings_value"] > \
        good_chip["secondary_findings_value"] * 5


def test_pharmacogenomics_does_not_dominate_the_prospective_estimate():
    # PGx has by far the highest prevalence of any category. Counting it
    # undiscounted would make sequencing look far more valuable than it is,
    # because consumer arrays already type the main CPIC star alleles.
    w = ed.wgs_marginal_value()
    assert w["pgx_incremental_value"] < w["secondary_findings_value"], (
        "pharmacogenomics is mostly already covered by an array and must not "
        "dominate the case for sequencing")


def test_higher_yield_raises_the_value_and_lowers_the_number_needed():
    low = ed.wgs_marginal_value()
    with ep.overridden({"wgs_yield_acmg_secondary": 0.04}):
        high = ed.wgs_marginal_value()
    assert high["gross_expected_value"] > low["gross_expected_value"]
    assert high["number_needed_to_sequence"] < low["number_needed_to_sequence"]


def test_sequencing_stops_being_worth_it_if_it_costs_enough():
    cheap = ed.wgs_marginal_value(chip_cost=100, wgs_cost=300)
    dear = ed.wgs_marginal_value(chip_cost=100, wgs_cost=5_000)
    assert cheap["worth_it"] is True
    assert dear["worth_it"] is False
    assert "is not worth it" in dear["plain"]


def test_yield_parameters_are_registered_with_provenance():
    for key in ("wgs_yield_acmg_secondary", "chip_detection_share_monogenic",
                "chip_pgx_coverage", "wgs_yield_carrier_expanded"):
        p = ep.get(key)
        assert p.tier in ep.TIERS
        if p.tier != "assumption":
            assert p.citation, f"{key} needs a citation"


def test_prospective_estimate_warns_that_the_average_hides_the_distribution():
    # An expected value is the right basis for a policy decision and a poor
    # guide to a personal one. Saying so is the difference between an honest
    # number and a misleading one.
    w = ed.wgs_marginal_value()
    assert "skewed" in w["caveat"]
    assert "policy decision" in w["caveat"]


def test_penetrance_used_is_ascertainment_corrected():
    # Using family-series penetrance here would inflate the value of every
    # secondary finding, which is the error this project already fixed once.
    w = ed.wgs_marginal_value()
    raw = (ep.value("brca2_penetrance_population")
           * ep.value("coi_pathogenic_generic")
           * ep.value("marginal_cost_fraction") * 0.45)
    assert w["value_per_finding"] < raw * 2, (
        "value per finding looks too high for a shrunk penetrance")
