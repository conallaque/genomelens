"""Correctness guarantees for the decision-analytic layer.

These outputs are easy to compute and easy to get subtly wrong in ways that
still look plausible — an EVPPI above its own EVPI, an efficiency frontier that
keeps an unchoosable option, a budget impact that credits decades of prevention
savings against a five-year window. Each of those happened while this module
was being written; the tests below are what caught them.
"""

import pytest

import econ_decision as ed
import econ_engine as ee
import econ_params as ep


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

def test_budget_impact_does_not_credit_a_lifetime_of_savings_to_five_years():
    # An early version credited the full modelled offset inside the budget
    # window and reported a programme saving $65m, which would have been read
    # as self-financing. Prevention savings accrue over decades.
    bi = ed.budget_impact(per_person_cost=2_100, per_person_offset=10_795,
                          population=10_000)
    assert bi["offset_realised_in_horizon"] < 0.5
    assert bi["total_net"] > 0, (
        "a prevention programme should cost a payer money over five years")


def test_budget_impact_is_undiscounted_and_cumulative():
    bi = ed.budget_impact(per_person_cost=100, per_person_offset=0,
                          population=1_000, uptake=(0.1, 0.1), years=2)
    rows = bi["rows"]
    assert rows[0]["net_budget_impact"] == rows[1]["net_budget_impact"], (
        "equal uptake in equal years must give equal nominal impact — any "
        "difference means discounting crept into a budget-impact analysis")
    assert rows[1]["cumulative"] == rows[0]["cumulative"] * 2


def test_budget_impact_scales_with_population():
    small = ed.budget_impact(per_person_cost=100, per_person_offset=10,
                             population=1_000)
    big = ed.budget_impact(per_person_cost=100, per_person_offset=10,
                           population=10_000)
    assert big["total_net"] == pytest.approx(small["total_net"] * 10, rel=1e-6)


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
    import value_of_information as voi
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
    for fn in (ed.budget_impact, ed.distributional_cea, ed.subgroup_analysis):
        for name, prm in inspect.signature(fn).parameters.items():
            if name in ("offset_realised_in_horizon", "inequality_aversion",
                        "annual_incidence"):
                assert prm.default is None, (
                    f"{fn.__name__}({name}=...) still carries a literal "
                    f"default; it should read from econ_params so it appears "
                    f"in the provenance count and the sensitivity analysis")


def test_registered_decision_parameters_drive_the_functions():
    with ep.overridden({"budget_offset_realised_in_horizon": 0.9}):
        generous = ed.budget_impact(per_person_cost=100,
                                    per_person_offset=1_000, population=1_000)
    with ep.overridden({"budget_offset_realised_in_horizon": 0.01}):
        stingy = ed.budget_impact(per_person_cost=100,
                                  per_person_offset=1_000, population=1_000)
    assert generous["total_net"] < stingy["total_net"]


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
