"""Correctness guarantees for the pooled cost-effectiveness engine.

The engine exists to fix one specific defect: eight of the report's finding
sources route onto the same cardiometabolic anchor, and the previous model
added their risk reductions together. These tests pin the fix, and pin the
arithmetic conventions (within-cycle correction, discounting, ICER quadrants)
against the published definitions rather than against the engine's own output.
"""

import math
import os
import sys

import pytest

import econ_engine as ee
import econ_params as ep


def _f(label="f", coi="CAD", p=0.20, rrr=0.30, hc=1.0, cost=500.0, src=""):
    return ee.Finding(label=label, coi_key=coi, p_event=p, rrr=rrr,
                      haircut=hc, intervention_cost=cost, source_category=src)


# ══════════════════════════════════════════════════════════════════════════
# The double-counting fix
# ══════════════════════════════════════════════════════════════════════════

def test_stacked_findings_do_not_sum_their_risk_reductions():
    # THE BUG THIS MODULE EXISTS FOR. Eight findings routed to CAD used to
    # produce eight independent risk reductions, added. Summed they exceed
    # any achievable effect; combined they must not.
    pool = ee.ConditionPool("CAD", [_f(label=f"f{i}") for i in range(8)])
    assert pool.naive_sum_rrr() == pytest.approx(2.4)      # 8 x 0.30, absurd
    assert pool.combined_rrr() <= ep.value("max_combined_rrr")
    assert pool.combined_rrr() < pool.naive_sum_rrr()


def test_combined_risk_reduction_is_always_a_probability():
    for n in (1, 2, 5, 20, 100):
        pool = ee.ConditionPool("CAD", [_f(label=f"f{i}", rrr=0.9)
                                        for i in range(n)])
        assert 0.0 <= pool.combined_rrr() <= 1.0


def test_each_additional_correlated_finding_adds_less_than_the_last():
    # A second polygenic signal for the same condition is mostly the same
    # information. Marginal contributions must be strictly decreasing.
    marginals = []
    prev = 0.0
    for n in range(1, 7):
        pool = ee.ConditionPool("CAD", [_f(label=f"f{i}") for i in range(n)])
        cur = pool.combined_rrr()
        marginals.append(cur - prev)
        prev = cur
    for a, b in zip(marginals, marginals[1:]):
        assert b < a + 1e-12, f"marginal contributions not decreasing: {marginals}"


def test_a_single_finding_is_unchanged_by_pooling():
    # The correction must not penalise the ordinary one-finding case.
    lone = _f(rrr=0.30, hc=0.5)
    pool = ee.ConditionPool("CAD", [lone])
    assert pool.combined_rrr() == pytest.approx(lone.effective_rrr)


def test_cost_of_illness_is_charged_once_per_condition():
    pools = ee.pool_findings([_f(label=f"f{i}") for i in range(8)])
    ev = ee.evaluate_pools(pools, test_cost=0.0)
    row = ev["conditions"][0]
    assert row["cost_averted"] <= row["coi_cost"], (
        "a condition cannot avert more than its own full cost of illness")


def test_intervention_cost_is_charged_once_not_once_per_finding():
    # The cost-side mirror of the benefit-side double count.
    pools = ee.pool_findings([_f(label=f"f{i}", cost=500.0) for i in range(8)])
    ev = ee.evaluate_pools(pools, test_cost=0.0)
    assert ev["conditions"][0]["intervention_cost"] == 500


def test_pooling_reports_the_size_of_its_own_correction():
    pools = ee.pool_findings([_f(label=f"f{i}") for i in range(8)])
    ev = ee.evaluate_pools(pools, test_cost=0.0)
    dc = ev["double_counting"]
    assert dc["pooled_cost_averted"] < dc["naive_cost_averted"]
    assert dc["inflation_removed"] > 0
    assert dc["pct_removed"] > 0


def test_pooling_never_inflates_relative_to_naive():
    # Whatever the finding mix, the correction may only ever reduce the total.
    import random
    rng = random.Random(20260822)
    for _ in range(200):
        fs = [_f(label=f"f{i}", coi=rng.choice(["CAD", "T2D", "Depression"]),
                 p=rng.uniform(0.01, 0.6), rrr=rng.uniform(0.05, 0.8),
                 hc=rng.choice([0.1, 0.25, 0.3, 0.5, 1.0]))
              for i in range(rng.randint(1, 10))]
        ev = ee.evaluate_pools(ee.pool_findings(fs), test_cost=0.0)
        dc = ev["double_counting"]
        assert dc["pooled_cost_averted"] <= dc["naive_cost_averted"] + 1


def test_baseline_risk_is_taken_once_not_summed():
    pool = ee.ConditionPool("CAD", [_f(p=0.10), _f(p=0.25), _f(p=0.15)])
    assert pool.baseline_risk() == pytest.approx(0.25)


def test_distinct_conditions_are_not_pooled_together():
    pools = ee.pool_findings([_f(coi="CAD"), _f(coi="T2D"), _f(coi="Alzheimer")])
    assert set(pools) == {"CAD", "T2D", "Alzheimer"}
    assert all(len(p.findings) == 1 for p in pools.values())


def test_evidence_haircut_and_risk_reduction_stay_separable():
    # How much to believe a source and how much acting helps are different
    # quantities; entangling them was part of what made the old numbers
    # impossible to audit.
    strong = _f(rrr=0.30, hc=1.0)
    weak = _f(rrr=0.30, hc=0.1)
    assert weak.effective_rrr == pytest.approx(strong.effective_rrr * 0.1)


# ══════════════════════════════════════════════════════════════════════════
# Disaggregated cost-effectiveness output
# ══════════════════════════════════════════════════════════════════════════

def test_icer_is_undefined_in_the_dominance_quadrants():
    # Reporting a ratio for a dominant or dominated strategy is a classic
    # error — the ratio is negative and reads like a bargain or a disaster
    # depending on sign, when the correct answer is "ratio does not apply".
    dominant = ee.CEAResult("s", inc_cost=-500.0, inc_qaly=0.5, wtp=100_000)
    dominated = ee.CEAResult("s", inc_cost=500.0, inc_qaly=-0.5, wtp=100_000)
    assert dominant.icer is None and "dominant" in dominant.verdict
    assert dominated.icer is None and "dominated" in dominated.verdict


def test_icer_is_reported_when_the_tradeoff_is_real():
    r = ee.CEAResult("s", inc_cost=20_000.0, inc_qaly=0.5, wtp=100_000)
    assert r.icer == pytest.approx(40_000.0)
    assert "cost-effective" in r.verdict


def test_inmb_follows_the_threshold():
    r = ee.CEAResult("s", inc_cost=20_000.0, inc_qaly=0.5, wtp=100_000)
    assert r.inmb == pytest.approx(30_000.0)
    cheap = ee.CEAResult("s", inc_cost=20_000.0, inc_qaly=0.5, wtp=20_000)
    assert cheap.inmb < 0 and "not cost-effective" in cheap.verdict


def test_costs_and_qalys_are_reported_separately():
    # The old headline added cash savings to monetised health under one label.
    ev = ee.evaluate_pools(ee.pool_findings([_f()]), test_cost=100.0)
    cea = ev["cea"]
    for key in ("incremental_cost", "incremental_qaly", "icer", "inmb", "wtp"):
        assert key in cea
    assert cea["incremental_qaly"] != cea["incremental_cost"]


def test_zero_qaly_gain_yields_no_icer_rather_than_infinity():
    r = ee.CEAResult("s", inc_cost=1_000.0, inc_qaly=0.0, wtp=100_000)
    assert r.icer is None


def test_no_findings_produces_an_unavailable_result_not_a_zero():
    ev = ee.evaluate_pools({}, test_cost=100.0)
    assert ev["available"] is False


def test_unpriced_condition_contributes_nothing_rather_than_guessing():
    # A condition with no registered cost anchor must not be given an invented
    # one — it should simply contribute zero and be visible as such.
    pools = ee.pool_findings([_f(coi="SomethingUnpriced")])
    ev = ee.evaluate_pools(pools, test_cost=0.0)
    assert ev["conditions"][0]["coi_cost"] == 0
    assert ev["conditions"][0]["cost_averted"] == 0


# ══════════════════════════════════════════════════════════════════════════
# Markov model arithmetic
# ══════════════════════════════════════════════════════════════════════════

def test_within_cycle_weights_match_the_published_implementation():
    # Transcribed from darthtools::gen_wcc, including its deliberate 1-indexed
    # asymmetry. Cross-checked against the independent implementation in the
    # heor-model-replication repository when that repository is present.
    for n in (2, 3, 10, 50):
        w = ee.simpson_weights(n)
        assert len(w) == n + 1
        assert w[0] == pytest.approx(1 / 3) and w[-1] == pytest.approx(1 / 3)

    repl = os.path.expanduser("~/heor-model-replication")
    if not os.path.isdir(repl):
        pytest.skip("heor-model-replication not available for cross-check")
    sys.path.insert(0, repl)
    try:
        from cstm.wcc import gen_wcc
    except Exception:
        pytest.skip("replication cstm package not importable")
    for n in (1, 2, 3, 5, 10, 49, 50, 75):
        assert ee.simpson_weights(n) == pytest.approx(list(gen_wcc(n)))


def test_discount_weights_match_the_closed_form():
    w = ee.discount_weights(10, 0.03)
    assert w[0] == 1.0
    assert w[5] == pytest.approx(1.0 / 1.03 ** 5)


def test_life_table_loads_and_mortality_rises_with_age():
    lt = ee.life_table("Total")
    assert lt, "vendored life table should be present"
    assert lt[80] > lt[50] > lt[30]


def test_life_table_is_sex_specific():
    assert ee.life_table("Male")[60] != ee.life_table("Female")[60]


def test_missing_life_table_degrades_rather_than_crashing(monkeypatch):
    monkeypatch.setattr(ee, "_LIFE_TABLE_PATH", "/nonexistent/table.csv")
    assert ee.life_table("Total") == {}
    r = ee.run_markov(start_age=50, annual_incidence=0.01, coi_cost=10_000,
                      disutility=0.15)
    assert r.qaly > 0, "model should lose precision, not disappear"


def test_markov_cohort_is_conserved():
    r = ee.run_markov(start_age=40, annual_incidence=0.02, coi_cost=50_000,
                      disutility=0.15)
    for well, sick, dead in r.trace:
        assert well + sick + dead == pytest.approx(1.0, abs=1e-9)
        assert min(well, sick, dead) >= -1e-12


def test_markov_everyone_eventually_dies():
    r = ee.run_markov(start_age=40, annual_incidence=0.01, coi_cost=1_000,
                      disutility=0.1, max_age=100)
    assert r.trace[-1][2] > 0.5, "most of a cohort followed to 100 should be dead"


def test_reducing_risk_gains_qalys_and_averts_cost():
    base = ee.run_markov(start_age=50, annual_incidence=0.02, coi_cost=80_000,
                         disutility=0.15, rrr=0.0)
    act = ee.run_markov(start_age=50, annual_incidence=0.02, coi_cost=80_000,
                        disutility=0.15, rrr=0.40)
    assert act.qaly > base.qaly
    assert act.cost < base.cost


def test_competing_mortality_shrinks_the_value_of_late_prevention():
    # The whole reason for the structural model: preventing a disease at 85 is
    # worth less than at 45 because fewer people survive to collect it.
    young = ee.incremental_analysis(start_age=45, annual_incidence=0.01,
                                    coi_key="CAD", rrr=0.3)
    old = ee.incremental_analysis(start_age=85, annual_incidence=0.01,
                                  coi_key="CAD", rrr=0.3)
    assert young["incremental_qaly"] > old["incremental_qaly"]


def test_higher_discount_rate_lowers_present_value():
    lo = ee.run_markov(start_age=50, annual_incidence=0.02, coi_cost=80_000,
                       disutility=0.15, rate=0.0)
    hi = ee.run_markov(start_age=50, annual_incidence=0.02, coi_cost=80_000,
                       disutility=0.15, rate=0.05)
    assert hi.qaly < lo.qaly


def test_incremental_analysis_refuses_unpriced_conditions():
    r = ee.incremental_analysis(start_age=50, annual_incidence=0.01,
                                coi_key="NotAThing", rrr=0.3)
    assert r["available"] is False and "no costed anchor" in r["reason"]


# ══════════════════════════════════════════════════════════════════════════
# Perspective, reporting, validation
# ══════════════════════════════════════════════════════════════════════════

def test_societal_perspective_adds_to_healthcare_and_itemises_the_difference():
    dp = ee.dual_perspective(10_000.0, 1.0, conditions=[], wtp=100_000)
    assert dp["societal"]["cost_averted"] >= dp["healthcare_sector"]["cost_averted"]
    assert dp["delta"] == sum(a["value"] for a in dp["societal_additions"])
    for a in dp["societal_additions"]:
        assert a["basis"], "every societal addition must state its basis"


def test_healthcare_perspective_excludes_productivity():
    # A reader who rejects productivity valuation must be able to use the
    # reference case untouched.
    dp = ee.dual_perspective(10_000.0, 1.0, conditions=[], wtp=100_000)
    assert dp["healthcare_sector"]["cost_averted"] == 10_000


def test_impact_inventory_declares_what_is_left_out():
    inv = ee.impact_inventory([])
    labels = {r["item"] for r in inv}
    assert "Reproductive decisions" in labels
    excluded = [r for r in inv if r["healthcare"] == "not counted"]
    assert excluded, "an inventory that counts everything is not an inventory"
    for row in inv:
        assert row["note"], "every inventory line needs a stated reason"


def test_cheers_checklist_admits_its_gaps():
    items = ee.cheers_checklist(wtp=100_000, rate=0.03, horizon=10)
    assert len(items) >= 15
    text = " ".join(i["response"] for i in items).lower()
    assert "not addressed" in text or "not applicable" in text, (
        "a reporting checklist with no gaps is not being filled in honestly")


def test_validation_checks_run_and_pass_on_a_normal_input():
    pools = ee.pool_findings([_f(label=f"f{i}") for i in range(6)]
                             + [_f(coi="T2D"), _f(coi="Alzheimer")])
    ev = ee.evaluate_pools(pools, test_cost=100.0)
    checks = ee.validate_model(pools, ev)
    assert checks
    failed = [c for c in checks if not c["pass"]]
    assert not failed, f"internal validation failed: {failed}"
    for c in checks:
        assert c["detail"], "a validation check must report what it saw"


def test_every_costed_condition_anchor_resolves_to_a_registered_parameter():
    for coi_key, (cost_key, qaly_key) in ee.COI_KEY_TO_PARAM.items():
        assert ep.get(cost_key).value > 0, f"{coi_key} cost anchor is zero"
        assert ep.get(qaly_key).value > 0, f"{coi_key} QALY anchor is zero"


def test_engine_pulls_constants_from_the_registry_not_literals(monkeypatch):
    # If the engine kept its own copies, changing a registered parameter would
    # leave the model unchanged and the provenance table would be decorative.
    import dataclasses
    pool = ee.ConditionPool("CAD", [_f(label=f"f{i}") for i in range(6)])
    before = pool.combined_rrr()
    tightened = dataclasses.replace(ep.PARAMS["correlated_signal_penalty"],
                                    value=0.1)
    monkeypatch.setitem(ep.PARAMS, "correlated_signal_penalty", tightened)
    assert pool.combined_rrr() < before, (
        "a harsher correlation penalty must reduce the combined effect; "
        "if it does not, the engine is not reading the registry")


def test_correlation_penalty_binds_before_the_cap():
    # Worth pinning explicitly: with the penalty compounding, the combined
    # effect converges well below the cap, so the cap is a backstop rather
    # than the operative constraint. If a future penalty change makes the cap
    # bind instead, that is a real behaviour change and should be noticed.
    many = ee.ConditionPool("CAD", [_f(label=f"f{i}") for i in range(30)])
    assert many.combined_rrr() < ep.value("max_combined_rrr")
