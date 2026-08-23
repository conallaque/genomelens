"""Correctness guarantees for the pooled cost-effectiveness engine.

The engine exists to fix one specific defect: eight of the report's finding
sources route onto the same cardiometabolic anchor, and the previous model
added their risk reductions together. These tests pin the fix, and pin the
arithmetic conventions (within-cycle correction, discounting, ICER quadrants)
against the published definitions rather than against the engine's own output.
"""

import os
import sys

import pytest

from econ import engine as ee
from econ import params as ep


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
    for a, b in zip(marginals, marginals[1:], strict=False):
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


# ══════════════════════════════════════════════════════════════════════════
# Duplicate-target pooling (the personal-economics path)
# ══════════════════════════════════════════════════════════════════════════

def _item(finding, net, category="Genomic"):
    return {"finding": finding, "category": category, "net": net,
            "avoided": net // 2, "qaly_value": net // 2, "qaly": 0.1,
            "intervention": 100}


def test_same_gene_in_two_panels_is_not_valued_twice():
    # THE BUG: COMT surfaced once via the neurochemistry panel and once via
    # pharmacogenomic prescribing guidance, and both were counted in full.
    items = [_item("COMT Val158Met — psychiatric medication response", 8_340),
             _item("COMT-guided SSRI/SNRI selection", 8_210)]
    out = ee.deduplicate_by_target(items)
    assert {i["pool_target"] for i in out} == {"COMT"}
    kept = sorted(i["net"] for i in out)
    assert kept[1] == 8_340, "the strongest line keeps its full value"
    assert kept[0] < 8_210, "the duplicate must be discounted"


def test_hyphenated_and_starred_gene_names_still_match():
    # "COMT-guided" and "HLA-B*58:01" must resolve to their gene, or the
    # duplicate they form with a plain mention goes uncaught.
    assert ee._extract_target("COMT-guided prescribing") == "COMT"
    assert ee._extract_target("HLA-B*58:01 allopurinol risk") == "HLA"


def test_unrelated_findings_are_never_pooled():
    items = [_item("Avoid clopidogrel non-response", 8_904),
             _item("1 actionable wellness variant(s)", 5_048),
             _item("Targeted supplementation (vitamin D/B12/folate)", 4_070)]
    out = ee.deduplicate_by_target(items)
    assert all(i["retained"] == 1.0 for i in out), (
        "independent findings must keep their full value")
    assert len({i["pool_target"] for i in out}) == 3


def test_shape_lookalikes_are_not_mistaken_for_genes():
    # A regex for gene-shaped words reads these as genes and invents targets.
    for text in ("Intensive cardiovascular risk reduction (MI/stroke)",
                 "Avoid MACE after stent thrombosis",
                 "Targeted supplementation (vitamin D/B12/folate)"):
        assert ee._extract_target(text) == "", f"{text!r} matched a gene"


def test_pooling_is_recorded_on_the_item_not_applied_silently():
    items = [_item("COMT Val158Met", 8_000), _item("COMT-guided dosing", 4_000)]
    out = ee.deduplicate_by_target(items)
    discounted = [i for i in out if i["pool_rank"] > 0]
    assert discounted and discounted[0]["pool_note"], (
        "a reduced line must say why it was reduced")


def test_pooling_scales_every_money_field_consistently():
    items = [_item("COMT a", 1_000), _item("COMT b", 1_000)]
    out = sorted(ee.deduplicate_by_target(items), key=lambda d: d["pool_rank"])
    second = out[1]
    keep = second["retained"]
    assert second["net"] == round(1_000 * keep)
    assert second["avoided"] == round(500 * keep)
    assert second["qaly_value"] == round(500 * keep)


def test_pooling_only_ever_reduces_a_total():
    items = [_item(f"COMT variant {i}", 1_000) for i in range(5)] \
        + [_item("Unrelated finding", 2_000)]
    before = sum(i["net"] for i in items)
    after = sum(i["net"] for i in ee.deduplicate_by_target(items))
    assert after <= before


def test_custom_vocabulary_is_honoured():
    items = [_item("ZZZ9 finding one", 100), _item("ZZZ9 finding two", 50)]
    assert all(i["retained"] == 1.0
               for i in ee.deduplicate_by_target(items))
    out = ee.deduplicate_by_target(items, vocabulary=frozenset({"ZZZ9"}))
    assert any(i["retained"] < 1.0 for i in out)


# ══════════════════════════════════════════════════════════════════════════
# Uncertainty propagation
# ══════════════════════════════════════════════════════════════════════════

def _pools_for_psa():
    return ee.pool_findings([_f(coi="CAD"), _f(coi="CAD", hc=0.3),
                             _f(coi="T2D"), _f(coi="Alzheimer")])


def test_psa_varies_parameters_and_reports_an_interval():
    r = ee.run_psa(_pools_for_psa(), n=300, test_cost=100)
    assert r["n_iterations"] == 300
    assert r["n_parameters_varied"] > 10
    assert r["inmb_ci_low"] < r["mean_inmb"] < r["inmb_ci_high"], (
        "a PSA that produces no spread is not varying anything")


def test_psa_is_reproducible_for_a_given_seed():
    a = ee.run_psa(_pools_for_psa(), n=200, seed=42, test_cost=100)
    b = ee.run_psa(_pools_for_psa(), n=200, seed=42, test_cost=100)
    assert a["mean_inmb"] == b["mean_inmb"]
    c = ee.run_psa(_pools_for_psa(), n=200, seed=43, test_cost=100)
    assert c["mean_inmb"] != a["mean_inmb"], "different seeds must differ"


def test_psa_restores_the_registry_afterwards():
    # The override context manager swaps registry entries in place. If it
    # leaked, every later calculation in the run would silently use a random
    # draw instead of the base value — a spectacular and near-undetectable bug.
    before = {k: p.value for k, p in ep.PARAMS.items()}
    ee.run_psa(_pools_for_psa(), n=50, test_cost=100)
    after = {k: p.value for k, p in ep.PARAMS.items()}
    assert before == after


def test_psa_respects_documented_bounds():
    import random
    rng = random.Random(3)
    for p in ep.sampleable():
        for _ in range(200):
            v = ep.draw(rng, p)
            if p.low is not None:
                assert v >= p.low - 1e-9, f"{p.key} drew below its stated low"
            if p.high is not None:
                assert v <= p.high + 1e-9, f"{p.key} drew above its stated high"


def test_parameters_without_a_published_spread_are_held_fixed():
    # Inventing uncertainty is the same error as inventing a value.
    import random
    rng = random.Random(1)
    fixed = [p for p in ep.PARAMS.values()
             if p.dist == "fixed" or (p.se is None and p.low is None)]
    assert fixed
    for p in fixed:
        if p.dist == "fixed":
            assert ep.draw(rng, p) == p.value


def test_ceac_is_monotonic_in_willingness_to_pay():
    # More willingness to pay can never make a strategy less likely to be
    # cost-effective. A non-monotonic curve means the thresholds were resampled
    # independently and the curve is showing Monte Carlo noise.
    curve = ee.ceac(_pools_for_psa(), n=200, test_cost=100)
    ps = [c["p_cost_effective"] for c in curve]
    assert ps == sorted(ps), f"CEAC not monotonic: {ps}"
    assert all(0.0 <= x <= 1.0 for x in ps)


def test_tornado_ranks_by_swing_and_names_assumption_tiers():
    rows = ee.tornado(_pools_for_psa(), test_cost=100)
    assert rows
    swings = [r["swing"] for r in rows]
    assert swings == sorted(swings, reverse=True)
    for r in rows:
        assert r["tier"] in ep.TIERS
        assert r["low_value"] < r["high_value"]


def test_tornado_swings_come_from_the_registry_range():
    rows = {r["parameter"]: r for r in ee.tornado(_pools_for_psa(), top=50)}
    for key, row in rows.items():
        p = ep.get(key)
        assert row["low_value"] == p.low and row["high_value"] == p.high


def test_override_context_manager_restores_on_exception():
    base = ep.value("discount_rate")
    try:
        with ep.overridden({"discount_rate": 0.99}):
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    assert ep.value("discount_rate") == base


def test_psa_without_rebuild_understates_uncertainty():
    # REGRESSION GUARD. Baseline risk, effect size and intervention cost are
    # read when a Finding is constructed. An earlier version ran PSA over
    # pools built from the base case, so those parameters stayed pinned while
    # only the cost-of-illness terms varied. The result was a strategy that
    # could not lose: cost-saving in 100% of simulations and a CEAC reading
    # 100% at every threshold, including $0/QALY.
    def build():
        return ee.pool_findings([
            ee.Finding(label="f", coi_key="CAD",
                       p_event=ep.value("baseline_event_probability"),
                       rrr=ep.value("actionable_rrr"),
                       intervention_cost=ep.value("intervention_cost_standard"))])

    pinned = ee.run_psa(build(), n=400, test_cost=100)
    live = ee.run_psa(build(), n=400, test_cost=100, rebuild=build)
    span_pinned = pinned["inmb_ci_high"] - pinned["inmb_ci_low"]
    span_live = live["inmb_ci_high"] - live["inmb_ci_low"]
    assert span_live > span_pinned, (
        "rebuilding pools inside each draw must widen the interval; if it "
        "does not, the finding-level parameters are still pinned")


def test_finding_defaults_are_registry_backed_not_literals():
    # The two numbers the whole benefit side rests on must be varyable.
    from econ import value_of_information as voi
    econ = {"findings_with_economics": [
        {"finding": "x", "category": "Polygenic Risk", "qaly_gain": 0.5}]}
    base = voi._collect(econ, None, None, [])[0]
    assert base["p_event"] == ep.value("baseline_event_probability")
    assert base["rrr"] == ep.value("actionable_rrr")
    with ep.overridden({"actionable_rrr": 0.05,
                        "baseline_event_probability": 0.4}):
        drawn = voi._collect(econ, None, None, [])[0]
    assert drawn["rrr"] == 0.05 and drawn["p_event"] == 0.4


def test_each_condition_has_its_own_quality_of_life_decrement():
    # Seven conditions once shared qaly_loss_mace — a cardiovascular event's
    # decrement standing in for dementia, depression and kidney stones. That
    # is wrong on its face, and it made the shared placeholder dominate the
    # tornado, so the report named it the key driver when it was really the
    # most overloaded constant.
    qaly_params = [q for _, q in ee.COI_KEY_TO_PARAM.values()]
    from collections import Counter
    overloaded = [k for k, n in Counter(qaly_params).items() if n > 2]
    assert not overloaded, (
        f"QALY anchors reused across more than two conditions: {overloaded}")
    # Conditions with genuinely different decrements must not share one.
    distinct = {ee.COI_KEY_TO_PARAM[c][1]
                for c in ("CAD", "Alzheimer", "Depression", "Urologic")}
    assert len(distinct) == 4, (
        "dementia, depression, kidney stones and heart attacks do not have "
        "the same quality-of-life decrement")


def test_dementia_costs_more_quality_of_life_than_kidney_stones():
    # A cheap ordering check that would have caught the shared-anchor bug.
    dementia = ep.value(ee.COI_KEY_TO_PARAM["Alzheimer"][1])
    urologic = ep.value(ee.COI_KEY_TO_PARAM["Urologic"][1])
    assert dementia > urologic * 3, (
        f"dementia decrement {dementia} vs urologic {urologic} — implausible")


# ══════════════════════════════════════════════════════════════════════════
# Adherence — efficacy vs. effectiveness
# ══════════════════════════════════════════════════════════════════════════

def _adh_findings(adherence, coi_key="CAD"):
    """Two CAD findings differing only in the adherence multiplier."""
    return [
        ee.Finding(label="statin", coi_key=coi_key, p_event=0.20, rrr=0.27,
                   haircut=1.0, intervention_cost=400.0, adherence=adherence),
        ee.Finding(label="prs", coi_key=coi_key, p_event=0.15, rrr=0.20,
                   haircut=0.8, intervention_cost=200.0, adherence=adherence),
    ]


def test_perfect_adherence_reproduces_the_pre_adherence_model():
    # The regression guard. Adherence must be a multiplier that vanishes at
    # 1.0, not a re-derivation that shifts the answer even when nobody stops.
    # Note this constructs Findings explicitly rather than wrapping a
    # pre-built pool in ep.overridden — adherence is fixed on the instance at
    # construction, so an override around an existing pool reaches nothing and
    # the test would pass without testing anything.
    pool = ee.pool_findings(_adh_findings(1.0))["CAD"]
    assert pool.combined_rrr() == pytest.approx(pool.pooled_efficacy_rrr())
    assert pool.adherence() == 1.0
    for f in pool.findings:
        assert f.effective_rrr == pytest.approx(f.efficacy_rrr)


def test_adherence_reduces_the_realised_benefit():
    half = ee.pool_findings(_adh_findings(0.5))["CAD"]
    full = ee.pool_findings(_adh_findings(1.0))["CAD"]
    assert half.combined_rrr() < full.combined_rrr()
    assert half.pooled_efficacy_rrr() == pytest.approx(full.pooled_efficacy_rrr()), \
        "efficacy is a property of the evidence and must not move with adherence"


def test_adherence_is_applied_before_the_product_not_after():
    # Adherence attenuates each intervention's own effect, and the attenuated
    # effects then combine — so the discount belongs inside the product. The
    # combination rule is concave, which means the two orderings genuinely
    # differ and the correct one is the slightly LESS conservative of the two.
    # Pinning the direction here so nobody "fixes" it toward the smaller
    # number on the assumption that smaller means safer.
    half = ee.pool_findings(_adh_findings(0.5))["CAD"]
    post_scaling = 0.5 * half.pooled_efficacy_rrr()
    assert half.combined_rrr() > post_scaling
    assert half.combined_rrr() == pytest.approx(post_scaling, rel=0.05), \
        "the two orderings should differ modestly, not by an order of magnitude"


def test_the_pooling_correction_does_not_absorb_the_adherence_discount():
    # Both are shrinkages of the same headline. If double_count_avoided were
    # measured against the adherence-discounted figure it would silently grow,
    # and the report's "size of the double-counting correction" banner would
    # be reporting two different things under one label.
    half = ee.pool_findings(_adh_findings(0.5))["CAD"]
    full = ee.pool_findings(_adh_findings(1.0))["CAD"]
    assert half.to_dict()["double_count_avoided"] == \
        pytest.approx(full.to_dict()["double_count_avoided"])
    assert half.to_dict()["adherence_drag_rrr"] > 0
    assert full.to_dict()["adherence_drag_rrr"] == pytest.approx(0.0)


def test_the_two_corrections_sum_to_the_total_shrinkage():
    d = ee.pool_findings(_adh_findings(0.5))["CAD"].to_dict()
    total = d["naive_additive_rrr"] - d["combined_rrr"]
    assert d["double_count_avoided"] + d["adherence_drag_rrr"] == \
        pytest.approx(total, abs=1e-4)


def test_intervention_cost_scales_with_adherence():
    # The cost side has to move with the benefit side. Charging the full
    # course while crediting half the effect is the mirror image of the
    # double-count this model was built to remove.
    half = ee.pool_findings(_adh_findings(0.5))["CAD"]
    full = ee.pool_findings(_adh_findings(1.0))["CAD"]
    assert half.intervention_cost() == pytest.approx(0.5 * full.intervention_cost())


def test_mixed_adherence_within_a_pool_is_rejected():
    # intervention_cost() scales the pool by one multiplier. If adherence ever
    # becomes per-finding that arithmetic is wrong, and this assertion is what
    # makes the change fail loudly instead of quietly mis-costing.
    fs = _adh_findings(0.5)
    fs[1].adherence = 0.9
    with pytest.raises(AssertionError):
        ee.pool_findings(fs)["CAD"].intervention_cost()


def test_every_valued_condition_has_an_adherence_archetype():
    # An unmapped condition falls back to adherence_default. That fallback
    # exists for safety, not as a resting place.
    missing = set(ee.COI_KEY_TO_PARAM) - set(ee.ADHERENCE_BY_COI_KEY)
    assert not missing, f"conditions with no adherence archetype: {sorted(missing)}"


def test_unmapped_conditions_do_not_default_to_perfect_adherence():
    assert ee.adherence_for("NotAConditionKey") == pytest.approx(
        ep.value("adherence_default"))
    assert ee.adherence_for("NotAConditionKey") < 1.0


def test_adherence_archetypes_all_resolve_to_registered_parameters():
    for key in set(ee.ADHERENCE_BY_COI_KEY.values()) | {"adherence_default"}:
        assert key in ep.PARAMS, f"{key} is mapped but not registered"
        assert 0.0 < ep.value(key) <= 1.0


def test_report_separates_efficacy_from_effectiveness():
    ev = ee.evaluate_pools(ee.pool_findings(_adh_findings(0.5)), test_cost=300.0)
    a = ev["adherence"]
    assert a["efficacy_qaly"] > a["effectiveness_qaly"] > 0
    assert a["value_lost_to_non_adherence"] > 0
    assert 0 < a["pct_of_benefit_lost"] < 100


def test_the_fixed_test_cost_is_what_moves_the_icer():
    # Adherence scales benefit and ongoing cost together, so on its own it
    # leaves cost-per-QALY roughly alone. The one-off test cost does not
    # scale, so it amortises over fewer realised QALYs. With no test cost the
    # ratio should barely move; with one it should get materially worse.
    def ratio(adh, test_cost):
        ev = ee.evaluate_pools(ee.pool_findings(_adh_findings(adh)),
                               test_cost=test_cost)
        c = ev["cea"]
        return c["incremental_cost"] / c["incremental_qaly"]

    free_gap = abs(ratio(0.5, 0.0) - ratio(1.0, 0.0))
    paid_gap = abs(ratio(0.5, 300.0) - ratio(1.0, 300.0))
    assert paid_gap > free_gap * 5, (
        "the adherence penalty to cost-effectiveness should come from the "
        "unscaled fixed test cost, not from the intervention itself")


def test_adherence_parameters_are_varied_in_sensitivity_analysis():
    # A parameter this influential cannot be a pinned point estimate.
    sampleable = {p.key for p in ep.sampleable()}
    for key in set(ee.ADHERENCE_BY_COI_KEY.values()):
        assert key in sampleable, f"{key} is held fixed in PSA"


def test_personal_economics_uses_the_same_adherence_archetypes():
    # The two economic pages describe one genome. When the pooled payer
    # analysis charged adherence and this one did not, the report said 62% of
    # the benefit was lost on one page and reported the undiscounted total as
    # the headline dollar figure on the other.
    from econ import health_economics as he
    assert set(he._CATEGORY_ADHERENCE.values()) <= (
        set(ee.ADHERENCE_BY_COI_KEY.values()) | {"adherence_default"}), \
        "the two pages must draw on one set of adherence parameters"
    for key in he._CATEGORY_ADHERENCE.values():
        assert key in ep.PARAMS


def test_every_personal_economics_category_is_mapped():
    # An unmapped category silently falls back to adherence_default. That is a
    # safe fallback, not a place to leave real categories sitting.
    import re

    from econ import health_economics as he
    src = open(he.__file__).read()
    used = set(re.findall(r'^\s+add\("([^"]+)"', src, re.M))
    assert used, "could not find the add() call sites to check"
    missing = used - set(he._CATEGORY_ADHERENCE)
    assert not missing, f"categories with no adherence archetype: {sorted(missing)}"


def test_personal_economics_categories_do_not_default_to_perfect_adherence():
    from econ import health_economics as he
    for cat in list(he._CATEGORY_ADHERENCE) + ["NoSuchCategory"]:
        assert he._adherence_for_category(cat) < 1.0


def test_the_efficacy_counterfactual_reconciles_with_the_discounted_total():
    # "At full adherence these findings would total X" has to be the same
    # findings pooled the same way, not a figure captured before the
    # correlated-target discount. Captured pre-pooling it overstated the
    # counterfactual by ~$4k and made the adherence drag look bigger than it is.
    from econ import health_economics as he
    items = [{"net": 1000, "adherence": 0.5}, {"net": 300, "adherence": 0.35}]
    eff = sum(i["net"] / i["adherence"] for i in items)
    assert eff == pytest.approx(2000 + 857.142857)
    assert hasattr(he, "_render_adherence_basis_html")
    html = he._render_adherence_basis_html({
        "adherence_applied": True, "efficacy_net": round(eff),
        "adherence_drag": round(eff - 1300), "mean_adherence": 0.425})
    assert "At full adherence" in html and "$2,857" in html


def test_the_adherence_basis_is_stated_whenever_it_is_applied():
    # A page that silently reports a smaller number than it used to is its own
    # kind of dishonesty.
    from econ import health_economics as he
    assert he._render_adherence_basis_html({"adherence_applied": False}) == ""
    assert "real-world figures" in he._render_adherence_basis_html({
        "adherence_applied": True, "efficacy_net": 100, "adherence_drag": 50,
        "mean_adherence": 0.5})


def test_the_vendored_life_table_is_actually_reachable():
    # THE BUG THIS CATCHES. life_table() degrades to {} on OSError so the model
    # loses precision rather than disappearing — which means a wrong path is
    # invisible. Moving this module into econ/ broke the __file__-relative
    # lookup exactly that way: every import still passed, every Markov run
    # still returned numbers, and the age-specific mortality had silently been
    # replaced by a constant hazard.
    tbl = ee.life_table()
    assert tbl, f"life table unreachable at {ee._LIFE_TABLE_PATH}"
    assert len(tbl) > 100, f"expected ~111 single-year ages, got {len(tbl)}"
    # Mortality must rise with age, or the file was parsed off the wrong column.
    assert tbl[80] > tbl[40] > tbl[20]


def test_the_life_table_stays_at_the_repository_root():
    # It cannot move into econ/ alongside this module: .gitignore protects the
    # repo from DNA with a blanket *.csv and one path-anchored exception for
    # data/LifeTable_USA_Mx_2015.csv. Relocate the file and that exception stops
    # matching, so the life table becomes untracked while the DNA guard looks
    # untouched.
    import os
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(ee.__file__))))
    gitignore = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(ee.__file__))), ".gitignore")
    if os.path.exists(gitignore):
        rules = open(gitignore).read()
        assert "*.csv" in rules, "the blanket DNA guard is gone"
        assert "!data/LifeTable_USA_Mx_2015.csv" in rules
        assert ee._LIFE_TABLE_PATH.endswith(
            os.path.join("data", "LifeTable_USA_Mx_2015.csv"))
    assert os.path.isdir(os.path.dirname(ee._LIFE_TABLE_PATH))
