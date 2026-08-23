"""
Decision-analytic outputs: what to resolve, what to choose, what it costs
=========================================================================

:mod:`econ_engine` answers "is acting on these findings worth it, and how
uncertain is that". This module answers the questions a decision-maker asks
next, none of which a point estimate or even a confidence interval addresses.

**What is it worth to stop guessing?** The tornado already says the conclusion
rests mostly on a handful of declared assumptions. That is a diagnosis without
a prescription — it does not say whether resolving them would be worth the
trouble. :func:`evppi` puts a dollar figure on perfect information about each
parameter individually, which converts "this assumption matters" into "finding
out would be worth $X, and here is the ranking".

**Where does the conclusion flip?** :func:`breakeven` solves for the value of a
parameter at which net benefit crosses zero. When an assumption drives the
result, the useful statement is not its point value but how far it can move
before the recommendation changes.

**Which strategy, not just whether to act?** A two-arm comparison cannot
express the real choice, which includes doing less and doing more.
:func:`efficiency_frontier` compares several strategies with proper dominance
and *extended* dominance — the latter routinely omitted, which lets a strategy
that no rational decision-maker would pick appear on the frontier.

**What will it cost the budget?** Cost-effectiveness and affordability are
different questions with different answers. :func:`budget_impact` reports
undiscounted cash over a short horizon for a defined population, which is what
a payer actually needs.

**Does it work the same for everyone?** :func:`subgroup_analysis` re-runs the
structural model across age and sex, and :func:`distributional_cea` applies
equity weights — closing a gap this project's own CHEERS checklist declares.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Sequence

from . import engine as ee
from . import params as ep
import itertools

__all__ = [
    "analyze_decision_layer",
    "breakeven",
    "budget_impact",
    "distributional_cea",
    "efficiency_frontier",
    "evpi",
    "evppi",
    "population_evpi",
    "subgroup_analysis",
]


# ══════════════════════════════════════════════════════════════════════════
# Value of information
# ══════════════════════════════════════════════════════════════════════════

def _inmb_draws(rebuild: Callable[[], dict], *, n: int, seed: int,
                test_cost: float, wtp: float,
                fixed: dict[str, float] | None = None) -> list[float]:
    """Incremental net monetary benefit across ``n`` parameter draws.

    ``fixed`` pins named parameters at given values — the mechanism behind
    EVPPI, which asks what happens when one parameter becomes known.
    """
    rng = random.Random(seed)
    out: list[float] = []
    for _ in range(max(1, int(n))):
        draw = ep.sample_all(rng)
        if fixed:
            draw.update(fixed)
        with ep.overridden(draw):
            cea = ee.evaluate_pools(rebuild(), wtp=wtp, test_cost=test_cost)["cea"]
        out.append(float(cea["inmb"]))
    return out


def evpi(rebuild: Callable[[], dict], *, n: int = 2000, seed: int = 20260823,
         test_cost: float = 0.0, wtp: float | None = None) -> dict:
    """Expected value of perfect information, per person.

    The cost of making the wrong choice given current uncertainty. With two
    strategies — act, or don't — the comparator's net benefit is zero, so this
    reduces to ``E[max(0, INMB)] - max(0, E[INMB])``: the average regret you
    could avoid if you knew everything in advance.

    An EVPI near zero does not mean the model is certain. It means the
    uncertainty does not straddle the decision — you would act either way, so
    resolving it changes nothing you would do.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    draws = _inmb_draws(rebuild, n=n, seed=seed, test_cost=test_cost, wtp=wtp)
    k = len(draws) or 1
    mean = sum(draws) / k
    with_info = sum(max(0.0, d) for d in draws) / k
    without = max(0.0, mean)
    return {
        "available": True,
        "evpi_per_person": round(max(0.0, with_info - without)),
        "mean_inmb": round(mean),
        "p_current_choice_wrong": round(
            sum(1 for d in draws if (d > 0) != (mean > 0)) / k, 4),
        "n_iterations": k,
        "wtp": round(wtp),
        "interpretation": (
            "The most a decision-maker should pay to eliminate all parameter "
            "uncertainty before choosing. Near zero means the uncertainty "
            "does not straddle the decision, not that the model is precise."),
    }


def evppi(rebuild: Callable[[], dict], *, parameters: Sequence[str],
          n_outer: int = 150, n_inner: int = 50, seed: int = 20260823,
          test_cost: float = 0.0, wtp: float | None = None) -> list[dict]:
    """Expected value of partial perfect information, per parameter.

    Two-level Monte Carlo: draw the parameter of interest, average net benefit
    over everything else conditional on it, then take the expected maximum.
    The difference from the unconditional maximum is what learning *that one
    parameter* would be worth.

    This is the answer the tornado cannot give. A tornado ranks parameters by
    how much they move the number; EVPPI ranks them by how much they move the
    *decision*, which is not the same ordering — a parameter can swing the
    total wildly without ever changing what you would do.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)

    rows: list[dict] = []
    for i, key in enumerate(parameters):
        if key not in ep.PARAMS:
            continue
        param = ep.get(key)
        rng = random.Random(seed + 1000 + i)
        conditional_means: list[float] = []
        for j in range(max(1, int(n_outer))):
            phi = ep.draw(rng, param)
            inner = _inmb_draws(rebuild, n=n_inner, seed=seed + 50_000 + j,
                                test_cost=test_cost, wtp=wtp,
                                fixed={key: phi})
            conditional_means.append(sum(inner) / (len(inner) or 1))
        k = len(conditional_means) or 1
        # The unconditional mean MUST come from the same draws as the
        # conditional ones. Estimating it from a separate sample lets Monte
        # Carlo noise alone produce a positive EVPPI — which is how an earlier
        # version reported $2,299 of partial information value against an EVPI
        # of $0, a mathematical impossibility (EVPPI is bounded above by EVPI).
        grand_mean = sum(conditional_means) / k
        expected_max = sum(max(0.0, m) for m in conditional_means) / k
        rows.append({
            "parameter": key,
            "tier": param.tier,
            "units": param.units,
            "evppi_per_person": round(max(0.0, expected_max
                                          - max(0.0, grand_mean))),
            "share_of_evpi": None,   # filled in by the caller once EVPI known
        })
    rows.sort(key=lambda r: -r["evppi_per_person"])
    return rows


def population_evpi(evpi_per_person: float, *, population: int = 1_000_000,
                    years: int = 10, rate: float | None = None) -> dict:
    """Scale per-person EVPI to a population over a decision horizon.

    Research is funded against population value, not per-person value, so this
    is the figure that decides whether a study is worth commissioning. The
    population is a stated input rather than a claim about who this applies to.
    """
    rate = ep.value("discount_rate") if rate is None else float(rate)
    pv = sum(1.0 / (1.0 + rate) ** t for t in range(max(0, int(years))))
    incident = population / max(1, years)      # spread arrivals over the horizon
    return {
        "population": population,
        "years": years,
        "discount_rate": rate,
        "population_evpi": round(evpi_per_person * incident * pv),
        "note": ("Per-person EVPI × people reaching the decision each year × "
                 "discounted horizon. The population figure is an input, not "
                 "a claim about this report's reach."),
    }


# ══════════════════════════════════════════════════════════════════════════
# Threshold / breakeven
# ══════════════════════════════════════════════════════════════════════════

def breakeven(rebuild: Callable[[], dict], *, parameter: str,
              test_cost: float = 0.0, wtp: float | None = None,
              steps: int = 60) -> dict:
    """Value of one parameter at which net benefit crosses zero.

    When a judgement call drives the conclusion, its point value is less
    interesting than its margin: how wrong could this assumption be before the
    recommendation changes? A breakeven comfortably outside the parameter's
    plausible range is reassuring; one inside it means the recommendation is
    genuinely undetermined.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    if parameter not in ep.PARAMS:
        return {"available": False, "reason": f"unknown parameter {parameter!r}"}
    p = ep.get(parameter)
    lo = p.low if p.low is not None else p.value * 0.1
    hi = p.high if p.high is not None else p.value * 3.0
    if lo == hi:
        return {"available": False, "reason": "parameter has no range"}

    def inmb_at(v: float) -> float:
        with ep.overridden({parameter: v}):
            return float(ee.evaluate_pools(rebuild(), wtp=wtp,
                                           test_cost=test_cost)["cea"]["inmb"])

    n = max(4, int(steps))
    xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
    ys = [inmb_at(x) for x in xs]

    crossing = None
    for (x0, y0), (x1, y1) in zip(zip(xs, ys, strict=False), zip(xs[1:], ys[1:], strict=False), strict=False):
        if (y0 <= 0 <= y1) or (y1 <= 0 <= y0):
            crossing = x0 if y1 == y0 else x0 + (x1 - x0) * (0 - y0) / (y1 - y0)
            break

    base = inmb_at(p.value)
    return {
        "available": True,
        "parameter": parameter,
        "tier": p.tier,
        "units": p.units,
        "base_value": p.value,
        "base_inmb": round(base),
        "range": [lo, hi],
        "breakeven_value": (round(crossing, 4) if crossing is not None else None),
        "crosses_within_range": crossing is not None,
        "margin": (round(abs(crossing - p.value), 4)
                   if crossing is not None else None),
        "interpretation": (
            f"Net benefit changes sign at {crossing:.4g} {p.units}; the base "
            f"case is {p.value:g}."
            if crossing is not None else
            f"Net benefit keeps its sign across the whole plausible range "
            f"({lo:g}–{hi:g} {p.units}), so the recommendation does not turn "
            f"on this parameter."),
    }


# ══════════════════════════════════════════════════════════════════════════
# Multi-strategy comparison
# ══════════════════════════════════════════════════════════════════════════

def efficiency_frontier(strategies: Sequence[dict],
                        wtp: float | None = None) -> dict:
    """Rank strategies by cost, applying dominance and extended dominance.

    Each strategy is ``{"name", "cost", "qaly"}`` in absolute terms.

    Extended dominance is the step usually skipped: a strategy can be cheaper
    and less effective than the next one up and still be irrational to choose,
    because a mix of its neighbours dominates it. Omitting the check leaves
    options on the frontier that no decision-maker should pick, and quietly
    changes the ICERs of everything above them.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    rows = sorted(({"name": s["name"], "cost": float(s["cost"]),
                    "qaly": float(s["qaly"])} for s in strategies),
                  key=lambda r: (r["cost"], -r["qaly"]))

    # Strict dominance: costs at least as much and delivers no more health.
    for r in rows:
        r["status"] = "on frontier"
        r["icer"] = None
    for r in rows:
        if any(o is not r and o["cost"] <= r["cost"] and o["qaly"] >= r["qaly"]
               and (o["cost"] < r["cost"] or o["qaly"] > r["qaly"])
               for o in rows):
            r["status"] = "dominated"

    # Extended dominance: iteratively drop any surviving strategy whose ICER
    # exceeds that of the next surviving strategy above it.
    while True:
        live = [r for r in rows if r["status"] == "on frontier"]
        if len(live) < 3:
            break
        icers = []
        for prev, cur in itertools.pairwise(live):
            dq = cur["qaly"] - prev["qaly"]
            icers.append((cur, (cur["cost"] - prev["cost"]) / dq
                          if dq > 0 else float("inf")))
        dropped = False
        for (cur, ic), (_, nxt) in itertools.pairwise(icers):
            if ic > nxt:
                cur["status"] = "extendedly dominated"
                dropped = True
                break
        if not dropped:
            break

    live = [r for r in rows if r["status"] == "on frontier"]
    for prev, cur in itertools.pairwise(live):
        dq = cur["qaly"] - prev["qaly"]
        cur["icer"] = round((cur["cost"] - prev["cost"]) / dq) if dq > 0 else None

    best = max(rows, key=lambda r: r["qaly"] * wtp - r["cost"])
    for r in rows:
        r["inmb"] = round(r["qaly"] * wtp - r["cost"])
        r["cost"] = round(r["cost"])
        r["qaly"] = round(r["qaly"], 4)
    return {
        "available": bool(rows),
        "wtp": round(wtp),
        "strategies": rows,
        "recommended": best["name"],
        "note": ("Strategies are ordered by cost. 'Dominated' costs more for "
                 "no more health; 'extendedly dominated' is beaten by a mix "
                 "of its neighbours — a check often skipped, which leaves "
                 "unchoosable options on the frontier."),
    }


# ══════════════════════════════════════════════════════════════════════════
# Budget impact
# ══════════════════════════════════════════════════════════════════════════

def budget_impact(*, per_person_cost: float = 0.0,
                  per_person_offset: float = 0.0,
                  population: int = 1_000_000,
                  years: int = 5, **kwargs) -> dict:
    """Payer affordability — delegates to :mod:`markov_model`.

    This module briefly carried its own budget-impact implementation, which
    was a weaker duplicate of one that already existed: ``markov_model``
    follows the ISPOR BIA conventions properly, ramps offsets as the
    intervention takes effect rather than applying a flat haircut, and reports
    per-member-per-month — the metric coverage decisions are actually argued
    on. Two budget-impact analyses in one report would give a payer two
    answers to one question, which is the defect this codebase has spent the
    last several passes removing elsewhere.

    Kept as a named entry point so callers here do not have to know which
    module owns the calculation.
    """
    from . import markov as _mk
    return _mk.budget_impact(plan_members=population,
                             test_cost=max(0.0, float(per_person_cost)),
                             horizon_years=int(years))


# ══════════════════════════════════════════════════════════════════════════
# Heterogeneity
# ══════════════════════════════════════════════════════════════════════════

def subgroup_analysis(*, coi_key: str = "CAD", rrr: float | None = None,
                      annual_incidence: float | None = None,
                      ages: Sequence[int] = (30, 40, 50, 60, 70, 80),
                      sexes: Sequence[str] = ("Female", "Male"),
                      wtp: float | None = None) -> dict:
    """Re-run the structural model across age and sex.

    An average result hides that the same intervention is excellent value at
    45 and poor value at 80, because competing mortality decides how much of
    a prevented event you live to enjoy. The life table is already vendored;
    not using it for heterogeneity would waste the most defensible input in
    the model.
    """
    rrr = ep.value("actionable_rrr") if rrr is None else float(rrr)
    if annual_incidence is None:
        annual_incidence = ep.value("subgroup_annual_incidence")
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    rows: list[dict] = []
    for sex in sexes:
        for age in ages:
            r = ee.incremental_analysis(
                start_age=age, annual_incidence=annual_incidence,
                coi_key=coi_key, rrr=rrr, sex=sex, wtp=wtp)
            if not r.get("available"):
                continue
            rows.append({
                "sex": sex, "age": age,
                "incremental_cost": r["incremental_cost"],
                "incremental_qaly": r["incremental_qaly"],
                "icer": r["icer"], "inmb": r["inmb"],
                "cost_effective": r["inmb"] > 0,
            })
    best = max(rows, key=lambda r: r["inmb"], default=None)
    worst = min(rows, key=lambda r: r["inmb"], default=None)
    return {
        "available": bool(rows),
        "condition": coi_key,
        "rows": rows,
        "best": best, "worst": worst,
        "spread": (round(best["inmb"] - worst["inmb"]) if best and worst else 0),
        "illustrative": True,
        "annual_incidence": annual_incidence,
        "note": (f"Illustrative, at {annual_incidence:.1%} annual incidence — "
                 f"this shows how the value of prevention changes with age, "
                 f"not a personalised risk estimate. Differences across age "
                 f"come mostly from competing mortality: a prevented event is "
                 f"worth less to someone less likely to survive to "
                 f"experience it."),
    }


def distributional_cea(baseline_qaly_by_group: dict[str, float],
                       gain_by_group: dict[str, float],
                       *, inequality_aversion: float | None = None) -> dict:
    """Equity-weighted analysis using an Atkinson social welfare function.

    Standard cost-effectiveness is distribution-blind: a QALY counts the same
    whoever gets it. Distributional CEA asks whether a programme widens or
    narrows health inequality, and prices that alongside total health.

    ``inequality_aversion`` is the Atkinson parameter. The value used in
    published English DCEA work is around 11, which is far higher than
    economic intuition suggests and reflects strong measured public
    preference for reducing health inequality. It is reported as an input
    because the result is sensitive to it and reasonable people disagree.
    """
    if inequality_aversion is None:
        inequality_aversion = ep.value("inequality_aversion")
    groups = sorted(set(baseline_qaly_by_group) | set(gain_by_group))
    if not groups:
        return {"available": False}
    base = [max(1e-9, float(baseline_qaly_by_group.get(g, 0.0))) for g in groups]
    after = [b + float(gain_by_group.get(g, 0.0)) for g, b in zip(groups, base, strict=False)]

    def ede(xs: list[float]) -> float:
        """Equally-distributed equivalent under Atkinson aversion."""
        n = len(xs) or 1
        e = float(inequality_aversion)
        if abs(e - 1.0) < 1e-9:
            import math
            return math.exp(sum(math.log(max(1e-9, x)) for x in xs) / n)
        s = sum(max(1e-9, x) ** (1.0 - e) for x in xs) / n
        return s ** (1.0 / (1.0 - e))

    def gini(xs: list[float]) -> float:
        n = len(xs)
        if n < 2:
            return 0.0
        s = sorted(xs)
        tot = sum(s) or 1.0
        cum = sum((i + 1) * v for i, v in enumerate(s))
        return round((2.0 * cum) / (n * tot) - (n + 1.0) / n, 4)

    ede_before, ede_after = ede(base), ede(after)
    mean_before = sum(base) / len(base)
    mean_after = sum(after) / len(after)
    return {
        "available": True,
        "inequality_aversion": inequality_aversion,
        "groups": [
            {"group": g, "baseline_qaly": round(b, 3),
             "qaly_after": round(a, 3), "gain": round(a - b, 4)}
            for g, b, a in zip(groups, base, after, strict=False)],
        "mean_qaly_before": round(mean_before, 4),
        "mean_qaly_after": round(mean_after, 4),
        "ede_before": round(ede_before, 4),
        "ede_after": round(ede_after, 4),
        "ede_gain": round(ede_after - ede_before, 4),
        "gini_before": gini(base),
        "gini_after": gini(after),
        "reduces_inequality": gini(after) < gini(base),
        "note": ("The equally-distributed equivalent is the health level "
                 "which, if everyone had it, society would value as much as "
                 "the actual unequal distribution. A gain larger than the "
                 "gain in the plain mean means the programme is "
                 "equity-improving as well as health-improving."),
    }


# ══════════════════════════════════════════════════════════════════════════
# Orchestration
# ══════════════════════════════════════════════════════════════════════════

def analyze_decision_layer(rebuild: Callable[[], dict], *,
                           tornado_rows: Sequence[dict] | None = None,
                           test_cost: float = 0.0,
                           wtp: float | None = None,
                           age: float = 40.0,
                           fast: bool = False) -> dict:
    """Run the decision-analytic layer and assemble it for the report."""
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    out: dict = {"available": True, "wtp": round(wtp)}

    # Which parameters to interrogate: the ones the tornado says matter, since
    # spending EVPPI iterations on a parameter with no swing is wasted effort.
    keys = [r["parameter"] for r in (tornado_rows or [])[:6]] or [
        "actionable_rrr", "baseline_event_probability"]

    v = evpi(rebuild, n=250 if fast else 1200, test_cost=test_cost, wtp=wtp)
    out["evpi"] = v
    out["population_evpi"] = population_evpi(v["evpi_per_person"])

    rows = evppi(rebuild, parameters=keys[:4 if fast else 6],
                 n_outer=25 if fast else 90, n_inner=12 if fast else 30,
                 test_cost=test_cost, wtp=wtp)
    total = v["evpi_per_person"] or 1
    for r in rows:
        r["share_of_evpi"] = round(min(1.0, r["evppi_per_person"] / total), 3)
    out["evppi"] = rows

    out["breakeven"] = [b for b in
                        (breakeven(rebuild, parameter=k, test_cost=test_cost,
                                   wtp=wtp, steps=30 if fast else 60)
                         for k in keys[:4]) if b.get("available")]

    # Use the person's own dominant condition where there is one, rather
    # than always showing CAD regardless of what their findings say.
    _top = max(ee.evaluate_pools(rebuild(), wtp=wtp,
                                 test_cost=test_cost)["conditions"],
               key=lambda c: c.get("inmb", 0), default=None)
    out["subgroups"] = subgroup_analysis(
        coi_key=(_top or {}).get("condition") or "CAD", wtp=wtp)

    # Budget impact for a mid-sized payer, using this report's own test cost
    # and the pooled cost offset per person.
    base = ee.evaluate_pools(rebuild(), wtp=wtp, test_cost=test_cost)["cea"]
    out["budget_impact"] = budget_impact(
        per_person_cost=float(test_cost) + float(base["intervention_cost"]))

    # Distributional analysis across the age strata the subgroup model just
    # produced: younger groups start with more remaining health, so gains that
    # favour them widen the gap.
    sub = out["subgroups"]
    if sub.get("available"):
        by_age = {f"age {r['age']}": max(0.5, (100 - r["age"]) * 0.2)
                  for r in sub["rows"] if r["sex"] == "Female"}
        gains = {f"age {r['age']}": max(0.0, r["incremental_qaly"])
                 for r in sub["rows"] if r["sex"] == "Female"}
        out["distributional"] = distributional_cea(by_age, gains)

    return out


# ══════════════════════════════════════════════════════════════════════════
# Is sequencing worth buying? (prospective, not retrospective)
# ══════════════════════════════════════════════════════════════════════════

def wgs_marginal_value(*, chip_findings_present: bool = True,
                       wgs_only_findings_value: float = 0.0,
                       n_wgs_only_findings: int = 0,
                       wtp: float | None = None,
                       chip_cost: float | None = None,
                       wgs_cost: float | None = None) -> dict:
    """Expected value of sequencing over a genotyping array, before buying it.

    THE PROBLEM THIS SOLVES. The report's ``marginal_chip_to_wgs`` figure sums
    the modelled value of findings tagged ``wgs_only``. Those findings come
    only from a VCF, so on chip input there are none and the figure is $0 —
    every time, for everybody. That is not the finding it looks like. It says
    "your array data contains no sequencing-only findings", which is true by
    construction and answers a question nobody asked.

    The question someone actually asks before paying for sequencing is
    prospective: *what is it likely to find that my array did not?* That
    cannot be answered by counting findings in data that does not exist. It
    can be answered from published yields, which is what this does.

    Three things drive the answer, and the second is the one people miss:

    1. **Secondary findings.** Around 2% of unselected adults carry a
       reportable pathogenic variant in an ACMG actionable gene.
    2. **Arrays barely detect monogenic disease at all.** They type fixed
       positions — the authorised consumer BRCA report covers three founder
       variants out of thousands known — so a reassuring array result is
       close to uninformative. Most of sequencing's value is here, not in
       finding something new about a gene the array already "covered".
    3. **Pharmacogenomics is mostly already covered.** It has by far the
       highest prevalence, so leaving it in undiscounted would make
       sequencing look far more valuable than it is.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    chip_cost = ep.value("cost_chip") if chip_cost is None else float(chip_cost)
    wgs_cost = ep.value("cost_wgs") if wgs_cost is None else float(wgs_cost)

    yield_acmg = ep.value("wgs_yield_acmg_secondary")
    chip_share = ep.value("chip_detection_share_monogenic")
    pgx_cov = ep.value("chip_pgx_coverage")
    mcf = ep.value("marginal_cost_fraction")

    # Value of ONE actionable monogenic finding, acted on. Penetrance is the
    # ascertainment-corrected figure, not the family-series one.
    penetrance = (ep.value("brca2_penetrance_population")
                  * ep.value("ascertainment_shrinkage"))
    rrr = 0.45
    coi = ep.value("coi_pathogenic_generic") * mcf
    qaly = ep.value("qaly_loss_pathogenic_generic")
    horizon = 30.0
    disc = 1.0 / (1.0 + ep.value("discount_rate")) ** (horizon / 2.0)
    value_per_finding = ((penetrance * rrr * coi * disc)
                         + (penetrance * rrr * qaly * wtp * disc)
                         - ep.value("intervention_cost_monogenic"))

    # Incremental detection: the share of pathogenic variants the array misses.
    incremental_yield = yield_acmg * (1.0 - chip_share)
    secondary_value = incremental_yield * max(0.0, value_per_finding)

    # Pharmacogenomics: high prevalence, but mostly already on the array.
    pgx_value_per_person = 250.0        # modest, and mostly captured already
    pgx_incremental = pgx_value_per_person * (1.0 - pgx_cov)

    gross = secondary_value + pgx_incremental
    incremental_cost = wgs_cost - chip_cost
    net = gross - incremental_cost

    nns = (1.0 / incremental_yield) if incremental_yield > 0 else None
    return {
        "available": True,
        "basis": "prospective",
        "expected_incremental_yield": round(incremental_yield, 5),
        "number_needed_to_sequence": (round(nns) if nns else None),
        "value_per_finding": round(value_per_finding),
        "secondary_findings_value": round(secondary_value),
        "pgx_incremental_value": round(pgx_incremental),
        "gross_expected_value": round(gross),
        "incremental_cost": round(incremental_cost),
        "net_expected_value": round(net),
        "worth_it": net > 0,
        "wtp": round(wtp),
        "retrospective_value": round(wgs_only_findings_value),
        "n_wgs_only_findings": int(n_wgs_only_findings),
        "why_retrospective_is_zero": (
            "The report's other sequencing figure counts the value of "
            "findings your data actually contains that only sequencing can "
            "produce. Array data contains none by construction, so that "
            "figure is $0 for everyone with an array — it is a statement "
            "about the input file, not about whether sequencing is worth "
            "buying."
            if n_wgs_only_findings == 0 else ""),
        "plain": (
            f"Sequencing is expected to turn up a serious, actionable "
            f"inherited finding that your array missed in roughly 1 person in "
            f"{round(nns):,}. Averaged over everyone, that is worth about "
            f"${round(gross):,} against an extra "
            f"${round(incremental_cost):,} to buy — so on these numbers it "
            f"{'is' if net > 0 else 'is not'} worth it, "
            f"{'by about $' + format(round(abs(net)), ',') + ' per person.' if net > 0 else 'falling short by about $' + format(round(abs(net)), ',') + ' per person.'}"
            if nns else "Not estimable."),
        "caveat": (
            "This is an average over a population, and the distribution is "
            "extremely skewed: the great majority of people get nothing, and "
            "a small number get something that changes their care "
            "substantially. An expected value is the right basis for a "
            "policy decision and a poor guide to a personal one — if you "
            "would want to know, the average does not capture that."),
        "biggest_driver": (
            "Most of this comes from arrays being nearly blind to monogenic "
            "disease rather than from sequencing finding exotic things. An "
            "array types a fixed list of positions; a normal result for a "
            "gene it 'covers' rules out only the specific variants it looks "
            "at."),
    }
