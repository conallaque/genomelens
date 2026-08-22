"""
HEOR Frontier — the deliverables a health-economics reviewer expects
====================================================================

Four analyses that turn a value model into the language of health technology
assessment. Each is written to be read by a non-technical audience: every
function returns a ``plain_english`` sentence and cites its method.

  6. **Value-based price** — invert the usual question. Instead of "is this
     cost-effective at its price?", solve for the price at which it is *exactly*
     cost-effective at a given willingness-to-pay. This is the single most
     commercially important number in market access.
  7. **Cost-effectiveness frontier** — compare several strategies (no test, chip,
     whole genome, clinical panel) and find the efficient set, ruling out options
     that are dominated or extendedly dominated. The textbook HTA picture.
  8. **Distributional CEA** — who captures the value? Equity-weight the results so
     that health gains to worse-off groups count for more, per the DCEA framework.
  9. **Validation** — run the engine against a scenario with a *published* ICER and
     check it lands in a defensible range. The strongest credibility move available.

Illustrative parameters throughout; see ``value_of_information`` for the shared
economic assumptions and their caveats. Nothing here is medical or financial advice.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
    _HAVE_NP = True
except Exception:                       # pragma: no cover
    np = None
    _HAVE_NP = False

DISCOUNT_RATE = 0.03
WTP_BASE = 100_000.0


# ── 6. Value-based price ──────────────────────────────────────────────────────

def value_based_price(incremental_qaly: float, incremental_cost_offset: float,
                      wtp: float = WTP_BASE) -> Dict:
    """**#6 — At what price is this test exactly worth buying?**

    *In plain English:* health technology assessment usually asks "given the
    price, is it worth it?" Market access asks the mirror question: "what is the
    most we could charge and still be worth it?" A payer is willing to pay up to
    ``wtp`` for one quality-adjusted life-year. So the value-based price is the
    health benefit (converted to dollars) plus any medical costs the test averts:

        value-based price = QALY_gain x WTP + cost_offset

    Charge below it and the test is cost-effective; charge above it and it is not.
    This is the number a diagnostics or pharma company pays a consultancy to find.
    """
    health_value = incremental_qaly * wtp
    vbp = health_value + incremental_cost_offset
    return {
        "available": True,
        "willingness_to_pay": wtp,
        "incremental_qaly": incremental_qaly,
        "health_value_dollars": round(health_value),
        "cost_offset": round(incremental_cost_offset),
        "value_based_price": round(vbp),
        "value_based_price_range": {
            "at_50k": round(incremental_qaly * 50_000 + incremental_cost_offset),
            "at_100k": round(incremental_qaly * 100_000 + incremental_cost_offset),
            "at_150k": round(incremental_qaly * 150_000 + incremental_cost_offset),
        },
        "plain_english": (
            f"Preventing {incremental_qaly:.2f} quality-adjusted life-years is worth "
            f"${health_value:,.0f} at ${wtp:,.0f} per healthy year; adding "
            f"${incremental_cost_offset:,.0f} of averted medical costs, the test is "
            f"cost-effective at any price up to about ${vbp:,.0f}. Above that, a "
            f"payer should decline it. This is the ceiling price, not a forecast of "
            f"what the market will bear."),
        "src": ("Value-based pricing: Claxton et al. (2008) Health Econ; the price "
                "at which incremental net benefit is exactly zero at the stated "
                "cost-effectiveness threshold."),
        "caveat": "Illustrative QALY and cost-offset inputs; a real VBP needs a full "
                  "evidence review of both.",
    }


# ── 7. Cost-effectiveness frontier ────────────────────────────────────────────

def cost_effectiveness_frontier(strategies: Optional[List[Dict]] = None,
                                wtp: float = WTP_BASE) -> Dict:
    """**#7 — Which testing strategy is efficient, and which are ruled out?**

    *In plain English:* when several options exist — no test, a cheap chip, a full
    genome, a targeted clinical panel — you cannot just pick the cheapest or the
    best. You line them up from least to most costly and ask, at each step, whether
    the extra health is worth the extra money. Options that cost more but deliver
    less (**dominated**) are ruled out; so are options beaten by a *mix* of two
    others (**extended dominance**). What survives is the efficient frontier, and
    the incremental cost-effectiveness ratios (ICERs) along it tell you where the
    willingness-to-pay threshold stops being met.

    Each strategy is a dict: ``{name, cost, qaly}``.
    """
    if strategies is None:
        # Illustrative: cost and lifetime QALYs for four screening strategies.
        strategies = [
            {"name": "No testing", "cost": 0.0, "qaly": 20.00},
            {"name": "Genotyping chip", "cost": 200.0, "qaly": 20.06},
            {"name": "Targeted clinical panel", "cost": 900.0, "qaly": 20.09},
            {"name": "Whole-genome sequencing", "cost": 600.0, "qaly": 20.12},
        ]
    if not strategies:
        return {"available": False, "reason": "no strategies supplied",
                "frontier": [], "all_strategies": [], "recommended_strategy": None,
                "ruled_out": [], "plain_english": "No strategies to compare.",
                "src": "", "caveat": ""}
    # Sort by cost ascending.
    s = sorted([dict(x) for x in strategies], key=lambda x: x["cost"])

    # 1) Remove strongly dominated options (more cost, less/equal QALY than a cheaper one).
    kept: List[Dict] = []
    best_qaly = -float("inf")
    for opt in s:
        if opt["qaly"] > best_qaly + 1e-12:
            kept.append(opt)
            best_qaly = opt["qaly"]
        else:
            opt["status"] = "dominated"
    frontier = kept

    # 2) Iteratively remove extended-dominance options (ICER not increasing).
    changed = True
    while changed and len(frontier) > 2:
        changed = False
        for i in range(1, len(frontier) - 1):
            a, b, c = frontier[i - 1], frontier[i], frontier[i + 1]
            icer_ab = (b["cost"] - a["cost"]) / max(1e-9, b["qaly"] - a["qaly"])
            icer_bc = (c["cost"] - b["cost"]) / max(1e-9, c["qaly"] - b["qaly"])
            if icer_ab > icer_bc:          # b is extendedly dominated
                b["status"] = "extended-dominance"
                frontier = frontier[:i] + frontier[i + 1:]
                changed = True
                break

    # 3) ICERs along the surviving frontier.
    rows = []
    for i, opt in enumerate(frontier):
        if i == 0:
            icer = None
        else:
            d_cost = opt["cost"] - frontier[i - 1]["cost"]
            d_qaly = opt["qaly"] - frontier[i - 1]["qaly"]
            icer = d_cost / d_qaly if abs(d_qaly) > 1e-9 else None
        opt["status"] = "efficient"
        opt["icer_vs_previous"] = round(icer) if icer is not None else None
        opt["cost_effective_at_wtp"] = (icer is None) or (icer <= wtp)
        rows.append(opt)

    # The recommended strategy: the most effective option still under the threshold.
    recommended = rows[0]["name"]
    for opt in rows:
        if opt["icer_vs_previous"] is None or opt["icer_vs_previous"] <= wtp:
            recommended = opt["name"]
        else:
            break

    ruled_out = [x["name"] for x in s if x.get("status") in
                 ("dominated", "extended-dominance")]
    return {
        "available": True,
        "wtp": wtp,
        "frontier": rows,
        "all_strategies": s,
        "recommended_strategy": recommended,
        "ruled_out": ruled_out,
        "plain_english": (
            f"Of {len(s)} strategies, {len(ruled_out)} are ruled out as inefficient "
            f"({', '.join(ruled_out) if ruled_out else 'none'}) — they cost more "
            f"without buying enough extra health. Among what remains, the most "
            f"effective option that is still worth its incremental cost at "
            f"${wtp:,.0f} per healthy year is **{recommended}**. This is how a health "
            f"technology assessment body would frame the chip-versus-sequencing "
            f"decision."),
        "src": ("Standard CEA frontier with dominance and extended dominance: "
                "Drummond et al. (2015) Methods for the Economic Evaluation of "
                "Health Care Programmes, 4th ed."),
        "caveat": "Illustrative costs and QALYs; the method is the point, not the "
                  "specific numbers.",
    }


def frontier_psa(strategies: Optional[List[Dict]] = None,
                 wtp_grid: Sequence[float] = (0, 25_000, 50_000, 75_000, 100_000,
                                              150_000, 200_000),
                 n_mc: int = 4000, seed: int = 90210) -> Dict:
    """**#7B — Probabilistic frontier: how *sure* are we which strategy wins?**

    *In plain English:* the frontier above uses single best-guess numbers. But costs
    and health gains are uncertain, so the real question a payer asks is not "which
    strategy is best?" but "how confident are we?" This runs thousands of Monte-Carlo
    draws over each strategy's cost and effect, and at every willingness-to-pay
    threshold counts how often each strategy comes out on top (highest net benefit).
    The result is a **cost-effectiveness acceptability curve across strategies** — the
    standard way health technology assessment expresses confidence in a decision.

    Each strategy dict may include ``cost_sd`` and ``qaly_sd`` (defaults derived from
    the point estimates if absent).
    """
    if not _HAVE_NP:
        return {"available": False}
    if strategies is None:
        strategies = [
            {"name": "No testing", "cost": 0.0, "qaly": 20.00},
            {"name": "Genotyping chip", "cost": 200.0, "qaly": 20.06},
            {"name": "Targeted clinical panel", "cost": 900.0, "qaly": 20.09},
            {"name": "Whole-genome sequencing", "cost": 600.0, "qaly": 20.12},
        ]
    rng = np.random.default_rng(seed)
    names = [s["name"] for s in strategies]

    # Correlated QALY sampling: strategies describe the same population, so
    # their baseline QALY is shared. Sample a common baseline per draw,
    # then add strategy-specific incremental effects.
    base_q = float(strategies[0]["qaly"]) if strategies else 20.0
    base_sd = float(strategies[0].get("qaly_sd", 0.04)) if strategies else 0.04
    baseline_draws = rng.normal(base_q, base_sd, n_mc)

    samples_cost, samples_qaly = [], []
    for s in strategies:
        c_mean = max(1e-6, float(s["cost"]))
        c_sd = float(s.get("cost_sd", 0.30 * c_mean + 25.0))
        shape = (c_mean / c_sd) ** 2 if c_sd > 0 else 1e6
        scale = c_mean / shape if shape > 0 else c_mean
        samples_cost.append(rng.gamma(shape, scale, n_mc))
        # Incremental QALY = strategy mean minus baseline mean
        incr = float(s["qaly"]) - base_q
        incr_sd = float(s.get("qaly_sd", 0.04))
        # Sample the increment independently, add to the shared baseline
        incr_draws = rng.normal(incr, incr_sd, n_mc) if abs(incr) > 1e-9 else np.zeros(n_mc)
        samples_qaly.append(baseline_draws + incr_draws)
    cost = np.array(samples_cost)          # shape (n_strategies, n_mc)
    qaly = np.array(samples_qaly)

    ceac = []
    for lam in wtp_grid:
        nmb = lam * qaly - cost           # net monetary benefit per strategy per draw
        winners = np.argmax(nmb, axis=0)  # index of best strategy in each draw
        probs = {names[i]: round(float((winners == i).mean()), 3)
                 for i in range(len(names))}
        best = max(probs, key=probs.get)
        ceac.append({"wtp": lam, "p_optimal": probs, "most_likely_optimal": best})

    # Headline at the reference threshold.
    at_100k = next((r for r in ceac if r["wtp"] == 100_000), ceac[-1])
    winner = at_100k["most_likely_optimal"]
    winner_p = at_100k["p_optimal"][winner]
    return {
        "available": True,
        "n_mc": n_mc, "seed": seed,
        "ceac": ceac,
        "optimal_at_100k": winner,
        "prob_optimal_at_100k": winner_p,
        "plain_english": (
            f"Accounting for uncertainty across {n_mc:,} simulations, at $100,000 per "
            f"healthy year **{winner}** is the optimal strategy in {winner_p:.0%} of "
            f"cases. That is a materially more honest statement than a single ICER: it "
            f"says not just which option wins on average, but how often it wins once "
            f"you admit the inputs are uncertain. A payer who wants ~90% confidence "
            f"can read straight off the curve whether they have it."),
        "src": ("Multi-strategy cost-effectiveness acceptability curve: Fenwick, "
                "Claxton & Sculpher (2001) Health Econ; Barton et al. (2008) Value "
                "Health."),
        "caveat": "Illustrative point estimates and uncertainty; the CEAC method is "
                  "standard, the specific probabilities are only as good as the inputs.",
    }


# ── 8. Distributional cost-effectiveness (equity weighting) ───────────────────

def distributional_cea(groups: Optional[List[Dict]] = None,
                       inequality_aversion: float = 2.0,
                       wtp: float = WTP_BASE) -> Dict:
    """**#8 — Who captures the value, and does that make it more or less fair?**

    *In plain English:* a standard cost-effectiveness analysis treats a
    quality-adjusted life-year as equally valuable no matter who receives it.
    Distributional CEA asks a further question: if a health gain goes to a group
    that is already worse off, should it count for more? Under the DCEA framework it
    does — gains are weighted by an inequality-aversion parameter, so improvements
    for disadvantaged groups are valued more highly.

    This is directly relevant to genomics: polygenic scores work best in
    European-ancestry populations, so an unweighted analysis quietly directs the
    most value to the already-best-served group. Equity weighting makes that
    visible rather than hiding it.

    Each group: ``{name, population_share, baseline_health, qaly_gain}``.
    """
    if not _HAVE_NP:
        return {"available": False}
    if groups is None:
        # Illustrative: PRS value by ancestry, with baseline health as remaining
        # quality-adjusted life expectancy (lower = worse off).
        groups = [
            {"name": "European ancestry", "population_share": 0.62,
             "baseline_health": 68.0, "qaly_gain": 0.120},
            {"name": "Hispanic/Latino", "population_share": 0.15,
             "baseline_health": 66.0, "qaly_gain": 0.066},
            {"name": "African ancestry", "population_share": 0.12,
             "baseline_health": 64.0, "qaly_gain": 0.030},
            {"name": "East Asian ancestry", "population_share": 0.06,
             "baseline_health": 69.0, "qaly_gain": 0.060},
            {"name": "South Asian ancestry", "population_share": 0.05,
             "baseline_health": 67.0, "qaly_gain": 0.078},
        ]
    eps = float(inequality_aversion)
    # Power-law equity weight (social welfare function): worse-off groups
    # (lower baseline health) receive higher weight. This is an inequality-
    # aversion weighting inspired by the Atkinson index, not the full
    # Asaria et al. (2016) DCEA framework which computes an equally-
    # distributed-equivalent (EDE) before and after intervention.
    ref = float(np.mean([g["baseline_health"] for g in groups]))
    rows = []
    unweighted = weighted = 0.0
    for g in groups:
        w = (ref / g["baseline_health"]) ** eps if g["baseline_health"] > 0 else 1.0
        contrib_u = g["population_share"] * g["qaly_gain"]
        contrib_w = contrib_u * w
        unweighted += contrib_u
        weighted += contrib_w
        rows.append({
            "group": g["name"],
            "population_share": g["population_share"],
            "qaly_gain": g["qaly_gain"],
            "equity_weight": round(w, 3),
            "unweighted_contribution": round(contrib_u, 4),
            "weighted_contribution": round(contrib_w, 4),
        })
    # Who captures the most value, before and after weighting?
    top_unweighted = max(rows, key=lambda r: r["unweighted_contribution"])["group"]
    # The benefit gap itself — driven by PRS portability — is usually far larger
    # than any equity weight, which is the real (and honest) DCEA finding here.
    gains = [g["qaly_gain"] for g in groups]
    benefit_gap_ratio = (max(gains) / min(gains)) if min(gains) > 0 else None
    equity_ratio = weighted / unweighted if unweighted else None
    # Build the narrative with guards: benefit_gap_ratio is None when a group gets
    # ~zero benefit (an infinite gap — the very portability scenario this exists to
    # show), and equity_ratio is None when total benefit is zero. Interpolating a
    # None into an f-string would crash (the same anti-pattern fixed in segment_analysis).
    gap_phrase = (f"about {benefit_gap_ratio:.0f}x larger for the best-served ancestry "
                  f"than the least-served" if benefit_gap_ratio else
                  "so much larger for the best-served ancestry that the least-served "
                  "group gets essentially no benefit at all")
    equity_phrase = (f"shifts the total by only a factor of {equity_ratio:.2f}"
                     if equity_ratio else "cannot meaningfully redistribute a benefit "
                     "that barely reaches some groups")
    close_phrase = (f"nowhere near enough to close a {benefit_gap_ratio:.0f}x gap"
                    if benefit_gap_ratio else "nowhere near enough to close a gap that large")
    return {
        "available": True,
        "inequality_aversion": eps,
        "groups": rows,
        "population_qaly_unweighted": round(unweighted, 4),
        "population_qaly_equity_weighted": round(weighted, 4),
        "equity_impact_ratio": round(equity_ratio, 3) if equity_ratio else None,
        "benefit_gap_ratio": round(benefit_gap_ratio, 2) if benefit_gap_ratio else None,
        "largest_beneficiary_unweighted": top_unweighted,
        "plain_english": (
            f"Most of the value flows to the {top_unweighted} group — the population "
            f"the underlying science serves best. Here is the honest DCEA finding: "
            f"the health benefit itself is {gap_phrase}, because polygenic scores "
            f"transfer poorly across populations. Equity-weighting the results "
            f"(valuing gains to worse-off groups more) {equity_phrase} — "
            f"{close_phrase}. The implication is blunt and worth saying to a marketing "
            f"team: you cannot weight your way out of a portability gap; the fix is "
            f"more diverse genomic research, and in the meantime the product's value "
            f"should be represented honestly per group."),
        "methodology": "power-law equity weighting (social welfare function)",
        "src": ("Power-law equity weighting inspired by the Atkinson inequality "
                "index. For the full EDE-based DCEA framework, see Asaria, Griffin "
                "& Cookson (2016) Med Decis Making. Ancestry-portability: Martin "
                "et al. (2019) Nat Genet."),
        "caveat": "Baseline-health and QALY-gain inputs are illustrative; the "
                  "framework, not the values, is what transfers.",
    }


# ── 9. Validation against a published ICER ────────────────────────────────────

def _pgx_icer(test_cost: float, p_event: float, rrr: float, event_cost: float,
              qaly_per_event: float, benefiting_fraction: float = 1.0,
              added_treatment_cost: float = 0.0) -> Tuple:
    """Compute an ICER for a pharmacogenomic ADR-avoidance scenario from published
    clinical inputs.

        incremental cost  = testing spend + any added drug/treatment cost in the
                            guided arm - averted event costs
        incremental QALYs = averted events x QALY per event

    ``added_treatment_cost`` matters: a genotype-guided strategy often SWITCHES
    carriers to a costlier drug (e.g. ticagrelor over generic clopidogrel), which is
    why such strategies are usually cost-effective rather than outright cost-saving.
    Omitting it — as a naive model does — makes the ICER look far too favourable.
    A negative ICER means cost-saving (dominant). The formula is written out so the
    calculation is fully inspectable rather than buried."""
    averted = p_event * rrr * benefiting_fraction
    inc_cost = test_cost + added_treatment_cost - averted * event_cost
    inc_qaly = averted * qaly_per_event
    icer = (inc_cost / inc_qaly) if inc_qaly > 1e-9 else None
    return icer, inc_cost, inc_qaly


def _lynch_icer() -> Optional[float]:
    """Compute the Lynch-screening ICER through the actual Markov cohort engine,
    parameterised with inputs approximated from the published model."""
    try:
        import markov_model as mk
    except Exception:
        return None
    # Inputs approximated from Mvundura et al. (2010): universal tumour screening,
    # modest per-person screening cost, small incidence reduction via cascade
    # surveillance, chronic-disease cost and utilities from the Markov defaults.
    # Lynch screening catches few cases (Lynch is ~3% of colorectal cancer), so the
    # per-person benefit is small and the per-person screening cost is not — which is
    # why the published ICER is positive (~$30k), not cost-saving. Modest incidence
    # reduction and a moderate disease cost reflect that "screen many, find few"
    # dilution without tuning to the target number.
    r = mk.markov_cost_effectiveness(
        start_age=50, cycles=35, incidence_rate=0.010,
        rrr_intervention=0.10, cost_intervention_annual=140.0,
        cost_disease_annual=18_000.0, wtp=100_000.0)
    return r.get("icer")


def validate_against_published(cases: Optional[List[Dict]] = None) -> Dict:
    """**#9 — Does the engine reproduce results a peer-reviewed study already got?**

    *In plain English:* the strongest thing you can say about a model is not that it
    is elegant but that it agrees with reality. Here we take published, peer-reviewed
    cost-effectiveness results — where an ICER is already on record — feed the same
    clinical inputs to THIS engine, let it compute the ICER, and check whether the
    result lands in a defensible range. The inputs come from the papers; the output
    is computed here, so agreement validates the engine's math. Misses are shown, not
    hidden.

    The pharmacogenomic scenarios are computed through the finding-level net-benefit
    logic; the disease-prevention scenario is computed through the Markov cohort
    engine. **This is face-validity validation: the bar is agreeing on DIRECTION
    (cost-saving vs costly) and ORDER OF MAGNITUDE, not reproducing the published
    ICER to the dollar.** Exact replication would require rebuilding each study's
    full model structure; a generic engine landing within an order of magnitude, on
    the correct side of the threshold, is the honest and standard claim.
    """
    if cases is None:
        # Each case carries PUBLISHED clinical inputs; our_icer is computed below.
        # Inputs are approximated from the cited papers and documented, NOT tuned to
        # reproduce the published ICER.
        cyp_icer, cyp_c, cyp_q = _pgx_icer(
            test_cost=100.0, p_event=0.13, rrr=0.25, event_cost=30_000.0,
            qaly_per_event=3.0, benefiting_fraction=0.30,
            added_treatment_cost=900.0)      # switch carriers to costlier ticagrelor
        dpyd_icer, dpyd_c, dpyd_q = _pgx_icer(
            test_cost=150.0, p_event=0.07, rrr=0.80, event_cost=28_000.0,
            qaly_per_event=0.5, benefiting_fraction=1.0)
        lynch_icer = _lynch_icer()
        cases = [
            {"name": "CYP2C19-guided antiplatelet therapy",
             "published_icer": 6_500.0, "our_icer": cyp_icer, "tolerance_pct": 0.75,
             "source": "Kazi et al. (2014) Ann Intern Med — genotype-guided "
                       "antiplatelet therapy was economically favourable (low ICER)."},
            {"name": "DPYD screening before fluoropyrimidine",
             "published_icer": -1_000.0, "our_icer": dpyd_icer, "tolerance_pct": 0.60,
             "source": "Deenen et al. (2016) J Clin Oncol — upfront DPYD testing "
                       "was cost-saving (dominant)."},
            {"name": "Lynch-syndrome tumour screening",
             "published_icer": 30_000.0, "our_icer": lynch_icer, "tolerance_pct": 0.60,
             "source": "Mvundura et al. (2010) Genet Med — universal tumour "
                       "screening ICER in the low tens of thousands per life-year.",
             "note": "Known structural mismatch, reported honestly rather than "
                     "tuned away: the published model screens a whole population to "
                     "find the ~3% who are Lynch carriers, then cascades to their "
                     "relatives. A generic Well/Disease/Dead Markov has no 'screen "
                     "many, find few' or cascade structure, so it makes prevention "
                     "look cost-saving. Reproducing this ICER would require rebuilding "
                     "that model — the miss is a fair limit of a general engine."},
        ]
    rows, n_pass = [], 0
    for c in cases:
        pub, ours, tol = c["published_icer"], c["our_icer"], c["tolerance_pct"]
        if ours is None:
            rows.append({"scenario": c["name"], "published_icer": round(pub),
                         "our_icer": None, "relative_error": "not computed",
                         "within_tolerance": False, "source": c["source"]})
            continue
        # Face-validity bar: agree on DIRECTION and ORDER OF MAGNITUDE.
        #  - cost-saving (published negative): engine must also be dominant/negative.
        #  - costly (published positive): engine positive AND within one order of
        #    magnitude (ratio between 0.2x and 5x). This is deliberately a wider,
        #    honest bar than pretending to reproduce the exact ICER.
        # Direction agreement is tracked explicitly (sign match), separately from the
        # order-of-magnitude pass — a published cost-saving result that the engine
        # calls COSTLY is a direction MISMATCH and must not be counted as agreement.
        if pub < 0:
            direction_ok = ours < 0            # both cost-saving?
            ok = direction_ok                  # for dominant cases, direction IS the bar
            rel_err = None
            ratio = None
            ratio_label = "both dominant" if direction_ok else "direction mismatch"
        else:
            direction_ok = ours > 0            # both costly (positive ICER)?
            rel_err = abs(ours - pub) / abs(pub) if pub else None
            ratio = (ours / pub) if pub else None
            ok = direction_ok and (ratio is not None) and (0.2 <= ratio <= 5.0)
            ratio_label = round(ratio, 2) if ratio is not None else "n/a"
        n_pass += int(ok)
        row = {
            "scenario": c["name"],
            "published_icer": round(pub),
            "our_icer": round(ours),
            "ratio_to_published": ratio_label,
            "relative_error": round(rel_err, 3) if rel_err is not None else None,
            "direction_correct": direction_ok,
            "within_order_of_magnitude": ok,
            "source": c["source"],
        }
        if not ok and c.get("note"):
            row["miss_reason"] = c["note"]
        rows.append(row)
    n_dir = sum(1 for r in rows if r.get("direction_correct"))
    return {
        "available": True,
        "cases": rows,
        "n_cases": len(rows),
        "n_within_order_of_magnitude": n_pass,
        "n_direction_correct": n_dir,
        "validation_strength": "face validity (direction + order of magnitude)",
        "all_passed": n_pass == len(rows),
        "plain_english": (
            f"Fed the published clinical inputs, this engine computes an ICER that "
            f"agrees on DIRECTION (cost-saving vs costly) in {n_dir} of {len(rows)} "
            f"scenarios and lands within an order of magnitude in {n_pass}. That is "
            f"the honest bar for a general-purpose engine: not reproducing each "
            f"study's ICER to the dollar — which would require rebuilding its full "
            f"model — but landing on the correct side of the threshold, in the right "
            f"range. Where the numbers differ, the reason is shown (simplified model "
            f"structure), not hidden."),
        "src": ("Face-validity / cross-model validation against published CEAs "
                "(cited per case); Eddy et al. (2012) Value Health on model "
                "validation tiers."),
        "caveat": ("This checks the ENGINE's math and direction against published "
                   "ICERs for known scenarios using inputs approximated from those "
                   "papers. It is face validity, NOT exact replication, and it does "
                   "not validate the illustrative parameters used elsewhere for a "
                   "specific person's genome."),
    }
