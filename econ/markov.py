"""
Markov State-Transition Cohort Model + Budget Impact Analysis
=============================================================

The two canonical HEOR deliverables that a static decision model cannot produce.

**1. Markov cohort model.** Health economics evaluates interventions whose effects
play out over decades, with people moving between health states and dying along the
way. The standard structure is a discrete-time state-transition ("Markov") cohort
model: a closed cohort is distributed across mutually exclusive health states, a
transition matrix moves them each cycle, and costs/QALYs accrue per state-cycle.
Implemented here with the conventions an HTA reviewer checks for:

  * annual cycles with an explicit **half-cycle correction** (trapezoidal),
  * **rate → probability** conversion  p = 1 − exp(−r·Δt)  (never r·Δt),
  * **age-dependent background mortality** (Gompertz), so competing death is real,
  * a validated transition matrix (rows sum to 1; cohort conserved every cycle),
  * discounting of **both** costs and QALYs at 3%,
  * incremental analysis → ΔCost, ΔQALY, **ICER**, and NMB.

**2. Budget impact analysis (BIA).** CEA answers "is it worth it?" (efficiency);
BIA answers "can we afford it?" (affordability) — the question a payer actually
asks. Different conventions by design (ISPOR Task Force): short horizon (1–5 years),
population-scaled to a real plan, **undiscounted**, uptake phased in over time, and
reported as **per-member-per-month (PMPM)** — the number that decides formulary
placement.

Both are illustrative teaching-grade models with transparent parameters, not
submissions. Every assumption is a named argument you can change and re-run.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

try:
    import numpy as np
    _HAVE_NP = True
except Exception:                       # pragma: no cover
    np = None
    _HAVE_NP = False

DISCOUNT_RATE = 0.03
STATES = ("Well", "Disease", "Dead")


# ── helpers ───────────────────────────────────────────────────────────────────

def rate_to_prob(rate: float, cycle_length: float = 1.0) -> float:
    """Convert a continuous-time rate to a per-cycle probability.

    p = 1 − exp(−r·Δt).  Using r·Δt directly is the most common error in applied
    modeling: it is only accurate for small r, and can exceed 1 for large r.
    """
    return 1.0 - math.exp(-max(0.0, float(rate)) * float(cycle_length))


def _gompertz_mortality(age: float, a: float = 0.0001, b: float = 0.085) -> float:
    """Age-specific all-cause mortality *rate* (Gompertz)."""
    return a * math.exp(b * max(0.0, age - 20.0))


def build_transition_matrix(p_well_to_disease: float, p_disease_death_excess: float,
                            p_background_death: float) -> list[list[float]]:
    """Assemble a validated 3×3 transition matrix for (Well, Disease, Dead).

    Competing risks are handled by construction: from **Well** you may develop the
    disease *or* die of something else, and the residual stays well. Rows are
    normalized so each sums to exactly 1.
    """
    p_wd = max(0.0, min(1.0, p_well_to_disease))
    p_bg = max(0.0, min(1.0, p_background_death))
    # Competing risks from Well: scale so the two exits cannot exceed 1.
    total_exit = p_wd + p_bg
    if total_exit > 1.0:
        p_wd, p_bg = p_wd / total_exit, p_bg / total_exit
    row_well = [1.0 - p_wd - p_bg, p_wd, p_bg]

    p_dd = max(0.0, min(1.0, p_background_death + p_disease_death_excess))
    row_dis = [0.0, 1.0 - p_dd, p_dd]

    row_dead = [0.0, 0.0, 1.0]                    # absorbing
    return [row_well, row_dis, row_dead]


# ── 1. Markov cohort model ────────────────────────────────────────────────────

def run_markov(strategy: str = "standard_care",
               start_age: float = 40.0, cycles: int = 45,
               incidence_rate: float = 0.010,
               rrr_intervention: float = 0.30,
               cost_intervention_annual: float = 250.0,
               cost_disease_annual: float = 12_000.0,
               cost_well_annual: float = 0.0,
               utility_well: float = 0.90, utility_disease: float = 0.68,
               excess_mortality_rate: float = 0.035,
               discount_rate: float = DISCOUNT_RATE,
               half_cycle: bool = True) -> dict:
    """Run one strategy through the cohort model and return the trace + totals.

    ``strategy='genomic_guided'`` applies ``rrr_intervention`` to the Well→Disease
    transition (prevention/early screening triggered by a genomic finding) and pays
    ``cost_intervention_annual`` while the patient remains well.
    """
    guided = strategy == "genomic_guided"
    cohort = [1.0, 0.0, 0.0]                       # everyone starts Well
    trace: list[dict] = []
    tot_cost = tot_qaly = tot_ly = 0.0

    for t in range(cycles):
        age = start_age + t
        # Effective incidence, reduced by the intervention if in the guided arm.
        inc = incidence_rate * (1.0 - rrr_intervention) if guided else incidence_rate
        p_wd = rate_to_prob(inc)
        p_bg = rate_to_prob(_gompertz_mortality(age))
        p_ex = rate_to_prob(excess_mortality_rate)
        P = build_transition_matrix(p_wd, p_ex, p_bg)

        # Per-cycle payoffs, valued on the state occupancy at the START of the cycle.
        c_cycle = (cohort[0] * (cost_well_annual + (cost_intervention_annual if guided else 0.0))
                   + cohort[1] * cost_disease_annual)
        q_cycle = cohort[0] * utility_well + cohort[1] * utility_disease
        l_cycle = cohort[0] + cohort[1]

        # Half-cycle correction: weight the first and last cycles at one half
        # (trapezoidal rule) — states are entered continuously, not at cycle bounds.
        w = 0.5 if (half_cycle and (t == 0 or t == cycles - 1)) else 1.0
        disc = 1.0 / ((1.0 + discount_rate) ** t)
        tot_cost += c_cycle * w * disc
        tot_qaly += q_cycle * w * disc
        tot_ly += l_cycle * w * disc

        trace.append({"cycle": t, "age": round(age, 1),
                      "well": round(cohort[0], 5), "disease": round(cohort[1], 5),
                      "dead": round(cohort[2], 5)})

        # Advance the cohort: new_state_j = Σ_i cohort_i · P[i][j]
        cohort = [sum(cohort[i] * P[i][j] for i in range(3)) for j in range(3)]
        total = sum(cohort)
        if abs(total - 1.0) > 1e-9:                # guard against leakage
            cohort = [x / total for x in cohort]

    return {
        "strategy": strategy,
        "cycles": cycles, "start_age": start_age,
        "total_cost": round(tot_cost, 2),
        "total_qaly": round(tot_qaly, 4),
        "total_life_years": round(tot_ly, 4),
        "final_distribution": {"well": round(cohort[0], 5),
                               "disease": round(cohort[1], 5),
                               "dead": round(cohort[2], 5)},
        # Unrounded totals, for the same reason cohort_sum_exact is unrounded:
        # display precision must not leak into arithmetic. Differencing the
        # rounded totals above and then multiplying the QALY delta by a
        # $100,000 threshold turns a 4dp rounding into ~$5 of net monetary
        # benefit — small, but it is noise reported as signal, and an HTA
        # reviewer recomputing NMB by hand will not reproduce the figure.
        # markov_cost_effectiveness differences THESE and rounds for display.
        "total_cost_exact": tot_cost,
        "total_qaly_exact": tot_qaly,
        "total_life_years_exact": tot_ly,
        # Unrounded sum, so validation tests the model rather than the display
        # precision (three values rounded to 5dp can sum to 1.00001).
        "cohort_sum_exact": sum(cohort),
        "trace": trace[::max(1, cycles // 15)],
        "half_cycle_correction": half_cycle,
        "discount_rate": discount_rate,
    }


def markov_cost_effectiveness(wtp: float = 100_000.0, **kwargs) -> dict:
    """Run both strategies and produce the incremental analysis (ΔC, ΔQ, ICER, NMB).

    This is the deliverable an HTA body reads: two arms, an incremental
    cost-effectiveness ratio, and the net monetary benefit at a stated threshold.
    """
    sc = run_markov(strategy="standard_care", **kwargs)
    gg = run_markov(strategy="genomic_guided", **kwargs)

    # Difference the UNROUNDED totals. Taking these from the rounded display
    # values (2dp cost, 4dp QALY) put the rounding inside the arithmetic: NMB
    # multiplies the QALY delta by the willingness-to-pay threshold, so a 5e-5
    # rounding became ~$5 at $100,000/QALY. The reported incrementals below are
    # still rounded — rounding for display is fine, rounding before the
    # multiply is not.
    d_cost = gg["total_cost_exact"] - sc["total_cost_exact"]
    d_qaly = gg["total_qaly_exact"] - sc["total_qaly_exact"]
    d_ly = gg["total_life_years_exact"] - sc["total_life_years_exact"]
    # An ICER is only meaningful when cost and effect move in the SAME
    # direction. In the dominance quadrants the ratio is negative, and a
    # negative ICER is ambiguous by construction — "-$6,054/QALY" reads as a
    # bargain whether the strategy saves money and adds health or costs money
    # and destroys it. HTA reporting convention is to state dominance and
    # suppress the ratio. The verdict below already did this correctly; the
    # ratio was reported alongside it anyway.
    icer = (d_cost / d_qaly) if abs(d_qaly) > 1e-9 else None
    dominant = d_cost < 0 and d_qaly > 0
    dominated = d_cost > 0 and d_qaly < 0
    if dominant or dominated:
        icer = None
    nmb = wtp * d_qaly - d_cost

    # Dominance language matters in HTA reporting.
    if dominant:
        verdict = "dominant (cheaper and more effective)"
    elif dominated:
        verdict = "dominated (costlier and less effective)"
    elif icer is not None and icer <= wtp:
        verdict = f"cost-effective at ${wtp:,.0f}/QALY"
    else:
        verdict = f"not cost-effective at ${wtp:,.0f}/QALY"

    return {
        "available": True,
        # Dominance is exposed as a flag, not left to be read off the sign of
        # a suppressed ratio. Downstream validation used to infer "cost-saving"
        # from a negative ICER; once the ratio is correctly withheld in the
        # dominance quadrants that inference has nothing to read, so the fact
        # has to be stated directly.
        "dominant": dominant,
        "dominated": dominated,
        "cost_saving": d_cost < 0,
        "standard_care": sc, "genomic_guided": gg,
        "incremental_cost": round(d_cost, 2),
        "incremental_qaly": round(d_qaly, 4),
        "incremental_life_years": round(d_ly, 4),
        "icer": round(icer, 2) if icer is not None else None,
        "nmb_at_wtp": round(nmb, 2),
        "wtp": wtp,
        "verdict": verdict,
        "note": ("Discrete-time Markov cohort model with annual cycles, half-cycle "
                 "correction, rate-to-probability conversion, age-dependent competing "
                 "mortality, and 3% discounting of both costs and QALYs."),
        "src": ("Sonnenberg & Beck (1993), Med Decis Making; Briggs, Sculpher & "
                "Claxton (2006), Decision Modeling for Health Economic Evaluation"),
    }


def validate_markov(result: dict) -> dict:
    """Structural validation an HTA reviewer would run: cohort conservation, an
    absorbing death state, and monotone accumulation of deaths."""
    checks, ok = [], True

    def add(name, passed, detail=""):
        nonlocal ok
        ok = ok and passed
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    for arm in ("standard_care", "genomic_guided"):
        r = result.get(arm) or {}
        total = r.get("cohort_sum_exact")
        if total is None:                          # fall back to the rounded display
            fd = r.get("final_distribution", {})
            total = sum(fd.values()) if fd else 0.0
        add(f"{arm}: cohort conserved (Σ states = 1)", abs(total - 1.0) < 1e-9,
            f"Σ = {total:.12f}")
        deaths = [row["dead"] for row in (r.get("trace") or [])]
        add(f"{arm}: death is absorbing (non-decreasing)",
            all(deaths[i] <= deaths[i + 1] + 1e-9 for i in range(len(deaths) - 1)),
            f"{deaths[0]:.4f} → {deaths[-1]:.4f}" if deaths else "")
    add("intervention reduces incidence (fewer diseased at end)",
        (result["genomic_guided"]["final_distribution"]["disease"]
         <= result["standard_care"]["final_distribution"]["disease"] + 1e-9))
    return {"all_passed": ok, "checks": checks}


# ── 2. Budget Impact Analysis ─────────────────────────────────────────────────

def budget_impact(plan_members: int = 1_000_000,
                  eligible_fraction: float = 0.08,
                  test_cost: float = 300.0,
                  uptake_curve: Sequence[float] = (0.05, 0.12, 0.20, 0.28, 0.35),
                  actionable_fraction: float = 0.20,
                  annual_offset_per_actionable: float = 900.0,
                  annual_intervention_cost: float = 250.0,
                  horizon_years: int = 5) -> dict:
    """**Budget impact analysis** — the payer's affordability question.

    Deliberately *not* a CEA. Per the ISPOR BIA Task Force the conventions differ:

      * short horizon (1–5 years), because budgets are annual;
      * scaled to a **real plan population**, not a hypothetical cohort;
      * **undiscounted** (these are actual cash outlays in the budget year);
      * **uptake phased in** rather than assuming instant 100% adoption;
      * headline metric is **per-member-per-month (PMPM)**, which is how formulary
        and coverage decisions are actually argued.

    Costs = testing + downstream preventive intervention for those who act.
    Offsets = averted disease costs among those managed earlier.
    """
    eligible = plan_members * eligible_fraction
    rows = []
    cum_net = 0.0
    for y in range(1, horizon_years + 1):
        uptake = uptake_curve[min(y - 1, len(uptake_curve) - 1)]
        tested_new = eligible * (uptake - (uptake_curve[y - 2] if y > 1 else 0.0))
        tested_cum = eligible * uptake
        actionable_cum = tested_cum * actionable_fraction

        cost_testing = tested_new * test_cost          # one-off per person tested
        cost_intervention = actionable_cum * annual_intervention_cost
        # Offsets ramp: benefits accrue only after the intervention has been running.
        offsets = actionable_cum * annual_offset_per_actionable * min(1.0, (y - 1) / 2.0)
        net = cost_testing + cost_intervention - offsets
        cum_net += net
        rows.append({
            "year": y,
            "uptake": round(uptake, 4),
            "tested_cumulative": round(tested_cum),
            "actionable_cumulative": round(actionable_cum),
            "cost_testing": round(cost_testing),
            "cost_intervention": round(cost_intervention),
            "offsets": round(offsets),
            "net_budget_impact": round(net),
            "pmpm": round(net / (plan_members * 12.0), 4),
            "cumulative_net": round(cum_net),
        })

    peak = max(rows, key=lambda r: r["net_budget_impact"])
    final = rows[-1]
    return {
        "available": True,
        "plan_members": plan_members,
        "eligible_population": round(eligible),
        "horizon_years": horizon_years,
        "rows": rows,
        # "Peak" alone does not say peak *what*, and the figure is not the
        # largest number in the series — once the program turns cost-saving
        # the last year is bigger in magnitude and smaller in burden. The name
        # states which: the worst year for the payer, selected on signed net
        # spend per ISPOR convention. Old keys retained; renderers migrate.
        "maximum_budget_burden_year": peak["year"],
        "maximum_budget_burden_net": peak["net_budget_impact"],
        "maximum_budget_burden_pmpm": peak["pmpm"],
        "peak_year": peak["year"],
        "peak_net_budget_impact": peak["net_budget_impact"],
        "peak_pmpm": peak["pmpm"],
        "year5_net": final["net_budget_impact"],
        "year5_pmpm": final["pmpm"],
        "cumulative_net": final["cumulative_net"],
        "becomes_cost_saving": any(r["net_budget_impact"] < 0 for r in rows),
        "note": ("Undiscounted by BIA convention (these are budget-year cash flows). "
                 "PMPM is the decision metric: payers typically treat well under "
                 "$1.00 PMPM as easily absorbable."),
        "src": "Sullivan et al. (2014), Value in Health — ISPOR BIA Task Force",
    }
