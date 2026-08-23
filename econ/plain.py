"""
Plain-language health economics: the same results, said in English
===================================================================

The rest of the economic layer is defensible and almost unreadable. It reports
incremental net monetary benefit, quality-adjusted life-years, acceptability
curves and partial expected value of perfect information — all correct, all
standard, and none of it meaning anything to someone who has not done a health
economics course. A result nobody can read is not a result.

This module does two things.

**It adds outcomes that are intuitive by construction.** Number needed to
screen is the plainest statistic in medicine — "about 1 in 40 people who do
this avoid a heart attack" needs no glossary. Quality-adjusted life-years
converted to months of healthy life are the same quantity a person can picture.
Payback period answers "when does this pay for itself" directly. These are not
simplifications of the technical results; they are standard health-economic
outputs that happen to be legible, and the model had the ingredients for all of
them without reporting any.

**It translates what is already there.** A 97% probability of cost-effectiveness
becomes "in 97 out of 100 versions of this model, doing it was worth it". An
ICER becomes "each extra healthy year costs about $X, against the $100,000 a
year of good health is usually valued at". A tornado becomes "the answer
depends most on a number nobody has measured".

Translation is where over-claiming creeps in, so two rules hold throughout:
never state a plain-language conclusion the technical result does not support,
and keep the uncertainty attached to the sentence rather than in a footnote.
Where the honest answer is "we don't know", these functions say so.
"""

from __future__ import annotations

from collections.abc import Sequence

from . import params as ep

__all__ = [
    "CONDITION_NAMES",
    "build_plain_summary",
    "healthy_time_gained",
    "number_needed_to_screen",
    "payback_period",
    "plain_actions",
    "plain_confidence",
    "plain_money",
    "plain_verdict",
]

# Human names for the internal condition keys. "coi_key: CAD" is not English.
CONDITION_NAMES: dict[str, str] = {
    "CAD": "heart attack or stroke",
    "T2D": "type 2 diabetes",
    "Alzheimer": "dementia",
    "Depression": "depression",
    "SubstanceUse": "a substance-use problem",
    "Autoimmune": "an autoimmune condition",
    "Urologic": "kidney stones or related urinary problems",
    "IronOverload": "iron overload",
    "Colorectal": "bowel cancer",
    "BreastOvarian": "breast or ovarian cancer",
    "Pathogenic": "a serious inherited condition",
}


def _name(coi_key: str) -> str:
    return CONDITION_NAMES.get(coi_key, coi_key or "a health problem")


# ══════════════════════════════════════════════════════════════════════════
# Outcomes that are intuitive by construction
# ══════════════════════════════════════════════════════════════════════════

def number_needed_to_screen(conditions: Sequence[dict],
                            *, assumed_baseline: bool = True) -> list[dict]:
    """"About 1 in N people who do this avoid X."

    Number needed to screen is the inverse of the absolute risk reduction, and
    it is the single most legible statistic health economics produces. The
    model already computes ``cases_averted`` per condition, which *is* the
    absolute risk reduction — it simply never inverted it.

    A large N is not a failure and should not be hidden. "1 in 900" is the
    honest description of most population screening, and a reader who sees it
    understands the trade-off better than one shown only a dollar figure.
    """
    out: list[dict] = []
    for c in conditions:
        arr = float(c.get("cases_averted", 0.0) or 0.0)
        if arr <= 0:
            continue
        nns = 1.0 / arr
        cond = _name(c.get("condition", ""))
        # MOOD MATTERS HERE. The baseline risk behind this number is, for most
        # findings, a generic registered assumption rather than a measured
        # risk for this person. Writing "1 in 17 people avoid a serious
        # inherited condition" in the indicative turns that assumption into a
        # claim — which is precisely how a plain-language layer can end up
        # less honest than the technical output it summarises. The
        # conditional mood costs nothing and keeps the sentence true.
        if assumed_baseline:
            plain = (f"If the model's assumptions hold, roughly 1 in "
                     f"{round(nns):,} people with this pattern of results "
                     f"would avoid {cond} by acting on it.")
        else:
            plain = (f"About 1 in {round(nns):,} people like you who act "
                     f"on this avoid {cond} that they would otherwise "
                     f"have had.")
        out.append({
            "condition": c.get("condition", ""),
            "condition_plain": cond,
            "absolute_risk_reduction": round(arr, 5),
            "number_needed_to_screen": round(nns),
            "baseline_is_assumed": bool(assumed_baseline),
            "plain": plain,
        })
    out.sort(key=lambda r: r["number_needed_to_screen"])
    return out


def healthy_time_gained(qalys: float) -> dict:
    """Quality-adjusted life-years, said as time a person can picture.

    A QALY is a year of life in full health; fractions of one are the usual
    output and are almost impossible to feel. Months and days are the same
    number in units people use.
    """
    q = float(qalys or 0.0)
    days = q * 365.25
    months = q * 12.0
    if abs(days) < 1:
        plain = ("Less than a day of extra healthy life, on average — the "
                 "health gain here is essentially a rounding error.")
    elif abs(days) < 60:
        plain = (f"About {days:.0f} extra days of healthy life, on average.")
    elif abs(months) < 24:
        plain = (f"About {months:.1f} extra months of healthy life, "
                 f"on average.")
    else:
        plain = (f"About {q:.1f} extra years of healthy life, on average.")
    return {
        "qalys": round(q, 4),
        "days": round(days, 1),
        "months": round(months, 2),
        "plain": plain,
        "caveat": ("An average across many people like you, not a promise to "
                   "you. Most people get nothing from any single preventive "
                   "action; a few get a great deal, and the average spreads "
                   "that across everyone."),
    }


def payback_period(*, upfront_cost: float, annual_saving: float,
                   max_years: int = 40) -> dict:
    """"It pays for itself after about N years" — or honestly, that it doesn't.

    Undiscounted deliberately: this answers a cash question in the terms
    someone actually asks it. The discounted version lives in the technical
    result and disagreeing slightly is expected.
    """
    cost = float(upfront_cost or 0.0)
    saving = float(annual_saving or 0.0)
    if cost <= 0:
        return {"pays_back": True, "years": 0.0,
                "plain": "There is no upfront cost to recover."}
    if saving <= 0:
        return {"pays_back": False, "years": None,
                "plain": (f"This does not pay for itself in money. It costs "
                          f"about {plain_money(cost)} and does not reduce "
                          f"medical spending by more than that — you are "
                          f"buying health, not savings.")}
    yrs = cost / saving
    if yrs > max_years:
        return {"pays_back": False, "years": round(yrs, 1),
                "plain": (f"On these numbers it would take about {yrs:.0f} "
                          f"years to pay for itself, which is long enough "
                          f"that it should be treated as not paying back.")}
    return {"pays_back": True, "years": round(yrs, 1),
            "plain": (f"The money spent comes back after about "
                      f"{yrs:.1f} years in medical costs you avoid.")}


# ══════════════════════════════════════════════════════════════════════════
# Translation
# ══════════════════════════════════════════════════════════════════════════

def plain_money(v) -> str:
    """Dollar amounts at a readable precision."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "an unknown amount"
    a = abs(v)
    if a >= 1_000_000:
        s = f"${a/1_000_000:.1f} million"
    elif a >= 10_000:
        s = f"${a/1_000:.0f},000".replace(",000", ",000")
        s = f"${round(a/1_000):,},000".replace(",000,000", " million")
    elif a >= 1_000:
        s = f"${a:,.0f}"
    else:
        s = f"${a:,.0f}"
    return f"minus {s}" if v < 0 else s


def plain_verdict(cea: dict, *, wtp: float | None = None) -> dict:
    """One sentence answering "is this worth doing?".

    The three cases a reader needs distinguished, and which the technical
    output blurs: it saves money, it costs money but buys health cheaply
    enough to be worth it, or it costs more than the health is usually valued
    at. Prevention is usually the middle one, and calling that "cost-saving"
    is the most common way these reports mislead.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    cost = float(cea.get("incremental_cost", 0) or 0)
    qaly = float(cea.get("incremental_qaly", 0) or 0)
    icer = cea.get("icer")

    if qaly <= 0 and cost > 0:
        return {"headline": "Not worth doing on these numbers",
                "plain": ("This costs money and the model does not find a "
                          "health gain to justify it."),
                "tone": "negative"}
    if cost < 0 and qaly > 0:
        return {"headline": "Worth doing — and it saves money",
                "plain": (f"Acting on these findings is modelled to improve "
                          f"health and reduce medical spending by about "
                          f"{plain_money(-cost)} over the period. That is "
                          f"unusual: most prevention improves health while "
                          f"costing money."),
                "tone": "positive"}
    if icer is not None and icer <= wtp:
        return {"headline": "Worth doing, but it costs money",
                "plain": (f"Each extra year of healthy life costs about "
                          f"{plain_money(icer)}. Health systems in the US "
                          f"generally treat anything under "
                          f"{plain_money(wtp)} per healthy year as good "
                          f"value, so this clears that bar — but it is a "
                          f"purchase, not a saving."),
                "tone": "positive"}
    if icer is not None:
        return {"headline": "Probably not good value at the usual threshold",
                "plain": (f"Each extra year of healthy life costs about "
                          f"{plain_money(icer)}, above the "
                          f"{plain_money(wtp)} that is usually treated as "
                          f"good value."),
                "tone": "negative"}
    return {"headline": "No clear answer",
            "plain": ("The model does not produce a clear enough result here "
                      "to call it either way."),
            "tone": "neutral"}


def plain_confidence(psa: dict | None) -> dict:
    """"In N out of 100 versions of this model, it was worth it."

    A probability of cost-effectiveness is a frequency, and frequencies are
    understood far better as counts out of a hundred than as percentages or
    as an interval.
    """
    if not psa or not psa.get("available"):
        return {"available": False}
    p = float(psa.get("p_cost_effective", 0) or 0)
    n = round(p * 100)
    lo, hi = psa.get("inmb_ci_low"), psa.get("inmb_ci_high")
    if n >= 95:
        strength = "The conclusion is not close to the line."
    elif n >= 80:
        strength = "The conclusion holds up in most versions, but not all."
    elif n >= 50:
        strength = "This is genuinely uncertain — it could go either way."
    else:
        strength = "The model more often finds this is not worth it."
    return {
        "available": True,
        "n_in_100": n,
        "plain": (f"The model was run {psa.get('n_iterations', 0):,} times, "
                  f"each time with slightly different assumptions drawn from "
                  f"the ranges the research supports. Acting came out worth "
                  f"it in about {n} of every 100 runs. {strength}"),
        "range_plain": (
            f"Across those runs the overall value ranged from about "
            f"{plain_money(lo)} to {plain_money(hi)}."
            if lo is not None and hi is not None else ""),
        "honesty": psa.get("note", ""),
    }


def plain_actions(conditions: Sequence[dict], *, top: int = 5) -> list[dict]:
    """The ranked list of what to actually do, in English.

    The technical table is sorted by net monetary benefit, which is the right
    ordering and an opaque label. This says what each row means.
    """
    rows = sorted(conditions, key=lambda c: c.get("inmb", 0), reverse=True)
    out: list[dict] = []
    for c in rows[:top]:
        if (c.get("inmb", 0) or 0) <= 0:
            continue
        arr = float(c.get("cases_averted", 0.0) or 0.0)
        nns = round(1.0 / arr) if arr > 0 else None
        rrr = float(c.get("combined_rrr", 0.0) or 0.0)
        n_src = int(c.get("n_findings", 0) or 0)
        flagged = ("One of your results flags this."
                   if n_src == 1 else
                   f"{n_src} separate results in your report point at this.")
        out.append({
            "condition": c.get("condition", ""),
            "what": f"Reduce your risk of {_name(c.get('condition', ''))}",
            "why": (f"{flagged} Acting on it is modelled to cut that risk by "
                    f"about {rrr:.0%}."),
            "scale": (f"On the model's assumptions, that is roughly 1 person "
                      f"in {nns:,} avoiding the condition altogether."
                      if nns else ""),
            "value": plain_money(c.get("inmb", 0)),
        })
    return out


# ══════════════════════════════════════════════════════════════════════════
# Assembly
# ══════════════════════════════════════════════════════════════════════════

def build_plain_summary(pooled: dict | None) -> dict:
    """Assemble the plain-language summary from the pooled economic result."""
    if not pooled or not pooled.get("available"):
        return {"available": False}

    cea = pooled.get("cea") or {}
    conds = pooled.get("conditions") or []
    psa = pooled.get("psa") or {}
    wtp = float(cea.get("wtp") or ep.value("wtp_per_qaly"))

    verdict = plain_verdict(cea, wtp=wtp)
    time_gain = healthy_time_gained(cea.get("incremental_qaly", 0))
    nns = number_needed_to_screen(conds)
    conf = plain_confidence(psa)

    cost_averted = float(cea.get("cost_averted", 0) or 0)
    upfront = float(cea.get("intervention_cost", 0) or 0) + \
        float(cea.get("test_cost", 0) or 0)
    horizon = float(cea.get("horizon_years") or
                    ep.value("horizon_years_personal"))
    payback = payback_period(upfront_cost=upfront,
                             annual_saving=cost_averted / max(1.0, horizon))

    # What would change the answer, in English.
    tornado = pooled.get("tornado") or []
    assumption_swing = sum(r.get("swing", 0) for r in tornado
                           if r.get("tier") == "assumption")
    total_swing = sum(r.get("swing", 0) for r in tornado) or 1
    share = 100.0 * assumption_swing / total_swing
    if share >= 50:
        weakest = ("Most of what drives this answer is guesswork, not "
                   "measurement. The biggest single driver is how much acting "
                   "on a genetic risk actually reduces it — a number nobody "
                   "has measured directly for most of these findings. Treat "
                   "the figures as a considered estimate, not a forecast.")
    elif share >= 25:
        weakest = ("A meaningful share of this answer rests on judgement "
                   "rather than published evidence. The direction is more "
                   "trustworthy than the exact numbers.")
    else:
        weakest = ("Most of what drives this answer comes from published "
                   "evidence rather than assumption.")

    return {
        "available": True,
        "verdict": verdict,
        "healthy_time": time_gain,
        "money": {
            "cost_averted": plain_money(cost_averted),
            "upfront": plain_money(upfront),
            "net": plain_money(-float(cea.get("incremental_cost", 0) or 0)),
            "plain": (
                f"Acting on these findings costs about {plain_money(upfront)} "
                f"and is modelled to avoid about {plain_money(cost_averted)} "
                f"of medical spending over {horizon:.0f} years."),
        },
        "payback": payback,
        "number_needed_to_screen": nns,
        "confidence": conf,
        "actions": plain_actions(conds),
        "what_would_change_it": weakest,
        "assumption_share": round(share, 1),
        "bottom_line": (
            f"{verdict['headline']}. {time_gain['plain']} "
            + (f"In about {conf['n_in_100']} of every 100 runs of the model, "
               f"acting was worth it." if conf.get("available") else "")),
        "disclaimer": (
            "These are model estimates for someone with your genetic results, "
            "using population-average figures. They are not medical advice, "
            "not a prediction about you specifically, and not a substitute "
            "for talking to a clinician."),
    }
