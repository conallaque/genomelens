"""
Cohort Simulator — from "what is MY genome worth" to "what is a MARKET worth"
=============================================================================

The value-of-information engine answers the question for one person. That is a
personal result, not a market analysis. This module simulates a synthetic
population of customers, runs each one through the same validated economics, and
reports the **distribution** of value — which is what a payer, an employer, or a
marketing team actually needs.

Four analyses live here:

  1. **Cohort simulation** — sample N synthetic people from published population
     prevalences, score each one, and return the distribution of value.
  2. **Segment analysis** — slice that distribution by age, family history and
     ancestry, so you can see who benefits most.
  3. **Demand curve** — at each price, what share of the population is "in the
     money"? Produces the chip-vs-whole-genome crossover point.
  4. **Adoption (Bass diffusion)** — how uptake spreads over time, combined with
     the behavioral adoption gap, to separate "no value" from "value not acted on".

**Where the heterogeneity comes from.** Each person is scored with point-estimate
parameters rather than a full probabilistic sensitivity analysis. That is
deliberate: here we are modeling *between-person* variation (do you carry an
actionable variant? how old are you?), not *parameter* uncertainty, which the PSA
in ``value_of_information`` already handles. Mixing the two would double-count.

Every prevalence below is sourced. The economic parameters are inherited from
``value_of_information`` and remain illustrative — see that module's caveats.
"""

from __future__ import annotations

from collections.abc import Sequence

try:
    import numpy as np
    _HAVE_NP = True
except Exception:                       # pragma: no cover
    np = None
    _HAVE_NP = False

from econ.value_of_information import (
    DISCOUNT_RATE,
    WTP,
    _finding_nmb,
)

try:
    from risk.genomic_statistics import PORTABILITY
except Exception:                       # pragma: no cover
    PORTABILITY = {"european": 1.0, "african": 0.25, "east_asian": 0.5,
                   "south_asian": 0.65, "hispanic_latino": 0.55}


# ── Population prevalences — every figure cited ───────────────────────────────
#
# These describe how common each *finding type* is in an unselected population.
# They are the engine of the whole simulation: they decide how many synthetic
# customers get an actionable result at all.

PREVALENCE = {
    "pgx_actionable": {
        "p": 0.95,
        "src": "Van Driest et al. (2014) Clin Pharmacol Ther; Dunnenberger et al. "
               "(2015) Annu Rev Pharmacol Toxicol — >90% of individuals carry at "
               "least one actionable pharmacogenomic variant.",
        "plain": "Almost everyone carries at least one gene variant that changes "
                 "how they respond to some common medication.",
    },
    "acmg_secondary": {
        "p": 0.030,
        "src": "Miller et al. (2023) Genet Med (ACMG SF v3.2); population "
               "sequencing programs (e.g. Geisinger MyCode) report ~2-3.5% of "
               "unselected adults carrying a reportable pathogenic variant.",
        "plain": "About 3 in 100 people carry a serious but medically actionable "
                 "variant — the kind worth screening or preventing against.",
    },
    "recessive_carrier": {
        "p": 0.60,
        "src": "Lazarin et al. (2013) Genet Med — carrier frequency rises with "
               "panel size; expanded panels identify a carrier state in the "
               "majority of individuals.",
        "plain": "Most people silently carry a recessive condition. It rarely "
                 "affects them, but it matters when planning a family.",
    },
    "apoe_e4": {
        "p": 0.25,
        "src": "Farrer et al. (1997) JAMA; consistently replicated — roughly a "
               "quarter of people of European ancestry carry at least one APOE e4 "
               "allele.",
        "plain": "About 1 in 4 people carry a version of the APOE gene linked to "
                 "higher dementia risk.",
    },
    "prs_high_decile": {
        "p": 0.10,
        "src": "True by construction — the top decile of any polygenic risk score "
               "distribution contains 10% of the population.",
        "plain": "By definition, 1 in 10 people sit in the top 10% of genetic risk "
                 "for any given common disease.",
    },
    "family_history": {
        "p": 0.30,
        "src": "Scheuner et al. (1997) Am J Med Genet; CDC family-history "
               "surveillance — roughly a third of adults report a significant "
               "family history of a common chronic disease.",
        "plain": "Around a third of people have a family history strong enough to "
                 "change how their own risk should be read.",
    },
}

# Age distribution of consumer-genomics purchasers. Deliberately skewed to
# 30-60: buyers cluster in mid-adulthood rather than matching the census.
AGE_DIST = {
    "mean": 45.0, "sd": 13.0, "min": 18.0, "max": 80.0,
    "src": "Approximate purchaser skew for direct-to-consumer genomics; treated "
           "as an assumption, not an observed distribution.",
}

# Ancestry mix for a US consumer base. Matters because polygenic scores were
# derived mostly in European-ancestry cohorts and transfer poorly (see
# genomic_statistics.PORTABILITY).
ANCESTRY_MIX = {
    "european": 0.62, "hispanic_latino": 0.15, "african": 0.12,
    "east_asian": 0.06, "south_asian": 0.05,
}
ANCESTRY_SRC = ("Approximate US population mix (US Census 2020 broad categories); "
                "consumer-genomics customer bases skew more European than the "
                "general population, so this is a conservative assumption.")


# ── persona sampling ──────────────────────────────────────────────────────────

def _sample_persona(rng) -> dict:
    """Draw one synthetic customer from the published prevalences above."""
    age = float(np.clip(rng.normal(AGE_DIST["mean"], AGE_DIST["sd"]),
                        AGE_DIST["min"], AGE_DIST["max"]))
    anc_keys = list(ANCESTRY_MIX)
    ancestry = str(rng.choice(anc_keys, p=[ANCESTRY_MIX[k] for k in anc_keys]))
    return {
        "age": round(age, 1),
        "ancestry": ancestry,
        "family_history": bool(rng.random() < PREVALENCE["family_history"]["p"]),
        "pgx_actionable": bool(rng.random() < PREVALENCE["pgx_actionable"]["p"]),
        "acmg_secondary": bool(rng.random() < PREVALENCE["acmg_secondary"]["p"]),
        "recessive_carrier": bool(rng.random() < PREVALENCE["recessive_carrier"]["p"]),
        "apoe_e4": bool(rng.random() < PREVALENCE["apoe_e4"]["p"]),
        "prs_high": bool(rng.random() < PREVALENCE["prs_high_decile"]["p"]),
    }


def _persona_findings(p: dict) -> list[dict]:
    """Convert a persona into the finding dicts the VOI engine already scores.

    ``wgs_only`` marks findings a genotyping chip cannot deliver — rare pathogenic
    variants and comprehensive carrier screening need sequencing. Chips *can* do
    pharmacogenomics, polygenic scores and APOE, so those are marked False. That
    split is what makes the chip-to-whole-genome comparison meaningful.
    """
    out: list[dict] = []
    horizon = max(5, int(85 - p["age"]))          # remaining years to benefit
    # Polygenic scores transfer poorly outside European-ancestry cohorts, so the
    # value they carry is attenuated rather than assumed constant.
    port = PORTABILITY.get(p["ancestry"], 0.5)

    # Model the health gain as an ANNUAL flow of quality-adjusted life-years over
    # the remaining horizon, not a fixed lump. Preventing a disease at 30 protects
    # more remaining life-years than preventing it at 70 — this is the Grossman
    # health-capital effect, and it is what makes value decline with age. (A fixed
    # QALY lump would do the reverse, because _finding_nmb applies the average
    # per-year discount, which shrinks as the horizon lengthens.)
    def q(annual_qaly: float) -> float:
        return annual_qaly * horizon

    if p["pgx_actionable"]:
        out.append({"label": "Actionable pharmacogenomic variant", "kind": "pgx",
                    "pgx_key": "PGx-generic", "intervention": 100.0,
                    "wgs_only": False, "haircut": 1.0, "confidence": "high"})
    if p["acmg_secondary"]:
        # Family history raises the posterior probability the variant matters.
        pen = 0.45 if p["family_history"] else 0.32
        out.append({"label": "Pathogenic variant (ACMG secondary finding)",
                    "kind": "coi", "coi_key": "BreastOvarian", "p_event": pen,
                    "rrr": 0.45, "qaly": q(0.040), "intervention": 1_500.0,
                    "horizon": horizon, "wgs_only": True, "haircut": 1.0,
                    "confidence": "high"})
    if p["recessive_carrier"]:
        # Reproductive-planning value is a one-time decision, not life-expectancy
        # scaled, so it keeps a fixed small QALY and short horizon.
        out.append({"label": "Recessive carrier status", "kind": "coi",
                    "coi_key": "Pathogenic", "p_event": 0.04, "rrr": 0.50,
                    "qaly": 0.2, "intervention": 400.0, "horizon": 5,
                    "wgs_only": True, "haircut": 1.0, "confidence": "moderate"})
    if p["prs_high"]:
        out.append({"label": "High polygenic risk (cardiometabolic)", "kind": "coi",
                    "coi_key": "CAD", "p_event": 0.20, "rrr": 0.30, "qaly": q(0.025),
                    "intervention": 500.0, "horizon": horizon, "wgs_only": False,
                    "haircut": port, "confidence": "moderate", "prevention": True})
    if p["apoe_e4"]:
        out.append({"label": "APOE e4 carrier", "kind": "coi",
                    "coi_key": "Alzheimer", "p_event": 0.15, "rrr": 0.20,
                    "qaly": q(0.035), "intervention": 500.0, "horizon": horizon,
                    "wgs_only": False, "haircut": 1.0, "confidence": "moderate",
                    "prevention": True})
    # Mark the life-expectancy-scaled disease-prevention finding (ACMG) too, so the
    # age gradient can be measured on prevention value separately from the flat
    # pharmacogenomic base that dominates most customers' totals.
    for f in out:
        if f.get("coi_key") == "BreastOvarian":
            f["prevention"] = True
    return out


def _score(findings: Sequence[dict], wtp: float) -> tuple[float, float, float]:
    """(total NMB, sequencing-only NMB, prevention quality-of-life value).

    The third term isolates the *quality-of-life* value of the life-expectancy-
    scaled prevention findings — QALYs gained x willingness-to-pay, before the
    averted-cost offset. This is where the Grossman effect is clean: preventing a
    disease earlier protects more remaining healthy years, so this component
    declines with age. (Net prevention value including averted medical cost is
    roughly flat, because a nearer disease is discounted less — the two forces
    partly cancel. Total value barely moves with age at all, because it is
    dominated by age-independent pharmacogenomics.)
    """
    total = wgs_only = prevention_qol = 0.0
    for f in findings:
        nmb, _dcost, dqaly, _interv = _finding_nmb(f, wtp, DISCOUNT_RATE)
        total += nmb
        if f.get("wgs_only"):
            wgs_only += nmb
        if f.get("prevention"):
            prevention_qol += dqaly * wtp
    return total, wgs_only, prevention_qol


# ── 1. cohort simulation ──────────────────────────────────────────────────────

def simulate_cohort(n: int = 10_000, seed: int = 20260803,
                    wtp: float | None = None) -> dict:
    """**Simulate a synthetic customer population and score every member.**

    *In plain English:* instead of asking "what is this one person's genome worth",
    we invent ten thousand plausible customers — with realistic ages, ancestries,
    family histories and genetic findings drawn from published rates — and run each
    one through the same economics. The result is a picture of the whole market
    rather than a single anecdote.

    The headline finding is almost always **concentration**: most people carry only
    a common pharmacogenomic result worth a modest amount, while a small minority
    with a serious actionable variant carry most of the total value. Averages hide
    that, which is why the percentiles matter more than the mean.
    """
    if not _HAVE_NP:
        return {"available": False, "reason": "numpy required for cohort simulation"}
    rng = np.random.default_rng(seed)
    lam = float(wtp if wtp is not None else WTP["base"])
    if n < 1:
        return {"available": False, "reason": "cohort size must be >= 1"}

    personas = []
    totals, wgs_parts, prev_parts = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        p = _sample_persona(rng)
        t, w, pv = _score(_persona_findings(p), lam)
        p["value"], p["value_wgs_only"], p["value_prevention_qol"] = \
            round(t), round(w), round(pv)
        personas.append(p)
        totals[i], wgs_parts[i], prev_parts[i] = t, w, pv

    order = np.argsort(totals)[::-1]
    top10 = totals[order[:max(1, n // 10)]].sum()
    concentration = float(top10 / totals.sum()) if totals.sum() > 0 else 0.0

    pct = {f"p{q}": round(float(np.percentile(totals, q)))
           for q in (5, 10, 25, 50, 75, 90, 95, 99)}
    return {
        "available": True,
        "n": n, "seed": seed, "wtp": lam,
        "mean_value": round(float(totals.mean())),
        "median_value": round(float(np.median(totals))),
        "sd_value": round(float(totals.std())),
        "percentiles": pct,
        "mean_wgs_only_value": round(float(wgs_parts.mean())),
        "share_of_value_in_top_decile": round(concentration, 3),
        "share_with_any_value": round(float((totals > 0).mean()), 3),
        "personas_sample": personas[:25],
        "_totals": totals, "_wgs": wgs_parts, "_personas": personas,
        "plain_english": (
            f"Across {n:,} simulated customers the average modeled value is "
            f"${totals.mean():,.0f}, but the median is ${np.median(totals):,.0f} — "
            f"and the top 10% of customers account for "
            f"{concentration:.0%} of all the value. In other words, a small number "
            f"of people get a great deal out of testing and most get a modest "
            f"amount. Marketing to the average customer misses this."),
        "sources": {k: v["src"] for k, v in PREVALENCE.items()},
        "caveat": ("Synthetic population built from published prevalences; the "
                   "economic parameters are illustrative. This estimates the shape "
                   "of a market, not a revenue forecast."),
    }


# ── 2. segment analysis ───────────────────────────────────────────────────────

def segment_analysis(cohort: dict) -> dict:
    """**Who gets the most value?** Slice the cohort by age, family history, ancestry.

    *In plain English:* the same test is worth more to some people than others.
    Younger customers have more years over which prevention pays off, people with a
    family history have a higher chance their variant matters, and polygenic scores
    are less reliable outside the populations they were developed in. This turns the
    economics into a targeting recommendation.
    """
    if not cohort.get("available"):
        return {"available": False}
    ps = cohort["_personas"]

    def _agg(rows):
        vals = [r["value"] for r in rows]
        prev = [r["value_prevention_qol"] for r in rows]
        return {"n": len(rows),
                "mean_value": round(float(np.mean(vals))) if rows else 0,
                "median_value": round(float(np.median(vals))) if rows else 0,
                "mean_prevention_qol_value": round(float(np.mean(prev))) if rows else 0}

    bands = [(18, 30), (30, 40), (40, 50), (50, 60), (60, 81)]
    by_age = [{"band": f"{a}-{b-1}", **_agg([r for r in ps if a <= r["age"] < b])}
              for a, b in bands]
    by_fh = {("with family history" if k else "no family history"):
             _agg([r for r in ps if r["family_history"] is k]) for k in (True, False)}
    by_anc = {a: _agg([r for r in ps if r["ancestry"] == a]) for a in ANCESTRY_MIX}

    # The age gradient is measured on PREVENTION value, where the Grossman effect
    # lives. Total value barely moves with age because age-independent
    # pharmacogenomics dominates every customer's total.
    young = next(x for x in by_age if x["band"] == "30-39")["mean_prevention_qol_value"]
    old = next(x for x in by_age if x["band"] == "60-80")["mean_prevention_qol_value"]
    gradient = round(young / old, 2) if old else None
    # Guard the narrative against an empty band (small cohort): only quote the
    # multiple when both bands are populated, else describe the gradient in words.
    gradient_phrase = (f"a customer in their thirties gets about {gradient:.1f}x the "
                       f"prevention benefit of one in their sixties"
                       if gradient else
                       "prevention benefit is markedly higher for younger customers")
    return {
        "available": True,
        "by_age_band": by_age,
        "by_family_history": by_fh,
        "by_ancestry": by_anc,
        "age_gradient_ratio": gradient,
        "plain_english": (
            f"The two products age very differently. Pharmacogenomic value — how you "
            f"respond to medications — is roughly constant across age, so that "
            f"product sells to everyone. Disease-prevention value, by contrast, "
            f"falls steadily with age: {gradient_phrase}, "
            f"because acting early protects more remaining years of health. So the "
            f"marketing split writes itself — pharmacogenomics to all ages, "
            f"prevention and whole-genome upgrades to younger customers. Family "
            f"history raises value further, and polygenic-score value is lower for "
            f"non-European ancestries — not because those customers matter less, but "
            f"because the underlying science was largely developed in European "
            f"cohorts and transfers imperfectly."),
        "sources": {
            "age gradient": "Grossman (1972) J Polit Econ — health capital "
                            "depreciates with age, so investment pays off over "
                            "fewer remaining years.",
            "family history": "Scheuner et al. (1997); family history raises the "
                              "posterior probability a variant is consequential.",
            "ancestry": "Martin et al. (2019) Nat Genet; Privé et al. (2022) AJHG — "
                        "polygenic scores lose accuracy with genetic distance from "
                        "the discovery cohort.",
        },
        "caveat": ("The ancestry result is a limitation of the science, not a "
                   "property of the customer, and should never be used to "
                   "deprioritise a group. It is an argument for more diverse "
                   "genomic research."),
    }


# ── 3. demand curve / price sensitivity ───────────────────────────────────────

def demand_curve(cohort: dict,
                 prices: Sequence[float] = (0, 99, 199, 299, 399, 599, 799, 999,
                                            1499, 1999, 2999, 4999)) -> dict:
    """**At each price, what share of customers get more value than they pay?**

    *In plain English:* a customer should rationally buy when the expected value of
    the test exceeds its price. Sweep the price and you get a demand curve — the
    share "in the money" at each point — plus modeled revenue per thousand
    customers, which peaks somewhere in the middle. This is a willingness-to-pay
    estimate grounded in modeled health value, not a survey.
    """
    if not cohort.get("available"):
        return {"available": False}
    totals = cohort["_totals"]
    wgs = cohort["_wgs"]

    rows = []
    for pr in prices:
        share = float((totals > pr).mean())
        rows.append({
            "price": pr,
            "share_in_the_money": round(share, 3),
            "revenue_per_1000": round(share * 1000 * pr),
        })
    best = max(rows, key=lambda r: r["revenue_per_1000"])

    # Chip-vs-sequencing crossover: the extra you can charge for sequencing is
    # bounded by the extra value only sequencing delivers.
    mean_wgs_premium = float(wgs.mean())
    median_wgs_premium = float(np.median(wgs))
    return {
        "available": True,
        "curve": rows,
        "revenue_maximising_price": best["price"],
        "revenue_at_that_price_per_1000": best["revenue_per_1000"],
        "mean_sequencing_premium": round(mean_wgs_premium),
        "median_sequencing_premium": round(median_wgs_premium),
        "plain_english": (
            f"Of the prices tested, ${best['price']:,} generates the most modeled "
            f"revenue per thousand customers (${best['revenue_per_1000']:,}), because "
            f"raising the price further loses more customers than it gains in margin. "
            f"Separately, the findings only sequencing can deliver are worth about "
            f"${mean_wgs_premium:,.0f} on average — that is the economic ceiling on "
            f"what a whole-genome upgrade can justifiably cost over a chip."),
        "sources": {
            "in-the-money rule": "Standard consumer-surplus logic: a rational buyer "
                                 "purchases when expected value exceeds price.",
            "caveat on behavior": "Real buyers systematically under-buy relative to "
                                   "this rule — see adoption_curve() and the "
                                   "prospect-theory adjustment in value_of_information.",
        },
        "caveat": ("This is a *normative* demand curve — what rational buyers would "
                   "do. Observed conversion is always lower. Use it as a ceiling, "
                   "not a forecast."),
    }


# ── 5. adoption / diffusion ───────────────────────────────────────────────────

def adoption_curve(periods: int = 10, p_innovate: float = 0.03,
                   q_imitate: float = 0.38, market_size: float = 1_000_000,
                   behavioral_drag: float = 0.30) -> dict:
    """**Bass diffusion — how adoption spreads, and why it lags the value.**

    *In plain English:* new products spread in a predictable S-curve. A few
    "innovators" buy immediately; everyone else buys because they see others doing
    it. The Bass model captures both forces with two numbers: ``p`` (innovation) and
    ``q`` (imitation, usually much larger — word of mouth dominates).

    On top of that we apply a **behavioral drag**. The economics may say a test is
    clearly worth buying, yet people still don't — because a certain cost today
    feels heavier than an uncertain benefit years away (loss aversion and present
    bias, quantified in ``value_of_information.analyze_behavioral``). The gap
    between the two curves is not a value problem; it is a *marketing* problem.
    """
    cum, rows = 0.0, []
    for t in range(1, periods + 1):
        frac = cum / market_size if market_size else 0.0
        new = (p_innovate + q_imitate * frac) * (market_size - cum)
        cum += new
        rows.append({
            "period": t,
            "new_adopters": round(new),
            "cumulative_adopters": round(cum),
            "penetration": round(cum / market_size, 4),
            "cumulative_with_behavioral_drag": round(cum * (1 - behavioral_drag)),
        })
    peak = max(rows, key=lambda r: r["new_adopters"])
    return {
        "available": True,
        "p_innovation": p_innovate, "q_imitation": q_imitate,
        "market_size": market_size, "behavioral_drag": behavioral_drag,
        "curve": rows,
        "peak_period": peak["period"],
        "final_penetration": rows[-1]["penetration"],
        "adopters_lost_to_behavioral_drag": round(rows[-1]["cumulative_adopters"]
                                                   * behavioral_drag),
        "plain_english": (
            f"Adoption peaks in period {peak['period']} and reaches "
            f"{rows[-1]['penetration']:.0%} of the market by period {periods}. "
            f"Word of mouth (q={q_imitate}) matters far more than advertising to "
            f"innovators (p={p_innovate}) — the single biggest lever is getting "
            f"early customers to talk. Applying the behavioral drag, roughly "
            f"{rows[-1]['cumulative_adopters'] * behavioral_drag:,.0f} people who "
            f"would rationally benefit still never buy."),
        "sources": {
            "diffusion model": "Bass, F.M. (1969) 'A New Product Growth for Model "
                               "Consumer Durables,' Management Science 15(5):215-227. "
                               "Typical values p~0.03, q~0.38 across consumer durables.",
            "behavioral drag": "Kahneman & Tversky (1979) prospect theory (loss "
                                "aversion ~2.25x); Laibson (1997) quasi-hyperbolic "
                                "discounting.",
        },
        "caveat": ("Bass parameters are category averages, not fitted to genomics. "
                   "Treat the shape as informative and the level as illustrative."),
    }


def data_asset_ltv(initial_value: float = 4_463.0, years: int = 10,
                   knowledge_growth: float = 0.07, discount_rate: float = DISCOUNT_RATE,
                   reanalysis_cost: float = 0.0) -> dict:
    """**#4 — The genome is bought once but appreciates forever.**

    *In plain English:* unlike almost any other purchase, a sequenced genome does
    not depreciate — it *appreciates*. You buy the data once, but the science that
    interprets it keeps growing: ClinVar adds variant classifications every month,
    new pharmacogenomic guidelines are published, and new polygenic scores appear.
    Because the file can be re-analyzed locally for free, all of that future
    knowledge accrues to the person who already owns their genome. This models the
    genome as an appreciating data asset.

    The value in year *t* grows with the stock of genomic knowledge and is then
    discounted back to today:

        V(t) = V0 * (1 + g)^t ,   PV = Σ_t V(t) * (1+r)^(-t)

    The headline output is the **appreciation share** — how much of the genome's
    lifetime value arrives *after* the initial purchase. That is the economic
    foundation of a subscription or re-contact model: you are not re-selling the
    test, you are delivering the value that accrued since.
    """
    rows, pv_total, pv_flat = [], 0.0, 0.0
    for t in range(0, years + 1):
        annual_value = initial_value * (1.0 + knowledge_growth) ** t
        disc = annual_value / ((1.0 + discount_rate) ** t)
        # Counterfactual: the same genome delivering FLAT value (no knowledge
        # growth) each year. The gap between the two isolates true appreciation
        # from mere durability.
        disc_flat = initial_value / ((1.0 + discount_rate) ** t)
        net = disc - reanalysis_cost
        pv_total += net
        pv_flat += disc_flat - reanalysis_cost
        rows.append({
            "year": t,
            "annual_value_nominal": round(annual_value),
            "annual_value_discounted": round(disc),
            "cumulative_pv": round(pv_total),
        })
    pv_initial = rows[0]["annual_value_discounted"]
    # Two DISTINCT effects, reported separately so neither is oversold:
    #  - durability: value arriving after year 0 even with NO growth (the asset lasts).
    #  - appreciation: the EXTRA lifetime value that knowledge growth adds on top.
    durability_share = 1.0 - (pv_initial / pv_flat) if pv_flat > 0 else 0.0
    appreciation_premium = (pv_total / pv_flat - 1.0) if pv_flat > 0 else 0.0
    return {
        "available": True,
        "initial_value": round(initial_value),
        "knowledge_growth_rate": knowledge_growth,
        "years": years,
        "lifetime_pv": round(pv_total),
        "lifetime_pv_no_growth": round(pv_flat),
        "value_at_purchase": round(pv_initial),
        "value_accrued_after_purchase": round(pv_total - pv_initial),
        "durability_share": round(durability_share, 3),
        "appreciation_premium": round(appreciation_premium, 3),
        "curve": rows,
        "plain_english": (
            f"Two separate things make a genome a good asset, and it is worth not "
            f"conflating them. First, durability: even if the science never improved, "
            f"a ${initial_value:,.0f} genome delivers value every year, so about "
            f"{durability_share:.0%} of its lifetime value arrives after the sale "
            f"simply because the data keeps working. Second, appreciation: because "
            f"the file is re-analyzed for free as knowledge grows, its lifetime value "
            f"is about {appreciation_premium:.0%} higher than a genome frozen at "
            f"today's knowledge. Together they are the economic case for staying in "
            f"contact with a customer rather than treating the test as a one-time "
            f"sale — but the honest headline is durability, with appreciation as the "
            f"upside."),
        "sources": {
            "ClinVar growth": "Landrum et al. (2018) Nucleic Acids Res — ClinVar "
                              "variant records grow rapidly year on year.",
            "PGx guideline expansion": "CPIC (cpicpgx.org) — the set of "
                                       "gene-drug guidelines expands continuously.",
            "PRS catalog growth": "Lambert et al. (2021) Nat Genet — the PGS "
                                    "Catalog has grown from a handful to thousands "
                                    "of published scores.",
            "free re-analysis / early ownership": "Mirrors the real-options result "
                                                  "in value_of_information: because "
                                                  "re-analysis is costless, knowledge "
                                                  "growth accrues to early owners.",
        },
        "caveat": ("The knowledge-growth rate is an assumption; it has been high "
                   "historically but is not guaranteed to continue. Illustrative."),
    }


_REFERENCE_COHORT = None


def _reference_totals():
    """A cached reference population value distribution, for percentile placement."""
    global _REFERENCE_COHORT
    if _REFERENCE_COHORT is None:
        _REFERENCE_COHORT = simulate_cohort(n=8000, seed=20260803)["_totals"]
    return _REFERENCE_COHORT


def personalize_for_report(voi_result: dict, age: float = 40.0,
                           wtp: float | None = None) -> dict:
    """Turn ONE person's value-of-information result into the individually-relevant
    health-economics panels for their personal report.

    Produces four things that depend on *this* genome (not the market):
      * a personal cost-effectiveness frontier over their real options (no test /
        chip / whole genome), so "is sequencing worth it *for me*?" is answered with
        dominance logic rather than a bare marginal number;
      * a personal CEAC — how confident that recommendation is given their own
        uncertainty interval;
      * their genome as an appreciating data asset (seeded with their own value);
      * where their value falls in the population distribution (a percentile).

    The market-level analyses (demand, adoption, distributional CEA, validation) are
    deliberately NOT included here — they describe a market, not a person.
    """
    if not _HAVE_NP or not voi_result or not voi_result.get("available"):
        return {"available": False}
    try:
        from econ import frontier as hf
    except Exception:
        return {"available": False}

    lam = float(wtp if wtp is not None else WTP["base"])
    if lam <= 0:                       # willingness-to-pay must be positive (divisor)
        return {"available": False, "reason": "willingness-to-pay must be > 0"}
    total = float(voi_result.get("voi_expost_mean",
                                 voi_result.get("voi_expost_point", 0.0)))
    wgs_marginal = float(voi_result.get("marginal_chip_to_wgs", 0.0) or 0.0)
    input_type = voi_result.get("input_type", "chip")
    # Split the person's value into chip-deliverable and sequencing-only parts.
    if input_type == "wgs":
        chip_value = max(0.0, total - wgs_marginal)
        wgs_value = total
    else:
        chip_value = total
        wgs_value = total + wgs_marginal
    # Uncertainty (SD) recovered from the reported credible interval.
    lo, hi = voi_result.get("voi_ci_low"), voi_result.get("voi_ci_high")
    sd = max(1.0, (float(hi) - float(lo)) / 3.92) if (lo is not None and hi is not None) \
        else 0.3 * max(1.0, total)

    base_qaly = 20.0
    # Dollars of health value → QALYs at the threshold.
    strategies = [
        {"name": "No testing", "cost": 0.0, "qaly": base_qaly, "qaly_sd": 0.001},
        {"name": "Genotyping chip", "cost": 100.0,
         "qaly": base_qaly + chip_value / lam, "qaly_sd": sd / lam},
        {"name": "Whole-genome sequencing", "cost": 600.0,
         "qaly": base_qaly + wgs_value / lam, "qaly_sd": sd / lam},
    ]
    frontier = hf.cost_effectiveness_frontier(strategies, wtp=lam)
    ceac = hf.frontier_psa(strategies, n_mc=3000)
    # The appreciation horizon depends on age: a younger person has more remaining
    # years over which free re-analysis keeps paying off. Capped at 20y because
    # knowledge-growth projections beyond that are too speculative to report.
    ltv_years = max(5, min(20, round(85.0 - float(age))))
    ltv = data_asset_ltv(initial_value=max(1.0, total), years=ltv_years)

    # Percentile placement against the reference population.
    ref = _reference_totals()
    pct = float((ref < total).mean()) * 100.0

    # Report the confidence in the RECOMMENDED strategy (the frontier's choice), not
    # the CEAC's own argmax — they can disagree (e.g. a genome with no sequencing-only
    # value makes chip and WGS a coin-flip), and attributing the wrong strategy's
    # probability would be misleading. Mirrors the renderer's consistency fix.
    _rec = frontier.get("recommended_strategy")
    _at100 = next((r for r in (ceac.get("ceac") or [])
                   if r.get("wtp") == 100_000), None)
    rec_prob = (_at100 or {}).get("p_optimal", {}).get(_rec,
                                                       ceac.get("prob_optimal_at_100k", 0.0))
    # Use the ACTUAL appreciation horizon (age-dependent, 5-20y), and compare against
    # the correct baseline. appreciation_premium is vs a genome frozen at today's
    # knowledge — so label it that way, not "than at purchase" (a different quantity).
    yrs = ltv.get("years", 10)
    return {
        "available": True,
        "expected_value": round(total),
        "chip_value": round(chip_value),
        "wgs_only_value": round(wgs_marginal),
        "frontier": frontier,
        "ceac": ceac,
        "ltv": ltv,
        "population_percentile": round(pct),
        "plain_english": (
            f"For this genome specifically: the efficient testing choice is "
            f"{_rec}, optimal in {rec_prob:.0%} of simulations at ${lam:,.0f} per "
            f"healthy year. The modeled value of ${total:,.0f} sits at about the "
            f"{pct:.0f}th percentile of the population, and — because the data is "
            f"re-analyzed for free as science advances — over the next {yrs} years it "
            f"is worth roughly {ltv.get('appreciation_premium', 0):.0%} more than a "
            f"genome frozen at today's knowledge."),
        "caveat": ("Personal decision-analytic estimate on illustrative parameters — "
                   "not medical or financial advice."),
    }


def explain_cohort(n: int = 10_000, seed: int = 20260803) -> str:
    """One plain-English narrative tying all four analyses together, suitable for
    a non-technical reader or a slide."""
    c = simulate_cohort(n=n, seed=seed)
    if not c.get("available"):
        return "Cohort simulation unavailable (numpy required)."
    s, d, a = segment_analysis(c), demand_curve(c), adoption_curve()
    return "\n\n".join([
        "MARKET VALUE OF GENOMIC TESTING — plain-language summary",
        "1. HOW MUCH IS IT WORTH?\n   " + c["plain_english"],
        "2. WHO BENEFITS MOST?\n   " + s["plain_english"],
        "3. WHAT SHOULD IT COST?\n   " + d["plain_english"],
        "4. HOW FAST DOES IT SPREAD?\n   " + a["plain_english"],
        "IMPORTANT: " + c["caveat"],
    ])
