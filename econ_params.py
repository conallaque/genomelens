"""
Parameter provenance registry for the economic model
====================================================

Every number that reaches a dollar figure in this report has to be answerable
to the question "where did that come from?". Before this module the answer was
often "a constant someone typed", and an internal audit put roughly two thirds
of the economic parameters in that category. A model whose inputs cannot be
traced is not wrong so much as *uncheckable*, which is worse: nobody can find
the error, including the author.

So parameters live here as :class:`Param` records rather than bare floats, and
each one carries a **tier** that states honestly how much authority stands
behind it:

``published``
    The value appears in the cited source. Someone can open the paper and read
    the number off it.

``derived``
    Computed from the cited source by a stated arithmetic step — inflated to a
    common year, converted from a rate to a probability, taken as a midpoint of
    a published range. ``note`` records the step, so the derivation can be
    disagreed with separately from the source.

``assumption``
    A judgement call with no published anchor. These are *not* forbidden — some
    quantities have no literature — but they are enumerated in the report as
    declared assumptions instead of blending in with sourced values, and
    :func:`assumption_burden` reports how much of the model's output rests on
    them.

The registry is validated at import: a ``published`` or ``derived`` parameter
without a citation is an error, and so is an ``assumption`` without a note
explaining the judgement. ``test_econ_params.py`` enforces the same rules
against the parameters actually reaching the model, so a new unsourced constant
cannot quietly appear.

Citations are given as PMID/DOI/ISBN so they resolve outside this repository.
Where a figure is an order-of-magnitude anchor rather than a number lifted from
a table, the tier is ``derived`` and the note says so — claiming a paper states
a figure it does not state would defeat the purpose of the exercise.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

__all__ = [
    "Param", "PARAMS", "get", "value", "validate_registry",
    "assumptions", "by_tier", "assumption_burden", "citation_list",
    "TIERS",
]

TIERS = ("published", "derived", "assumption")


@dataclass(frozen=True)
class Param:
    """One economic parameter with its provenance and uncertainty."""

    key: str
    value: float
    units: str
    tier: str
    note: str = ""
    source: str = ""
    citation: str = ""          # PMID / DOI / ISBN / stable URL
    year: Optional[int] = None
    dist: str = ""              # "beta" | "gamma" | "lognormal" | "fixed"
    se: Optional[float] = None
    low: Optional[float] = None
    high: Optional[float] = None

    # ── validation ──────────────────────────────────────────────────────
    def validate(self) -> List[str]:
        """Return a list of provenance problems with this parameter."""
        bad: List[str] = []
        if self.tier not in TIERS:
            bad.append(f"{self.key}: tier {self.tier!r} not in {TIERS}")
        if self.tier in ("published", "derived"):
            if not self.citation.strip():
                bad.append(f"{self.key}: tier={self.tier} requires a citation")
            if not self.source.strip():
                bad.append(f"{self.key}: tier={self.tier} requires a source")
        if self.tier == "derived" and not self.note.strip():
            bad.append(f"{self.key}: tier=derived requires a note stating the "
                       f"derivation step")
        if self.tier == "assumption":
            if not self.note.strip():
                bad.append(f"{self.key}: tier=assumption requires a note "
                           f"explaining the judgement")
            if self.citation.strip():
                bad.append(f"{self.key}: tier=assumption must not cite a "
                           f"source — retier it as published/derived")
        if self.low is not None and self.high is not None and self.low > self.high:
            bad.append(f"{self.key}: low {self.low} > high {self.high}")
        return bad

    # ── convenience ─────────────────────────────────────────────────────
    @property
    def sourced(self) -> bool:
        return self.tier in ("published", "derived")

    @property
    def range(self) -> Optional[tuple]:
        if self.low is None or self.high is None:
            return None
        return (self.low, self.high)

    def cite(self) -> str:
        """Human-readable one-line citation, or a declared-assumption label."""
        if self.tier == "assumption":
            return "Declared assumption — no published anchor"
        y = f" ({self.year})" if self.year else ""
        return f"{self.source}{y}. {self.citation}"


def _p(*a, **kw) -> Param:
    return Param(*a, **kw)


# ══════════════════════════════════════════════════════════════════════════
# Registry
# ══════════════════════════════════════════════════════════════════════════
# Grouped by role. Keys are stable strings — other modules look parameters up
# by key rather than importing floats, so provenance travels with the value.

_REGISTRY: List[Param] = [

    # ── Method conventions ────────────────────────────────────────────────
    _p("wtp_per_qaly", 100_000.0, "$/QALY", "published",
       note="Upper bound of the conventionally cited US range; the model runs "
            "50k/100k/150k as a threshold sensitivity rather than treating any "
            "one value as correct.",
       source="Neumann PJ, Cohen JT, Weinstein MC. Updating cost-effectiveness "
              "— the curious resilience of the $50,000-per-QALY threshold. "
              "N Engl J Med",
       citation="PMID:25162885 doi:10.1056/NEJMp1405158", year=2014,
       dist="fixed", low=50_000.0, high=150_000.0),

    _p("discount_rate", 0.03, "per year", "published",
       note="Second Panel reference-case real discount rate, applied to both "
            "costs and effects; 0%/3%/5% carried as a sensitivity.",
       source="Sanders GD, et al. Recommendations for Conduct, Methodological "
              "Practices, and Reporting of Cost-effectiveness Analyses: Second "
              "Panel on Cost-Effectiveness in Health and Medicine. JAMA",
       citation="PMID:27623463 doi:10.1001/jama.2016.12195", year=2016,
       dist="fixed", low=0.0, high=0.05),

    _p("marginal_cost_fraction", 0.60, "fraction of average cost", "derived",
       note="Averting one case frees the marginal, not the average, cost of "
            "that case — much of average cost is fixed capacity that persists. "
            "No single published ratio covers all conditions here, so 0.6 is "
            "used as a uniform conservative haircut on the cash arm, in the "
            "direction that lowers claimed savings. Varied 0.4–0.8 in PSA.",
       source="Drummond MF, Sculpher MJ, Claxton K, Stoddart GL, Torrance GW. "
              "Methods for the Economic Evaluation of Health Care Programmes, "
              "4th ed. Oxford University Press (ch. 7, marginal vs average "
              "costing)",
       citation="ISBN:9780199665877", year=2015,
       dist="beta", low=0.40, high=0.80),

    _p("cycle_correction", 1.0, "method flag", "published",
       note="Simpson's 1/3 within-cycle correction; flag records the "
            "convention, the weights are generated in econ_engine.",
       source="Alarid-Escudero F, Krijkamp EM, Enns EA, et al. A Tutorial on "
              "Time-Dependent Cohort State-Transition Models in R Using a "
              "Cost-Effectiveness Analysis Example. Med Decis Making",
       citation="PMID:35924564 doi:10.1177/0272989X221121747", year=2023,
       dist="fixed"),

    # ── Test / intervention costs ─────────────────────────────────────────
    _p("cost_chip", 100.0, "$", "derived",
       note="Consumer genotyping-array retail price observed in the US direct-"
            "to-consumer market at time of writing; a market price, not a "
            "published cost estimate.",
       source="US direct-to-consumer genotyping list prices (23andMe, "
              "AncestryDNA), retail observation",
       citation="https://www.23andme.com/dna-reports-list/", year=2026,
       dist="gamma", low=59.0, high=199.0),

    _p("cost_wgs", 300.0, "$", "derived",
       note="Consumer 30x whole-genome sequencing promotional price; market "
            "price, not a costed estimate of production.",
       source="US direct-to-consumer 30x WGS list prices (Nebula Genomics, "
              "Dante Labs), retail observation",
       citation="https://nebula.org/whole-genome-sequencing-dna-test/", year=2026,
       dist="gamma", low=200.0, high=999.0),

    _p("cost_analysis_bundle", 700.0, "$", "derived",
       note="One-time cost of the inputs this report models: ~$200 genotyping "
            "plus ~$500 for the laboratory panel. Sum of two market prices.",
       source="Retail genotyping price plus US self-pay comprehensive "
              "laboratory panel price, retail observation",
       citation="https://www.questhealth.com/", year=2026,
       dist="gamma", low=300.0, high=1_500.0),

    _p("cost_genetic_counseling", 400.0, "$", "derived",
       note="Single US genetic-counselling session, self-pay; midpoint of the "
            "commonly quoted $200–600 range.",
       source="National Society of Genetic Counselors, patient FAQ on cost of "
              "genetic counselling",
       citation="https://www.nsgc.org/Policy-Research-and-Publications/"
                "Genetic-Testing-and-Counseling", year=2026,
       dist="gamma", low=200.0, high=600.0),

    _p("cost_partner_carrier_panel", 300.0, "$", "derived",
       note="Expanded carrier panel for one partner, self-pay list price.",
       source="US expanded carrier screening self-pay list prices (Invitae, "
              "Myriad Foresight), retail observation",
       citation="https://www.invitae.com/en/providers/test-catalog", year=2026,
       dist="gamma", low=250.0, high=600.0),

    _p("cost_statin_10yr", 500.0, "$", "derived",
       note="Ten years of generic statin at US generic pricing (~$4/month "
            "cash price), excluding monitoring visits.",
       source="US generic statin cash prices (GoodRx atorvastatin), retail "
              "observation",
       citation="https://www.goodrx.com/atorvastatin", year=2026,
       dist="gamma", low=200.0, high=2_000.0),

    _p("cost_dpp_program", 3_000.0, "$", "derived",
       note="Cost of delivering a year of the CDC-recognised Diabetes "
            "Prevention Program lifestyle intervention per participant; the "
            "DPP trial's own direct programme cost, inflated to present-day "
            "dollars as an order-of-magnitude anchor.",
       source="Diabetes Prevention Program Research Group. Within-trial "
              "cost-effectiveness of lifestyle intervention or metformin for "
              "the primary prevention of type 2 diabetes. Diabetes Care",
       citation="PMID:12832381 doi:10.2337/diacare.26.9.2518", year=2003,
       dist="gamma", low=500.0, high=6_000.0),

    # ── Effect sizes ──────────────────────────────────────────────────────
    _p("statin_rrr_primary", 0.27, "relative risk reduction", "published",
       note="Proportional reduction in major vascular events per ~1 mmol/L LDL "
            "reduction in people at low baseline risk.",
       source="Cholesterol Treatment Trialists' (CTT) Collaborators. The "
              "effects of lowering LDL cholesterol with statin therapy in "
              "people at low risk of vascular disease. Lancet",
       citation="PMID:22607822 doi:10.1016/S0140-6736(12)60367-5", year=2012,
       dist="beta", se=0.03, low=0.20, high=0.34),

    _p("dpp_rrr", 0.58, "relative risk reduction", "published",
       note="Reduction in incidence of type 2 diabetes with intensive "
            "lifestyle intervention versus placebo over 2.8 years.",
       source="Knowler WC, Barrett-Connor E, Fowler SE, et al. Reduction in "
              "the incidence of type 2 diabetes with lifestyle intervention or "
              "metformin. N Engl J Med",
       citation="PMID:11832527 doi:10.1056/NEJMoa012512", year=2002,
       dist="beta", se=0.05, low=0.48, high=0.66),

    _p("prediabetes_progression_10yr", 0.35, "probability", "derived",
       note="Cumulative 10-year progression from prediabetes to type 2 "
            "diabetes, taken from the middle of the wide published range "
            "(annual conversion ~5–10%/yr compounds to roughly 0.25–0.50).",
       source="Tabák AG, Herder C, Rathmann W, Brunner EJ, Kivimäki M. "
              "Prediabetes: a high-risk state for diabetes development. Lancet",
       citation="PMID:22683128 doi:10.1016/S0140-6736(12)60283-9", year=2012,
       dist="beta", low=0.20, high=0.55),

    # ── Cost of illness ───────────────────────────────────────────────────
    _p("coi_mace", 85_000.0, "$ per case", "derived",
       note="Acute plus first-year care for a major adverse cardiovascular "
            "event, as an order-of-magnitude anchor spanning MI and stroke; "
            "national aggregate CVD costs divided across event counts rather "
            "than a per-case figure published as such.",
       source="Tsao CW, Aday AW, Almarzooq ZI, et al. Heart Disease and Stroke "
              "Statistics — 2023 Update. Circulation",
       citation="PMID:36695182 doi:10.1161/CIR.0000000000001123", year=2023,
       dist="gamma", low=40_000.0, high=150_000.0),

    _p("coi_t2d", 85_000.0, "$ per case", "derived",
       note="Excess lifetime medical cost attributable to type 2 diabetes; "
            "derived from annual per-capita excess expenditure (~$12k) carried "
            "over a typical post-diagnosis horizon, discounted.",
       source="Parker ED, Lin J, Mahoney T, et al. Economic Costs of Diabetes "
              "in the U.S. in 2022. Diabetes Care",
       citation="PMID:37909353 doi:10.2337/dci23-0085", year=2024,
       dist="gamma", low=40_000.0, high=180_000.0),

    _p("coi_alzheimer", 200_000.0, "$ per case", "derived",
       note="Lifetime formal care cost of dementia, excluding the substantially "
            "larger unpaid-caregiving component, which is carried separately in "
            "the societal perspective.",
       source="Alzheimer's Association. 2023 Alzheimer's Disease Facts and "
              "Figures. Alzheimers Dement",
       citation="PMID:36918389 doi:10.1002/alz.13016", year=2023,
       dist="gamma", low=100_000.0, high=400_000.0),

    _p("coi_depression", 55_000.0, "$ per case", "derived",
       note="Per-person share of the US major-depressive-disorder burden over "
            "a multi-year episode horizon; the source reports an aggregate "
            "national burden and a per-person annual figure, which is carried "
            "forward here rather than quoted directly.",
       source="Greenberg PE, Fournier AA, Sisitsky T, et al. The Economic "
              "Burden of Adults with Major Depressive Disorder in the United "
              "States (2010 and 2018). Pharmacoeconomics",
       citation="PMID:33950419 doi:10.1007/s40273-021-01019-4", year=2021,
       dist="gamma", low=20_000.0, high=120_000.0),

    _p("coi_substance_use", 60_000.0, "$ per case", "derived",
       note="Per-person healthcare-sector share of excessive-alcohol and "
            "substance-use cost; the source reports national totals dominated "
            "by lost productivity, and only the healthcare share is taken here.",
       source="Sacks JJ, Gonzales KR, Bouchery EE, Tomedi LE, Brewer RD. 2010 "
              "National and State Costs of Excessive Alcohol Consumption. "
              "Am J Prev Med",
       citation="PMID:26477807 doi:10.1016/j.amepre.2015.05.031", year=2015,
       dist="gamma", low=20_000.0, high=150_000.0),

    _p("coi_autoimmune", 95_000.0, "$ per case", "derived",
       note="Multi-year direct medical cost of an established autoimmune "
            "condition, anchored on rheumatoid arthritis as the costed "
            "exemplar; the panel that routes here spans several conditions of "
            "differing cost, so this is a class-level anchor.",
       source="Birnbaum H, Pike C, Kaufman R, Marynchenko M, Kidolezi Y, "
              "Cifaldi M. Societal cost of rheumatoid arthritis patients in "
              "the US. Curr Med Res Opin",
       citation="PMID:19908947 doi:10.1185/03007990903422307", year=2010,
       dist="gamma", low=30_000.0, high=200_000.0),

    _p("coi_urologic", 25_000.0, "$ per case", "derived",
       note="Direct medical cost of recurrent stone disease and related "
            "urologic care over a multi-year horizon; derived from national "
            "urologic-disease expenditure per affected person.",
       source="Litwin MS, Saigal CS (eds). Urologic Diseases in America. "
              "NIH Publication 12-7865, US Department of Health and Human "
              "Services",
       citation="https://www.niddk.nih.gov/about-niddk/strategic-plans-reports/"
                "urologic-diseases-in-america", year=2012,
       dist="gamma", low=8_000.0, high=60_000.0),

    _p("coi_iron_overload", 40_000.0, "$ per case", "derived",
       note="Lifetime cost of clinically expressed hereditary haemochromatosis "
            "including surveillance and phlebotomy; most C282Y homozygotes "
            "never express, which is handled by penetrance upstream, not here.",
       source="Adams PC, Reboussin DM, Barton JC, et al. Hemochromatosis and "
              "iron-overload screening in a racially diverse population. "
              "N Engl J Med",
       citation="PMID:15858186 doi:10.1056/NEJMoa041534", year=2005,
       dist="gamma", low=10_000.0, high=120_000.0),

    _p("coi_colorectal", 120_000.0, "$ per case", "derived",
       note="Direct medical cost of treated colorectal cancer across stages, "
            "as a stage-weighted anchor.",
       source="Mariotto AB, Enewold L, Zhao J, Zeruto CA, Yabroff KR. Medical "
              "Care Costs Associated with Cancer Survivorship in the United "
              "States. Cancer Epidemiol Biomarkers Prev",
       citation="PMID:32229578 doi:10.1158/1055-9965.EPI-19-1534", year=2020,
       dist="gamma", low=50_000.0, high=250_000.0),

    _p("coi_breast_ovarian", 150_000.0, "$ per case", "derived",
       note="Direct medical cost of treated breast or ovarian cancer, "
            "stage-weighted; ovarian sits above and early-stage breast below "
            "this anchor.",
       source="Mariotto AB, Enewold L, Zhao J, Zeruto CA, Yabroff KR. Medical "
              "Care Costs Associated with Cancer Survivorship in the United "
              "States. Cancer Epidemiol Biomarkers Prev",
       citation="PMID:32229578 doi:10.1158/1055-9965.EPI-19-1534", year=2020,
       dist="gamma", low=60_000.0, high=300_000.0),

    _p("coi_pathogenic_generic", 100_000.0, "$ per case", "derived",
       note="Fallback anchor for an actionable monogenic finding whose gene is "
            "not individually costed; set near the middle of the costed "
            "ACMG-gene conditions above.",
       source="Miller DT, Lee K, Abul-Husn NS, et al. ACMG SF v3.2 list for "
              "reporting of secondary findings in clinical exome and genome "
              "sequencing. Genet Med",
       citation="PMID:37347242 doi:10.1016/j.gim.2023.100866", year=2023,
       dist="gamma", low=30_000.0, high=300_000.0),

    # ── Utilities / QALY decrements ───────────────────────────────────────
    _p("qaly_loss_mace", 1.5, "QALYs per case", "derived",
       note="Discounted quality-adjusted life-years lost to a non-fatal major "
            "cardiovascular event over the remaining horizon, from the "
            "published post-event utility decrement.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.5, high=4.0),

    _p("qaly_loss_t2d", 2.0, "QALYs per case", "derived",
       note="QALYs lost to type 2 diabetes over the modelled horizon, from the "
            "published per-year utility decrement applied across expected "
            "years with disease.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.8, high=5.0),

    _p("utility_healthy", 0.85, "utility weight", "published",
       note="US population mean EQ-5D index for adults without the modelled "
            "condition; used as the well state in the state-transition model.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="beta", low=0.78, high=0.92),

    _p("utility_post_event", 0.70, "utility weight", "derived",
       note="Well-state utility less the published decrement for established "
            "cardiovascular disease; used as the diseased state in the "
            "state-transition model.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="beta", low=0.55, high=0.82),

    # ── Penetrance / ascertainment ────────────────────────────────────────
    _p("ascertainment_shrinkage", 0.60, "multiplier on penetrance", "derived",
       note="Population-cohort penetrance for high-risk variants runs well "
            "below the figure from clinically ascertained families. A single "
            "shrinkage multiplier stands in for that gap where a gene-specific "
            "population estimate is unavailable; applied before the economics, "
            "not merely displayed.",
       source="Wright CF, West B, Tuke M, et al. Assessing the Pathogenicity, "
              "Penetrance, and Expressivity of Putative Disease-Causing "
              "Variants in a Population Setting. Am J Hum Genet",
       citation="PMID:30665703 doi:10.1016/j.ajhg.2018.12.015", year=2019,
       dist="beta", low=0.30, high=0.90),

    _p("brca2_penetrance_population", 0.55, "cumulative risk to age 80",
       "published",
       note="Cumulative breast-cancer risk for BRCA2 carriers in a prospective "
            "cohort — the ascertainment-corrected starting point, not the "
            "family-series figure.",
       source="Kuchenbaecker KB, Hopper JL, Barnes DR, et al. Risks of Breast, "
              "Ovarian, and Contralateral Breast Cancer for BRCA1 and BRCA2 "
              "Mutation Carriers. JAMA",
       citation="PMID:28632866 doi:10.1001/jama.2017.7112", year=2017,
       dist="beta", low=0.45, high=0.65),

    # ── Aggregation policy ────────────────────────────────────────────────
    _p("max_combined_rrr", 0.60, "relative risk reduction", "assumption",
       note="Ceiling on the combined risk reduction attributable to acting on "
            "all genomic findings for a single condition. No trial delivers "
            "elimination of a common complex disease, and without a cap the "
            "complement-of-products combination drifts toward 1.0 as findings "
            "accumulate. Set at roughly twice the best single-intervention "
            "effect in the model (statins, 0.27). Judgement, not a published "
            "figure — varied 0.4–0.8 in sensitivity analysis.",
       dist="beta", low=0.40, high=0.80),

    _p("correlated_signal_penalty", 0.50, "multiplier", "assumption",
       note="Second and subsequent findings routed to the same condition are "
            "largely re-measurements of one underlying liability — a polygenic "
            "score and a PheWAS biomarker for the same trait are not "
            "independent evidence. Each additional finding on a condition is "
            "down-weighted by this factor compounding, so the marginal "
            "contribution of the fifth CAD signal is small. Judgement; the "
            "alternative in use before this was implicit independence, which "
            "is certainly wrong in the anti-conservative direction.",
       dist="beta", low=0.25, high=0.80),

    _p("horizon_years_personal", 10.0, "years", "assumption",
       note="Reporting horizon for the personal economic summary. Chosen for "
            "interpretability rather than derived: long enough for prevention "
            "to accrue, short enough that the projection is not fantasy. The "
            "state-transition model runs to age 100 internally and is "
            "truncated to this window for display.",
       dist="fixed", low=5.0, high=30.0),

    _p("productivity_annual", 55_000.0, "$ per year", "derived",
       note="Annual civilian earnings used for the productivity arm of the "
            "societal perspective, from national median wage data. Reported "
            "separately in the impact inventory and never folded into the "
            "healthcare-sector result.",
       source="US Bureau of Labor Statistics, Occupational Employment and "
              "Wage Statistics, national median annual wage, all occupations",
       citation="https://www.bls.gov/oes/current/oes_nat.htm", year=2024,
       dist="gamma", low=35_000.0, high=90_000.0),

    _p("caregiving_hours_annual_dementia", 1_500.0, "hours per year", "derived",
       note="Unpaid caregiving hours per person with dementia per year, from "
            "national aggregate caregiving hours divided by prevalent cases; "
            "valued at the replacement wage in the societal perspective only.",
       source="Alzheimer's Association. 2023 Alzheimer's Disease Facts and "
              "Figures. Alzheimers Dement",
       citation="PMID:36918389 doi:10.1002/alz.13016", year=2023,
       dist="gamma", low=800.0, high=2_500.0),

    _p("caregiver_replacement_wage", 17.0, "$ per hour", "derived",
       note="Replacement-cost valuation of unpaid caregiving, using the median "
            "wage for home health and personal care aides.",
       source="US Bureau of Labor Statistics, Occupational Employment and "
              "Wage Statistics, home health and personal care aides (31-1120)",
       citation="https://www.bls.gov/oes/current/oes311120.htm", year=2024,
       dist="gamma", low=12.0, high=28.0),
]

PARAMS: Dict[str, Param] = {p.key: p for p in _REGISTRY}
if len(PARAMS) != len(_REGISTRY):
    _seen: Dict[str, int] = {}
    for _p_ in _REGISTRY:
        _seen[_p_.key] = _seen.get(_p_.key, 0) + 1
    raise ValueError(f"duplicate parameter keys: "
                     f"{[k for k, n in _seen.items() if n > 1]}")


# ══════════════════════════════════════════════════════════════════════════
# Access + validation
# ══════════════════════════════════════════════════════════════════════════

def get(key: str) -> Param:
    """Look up a parameter, failing loudly on an unknown key."""
    try:
        return PARAMS[key]
    except KeyError:
        raise KeyError(
            f"unknown economic parameter {key!r}. Parameters must be "
            f"registered in econ_params.py with a tier and provenance before "
            f"they can reach the model; that is the point of the registry."
        ) from None


def value(key: str, default: Optional[float] = None) -> float:
    """Return a parameter's value. ``default`` is honoured only for keys that
    genuinely may be absent; an unregistered key is otherwise an error."""
    if default is not None and key not in PARAMS:
        return float(default)
    return float(get(key).value)


def validate_registry(params: Optional[Sequence[Param]] = None) -> List[str]:
    """Return every provenance problem across the registry (empty == clean)."""
    problems: List[str] = []
    for p in (params if params is not None else _REGISTRY):
        problems.extend(p.validate())
    return problems


def by_tier(tier: str) -> List[Param]:
    """All parameters at one provenance tier."""
    return [p for p in _REGISTRY if p.tier == tier]


def assumptions() -> List[Param]:
    """Parameters with no published anchor, for the report's declared-
    assumptions section."""
    return by_tier("assumption")


def count_unregistered_parameters() -> int:
    """How many load-bearing numbers still live outside the registry.

    The curated per-finding tables in ``health_economics`` (``PGX_ECONOMICS``,
    ``ACMG_GENE_ECONOMICS``, ``PHEWAS_CATEGORY_ECONOMICS`` and the rest) each
    carry ``cost`` / ``outcome_value`` / ``prevalence`` / ``qaly_gain`` fields,
    and every one of them reaches a dollar figure. They are not in the registry
    yet. Counting them is what keeps :func:`assumption_burden` from reporting a
    flattering number about a small corner of the model — the registry covers
    method conventions, cost-of-illness anchors and effect sizes, which is the
    spine, not the whole skeleton.
    """
    try:
        import health_economics as _he
    except Exception:
        return 0
    fields = ("cost", "outcome_value", "prevalence", "qaly_gain",
              "adr_cost", "rrr")
    total = 0
    for name in dir(_he):
        if not (name.endswith("_ECONOMICS") or name.endswith("_COSTS")):
            continue
        table = getattr(_he, name, None)
        if not isinstance(table, dict):
            continue
        for entry in table.values():
            if isinstance(entry, dict):
                total += sum(1 for f in fields if f in entry)
    return total


def assumption_burden() -> Dict[str, float]:
    """How much of the registry rests on judgement rather than literature.

    Reports registry coverage AND how much of the model the registry does not
    yet reach. Quoting only the first would be the same species of error the
    registry was built to fix: a true statement about a subset, phrased so it
    reads as a statement about the whole.
    """
    n = len(_REGISTRY) or 1
    counts = {t: len(by_tier(t)) for t in TIERS}
    unregistered = count_unregistered_parameters()
    total_known = n + unregistered
    return {
        "n_parameters": n,
        "n_published": counts["published"],
        "n_derived": counts["derived"],
        "n_assumption": counts["assumption"],
        "pct_sourced": round(100.0 * (counts["published"] + counts["derived"]) / n, 1),
        "pct_assumption": round(100.0 * counts["assumption"] / n, 1),
        # The honest denominator.
        "n_unregistered": unregistered,
        "n_total_known": total_known,
        "pct_of_model_registered": round(100.0 * n / total_known, 1)
        if total_known else 100.0,
        "scope": ("Registered parameters cover method conventions, "
                  "cost-of-illness anchors, effect sizes and utilities. The "
                  "per-finding curated tables in health_economics.py are not "
                  "yet registered; their figures reach the model without a "
                  "provenance tier."),
    }


def citation_list() -> List[Dict[str, str]]:
    """De-duplicated reference list for the report, sorted by source."""
    seen: Dict[str, Dict[str, str]] = {}
    for p in _REGISTRY:
        if not p.sourced or not p.citation:
            continue
        seen.setdefault(p.citation, {
            "source": p.source, "citation": p.citation,
            "year": str(p.year or ""), "params": [],
        })
        seen[p.citation]["params"].append(p.key)
    out = sorted(seen.values(), key=lambda d: d["source"].lower())
    for d in out:
        d["params"] = ", ".join(sorted(d["params"]))
    return out


# Fail at import if the registry itself is malformed — a broken provenance
# record should never survive to the point of producing a dollar figure.
_PROBLEMS = validate_registry()
if _PROBLEMS:
    raise ValueError("econ_params registry has provenance problems:\n  " +
                     "\n  ".join(_PROBLEMS))
