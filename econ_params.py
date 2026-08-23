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

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

__all__ = [
    "Param", "PARAMS", "get", "value", "validate_registry",
    "assumptions", "by_tier", "assumption_burden", "citation_list",
    "overridden", "sampleable", "draw", "sample_all",
    "CURATED_SOURCE_IDS", "resolve_curated_source", "audit_curated_tables",
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

    # ── Per-condition QALY decrements ─────────────────────────────────────
    # These were briefly all pointing at qaly_loss_mace, which claimed a
    # non-fatal cardiovascular event's decrement for dementia, depression and
    # kidney stones alike. That is wrong on its face, and it also corrupted the
    # sensitivity analysis: one placeholder doing the work of seven conditions
    # dominated the tornado, so the report named it as the model's key driver
    # when it was really just the most overloaded constant.
    _p("qaly_loss_dementia", 4.0, "QALYs per case", "derived",
       note="Dementia carries a large utility decrement over a long course. "
            "Derived from the published EQ-5D decrement for dementia applied "
            "across expected years with disease — substantially larger than a "
            "cardiovascular event's, which is why it needs its own parameter.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=1.5, high=8.0),

    _p("qaly_loss_depression", 1.2, "QALYs per case", "derived",
       note="Recurrent major depression, from the published per-year utility "
            "decrement across an expected number of episode-years within the "
            "horizon rather than a lifetime total.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.4, high=3.5),

    _p("qaly_loss_autoimmune", 1.8, "QALYs per case", "derived",
       note="Established autoimmune disease, anchored on the rheumatoid-"
            "arthritis decrement as the costed exemplar; the panel routing "
            "here spans conditions of differing severity.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.5, high=4.5),

    _p("qaly_loss_substance_use", 1.5, "QALYs per case", "derived",
       note="Substance-use disorder over a multi-year course, from the "
            "published decrement for the condition class.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.4, high=4.0),

    _p("qaly_loss_urologic", 0.4, "QALYs per case", "derived",
       note="Recurrent stone disease and related urologic care: episodic and "
            "largely reversible, so the decrement is far smaller than the "
            "chronic conditions above. Kept distinct precisely because "
            "borrowing a cardiovascular figure here would overstate it "
            "several-fold.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.1, high=1.5),

    _p("qaly_loss_iron_overload", 1.0, "QALYs per case", "derived",
       note="Clinically expressed haemochromatosis with organ involvement; "
            "most carriers never express, which is handled by penetrance "
            "upstream rather than by shrinking this figure.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.2, high=3.0),

    _p("qaly_loss_cancer", 2.5, "QALYs per case", "derived",
       note="Treated solid-tumour cancer across stages, stage-weighted; "
            "applied to the colorectal and breast/ovarian anchors.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.8, high=6.0),

    _p("qaly_loss_pathogenic_generic", 2.2, "QALYs per case", "derived",
       note="Fallback decrement for an actionable monogenic finding whose "
            "gene is not individually costed. The ACMG secondary-findings "
            "list is dominated by hereditary cancer and inherited cardiac "
            "conditions, so this sits between the two rather than borrowing "
            "either — a generic bucket should not silently inherit a "
            "specific condition's quality-of-life loss.",
       source="Sullivan PW, Ghushchyan V. Preference-Based EQ-5D Index Scores "
              "for Chronic Conditions in the United States. Med Decis Making",
       citation="PMID:16855125 doi:10.1177/0272989X06290495", year=2006,
       dist="gamma", low=0.5, high=6.0),

    # ── Prospective yield of sequencing over a genotyping chip ────────────
    # These answer "what is sequencing likely to find that my chip did not",
    # which is a different question from "what did sequencing find" and the
    # only one that can be answered before buying it.
    _p("wgs_yield_acmg_secondary", 0.020, "probability", "derived",
       note="Probability that an unselected adult carries a reportable "
            "pathogenic variant in an ACMG secondary-findings gene. "
            "Population screening cohorts report roughly 1–3% depending on "
            "how many genes are returned; 2% sits in the middle of that "
            "range. Genotyping arrays detect almost none of these, because "
            "they type specific known positions rather than reading the gene.",
       source="Grzymski JJ, Elhanan G, Morales Rosado JA, et al. Population "
              "genetic screening efficiently identifies carriers of "
              "autosomal dominant diseases. Nat Med",
       citation="PMID:33020650 doi:10.1038/s41591-020-1093-z", year=2020,
       dist="beta", low=0.008, high=0.040),

    _p("chip_detection_share_monogenic", 0.10, "fraction", "derived",
       note="Share of pathogenic variants in a given actionable gene that a "
            "consumer genotyping array actually detects. Arrays type a fixed "
            "set of positions — the authorised consumer BRCA report covers "
            "three founder variants out of thousands of known pathogenic "
            "ones — so a negative chip result is close to uninformative for "
            "monogenic disease. This is the single biggest reason sequencing "
            "adds anything.",
       source="US Food and Drug Administration. FDA authorises, with special "
              "controls, direct-to-consumer test for three BRCA1/BRCA2 "
              "breast cancer mutations (limitations of array-based reporting)",
       citation="https://www.fda.gov/news-events/press-announcements/"
                "fda-authorizes-special-controls-direct-consumer-test-"
                "reporting-three-mutations-brca-breast-cancer", year=2018,
       dist="beta", low=0.02, high=0.35),

    _p("chip_pgx_coverage", 0.85, "fraction", "assumption",
       note="Share of the clinically actionable pharmacogenomic variation a "
            "consumer array already captures. The main CPIC star alleles sit "
            "at common typed positions, so sequencing adds comparatively "
            "little here — which matters, because pharmacogenomics is the "
            "highest-prevalence category and would otherwise dominate a "
            "prospective estimate of what sequencing is worth. Judgement: "
            "coverage varies by array and by gene, and CYP2D6 in particular "
            "is poorly resolved by both technologies.",
       dist="beta", low=0.50, high=0.98),

    _p("wgs_yield_carrier_expanded", 0.60, "probability", "derived",
       note="Probability that an individual carries at least one recessive "
            "condition on an expanded carrier panel. Reported for "
            "completeness only — reproductive findings are deliberately not "
            "monetised anywhere in this model, so this contributes nothing "
            "to the dollar figure.",
       source="Haque IS, Lazarin GA, Kang HP, Evans EA, Goldberg JD, Wapner "
              "RJ. Modeled Fetal Risk of Genetic Diseases Identified by "
              "Expanded Carrier Screening. JAMA",
       citation="PMID:27532916 doi:10.1001/jama.2016.11139", year=2016,
       dist="beta", low=0.30, high=0.90),

    # ── Decision-layer parameters ─────────────────────────────────────────
    # These sat as default arguments in econ_decision.py, where each one
    # decided the sign or the substance of a reported result while being
    # invisible to the provenance count, the tornado and EVPPI. By the standard
    # applied to the rest of the model they belong here.
    _p("inequality_aversion", 11.0, "Atkinson parameter", "assumption",
       note="Strength of preference for reducing health inequality in the "
            "distributional analysis. Values around this magnitude appear in "
            "the English distributional cost-effectiveness literature, which "
            "elicits far stronger aversion than economic intuition suggests; "
            "this is recorded as a judgement rather than cited because the "
            "figure is carried here from memory of that literature and has "
            "not been verified against a specific paper. The result is "
            "sensitive to it, so it is reported as an input.",
       dist="gamma", low=0.0, high=30.0),

    _p("subgroup_annual_incidence", 0.01, "probability per year", "assumption",
       note="Annual incidence used for the illustrative age-and-sex subgroup "
            "table. The table exists to show how competing mortality changes "
            "the value of prevention with age, which it does at any plausible "
            "incidence; it is not a personalised risk estimate and is labelled "
            "as illustrative in the report.",
       dist="beta", low=0.001, high=0.05),

    # ── Default finding parameters ────────────────────────────────────────
    # The entire benefit side of the pooled model runs on these two numbers.
    # They were bare literals in value_of_information._collect, which meant
    # they could not be varied in sensitivity analysis and never appeared in
    # any provenance count — the most load-bearing figures in the model were
    # also the least visible. Registering them as declared assumptions is the
    # honest description, and it lets the tornado report their influence.
    _p("baseline_event_probability", 0.20, "probability", "assumption",
       note="Probability, over the modelled horizon, that a person with an "
            "elevated-risk finding experiences the condition it points to, "
            "where no condition-specific estimate is available. A generic "
            "stand-in for a quantity that genuinely varies by condition, age "
            "and sex; it drives every dollar on the benefit side, so it is "
            "varied widely in sensitivity analysis rather than defended.",
       dist="beta", low=0.05, high=0.45),

    _p("baseline_event_probability_dementia", 0.15, "probability", "assumption",
       note="As above, for dementia specifically, set lower than the generic "
            "default to reflect that the horizon ends before most of the "
            "lifetime risk is realised. Judgement, not an epidemiological "
            "estimate.",
       dist="beta", low=0.03, high=0.35),

    _p("actionable_rrr", 0.30, "relative risk reduction", "assumption",
       note="Risk reduction achievable by acting on a genomic finding where "
            "no trial-specific effect is available. Chosen to sit near the "
            "measured statin primary-prevention effect (0.27) as a plausible "
            "order of magnitude for a well-executed preventive response. This "
            "is the single most influential assumption in the model; the "
            "range spans doing appreciably less and appreciably better than a "
            "statin.",
       dist="beta", low=0.05, high=0.55),

    _p("intervention_cost_standard", 500.0, "$", "derived",
       note="Cost of acting on a typical risk finding: additional clinical "
            "visits and monitoring over the horizon, at US self-pay office-"
            "visit prices. An order-of-magnitude anchor, not a costed pathway.",
       source="Centers for Medicare & Medicaid Services, Physician Fee "
              "Schedule (office/outpatient evaluation and management codes)",
       citation="https://www.cms.gov/medicare/physician-fee-schedule/search",
       year=2025, dist="gamma", low=100.0, high=2_000.0),

    _p("intervention_cost_monogenic", 1_500.0, "$", "derived",
       note="Cost of acting on an actionable monogenic finding: confirmatory "
            "testing, specialist referral and enhanced surveillance. Larger "
            "than the standard anchor because the clinical response is "
            "specified by guideline rather than discretionary.",
       source="Miller DT, Lee K, Abul-Husn NS, et al. ACMG SF v3.2 list for "
              "reporting of secondary findings in clinical exome and genome "
              "sequencing. Genet Med",
       citation="PMID:37347242 doi:10.1016/j.gim.2023.100866", year=2023,
       dist="gamma", low=400.0, high=5_000.0),

    _p("intervention_cost_pgx", 100.0, "$", "derived",
       note="Cost of acting on a pharmacogenomic result: substituting or "
            "dose-adjusting a prescription, which is usually a change of drug "
            "rather than an added service.",
       source="US generic drug cash prices (GoodRx), retail observation",
       citation="https://www.goodrx.com/", year=2026,
       dist="gamma", low=0.0, high=600.0),

    _p("intervention_cost_predicted_variant", 800.0, "$", "assumption",
       note="Cost of following up a computationally predicted pathogenic "
            "variant. No standard pathway exists — that is the point of the "
            "category — so this is a judgement about what a cautious "
            "confirmatory workup would cost.",
       dist="gamma", low=200.0, high=3_000.0),

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

    # ── Adherence (efficacy → effectiveness) ─────────────────────────────
    # Every effect size in this model is trial efficacy: what happens when a
    # protocol is followed. What a payer buys is effectiveness: what happens
    # in a population where roughly half of people stop. Running the model at
    # implicit 100% adherence overstated every benefit by the size of that gap.
    #
    # These scale BOTH the effect and the ongoing intervention cost (real-world
    # effectiveness framing), not the effect alone (intention-to-treat). People
    # who stop taking a statin stop paying for it. The consequence is the
    # interesting one: the fixed test cost does NOT scale, so it amortises over
    # fewer realised QALYs and the ICER worsens even though the intervention's
    # own cost per unit of benefit is roughly unchanged.

    _p("adherence_pharmacological", 0.50, "proportion", "published",
       note="Long-term persistence with preventive medication. Applied to "
            "conditions whose modelled intervention is chronic drug therapy "
            "(statins, metformin, antidepressants). The 50% figure is the "
            "headline finding of the WHO review and is consistent with the "
            "statin discontinuation literature; varied 0.35-0.70.",
       source="World Health Organization. Adherence to Long-Term Therapies: "
              "Evidence for Action. Geneva: WHO",
       citation="https://iris.who.int/handle/10665/42682", year=2003,
       dist="beta", low=0.35, high=0.70),

    _p("adherence_screening", 0.65, "proportion", "assumption",
       note="Uptake of a recommended screening or surveillance programme "
            "(colonoscopy, mammography, ferritin monitoring). Higher than "
            "chronic drug therapy because the ask is episodic rather than "
            "daily, but well short of complete: US population screening "
            "uptake sits in the 60-70% band. Anchored on that band rather "
            "than on a single study, so it is declared as judgement.",
       dist="beta", low=0.45, high=0.85),

    _p("adherence_lifestyle", 0.35, "proportion", "assumption",
       note="Sustained behaviour change — diet, exercise, alcohol reduction — "
            "maintained over the horizon. The lowest of the three archetypes "
            "because maintenance, not initiation, is what the model needs, "
            "and maintenance attrition in behavioural trials is severe. "
            "Judgement; varied 0.20-0.55.",
       dist="beta", low=0.20, high=0.55),

    _p("adherence_default", 0.50, "proportion", "assumption",
       note="Fallback for a condition with no archetype assigned. Set equal "
            "to the pharmacological figure so an unmapped condition is "
            "treated no more optimistically than a mapped one. A condition "
            "reaching this parameter is a gap in the mapping, not a finding.",
       dist="beta", low=0.30, high=0.75),
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


class overridden:
    """Temporarily replace registry values, for sensitivity analysis.

    Probabilistic and one-way sensitivity analyses need to re-run the model
    with different parameter values. Doing that by passing overrides down
    through every call signature would be invasive and easy to get partially
    wrong — a parameter read from the registry somewhere deep would silently
    keep its base value and quietly narrow the reported uncertainty.

    Swapping the registry itself means every reader sees the drawn value, by
    construction. Restores on exit, including on exception.

        with overridden({"discount_rate": 0.05}):
            ...
    """

    def __init__(self, values: Dict[str, float]):
        self._values = dict(values or {})
        self._saved: Dict[str, Param] = {}

    def __enter__(self):
        import dataclasses
        for key, val in self._values.items():
            if key not in PARAMS:
                continue
            self._saved[key] = PARAMS[key]
            PARAMS[key] = dataclasses.replace(PARAMS[key], value=float(val))
        return self

    def __exit__(self, *exc):
        for key, param in self._saved.items():
            PARAMS[key] = param
        self._saved.clear()
        return False


def sampleable() -> List[Param]:
    """Parameters carrying enough information to be drawn from.

    A parameter with no distribution and no range is held fixed rather than
    given invented uncertainty — pretending to know a spread is the same error
    as pretending to know a value.
    """
    return [p for p in _REGISTRY
            if p.dist in ("beta", "gamma", "lognormal")
            and (p.se is not None or (p.low is not None and p.high is not None))]


def draw(rng, param: Param) -> float:
    """One draw from a parameter's documented uncertainty.

    Distributions follow the usual health-economic conventions: beta for
    quantities bounded on [0, 1], gamma for costs and other non-negative
    unbounded quantities. Where only a range is published, it is read as an
    approximate 95% interval and converted to a standard error.
    """
    v = float(param.value)
    se = param.se
    if se is None and param.low is not None and param.high is not None:
        se = (float(param.high) - float(param.low)) / 3.92   # 95% CI -> SE
    if not se or se <= 0:
        return v

    if param.dist == "beta":
        if not (0.0 < v < 1.0):
            return v
        var = min(se ** 2, v * (1 - v) * 0.999)
        if var <= 0:
            return v
        common = v * (1 - v) / var - 1.0
        a, b = max(1e-6, v * common), max(1e-6, (1 - v) * common)
        out = rng.betavariate(a, b)
    elif param.dist == "gamma":
        if v <= 0:
            return v
        shape = max(1e-6, (v / se) ** 2)
        out = rng.gammavariate(shape, v / shape)
    elif param.dist == "lognormal":
        if v <= 0:
            return v
        sigma = (se / v)
        out = math.exp(rng.gauss(math.log(v), max(1e-9, sigma)))
    else:
        return v

    if param.low is not None:
        out = max(float(param.low), out)
    if param.high is not None:
        out = min(float(param.high), out)
    return out


def sample_all(rng, keys: Optional[Sequence[str]] = None) -> Dict[str, float]:
    """A full parameter draw, for one probabilistic-sensitivity iteration."""
    pool = ([get(k) for k in keys] if keys is not None else sampleable())
    return {p.key: draw(rng, p) for p in pool}


def by_tier(tier: str) -> List[Param]:
    """All parameters at one provenance tier."""
    return [p for p in _REGISTRY if p.tier == tier]


def assumptions() -> List[Param]:
    """Parameters with no published anchor, for the report's declared-
    assumptions section."""
    return by_tier("assumption")


# ══════════════════════════════════════════════════════════════════════════
# Curated-table source resolution
# ══════════════════════════════════════════════════════════════════════════
# The per-finding tables in health_economics.py each carry a free-text ``src``
# like "Kazi et al. (2014) Ann Intern Med". Those are real attributions — the
# figures were not invented — but an author-year string cannot be resolved
# automatically, checked for existence, or followed by a reader without a
# search. This map upgrades them to PMIDs/DOIs in ONE place, keyed by the
# ``src`` string, so 300-odd table entries gain resolvable citations without
# editing 300-odd dictionaries.
#
# Only identifiers verified as matching the cited work appear here. An
# unlisted source is reported as "attributed, identifier not yet verified",
# which is the honest state and doubles as the work queue. A WRONG identifier
# would be worse than none: it sends a reader to the wrong paper while looking
# more rigorous, so guessing is not an option.
CURATED_SOURCE_IDS: Dict[str, str] = {
    # ── Cost-effectiveness and trial evidence ──
    "Ladabaum et al. (2011) Ann Intern Med":
        "PMID:21768580 doi:10.7326/0003-4819-155-2-201107190-00002",
    "Kazi et al. (2014) Ann Intern Med — genotype-guided antiplatelet":
        "PMID:25089860 doi:10.7326/M13-1999",
    "CTT Collaboration (2010) Lancet":
        "PMID:21067804 doi:10.1016/S0140-6736(10)61350-5",
    "DPP Research Group (2002) NEJM": "PMID:11832527 doi:10.1056/NEJMoa012512",
    "Knowler (2002) NEJM — DPP; Khera (2016) NEJM — PRS × lifestyle":
        "PMID:11832527 doi:10.1056/NEJMoa012512; "
        "PMID:27959714 doi:10.1056/NEJMoa1605086",
    "Yusuf et al. (2004) Lancet — INTERHEART":
        "PMID:15364185 doi:10.1016/S0140-6736(04)17018-9",
    "Mega et al. (2015) Lancet — PRS-guided statin benefit":
        "PMID:25748612 doi:10.1016/S0140-6736(14)61730-X",
    "Ngandu et al. (2015) Lancet — FINGER trial; Livingston (2020) Lancet "
    "— dementia prevention":
        "PMID:25771249 doi:10.1016/S0140-6736(15)60461-5; "
        "PMID:32738937 doi:10.1016/S0140-6736(20)30367-6",
    "Pashayan et al. (2018) Genet Med — PRS-stratified screening CEA":
        "PMID:29236091 doi:10.1038/gim.2017.246",

    # ── CPIC pharmacogenomic guidelines ──
    "Johnson et al. (2017) Clin Pharmacol Ther — CPIC warfarin guideline":
        "PMID:28198005 doi:10.1002/cpt.668",
    "Relling et al. (2019) Clin Pharmacol Ther — CPIC thiopurines":
        "PMID:30447069 doi:10.1002/cpt.1304",
    "Ramsey et al. (2014) Clin Pharmacol Ther — CPIC simvastatin/SLCO1B1":
        "PMID:24918167 doi:10.1038/clpt.2014.125",
    "Birdwell et al. (2015) Clin Pharmacol Ther — CPIC tacrolimus":
        "PMID:25801146 doi:10.1002/cpt.113",

    # ── GWAS and genetic-epidemiology anchors ──
    "Locke et al. (2015) Nature — BMI GWAS; NICE CG189 weight management":
        "PMID:25673413 doi:10.1038/nature14177",
    "Locke et al. (2015) Nature — BMI": "PMID:25673413 doi:10.1038/nature14177",
    "Okada et al. (2014) Nature — RA GWAS; Finckh (2006) Arthritis Rheum "
    "— early DMARD CEA": "PMID:24390342 doi:10.1038/nature12873",
    "Demenais et al. (2018) Nat Genet — asthma GWAS; GINA 2023 guidelines":
        "PMID:29273806 doi:10.1038/s41588-017-0014-7",
    "Howard et al. (2019) Nat Neurosci — MDD GWAS; Chisholm (2016) Lancet "
    "Psych — CBT CEA": "PMID:30718901 doi:10.1038/s41593-018-0326-7",
    "de Lange et al. (2017) Nat Genet — IBD GWAS; van der Valk (2016) IBD "
    "— IBD cost burden": "PMID:28067908 doi:10.1038/ng.3760",
    "Jones et al. (2019) Nat Commun — chronotype GWAS":
        "PMID:30696823 doi:10.1038/s41467-018-08259-7",
    "Samson et al. (1996) Nature — CCR5-Δ32": "PMID:8757135 doi:10.1038/382722a0",
    "Ge et al. (2009) Nature — IL28B × HCV treatment":
        "PMID:19684573 doi:10.1038/nature08309",
    "Garcia-Closas et al. (2005) Lancet — NAT2 × bladder cancer":
        "PMID:16144894 doi:10.1016/S0140-6736(05)67137-1",
    "Binder et al. (2008) JAMA Psych — FKBP5 × PTSD treatment response":
        "PMID:18349090 doi:10.1001/jama.299.11.1291",
    "Healy et al. (2008) Lancet Neurol — LRRK2 penetrance + PD trials":
        "PMID:18539534 doi:10.1016/S1474-4422(08)70117-0",
    "Adams et al. (2005) NEJM; Allen et al. (2008) NEJM — C282Y penetrance":
        "PMID:15858186 doi:10.1056/NEJMoa041534; "
        "PMID:18199861 doi:10.1056/NEJMoa073286",
    "Holick (2007) NEJM — vitamin D; WHO nutrition guidelines":
        "PMID:17634462 doi:10.1056/NEJMra070553",
    "Ridker (2003) Circulation — CRP":
        "PMID:12551878 doi:10.1161/01.CIR.0000053730.47739.3C",
}

_RESOLVABLE_RE = None


def resolve_curated_source(src: str) -> Dict[str, str]:
    """Classify one curated-table ``src`` string and attach an identifier.

    Returns ``{"text", "identifier", "state"}`` where ``state`` is one of
    ``resolvable`` (carries a PMID/DOI/URL a reader can follow),
    ``attributed`` (a real citation whose identifier is not yet verified), or
    ``missing`` (no attribution at all — the only genuinely bad case).
    """
    global _RESOLVABLE_RE
    if _RESOLVABLE_RE is None:
        import re
        _RESOLVABLE_RE = re.compile(
            r"(PMID:\s*\d+|doi:\s*10\.\S+|ISBN:[\dX-]+|https?://\S+)", re.I)
    text = (src or "").strip()
    if not text:
        return {"text": "", "identifier": "", "state": "missing"}
    if _RESOLVABLE_RE.search(text):
        return {"text": text, "identifier": text, "state": "resolvable"}
    ident = CURATED_SOURCE_IDS.get(text, "")
    if ident:
        return {"text": text, "identifier": ident, "state": "resolvable"}
    return {"text": text, "identifier": "", "state": "attributed"}


def audit_curated_tables() -> Dict:
    """Provenance state of every numeric field in the curated econ tables.

    The registry covers the model's spine — method conventions, cost-of-illness
    anchors, effect sizes. This covers the other several hundred numbers, and
    reporting both is what stops the provenance section from being a flattering
    statement about a small corner of the model.
    """
    try:
        import health_economics as _he
    except Exception:
        return {"available": False, "n_params": 0}
    fields = ("cost", "outcome_value", "prevalence", "qaly_gain",
              "adr_cost", "rrr")
    counts = {"resolvable": 0, "attributed": 0, "missing": 0}
    unresolved: Dict[str, int] = {}
    tables: List[Dict] = []
    for name in sorted(dir(_he)):
        if not (name.endswith("_ECONOMICS") or name.endswith("_COSTS")):
            continue
        table = getattr(_he, name, None)
        if not isinstance(table, dict):
            continue
        t_counts = {"resolvable": 0, "attributed": 0, "missing": 0}
        for entry in table.values():
            if not isinstance(entry, dict):
                continue
            n = sum(1 for f in fields if f in entry)
            if not n:
                continue
            src = (entry.get("src") or entry.get("source")
                   or entry.get("evidence") or "")
            state = resolve_curated_source(src)["state"]
            counts[state] += n
            t_counts[state] += n
            if state == "attributed":
                unresolved[src] = unresolved.get(src, 0) + n
        if sum(t_counts.values()):
            tables.append({"table": name, **t_counts})
    total = sum(counts.values()) or 1
    return {
        "available": True,
        "n_params": sum(counts.values()),
        "n_resolvable": counts["resolvable"],
        "n_attributed": counts["attributed"],
        "n_missing": counts["missing"],
        "pct_resolvable": round(100.0 * counts["resolvable"] / total, 1),
        "pct_attributed_or_better": round(
            100.0 * (counts["resolvable"] + counts["attributed"]) / total, 1),
        "tables": tables,
        # Work queue, largest first: which unresolved sources would buy the
        # most coverage if someone verified their identifier next.
        "unresolved_sources": [
            {"source": s, "n_params": n}
            for s, n in sorted(unresolved.items(), key=lambda kv: -kv[1])],
    }


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
    curated = audit_curated_tables()
    unregistered = curated.get("n_params", 0) if curated.get("available") \
        else count_unregistered_parameters()
    total_known = n + unregistered

    # Whole-model view. Registry parameters are all sourced or explicitly
    # declared; curated-table parameters split by whether their attribution
    # carries a resolvable identifier.
    model_resolvable = (counts["published"] + counts["derived"]
                        + curated.get("n_resolvable", 0))
    model_attributed = curated.get("n_attributed", 0)
    model_unsourced = counts["assumption"] + curated.get("n_missing", 0)
    denom = total_known or 1
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
        "pct_of_model_registered": round(100.0 * n / denom, 1),
        # Whole-model provenance, which is the figure worth quoting.
        "n_curated_resolvable": curated.get("n_resolvable", 0),
        "n_curated_attributed": curated.get("n_attributed", 0),
        "n_curated_missing": curated.get("n_missing", 0),
        "model_pct_resolvable": round(100.0 * model_resolvable / denom, 1),
        "model_pct_attributed_or_better": round(
            100.0 * (model_resolvable + model_attributed) / denom, 1),
        "model_pct_unsourced": round(100.0 * model_unsourced / denom, 1),
        "n_unresolved_sources": len(curated.get("unresolved_sources", [])),
        "unresolved_sources": curated.get("unresolved_sources", [])[:12],
        "scope": ("Registered parameters cover method conventions, "
                  "cost-of-illness anchors, effect sizes and utilities. The "
                  "per-finding curated tables in health_economics.py are not "
                  "in the registry, but every one carries a literature "
                  "attribution; those whose identifier has been verified are "
                  "counted as resolvable, the rest as attributed."),
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
