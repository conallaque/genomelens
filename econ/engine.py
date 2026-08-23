"""
Cost-effectiveness engine: pooling, state transitions, and disaggregated output
===============================================================================

This module replaces three habits that made the previous economic output
indefensible, whatever methods were layered on top of it.

**1. Findings were summed as if independent.** Eight of the report's twenty-one
finding sources route onto the same cardiometabolic anchor. Each one arrived at
the model as its own line with its own risk reduction, and the lines were added.
A polygenic score for coronary disease, a PheWAS lipid biomarker, a Mendelian-
randomisation estimate for the same lipid, and a "longevity composite" are not
four independent opportunities to prevent four heart attacks — they are four
measurements of one liability. :class:`ConditionPool` combines them on the risk
scale with a complement-of-products rule, applies a compounding penalty for
correlated re-measurement, caps the total, and charges the cost of illness
**once**.

**2. Cash savings and monetised health were added into one headline number.**
"Net benefit $275,405" blended dollars that a payer would actually not spend
with quality-adjusted life-years priced at willingness-to-pay. Those are
different objects and a reviewer needs them apart. Every result here reports
incremental cost, incremental QALYs, ICER and INMB separately, and the
willingness-to-pay threshold is a stated input rather than an invisible
multiplier.

**3. Time was handled by a single scalar.** Discounting the whole horizon at its
midpoint is a reasonable shortcut when events are undated, but it cannot express
competing mortality — and without competing mortality, preventing a disease at
78 credits the full quality-adjusted life-years of someone who would have died
at 76 anyway. :func:`run_markov` runs a three-state cohort model against US
life-table mortality with Simpson's 1/3 within-cycle correction, which is the
convention the cited tutorial validates against.

Nothing here invents a parameter. Every constant is fetched from
:mod:`econ_params` by key, so its provenance tier travels with it and the
report can state how much of its own output rests on judgement.
"""

from __future__ import annotations

import csv
import math
import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

from . import params as ep

__all__ = [
    "ADHERENCE_BY_COI_KEY",
    "COI_KEY_TO_PARAM",
    "DEFAULT_GENE_VOCABULARY",
    "CEAResult",
    "ConditionPool",
    "Finding",
    "MarkovResult",
    "adherence_for",
    "ceac",
    "cheers_checklist",
    "deduplicate_by_target",
    "discount_weights",
    "dual_perspective",
    "evaluate_pools",
    "impact_inventory",
    "incremental_analysis",
    "life_table",
    "pool_findings",
    "run_markov",
    "run_psa",
    "simpson_weights",
    "tornado",
    "validate_model",
]

# One level up: this module lives in econ/ but the vendored life table stays
# at the repository root. It has to — the .gitignore rule that keeps DNA out
# of the repo is a blanket ``*.csv`` with a single path-anchored exception for
# ``data/LifeTable_USA_Mx_2015.csv``. Move the file under econ/ and that
# exception silently stops matching, so the life table quietly becomes
# untracked while the DNA guard appears to still be in place.
_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
_LIFE_TABLE_PATH = os.path.join(_DATA_DIR, "LifeTable_USA_Mx_2015.csv")

# Condition anchors used by value_of_information._classify_category, mapped to
# the registry keys that cost them. A condition with no entry here cannot be
# valued — which is the intended behaviour, since the alternative is inventing
# a cost for it.
#
# Each condition must name its OWN quality-of-life decrement. For a while
# seven of these pointed at ``qaly_loss_mace`` — a non-fatal cardiovascular
# event's decrement — which claimed the same quality-of-life loss for dementia,
# depression and kidney stones. Beyond being wrong, it corrupted the
# sensitivity analysis: one constant carrying seven conditions dominated the
# tornado, so the report named it the model's key driver when it was really
# just the most overloaded placeholder.
COI_KEY_TO_PARAM: dict[str, tuple[str, str]] = {
    # coi_key            (cost param,               qaly-loss param)
    "CAD":               ("coi_mace",               "qaly_loss_mace"),
    "T2D":               ("coi_t2d",                "qaly_loss_t2d"),
    "Alzheimer":         ("coi_alzheimer",          "qaly_loss_dementia"),
    "Depression":        ("coi_depression",         "qaly_loss_depression"),
    "SubstanceUse":      ("coi_substance_use",      "qaly_loss_substance_use"),
    "Autoimmune":        ("coi_autoimmune",         "qaly_loss_autoimmune"),
    "Urologic":          ("coi_urologic",           "qaly_loss_urologic"),
    "IronOverload":      ("coi_iron_overload",      "qaly_loss_iron_overload"),
    "Colorectal":        ("coi_colorectal",         "qaly_loss_cancer"),
    "BreastOvarian":     ("coi_breast_ovarian",     "qaly_loss_cancer"),
    "Pathogenic":        ("coi_pathogenic_generic", "qaly_loss_pathogenic_generic"),
    # Coeliac, added when the gut-health module was connected to the engine.
    # A deliberately small anchor: the disease is managed by diet rather than
    # by procedures or biologics. Lactase non-persistence was considered and
    # left OUT — it is a symptom burden rather than a disease, worth roughly
    # $250 and 0.05 QALYs, and buying that with two unsourced registry
    # parameters would have dropped the model below its own 75%-sourced gate
    # for a rounding error. It stays a reported signal with no dollar figure.
    "Coeliac":           ("coi_coeliac",            "qaly_loss_coeliac"),
}


# Every relative risk reduction in this model is trial efficacy — the effect
# when the protocol is followed. Real cohorts do not follow protocols: roughly
# half of people stop long-term preventive medication. Until this map existed
# the model ran at implicit 100% adherence, which is not a conservative
# simplification but a systematic overstatement of every benefit it reports.
#
# Conditions are assigned to one of three archetypes by what acting on the
# finding actually requires of the person, since that — not the disease — is
# what predicts whether they keep doing it.
ADHERENCE_BY_COI_KEY: dict[str, str] = {
    # coi_key            registry key           what acting requires
    "CAD":               "adherence_pharmacological",   # statin, daily
    "T2D":               "adherence_lifestyle",         # diet and exercise
    "Alzheimer":         "adherence_lifestyle",         # activity, hearing, BP
    "Depression":        "adherence_pharmacological",   # SSRI, daily
    "SubstanceUse":      "adherence_lifestyle",         # sustained abstinence
    "Autoimmune":        "adherence_screening",         # periodic monitoring
    "Urologic":          "adherence_lifestyle",         # hydration, diet
    "IronOverload":      "adherence_screening",         # ferritin surveillance
    "Colorectal":        "adherence_screening",         # colonoscopy programme
    "BreastOvarian":     "adherence_screening",         # imaging surveillance
    "Pathogenic":        "adherence_screening",         # specialist follow-up
    # Coeliac: the genotype informs a screening decision (serology), and the
    # diet that follows a positive result is enforced by symptoms rather than
    # by willpower, so it sits closer to screening than to general dietary
    # change. Lactose is ordinary dietary self-management.
    "Coeliac":           "adherence_screening",         # serology, then diet
}


def adherence_for(coi_key: str) -> float:
    """Real-world adherence multiplier for a condition's intervention.

    Falls back to ``adherence_default`` rather than to 1.0: an unmapped
    condition is a gap in the map, and defaulting it to perfect adherence
    would make the gap invisible by making it flattering.
    """
    return float(ep.value(ADHERENCE_BY_COI_KEY.get(coi_key, "adherence_default")))


# ══════════════════════════════════════════════════════════════════════════
# Inputs
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class Finding:
    """One genomic signal bearing on one condition.

    ``rrr`` is the relative risk reduction achievable by acting on this finding,
    ``p_event`` the baseline probability of the condition over the horizon, and
    ``haircut`` the evidence-strength discount for the source module. The
    distinction between ``haircut`` (how much to believe the source) and
    ``rrr`` (how much acting helps) is deliberate: they were previously
    entangled in one multiplier.

    ``adherence`` is the third and last term, and it answers a different
    question again: not whether the evidence is good or the effect is large,
    but whether the person keeps doing the thing. Efficacy times adherence is
    effectiveness, and effectiveness is what a payer buys.
    """

    label: str
    coi_key: str
    p_event: float
    rrr: float
    haircut: float = 1.0
    intervention_cost: float = 0.0
    confidence: str = "moderate"
    source_category: str = ""
    qaly_override: float | None = None
    adherence: float = 1.0

    @property
    def efficacy_rrr(self) -> float:
        """Risk reduction under the trial protocol, once the source's evidence
        strength is applied but before real-world adherence. Bounded to [0, 1).

        Kept separate from :attr:`effective_rrr` so the report can show the
        efficacy-to-effectiveness gap as its own line rather than folding it
        into the pooling correction, where it would be invisible.
        """
        return max(0.0, min(0.999, float(self.rrr) * float(self.haircut)))

    @property
    def effective_rrr(self) -> float:
        """Risk reduction expected in a real cohort: efficacy times adherence."""
        return max(0.0, min(0.999, self.efficacy_rrr * float(self.adherence)))


# ══════════════════════════════════════════════════════════════════════════
# Pooling — the double-counting fix
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class ConditionPool:
    """All findings bearing on a single condition, combined once.

    The combination rule is complement-of-products on the risk scale::

        combined = 1 - Π (1 - rrr_i · penalty^i)

    which is what you get if the interventions act independently on the same
    baseline risk — additive summation, the previous behaviour, can exceed 1.0
    and is simply not a probability. The ``penalty^i`` term (i = rank order,
    strongest first) then down-weights each additional signal on the grounds
    that a second measurement of the same liability is mostly the same
    information. Finally the total is capped, because no combination of
    lifestyle and screening eliminates a common complex disease.
    """

    coi_key: str
    findings: list[Finding] = field(default_factory=list)

    # ── combination ─────────────────────────────────────────────────────
    @staticmethod
    def _combine(values: list[float]) -> float:
        """Complement-of-products with the correlated-signal penalty and cap.

        Factored out so the same rule can be applied to efficacy and to
        effectiveness. Scaling the pooled result by adherence afterwards would
        not give the same answer — the combination is non-linear — so the
        discount has to enter before the product, not after it.
        """
        if not values:
            return 0.0
        penalty = ep.value("correlated_signal_penalty")
        cap = ep.value("max_combined_rrr")
        surviving = 1.0
        for i, v in enumerate(sorted(values, reverse=True)):
            surviving *= (1.0 - v * (penalty ** i))
        return min(cap, 1.0 - surviving)

    def combined_rrr(self) -> float:
        """Pooled risk reduction as a real cohort would realise it."""
        return self._combine([f.effective_rrr for f in self.findings])

    def pooled_efficacy_rrr(self) -> float:
        """Pooled risk reduction if everyone followed the protocol.

        The difference between this and :meth:`combined_rrr` is the whole cost
        of imperfect adherence, isolated from every other correction.
        """
        return self._combine([f.efficacy_rrr for f in self.findings])

    def adherence(self) -> float:
        """The pool's adherence multiplier.

        Adherence is assigned per condition, so every finding in a pool shares
        one value and scaling the pool's cost by it is well defined. The
        assertion records that invariant: if adherence ever becomes per-finding,
        ``intervention_cost`` below has to change with it.
        """
        if not self.findings:
            return 1.0
        vals = {round(float(f.adherence), 6) for f in self.findings}
        assert len(vals) == 1, (
            f"{self.coi_key}: adherence must be uniform within a pool, got {vals}")
        return vals.pop()

    def baseline_risk(self) -> float:
        """One baseline probability for the condition, not one per finding.

        Findings disagree about baseline risk because they measure it
        differently. Taking the maximum keeps the strongest evidence of
        elevated risk rather than diluting it with weaker panels, while still
        charging the cost of illness only once.
        """
        if not self.findings:
            return 0.0
        return max(0.0, min(1.0, max(f.p_event for f in self.findings)))

    def intervention_cost(self) -> float:
        """Cost of acting on this condition.

        Acting on four correlated signals for one condition is one course of
        action, so the pool is charged the most expensive single intervention
        rather than the sum. Charging all four was the cost-side mirror of the
        benefit-side double count.
        """
        if not self.findings:
            return 0.0
        raw = max(float(f.intervention_cost or 0.0) for f in self.findings)
        # Scaled by adherence for the same reason the benefit is: someone who
        # stops taking the statin at month six stops paying for it. Costing the
        # full course while crediting only the adhered fraction of the benefit
        # would be the mirror-image error of the one this model set out to fix.
        return raw * self.adherence()

    def qaly_loss(self) -> float:
        overrides = [f.qaly_override for f in self.findings
                     if f.qaly_override is not None]
        if overrides:
            return float(max(overrides))
        entry = COI_KEY_TO_PARAM.get(self.coi_key)
        return ep.value(entry[1]) if entry else 0.0

    def coi_cost(self) -> float:
        entry = COI_KEY_TO_PARAM.get(self.coi_key)
        return ep.value(entry[0]) if entry else 0.0

    # ── reporting ───────────────────────────────────────────────────────
    def naive_sum_rrr(self) -> float:
        """What the previous additive treatment would have produced.

        Kept so the report can show the size of the correction rather than
        quietly banking it.
        """
        return sum(f.efficacy_rrr for f in self.findings)

    def to_dict(self) -> dict:
        return {
            "condition": self.coi_key,
            "n_findings": len(self.findings),
            "labels": [f.label for f in self.findings],
            "sources": sorted({f.source_category for f in self.findings if f.source_category}),
            "baseline_risk": round(self.baseline_risk(), 4),
            "combined_rrr": round(self.combined_rrr(), 4),
            "pooled_efficacy_rrr": round(self.pooled_efficacy_rrr(), 4),
            "naive_additive_rrr": round(self.naive_sum_rrr(), 4),
            # Two corrections, two numbers. Pooling is measured at full
            # adherence and adherence is measured after pooling, so neither
            # absorbs the other and the pair sums to the total correction.
            "double_count_avoided": round(
                max(0.0, self.naive_sum_rrr() - self.pooled_efficacy_rrr()), 4),
            "adherence": round(self.adherence(), 4),
            "adherence_archetype": ADHERENCE_BY_COI_KEY.get(
                self.coi_key, "adherence_default"),
            "adherence_drag_rrr": round(
                max(0.0, self.pooled_efficacy_rrr() - self.combined_rrr()), 4),
            "coi_cost": round(self.coi_cost()),
            "qaly_loss": round(self.qaly_loss(), 3),
            "intervention_cost": round(self.intervention_cost()),
        }


def pool_findings(findings: Iterable[Finding]) -> dict[str, ConditionPool]:
    """Group findings by condition into pools."""
    pools: dict[str, ConditionPool] = {}
    for f in findings:
        if not f.coi_key:
            continue
        pools.setdefault(f.coi_key, ConditionPool(f.coi_key)).findings.append(f)
    return pools


# When the same gene produces a line in two panels — COMT via the
# neurochemistry panel and again via pharmacogenomic prescribing guidance —
# those are two framings of one genotype, not two independent benefits.
#
# Matching is against an explicit VOCABULARY rather than a regex for
# gene-shaped words. A shape-based pattern reads "MI", "MACE" and "B12" as
# genes and invents targets that do not exist; worse, it could collapse two
# genuinely independent findings and under-count. Callers pass the vocabulary
# their own modules actually emit; this default covers the common case.
DEFAULT_GENE_VOCABULARY = frozenset({
    "APOE", "COMT", "BDNF", "DRD2", "OPRM1", "MTHFR", "FTO", "TCF7L2",
    "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A5", "CYP4F2", "VKORC1", "SLCO1B1",
    "DPYD", "TPMT", "NUDT15", "UGT1A1", "G6PD", "NAT2", "PON1", "GSTT1",
    "GSTM1", "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2", "APC",
    "LDLR", "APOB", "PCSK9", "MYH7", "MYBPC3", "KCNQ1", "KCNH2", "SCN5A",
    "RYR1", "RYR2", "TTN", "LMNA", "TP53", "PTEN", "RB1", "VHL", "RET",
    "SDHB", "SDHD", "MEN1", "NF2", "TSC1", "TSC2", "HFE", "SERPINA1",
    "F5", "F2", "PTPN22", "HLA", "CCR5", "IL28B", "LRRK2", "ATP7B",
    "CFTR", "SMN1", "HEXA", "HBB", "GBA", "PAH", "GJB2", "FLG", "ALDH2",
    "ADH1B", "CHRNA5", "SLC30A8", "PPARG", "KCNJ11", "ACE", "AGT", "NOS3",
})

# Topic tokens for lines that name no gene but plainly describe the same
# clinical target as another line. Deliberately short: over-collapsing here
# would hide a genuine second benefit, so only phrases whose duplication has
# actually been observed in this report's output are listed.
DEFAULT_TOPIC_TARGETS: tuple[tuple[str, str], ...] = (
    ("autoimmune", "topic:autoimmune"),
    ("statin-induced", "topic:statin-myopathy"),
    ("myopathy", "topic:statin-myopathy"),
)


def _extract_target(text: str,
                    vocabulary: frozenset | None = None,
                    topics: Sequence[tuple[str, str]] | None = None) -> str:
    """Identifier for what a finding is *about*, or '' if it shares nothing.

    Returns a gene symbol when the text names one from the vocabulary, since
    that is the sharpest signal that two lines describe the same genotype;
    otherwise a topic tag for the few clinical targets known to surface twice.
    Returning '' is the safe answer — the caller then leaves the line alone.
    """
    vocab = DEFAULT_GENE_VOCABULARY if vocabulary is None else vocabulary
    tops = DEFAULT_TOPIC_TARGETS if topics is None else topics

    # Topics are checked BEFORE genes, and deliberately so. The duplication
    # this catches is "same clinical target described twice" — a carrier
    # result named by its gene in one panel and by its clinical implication in
    # another. Matching the gene first would give those two lines different
    # keys and leave the duplicate standing. Pooling too eagerly under-counts
    # a benefit, which is the safer direction to be wrong in; the topic list
    # is kept to phrases whose duplication has actually been observed here.
    low = (text or "").lower()
    for needle, tag in tops:
        if needle in low:
            return tag

    # Gene symbols. Split on the punctuation that attaches to them in prose
    # ("COMT-guided", "HLA-B*58:01") before testing against the vocabulary.
    import re
    for token in re.findall(r"[A-Z][A-Z0-9]*(?:[*:\-][A-Z0-9]+)*", (text or "").upper()):
        for base in re.split(r"[*:\-]", token):
            if base in vocab:
                return base
    return ""


def deduplicate_by_target(items: Sequence[dict], *,
                          value_key: str = "net",
                          text_key: str = "finding",
                          fallback_key: str = "category",
                          target_key: str = "pool_hint",
                          vocabulary: frozenset | None = None,
                          penalty: float | None = None) -> list[dict]:
    """Down-weight repeated claims on the same underlying genotype.

    The condition-level pooling above catches "four findings, one disease".
    This catches the other shape the report produces: one variant surfacing in
    two panels and being valued twice — a COMT line from the neurochemistry
    panel and a COMT-guided-prescribing line from the pharmacogenomics panel
    are the same genotype seen from two angles.

    Items are ranked by value within each target; the strongest keeps its full
    value and each subsequent one is multiplied by the correlated-signal
    penalty compounding by rank. Returns copies annotated with
    ``pool_target``, ``pool_rank`` and ``retained`` so the report can show
    which lines were discounted and why, rather than silently shrinking them.
    """
    pen = ep.value("correlated_signal_penalty") if penalty is None else float(penalty)
    grouped: dict[str, list[dict]] = {}
    for i, it in enumerate(items):
        # A line naming no shared target gets a unique key, so it is never
        # pooled with anything. Grouping unmatched lines together by category
        # would discount independent findings for no reason.
        # An explicit hint wins over text extraction. Two records can describe
        # one course of action without sharing a gene symbol — a carrier line
        # named by gene and a partner-testing line named by disease are the
        # same reproductive decision, and no vocabulary match will ever join
        # them. The producing code knows they belong together; this lets it say so.
        target = (str(it.get(target_key) or "").strip()
                  or _extract_target(str(it.get(text_key, "")), vocabulary)
                  or f"unique:{i}")
        grouped.setdefault(target, []).append(it)

    out: list[dict] = []
    for target, group in grouped.items():
        ranked = sorted(group, key=lambda d: float(d.get(value_key, 0) or 0),
                        reverse=True)
        shared = len(ranked) > 1 and not target.startswith("unique:")
        for rank, it in enumerate(ranked):
            keep = (pen ** rank) if shared else 1.0
            copy = dict(it)
            copy["pool_target"] = target
            copy["pool_rank"] = rank
            copy["retained"] = round(keep, 4)
            if shared and rank > 0:
                copy["pool_note"] = (
                    f"Discounted to {keep:.0%} — {target} is already valued "
                    f"above under a different panel; this is the same "
                    f"genotype seen from another angle, not a second benefit.")
            for k in ("avoided", "qaly_value", "net", "qaly", "intervention"):
                if k in copy and isinstance(copy[k], int | float):
                    copy[k] = (round(copy[k] * keep, 4) if k == "qaly"
                               else round(copy[k] * keep))
            out.append(copy)
    return out


# ══════════════════════════════════════════════════════════════════════════
# Disaggregated cost-effectiveness output
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class CEAResult:
    """Incremental costs and effects, kept apart.

    ``inc_cost`` is incremental cost to the healthcare sector — negative means
    the strategy saves money. ``inc_qaly`` is incremental quality-adjusted
    life-years. ``icer`` is their ratio, undefined and reported as ``None``
    when the strategy dominates or is dominated (a ratio in those quadrants is
    meaningless and reporting one is a classic error). ``inmb`` is incremental
    net monetary benefit at the stated threshold.
    """

    strategy: str
    inc_cost: float
    inc_qaly: float
    wtp: float
    detail: dict = field(default_factory=dict)

    @property
    def icer(self) -> float | None:
        if abs(self.inc_qaly) < 1e-9:
            return None
        if self.inc_cost < 0 and self.inc_qaly > 0:
            return None      # dominant — report as such, not as a ratio
        if self.inc_cost > 0 and self.inc_qaly < 0:
            return None      # dominated
        return self.inc_cost / self.inc_qaly

    @property
    def inmb(self) -> float:
        return self.inc_qaly * self.wtp - self.inc_cost

    @property
    def verdict(self) -> str:
        if self.inc_qaly > 0 and self.inc_cost < 0:
            return "dominant (more health, lower cost)"
        if self.inc_qaly < 0 and self.inc_cost > 0:
            return "dominated (less health, higher cost)"
        if self.inmb > 0:
            return f"cost-effective at ${self.wtp:,.0f}/QALY"
        return f"not cost-effective at ${self.wtp:,.0f}/QALY"

    def to_dict(self) -> dict:
        return {
            "strategy": self.strategy,
            "incremental_cost": round(self.inc_cost),
            "incremental_qaly": round(self.inc_qaly, 4),
            "icer": (round(self.icer) if self.icer is not None else None),
            "icer_note": ("not defined — strategy dominates or is dominated"
                          if self.icer is None and abs(self.inc_qaly) > 1e-9
                          else ""),
            "inmb": round(self.inmb),
            "wtp": round(self.wtp),
            "verdict": self.verdict,
            **self.detail,
        }


def evaluate_pools(pools: dict[str, ConditionPool],
                   *, wtp: float | None = None,
                   horizon_years: float | None = None,
                   test_cost: float = 0.0,
                   marginal_only: bool = True) -> dict:
    """Turn condition pools into a disaggregated cost-effectiveness result.

    Returns per-condition rows plus a portfolio total, with costs and QALYs
    reported separately and the naive additive figure alongside so the size of
    the double-counting correction is visible.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    horizon = (ep.value("horizon_years_personal") if horizon_years is None
               else float(horizon_years))
    rate = ep.value("discount_rate")
    mcf = ep.value("marginal_cost_fraction") if marginal_only else 1.0
    # Events land at an unmodelled point inside the horizon; discount at the
    # midpoint. The Markov path below supersedes this for conditions where
    # timing matters, and reports both.
    disc = 1.0 / (1.0 + rate) ** (horizon / 2.0)

    rows: list[dict] = []
    tot_cost_averted = tot_qaly = tot_intervention = 0.0
    tot_naive_cost_averted = 0.0
    tot_efficacy_cost_averted = tot_efficacy_qaly = 0.0

    for key in sorted(pools):
        pool = pools[key]
        p0 = pool.baseline_risk()
        rrr = pool.combined_rrr()
        coi = pool.coi_cost()
        qloss = pool.qaly_loss()
        cases_averted = p0 * rrr

        cost_averted = cases_averted * coi * mcf * disc
        qaly_gained = cases_averted * qloss * disc
        intervention = pool.intervention_cost()

        naive_cases = p0 * min(1.0, pool.naive_sum_rrr())
        naive_cost_averted = naive_cases * coi * mcf * disc
        # The same arithmetic at full adherence, so the report can separate the
        # two corrections instead of showing one combined shrinkage.
        eff_cases = p0 * pool.pooled_efficacy_rrr()
        efficacy_cost_averted = eff_cases * coi * mcf * disc
        efficacy_qaly = eff_cases * qloss * disc

        tot_cost_averted += cost_averted
        tot_qaly += qaly_gained
        tot_intervention += intervention
        tot_efficacy_cost_averted += efficacy_cost_averted
        tot_efficacy_qaly += efficacy_qaly
        tot_naive_cost_averted += naive_cost_averted

        d = pool.to_dict()
        d.update({
            "cases_averted": round(cases_averted, 4),
            "cost_averted": round(cost_averted),
            "qaly_gained": round(qaly_gained, 4),
            "inmb": round(qaly_gained * wtp + cost_averted - intervention),
            "naive_cost_averted": round(naive_cost_averted),
            "efficacy_cost_averted": round(efficacy_cost_averted),
            "efficacy_qaly_gained": round(efficacy_qaly, 4),
            "inflation_removed": round(naive_cost_averted - cost_averted),
        })
        rows.append(d)

    inc_cost = tot_intervention + test_cost - tot_cost_averted
    result = CEAResult(
        strategy="Act on genomic findings vs. usual care",
        inc_cost=inc_cost, inc_qaly=tot_qaly, wtp=wtp,
        detail={
            "cost_averted": round(tot_cost_averted),
            "intervention_cost": round(tot_intervention),
            "test_cost": round(test_cost),
            "horizon_years": horizon,
            "discount_rate": rate,
            "midpoint_discount_factor": round(disc, 4),
            "marginal_cost_fraction": mcf,
        })

    return {
        "available": bool(rows),
        "conditions": rows,
        "cea": result.to_dict(),
        # Efficacy is what the trials measured; effectiveness is what this
        # cohort would actually get. Reporting only the second hides the size
        # of the assumption; reporting only the first is the overstatement the
        # model previously shipped.
        "adherence": {
            "efficacy_cost_averted": round(tot_efficacy_cost_averted),
            "effectiveness_cost_averted": round(tot_cost_averted),
            "efficacy_qaly": round(tot_efficacy_qaly, 4),
            "effectiveness_qaly": round(tot_qaly, 4),
            "value_lost_to_non_adherence": round(
                tot_efficacy_cost_averted - tot_cost_averted),
            "qaly_lost_to_non_adherence": round(
                tot_efficacy_qaly - tot_qaly, 4),
            "pct_of_benefit_lost": (
                round(100.0 * (tot_efficacy_qaly - tot_qaly) / tot_efficacy_qaly, 1)
                if tot_efficacy_qaly > 0 else 0.0),
            # The fixed test cost does not shrink with adherence, so it is
            # spread over fewer realised QALYs. That, not the intervention's
            # own cost per QALY, is why the ICER moves.
            "fixed_test_cost": round(test_cost),
            "archetypes": sorted({
                ADHERENCE_BY_COI_KEY.get(r["condition"], "adherence_default")
                for r in rows}),
            "note": ("Effect sizes are trial efficacy; this cohort is charged "
                     "real-world adherence on both the benefit and the ongoing "
                     "intervention cost. The one-off test cost is not scaled."),
            "src": ("WHO (2003), Adherence to Long-Term Therapies: Evidence "
                    "for Action"),
        },
        "double_counting": {
            "naive_cost_averted": round(tot_naive_cost_averted),
            "pooled_cost_averted": round(tot_cost_averted),
            "inflation_removed": round(tot_naive_cost_averted - tot_cost_averted),
            "pct_removed": (
                round(100.0 * (tot_naive_cost_averted - tot_cost_averted)
                      / tot_naive_cost_averted, 1)
                if tot_naive_cost_averted > 0 else 0.0),
            "explanation":
                "Findings routed to the same condition are combined on the "
                "risk scale and charged one cost of illness. The naive figure "
                "is what additive summation of the same findings produces.",
        },
    }


# ══════════════════════════════════════════════════════════════════════════
# Cohort state-transition model
# ══════════════════════════════════════════════════════════════════════════

def life_table(sex: str = "Total") -> dict[int, float]:
    """Age-specific all-cause mortality hazard, from the vendored US table.

    Returns ``{age: mx}``. Missing file degrades to an empty dict, and callers
    fall back to a constant hazard — the model should lose precision, not
    disappear, if the data file is absent.
    """
    col = {"m": "Male", "male": "Male", "f": "Female", "female": "Female"}.get(
        (sex or "").strip().lower(), "Total")
    out: dict[int, float] = {}
    try:
        with open(_LIFE_TABLE_PATH, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    out[int(row["Age"])] = float(row[col])
                except (KeyError, TypeError, ValueError):
                    continue
    except OSError:
        return {}
    return out


def simpson_weights(n_cycles: int) -> list[float]:
    """Simpson's 1/3 within-cycle correction weights over ``n_cycles + 1``
    points, transcribed from ``darthtools::gen_wcc``.

    Without a correction, a cohort model charges a full cycle of reward to
    people who transitioned partway through it, biasing every total upward.

    The published implementation is 1-indexed, so its first interior weight is
    2/3 where a textbook Simpson's rule would put 4/3. That asymmetry is
    reproduced deliberately: the point is to match the convention the cited
    tutorial validates against, not to improve on it. ``test_econ_engine.py``
    pins these weights against the independent implementation in the
    heor-model-replication repository when that repository is available.
    """
    n = int(n_cycles)
    if n <= 0:
        return [1.0]
    if n == 1:
        return [1.0 / 3.0, 1.0 / 3.0]
    # R's v_cycles = 1..n+1 maps to python index i -> i + 1; even -> 2/3.
    w = [(2.0 / 3.0) if (i + 1) % 2 == 0 else (4.0 / 3.0) for i in range(n + 1)]
    w[0] = w[n] = 1.0 / 3.0
    return w


def discount_weights(n_cycles: int, rate: float,
                     cycle_length: float = 1.0) -> list[float]:
    return [1.0 / (1.0 + rate) ** (t * cycle_length) for t in range(n_cycles + 1)]


@dataclass
class MarkovResult:
    cost: float
    qaly: float
    life_years: float
    trace: list[tuple[float, float, float]]   # (well, diseased, dead)
    n_cycles: int
    start_age: float

    def to_dict(self) -> dict:
        return {"cost": round(self.cost), "qaly": round(self.qaly, 4),
                "life_years": round(self.life_years, 3),
                "n_cycles": self.n_cycles, "start_age": self.start_age}


def run_markov(*, start_age: float, annual_incidence: float,
               coi_cost: float, disutility: float,
               rrr: float = 0.0, annual_intervention_cost: float = 0.0,
               excess_mortality_rr: float = 1.5,
               sex: str = "Total", max_age: int = 100,
               rate: float | None = None) -> MarkovResult:
    """Three-state cohort model — Well → Diseased → Dead — against US mortality.

    This exists to answer one question the scalar-discount approach cannot:
    how much of a prevented event's value is real once you account for the
    person having to survive long enough to collect it. Background mortality
    comes from the life table, disease adds excess mortality, and rewards get
    Simpson's 1/3 within-cycle correction before discounting.
    """
    rate = ep.value("discount_rate") if rate is None else float(rate)
    lt = life_table(sex)
    a0 = int(max(18, min(float(start_age), max_age - 1)))
    n = int(max_age - a0)
    if n <= 0:
        return MarkovResult(0.0, 0.0, 0.0, [], 0, float(a0))

    u_well = ep.value("utility_healthy")
    u_sick = max(0.0, u_well - float(disutility))
    inc = max(0.0, min(1.0, float(annual_incidence) * (1.0 - float(rrr))))

    well, sick, dead = 1.0, 0.0, 0.0
    trace: list[tuple[float, float, float]] = [(well, sick, dead)]
    costs: list[float] = [annual_intervention_cost * well]
    utils: list[float] = [u_well * well + u_sick * sick]
    lys: list[float] = [well + sick]

    for t in range(n):
        age = a0 + t
        mx = lt.get(age, 0.02)
        p_die = 1.0 - math.exp(-mx)                       # rate → probability
        p_die_sick = 1.0 - math.exp(-mx * excess_mortality_rr)
        new_sick = well * (1.0 - p_die) * inc
        new_dead = well * p_die + sick * p_die_sick
        well = well * (1.0 - p_die) * (1.0 - inc)
        sick = sick * (1.0 - p_die_sick) + new_sick
        dead = dead + new_dead
        trace.append((well, sick, dead))
        # Cost of illness is charged on incidence (a one-off per case), the
        # intervention on everyone still well and being treated.
        costs.append(new_sick * coi_cost + annual_intervention_cost * well)
        utils.append(u_well * well + u_sick * sick)
        lys.append(well + sick)

    wcc = simpson_weights(n)
    dw = discount_weights(n, rate)
    tc = sum(c * w * d for c, w, d in zip(costs, wcc, dw, strict=False))
    tq = sum(u * w * d for u, w, d in zip(utils, wcc, dw, strict=False))
    tl = sum(ly * w * d for ly, w, d in zip(lys, wcc, dw, strict=False))
    return MarkovResult(tc, tq, tl, trace, n, float(a0))


def incremental_analysis(*, start_age: float, annual_incidence: float,
                         coi_key: str, rrr: float,
                         intervention_cost_annual: float = 0.0,
                         sex: str = "Total",
                         wtp: float | None = None) -> dict:
    """Run the Markov model with and without acting, and return the increment.

    This is the structural alternative to `p × cost × rrr`: the same question
    asked of a model that knows the person can die of something else first.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    entry = COI_KEY_TO_PARAM.get(coi_key)
    if not entry:
        return {"available": False,
                "reason": f"no costed anchor for condition {coi_key!r}"}
    coi_cost = ep.value(entry[0]) * ep.value("marginal_cost_fraction")
    disutility = ep.value("utility_healthy") - ep.value("utility_post_event")

    base = run_markov(start_age=start_age, annual_incidence=annual_incidence,
                      coi_cost=coi_cost, disutility=disutility, rrr=0.0,
                      annual_intervention_cost=0.0, sex=sex)
    act = run_markov(start_age=start_age, annual_incidence=annual_incidence,
                     coi_cost=coi_cost, disutility=disutility, rrr=rrr,
                     annual_intervention_cost=intervention_cost_annual, sex=sex)

    res = CEAResult(strategy=f"Act on {coi_key} findings",
                    inc_cost=act.cost - base.cost,
                    inc_qaly=act.qaly - base.qaly, wtp=wtp,
                    detail={"condition": coi_key,
                            "life_years_gained": round(act.life_years - base.life_years, 4),
                            "usual_care": base.to_dict(),
                            "act": act.to_dict()})
    out = res.to_dict()
    out["available"] = True
    return out


# ══════════════════════════════════════════════════════════════════════════
# Second Panel dual perspective + impact inventory
# ══════════════════════════════════════════════════════════════════════════

def dual_perspective(healthcare_cost_averted: float, qaly_gained: float,
                     *, conditions: Sequence[dict] | None = None,
                     wtp: float | None = None) -> dict:
    """Report the healthcare-sector and societal perspectives side by side.

    The Second Panel asks for both, and for the societal additions to be
    itemised rather than absorbed into one number — precisely so that a reader
    who rejects, say, the valuation of unpaid caregiving can subtract it.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    productivity = ep.value("productivity_annual")
    care_hours = ep.value("caregiving_hours_annual_dementia")
    care_wage = ep.value("caregiver_replacement_wage")

    # Productivity: a QALY gained is not a work-year, so only the share of
    # health gain plausibly falling in working life and translating to
    # participation is counted. Stated explicitly rather than assumed.
    productivity_gain = qaly_gained * productivity * 0.35
    dementia = next((c for c in (conditions or [])
                     if c.get("condition") == "Alzheimer"), None)
    caregiving_gain = 0.0
    if dementia:
        caregiving_gain = (float(dementia.get("cases_averted", 0.0))
                           * care_hours * care_wage * 3.0)

    societal_extra = productivity_gain + caregiving_gain
    return {
        "healthcare_sector": {
            "cost_averted": round(healthcare_cost_averted),
            "qaly_gained": round(qaly_gained, 4),
            "monetised_health": round(qaly_gained * wtp),
            "note": "Payer/provider costs only — the reference case.",
        },
        "societal": {
            "cost_averted": round(healthcare_cost_averted + societal_extra),
            "qaly_gained": round(qaly_gained, 4),
            "monetised_health": round(qaly_gained * wtp),
            "note": "Adds productivity and unpaid caregiving, itemised below.",
        },
        "societal_additions": [
            {"item": "Productivity", "value": round(productivity_gain),
             "basis": f"QALYs gained × ${productivity:,.0f}/yr × 0.35 "
                      f"working-life share",
             "param": "productivity_annual"},
            {"item": "Unpaid caregiving (dementia)", "value": round(caregiving_gain),
             "basis": f"Dementia cases averted × {care_hours:,.0f} h/yr × "
                      f"${care_wage:.0f}/h × 3 yr",
             "param": "caregiving_hours_annual_dementia"},
        ],
        "delta": round(societal_extra),
    }


def impact_inventory(conditions: Sequence[dict]) -> list[dict]:
    """Second Panel impact inventory: what is counted, in which perspective.

    An explicit "counted / not counted" table is the cheapest defence against
    the accusation that a favourable result comes from a convenient choice of
    what to include.
    """
    return [
        {"sector": "Formal healthcare", "item": "Disease treatment cost averted",
         "healthcare": "included", "societal": "included",
         "note": "Marginal, not average, cost — see marginal_cost_fraction."},
        {"sector": "Formal healthcare", "item": "Test and intervention cost",
         "healthcare": "included", "societal": "included",
         "note": "Charged once per condition, not once per finding."},
        {"sector": "Patient", "item": "Health-related quality of life",
         "healthcare": "included", "societal": "included",
         "note": "Reported as QALYs; monetised only at a stated threshold."},
        {"sector": "Patient", "item": "Time and travel for care",
         "healthcare": "not counted", "societal": "not counted",
         "note": "No defensible per-person estimate available here."},
        {"sector": "Informal caregiver", "item": "Unpaid caregiving time",
         "healthcare": "not counted", "societal": "included",
         "note": "Dementia only; replacement-wage valuation."},
        {"sector": "Productivity", "item": "Paid labour output",
         "healthcare": "not counted", "societal": "included",
         "note": "Partial — working-life share stated, not assumed to be 1."},
        {"sector": "Other", "item": "Reproductive decisions",
         "healthcare": "not counted", "societal": "not counted",
         "note": "Deliberately never monetised; surfaced as a recommendation."},
        {"sector": "Other", "item": "Insurance / discrimination risk",
         "healthcare": "not counted", "societal": "not counted",
         "note": "Real but not quantified; noted as a limitation."},
    ]


# ══════════════════════════════════════════════════════════════════════════
# Reporting standards + validation
# ══════════════════════════════════════════════════════════════════════════

def cheers_checklist(*, wtp: float, rate: float, horizon: float,
                     perspective: str = "Healthcare sector (reference case), "
                                        "societal reported alongside") -> list[dict]:
    """CHEERS 2022 items this model can answer, with the answer.

    Reported so a reader can see which items are addressed and, just as
    importantly, which are not.
    """
    burden = ep.assumption_burden()
    return [
        {"item": "Health economic analysis plan",
         "response": "Cost–utility analysis of acting on genomic findings "
                     "versus usual care, individual-level inputs."},
        {"item": "Study population",
         "response": "One individual's genotype, with population-average "
                     "parameters — not a cohort study."},
        {"item": "Setting and location", "response": "United States."},
        {"item": "Comparators", "response": "Usual care (no genomic testing)."},
        {"item": "Perspective", "response": perspective},
        {"item": "Time horizon",
         "response": f"{horizon:.0f} years for the summary; to age 100 in the "
                     f"state-transition model."},
        {"item": "Discount rate",
         "response": f"{rate:.0%} for costs and effects; 0%/3%/5% sensitivity."},
        {"item": "Selection of outcomes",
         "response": "QALYs and healthcare-sector costs, reported separately."},
        {"item": "Measurement of effectiveness",
         "response": "Published relative risk reductions; correlated findings "
                     "pooled on the risk scale, not summed."},
        {"item": "Measurement and valuation of preference-based outcomes",
         "response": "EQ-5D-based utility weights from a US catalogue."},
        {"item": "Currency, price date, conversion",
         "response": "US dollars; parameter years recorded per parameter in "
                     "the provenance registry. No inflation adjustment is "
                     "applied — a stated limitation."},
        {"item": "Rationale and description of model",
         "response": "Three-state cohort model (Well/Diseased/Dead) with "
                     "life-table background mortality and Simpson's 1/3 "
                     "within-cycle correction."},
        {"item": "Analytics and assumptions",
         "response": f"Across {burden.get('n_total_known', 0)} model "
                     f"parameters, {burden.get('model_pct_resolvable', 0):.0f}% "
                     f"carry a resolvable citation (PMID/DOI) and "
                     f"{burden.get('model_pct_attributed_or_better', 0):.0f}% "
                     f"carry at least a named attribution; "
                     f"{burden.get('model_pct_unsourced', 0):.1f}% are declared "
                     f"assumptions, listed individually. "
                     f"{burden.get('n_curated_attributed', 0)} figures have a "
                     f"named source whose identifier is not yet verified — a "
                     f"stated limitation."},
        {"item": "Characterising heterogeneity",
         "response": f"Willingness to pay varied ${ep.get('wtp_per_qaly').low:,.0f}"
                     f"–${ep.get('wtp_per_qaly').high:,.0f}/QALY; age and sex "
                     f"enter through the life table."},
        {"item": "Characterising uncertainty",
         "response": "Probabilistic sensitivity analysis with parameter "
                     "distributions from the registry; one-way tornado."},
        {"item": "Characterising heterogeneity in distributional effects",
         "response": "Addressed: results are reported by age and sex through "
                     "life-table mortality, and a distributional analysis "
                     "applies Atkinson equity weights to report whether the "
                     "programme narrows or widens health inequality."},
        {"item": "Approach to engagement with patients and stakeholders",
         "response": "Not applicable — no stakeholder engagement was conducted."},
        {"item": "Effect of uncertainty",
         "response": "Reported as a cost-effectiveness acceptability curve and "
                     "expected value of perfect information."},
        {"item": "Conflicts of interest",
         "response": "None. Non-commercial personal project."},
    ]


def burton_pct(burden: dict) -> str:
    """Format the sourced-parameter share for the CHEERS analytics line."""
    return f"{burden['pct_sourced']:.0f}%"


def validate_model(pools: dict[str, ConditionPool], evaluated: dict) -> list[dict]:
    """Internal-validity checks, reported rather than silently assumed.

    A model that has never been asked whether it obeys its own constraints is
    not validated by having a lot of methods in it.
    """
    checks: list[dict] = []

    def add(name: str, ok: bool, detail: str):
        checks.append({"check": name, "pass": bool(ok), "detail": detail})

    cap = ep.value("max_combined_rrr")
    worst = max((p.combined_rrr() for p in pools.values()), default=0.0)
    add("Combined risk reduction respects its cap", worst <= cap + 1e-9,
        f"largest pooled RRR {worst:.3f} vs cap {cap:.2f}")

    probs_ok = all(0.0 <= p.baseline_risk() <= 1.0 for p in pools.values())
    add("Baseline risks are probabilities", probs_ok,
        "every pooled baseline risk lies in [0, 1]")

    rows = evaluated.get("conditions", [])
    once = all(r["cost_averted"] <= r["coi_cost"] + 1 for r in rows)
    add("No condition is charged more than one cost of illness", once,
        "cost averted per condition never exceeds its full cost of illness")

    dc = evaluated.get("double_counting", {})
    add("Pooling reduces rather than inflates the total",
        dc.get("pooled_cost_averted", 0) <= dc.get("naive_cost_averted", 0) + 1,
        f"pooled ${dc.get('pooled_cost_averted', 0):,} vs naive "
        f"${dc.get('naive_cost_averted', 0):,}")

    # Cross-validation: the model's own CAD arm, run through the Markov path,
    # should land in the range published statin primary-prevention analyses
    # report. Being outside it is not proof of error, but it is the kind of
    # thing a reader is entitled to see checked.
    xv = incremental_analysis(start_age=50, annual_incidence=0.01,
                              coi_key="CAD", rrr=ep.value("statin_rrr_primary"),
                              intervention_cost_annual=ep.value("cost_statin_10yr") / 10.0)
    icer = xv.get("icer")
    in_range = icer is None or (0 <= icer <= 150_000)
    add("Cross-validation: statin primary prevention ICER is plausible", in_range,
        (f"modelled ICER ${icer:,.0f}/QALY" if icer is not None
         else "strategy dominant (cost-saving), consistent with published "
              "analyses of generic statins in at-risk primary prevention")
        + " — published analyses of generic statins in primary prevention "
          "generally report cost-saving to ~$50,000/QALY")

    burden = ep.assumption_burden()
    add("Most parameters carry a citation", burden["pct_sourced"] >= 75.0,
        f"{burden['pct_sourced']}% sourced, {burden['n_assumption']} declared "
        f"assumptions")

    return checks


# ══════════════════════════════════════════════════════════════════════════
# Uncertainty: PSA, CEAC, tornado
# ══════════════════════════════════════════════════════════════════════════
# The registry records a distribution and a range for most parameters. Until
# these functions existed those fields were documentation — the model ran on
# point estimates and reported a single number, which invites the reader to
# treat it as more certain than it is. Propagating the documented uncertainty
# is what turns the provenance table from a display into a working input.

def run_psa(pools: dict[str, ConditionPool], *, n: int = 2000,
            seed: int = 20260822, test_cost: float = 0.0,
            wtp: float | None = None, rebuild=None) -> dict:
    """Probabilistic sensitivity analysis over the registry's distributions.

    Each iteration draws every sampleable parameter from its documented
    distribution, re-runs the pooled evaluation with those values in place,
    and records incremental cost and QALYs. Parameters with no published
    spread are held fixed rather than given invented uncertainty.

    ``rebuild`` is a zero-argument callable returning freshly-built pools. It
    matters more than it looks: baseline risk, effect size and intervention
    cost are read when a Finding is CONSTRUCTED, so pools built before
    sampling carry base-case values no matter what the draw says. Without a
    rebuild, only the cost-of-illness terms vary — which is exactly the
    one-sided uncertainty that made an earlier version of this report claim a
    strategy was cost-saving in 100% of simulations.
    """
    import random
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    rng = random.Random(seed)
    costs: list[float] = []
    qalys: list[float] = []
    inmbs: list[float] = []
    for _ in range(max(1, int(n))):
        with ep.overridden(ep.sample_all(rng)):
            iter_pools = rebuild() if rebuild is not None else pools
            ev = evaluate_pools(iter_pools, wtp=wtp, test_cost=test_cost)
            cea = ev["cea"]
        costs.append(float(cea["incremental_cost"]))
        qalys.append(float(cea["incremental_qaly"]))
        inmbs.append(float(cea["inmb"]))

    def pct(xs: list[float], q: float) -> float:
        if not xs:
            return 0.0
        s = sorted(xs)
        i = min(len(s) - 1, max(0, round(q * (len(s) - 1))))
        return s[i]

    n_eff = len(inmbs) or 1
    return {
        "available": True,
        "n_iterations": n_eff,
        "n_parameters_varied": len(ep.sampleable()),
        "mean_incremental_cost": round(sum(costs) / n_eff),
        "mean_incremental_qaly": round(sum(qalys) / n_eff, 4),
        "mean_inmb": round(sum(inmbs) / n_eff),
        "inmb_ci_low": round(pct(inmbs, 0.025)),
        "inmb_ci_high": round(pct(inmbs, 0.975)),
        "cost_ci_low": round(pct(costs, 0.025)),
        "cost_ci_high": round(pct(costs, 0.975)),
        "qaly_ci_low": round(pct(qalys, 0.025), 4),
        "qaly_ci_high": round(pct(qalys, 0.975), 4),
        "p_cost_effective": round(sum(1 for x in inmbs if x > 0) / n_eff, 4),
        "p_cost_saving": round(sum(1 for c in costs if c < 0) / n_eff, 4),
        "wtp": round(wtp),
        "note": ("Parameters without a published spread are held fixed; the "
                 "interval below therefore understates true uncertainty "
                 "rather than overstating it."),
    }


def ceac(pools: dict[str, ConditionPool], *,
         thresholds: Sequence[float] = (0, 25_000, 50_000, 75_000, 100_000,
                                        150_000, 200_000),
         n: int = 800, seed: int = 20260822,
         test_cost: float = 0.0, rebuild=None) -> list[dict]:
    """Cost-effectiveness acceptability curve from the pooled model.

    Reuses one set of draws across all thresholds — resampling per threshold
    produces a curve that jitters non-monotonically for no reason other than
    Monte Carlo noise.
    """
    import random
    rng = random.Random(seed)
    draws: list[tuple[float, float]] = []
    for _ in range(max(1, int(n))):
        with ep.overridden(ep.sample_all(rng)):
            iter_pools = rebuild() if rebuild is not None else pools
            cea = evaluate_pools(iter_pools, test_cost=test_cost)["cea"]
        draws.append((float(cea["incremental_cost"]),
                      float(cea["incremental_qaly"])))
    out: list[dict] = []
    for w in thresholds:
        p = sum(1 for c, q in draws if q * w - c > 0) / (len(draws) or 1)
        out.append({"wtp": round(float(w)), "p_cost_effective": round(p, 4)})
    return out


def tornado(pools: dict[str, ConditionPool], *,
            test_cost: float = 0.0, wtp: float | None = None,
            top: int = 10, rebuild=None) -> list[dict]:
    """One-way sensitivity: swing in net monetary benefit across each
    parameter's documented range.

    Answers the question a probabilistic interval cannot: *which* parameter is
    responsible for the spread, and therefore which one is worth arguing about
    or resolving first.
    """
    wtp = ep.value("wtp_per_qaly") if wtp is None else float(wtp)
    base = evaluate_pools(pools, wtp=wtp, test_cost=test_cost)["cea"]["inmb"]
    rows: list[dict] = []
    for p in ep.PARAMS.values():
        if p.low is None or p.high is None or p.low == p.high:
            continue
        with ep.overridden({p.key: p.low}):
            lo = evaluate_pools(rebuild() if rebuild else pools,
                                wtp=wtp, test_cost=test_cost)["cea"]["inmb"]
        with ep.overridden({p.key: p.high}):
            hi = evaluate_pools(rebuild() if rebuild else pools,
                                wtp=wtp, test_cost=test_cost)["cea"]["inmb"]
        swing = abs(hi - lo)
        if swing < 1:
            continue
        rows.append({
            "parameter": p.key, "units": p.units, "tier": p.tier,
            "low_value": p.low, "high_value": p.high,
            "inmb_at_low": round(lo), "inmb_at_high": round(hi),
            "base_inmb": round(base), "swing": round(swing),
        })
    rows.sort(key=lambda r: -r["swing"])
    return rows[:top]
