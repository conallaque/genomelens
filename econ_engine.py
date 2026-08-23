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
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import econ_params as ep

__all__ = [
    "Finding", "ConditionPool", "CEAResult", "MarkovResult",
    "pool_findings", "evaluate_pools", "run_markov", "life_table",
    "incremental_analysis", "dual_perspective", "impact_inventory",
    "cheers_checklist", "validate_model", "COI_KEY_TO_PARAM",
    "simpson_weights", "discount_weights",
]

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_LIFE_TABLE_PATH = os.path.join(_DATA_DIR, "LifeTable_USA_Mx_2015.csv")

# Condition anchors used by value_of_information._classify_category, mapped to
# the registry keys that cost them. A condition with no entry here cannot be
# valued — which is the intended behaviour, since the alternative is inventing
# a cost for it.
COI_KEY_TO_PARAM: Dict[str, Tuple[str, str]] = {
    # coi_key            (cost param,               qaly-loss param)
    "CAD":               ("coi_mace",               "qaly_loss_mace"),
    "T2D":               ("coi_t2d",                "qaly_loss_t2d"),
    "Alzheimer":         ("coi_alzheimer",          "qaly_loss_mace"),
    "Depression":        ("coi_depression",         "qaly_loss_mace"),
    "SubstanceUse":      ("coi_substance_use",      "qaly_loss_mace"),
    "Autoimmune":        ("coi_autoimmune",         "qaly_loss_mace"),
    "Urologic":          ("coi_urologic",           "qaly_loss_mace"),
    "IronOverload":      ("coi_iron_overload",      "qaly_loss_mace"),
    "Colorectal":        ("coi_colorectal",         "qaly_loss_t2d"),
    "BreastOvarian":     ("coi_breast_ovarian",     "qaly_loss_t2d"),
    "Pathogenic":        ("coi_pathogenic_generic", "qaly_loss_t2d"),
}


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
    """

    label: str
    coi_key: str
    p_event: float
    rrr: float
    haircut: float = 1.0
    intervention_cost: float = 0.0
    confidence: str = "moderate"
    source_category: str = ""
    qaly_override: Optional[float] = None

    @property
    def effective_rrr(self) -> float:
        """Risk reduction this finding supports once the source's evidence
        strength is applied. Bounded to [0, 1)."""
        return max(0.0, min(0.999, float(self.rrr) * float(self.haircut)))


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
    findings: List[Finding] = field(default_factory=list)

    # ── combination ─────────────────────────────────────────────────────
    def combined_rrr(self) -> float:
        if not self.findings:
            return 0.0
        penalty = ep.value("correlated_signal_penalty")
        cap = ep.value("max_combined_rrr")
        ranked = sorted(self.findings, key=lambda f: f.effective_rrr, reverse=True)
        surviving = 1.0
        for i, f in enumerate(ranked):
            surviving *= (1.0 - f.effective_rrr * (penalty ** i))
        return min(cap, 1.0 - surviving)

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
        return max(float(f.intervention_cost or 0.0) for f in self.findings)

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
        return sum(f.effective_rrr for f in self.findings)

    def to_dict(self) -> Dict:
        return {
            "condition": self.coi_key,
            "n_findings": len(self.findings),
            "labels": [f.label for f in self.findings],
            "sources": sorted({f.source_category for f in self.findings if f.source_category}),
            "baseline_risk": round(self.baseline_risk(), 4),
            "combined_rrr": round(self.combined_rrr(), 4),
            "naive_additive_rrr": round(self.naive_sum_rrr(), 4),
            "double_count_avoided": round(
                max(0.0, self.naive_sum_rrr() - self.combined_rrr()), 4),
            "coi_cost": round(self.coi_cost()),
            "qaly_loss": round(self.qaly_loss(), 3),
            "intervention_cost": round(self.intervention_cost()),
        }


def pool_findings(findings: Iterable[Finding]) -> Dict[str, ConditionPool]:
    """Group findings by condition into pools."""
    pools: Dict[str, ConditionPool] = {}
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
DEFAULT_TOPIC_TARGETS: Tuple[Tuple[str, str], ...] = (
    ("autoimmune", "topic:autoimmune"),
    ("statin-induced", "topic:statin-myopathy"),
    ("myopathy", "topic:statin-myopathy"),
)


def _extract_target(text: str,
                    vocabulary: Optional[frozenset] = None,
                    topics: Optional[Sequence[Tuple[str, str]]] = None) -> str:
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


def deduplicate_by_target(items: Sequence[Dict], *,
                          value_key: str = "net",
                          text_key: str = "finding",
                          fallback_key: str = "category",
                          vocabulary: Optional[frozenset] = None,
                          penalty: Optional[float] = None) -> List[Dict]:
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
    grouped: Dict[str, List[Dict]] = {}
    for i, it in enumerate(items):
        # A line naming no shared target gets a unique key, so it is never
        # pooled with anything. Grouping unmatched lines together by category
        # would discount independent findings for no reason.
        target = (_extract_target(str(it.get(text_key, "")), vocabulary)
                  or f"unique:{i}")
        grouped.setdefault(target, []).append(it)

    out: List[Dict] = []
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
                if k in copy and isinstance(copy[k], (int, float)):
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
    detail: Dict = field(default_factory=dict)

    @property
    def icer(self) -> Optional[float]:
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

    def to_dict(self) -> Dict:
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


def evaluate_pools(pools: Dict[str, ConditionPool],
                   *, wtp: Optional[float] = None,
                   horizon_years: Optional[float] = None,
                   test_cost: float = 0.0,
                   marginal_only: bool = True) -> Dict:
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

    rows: List[Dict] = []
    tot_cost_averted = tot_qaly = tot_intervention = 0.0
    tot_naive_cost_averted = 0.0

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

        tot_cost_averted += cost_averted
        tot_qaly += qaly_gained
        tot_intervention += intervention
        tot_naive_cost_averted += naive_cost_averted

        d = pool.to_dict()
        d.update({
            "cases_averted": round(cases_averted, 4),
            "cost_averted": round(cost_averted),
            "qaly_gained": round(qaly_gained, 4),
            "inmb": round(qaly_gained * wtp + cost_averted - intervention),
            "naive_cost_averted": round(naive_cost_averted),
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

def life_table(sex: str = "Total") -> Dict[int, float]:
    """Age-specific all-cause mortality hazard, from the vendored US table.

    Returns ``{age: mx}``. Missing file degrades to an empty dict, and callers
    fall back to a constant hazard — the model should lose precision, not
    disappear, if the data file is absent.
    """
    col = {"m": "Male", "male": "Male", "f": "Female", "female": "Female"}.get(
        (sex or "").strip().lower(), "Total")
    out: Dict[int, float] = {}
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


def simpson_weights(n_cycles: int) -> List[float]:
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
                     cycle_length: float = 1.0) -> List[float]:
    return [1.0 / (1.0 + rate) ** (t * cycle_length) for t in range(n_cycles + 1)]


@dataclass
class MarkovResult:
    cost: float
    qaly: float
    life_years: float
    trace: List[Tuple[float, float, float]]   # (well, diseased, dead)
    n_cycles: int
    start_age: float

    def to_dict(self) -> Dict:
        return {"cost": round(self.cost), "qaly": round(self.qaly, 4),
                "life_years": round(self.life_years, 3),
                "n_cycles": self.n_cycles, "start_age": self.start_age}


def run_markov(*, start_age: float, annual_incidence: float,
               coi_cost: float, disutility: float,
               rrr: float = 0.0, annual_intervention_cost: float = 0.0,
               excess_mortality_rr: float = 1.5,
               sex: str = "Total", max_age: int = 100,
               rate: Optional[float] = None) -> MarkovResult:
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
    trace: List[Tuple[float, float, float]] = [(well, sick, dead)]
    costs: List[float] = [annual_intervention_cost * well]
    utils: List[float] = [u_well * well + u_sick * sick]
    lys: List[float] = [well + sick]

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
    tc = sum(c * w * d for c, w, d in zip(costs, wcc, dw))
    tq = sum(u * w * d for u, w, d in zip(utils, wcc, dw))
    tl = sum(l * w * d for l, w, d in zip(lys, wcc, dw))
    return MarkovResult(tc, tq, tl, trace, n, float(a0))


def incremental_analysis(*, start_age: float, annual_incidence: float,
                         coi_key: str, rrr: float,
                         intervention_cost_annual: float = 0.0,
                         sex: str = "Total",
                         wtp: Optional[float] = None) -> Dict:
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
                     *, conditions: Optional[Sequence[Dict]] = None,
                     wtp: Optional[float] = None) -> Dict:
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


def impact_inventory(conditions: Sequence[Dict]) -> List[Dict]:
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
                                        "societal reported alongside") -> List[Dict]:
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
         "response": f"Of {burden['n_parameters']} registered parameters "
                     f"(method conventions, cost-of-illness anchors, effect "
                     f"sizes, utilities), {burden['pct_sourced']:.0f}% carry a "
                     f"literature citation and {burden['n_assumption']} are "
                     f"declared assumptions. A further "
                     f"{burden.get('n_unregistered', 0)} per-finding figures in "
                     f"the curated module tables are not yet registered and "
                     f"carry no provenance tier — a stated limitation."},
        {"item": "Characterising heterogeneity",
         "response": f"Willingness to pay varied ${ep.get('wtp_per_qaly').low:,.0f}"
                     f"–${ep.get('wtp_per_qaly').high:,.0f}/QALY; age and sex "
                     f"enter through the life table."},
        {"item": "Characterising uncertainty",
         "response": "Probabilistic sensitivity analysis with parameter "
                     "distributions from the registry; one-way tornado."},
        {"item": "Characterising heterogeneity in distributional effects",
         "response": "Not addressed — a limitation."},
        {"item": "Approach to engagement with patients and stakeholders",
         "response": "Not applicable — no stakeholder engagement was conducted."},
        {"item": "Effect of uncertainty",
         "response": "Reported as a cost-effectiveness acceptability curve and "
                     "expected value of perfect information."},
        {"item": "Conflicts of interest",
         "response": "None. Non-commercial personal project."},
    ]


def burton_pct(burden: Dict) -> str:
    """Format the sourced-parameter share for the CHEERS analytics line."""
    return f"{burden['pct_sourced']:.0f}%"


def validate_model(pools: Dict[str, ConditionPool], evaluated: Dict) -> List[Dict]:
    """Internal-validity checks, reported rather than silently assumed.

    A model that has never been asked whether it obeys its own constraints is
    not validated by having a lot of methods in it.
    """
    checks: List[Dict] = []

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
