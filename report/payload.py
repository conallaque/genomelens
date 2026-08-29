"""The canonical report payload — one normalised view of every economic result.

WHY THIS EXISTS. The economics engine returns nested dictionaries whose keys
carry no semantics: ``inmb`` appears in the deterministic reference case, in the
probabilistic analysis, in the Markov cross-check and in the per-finding tables,
and it means something different in each. A renderer reading ``x["inmb"]`` has
no way to know which question that number answers, and the report has printed
two of them beside each other under one heading.

This module does *not* recompute anything. Every field is lifted from a payload
the engine already returned. What it adds is **naming**: each quantity arrives
at the renderer with its identity attached, so the renderer never has to infer
which economic concept it is holding.

The rules the schema enforces by construction:

* **Deterministic and probabilistic results are different fields.**
  ``reference_case.nmb`` and ``uncertainty.psa_mean_nmb`` are both correct and
  are not equal. They cannot be swapped, and neither can be labelled "net
  monetary benefit" without a qualifier.
* **The Markov model is a cross-check, not the answer.** It lives in
  ``structural_crosscheck`` with its own field names, so it cannot be rendered
  as the genome-specific reference case.
* **Per-finding and pooled values are separate collections.** The reference case
  pools correlated findings before crediting value; per-finding figures do not
  sum to it, and the schema does not let a renderer pretend otherwise.
* **Every whole-genome-sequencing quantity declares its basis.** Either it
  counts findings this genome actually contains
  (``AnalysisBasis.OBSERVED_FINDINGS``) or it is a population-average
  expectation for someone who has not been tested yet
  (``AnalysisBasis.PROSPECTIVE_EXPECTED_YIELD``). Mixing them in one figure is
  the defect the testing-decision page exists to avoid.

Calculation stays where it is. This is a boundary, not a replacement.
"""
from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from econ.identity import is_legacy as _identity_is_legacy

__all__ = [
    "AdvancedAnalyses",
    "AnalysisBasis",
    "ConditionResult",
    "Corrections",
    "EconomicsReportPayload",
    "FindingEconomics",
    "PersonalView",
    "PricingPath",
    "Provenance",
    "ReferenceCase",
    "ReportMetadata",
    "StructuralCrossCheck",
    "TestingDecision",
    "Uncertainty",
    "build_report_payload",
    "payload_to_json",
]


class AnalysisBasis(str, Enum):
    """What kind of question a number answers. Never mix these in one figure."""

    OBSERVED_FINDINGS = "observed_findings"
    """Counted from findings this genome actually contains. Retrospective:
    the sequencing already happened and either found something or did not."""

    PROSPECTIVE_EXPECTED_YIELD = "prospective_expected_yield"
    """A population-average expectation for someone deciding whether to test.
    Says nothing about what any particular genome contains."""


class PricingPath(str, Enum):
    """Which code path priced a finding.

    The repository prices individual findings twice, through parameterisations
    that were never reconciled, and they disagree — see ``C3`` in
    ``docs/dev/econ-report-architecture.md``. Recording the path keeps the
    disagreement visible instead of letting whichever table renders last win.
    """

    CURATED_TABLE = "curated_per_finding_table"
    """``econ/health_economics.py`` — curated cost/value/QALY per category."""

    VOI_PARAMETRIC = "voi_parametric"
    """``econ/value_of_information.py`` — p_event x RRR x cost-of-illness."""

    CONDITION_POOLED = "condition_pooled"
    """``econ/engine.py`` — the reference case, after pooling."""


def _f(v, default: float = 0.0) -> float:
    try:
        if v is None:
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _i(v, default: int = 0) -> int:
    try:
        if v is None:
            return default
        return int(v)
    except (TypeError, ValueError):
        return default


def _s(v, default: str = "") -> str:
    return default if v is None else str(v)


# Vocabulary the economics layer understands. Upstream modules emit strings
# outside it — the novel-variants module reports "higher" for its most confident
# predictions — and `_CONFIDENCE_CONC` in value_of_information.py silently
# defaults those to the moderate concentration. Normalising here fixes the
# *display* and grouping only; the PSA spread is a calculation concern and is
# reported by the validator rather than changed here.
_CONFIDENCE_VOCABULARY = {"high", "moderate", "low"}
_CONFIDENCE_ALIASES = {
    "higher": "high", "very high": "high", "strong": "high",
    "medium": "moderate", "intermediate": "moderate",
    "lower": "low", "weak": "low", "very low": "low",
}


def normalise_confidence(v) -> str:
    raw = _s(v).strip().lower()
    if raw in _CONFIDENCE_VOCABULARY:
        return raw
    return _CONFIDENCE_ALIASES.get(raw, raw)


# ── metadata ──────────────────────────────────────────────────────────────────

@dataclass
class ReportMetadata:
    build_id: str = ""
    generated_at: str = ""
    input_type: str = ""
    input_label: str = ""
    is_synthetic: bool = False
    """Drives the synthetic-data banner. Read from metadata, never hard-coded
    into report prose."""
    analysis_horizon_years: float = 0.0
    perspective: str = "healthcare sector"
    currency: str = "USD"
    discount_rate: float = 0.0
    willingness_to_pay: float = 0.0


# ── findings ──────────────────────────────────────────────────────────────────

@dataclass
class FindingEconomics:
    """One finding as priced by one pricing path.

    The same genomic finding may appear more than once in a payload with
    different ``pricing_path`` values and different numbers. That is a faithful
    representation of the current engine, not an error in this schema.
    """

    finding_id: str = ""
    economic_pathway_id: str = ""
    """The join key: one finding acting on one condition through one
    intervention. Built from semantic components, not display text."""
    pathway_id_is_legacy: bool = True
    """True when the id was slugified from display text because the producing
    extractor supplies no semantic components yet. Such ids change if the
    wording changes; no canonical monetised record may rely on one."""
    action_id: str = ""
    condition_id: str = ""
    display_name: str = ""
    gene: str = ""

    @property
    def short_name(self) -> str:
        """A compact label that does not truncate mid-word.

        Finding names carry qualifiers in brackets and after dashes — "CYP2C19
        Intermediate Metabolizer (IM) (clopidogrel)". Chips and narrow columns
        cut those to "CYP2C19 Intermediate Metabolizer (IM) (clopi", which reads
        as a rendering fault. Cut at the first structural boundary instead.
        """
        s = (self.display_name or "").strip()
        for sep in (" \u2014 ", " - ", " ("):
            if sep in s:
                head = s.split(sep, 1)[0].strip()
                if len(head) >= 6:
                    s = head
                    break
        if len(s) > 34:
            s = s[:34].rsplit(" ", 1)[0] + "\u2026"
        return s
    variant: str = ""
    category: str = ""
    condition: str = ""
    action_summary: str = ""
    action_caveat: str = ""
    """A qualifier the pathway attaches to its action — e.g. that a predicted
    variant needs confirming before anything follows from it."""
    evidence_confidence: str = ""
    """Normalised to high / moderate / low for display and grouping."""
    raw_evidence_confidence: str = ""
    """Exactly what the producing module said, before normalisation. Retained so
    a vocabulary mismatch stays visible instead of being smoothed away."""
    provenance_tier: str = ""
    pricing_path: PricingPath = PricingPath.CURATED_TABLE

    expected_qaly_gain: float = 0.0
    medical_cost_averted: float = 0.0
    intervention_cost: float = 0.0
    net_cash: float = 0.0
    """medical_cost_averted - intervention_cost. Actual modelled money."""

    canonical_expected_nmb: float | None = None
    """THE ONE DOLLAR FIGURE THE MAIN REPORT SHOWS.

    Registry-backed parametric expected net monetary benefit: decomposable into
    explicit probabilities, effect sizes, costs and utilities; able to
    participate in probabilistic sensitivity analysis; responsive to
    willingness-to-pay, adherence and horizon. ``None`` when no parametric
    pathway exists for this finding yet — in which case the finding still
    renders, with "not yet standardised" in place of a value. A missing
    standardised figure must never hide a genomic result."""

    economic_value_basis: str = ""
    """``parametric_expected_nmb`` when canonical_expected_nmb is populated."""

    legacy_curated_value: float | None = None
    """AUDIT ONLY — never rendered beside the canonical figure.

    The per-finding total from the curated tables. Tracing it (see C3 in
    docs/dev/econ-report-architecture.md) showed it is not an alternative NMB
    for the same quantity: ``outcome_value`` is already an expected cost averted
    (p_rx x p_adr x rrr x adr_cost), and the curated path then multiplies it —
    and the QALY gain — by the *genotype prevalence*, the probability that a
    random person carries a variant this person is known to carry. So it is
    neither a gross event value nor an expected NMB but a mixture containing a
    conditioning error, which is why it may not carry a generic
    ``economic_value`` label."""

    legacy_curated_value_basis: str = ""
    """``curated_prevalence_weighted_mixture``."""

    health_economic_value: float = 0.0
    """Deprecated in the main report. Retained so existing renderers keep
    working during migration; the findings-first pages read
    ``canonical_expected_nmb``."""

    monetized_qaly_value: float = 0.0
    """expected_qaly_gain x willingness-to-pay. The health half of
    health_economic_value, separated so a renderer can keep cash and monetised
    health visually distinct."""

    event_probability: float | None = None
    effect_size: float | None = None
    adherence: float | None = None
    horizon_years: float | None = None
    pool_target: str = ""
    """Which condition pool this finding was routed to, when the pricing path
    records one. The only key that links a finding to ``ConditionResult``."""

    is_monetized: bool = True
    reason_not_monetized: str = ""
    is_hypothetical_or_awareness: bool = False
    is_wgs_only: bool = False
    sources: list[str] = field(default_factory=list)


@dataclass
class ConditionResult:
    condition: str = ""
    n_contributing_findings: int = 0
    contributing_findings: list[str] = field(default_factory=list)
    baseline_risk: float = 0.0
    naive_additive_rrr: float = 0.0
    """What summing the findings independently would have claimed."""
    pooled_efficacy_rrr: float = 0.0
    """After correlated-signal pooling, still at trial efficacy."""
    adherence_adjusted_rrr: float = 0.0
    """After real-world adherence. This is what the reference case uses."""
    adherence: float = 1.0
    adherence_archetype: str = ""
    cost_averted: float = 0.0
    qaly_gain: float = 0.0
    nmb: float = 0.0
    sources: list[str] = field(default_factory=list)


# ── the headline results ──────────────────────────────────────────────────────

@dataclass
class ReferenceCase:
    """The canonical genome-specific deterministic result. One per report."""

    strategy: str = ""
    incremental_cost: float = 0.0
    incremental_qalys: float = 0.0
    medical_cost_averted: float = 0.0
    intervention_cost: float = 0.0
    test_cost: float = 0.0
    nmb: float = 0.0
    """Deterministic point estimate. NOT the PSA mean — see
    ``Uncertainty.psa_mean_nmb``, which is a different quantity."""
    icer: float | None = None
    icer_note: str = ""
    dominance_status: str = ""
    wtp: float = 0.0
    horizon_years: float = 0.0
    discount_rate: float = 0.0
    marginal_cost_fraction: float = 0.0


@dataclass
class PersonalView:
    """The individual's own sheet. Answers a different question from the
    reference case: what the listed findings are worth one by one, before
    pooling. These totals do not reconcile to ``ReferenceCase`` and are not
    meant to."""

    medical_cost_avoided: float = 0.0
    intervention_cost: float = 0.0
    net_cash: float = 0.0
    qaly_gain: float = 0.0
    monetized_qaly_value: float = 0.0
    total_health_economic_value: float = 0.0
    horizon_years: float = 0.0
    mean_adherence: float = 0.0
    n_items: int = 0
    n_not_monetized: int = 0
    verdict: str = ""
    canonical_total_expected_nmb: float = 0.0
    """Sum of the canonical per-finding expected NMB. Reported alongside the
    legacy curated total rather than replacing it silently, because the two are
    different quantities: this one is registry-backed and standalone, that one
    is a prevalence-weighted mixture. NEITHER is the reference case — see
    ``ReferenceCase.nmb``, which pools overlapping findings first."""
    n_findings_standardised: int = 0
    n_findings_not_standardised: int = 0
    legacy_curated_total: float = 0.0
    """The pre-canonicalisation total, retained for audit."""

    value_per_dollar_of_testing: float | None = None
    """Formerly rendered as "ROI" and "26.7:1". It divides monetised health by
    test cost, so it is not a financial return. Any renderer using it must say
    so."""


@dataclass
class TestingDecision:
    """Chip vs whole-genome sequencing.

    Every field is tagged with the basis it was computed on. The two bases
    answer different questions and are never summed.
    """

    __test__ = False    # not a pytest class; the name only looks like one

    no_testing_cost: float = 0.0
    chip_cost: float = 0.0
    wgs_cost: float = 0.0
    incremental_chip_to_wgs_cost: float = 0.0

    # ── prospective: before testing, averaged over a population ──
    prospective_expected_yield: float = 0.0
    prospective_number_needed_to_sequence: int | None = None
    prospective_value_per_finding: float = 0.0
    prospective_gross_expected_value: float = 0.0
    prospective_net_expected_value: float = 0.0
    prospective_pgx_incremental_value: float = 0.0
    prospective_basis: AnalysisBasis = AnalysisBasis.PROSPECTIVE_EXPECTED_YIELD

    # ── observed: what this genome actually contains ──
    observed_wgs_only_findings: int = 0
    observed_wgs_only_value: float = 0.0
    observed_basis: AnalysisBasis = AnalysisBasis.OBSERVED_FINDINGS

    worth_it: bool = False
    strategies: list[dict[str, Any]] = field(default_factory=list)
    caveat: str = ""


@dataclass
class Uncertainty:
    psa_available: bool = False
    psa_iterations: int = 0
    n_parameters_varied: int = 0
    psa_mean_nmb: float = 0.0
    """Mean across simulations. Differs legitimately from
    ``ReferenceCase.nmb``; both are correct."""
    psa_mean_incremental_cost: float = 0.0
    psa_mean_incremental_qalys: float = 0.0
    nmb_ci_low: float = 0.0
    nmb_ci_high: float = 0.0
    probability_cost_effective: float = 0.0
    probability_cost_saving: float = 0.0
    ceac: list[dict[str, Any]] = field(default_factory=list)
    tornado: list[dict[str, Any]] = field(default_factory=list)
    evpi: dict[str, Any] = field(default_factory=dict)
    evppi: dict[str, Any] = field(default_factory=dict)
    note: str = ""


@dataclass
class Corrections:
    """The chain that turns an optimistic raw estimate into the reference case.
    Every step here is a documented reduction."""

    naive_cost_averted: float = 0.0
    pooled_cost_averted: float = 0.0
    duplicate_signal_removed: float = 0.0
    duplicate_signal_pct: float = 0.0
    efficacy_cost_averted: float = 0.0
    effectiveness_cost_averted: float = 0.0
    adherence_reduction: float = 0.0
    adherence_qaly_reduction: float = 0.0
    adherence_pct_of_benefit_lost: float = 0.0
    marginal_cost_fraction: float = 0.0
    midpoint_discount_factor: float = 0.0
    explanation: str = ""


@dataclass
class Provenance:
    registry_n_parameters: int = 0
    registry_n_published: int = 0
    registry_n_derived: int = 0
    registry_n_assumption: int = 0
    registry_pct_sourced: float = 0.0
    """Denominator is the registry only. Do not present beside
    ``model_pct_resolvable`` without labelling both denominators."""
    model_n_total_known: int = 0
    model_pct_resolvable: float = 0.0
    model_pct_attributed_or_better: float = 0.0
    model_pct_unsourced: float = 0.0
    unresolved_sources: list[dict[str, Any]] = field(default_factory=list)
    declared_assumptions: list[dict[str, Any]] = field(default_factory=list)
    scope_note: str = ""
    weighted_provenance: None = None
    """TODO(design): a sensitivity-weighted evidence score, so parameters that
    drive large NMB swings count for more than parameters that barely move the
    result. Deliberately not implemented here: any defensible weighting needs
    new registry parameters, and ``tests/unit/test_econ_params.py`` asserts
    ``pct_sourced >= 75.0`` with the live value at 75.4% — one unsourced
    addition fails CI. Raw provenance categories are retained instead."""


@dataclass
class StructuralCrossCheck:
    """The Markov cohort model.

    NOT the genome-specific reference case. It runs generic cohort parameters
    through a state-transition structure to check that the reference case's
    midpoint discounting is not badly wrong about timing. Its incremental QALYs
    are an order of magnitude larger because it models a different cohort over a
    different horizon. Belongs in the appendix.
    """

    available: bool = False
    is_reference_case: bool = False
    markov_incremental_cost: float = 0.0
    markov_qaly_gain: float = 0.0
    markov_incremental_life_years: float = 0.0
    markov_nmb: float = 0.0
    markov_icer: float | None = None
    markov_verdict: str = ""
    standard_care: dict[str, Any] = field(default_factory=dict)
    genomic_guided: dict[str, Any] = field(default_factory=dict)
    validation: list[dict[str, Any]] = field(default_factory=list)
    note: str = ""
    source: str = ""


@dataclass
class AdvancedAnalyses:
    budget_impact: dict[str, Any] = field(default_factory=dict)
    societal: dict[str, Any] = field(default_factory=dict)
    equity: dict[str, Any] = field(default_factory=dict)
    health_capital: dict[str, Any] = field(default_factory=dict)
    real_options: dict[str, Any] = field(default_factory=dict)
    risk_preferences: dict[str, Any] = field(default_factory=dict)
    behavioural: dict[str, Any] = field(default_factory=dict)
    privacy_welfare: dict[str, Any] = field(default_factory=dict)
    genomic_corrections: dict[str, Any] = field(default_factory=dict)
    cohort_projection: dict[str, Any] = field(default_factory=dict)


@dataclass
class EconomicsReportPayload:
    metadata: ReportMetadata = field(default_factory=ReportMetadata)
    reference_case: ReferenceCase = field(default_factory=ReferenceCase)
    personal_view: PersonalView = field(default_factory=PersonalView)
    findings: list[FindingEconomics] = field(default_factory=list)
    condition_results: list[ConditionResult] = field(default_factory=list)
    testing_decision: TestingDecision = field(default_factory=TestingDecision)
    uncertainty: Uncertainty = field(default_factory=Uncertainty)
    corrections: Corrections = field(default_factory=Corrections)
    provenance: Provenance = field(default_factory=Provenance)
    structural_crosscheck: StructuralCrossCheck = field(
        default_factory=StructuralCrossCheck)
    advanced: AdvancedAnalyses = field(default_factory=AdvancedAnalyses)
    plain_language: dict[str, Any] = field(default_factory=dict)
    """The engine's plain-language translation layer: number-needed-to-screen
    per condition, what would change the answer, the share of the result that
    rests on judgement. Carried through as-is; its own probability wording is
    NOT used, because it renders 0.9987 as "about 100 of every 100 runs"."""
    methods: list[str] = field(default_factory=list)
    """Method statements the engine emits for the appendix."""
    engine_validation: list[dict[str, Any]] = field(default_factory=list)
    """Checks the engine ran on itself, carried through unchanged."""
    report_validation: list[dict[str, Any]] = field(default_factory=list)
    """Filled by ``report.validate.validate_payload``."""

    # ── convenience views the renderer would otherwise have to derive ──

    def findings_by_path(self, path: PricingPath) -> list[FindingEconomics]:
        return [f for f in self.findings if f.pricing_path == path]

    @property
    def monetized_findings(self) -> list[FindingEconomics]:
        return [f for f in self.findings if f.is_monetized]

    @property
    def unmonetized_findings(self) -> list[FindingEconomics]:
        return [f for f in self.findings if not f.is_monetized]

    def findings_page_groups(self) -> list[tuple[str, list[FindingEconomics]]]:
        """Findings ordered for the findings-first page.

        CLINICAL SEMANTICS OUTRANK DOLLARS. Grouping comes first and sorting
        second, so a low-confidence wellness result cannot appear above a
        high-confidence prescribing finding on the strength of a modelled dollar
        estimate. Within a group, ordering is by canonical expected NMB
        descending; findings with no standardised value sort last but are never
        dropped, because a missing dollar figure must not hide a genomic result.
        """
        groups: dict[str, list[FindingEconomics]] = {
            "Medication & prescribing": [],
            "Risk & prevention": [],
            "Lower-confidence & exploratory": [],
            "Awareness — not monetised": [],
        }
        for f in self.findings:
            conf = (f.evidence_confidence or "").lower()
            cat = (f.category or "").lower()
            name = (f.display_name or "").lower()
            is_pgx = ("pharmacogenom" in cat or "pgx" in cat
                      or f.economic_pathway_id.startswith("pgx:")
                      or "metabolizer" in name or "function" in name)
            if not f.is_monetized:
                groups["Awareness — not monetised"].append(f)
            elif is_pgx and conf == "high":
                groups["Medication & prescribing"].append(f)
            elif conf in ("high", "moderate"):
                groups["Risk & prevention"].append(f)
            else:
                groups["Lower-confidence & exploratory"].append(f)

        def key(f: FindingEconomics) -> tuple[int, float]:
            # None sorts after every number, regardless of sign.
            if f.canonical_expected_nmb is None:
                return (1, 0.0)
            return (0, -f.canonical_expected_nmb)

        return [(name, sorted(items, key=key))
                for name, items in groups.items() if items]


# ── builder ───────────────────────────────────────────────────────────────────

def build_report_payload(
    economics_result: dict | None = None,
    voi_result: dict | None = None,
    personal_econ: dict | None = None,
    frontier_result: dict | None = None,
    *,
    metadata: dict | None = None,
) -> EconomicsReportPayload:
    """Normalise the engine's outputs into one typed payload.

    Reads only. Recomputes nothing. Where the engine reports a quantity, that
    quantity is carried through unchanged; where the engine reports two
    quantities that a renderer might confuse, they land in differently named
    fields.
    """
    economics_result = economics_result or {}
    voi_result = voi_result or {}
    personal_econ = personal_econ or {}
    frontier_result = frontier_result or {}
    meta_in = metadata or {}

    pooled = (voi_result.get("pooled_economics") or {})
    cea = pooled.get("cea") or {}
    dc = pooled.get("double_counting") or {}
    adh = pooled.get("adherence") or {}
    prov = pooled.get("provenance") or {}
    psa = pooled.get("psa") or {}
    wgs = pooled.get("wgs_decision") or {}
    markov = voi_result.get("markov") or {}

    wtp = _f(cea.get("wtp") or voi_result.get("wtp_base"), 100_000.0)

    # build_stamp() returns a dict; accept either it or a plain string so the
    # caller does not have to know which.
    stamp = meta_in.get("build_id")
    if isinstance(stamp, dict):
        build_id = _s(stamp.get("marker") or stamp.get("commit"))
        stamped_at = _s(stamp.get("generated_at"))
    else:
        build_id, stamped_at = _s(stamp), ""

    md = ReportMetadata(
        build_id=build_id,
        generated_at=_s(meta_in.get("generated_at")) or stamped_at,
        input_type=_s(voi_result.get("input_type") or meta_in.get("input_type")),
        input_label=_s(meta_in.get("input_label")),
        is_synthetic=bool(meta_in.get("is_synthetic", False)),
        analysis_horizon_years=_f(cea.get("horizon_years")),
        perspective=_s(meta_in.get("perspective"), "healthcare sector"),
        discount_rate=_f(cea.get("discount_rate") or voi_result.get("discount_rate")),
        willingness_to_pay=wtp,
    )

    ref = ReferenceCase(
        strategy=_s(cea.get("strategy")),
        incremental_cost=_f(cea.get("incremental_cost")),
        incremental_qalys=_f(cea.get("incremental_qaly")),
        medical_cost_averted=_f(cea.get("cost_averted")),
        intervention_cost=_f(cea.get("intervention_cost")),
        test_cost=_f(cea.get("test_cost")),
        nmb=_f(cea.get("inmb")),
        icer=None if cea.get("icer") is None else _f(cea.get("icer")),
        icer_note=_s(cea.get("icer_note")),
        dominance_status=_s(cea.get("verdict")),
        wtp=wtp,
        horizon_years=_f(cea.get("horizon_years")),
        discount_rate=_f(cea.get("discount_rate")),
        marginal_cost_fraction=_f(cea.get("marginal_cost_fraction")),
    )

    pv = PersonalView(
        medical_cost_avoided=_f(personal_econ.get("total_avoided")),
        intervention_cost=_f(personal_econ.get("total_intervention")),
        net_cash=_f(personal_econ.get("net_cash")),
        qaly_gain=_f(personal_econ.get("total_qaly")),
        monetized_qaly_value=_f(personal_econ.get("total_qaly_value")),
        total_health_economic_value=_f(personal_econ.get("total_net")),
        horizon_years=_f(personal_econ.get("horizon_years")),
        mean_adherence=_f(personal_econ.get("mean_adherence")),
        n_items=_i(personal_econ.get("n_items")),
        n_not_monetized=_i(personal_econ.get("n_not_monetised")),
        verdict=_s(personal_econ.get("verdict")),
        legacy_curated_total=_f(personal_econ.get("total_net")),
        value_per_dollar_of_testing=(
            None if personal_econ.get("value_to_cost_ratio") is None
            else _f(personal_econ.get("value_to_cost_ratio"))),
    )

    # ── findings: ONE record per economic pathway ────────────────────────────
    #
    # The report used to carry two dollar figures for the same finding —
    # $4,452 and $490 for CYP2C19 — because two pricing paths built separate
    # lists and both rendered. Merging on `economic_pathway_id` makes that
    # impossible to reproduce: a pathway priced by both paths becomes one
    # record whose canonical value is the parametric one and whose curated
    # figure is demoted to an audit field.
    #
    # Registry-backed parametric expected NMB is canonical because it
    # decomposes into explicit probabilities, effect sizes, costs and
    # utilities; participates in the probabilistic sensitivity analysis; and
    # responds correctly to willingness-to-pay, adherence and horizon. The
    # curated figure does none of those things and additionally conditions on
    # genotype prevalence — the chance a random person carries a variant this
    # person is known to carry.

    parametric: dict[str, dict] = {}
    parametric_order: list[str] = []
    for r in (voi_result.get("nmb_rows") or []):
        pid = _s(r.get("economic_pathway_id")) or f"unkeyed:{_s(r.get('label'))}"
        parametric[pid] = r
        parametric_order.append(pid)

    curated: dict[str, list[dict]] = {}
    for it in (personal_econ.get("items") or []):
        pid = _s(it.get("economic_pathway_id"))
        if pid:
            curated.setdefault(pid, []).append(it)

    findings: list[FindingEconomics] = []
    seen: set[str] = set()

    def _curated_totals(items: list[dict]) -> tuple[float, float, float, float]:
        return (sum(_f(i.get("avoided")) for i in items),
                sum(_f(i.get("intervention")) for i in items),
                sum(_f(i.get("qaly")) for i in items),
                sum(_f(i.get("net")) for i in items))

    # 1. Pathways the parametric model prices. Canonical.
    for pid in parametric_order:
        r = parametric[pid]
        seen.add(pid)
        cur = curated.get(pid) or []
        c_av, c_iv, _c_q, c_net = _curated_totals(cur)
        nmb = _f(r.get("nmb"))
        findings.append(FindingEconomics(
            finding_id=_s(r.get("label")),
            economic_pathway_id=pid,
            pathway_id_is_legacy=_identity_is_legacy(pid),
            display_name=_s(r.get("label")),
            evidence_confidence=normalise_confidence(r.get("confidence")),
            raw_evidence_confidence=_s(r.get("confidence")),
            pricing_path=PricingPath.VOI_PARAMETRIC,
            expected_qaly_gain=_f(r.get("dqaly")),
            medical_cost_averted=_f(r.get("dcost_averted")),
            intervention_cost=c_iv,
            net_cash=_f(r.get("dcost_averted")) - c_iv,
            canonical_expected_nmb=nmb,
            economic_value_basis="parametric_expected_nmb",
            legacy_curated_value=(c_net if cur else None),
            legacy_curated_value_basis=(
                "curated_prevalence_weighted_mixture" if cur else ""),
            health_economic_value=nmb,
            is_wgs_only=bool(r.get("wgs_only", False)),
            # The curated record's display name IS the record's
            # `clinical_benefit` — "Avoid clopidogrel non-response / stent
            # thrombosis (MACE)" — which is the decision this finding bears on.
            # Its `basis` describes the *method* ("Avoided adverse event x
            # probability of relevant exposure"), which is not what belongs in
            # an action column.
            # Prefer the curated record's clinical benefit; fall back to the
            # action the parametric pathway itself declares. The renderer never
            # supplies one.
            action_summary=(_s(cur[0].get("finding")) if cur
                            else _s(r.get("action"))),
            action_caveat=_s(r.get("action_caveat")),
            category=_s(cur[0].get("category") if cur else ""),
            adherence=(_f(cur[0].get("adherence"))
                       if cur and cur[0].get("adherence") is not None else None),
        ))

    # 2. Curated-only pathways. No standardised value exists, so none is
    #    invented — but the genomic finding still renders.
    # MIGRATION FALLBACK. Some curated items are synthesised inside
    # analyze_personal_economics from module results that never passed through
    # _econ_record, so they carry no pathway id — the carrier-derived "PTPN22
    # R620W carrier" row is one. Without a fallback the page renders PTPN22
    # twice: once canonical, once unstandardised. Matching on display name is
    # display-dependent and therefore explicitly a migration measure, not an
    # identity scheme; it folds the orphan in as a legacy value rather than
    # emitting a second row for one finding.
    by_name = {f.display_name.strip().lower(): f for f in findings
               if f.display_name}
    for it in (personal_econ.get("items") or []):
        pid = _s(it.get("economic_pathway_id"))
        if pid and pid in seen:
            continue
        name_key = _s(it.get("finding")).strip().lower()
        twin = by_name.get(name_key)
        if twin is not None and twin.canonical_expected_nmb is not None:
            # ACCUMULATE, do not overwrite or skip. PTPN22 carries two curated
            # records — the carrier row and the symptom-awareness row — and the
            # pathway's legacy total is both of them. Setting only the first
            # reported $174 against a pathway the curated sheet valued at
            # $5,800, which understates exactly the gap this table exists to
            # show.
            twin.legacy_curated_value = (
                (twin.legacy_curated_value or 0.0) + _f(it.get("net")))
            twin.legacy_curated_value_basis = (
                "curated_prevalence_weighted_mixture")
            continue
        key = pid or f"unkeyed:{_s(it.get('finding'))}"
        if key in seen:
            continue
        seen.add(key)
        c_av, c_iv = _f(it.get("avoided")), _f(it.get("intervention"))
        findings.append(FindingEconomics(
            finding_id=_s(it.get("finding")),
            economic_pathway_id=key,
            pathway_id_is_legacy=bool(it.get("pathway_id_is_legacy", True)),
            display_name=_s(it.get("finding")),
            category=_s(it.get("category")),
            action_summary=_s(it.get("basis")),
            evidence_confidence=normalise_confidence(it.get("confidence")),
            raw_evidence_confidence=_s(it.get("confidence")),
            pricing_path=PricingPath.CURATED_TABLE,
            expected_qaly_gain=_f(it.get("qaly")),
            medical_cost_averted=c_av,
            intervention_cost=c_iv,
            net_cash=c_av - c_iv,
            monetized_qaly_value=_f(it.get("qaly_value")),
            canonical_expected_nmb=None,
            economic_value_basis="not_yet_standardised",
            legacy_curated_value=_f(it.get("net")),
            legacy_curated_value_basis="curated_prevalence_weighted_mixture",
            health_economic_value=_f(it.get("net")),
            adherence=(None if it.get("adherence") is None
                       else _f(it.get("adherence"))),
            pool_target=_s(it.get("pool_target")),
            is_monetized=True,
        ))

    # 3. Findings deliberately carrying no dollar figure.
    for nm in (personal_econ.get("not_monetised") or []):
        findings.append(FindingEconomics(
            finding_id=_s(nm.get("id") or nm.get("finding")),
            display_name=_s(nm.get("finding") or nm.get("label")),
            category=_s(nm.get("category")),
            action_summary=_s(nm.get("action")),
            evidence_confidence=_s(nm.get("confidence")),
            pricing_path=PricingPath.CURATED_TABLE,
            canonical_expected_nmb=None,
            economic_value_basis="not_monetised",
            is_monetized=False,
            reason_not_monetized=_s(nm.get("reason")),
            is_hypothetical_or_awareness=True,
        ))

    # Canonical personal totals, derived from the merged findings rather than
    # from the curated sheet. Reported beside the legacy total, not instead of
    # it, because they are different quantities.
    _std = [f for f in findings if f.canonical_expected_nmb is not None]
    _unstd = [f for f in findings
              if f.canonical_expected_nmb is None and f.is_monetized]
    pv.canonical_total_expected_nmb = round(
        sum(f.canonical_expected_nmb or 0.0 for f in _std), 2)
    pv.n_findings_standardised = len(_std)
    pv.n_findings_not_standardised = len(_unstd)

    conditions = [
        ConditionResult(
            condition=_s(c.get("condition")),
            n_contributing_findings=_i(c.get("n_findings")),
            contributing_findings=list(c.get("findings") or []),
            baseline_risk=_f(c.get("baseline_risk")),
            naive_additive_rrr=_f(c.get("naive_additive_rrr")),
            pooled_efficacy_rrr=_f(c.get("pooled_efficacy_rrr",
                                          c.get("combined_rrr"))),
            adherence_adjusted_rrr=_f(c.get("combined_rrr")),
            adherence=_f(c.get("adherence"), 1.0),
            adherence_archetype=_s(c.get("adherence_archetype")),
            cost_averted=_f(c.get("cost_averted")),
            qaly_gain=_f(c.get("qaly_gained")),
            nmb=_f(c.get("inmb")),
            sources=list(c.get("sources") or []),
        )
        for c in (pooled.get("conditions") or [])
    ]

    # Comparator prices come from the strategy frontier. `voi_result["price"]`
    # is a different thing entirely — it prices GenomeLens's own panels against
    # buying the equivalent tests a la carte — so reading chip/WGS cost from it
    # silently produced $0.
    strategies = list(frontier_result.get("all_strategies")
                      or frontier_result.get("strategies") or [])

    def _strategy_cost(*names: str) -> float:
        for s in strategies:
            label = str(s.get("name") or s.get("strategy") or "").lower()
            if any(n in label for n in names):
                return _f(s.get("cost"))
        return 0.0

    td = TestingDecision(
        no_testing_cost=_strategy_cost("no test", "no-test", "baseline"),
        chip_cost=_strategy_cost("chip", "array", "genotyp"),
        wgs_cost=_strategy_cost("whole-genome", "whole genome", "wgs", "sequenc"),
        strategies=strategies,
        incremental_chip_to_wgs_cost=_f(wgs.get("incremental_cost")),
        prospective_expected_yield=_f(wgs.get("expected_incremental_yield")),
        prospective_number_needed_to_sequence=(
            None if wgs.get("number_needed_to_sequence") is None
            else _i(wgs.get("number_needed_to_sequence"))),
        prospective_value_per_finding=_f(wgs.get("value_per_finding")),
        prospective_gross_expected_value=_f(wgs.get("gross_expected_value")),
        prospective_net_expected_value=_f(wgs.get("net_expected_value")),
        prospective_pgx_incremental_value=_f(wgs.get("pgx_incremental_value")),
        observed_wgs_only_findings=_i(wgs.get("n_wgs_only_findings")),
        observed_wgs_only_value=_f(wgs.get("retrospective_value")),
        worth_it=bool(wgs.get("worth_it", False)),
        caveat=_s(wgs.get("caveat")),
    )

    unc = Uncertainty(
        psa_available=bool(psa.get("available", False)),
        psa_iterations=_i(psa.get("n_iterations")),
        n_parameters_varied=_i(psa.get("n_parameters_varied")),
        psa_mean_nmb=_f(psa.get("mean_inmb")),
        psa_mean_incremental_cost=_f(psa.get("mean_incremental_cost")),
        psa_mean_incremental_qalys=_f(psa.get("mean_incremental_qaly")),
        nmb_ci_low=_f(psa.get("inmb_ci_low")),
        nmb_ci_high=_f(psa.get("inmb_ci_high")),
        probability_cost_effective=_f(psa.get("p_cost_effective")),
        probability_cost_saving=_f(psa.get("p_cost_saving")),
        ceac=list(pooled.get("ceac") or voi_result.get("ceac") or []),
        tornado=list(pooled.get("tornado") or voi_result.get("tornado") or []),
        evpi=dict(voi_result.get("evpi") or {}),
        evppi=dict(voi_result.get("evppi") or {}),
        note=_s(psa.get("note")),
    )

    corr = Corrections(
        naive_cost_averted=_f(dc.get("naive_cost_averted")),
        pooled_cost_averted=_f(dc.get("pooled_cost_averted")),
        duplicate_signal_removed=_f(dc.get("inflation_removed")),
        duplicate_signal_pct=_f(dc.get("pct_removed")),
        efficacy_cost_averted=_f(adh.get("efficacy_cost_averted")),
        effectiveness_cost_averted=_f(adh.get("effectiveness_cost_averted")),
        adherence_reduction=_f(adh.get("value_lost_to_non_adherence")),
        adherence_qaly_reduction=_f(adh.get("qaly_lost_to_non_adherence")),
        adherence_pct_of_benefit_lost=_f(adh.get("pct_of_benefit_lost")),
        marginal_cost_fraction=_f(cea.get("marginal_cost_fraction")),
        midpoint_discount_factor=_f(cea.get("midpoint_discount_factor")),
        explanation=_s(dc.get("explanation")),
    )

    pr = Provenance(
        registry_n_parameters=_i(prov.get("n_parameters")),
        registry_n_published=_i(prov.get("n_published")),
        registry_n_derived=_i(prov.get("n_derived")),
        registry_n_assumption=_i(prov.get("n_assumption")),
        registry_pct_sourced=_f(prov.get("pct_sourced")),
        model_n_total_known=_i(prov.get("n_total_known")),
        model_pct_resolvable=_f(prov.get("model_pct_resolvable")),
        model_pct_attributed_or_better=_f(prov.get("model_pct_attributed_or_better")),
        model_pct_unsourced=_f(prov.get("model_pct_unsourced")),
        unresolved_sources=list(prov.get("unresolved_sources") or []),
        declared_assumptions=list(pooled.get("declared_assumptions") or []),
        scope_note=_s(prov.get("scope")),
    )

    sx = StructuralCrossCheck(
        available=bool(markov.get("available", False)),
        is_reference_case=False,
        markov_incremental_cost=_f(markov.get("incremental_cost")),
        markov_qaly_gain=_f(markov.get("incremental_qaly")),
        markov_incremental_life_years=_f(markov.get("incremental_life_years")),
        markov_nmb=_f(markov.get("nmb_at_wtp")),
        markov_icer=(None if markov.get("icer") is None
                     else _f(markov.get("icer"))),
        markov_verdict=_s(markov.get("verdict")),
        standard_care=dict(markov.get("standard_care") or {}),
        genomic_guided=dict(markov.get("genomic_guided") or {}),
        validation=list(markov.get("validation") or []),
        note=_s(markov.get("note")),
        source=_s(markov.get("src")),
    )

    adv = AdvancedAnalyses(
        budget_impact=dict(voi_result.get("budget_impact") or {}),
        societal=dict(pooled.get("dual_perspective") or {}),
        equity=dict(voi_result.get("longevity") or {}),
        health_capital=dict(voi_result.get("health_capital") or {}),
        real_options=dict(voi_result.get("real_option") or {}),
        risk_preferences=dict(voi_result.get("utility") or {}),
        behavioural=dict(voi_result.get("behavioural") or {}),
        privacy_welfare=dict(voi_result.get("information_economics") or {}),
        genomic_corrections=dict(voi_result.get("genomic_corrections") or {}),
        cohort_projection=dict(economics_result.get("payer_impact") or {}),
    )

    return EconomicsReportPayload(
        metadata=md, reference_case=ref, personal_view=pv, findings=findings,
        condition_results=conditions, testing_decision=td, uncertainty=unc,
        corrections=corr, provenance=pr, structural_crosscheck=sx, advanced=adv,
        plain_language=dict(pooled.get("plain") or {}),
        methods=list(voi_result.get("methods") or []),
        engine_validation=list(pooled.get("validation") or []),
    )


# ── serialisation ─────────────────────────────────────────────────────────────

def _default(o: Any) -> Any:
    if isinstance(o, Enum):
        return o.value
    if dataclasses.is_dataclass(o) and not isinstance(o, type):
        return dataclasses.asdict(o)
    return str(o)


def payload_to_json(payload: EconomicsReportPayload, *, indent: int = 2) -> str:
    """Serialise the payload. Enums become their string values, so the JSON is
    readable without the schema in hand."""
    return json.dumps(dataclasses.asdict(payload), indent=indent, default=_default)
