"""Consistency checks that run before a report is rendered.

WHAT THIS IS FOR. A report can be arithmetically correct in every module and
still publish a contradiction, because contradictions live between modules. The
engine already validates itself — pooling caps, probability ranges, cost-of-
illness ceilings — and those results are carried through as
``payload.engine_validation``. This layer checks the things only the *report*
can get wrong: identities that must hold across sections, quantities that must
not be conflated, and prose that must agree with the numbers beside it.

THREE SEVERITIES, AND THE DISTINCTION MATTERS.

``ERROR``
    The report would publish something false: a broken arithmetic identity, an
    impossible dominance classification, a probability outside [0, 1], a
    budget-impact peak that is not the maximum, or observed and prospective
    sequencing value contaminating one another. Callers should refuse to render.

``WARNING``
    The report would publish something true but weakly supported, or two results
    that diverge enough that a reader deserves to be told. Render, but say so.

``INFO``
    Two numbers differ and both are correct. The deterministic reference case
    and the probabilistic mean are not equal and never will be; per-finding
    values do not sum to the pooled total by design. These are recorded so a
    renderer can *explain* them, not suppress them.

A valid but conceptually different number is not an error. Most of the work here
is refusing to treat one as though it were.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from report.payload import (
    AnalysisBasis,
    EconomicsReportPayload,
    PricingPath,
)

__all__ = ["Finding", "Severity", "errors_in", "format_report", "pdf_is_blocked",
           "validate_payload"]


class Severity:
    ERROR = "ERROR"
    WARNING = "WARNING"
    INFO = "INFO"


@dataclass
class Finding:
    severity: str
    check: str
    detail: str
    field: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {"severity": self.severity, "check": self.check,
                "detail": self.detail, "field": self.field}


def _money_tolerance(magnitude: float) -> float:
    """How far an identity may drift before it counts as broken.

    The engine rounds incremental QALYs to four decimals and costs to whole
    dollars before reporting them. At a $100,000 threshold a four-decimal QALY
    carries up to $5 of rounding on its own, so a strict equality test would fire
    on correct arithmetic. The floor absorbs that; the proportional term keeps
    the check meaningful as the numbers grow.
    """
    return max(10.0, abs(magnitude) * 0.005)


def validate_payload(p: EconomicsReportPayload) -> list[dict[str, Any]]:
    """Run every report-level check. Returns findings, most severe first."""
    out: list[Finding] = []
    out += _check_nmb_identity(p)
    out += _check_net_cash_identity(p)
    out += _check_dominance(p)
    out += _check_probabilities(p)
    out += _check_pooling(p)
    out += _check_wgs_separation(p)
    out += _check_budget_impact(p)
    out += _check_structural_crosscheck(p)
    out += _check_pricing_path_divergence(p)
    out += _check_provenance(p)
    out += _check_monetised_flag_matches_the_value(p)
    out += _check_withheld_findings_carry_no_value(p)
    out += _check_path_reconciliation(p)
    out += _check_one_gene_one_finding(p)
    out += _check_render_identity(p)
    out += _note_legitimate_differences(p)

    order = {Severity.ERROR: 0, Severity.WARNING: 1, Severity.INFO: 2}
    out.sort(key=lambda f: order.get(f.severity, 3))
    return [f.as_dict() for f in out]


def errors_in(findings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [f for f in findings if f.get("severity") == Severity.ERROR]


def pdf_is_blocked(findings: list[dict[str, Any]]) -> bool:
    """Whether these validation findings must stop the PDF from being written.

    The policy — any ERROR blocks, warnings never do — used to live as an
    inline comprehension inside the pipeline, which meant the README could
    claim "the renderer blocks on broken arithmetic identities" with nothing
    asserting it. The behaviour was real; it was simply unreachable from a
    test without running the whole pipeline. Naming it here gives the claim
    something to point at, and `test_pdf_blocked_on_arithmetic_violation`
    exercises the same predicate the pipeline branches on.
    """
    return bool(errors_in(findings))


# ── identities ────────────────────────────────────────────────────────────────

def _check_nmb_identity(p: EconomicsReportPayload) -> list[Finding]:
    r = p.reference_case
    if not r.wtp:
        return [Finding(Severity.WARNING, "NMB identity",
                        "no willingness-to-pay threshold on the reference case; "
                        "identity not checkable", "reference_case.wtp")]
    expected = r.wtp * r.incremental_qalys - r.incremental_cost
    drift = abs(expected - r.nmb)
    tol = _money_tolerance(expected)
    if drift > tol:
        return [Finding(
            Severity.ERROR, "NMB identity",
            f"reference_case.nmb is ${r.nmb:,.0f} but "
            f"lambda*dQALY - dCost = ${expected:,.0f} "
            f"({r.wtp:,.0f} x {r.incremental_qalys:.4f} - "
            f"{r.incremental_cost:,.0f}); drift ${drift:,.2f} exceeds "
            f"tolerance ${tol:,.2f}", "reference_case.nmb")]
    return []


def _check_net_cash_identity(p: EconomicsReportPayload) -> list[Finding]:
    out: list[Finding] = []
    pv = p.personal_view
    expected = pv.medical_cost_avoided - pv.intervention_cost
    if abs(expected - pv.net_cash) > _money_tolerance(expected):
        out.append(Finding(
            Severity.ERROR, "Net-cash identity",
            f"personal_view.net_cash is ${pv.net_cash:,.0f} but "
            f"cost avoided - intervention = ${expected:,.0f}",
            "personal_view.net_cash"))

    for f in p.findings:
        if not f.is_monetized or f.pricing_path is not PricingPath.CURATED_TABLE:
            continue
        exp = f.medical_cost_averted - f.intervention_cost
        if abs(exp - f.net_cash) > _money_tolerance(exp):
            out.append(Finding(
                Severity.ERROR, "Net-cash identity (per finding)",
                f"{f.display_name}: net_cash ${f.net_cash:,.0f} but "
                f"averted - intervention = ${exp:,.0f}",
                f"findings[{f.finding_id}].net_cash"))
    return out


def _check_dominance(p: EconomicsReportPayload) -> list[Finding]:
    r = p.reference_case
    status = (r.dominance_status or "").lower()
    if "dominant" in status and "dominated" not in status:
        bad = []
        if r.incremental_cost >= 0:
            bad.append(f"incremental cost is ${r.incremental_cost:,.0f}, "
                       f"which is not negative")
        if r.incremental_qalys <= 0:
            bad.append(f"incremental QALYs are {r.incremental_qalys:.4f}, "
                       f"which is not positive")
        if bad:
            return [Finding(
                Severity.ERROR, "Dominance classification",
                f"reference case is labelled '{r.dominance_status}' but "
                + "; ".join(bad), "reference_case.dominance_status")]
    if ("dominated" in status and "not dominated" not in status
            and r.incremental_cost <= 0 and r.incremental_qalys >= 0):
        return [Finding(
            Severity.ERROR, "Dominance classification",
            f"reference case is labelled '{r.dominance_status}' but costs "
            f"less (${r.incremental_cost:,.0f}) and gains health "
            f"({r.incremental_qalys:.4f} QALYs)",
            "reference_case.dominance_status")]
    if r.icer is None and not (r.icer_note or "").strip():
        return [Finding(
            Severity.WARNING, "ICER note",
            "ICER is undefined and carries no explanatory note; a blank ICER "
            "reads as a missing calculation rather than a dominant strategy",
            "reference_case.icer_note")]
    return []


def _check_probabilities(p: EconomicsReportPayload) -> list[Finding]:
    out: list[Finding] = []
    u = p.uncertainty
    for name, v in (("probability_cost_effective", u.probability_cost_effective),
                    ("probability_cost_saving", u.probability_cost_saving)):
        if not (0.0 <= v <= 1.0):
            out.append(Finding(Severity.ERROR, "Probability range",
                               f"uncertainty.{name} is {v}, outside [0, 1]",
                               f"uncertainty.{name}"))
    for c in p.condition_results:
        if not (0.0 <= c.baseline_risk <= 1.0):
            out.append(Finding(Severity.ERROR, "Probability range",
                               f"{c.condition}: baseline risk {c.baseline_risk} "
                               f"outside [0, 1]",
                               f"condition_results[{c.condition}].baseline_risk"))
        if not (0.0 <= c.adherence <= 1.0):
            out.append(Finding(Severity.ERROR, "Probability range",
                               f"{c.condition}: adherence {c.adherence} outside "
                               f"[0, 1]",
                               f"condition_results[{c.condition}].adherence"))
    if u.psa_available and u.nmb_ci_low > u.nmb_ci_high:
        out.append(Finding(Severity.ERROR, "Interval ordering",
                           f"PSA interval is inverted: low ${u.nmb_ci_low:,.0f} "
                           f"> high ${u.nmb_ci_high:,.0f}", "uncertainty"))
    return out


def _check_pooling(p: EconomicsReportPayload) -> list[Finding]:
    out: list[Finding] = []
    c = p.corrections
    if c.naive_cost_averted and c.pooled_cost_averted > c.naive_cost_averted:
        out.append(Finding(
            Severity.ERROR, "Pooling direction",
            f"pooled cost averted ${c.pooled_cost_averted:,.0f} exceeds the "
            f"naive additive figure ${c.naive_cost_averted:,.0f}; pooling must "
            f"reduce, not inflate", "corrections.pooled_cost_averted"))

    for cond in p.condition_results:
        if cond.n_contributing_findings > 1:
            if cond.adherence_adjusted_rrr > cond.pooled_efficacy_rrr + 1e-9:
                out.append(Finding(
                    Severity.ERROR, "Adherence direction",
                    f"{cond.condition}: adherence-adjusted risk reduction "
                    f"{cond.adherence_adjusted_rrr:.3f} exceeds trial efficacy "
                    f"{cond.pooled_efficacy_rrr:.3f}",
                    f"condition_results[{cond.condition}]"))
            if (cond.naive_additive_rrr
                    and cond.pooled_efficacy_rrr > cond.naive_additive_rrr + 1e-9):
                out.append(Finding(
                    Severity.ERROR, "Pooling direction",
                    f"{cond.condition}: pooled efficacy "
                    f"{cond.pooled_efficacy_rrr:.3f} exceeds naive additive "
                    f"{cond.naive_additive_rrr:.3f}",
                    f"condition_results[{cond.condition}]"))
    return out


def _check_wgs_separation(p: EconomicsReportPayload) -> list[Finding]:
    """The check the testing-decision page exists for.

    Two quantities, two bases, and they must not leak into each other. The
    observed figure is a sum over findings this genome contains; the prospective
    figure is a population expectation. This reconstructs the observed sum from
    the findings themselves, so a future change that quietly sources it from the
    population prior fails here.
    """
    out: list[Finding] = []
    t = p.testing_decision

    if t.prospective_basis is not AnalysisBasis.PROSPECTIVE_EXPECTED_YIELD:
        out.append(Finding(Severity.ERROR, "WGS basis tagging",
                           "prospective fields are not tagged as prospective",
                           "testing_decision.prospective_basis"))
    if t.observed_basis is not AnalysisBasis.OBSERVED_FINDINGS:
        out.append(Finding(Severity.ERROR, "WGS basis tagging",
                           "observed fields are not tagged as observed",
                           "testing_decision.observed_basis"))

    wgs_only = [f for f in p.findings
                if f.is_wgs_only and f.pricing_path is PricingPath.VOI_PARAMETRIC]
    n_observed = len(wgs_only)
    sum_observed = sum(f.health_economic_value for f in wgs_only)

    if n_observed != t.observed_wgs_only_findings:
        out.append(Finding(
            Severity.ERROR, "WGS observed count",
            f"testing_decision reports {t.observed_wgs_only_findings} "
            f"sequencing-only finding(s) but the payload carries {n_observed}",
            "testing_decision.observed_wgs_only_findings"))

    if abs(sum_observed - t.observed_wgs_only_value) > _money_tolerance(sum_observed):
        out.append(Finding(
            Severity.ERROR, "WGS observed value contamination",
            f"observed_wgs_only_value is ${t.observed_wgs_only_value:,.0f} but "
            f"the sequencing-only findings in this payload sum to "
            f"${sum_observed:,.0f}; the observed figure must be a sum over "
            f"findings actually present, never the population prior",
            "testing_decision.observed_wgs_only_value"))

    if t.observed_wgs_only_findings == 0 and t.observed_wgs_only_value != 0:
        out.append(Finding(
            Severity.ERROR, "WGS observed value contamination",
            f"no sequencing-only findings were observed, yet "
            f"observed_wgs_only_value is ${t.observed_wgs_only_value:,.0f}",
            "testing_decision.observed_wgs_only_value"))

    expected_net = (t.prospective_gross_expected_value
                    - t.incremental_chip_to_wgs_cost)
    if abs(expected_net - t.prospective_net_expected_value) > _money_tolerance(expected_net):
        out.append(Finding(
            Severity.ERROR, "Prospective WGS identity",
            f"prospective net ${t.prospective_net_expected_value:,.0f} != gross "
            f"${t.prospective_gross_expected_value:,.0f} - incremental cost "
            f"${t.incremental_chip_to_wgs_cost:,.0f}",
            "testing_decision.prospective_net_expected_value"))

    if (t.prospective_number_needed_to_sequence
            and t.prospective_expected_yield > 0):
        implied = round(1.0 / t.prospective_expected_yield)
        if abs(implied - t.prospective_number_needed_to_sequence) > 1:
            out.append(Finding(
                Severity.ERROR, "Number needed to sequence",
                f"reported 1 in {t.prospective_number_needed_to_sequence} but "
                f"the yield {t.prospective_expected_yield:.5f} implies "
                f"1 in {implied}",
                "testing_decision.prospective_number_needed_to_sequence"))

    if t.observed_wgs_only_findings > 0:
        out.append(Finding(
            Severity.INFO, "WGS basis coexistence",
            f"this genome carries {t.observed_wgs_only_findings} observed "
            f"sequencing-only finding(s) worth "
            f"${t.observed_wgs_only_value:,.0f}, and the prospective model "
            f"separately reports 1 in "
            f"{t.prospective_number_needed_to_sequence} with "
            f"${t.prospective_gross_expected_value:,.0f} expected value. Both "
            f"are correct; the renderer must not add them or present the "
            f"population probability as a forecast about this genome",
            "testing_decision"))
    return out


def _check_budget_impact(p: EconomicsReportPayload) -> list[Finding]:
    bi = p.advanced.budget_impact
    if not bi or not bi.get("available"):
        return []
    rows = bi.get("rows") or []
    if not rows:
        return []

    def pmpm_of(r: dict) -> float | None:
        for k in ("pmpm", "net_pmpm", "pmpm_net"):
            if r.get(k) is not None:
                try:
                    return float(r[k])
                except (TypeError, ValueError):
                    return None
        return None

    vals: list[tuple[float, dict]] = [
        (pv, r) for pv, r in ((pmpm_of(r), r) for r in rows) if pv is not None]
    if not vals:
        return [Finding(Severity.WARNING, "Budget-impact peak",
                        "no per-year PMPM values found; peak cannot be verified",
                        "advanced.budget_impact")]

    reported_peak = bi.get("maximum_budget_burden_pmpm",
                           bi.get("peak_pmpm"))
    if reported_peak is None:
        return []
    # PEAK MEANS HIGHEST NET SPEND, not largest magnitude. A budget-impact
    # analysis answers "can the plan afford this?", so the peak is the worst
    # year for the payer — the maximum of the signed net figure. Once the
    # programme turns cost-saving the later years are large *negative* numbers,
    # and selecting on absolute value would report a year of maximum savings as
    # the peak budget impact. This check originally did exactly that and fired
    # against correct engine output; the engine selects with
    # max(rows, key=net_budget_impact), which is the ISPOR convention.
    top_v, top_r = max(vals, key=lambda t: t[0])
    out: list[Finding] = []
    if abs(float(reported_peak) - top_v) > 1e-6:
        out.append(Finding(
            Severity.ERROR, "Budget-impact peak",
            f"reported peak PMPM {reported_peak} is not the maximum of the "
            f"per-year values (max {top_v})", "advanced.budget_impact.peak_pmpm"))
    reported_year = bi.get("maximum_budget_burden_year",
                           bi.get("peak_year"))
    actual_year = top_r.get("year")
    if (reported_year is not None and actual_year is not None
            and str(reported_year) != str(actual_year)):
        out.append(Finding(
            Severity.ERROR, "Budget-impact peak year",
            f"reported peak year {reported_year} does not match the year of "
            f"the maximum PMPM ({actual_year})",
            "advanced.budget_impact.peak_year"))
    return out


def _check_monetised_flag_matches_the_value(
        p: EconomicsReportPayload) -> list[Finding]:
    """`is_monetized` must be false when there is no figure to show.

    The flag was set to a literal True at construction, so a finding whose value
    had been withheld upstream still arrived claiming to be monetised. HFE
    C282Y reported is_monetized=True with canonical_expected_nmb=None after
    NOT_VALUED correctly withheld its figure — the amount was gone and the flag
    still promised one.

    The general shape, which this file now has three instances of: a policy
    enforced at a call site rather than at the point of monetisation is bypassed
    by the next path that reaches monetisation. Duplicate valuation, NOT_VALUED,
    and the curated/parametric split were all the same defect wearing different
    labels. An invariant checked on the assembled payload is enforced once,
    wherever the value came from.
    """
    out: list[Finding] = []
    for f in p.findings:
        if not f.is_monetized:
            continue
        has = any(v for v in (f.canonical_expected_nmb, f.medical_cost_averted,
                              f.monetized_qaly_value, f.legacy_curated_value,
                              f.health_economic_value, f.net_cash))
        if not has:
            out.append(Finding(
                Severity.ERROR, "Monetised flag without a value",
                f"{f.display_name or f.finding_id or '?'} is flagged monetised "
                f"but carries no figure in any value field; the report would "
                f"present it as priced and have nothing to print",
                "findings.is_monetized"))
    return out


# Divergence between the curated and parametric pricing paths. The two model
# different things — curated is a prevalence-weighted lifetime cost of illness,
# parametric a discounted expected NMB over a stated horizon — so a small gap is
# structural and expected. Two orders of magnitude is not.
#
# 3x matches the tolerance _check_structural_crosscheck already applies to the
# Markov cross-check in this file, so the two divergence checks agree on what
# counts as "worth telling the reader". 10x is an order of magnitude above that
# and an order below what was actually observed: the reconciliation file has
# been reporting 97.7x, 93.6x, 80.5x and 59.0x for the same variants, computed
# correctly, written to disk, and recorded as informational. A check that
# measures a discrepancy and does not act on it is worse than no check, because
# it produces the appearance of verification without the substance.
_RECONCILE_ERROR_RATIO = 10.0
_RECONCILE_WARN_RATIO = 3.0


def _check_path_reconciliation(p: EconomicsReportPayload) -> list[Finding]:
    """Curated and parametric prices for one pathway must not diverge wildly."""
    try:
        from report.reconcile import reconcile_paths
    except Exception:
        return []
    out: list[Finding] = []
    for r in reconcile_paths(p):
        if not r.ratio:
            continue
        name = (r.curated_name or r.pathway_id or "?")[:52]
        if r.ratio >= _RECONCILE_ERROR_RATIO:
            out.append(Finding(
                Severity.ERROR, "Pricing paths diverge beyond reconciliation",
                f"{name}: curated ${r.curated_value:,.0f} against parametric "
                f"${r.parametric_value:,.0f} ({r.ratio:.1f}x apart). The two "
                f"paths model different things, but not by this much — one of "
                f"them is wrong about this finding",
                "findings.legacy_curated_value"))
        elif r.ratio >= _RECONCILE_WARN_RATIO:
            out.append(Finding(
                Severity.WARNING, "Pricing paths diverge",
                f"{name}: curated ${r.curated_value:,.0f} against parametric "
                f"${r.parametric_value:,.0f} ({r.ratio:.1f}x apart); reported "
                f"rather than reconciled, because the paths weight prevalence "
                f"and horizon differently",
                "findings.legacy_curated_value"))
    return out


def _check_withheld_findings_carry_no_value(
        p: EconomicsReportPayload) -> list[Finding]:
    """A finding the model declined to price must not be priced anywhere.

    NOT_VALUED was enforced at each pricing path in turn, and each time the next
    path reached monetisation without it. HFE C282Y was withheld by the carrier
    extractor — partner testing is a reproductive action, and pricing one
    prices a prospective child — and simultaneously valued at $5,576 by a second
    carrier block that passed no identity, so the two records could not even be
    linked to notice the contradiction.

    Checked here because the assembled payload is the one place every pricing
    path has already converged. Enforcement at a call site protects that call
    site; enforcement at the point of monetisation protects the report.

    Matching is on the pooling target and the display name as well as the gene,
    because a withheld record and its priced twin frequently carry neither the
    same name nor any gene at all — which is precisely how this survived.
    """
    def _keys(f) -> set[str]:
        out = set()
        for v in ((f.gene or ""), (f.pool_target or ""),
                  (f.condition_id or "")):
            v = str(v).strip().upper()
            if v:
                out.add(v)
        return out

    withheld = [f for f in p.findings if (f.reason_not_monetized or "").strip()]
    if not withheld:
        return []
    withheld_keys: set[str] = set()
    for f in withheld:
        withheld_keys |= _keys(f)

    out: list[Finding] = []
    for f in p.findings:
        if (f.reason_not_monetized or "").strip():
            # The withheld record itself must carry no figure either.
            vals = [("canonical_expected_nmb", f.canonical_expected_nmb),
                    ("legacy_curated_value", f.legacy_curated_value),
                    ("health_economic_value", f.health_economic_value)]
            live = [n for n, v in vals if v]
            if live:
                out.append(Finding(
                    Severity.ERROR, "Withheld finding carries a value",
                    f"{f.display_name or '?'} states a reason for not being "
                    f"monetised but still reports {', '.join(live)}",
                    "findings.reason_not_monetized"))
            continue
        shared = _keys(f) & withheld_keys
        if shared and (f.canonical_expected_nmb or f.legacy_curated_value
                       or f.health_economic_value):
            out.append(Finding(
                Severity.ERROR, "Withheld finding priced by another path",
                f"{f.display_name or '?'} is priced while a finding sharing "
                f"{'/'.join(sorted(shared))} was withheld from monetisation; "
                f"one path declined to value this and another valued it anyway",
                "findings.health_economic_value"))
    return out


def _check_one_gene_one_finding(p: EconomicsReportPayload) -> list[Finding]:
    """A gene must not resolve to more than one monetised finding.

    Two valuation paths used to price ACMG genes independently and disagree by
    200x on the same gene in the same run — one passed the curated QALY figure
    through raw, the other multiplied it by a 0.005 population frequency. Both
    rendered, on different pages, each internally consistent, so nothing
    compared them. The report said two different things about one variant and
    the gate had no opinion.

    Keyed on gene rather than on pathway id deliberately: the ids were derived
    from display text, so rewording a finding produced a new id and two records
    for one gene looked like two findings.
    """
    # A GATE KEYED ON IDENTITY IS BLIND TO RECORDS THAT LACK IT, which is the
    # population it most needs to see. When this check first shipped, gene was
    # populated on 3 of 49 findings, so the duplicate it was written to catch —
    # one threaded copy of LDLR and one untraced copy — sat directly underneath
    # it and did not fire. The six ERRORs it raised against reconstructed
    # pre-fix state gave a misleading impression of its reach.
    #
    # So identity falls back to the gene symbol recoverable from the display
    # name. That is the same text-sniffing this codebase is removing elsewhere,
    # and it is the right call *here* specifically because a validator's job is
    # to catch records that are malformed. Refusing to look at a record until it
    # is well-formed makes the check agree with the bug.
    def _identity(f) -> str:
        g = (f.gene or "").strip().upper()
        if g:
            return g
        name = (f.display_name or f.finding_id or "")
        for tok in name.replace("/", " ").replace("-", " ").split():
            t = tok.strip("().,").upper()
            if len(t) >= 3 and t.isalnum() and any(c.isdigit() for c in t) \
                    and t[0].isalpha():
                return t          # gene-shaped: letters then digits, e.g. MLH1
        return ""

    by_gene: dict[str, list[str]] = {}
    for f in p.findings:
        g = _identity(f)
        if not g or not f.is_monetized:
            continue
        by_gene.setdefault(g, []).append(f.display_name or f.finding_id or "?")
    out: list[Finding] = []
    for gene, names in sorted(by_gene.items()):
        if len(names) > 1:
            out.append(Finding(
                Severity.ERROR, "Gene valued more than once",
                f"{gene} resolves to {len(names)} monetised findings "
                f"({'; '.join(n[:44] for n in names)}). One variant cannot "
                f"carry two dollar figures — one of the paths is "
                f"double-counting it",
                "findings.gene"))
    return out


def _check_render_identity(p: EconomicsReportPayload) -> list[Finding]:
    """A finding rendered in more than one place must render identical values.

    Page 3 and page 4 of one report showed $2,579 vs $6,891 averted, 0.004 vs
    0.124 QALYs, and $100 vs $0 to act for what a reader would read as the same
    finding — while the model had charged $500. Three values for one field, none
    of them the one used. The $0 was the damaging one: it reads as "free to act"
    and inflates net benefit for anyone checking the arithmetic themselves.

    Enforced on identity rather than on position, so it holds however the pages
    are laid out: two records claiming the same gene AND the same condition must
    agree on every reported figure, including intervention cost.
    """
    seen: dict[tuple[str, str], object] = {}
    out: list[Finding] = []
    for f in p.findings:
        key = ((f.gene or "").strip().upper(),
               (f.condition_id or "").strip().lower())
        if not key[0]:
            continue
        prev = seen.get(key)
        if prev is None:
            seen[key] = f
            continue
        for field, label in (("medical_cost_averted", "averted cost"),
                             ("expected_qaly_gain", "QALY gain"),
                             ("intervention_cost", "cost to act"),
                             ("canonical_expected_nmb", "expected NMB")):
            a, b = getattr(prev, field, None), getattr(f, field, None)
            if a is None or b is None or a == b:
                continue
            out.append(Finding(
                Severity.ERROR, "Same finding rendered with different values",
                f"{key[0]} reports {label} as {a} and {b} in the same report; "
                f"a reader comparing two pages would see two answers for one "
                f"finding",
                f"findings.{field}"))
    return out


# ── divergences that are real but not errors ──────────────────────────────────

def _check_structural_crosscheck(p: EconomicsReportPayload) -> list[Finding]:
    sx = p.structural_crosscheck
    if not sx.available:
        return []
    out: list[Finding] = []
    if sx.is_reference_case:
        out.append(Finding(
            Severity.ERROR, "Structural cross-check misclassified",
            "the Markov model is flagged as the reference case; it uses generic "
            "cohort parameters and is not genome-specific",
            "structural_crosscheck.is_reference_case"))

    ref_q = p.reference_case.incremental_qalys
    if ref_q and sx.markov_qaly_gain:
        factor = sx.markov_qaly_gain / ref_q
        if factor > 3.0 or factor < (1 / 3.0):
            out.append(Finding(
                Severity.WARNING, "Structural cross-check divergence",
                f"the Markov cross-check reports {sx.markov_qaly_gain:.3f} "
                f"incremental QALYs against the reference case's {ref_q:.4f} "
                f"({factor:.1f}x). Both may be correct — they model different "
                f"cohorts over different horizons — but the two must not be "
                f"rendered at equal visual weight, and the cross-check must "
                f"not be presented as this genome's result",
                "structural_crosscheck.markov_qaly_gain"))
    return out


def _check_pricing_path_divergence(p: EconomicsReportPayload) -> list[Finding]:
    """What remains of C3 after canonicalisation.

    Findings are now merged one record per economic pathway, with the
    registry-backed parametric expected NMB canonical and the curated figure
    demoted to an audit field. A pathway therefore cannot carry two competing
    values any more, and the old "which of these two numbers is it" warning
    describes a world that no longer exists.

    Three things are still worth saying: which findings have no standardised
    value at all, which canonical records still hang off an identifier
    slugified from display text, and how far the retained curated figures sit
    from the canonical ones.
    """
    out: list[Finding] = []

    unstandardised = [f for f in p.findings
                      if f.is_monetized and f.canonical_expected_nmb is None]
    if unstandardised:
        out.append(Finding(
            Severity.WARNING, "Findings without a standardised value",
            f"{len(unstandardised)} monetised finding(s) have no registry-backed "
            f"parametric pathway and therefore no canonical expected NMB: "
            + ", ".join(sorted(f.display_name for f in unstandardised)[:4])
            + '. They render with finding, action and evidence intact and '
              '"not yet standardised" in place of a value; no figure is '
              'invented for them', "findings.canonical_expected_nmb"))

    legacy_ids = [f for f in p.findings
                  if f.canonical_expected_nmb is not None
                  and f.pathway_id_is_legacy]
    if legacy_ids:
        out.append(Finding(
            Severity.WARNING, "Canonical records on legacy identifiers",
            f"{len(legacy_ids)} canonical finding(s) still key on an id derived "
            f"from display text, so rewording the finding would move the id: "
            + ", ".join(sorted(f.display_name[:40] for f in legacy_ids)[:4])
            + ". Those extractors need semantic components, as the "
              "pharmacogenomics extractor now has",
            "findings.economic_pathway_id"))

    for f in p.findings:
        if f.canonical_expected_nmb is None or f.legacy_curated_value is None:
            continue
        a, b = f.legacy_curated_value, f.canonical_expected_nmb
        if a == 0 and b == 0:
            continue
        sign_flip = (a >= 0) != (b >= 0)
        if sign_flip or abs(a - b) > max(1000.0, abs(b) * 0.5):
            out.append(Finding(
                Severity.INFO, "Legacy curated value diverges",
                f"'{f.display_name}': canonical expected NMB ${b:,.0f}, legacy "
                f"curated figure ${a:,.0f}"
                + (" — opposite signs" if sign_flip else "")
                + ". The curated figure is retained for audit only and is not "
                  "an alternative NMB: it conditions on genotype prevalence "
                  "rather than event probability. See report/reconcile.py",
                f"findings[{f.economic_pathway_id}]"))
    return out


def _check_provenance(p: EconomicsReportPayload) -> list[Finding]:
    out: list[Finding] = []
    pr = p.provenance
    if pr.registry_pct_sourced and pr.registry_pct_sourced < 75.0:
        out.append(Finding(
            Severity.ERROR, "Registry provenance floor",
            f"registry is {pr.registry_pct_sourced:.1f}% sourced, below the "
            f"75% floor asserted in tests/unit/test_econ_params.py",
            "provenance.registry_pct_sourced"))
    if pr.model_pct_unsourced and pr.model_pct_unsourced > 10.0:
        out.append(Finding(
            Severity.WARNING, "Unsourced model parameters",
            f"{pr.model_pct_unsourced:.1f}% of known model parameters carry no "
            f"attribution at all", "provenance.model_pct_unsourced"))
    if pr.registry_pct_sourced and pr.model_pct_resolvable:
        out.append(Finding(
            Severity.INFO, "Two provenance denominators",
            f"registry coverage is {pr.registry_pct_sourced:.1f}% of "
            f"{pr.registry_n_parameters} registered parameters; whole-model "
            f"coverage is {pr.model_pct_resolvable:.1f}% of "
            f"{pr.model_n_total_known} known parameters. These are different "
            f"denominators and must each be labelled as such — collapsing them "
            f"into one figure would create an inconsistency, not remove one",
            "provenance"))
    if pr.weighted_provenance is None:
        out.append(Finding(
            Severity.INFO, "Weighted provenance not implemented",
            "sensitivity-weighted evidence scoring is deferred; raw provenance "
            "categories are reported instead. See the TODO on "
            "Provenance.weighted_provenance", "provenance.weighted_provenance"))

    odd = sorted({f.raw_evidence_confidence for f in p.findings
                  if f.raw_evidence_confidence
                  and f.raw_evidence_confidence.strip().lower()
                  not in ("high", "moderate", "low")})
    if odd:
        out.append(Finding(
            Severity.WARNING, "Unrecognised confidence vocabulary",
            f"producing modules emitted confidence value(s) {odd} that the "
            f"economics layer does not know. _CONFIDENCE_CONC in "
            f"value_of_information.py has keys high/moderate/low only, so "
            f"these silently take the moderate concentration in the "
            f"probabilistic analysis. Display is normalised here; the PSA "
            f"spread is not", "findings.raw_evidence_confidence"))

    low = [f for f in p.findings
           if (f.evidence_confidence or "").lower() == "low" and f.is_monetized]
    if low:
        out.append(Finding(
            Severity.WARNING, "Low-confidence findings carry value",
            f"{len(low)} monetised finding(s) rest on low-confidence evidence: "
            + ", ".join(sorted({f.display_name for f in low})[:5]),
            "findings"))
    return out


def _note_legitimate_differences(p: EconomicsReportPayload) -> list[Finding]:
    """Differences that are correct and must be explained, never equalised."""
    out: list[Finding] = []
    r, u = p.reference_case, p.uncertainty
    if u.psa_available and r.nmb and u.psa_mean_nmb:
        out.append(Finding(
            Severity.INFO, "Deterministic vs probabilistic NMB",
            f"reference-case NMB is ${r.nmb:,.0f}; the mean across "
            f"{u.psa_iterations:,} probabilistic simulations is "
            f"${u.psa_mean_nmb:,.0f}. Both are correct. Neither may be "
            f"labelled 'Net monetary benefit' without a qualifier",
            "reference_case.nmb vs uncertainty.psa_mean_nmb"))

    pooled_total = sum(c.nmb for c in p.condition_results)
    per_finding_total = sum(f.canonical_expected_nmb or 0.0
                            for f in p.findings)
    if pooled_total and per_finding_total:
        out.append(Finding(
            Severity.INFO, "Per-finding vs pooled totals",
            f"per-finding values sum to ${per_finding_total:,.0f}; the pooled "
            f"condition results sum to ${pooled_total:,.0f}. Correlated "
            f"findings are combined on the risk scale and each condition is "
            f"charged one cost of illness, so per-finding values are not "
            f"additive and must not be presented as though they were",
            "findings vs condition_results"))
    return out


# ── presentation ──────────────────────────────────────────────────────────────

def format_report(findings: list[dict[str, Any]]) -> str:
    """A plain-text summary, for CLI output and test failure messages."""
    if not findings:
        return "No report-consistency findings."
    lines = []
    for sev in (Severity.ERROR, Severity.WARNING, Severity.INFO):
        group = [f for f in findings if f.get("severity") == sev]
        if not group:
            continue
        lines.append(f"{sev} ({len(group)})")
        for f in group:
            lines.append(f"  - [{f.get('check')}] {f.get('detail')}")
    return "\n".join(lines)
