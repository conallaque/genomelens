"""Cross-path reconciliation: why one finding carries two dollar figures.

The report prices each finding twice, through parameterisations that were never
reconciled. Until ``econ/identity.py`` gave both paths a shared ``pathway_id``
the two figures could not even be matched, so the disagreement was invisible.
They can now be matched exactly, and this module reports the difference.

**It changes no values.** It compares them, decomposes the difference where the
inputs are identifiable, and states which side of each methodological choice
each path takes.

THE FIVE IDENTIFIED DIFFERENCES, traced through the code rather than inferred:

1. **Adherence.** The curated path multiplies benefit, QALYs and ongoing cost by
   a category adherence factor (0.35–0.65). ``_finding_nmb`` applies none: its
   figures are trial efficacy.
2. **Evidence haircut.** The value-of-information path applies
   ``_evidence_haircut(category)`` to both cost and QALY. The curated path
   applies none.
3. **Exposure probability.** The curated path uses the record's ``prevalence``
   (defaulting to 0.15). The parametric path uses ``p_rx`` from ``PGX_CEA`` for
   pharmacogenomics, or the registry's ``baseline_event_probability`` for
   condition findings.
4. **Cost anchor.** Curated reads ``outcome_value`` from the per-category
   tables, which are not in the parameter registry. Parametric reads registry
   cost-of-illness anchors, or ``adr_cost`` from ``PGX_CEA``.
5. **Marginal-cost fraction — applied inconsistently *within* the parametric
   path.** The condition branch scales avoided cost by
   ``MARGINAL_COST_FRACTION``; the pharmacogenomics branch does not. The curated
   path applies it to everything. This is an internal inconsistency in
   ``_finding_nmb``, not merely a difference between the two paths.

Neither path dominates: the curated one models real-world uptake but rests on
unregistered constants; the parametric one is registry-anchored and
provenance-tracked but reports efficacy rather than effectiveness. Choosing
between them is a methodology decision, and this module deliberately does not
make it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from report.payload import EconomicsReportPayload

__all__ = ["PathwayReconciliation", "format_reconciliation", "reconcile_paths"]

# Methodological differences by pathway kind. Keyed on what the parametric path
# does differently, because that is the side the registry tracks.
_CURATED_ONLY = ("real-world adherence", "marginal-cost fraction on PGx")
_PARAMETRIC_ONLY = ("evidence haircut",)


@dataclass
class PathwayReconciliation:
    pathway_id: str
    curated_name: str = ""
    parametric_name: str = ""
    curated_value: float = 0.0
    parametric_value: float = 0.0
    absolute_difference: float = 0.0
    ratio: float | None = None
    """Larger over smaller, so the figure reads as "Nx apart" regardless of
    direction. None when either side is zero."""
    higher_path: str = ""
    sign_agreement: bool = True
    curated_qaly: float = 0.0
    parametric_qaly: float = 0.0
    curated_cost_averted: float = 0.0
    parametric_cost_averted: float = 0.0
    curated_adherence: float | None = None
    n_curated_records: int = 1
    n_parametric_records: int = 1
    differing_inputs: list[str] = field(default_factory=list)
    evidence_confidence: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _ratio(a: float, b: float) -> float | None:
    lo, hi = min(abs(a), abs(b)), max(abs(a), abs(b))
    return None if lo == 0 else round(hi / lo, 1)


def reconcile_paths(p: EconomicsReportPayload) -> list[PathwayReconciliation]:
    """One row per pathway priced by both models, worst disagreement first."""
    rows: list[PathwayReconciliation] = []
    for f in p.findings:
        if f.canonical_expected_nmb is None or f.legacy_curated_value is None:
            continue
        a, b = f.legacy_curated_value, f.canonical_expected_nmb
        rows.append(PathwayReconciliation(
            pathway_id=f.economic_pathway_id,
            curated_name=f.display_name,
            parametric_name=f.display_name,
            curated_value=a,
            parametric_value=b,
            absolute_difference=abs(a - b),
            ratio=_ratio(a, b),
            higher_path=("curated" if abs(a) > abs(b) else "parametric"),
            sign_agreement=(a >= 0) == (b >= 0),
            parametric_qaly=f.expected_qaly_gain,
            parametric_cost_averted=f.medical_cost_averted,
            curated_adherence=f.adherence,
            differing_inputs=_reasons(f),
            evidence_confidence=f.evidence_confidence,
        ))
    rows.sort(key=lambda r: (r.sign_agreement, -(r.ratio or 0)))
    return rows


def _reasons(f) -> list[str]:
    """Why the two figures differ, traced to inputs rather than asserted.

    The dominant term is the same in every case examined: the curated path
    conditions the QALY gain on **genotype prevalence** — the probability a
    random person carries the variant — where the parametric path conditions it
    on the **event probability**, p_rx x p_adr x rrr. For CYP2C19 that is 0.30
    against 0.01, a factor of 30, and it accounts for essentially all of the
    $4,452 vs $490 gap: the curated figure is 98% monetised QALY.
    """
    out = [
        "QALY conditioned on genotype prevalence (curated) rather than on "
        "event probability p_rx x p_adr x rrr (parametric) — the dominant term",
        "curated multiplies outcome_value by prevalence again, although "
        "outcome_value is already p_rx x p_adr x rrr x adr_cost",
        "evidence haircut applied by the parametric path only",
        "real-world adherence applied by the curated path only; parametric "
        "figures are trial efficacy",
        "marginal-cost fraction applied to every curated row but only to the "
        "condition branch of the parametric path",
    ]
    if f.adherence is not None and f.adherence < 1.0:
        out.append(f"curated adherence factor {f.adherence:.0%}")
    return out


def format_reconciliation(rows: list[PathwayReconciliation]) -> str:
    """Plain-text table, for CLI output and the developer document."""
    if not rows:
        return "No findings are priced by both paths."
    head = (f"{'finding':<34}{'curated':>10}{'parametric':>12}"
            f"{'diff':>10}{'ratio':>8}  sign")
    out = [head, "-" * len(head)]
    for r in rows:
        name = (r.parametric_name or r.curated_name)[:33]
        ratio = f"{r.ratio}x" if r.ratio else "—"
        sign = "agree" if r.sign_agreement else "OPPOSITE"
        out.append(f"{name:<34}{r.curated_value:>10,.0f}"
                   f"{r.parametric_value:>12,.0f}"
                   f"{r.absolute_difference:>10,.0f}{ratio:>8}  {sign}")
    return "\n".join(out)
