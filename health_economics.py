"""
Health Economics — Clinical & Payer ROI for Genomic Interventions
=================================================================

Translates the genomic findings already produced by the pipeline (PGx
phenotypes, polygenic-risk tiers, APOE genotype) into a health-economic
view: what does acting on each finding cost, what adverse outcome does it
avert, and what is the return on that spend — for an individual, for a
clinic, and for a payer (insurer) covering a large member population.

This module is **decision-support economics, not accounting**. Every
dollar figure is an order-of-magnitude estimate drawn from published
pharmacoeconomic and screening literature. Actual ROI varies enormously by
clinic, payer, patient adherence, drug mix, and local costs. The numbers
here are meant to size opportunity, not to bill against.

Input shape
-----------
``analyze_health_economics(findings, snps_df)`` consumes the same dict the
pipeline already persists to ``tier1_results.json`` — specifically the
``pgx_summary``, ``prs_summary`` and ``apoe_genotype`` keys — so it can be
run either inside the pipeline or, for testing, straight against that JSON.
``snps_df`` is accepted for SNP-level confirmation (APOE / ACTN3 via the
unified ``snp_registry``); it is optional and the module degrades
gracefully when a SNP is absent.

Output shape
------------
``{status, findings_with_economics, clinic_dashboard, payer_impact,
disclaimer, ...}`` — see ``analyze_health_economics`` for the contract.

A finding is suppressed (not an error) when its economics are unknown, when
the genotype is normal/indeterminate, or when the polygenic tier is not
elevated. Missing data never raises; it just yields fewer findings.
"""

from __future__ import annotations

from typing import Dict, List, Optional
import pandas as pd

try:
    import snp_registry  # optional SNP-level confirmation (APOE / ACTN3)
except Exception:  # pragma: no cover - registry should always import
    snp_registry = None  # type: ignore


# ─── Economic model constants (documented assumptions) ─────────────────────
# All are estimates; surfaced in the disclaimer and editable in one place.

DISCOUNT_RATE = 0.03          # standard health-economic real discount rate
NPV_HORIZON_YEARS = 3         # NPV evaluated over a 3-year window

# Clinic revenue model — a genomics-enabled clinic offering this as a
# subscription service. These are business-model assumptions, not clinical.
CLINIC_REVENUE_MONTHLY = 500.0   # $/patient/month subscription
CLINIC_GROSS_MARGIN = 0.85       # gross margin on that subscription

DEFAULT_CLINIC_PATIENTS = 100
DEFAULT_PAYER_MEMBERS = 100_000


# ─── Pharmacogenomic intervention economics (one-time test cost) ───────────
# Keyed by the CPIC gene as it appears in pgx_summary. Each entry:
#   cost            — $ for the PGx test + consult (one-time, upfront)
#   outcome_value   — $ health-economic value of the averted adverse event
#   clinical_benefit— plain-language description of what is prevented
#   prevalence      — fraction of a general population with an *actionable*
#                     (atypical) phenotype for this gene (payer scaling)
#   qaly_gain       — QALYs gained per averted event (cost-per-QALY)
# Sources: CPIC guidelines + published PGx cost-effectiveness analyses;
# rounded to order-of-magnitude. ROI is derived (outcome_value / cost).
PGX_ECONOMICS: Dict[str, Dict] = {
    "CYP2C9": {
        "drug": "warfarin",
        "clinical_benefit": "Prevent major bleeding event from warfarin overdosing",
        "cost": 300, "outcome_value": 15_000,
        "prevalence": 0.35, "qaly_gain": 0.30,
    },
    "VKORC1": {
        "drug": "warfarin",
        "clinical_benefit": "Improve warfarin dosing accuracy, fewer INR excursions",
        "cost": 150, "outcome_value": 8_000,
        "prevalence": 0.37, "qaly_gain": 0.15,
    },
    "CYP2C19": {
        "drug": "clopidogrel",
        "clinical_benefit": "Avoid clopidogrel non-response / stent thrombosis (MACE)",
        "cost": 250, "outcome_value": 18_000,
        "prevalence": 0.30, "qaly_gain": 0.35,
    },
    "CYP2D6": {
        "drug": "codeine / antidepressants",
        "clinical_benefit": "Avoid opioid toxicity or ineffective analgesia/SSRI dosing",
        "cost": 250, "outcome_value": 12_000,
        "prevalence": 0.30, "qaly_gain": 0.20,
    },
    "CYP3A5": {
        "drug": "tacrolimus",
        "clinical_benefit": "Hit tacrolimus target faster, avoid rejection/toxicity",
        "cost": 250, "outcome_value": 10_000,
        "prevalence": 0.40, "qaly_gain": 0.20,
    },
    "TPMT": {
        "drug": "thiopurines (azathioprine / 6-MP)",
        "clinical_benefit": "Prevent serious myelosuppression from thiopurines",
        "cost": 200, "outcome_value": 25_000,
        "prevalence": 0.10, "qaly_gain": 0.40,
    },
    "NUDT15": {
        "drug": "thiopurines (azathioprine / 6-MP)",
        "clinical_benefit": "Prevent thiopurine-induced leukopenia (esp. East Asian ancestry)",
        "cost": 200, "outcome_value": 25_000,
        "prevalence": 0.10, "qaly_gain": 0.40,
    },
    "SLCO1B1": {
        "drug": "simvastatin",
        "clinical_benefit": "Prevent statin-induced myopathy / rhabdomyolysis",
        "cost": 150, "outcome_value": 9_000,
        "prevalence": 0.25, "qaly_gain": 0.15,
    },
    "UGT1A1": {
        "drug": "irinotecan",
        "clinical_benefit": "Prevent severe irinotecan neutropenia / diarrhea",
        "cost": 200, "outcome_value": 14_000,
        "prevalence": 0.15, "qaly_gain": 0.30,
    },
    "HLA-B*57:01": {
        "drug": "abacavir",
        "clinical_benefit": "Prevent abacavir hypersensitivity reaction",
        "cost": 150, "outcome_value": 50_000,
        "prevalence": 0.06, "qaly_gain": 0.50,
    },
}

# Phenotype keywords that make a PGx result *actionable* (i.e. it would
# change prescribing). Normal-function / indeterminate results are not
# counted as realized economic findings.
_ACTIONABLE_PGX_KEYWORDS = (
    "poor", "intermediate", "rapid", "ultrarapid", "decreased",
    "increased", "non-expressor", "high warfarin sensitivity",
    "positive", "deficien",
)


# ─── Polygenic / lifestyle intervention economics (recurring cost) ─────────
# Keyed by the PRS panel name in prs_summary (and APOE handled separately).
# Cost here is an *annual* program cost; outcome_value is the averted
# event's value. Triggered when the polygenic tier is Elevated or High.
PRS_ECONOMICS: Dict[str, Dict] = {
    "Coronary Artery Disease": {
        "finding": "Elevated coronary-artery-disease polygenic risk",
        "clinical_benefit": "Intensive lipid management (statin) to prevent MI",
        "intervention": "Intensive statin therapy + lipid monitoring",
        "cost": 500, "outcome_value": 250_000,
        "prevalence": 0.20, "qaly_gain": 1.50, "recurring": True,
    },
    "Type 2 Diabetes": {
        "finding": "Elevated type-2-diabetes polygenic risk",
        "clinical_benefit": "Prevent / delay T2D onset",
        "intervention": "CGM + lifestyle coaching",
        "cost": 3_600, "outcome_value": 50_000,
        "prevalence": 0.20, "qaly_gain": 0.80, "recurring": True,
    },
    "BMI / Obesity Tendency": {
        "finding": "Elevated obesity-tendency polygenic risk",
        "clinical_benefit": "Prevent obesity-driven T2D via structured program",
        "intervention": "Structured weight-management program",
        "cost": 1_200, "outcome_value": 30_000,
        "prevalence": 0.20, "qaly_gain": 0.40, "recurring": True,
    },
}

# APOE ε4 carriers — high CAD/AD risk; same statin-style intervention as the
# CAD PRS path but driven by genotype. Kept separate so it can fire even
# when no CAD PRS panel is present.
APOE_E4_ECONOMICS = {
    "finding": "APOE ε4 carrier (elevated cardiovascular / Alzheimer's risk)",
    "clinical_benefit": "Intensive cardiovascular risk reduction (statin) to prevent MI",
    "intervention": "Intensive statin therapy + lipid monitoring",
    "cost": 500, "outcome_value": 250_000,
    "prevalence": 0.25, "qaly_gain": 1.50, "recurring": True,
}

# Optional exercise / longevity inputs (only used if the caller supplies
# them in `findings`). Kept here so the economic assumptions live together.
VO2MAX_ECONOMICS = {
    "finding": "Low cardiorespiratory fitness (VO2max) capacity",
    "clinical_benefit": "Structured training to raise VO2max (all-cause mortality)",
    "intervention": "Supervised aerobic training program",
    "cost": 2_000, "outcome_value": 40_000,
    "prevalence": 0.30, "qaly_gain": 0.50, "recurring": True,
}
LONGEVITY_VALUE_PER_PERCENTILE = 10_000  # $ per 1-percentile longevity gain


DISCLAIMER = (
    "These are estimates based on published pharmacoeconomic literature, not "
    "guarantees. Actual ROI varies by clinic, payer, patient adherence, drug "
    "mix, and local costs. Intervention costs are one-time for pharmacogenomic "
    "tests and annual for lifestyle programs; outcome values are the modeled "
    "health-economic value of an averted adverse event. Figures size "
    "opportunity for decision support — they are not billing or actuarial values."
)


# ─── Core economic math ────────────────────────────────────────────────────

def calculate_roi(cost: float, outcome_value: float) -> Optional[float]:
    """ROI as a simple ratio (outcome value / cost). None if cost is zero."""
    if not cost or cost <= 0:
        return None
    return round(outcome_value / cost, 1)


def calculate_payback_months(cost: float, outcome_value: float) -> Optional[float]:
    """Months for the averted-outcome value to repay the intervention cost,
    treating ``outcome_value`` as an annualized benefit."""
    if not outcome_value or outcome_value <= 0:
        return None
    return round((cost / outcome_value) * 12.0, 2)


def calculate_npv(
    cost: float,
    outcome_value: float,
    recurring_cost: bool = False,
    horizon: int = NPV_HORIZON_YEARS,
    rate: float = DISCOUNT_RATE,
) -> float:
    """Net present value over ``horizon`` years at ``rate``.

    Benefit (``outcome_value``) accrues each year and is discounted. Cost is
    either a one-time upfront spend at t=0 (pharmacogenomic tests) or a
    recurring annual spend discounted over the horizon (lifestyle programs).
    """
    benefit = sum(outcome_value / (1 + rate) ** t for t in range(1, horizon + 1))
    if recurring_cost:
        spend = sum(cost / (1 + rate) ** t for t in range(1, horizon + 1))
    else:
        spend = cost  # upfront, undiscounted at t=0
    return round(benefit - spend, 2)


def _econ_record(
    finding: str,
    clinical_benefit: str,
    cost: float,
    outcome_value: float,
    confidence: str,
    *,
    recurring: bool = False,
    source: str = "",
    prevalence: float = 0.0,
    qaly_gain: float = 0.0,
    evidence: str = "",
) -> Dict:
    """Assemble one finding's economics record with derived metrics."""
    return {
        "finding": finding,
        "clinical_benefit": clinical_benefit,
        "intervention_cost": cost,
        "cost_basis": "annual" if recurring else "one-time",
        "outcome_value": outcome_value,
        "roi": calculate_roi(cost, outcome_value),
        "payback_months": calculate_payback_months(cost, outcome_value),
        "npv_3year": calculate_npv(cost, outcome_value, recurring_cost=recurring),
        "confidence": confidence,
        "category": source,
        "prevalence": prevalence,
        "qaly_gain": qaly_gain,
        "evidence": evidence,
    }


# ─── Finding extractors ────────────────────────────────────────────────────

def _is_actionable_pgx(phenotype: str) -> bool:
    p = (phenotype or "").lower()
    if not p or "indeterminate" in p or "normal" in p:
        # "Normal Warfarin Sensitivity" etc. — not actionable.
        return False
    return any(k in p for k in _ACTIONABLE_PGX_KEYWORDS)


def _pgx_findings(pgx_summary: Dict) -> List[Dict]:
    out: List[Dict] = []
    for gene, info in (pgx_summary or {}).items():
        econ = PGX_ECONOMICS.get(gene)
        if econ is None:
            continue  # no cost data → suppress (graceful)
        phenotype = (info or {}).get("phenotype", "")
        if not _is_actionable_pgx(phenotype):
            continue
        out.append(_econ_record(
            finding=f"{gene} {phenotype} ({econ['drug']})",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="high", source="Pharmacogenomics",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{gene} phenotype: {phenotype}",
        ))
    return out


def _prs_findings(prs_summary: Dict) -> List[Dict]:
    out: List[Dict] = []
    for panel, info in (prs_summary or {}).items():
        econ = PRS_ECONOMICS.get(panel)
        if econ is None:
            continue
        tier = (info or {}).get("tier")
        if tier not in ("Elevated", "High"):
            continue
        pct = (info or {}).get("percentile")
        ev = f"{panel} polygenic tier: {tier}"
        if pct is not None:
            ev += f" ({pct:g}th percentile)"
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=f"{econ['clinical_benefit']} — {econ['intervention']}",
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate", source="Polygenic Risk", recurring=econ["recurring"],
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=ev,
        ))
    return out


def _apoe_genotype_from_snps(snps_df: Optional[pd.DataFrame]) -> Optional[str]:
    """Best-effort APOE ε4 confirmation via the registry SNP, used only when
    the caller did not supply an apoe_genotype. Returns None if undetermined."""
    if snps_df is None or snp_registry is None:
        return None
    try:
        if "rs429358" not in snps_df.index:
            return None
        gt = str(snps_df.loc["rs429358"].get("genotype", "")).upper()
        gt = gt.replace(" ", "").replace("-", "")
        if "C" in gt:  # C allele at rs429358 marks an ε4 haplotype
            return "e4-carrier (rs429358 C)"
    except Exception:
        return None
    return None


def _apoe_findings(apoe_genotype: Optional[str], snps_df) -> List[Dict]:
    geno = apoe_genotype
    confirm = ""
    if not geno:
        geno = _apoe_genotype_from_snps(snps_df)
        confirm = " (rs429358-derived)"
    if not geno or "4" not in str(geno) and "e4" not in str(geno).lower():
        return []
    e = APOE_E4_ECONOMICS
    return [_econ_record(
        finding=e["finding"],
        clinical_benefit=f"{e['clinical_benefit']} — {e['intervention']}",
        cost=e["cost"], outcome_value=e["outcome_value"],
        confidence="moderate", source="Genotype", recurring=e["recurring"],
        prevalence=e["prevalence"], qaly_gain=e["qaly_gain"],
        evidence=f"APOE genotype: {geno}{confirm}",
    )]


def _exercise_longevity_findings(findings: Dict) -> List[Dict]:
    """Optional findings the caller may supply: low VO2max tier, and a
    longevity percentile to quantify improvement headroom."""
    out: List[Dict] = []
    vo2 = str(findings.get("vo2max_tier", "")).lower()
    if vo2 in ("low", "poor", "below average"):
        e = VO2MAX_ECONOMICS
        out.append(_econ_record(
            finding=e["finding"],
            clinical_benefit=f"{e['clinical_benefit']} — {e['intervention']}",
            cost=e["cost"], outcome_value=e["outcome_value"],
            confidence="moderate", source="Exercise / Lifestyle",
            recurring=e["recurring"], prevalence=e["prevalence"],
            qaly_gain=e["qaly_gain"], evidence=f"VO2max tier: {vo2}",
        ))
    pct = findings.get("longevity_percentile")
    if isinstance(pct, (int, float)) and 0 <= pct < 50:
        headroom = 50 - pct  # percentiles of achievable improvement to median
        value = round(headroom * LONGEVITY_VALUE_PER_PERCENTILE)
        cost = 1_500  # representative annual longevity-program cost
        out.append(_econ_record(
            finding="Below-median longevity composite (improvement headroom)",
            clinical_benefit=(
                f"~{headroom:g} percentile points of achievable longevity gain "
                f"(~${LONGEVITY_VALUE_PER_PERCENTILE:,}/percentile)"),
            cost=cost, outcome_value=value,
            confidence="low", source="Longevity", recurring=True,
            prevalence=0.50, qaly_gain=round(headroom * 0.02, 2),
            evidence=f"Longevity percentile: {pct:g}",
        ))
    return out


# ─── Scaling ───────────────────────────────────────────────────────────────

def scale_to_clinic(findings_econ: Dict, patient_count: int = DEFAULT_CLINIC_PATIENTS) -> Dict:
    """Scale per-finding economics to a clinic population and overlay a
    subscription revenue model."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return {"patient_count": patient_count, "n_findings": 0,
                "note": "No actionable findings with economics for this profile."}

    costs = [f["intervention_cost"] for f in findings]
    benefits = [f["outcome_value"] for f in findings]
    avg_cost = round(sum(costs) / len(costs), 2)
    avg_benefit = round(sum(benefits) / len(benefits), 2)

    monthly_margin = CLINIC_REVENUE_MONTHLY * CLINIC_GROSS_MARGIN
    payback = round(avg_cost / monthly_margin, 1) if monthly_margin else None

    return {
        "patient_count": patient_count,
        "n_findings": len(findings),
        "avg_cost_per_patient": avg_cost,
        "avg_benefit_per_patient": avg_benefit,
        "avg_roi": calculate_roi(avg_cost, avg_benefit),
        "total_cost": round(avg_cost * patient_count, 2),
        "total_benefit": round(avg_benefit * patient_count, 2),
        "revenue_model_monthly": CLINIC_REVENUE_MONTHLY,
        "gross_margin": CLINIC_GROSS_MARGIN,
        "payback_period_months": payback,
        "summary": (
            f"Applied to {patient_count} patients: "
            f"cost ${round(avg_cost * patient_count):,}, "
            f"modeled benefit ${round(avg_benefit * patient_count):,}, "
            f"ROI {calculate_roi(avg_cost, avg_benefit)}:1"
        ),
    }


def scale_to_payer(findings_econ: Dict, member_population: int = DEFAULT_PAYER_MEMBERS) -> Dict:
    """Scale per-finding economics to a payer's member population using each
    finding's population prevalence."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return {"member_population": member_population, "affected_members": 0,
                "note": "No actionable findings with economics for this profile."}

    affected_total = 0
    total_cost = 0.0
    total_benefit = 0.0
    total_qalys = 0.0
    per_finding: List[Dict] = []
    for f in findings:
        affected = round(member_population * f.get("prevalence", 0.0))
        cost = affected * f["intervention_cost"]
        benefit = affected * f["outcome_value"]
        qalys = affected * f.get("qaly_gain", 0.0)
        affected_total += affected
        total_cost += cost
        total_benefit += benefit
        total_qalys += qalys
        per_finding.append({
            "finding": f["finding"], "affected_members": affected,
            "total_cost": round(cost), "total_benefit": round(benefit),
        })

    cost_per_qaly = round(total_cost / total_qalys) if total_qalys else None
    return {
        "member_population": member_population,
        "affected_members": affected_total,
        "total_cost": round(total_cost),
        "total_benefit": round(total_benefit),
        "roi": calculate_roi(total_cost, total_benefit),
        "cost_per_qaly": cost_per_qaly,
        "net_savings": round(total_benefit - total_cost),
        "per_finding": per_finding,
        "summary": (
            f"Applied to {member_population:,} members: "
            f"{affected_total:,} interventions, cost ${round(total_cost):,}, "
            f"modeled savings ${round(total_benefit):,}, "
            f"ROI {calculate_roi(total_cost, total_benefit)}:1"
        ),
    }


# ─── Markdown summary ──────────────────────────────────────────────────────

def generate_economics_summary(findings_econ: Dict) -> str:
    """Human-readable markdown for the HTML report."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return ("## Health Economics\n\nNo actionable genomic findings with "
                "modeled economics for this profile.\n")

    # Main summary shows only high-confidence findings; if none are
    # high-confidence we fall back to the best available with a note so the
    # summary is never silently empty.
    high_conf = findings_econ.get("high_confidence") or [
        f for f in findings if f.get("confidence") == "high"
    ]
    main = high_conf if high_conf else findings
    top = sorted(main, key=lambda f: (f.get("roi") or 0), reverse=True)[:3]
    lines = ["## Health Economics — Clinical & Payer ROI", ""]
    heading = ("### Top high-confidence interventions by ROI" if high_conf
               else "### Top interventions by ROI (no high-confidence findings)")
    lines.append(heading)
    for f in top:
        lines.append(
            f"- **{f['finding']}** — {f['clinical_benefit']}. "
            f"Cost ${f['intervention_cost']:,} ({f['cost_basis']}), "
            f"value ${f['outcome_value']:,}, **ROI {f['roi']}:1**, "
            f"payback {f['payback_months']} mo, NPV(3y) ${f['npv_3year']:,.0f}."
        )

    clinic = findings_econ.get("clinic_dashboard", {})
    if clinic.get("n_findings"):
        lines += ["", "### Clinic dashboard", clinic.get("summary", ""),
                  f"- Avg cost/patient: ${clinic['avg_cost_per_patient']:,} · "
                  f"avg benefit/patient: ${clinic['avg_benefit_per_patient']:,} · "
                  f"avg ROI {clinic['avg_roi']}:1",
                  f"- Subscription revenue ${clinic['revenue_model_monthly']:,.0f}/mo "
                  f"at {int(clinic['gross_margin']*100)}% margin · "
                  f"payback {clinic['payback_period_months']} mo"]

    payer = findings_econ.get("payer_impact", {})
    if payer.get("affected_members"):
        lines += ["", "### Payer impact", payer.get("summary", ""),
                  f"- Cost per QALY: "
                  f"${payer['cost_per_qaly']:,}" if payer.get("cost_per_qaly")
                  else "- Cost per QALY: n/a",
                  f"- Net modeled savings: ${payer['net_savings']:,}"]

    lines += ["", f"_{DISCLAIMER}_"]
    return "\n".join(lines)


# ─── Master analyzer ───────────────────────────────────────────────────────

def analyze_health_economics(findings: Dict, snps_df: pd.DataFrame) -> Dict:
    """Compute clinical & payer ROI for genomic interventions.

    Parameters
    ----------
    findings : dict
        Pipeline findings in ``tier1_results.json`` shape — uses
        ``pgx_summary``, ``prs_summary``, ``apoe_genotype`` and optionally
        ``vo2max_tier`` / ``longevity_percentile``.
    snps_df : DataFrame
        Genotypes indexed by rsID (for optional APOE confirmation).

    Returns
    -------
    dict with ``status``, ``findings_with_economics`` (ranked by ROI),
    ``clinic_dashboard``, ``payer_impact``, ``high_confidence`` and
    ``disclaimer``.
    """
    findings = findings or {}
    econ_findings: List[Dict] = []
    try:
        econ_findings += _pgx_findings(findings.get("pgx_summary", {}))
        econ_findings += _prs_findings(findings.get("prs_summary", {}))
        econ_findings += _apoe_findings(findings.get("apoe_genotype"), snps_df)
        econ_findings += _exercise_longevity_findings(findings)
    except Exception as e:  # never raise from the pipeline
        return {"status": "error", "error": str(e),
                "findings_with_economics": [], "disclaimer": DISCLAIMER}

    # Rank by ROI (highest first); keep None ROIs last.
    econ_findings.sort(key=lambda f: (f.get("roi") or 0), reverse=True)

    result: Dict = {
        "status": "success" if econ_findings else "no_findings",
        "n_findings": len(econ_findings),
        "findings_with_economics": econ_findings,
        "high_confidence": [f for f in econ_findings if f["confidence"] == "high"],
        "disclaimer": DISCLAIMER,
    }
    result["clinic_dashboard"] = scale_to_clinic(result)
    result["payer_impact"] = scale_to_payer(result)
    return result
