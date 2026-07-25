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


# ══════════════════════════════════════════════════════════════════════════
# PERSONAL ECONOMIC-IMPACT ANALYSIS (standalone economic_analysis.html)
# ══════════════════════════════════════════════════════════════════════════
#
# Models the individual's 10-year modeled economic impact of ACTING on their
# results — expected medical-cost avoidance + monetised quality-of-life (QALY)
# gains, net of intervention cost — across genomic (PGx / carrier / PRS) and
# blood-work (PREVENT cardiovascular risk, prediabetes, biological age)
# findings. All figures are population-average model estimates, not a promise
# for any individual; the disclaimer is surfaced prominently in the report.

VALUE_PER_QALY = 100_000          # standard US cost-effectiveness threshold
PERSONAL_HORIZON_YEARS = 10
_ANALYSIS_COST = 700              # ~$200 genome + ~$500 lab panel, one-time

# Cost-of-illness / intervention assumptions (US, order-of-magnitude, 10-yr).
_MACE_COST = 85_000               # acute + 1-yr MI/stroke care
_MACE_QALY = 1.5
_STATIN_RRR = 0.27                # relative risk reduction, primary prevention
_STATIN_COST_10YR = 500
_T2D_COST = 85_000                # lifetime excess medical cost of T2D
_T2D_QALY = 2.0
_DPP_RRR = 0.58                   # lifestyle-program incidence reduction (DPP)
_DPP_COST_10YR = 3_000
_PREDIAB_PROGRESSION_10YR = 0.35  # prediabetes → diabetes over ~10 yr


def _money(x) -> str:
    x = round(x)
    return f"-${abs(x):,}" if x < 0 else f"${x:,}"


def analyze_personal_economics(economics_result: Optional[Dict] = None,
                               bloodwork_result: Optional[Dict] = None,
                               genetic_age_result: Optional[Dict] = None,
                               meta: Optional[Dict] = None) -> Dict:
    """Build the personal 10-year economic-impact model from the run's results."""
    items: List[Dict] = []

    def add(category, finding, avoided, qaly, intervention, confidence, basis):
        qv = qaly * VALUE_PER_QALY
        items.append({
            "category": category, "finding": finding,
            "avoided": round(avoided), "qaly": round(qaly, 2),
            "qaly_value": round(qv), "intervention": round(intervention),
            "net": round(avoided + qv - intervention),
            "confidence": confidence, "basis": basis,
        })

    # ── Genomic actionable findings (reuse the curated per-condition econ) ──
    if economics_result:
        for f in economics_result.get("findings_with_economics", []):
            outcome = f.get("outcome_value") or f.get("benefit") or 0
            prev = f.get("prevalence", 0.15)
            qaly = (f.get("qaly_gain") or f.get("qaly") or 0) * prev
            cost = f.get("cost", 200)
            avoided = outcome * prev * 0.30          # modeled realised benefit
            label = (f.get("clinical_benefit") or f.get("finding")
                     or f.get("drug") or "Genomic finding")
            if avoided <= 0 and qaly <= 0:
                continue
            add("Pharmacogenomic / genomic", label, avoided, qaly, cost,
                f.get("confidence", "moderate"),
                "Avoided adverse event × probability of relevant exposure (curated per-condition model).")

    # ── Blood-work derived ──
    adv = ((bloodwork_result or {}).get("clinical") or {}).get("advanced") or {}
    prevent = next((i["value"] for i in adv.get("indices", [])
                    if i.get("id") == "prevent_ascvd"), None)
    if prevent is not None:
        p = prevent / 100.0
        avoided = p * _MACE_COST * _STATIN_RRR
        qaly = p * _MACE_QALY * _STATIN_RRR
        add("Cardiovascular", f"10-yr ASCVD risk {prevent:.1f}% — risk-factor management",
            avoided, qaly, _STATIN_COST_10YR, "high" if prevent >= 7.5 else "moderate",
            f"PREVENT 10-yr ASCVD probability × ${_MACE_COST:,} event cost × {_STATIN_RRR:.0%} "
            "relative risk reduction from statin/lifestyle.")

    # Prediabetes → T2D (from biological-age input glucose or an HbA1c flag)
    bio = adv.get("biological_age") or {}
    glucose = (bio.get("inputs") or {}).get("glucose")
    a1c_flag = any("HbA1c" in (fl.get("name", "")) or "Glucose" in (fl.get("name", ""))
                   for fl in (((bloodwork_result or {}).get("clinical") or {}).get("flags") or []))
    prediabetic = (glucose is not None and 100 <= glucose < 126) or a1c_flag
    if prediabetic:
        avoided = _PREDIAB_PROGRESSION_10YR * _T2D_COST * _DPP_RRR
        qaly = _PREDIAB_PROGRESSION_10YR * _T2D_QALY * _DPP_RRR
        add("Metabolic", "Prediabetes flagged — intensive lifestyle (DPP)",
            avoided, qaly, _DPP_COST_10YR, "high",
            f"~{_PREDIAB_PROGRESSION_10YR:.0%} 10-yr progression risk × ${_T2D_COST:,} lifetime "
            f"T2D cost × {_DPP_RRR:.0%} reduction from a diabetes-prevention program.")

    # Biological-age "aging tax" (illustrative)
    accel = bio.get("accel")
    if accel is not None and abs(accel) >= 0.5:
        # ~$500/yr excess (or saved) healthcare cost per year of acceleration.
        val = -accel * 500 * PERSONAL_HORIZON_YEARS   # positive value if younger
        add("Biological aging", f"Biological age {accel:+.1f} yr vs chronological",
            max(0, val), 0.0, max(0, -val), "low",
            "Illustrative: accelerated biological age tracks higher healthcare "
            "utilisation; younger biological age tracks lower.")

    total_avoided = sum(i["avoided"] for i in items)
    total_qaly = sum(i["qaly"] for i in items)
    total_qaly_value = sum(i["qaly_value"] for i in items)
    total_intervention = sum(i["intervention"] for i in items)
    total_net = sum(i["net"] for i in items)
    gross_benefit = total_avoided + total_qaly_value
    roi = round((total_net) / _ANALYSIS_COST, 1) if _ANALYSIS_COST else None

    items.sort(key=lambda i: -i["net"])
    return {
        "available": bool(items),
        "n_items": len(items),
        "horizon_years": PERSONAL_HORIZON_YEARS,
        "items": items,
        "total_avoided": round(total_avoided),
        "total_qaly": round(total_qaly, 2),
        "total_qaly_value": round(total_qaly_value),
        "total_intervention": round(total_intervention),
        "total_net": round(total_net),
        "net_low": round(total_net * 0.5),          # ±50% sensitivity band
        "net_high": round(total_net * 1.5),
        "gross_benefit": round(gross_benefit),
        "top_preventable": (max(items, key=lambda i: i["avoided"])["finding"]
                            if items else None),
        "analysis_cost": _ANALYSIS_COST,
        "roi": roi,
        "value_per_qaly": VALUE_PER_QALY,
        "disclaimer": DISCLAIMER,
    }


def _esc_econ(s) -> str:
    s = "" if s is None else str(s)
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_economic_analysis_html(econ: Dict, file_label: str = "") -> str:
    """Standalone economic-impact sheet (economic_analysis.html)."""
    if not econ or not econ.get("available"):
        body = ("<p>No modeled economic-impact items for this profile — this needs "
                "actionable genomic findings and/or a supplied blood-work panel "
                "(<code>--bloodwork labs.json</code>).</p>")
    else:
        conf_col = {"high": "#1a7f37", "moderate": "#8a6100", "low": "#8a94a3"}
        rows = ""
        for i in econ["items"]:
            c = conf_col.get(i["confidence"], "#8a94a3")
            net_col = "#1a7f37" if i["net"] >= 0 else "#b3261e"
            rows += f"""
            <tr>
              <td><strong>{_esc_econ(i['finding'])}</strong>
                  <div style="color:#8a94a3;font-size:.82em">{_esc_econ(i['category'])} ·
                  <span style="color:{c}">{_esc_econ(i['confidence'])} confidence</span></div>
                  <div style="color:#9aa4b0;font-size:.78em">{_esc_econ(i['basis'])}</div></td>
              <td style="text-align:right;font-variant-numeric:tabular-nums">{_money(i['avoided'])}</td>
              <td style="text-align:right;font-variant-numeric:tabular-nums">{i['qaly']:g}<br>
                  <span style="color:#8a94a3;font-size:.8em">{_money(i['qaly_value'])}</span></td>
              <td style="text-align:right;font-variant-numeric:tabular-nums;color:#8a6100">{_money(i['intervention'])}</td>
              <td style="text-align:right;font-variant-numeric:tabular-nums;font-weight:700;color:{net_col}">{_money(i['net'])}</td>
            </tr>"""
        roi = econ.get("roi")
        net = econ["total_net"]
        net_col = "#1a7f37" if net >= 0 else "#b3261e"

        # SVG contribution chart — net value per finding (green +, red −)
        chart = ""
        ch_items = [i for i in econ["items"] if i["net"]]
        if ch_items:
            maxabs = max(abs(i["net"]) for i in ch_items) or 1
            barW, rowH, midX = 260, 26, 250
            svg_rows = ""
            for idx, i in enumerate(ch_items):
                y = idx * rowH + 4
                w = barW * abs(i["net"]) / maxabs
                col = "#2fae57" if i["net"] >= 0 else "#e0524a"
                x = midX if i["net"] >= 0 else midX - w
                lab = (i["finding"][:34] + "…") if len(i["finding"]) > 35 else i["finding"]
                svg_rows += (
                    f'<text x="{midX-8}" y="{y+15}" font-size="10.5" fill="#33404d" '
                    f'text-anchor="end">{_esc_econ(lab)}</text>'
                    f'<rect x="{x:.0f}" y="{y}" width="{w:.0f}" height="16" rx="3" fill="{col}"/>'
                    f'<text x="{(x+w+5) if i["net"]>=0 else (x-5):.0f}" y="{y+13}" font-size="10" '
                    f'fill="{col}" text-anchor="{"start" if i["net"]>=0 else "end"}">{_money(i["net"])}</text>')
            H = len(ch_items) * rowH + 10
            chart = (
                f'<div style="margin:14px 0"><div style="font-weight:700;color:#12467a;'
                f'margin-bottom:4px">Net value by finding</div>'
                f'<svg viewBox="0 0 560 {H}" width="100%" style="max-width:640px" '
                f'xmlns="http://www.w3.org/2000/svg">'
                f'<line x1="{midX}" y1="0" x2="{midX}" y2="{H}" stroke="#dbe3ec"/>{svg_rows}</svg></div>')

        rng = (f'<div style="color:#5b6673;font-size:.85em;margin-top:2px">'
               f'plausible range {_money(econ["net_low"])} – {_money(econ["net_high"])} (±50%)</div>')
        top_prev = (f'<div style="color:#8a94a3;font-size:.82em;margin-top:6px">'
                    f'Largest avoidable cost: <strong>{_esc_econ(econ.get("top_preventable",""))}</strong></div>'
                    if econ.get("top_preventable") else "")
        body = f"""
        {chart}
        <div style="display:flex;gap:18px;flex-wrap:wrap;align-items:center;
             background:linear-gradient(135deg,#f2f9f4,#eef4fb);border:1px solid #dbe8dd;
             border-radius:14px;padding:18px 22px;margin:14px 0">
          <div style="text-align:center;min-width:200px">
            <div style="font-size:.8em;color:#5b6673;text-transform:uppercase;letter-spacing:.04em">
              Modeled {econ['horizon_years']}-yr net value</div>
            <div style="font-size:2.8em;font-weight:800;color:{net_col}">{_money(net)}</div>
            <div style="color:#8a94a3;font-size:.82em">of acting on your results</div>
            {rng}
          </div>
          <div style="flex:1;min-width:240px;display:grid;grid-template-columns:1fr 1fr;gap:10px">
            <div><div style="font-size:.78em;color:#8a94a3">Medical cost avoided</div>
              <div style="font-size:1.4em;font-weight:700">{_money(econ['total_avoided'])}</div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Quality-of-life value</div>
              <div style="font-size:1.4em;font-weight:700">{_money(econ['total_qaly_value'])}
              <span style="font-size:.55em;color:#8a94a3">({econ['total_qaly']:g} QALYs)</span></div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Intervention cost</div>
              <div style="font-size:1.4em;font-weight:700;color:#8a6100">{_money(econ['total_intervention'])}</div></div>
            <div><div style="font-size:.78em;color:#8a94a3">ROI vs ~{_money(econ['analysis_cost'])} analysis</div>
              <div style="font-size:1.4em;font-weight:700;color:{net_col}">{roi}×</div></div>
          </div>
        </div>
        {top_prev}
        <table style="width:100%;border-collapse:collapse;font-size:.9em;margin-top:8px">
          <thead><tr style="background:#f7f9fb;text-align:right">
            <th style="text-align:left;padding:8px 10px">Finding</th>
            <th style="padding:8px 10px">Cost avoided</th>
            <th style="padding:8px 10px">QALY (value)</th>
            <th style="padding:8px 10px">Intervention</th>
            <th style="padding:8px 10px">Net value</th></tr></thead>
          <tbody>{rows}</tbody>
        </table>"""

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Economic-Impact Analysis{(' — ' + _esc_econ(file_label)) if file_label else ''}</title>
<style>
body{{font-family:-apple-system,"Segoe UI",Roboto,sans-serif;color:#1a1f26;max-width:1000px;
 margin:24px auto;padding:0 16px;line-height:1.5}}
h1{{font-size:1.7em;border-bottom:2px solid #12467a;padding-bottom:6px;color:#12467a}}
td,th{{padding:8px 10px;border-bottom:1px solid #eef1f4;vertical-align:top}}
.disc{{background:#fff8e6;border:1px solid #f0e0b0;border-radius:8px;padding:10px 14px;
 font-size:.85em;color:#6a5b2a;margin:12px 0}}
</style></head><body>
<h1>Economic-Impact Analysis</h1>
<p style="color:#667">A model of the 10-year economic impact of acting on your genomic and
blood-work results — expected medical-cost avoidance plus monetised quality-of-life gains,
net of intervention cost.</p>
<div class="disc">⚠️ <strong>Illustrative model, not financial or medical advice.</strong>
Figures use population-average cost-of-illness and risk-reduction estimates applied to your
findings; individual outcomes vary widely. Quality-of-life gains are monetised at the standard
${econ.get('value_per_qaly', 100000):,}/QALY threshold. Not a guarantee of savings.</div>
{body}
<p style="color:#9aa4b0;font-size:.8em;margin-top:16px">{_esc_econ(econ.get('disclaimer',''))}</p>
</body></html>"""
