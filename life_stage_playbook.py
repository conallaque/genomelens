"""
Life-Stage Playbook — decade-by-decade priorities from your genome
==================================================================

Synthesises every upstream result into a **decade-by-decade** priority list:
what to focus on in your 20s vs 30s vs 40s vs 50s vs 60s+, given your specific
genetic risk profile and (if known) your current labs and age.

This module deliberately runs LATE in the pipeline — after `holistic_synthesis`
— because its most important input is the **Genome Leverage Score**, which
itself aggregates APOE, PRS, longevity variants, immunogenetics, neurochemistry,
and PhenoAge. A decade playbook without that context would be generic.

Each decade starts from an evidence-based preventive-medicine baseline and is
then *modulated* by genome-specific inputs, with every genome-driven item
tagged to its source module so the reasoning is auditable.

Age handling
------------
Age is optional. If provided (via --age or --bloodwork labs), the current decade
is highlighted and framed as "you are here." If unknown, all decades are shown
without a highlight and a note explains why. There is no silent default — a
55-year-old is never told "your 20s" because age defaulted.
"""

from __future__ import annotations
from typing import Dict, List, Optional


def _get(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


# ── Baseline preventive-medicine priorities per decade ────────────────────────

_DECADE_BASE = {
    "20s": {
        "label": "Your 20s",
        "range": (18, 29),
        "theme": "Build the baseline. Habits set now compound for 60 years.",
        "base": [
            "Establish lifelong habits: regular exercise (strength + zone-2 cardio), "
            "whole-food diet, 7-9 h sleep, no smoking, moderate alcohol.",
            "Get a baseline lipid panel, fasting glucose/HbA1c, blood pressure, "
            "and vitamin D — you want a young-adult reference point to track from.",
            "Protect skin from UV; establish dental and vision care cadence.",
            "Build cardiovascular and bone-density capital while it's cheapest to add.",
        ],
    },
    "30s": {
        "label": "Your 30s",
        "range": (30, 39),
        "theme": "Hold the line as metabolism and time pressure rise.",
        "base": [
            "Re-check lipids, glucose/HbA1c, blood pressure every 1-3 years; watch "
            "for the first drift toward metabolic syndrome.",
            "Preserve muscle mass and VO2max — both start their natural decline; "
            "resistance training is the counter-lever.",
            "Manage stress and sleep deliberately as career/family load increases.",
            "If planning children, this is the window for carrier screening and "
            "fertility awareness.",
        ],
    },
    "40s": {
        "label": "Your 40s",
        "range": (40, 49),
        "theme": "Early-detection decade. Most preventable disease is silent now.",
        "base": [
            "Formal cardiovascular risk assessment (consider ApoB, Lp(a) once, "
            "and a coronary artery calcium score if risk is borderline).",
            "Begin colorectal cancer screening at 45 (earlier with family history).",
            "Baseline cognitive, hearing and vision checks; manage blood pressure "
            "tightly — midlife BP is a leading dementia risk factor.",
            "Prioritise sleep quality and muscle mass; both protect late-life "
            "independence.",
        ],
    },
    "50s": {
        "label": "Your 50s",
        "range": (50, 59),
        "theme": "Screening ramps up; protect the healthspan you've built.",
        "base": [
            "Full cancer-screening cadence (colorectal, plus sex-specific: "
            "mammography / prostate discussion).",
            "Aggressive cardiovascular and metabolic management — this is where "
            "risk curves steepen.",
            "Bone-density (DEXA) baseline, especially for women near menopause.",
            "Maintain strength training to fight sarcopenia; guard cognitive "
            "and cardiovascular reserve.",
        ],
    },
    "60s+": {
        "label": "Your 60s and beyond",
        "range": (60, 200),
        "theme": "Compress morbidity — spend the last years healthy, not declining.",
        "base": [
            "Fall prevention: strength, balance, vision, and medication review.",
            "Cognitive engagement, social connection, and hearing correction "
            "(untreated hearing loss accelerates cognitive decline).",
            "Continue cancer and cardiovascular screening per guidelines; "
            "vaccinations (influenza, pneumococcal, shingles, COVID boosters).",
            "Preserve muscle and protein intake to maintain independence.",
        ],
    },
}

_DECADE_ORDER = ["20s", "30s", "40s", "50s", "60s+"]


def resolve_age(explicit_age: Optional[int],
                bloodwork_result: Optional[Dict]) -> Optional[int]:
    """Resolve the subject's age for the playbook. Priority: an explicit age
    (from --age), else the age nested in a compare_bloodwork result at
    ``result["clinical"]["age_used"]`` (its true location — the top level does
    NOT carry it), else None. Coerced to int."""
    if explicit_age is not None:
        return int(explicit_age)
    if bloodwork_result:
        a = (bloodwork_result.get("clinical") or {}).get("age_used")
        if a is not None:
            return int(a)
    return None


def _decade_for_age(age: Optional[int]) -> Optional[str]:
    if age is None:
        return None
    for key in _DECADE_ORDER:
        lo, hi = _DECADE_BASE[key]["range"]
        if lo <= age <= hi:
            return key
    return "60s+" if age >= 60 else "20s"


def analyze_life_stage_playbook(
    age: Optional[int] = None,
    holistic_synthesis_result: Optional[Dict] = None,
    immunogenetics_result: Optional[Dict] = None,
    addiction_genetics_result: Optional[Dict] = None,
    neurochemistry_result: Optional[Dict] = None,
    family_planning_result: Optional[Dict] = None,
    tier1_results: Optional[Dict] = None,
) -> Dict:
    """Build the decade-by-decade playbook, modulated by genome context."""

    # ── Gather genome-specific modulators ────
    leverage = _get(holistic_synthesis_result, "genome_leverage", default={}) or {}
    leverage_tier = leverage.get("tier")
    leverage_score = leverage.get("score")

    apoe = (tier1_results or {}).get("apoe_genotype", "") if tier1_results else ""
    apoe_e4 = "4" in (apoe or "")

    # smoking / CHRNA5 flag
    chrna5_flag = False
    for f in _get(addiction_genetics_result, "composite", "clinical_flags", default=[]) or []:
        if "smok" in (f.get("title", "") + f.get("text", "")).lower():
            chrna5_flag = True

    # neuroplasticity window
    nc = _get(neurochemistry_result, "composite", default={}) or {}
    high_plasticity = str(nc.get("bdnf_class", "")).startswith("Val/Val")

    # immunogenetics susceptibilities (e.g. severe flu) & headline resistances
    immuno_susc = [f for f in _get(immunogenetics_result, "findings", default=[]) or []
                   if f.get("impact") == "susceptible"]
    flu_risk = any("influenza" in f.get("name", "").lower() or "IFITM3" in f.get("gene", "")
                   for f in immuno_susc)

    # holistic insights (early metabolic signals etc.)
    insight_ids = {i.get("id") for i in _get(holistic_synthesis_result, "insights", default=[]) or []}
    early_metabolic = ("glucose_hba1c_stress_discordance" in insight_ids
                       or "apoe_e4_lipid_amplification" in insight_ids)

    # family planning relevance
    has_repro_findings = bool(_get(family_planning_result, "n_recessive", default=0)
                              or _get(family_planning_result, "n_dominant", default=0))

    # ── Build decades ────
    current = _decade_for_age(age)
    decades: List[Dict] = []
    for key in _DECADE_ORDER:
        base = _DECADE_BASE[key]
        genome_items: List[Dict] = []

        def add(text, source):
            genome_items.append({"text": text, "source": source})

        # Leverage framing appears in every decade at/after the current one
        if leverage_tier in ("Very favorable", "Favorable") and key in ("20s", "30s"):
            add(f"Your Genome Leverage Score ({leverage_score}/100, "
                f"'{leverage_tier}') means environment dominates your trajectory "
                "— early investment compounds unusually well. Treat health as a "
                "compounding asset now.", "holistic_synthesis")
        elif leverage_tier == "Actionable risk" and key in ("20s", "30s", "40s"):
            add("Your genome carries specific risk anchors — starting targeted "
                "prevention and monitoring earlier than average has strong "
                "evidence for you.", "holistic_synthesis")

        # CHRNA5 / smoking — heaviest weight in 20s
        if chrna5_flag and key == "20s":
            add("Never start smoking or vaping — you carry the CHRNA5 risk allele "
                "that makes smoking both more addictive and more dangerous. This "
                "is your single highest-value behavioural decision.", "addiction_genetics")

        # Plasticity window
        if high_plasticity and key in ("20s", "30s"):
            add("Peak neuroplasticity window: your BDNF Val/Val genotype makes "
                "deliberate skill-building compound faster than average. Pick "
                "1-2 domains and invest the reps now.", "neurochemistry")

        # APOE e4 — cardio/cognitive prevention early
        if apoe_e4 and key in ("30s", "40s", "50s"):
            add("APOE ε4 carrier: begin lipid, blood-pressure and cognitive-"
                "reserve protection earlier and more aggressively — the benefit "
                "of midlife control is larger for you.", "APOE")

        # Early metabolic signals
        if early_metabolic and key in ("20s", "30s"):
            add("Early metabolic signals were flagged in your labs — address diet, "
                "activity and sleep now, and re-test; these are far cheaper to "
                "reverse in this decade than later.", "holistic_synthesis / bloodwork")

        # Flu susceptibility
        if flu_risk and key in ("20s", "30s", "40s", "50s", "60s+"):
            add("Annual influenza vaccination matters more for you — your IFITM3 "
                "genotype is associated with more severe influenza.", "immunogenetics")

        # Family planning window
        if has_repro_findings and key in ("20s", "30s"):
            add("You carry reproductive-relevant variants — if planning children, "
                "do partner carrier screening and review the Family Planning "
                "section in this window.", "family_planning")
        elif key in ("30s",):
            add("If planning children, fertility naturally declines through this "
                "decade — factor timing and consider carrier screening.", "general")

        decades.append({
            "key": key,
            "label": base["label"],
            "theme": base["theme"],
            "range": base["range"],
            "base": base["base"],
            "genome_items": genome_items,
            "is_current": key == current,
        })

    return {
        "available": True,
        "age": age,
        "current_decade": current,
        "age_known": age is not None,
        "leverage_tier": leverage_tier,
        "leverage_score": leverage_score,
        "decades": decades,
        "note": (
            f"You are currently in {_DECADE_BASE[current]['label'].lower()} — "
            "that decade is highlighted below."
            if current else
            "Age was not provided (pass --age or --bloodwork labs), so all "
            "decades are shown without a 'you are here' highlight."),
    }
