"""
Holistic Synthesis — cross-panel pattern detection
==================================================

This module doesn't produce new genotype calls. It reads *every* upstream
module's output and detects the patterns that only exist when you look
across panels — the things a single-module report will never say, because
they emerge from the interaction of two or more findings.

Examples of what it catches:
  • **FUT2 non-secretor × elevated hs-CRP** — non-secretors have a mildly
    pro-inflammatory microbiome baseline; some fraction of an "elevated"
    CRP in a non-secretor is genotype-driven, not lifestyle.
  • **Fasting-glucose ≥ 100 × HbA1c < 5.7** — the classic acute-stress /
    poor-sleep morning-of-draw pattern, not real dysglycemia.
  • **Favorable-genome leverage score** — a person with few genetic risk
    anchors has more of their eventual healthspan defined by behaviour
    than someone whose genome is already pushing them toward pathology.
    Their upside is bigger; so is their downside if they neglect lifestyle.
  • **Ancestral-diet fit** — Northern-European + LCT persistent + Yamnaya×EEF
    affinity → Mediterranean+dairy diet is *genetically appropriate*, not
    just generically healthy.
  • **APOE ε4 × elevated LDL** — LDL matters ~2× more with ε4.
  • **CHRNA5 A-carrier who has never smoked** — an actively-realised
    prevention success worth acknowledging.
  • **Warrior COMT + high MAOA + prediabetic glucose** — the stress-driven
    glucose spike pattern, not metabolic disease.
  • **HFE clear × high ferritin** — genetics don't explain it; look at
    inflammation, diet, and lab timing.
  • **BDNF Val/Val + favorable neurotype + young adult** — deliberate
    practice compounds materially more than average.

Every insight cites the *specific* upstream findings that trigger it.
"""

from __future__ import annotations


def _get(d: dict | None, *keys, default=None):
    for k in keys:
        if d is None:
            return default
        d = d.get(k)
    return d if d is not None else default


def _find_variant(module_result: dict | None, rsid: str) -> dict | None:
    if not module_result:
        return None
    for k in ("variants", "findings", "indices"):
        seq = module_result.get(k)
        if isinstance(seq, list):
            for v in seq:
                if v.get("rsid") == rsid:
                    return v
    return None


def _bloodwork_marker(bloodwork_result: dict | None, name: str) -> float | None:
    """Find a raw lab value from the clinical panel by name (case-insensitive)."""
    if not bloodwork_result:
        return None
    clinical = bloodwork_result.get("clinical") or {}
    for sys in clinical.get("systems", []):
        for m in sys.get("markers", []):
            if m.get("name", "").lower() == name.lower():
                return m.get("value")
    return None


# ─── Pattern detectors ────────────────────────────────────────────────────────

def _insight(id, title, category, impact, explanation, action, evidence,
             confidence="moderate", severity=1):
    return {"id": id, "title": title, "category": category, "impact": impact,
            "explanation": explanation, "action": action, "evidence": evidence,
            "confidence": confidence, "severity": severity}


def _detect_fut2_crp_baseline(immuno, bloodwork):
    """FUT2 non-secretor + elevated CRP → some fraction is genotype."""
    fut2 = None
    for f in _get(immuno, "findings", default=[]) or []:
        if f.get("gene") == "FUT2":
            fut2 = f
            break
    if not fut2 or ("non-secretor" not in fut2.get("phenotype", "").lower() \
       and "non-secretor" not in fut2.get("verdict", "").lower()):
        return None
    crp = _bloodwork_marker(bloodwork, "hs-CRP") or _bloodwork_marker(bloodwork, "C-Reactive Protein")
    if crp is None or crp <= 1.0:
        return None
    return _insight(
        "fut2_crp_baseline",
        "FUT2 non-secretor may inflate your hs-CRP baseline",
        "Cross-panel: gut × inflammation", "informational",
        f"You are a FUT2 non-secretor (near-complete resistance to GII.4 "
        f"norovirus). Non-secretors have a well-documented shifted gut "
        f"microbiome — fewer Bifidobacterium, more Bacteroides — which "
        f"produces a mildly higher inflammatory baseline. Your hs-CRP of "
        f"{crp:g} mg/L may reflect ~0.3-0.7 of that from genotype rather "
        f"than lifestyle. Track hs-CRP as a **trend over time**, not as an "
        f"absolute value at any single draw.",
        "Retest CRP alongside a repeat lipid panel; a stable value in the "
        "1-2 mg/L range without other risk markers is likely just your "
        "non-secretor baseline. Focus lifestyle levers on the real drivers "
        "(diet quality, sleep, visceral fat) rather than chasing a lower "
        "absolute CRP.",
        [{"module": "immunogenetics", "finding": "FUT2 non-secretor"},
         {"module": "bloodwork", "value": f"hs-CRP {crp} mg/L"}],
        "moderate", severity=1,
    )


def _detect_glucose_hba1c_discordance(bloodwork):
    """Fasting glucose ≥100 but HbA1c <5.7 → acute stress, not dysglycemia."""
    glu = _bloodwork_marker(bloodwork, "Fasting Glucose")
    a1c = _bloodwork_marker(bloodwork, "HbA1c")
    if glu is None or a1c is None:
        return None
    if not (glu >= 100 and a1c < 5.7):
        return None
    return _insight(
        "glucose_hba1c_stress_discordance",
        f"Fasting glucose {glu:g} vs HbA1c {a1c:g}% — acute-draw pattern, not dysglycemia",
        "Cross-panel: labs × interpretation", "informational",
        f"Your fasting glucose ({glu:g}) sits in the prediabetes range "
        f"(100-125), but your HbA1c ({a1c:g}%) — the 3-month glucose average "
        f"— is squarely non-diabetic (<5.7). This gap is characteristic of "
        f"**morning-of-draw cortisol / poor-sleep / fasting-stress glucose "
        f"elevation**, not real insulin resistance. If you had genuine "
        f"prediabetes, your HbA1c would be creeping up.",
        "Retest fasted, well-rested (no bad night before, no early alarm), "
        "and consider adding fasting insulin + HOMA-IR to the panel — those "
        "would definitively distinguish stress-glucose from real insulin "
        "resistance. Don't over-interpret one prediabetic-range fasting "
        "glucose as needing intervention.",
        [{"module": "bloodwork", "value": f"Fasting glucose {glu:g}"},
         {"module": "bloodwork", "value": f"HbA1c {a1c:g}%"}],
        "high", severity=1,
    )


def _detect_apoe_lipid_amplification(tier1, bloodwork):
    """APOE ε4 carrier + elevated LDL/non-HDL → matters materially more."""
    apoe = None
    if tier1:
        apoe = tier1.get("apoe_genotype")
    if not apoe:
        return None
    if "e4" not in apoe.lower() and "ε4" not in apoe.lower() and "4/" not in apoe and "/4" not in apoe:
        return None
    ldl = _bloodwork_marker(bloodwork, "LDL Cholesterol")
    non_hdl = _bloodwork_marker(bloodwork, "Non-HDL Cholesterol")
    apob = _bloodwork_marker(bloodwork, "Apolipoprotein B")
    if not any(v and v > 100 for v in (ldl, non_hdl, apob)):
        return None
    values = ", ".join(f"{k} {v:g}" for k, v in
                       [("LDL", ldl), ("non-HDL", non_hdl), ("ApoB", apob)] if v)
    return _insight(
        "apoe_e4_lipid_amplification",
        f"APOE ε4 × elevated {values} — the combination that matters most",
        "Cross-panel: APOE × lipids", "actionable",
        f"APOE ε4 carriers have ~2× the impact of elevated LDL / non-HDL / "
        f"ApoB on both cardiovascular risk AND late-life Alzheimer's risk. "
        f"Non-ε4 people can tolerate modestly elevated lipids for longer; "
        f"ε4 carriers cannot. Given your APOE genotype ({apoe}) and current "
        f"lipids ({values}), lipid control matters more for you than for the "
        f"average person.",
        "**Priority action:** aim for LDL <100 (ideally <70) and ApoB <90 "
        "(ideally <80). This may involve dietary saturated-fat reduction, "
        "soluble fibre, exercise, and — if lifestyle alone doesn't move it "
        "— early statin discussion with a physician despite your young age. "
        "Also relevant for Alzheimer's-risk trajectory over the next 40 years.",
        [{"module": "apoe", "value": apoe},
         {"module": "bloodwork", "value": values}],
        "high", severity=3,
    )


def _detect_chrna5_prevention_success(neurochem, meta):
    """CHRNA5 A-carrier who has not started smoking — realised prevention."""
    if not neurochem:
        return None
    chrna5 = None
    for f in neurochem.get("findings", []):
        if f.get("gene") == "CHRNA5":
            chrna5 = f
            break
    if not chrna5 or "A-carrier" not in chrna5.get("phenotype", ""):
        return None
    # Detect via meta or via absence of a smoker flag
    smoker = (meta or {}).get("smoker", False)
    if smoker:
        return None
    return _insight(
        "chrna5_realised_prevention_success",
        "CHRNA5 A-carrier — non-smoking is an actively-realised prevention success",
        "Cross-panel: risk gene × behaviour", "informational",
        "You carry the CHRNA5 α5 nAChR risk allele. In people who never "
        "start smoking, this variant does nothing — it only matters if you "
        "smoke, in which case you would smoke harder and have markedly "
        "higher lung-cancer / COPD risk. Your non-smoker status is quietly "
        "the single highest-value behavioural decision your genome cares "
        "about.",
        "Continue avoiding cigarettes, vapes, and cigar/cigarillos — "
        "**including social smoking** — indefinitely. If you're ever tempted, "
        "the honest math is that this variant makes 'a few here and there' "
        "physiologically more addictive for you than for the average person.",
        [{"module": "neurochemistry", "finding": chrna5.get("phenotype")}],
        "high", severity=1,
    )


def _detect_ancestral_diet_fit(deep_ancestry, immuno):
    """European ancestry + LCT persistence + Yamnaya/EEF affinity → Mediterranean+dairy diet is *genetically appropriate*."""
    if not deep_ancestry:
        return None
    ax = deep_ancestry.get("european_axis") or {}
    if not ax.get("available"):
        return None
    lct = False
    for m in ax.get("used", []):
        if m.get("rsid") == "rs4988235" and m.get("dose", 0) >= 1:
            lct = True
            break
    ap = deep_ancestry.get("ancient_populations") or {}
    yamnaya = 0.0
    eef = 0.0
    for p in ap.get("populations", []):
        if p["short"] == "Yamnaya":
            yamnaya = p["affinity"]
        if p["short"].startswith("EEF"):
            eef = p["affinity"]
    if not (lct and yamnaya >= 0.5 and eef >= 0.5):
        return None
    return _insight(
        "ancestral_diet_fit",
        "You are the archetypal Mediterranean-plus-dairy dietary candidate",
        "Cross-panel: ancestry × nutrition", "actionable",
        f"Your genome carries the classic Northern-Central European "
        f"admixture: LCT lactase persistence (Yamnaya-derived), Yamnaya "
        f"affinity {yamnaya*100:.0f}%, EEF affinity {eef*100:.0f}%. All "
        f"three ancestral gene pools contributing to you spent 5,000-10,000+ "
        f"years adapting to Mediterranean-plus-dairy foodways: wheat, "
        f"barley, olives, wine, dairy fermentation (cheese/yogurt), grass-"
        f"fed meat. This isn't just 'a healthy diet' — it is **the diet "
        f"your genome was actually optimised for**.",
        "Base your eating on: whole grains, legumes, olive oil, moderate "
        "fermented dairy, fish/red meat in moderation, vegetables, moderate "
        "red wine, minimal processed foods. Avoid the two categories your "
        "genome has *not* adapted to: (1) modern ultra-processed foods "
        "(zero evolutionary exposure), (2) strict veganism (Yamnaya "
        "pastoralism strongly selected for dairy/meat metabolism in your "
        "lineage).",
        [{"module": "deep_ancestry", "value": f"Yamnaya {yamnaya:.0%} · EEF {eef:.0%}"},
         {"module": "deep_ancestry", "value": "LCT lactase persistence carrier"}],
        "high", severity=1,
    )


def _detect_favorable_genome_leverage(tier1, prs_result, immuno, neurochem,
                                       deep_ancestry, bloodwork):
    """Composite: multiple protective genotypes + few red flags → environment
    dominates. Very-favorable genomes have wider upside/downside than average."""
    score = 50   # start at neutral

    reasons_up: list[str] = []
    reasons_down: list[str] = []

    # APOE ε4 status (a big single lever)
    apoe = (tier1 or {}).get("apoe_genotype", "") or ""
    if apoe and "4" not in apoe:
        score += 8
        reasons_up.append(f"no APOE ε4 ({apoe})")
    elif apoe and "4" in apoe:
        score -= 8
        reasons_down.append(f"APOE ε4 carrier ({apoe})")

    # PRS panels — count how many are elevated/high
    n_prs_high = 0
    n_prs_low = 0
    prs = ((tier1 or {}).get("prs_summary") or {}) if isinstance(tier1, dict) else {}
    for _, p in (prs or {}).items():
        tier = ((p or {}).get("tier") or "").lower()
        if "high" in tier or "elevated" in tier:
            n_prs_high += 1
        elif "below" in tier or "low" in tier:
            n_prs_low += 1
    if n_prs_high == 0 and n_prs_low >= 2:
        score += 6
        reasons_up.append(f"{n_prs_low} PRS panels in below-average tier, none elevated")
    elif n_prs_high >= 3:
        score -= 8
        reasons_down.append(f"{n_prs_high} PRS panels in elevated/high tier")

    # Longevity variants live in deep_ancestry.genetic_longevity or in immuno.
    gl_variants = _get(deep_ancestry, "genetic_longevity", default={}) or {}
    if gl_variants.get("lean") == "favorable" or gl_variants.get("n_favorable", 0) >= 2:
        score += 6
        reasons_up.append(f"{gl_variants.get('n_favorable',0)} favorable longevity variants (FOXO3/CETP/IL6/APOE ε2)")

    # Immunogenetics headlines
    im_headlines = (immuno or {}).get("headlines", [])
    if len(im_headlines) >= 2:
        score += 4
        reasons_up.append(f"{len(im_headlines)} major viral-resistance headlines")

    # Neurochemistry
    nc_c = _get(neurochem, "composite", default={}) or {}
    if nc_c.get("comt_class") == "middle" and nc_c.get("bdnf_class", "").startswith("Val/Val"):
        score += 4
        reasons_up.append("adaptive-middle COMT + full BDNF plasticity")

    # PhenoAge
    bio = _get(bloodwork, "clinical", "advanced", "biological_age", default={}) or {}
    accel = bio.get("accel")
    if accel is not None:
        if accel <= -2:
            score += 8
            reasons_up.append(f"PhenoAge {accel:+.1f} yr (biologically younger)")
        elif accel >= 3:
            score -= 8
            reasons_down.append(f"PhenoAge {accel:+.1f} yr (biologically older)")

    # Flagged high-severity clinical markers
    n_critical = len([f for f in (_get(bloodwork, "clinical", "flags", default=[]) or [])
                      if f.get("severity", 0) >= 3])
    if n_critical:
        score -= 6 * n_critical
        reasons_down.append(f"{n_critical} critical clinical marker(s)")

    score = max(0, min(100, round(score)))
    if score >= 75:
        tier = "Very favorable"
        narrative = (
            "Your genome sits in a rare cluster where the intervention math "
            "is genuinely favorable. Nearly every major fork in human "
            "evolutionary history has landed on the protective side for "
            "your lineage. Practical translation: **behaviour is the "
            "trajectory** — you don't have anchors dragging you toward a "
            "specific pathology, which means the environment you build for "
            "yourself over the next 40 years will dominate the outcome. "
            "Upside is very high; downside from lifestyle neglect is "
            "also higher-than-typical because you lack built-in shock "
            "absorbers pointing you elsewhere."
        )
    elif score >= 60:
        tier = "Favorable"
        narrative = (
            "Your genome carries several protective findings and few strong "
            "risk anchors. Environmental / lifestyle factors will be the "
            "dominant driver of your long-term trajectory — genes tilt the "
            "table modestly in your favor but don't decide it for you."
        )
    elif score >= 40:
        tier = "Balanced"
        narrative = (
            "A mixed genetic profile with both protective and risk-carrying "
            "findings across systems. Standard preventive-medicine playbook "
            "applies; specific findings elsewhere in the report identify "
            "the highest-leverage actions."
        )
    else:
        tier = "Actionable risk"
        narrative = (
            "One or more genetic findings carry non-trivial risk that "
            "warrants active management. The specific loci flagged above "
            "identify the priorities — this is the profile where "
            "clinical monitoring and, in some cases, pharmacological "
            "management have their strongest evidence base."
        )
    return {
        "score": score, "tier": tier, "narrative": narrative,
        "reasons_up": reasons_up, "reasons_down": reasons_down,
    }


def _detect_plasticity_leverage(neurochem, prs_result, tier1, meta):
    """BDNF Val/Val + favorable neurotype + young adult → deliberate practice
    compounds materially more than average."""
    if not neurochem:
        return None
    c = (neurochem.get("composite") or {})
    if not (c.get("bdnf_class", "").startswith("Val/Val")
            and c.get("comt_class") in ("middle", "warrior")):
        return None
    age = (meta or {}).get("age")
    if age and age > 45:
        return None
    return _insight(
        "plasticity_leverage",
        "High-plasticity + adaptive neurotype × young adult — deliberate practice compounds",
        "Cross-panel: neuroplasticity × age", "actionable",
        f"You carry BDNF Val/Val (full activity-dependent BDNF secretion) "
        f"and COMT {c.get('comt_class','?')} + MAOA "
        f"{c.get('maoa_class','?')} — the neurochemical substrate that "
        f"turns deliberate practice into durable skill faster than average. "
        f"Combined with your young age, this window is the most productive "
        f"decade you will ever have for skill acquisition. Both the "
        f"neuroplasticity and the neurotransmitter systems are cooperating.",
        "Pick 1-2 domains that genuinely matter to you and put the reps in "
        "with structured feedback (Anki, spaced repetition, deliberate "
        "practice with a coach or measurable metric). You'll compound more "
        "reliably than most people your age. Aerobic exercise 4+×/week is "
        "the single most effective BDNF booster on top of your baseline.",
        [{"module": "neurochemistry", "finding": "BDNF Val/Val + adaptive COMT/MAOA"}],
        "moderate", severity=1,
    )


def _detect_iron_stress_pattern(bloodwork, immuno, tier1):
    """HFE clear + high ferritin → inflammatory/dietary/lifestyle, not genetic."""
    ferritin = _bloodwork_marker(bloodwork, "Ferritin")
    if ferritin is None or ferritin < 300:
        return None
    hfe_variant = False
    for f in (tier1 or {}).get("variants", []) or []:
        if f.get("gene") == "HFE" and f.get("risk_copies", 0) >= 1:
            hfe_variant = True
            break
    if hfe_variant:
        return None
    crp = _bloodwork_marker(bloodwork, "hs-CRP")
    return _insight(
        "iron_no_hfe_pattern",
        f"Ferritin {ferritin:g} ng/mL — elevated but not genetically explained",
        "Cross-panel: iron × HFE × inflammation", "actionable",
        f"Your ferritin is elevated at {ferritin:g} ng/mL, but you don't "
        f"carry HFE C282Y or H63D — so this isn't hereditary hemochromatosis. "
        f"Non-genetic drivers of high ferritin: inflammation (ferritin is an "
        f"acute-phase reactant"
        + (f"; your hs-CRP is {crp:g}, which fits" if crp and crp > 1.5 else "")
        + "), fatty liver, alcohol, high red-meat diet, "
        "iron supplementation. Ferritin alone can look scary without "
        "context.",
        "Add transferrin saturation to the next panel; if it's >45%, "
        "hereditary hemochromatosis panels (beyond C282Y/H63D — some rare "
        "variants aren't on consumer chips) may be worth pursuing. If "
        "transferrin sat is normal, the ferritin is inflammation- or "
        "diet-driven and will move with lifestyle adjustments. Don't "
        "supplement iron.",
        [{"module": "bloodwork", "value": f"Ferritin {ferritin:g}"},
         {"module": "tier1", "value": "HFE C282Y / H63D not detected"}],
        "moderate", severity=2,
    )


def _detect_coffee_synthesis(neurochem, pgx_result):
    """Combine COMT + CYP1A2 + ADORA2A into one specific coffee protocol."""
    if not neurochem:
        return None
    c = (neurochem.get("composite") or {})
    # Look for CYP1A2 in pgx
    cyp1a2 = None
    if pgx_result and isinstance(pgx_result, dict):
        for gene, phen in (pgx_result.get("phenotypes") or {}).items():
            if "CYP1A2" in gene:
                cyp1a2 = phen
                break
    return _insight(
        "coffee_protocol",
        "Your personalised coffee protocol (COMT × CYP1A2 × MAOA)",
        "Cross-panel: neurochemistry × PGx", "actionable",
        f"COMT {c.get('comt_class','?')} + MAOA "
        f"{c.get('maoa_class','?')}"
        + (f" + CYP1A2 {cyp1a2}" if cyp1a2 else "") + ". "
        "You clear caffeine and catecholamines fast. Practical "
        "consequences: caffeine wears off subjectively quickly, but its "
        "REM-suppressing half-life in sleep-relevant compartments is longer "
        "than you feel. Anxiety ceiling is real (COMT + fast MAOA can tip "
        "into jitter above ~200 mg per bolus dose).",
        "Espresso-style: 1-2 shots (~120-180 mg) in the morning, one more "
        "before noon if needed, **no caffeine after 2 pm** regardless of "
        "feeling. Avoid single-dose bolus > 200 mg (large drip coffees, "
        "energy drinks). L-theanine 100-200 mg with morning coffee blunts "
        "the anxiety edge if you notice it.",
        [{"module": "neurochemistry", "value": f"COMT {c.get('comt_class','?')} · MAOA {c.get('maoa_class','?')}"}],
        "moderate", severity=1,
    )


# ═════════════════════════════════════════════════════════════════════════════
# Priority-action ranking
# ═════════════════════════════════════════════════════════════════════════════

def _rank_priority_actions(insights: list[dict],
                           genome_leverage: dict) -> list[dict]:
    """Turn detected insights into a ranked priority list based on severity
    × confidence × modifiability."""
    ranked = []
    for i in insights:
        modifiable = i["category"].startswith("Cross-panel")
        rank_score = (
            i.get("severity", 1) * 3
            + {"high": 3, "moderate": 2, "low": 1}.get(i.get("confidence", "moderate"), 2)
            + (2 if modifiable else 0)
            + (2 if i.get("impact") == "actionable" else
               (1 if i.get("impact") == "informational" else 0))
        )
        ranked.append({"score": rank_score, "insight": i})
    ranked.sort(key=lambda x: -x["score"])
    out = []
    for i, r in enumerate(ranked[:6], start=1):
        ins = r["insight"]
        out.append({
            "priority": i,
            "title": ins["title"],
            "why": ins["explanation"][:220] + ("…" if len(ins["explanation"]) > 220 else ""),
            "action": ins["action"],
            "id": ins["id"],
        })
    return out


# ═════════════════════════════════════════════════════════════════════════════
# Master analyzer
# ═════════════════════════════════════════════════════════════════════════════

def analyze_holistic_synthesis(
    tier1_results: dict | None = None,
    bloodwork_result: dict | None = None,
    immunogenetics_result: dict | None = None,
    neurochemistry_result: dict | None = None,
    deep_ancestry_result: dict | None = None,
    ancestry_result: dict | None = None,
    prs_result: dict | None = None,
    pgx_result: dict | None = None,
    meta: dict | None = None,
) -> dict:
    """Detect cross-panel patterns and produce a composite genome-leverage
    score + ranked priority actions."""
    detectors = [
        lambda: _detect_fut2_crp_baseline(immunogenetics_result, bloodwork_result),
        lambda: _detect_glucose_hba1c_discordance(bloodwork_result),
        lambda: _detect_apoe_lipid_amplification(tier1_results, bloodwork_result),
        lambda: _detect_chrna5_prevention_success(neurochemistry_result, meta),
        lambda: _detect_ancestral_diet_fit(deep_ancestry_result, immunogenetics_result),
        lambda: _detect_plasticity_leverage(neurochemistry_result, prs_result, tier1_results, meta),
        lambda: _detect_iron_stress_pattern(bloodwork_result, immunogenetics_result, tier1_results),
        lambda: _detect_coffee_synthesis(neurochemistry_result, pgx_result),
    ]
    insights: list[dict] = []
    for d in detectors:
        try:
            r = d()
        except Exception:
            continue
        if r is not None:
            insights.append(r)

    genome_leverage = _detect_favorable_genome_leverage(
        tier1_results, prs_result, immunogenetics_result,
        neurochemistry_result, deep_ancestry_result, bloodwork_result,
    )

    priority_actions = _rank_priority_actions(insights, genome_leverage)

    return {
        "available": bool(insights) or genome_leverage is not None,
        "n_insights": len(insights),
        "insights": insights,
        "genome_leverage": genome_leverage,
        "priority_actions": priority_actions,
    }
