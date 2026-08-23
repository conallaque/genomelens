"""
Longevity Composite & Integrated Year-Long Plan
================================================

Cross-cutting synthesis layer that combines nutrition + exercise analyses into:

  • Longevity composite score (0-100) with component breakdown and biggest levers
  • Year-long integrated periodisation (4 mesocycles × nutrition + exercise)
  • Body-composition trajectory model (12-week realistic-outcome estimator)
  • Print-friendly executive summary (one-page synthesis)
  • JSON export envelope
"""

from __future__ import annotations


def _safe_pct(d: dict, key: str, default: float = 50.0) -> float:
    s = d.get(key) if d else None
    if not s:
        return default
    p = s.get("percentile")
    return float(p) if p is not None else default


def _inv(p: float) -> float:
    """Invert a percentile (higher = better) for risk-direction traits."""
    return 100 - p


# ── Longevity Composite Score ──────────────────────────────────────────────

def longevity_composite(nutrition: dict, exercise: dict) -> dict:
    pgs = nutrition.get("polygenic_scores", {})
    inflam = nutrition.get("inflammation", {})
    sat_fat = nutrition.get("saturated_fat", {})
    profile = exercise.get("composite_profile", {})
    injury_map = exercise.get("injury_risk_map", {})

    components: list[dict] = []

    # Cardiometabolic — invert LDL, TG, BP, T2D so higher score = better
    cmd_components = []
    for k in ("LDL_cholesterol", "Triglycerides", "BP_systolic", "T2D"):
        cmd_components.append(_inv(_safe_pct(pgs, k)))
    hdl = _safe_pct(pgs, "HDL_cholesterol")  # higher percentile here = lower HDL (we inverted in dashboard logic)
    cmd_components.append(_inv(hdl))
    cardiometabolic = sum(cmd_components) / len(cmd_components)
    components.append({
        "component": "Cardiometabolic genetic baseline",
        "score": round(cardiometabolic, 1),
        "weight": 0.25,
    })

    # Inflammatory tone (score 0-3+ from inflammation index; invert)
    inflam_score = inflam.get("score", 0) if inflam else 0
    inflam_norm = max(20, 100 - 20 * inflam_score)
    components.append({"component": "Inflammatory tone", "score": round(inflam_norm, 1), "weight": 0.15})

    # APOE / brain — ε4 carrier penalty
    apoe_geno = sat_fat.get("apoe_genotype", "") if sat_fat else ""
    if apoe_geno.startswith("ε4/ε4"):
        apoe_score = 30
    elif apoe_geno.startswith("ε4"):
        apoe_score = 55
    elif apoe_geno.startswith("ε2"):
        apoe_score = 90
    else:
        apoe_score = 75
    components.append({"component": "APOE / cognitive ageing", "score": apoe_score, "weight": 0.10})

    # Athletic capacity (predicts all-cause mortality strongly)
    athl = profile.get("overall_score", 50)
    components.append({"component": "Athletic capacity index", "score": round(athl, 1), "weight": 0.20})

    # Injury / musculoskeletal resilience
    injury_idx = injury_map.get("overall_index", 30) if injury_map else 30
    musculo = max(20, 100 - injury_idx)
    components.append({"component": "Musculoskeletal resilience", "score": round(musculo, 1), "weight": 0.10})

    # Recovery (recovery speed)
    rec = exercise.get("recovery", {}).get("speed", "Moderate")
    rec_score = {"Fast": 85, "Moderate": 60, "Slow": 35}.get(rec, 50)
    components.append({"component": "Recovery capacity", "score": rec_score, "weight": 0.10})

    # Sleep / circadian
    chrono = exercise.get("chronotype", {}).get("chronotype", "Neutral")
    sleep_score = 75 if "Morning" in chrono else 60
    components.append({"component": "Sleep / circadian alignment", "score": sleep_score, "weight": 0.10})

    total = sum(c["score"] * c["weight"] for c in components)

    # Identify the top three improvable levers (lowest scoring × highest weight)
    leverable = sorted(components, key=lambda c: c["score"] * (1 / max(0.05, c["weight"])))
    biggest_levers = [
        {
            "lever": leverable[0]["component"],
            "current_score": leverable[0]["score"],
            "improvement_action": _improvement_lever(leverable[0]["component"], nutrition, exercise),
        },
        {
            "lever": leverable[1]["component"],
            "current_score": leverable[1]["score"],
            "improvement_action": _improvement_lever(leverable[1]["component"], nutrition, exercise),
        },
        {
            "lever": leverable[2]["component"],
            "current_score": leverable[2]["score"],
            "improvement_action": _improvement_lever(leverable[2]["component"], nutrition, exercise),
        },
    ]

    return {
        "composite_score": round(total, 1),
        "tier": _longevity_tier(total),
        "components": components,
        "biggest_levers": biggest_levers,
        "interpretation": (
            f"Composite longevity index: {round(total, 1)}/100. This score integrates "
            f"genetic cardiometabolic risk, inflammatory tone, APOE status, athletic "
            f"capacity, musculoskeletal resilience, recovery, and circadian alignment. "
            f"It is a snapshot of HEADWINDS and TAILWINDS — actual lifespan outcomes are "
            f"~70% behavioural, so even a low genetic score can be substantially offset "
            f"by lifestyle, and a high score can be squandered."
        ),
    }


def _longevity_tier(score: float) -> str:
    if score >= 80:
        return "Tailwind — strong genetic + capacity profile"
    if score >= 65:
        return "Above average — modest headwinds, easily counterable"
    if score >= 50:
        return "Average — usual mix of strengths and risks"
    if score >= 35:
        return "Headwind — multiple risk axes; behavioural levers matter more"
    return "Strong headwind — lifestyle interventions are mandatory, not optional"


def _improvement_lever(component: str, nutrition: dict, exercise: dict) -> str:
    if "Cardiometabolic" in component:
        return ("Adopt Mediterranean pattern strictly, eliminate sugar-sweetened beverages, "
                "150+ min Z2/wk, achieve 25 g oat β-glucan + 2 g plant sterols/day. "
                "Re-check lipid panel in 12 weeks — expect 15-30% LDL drop.")
    if "Inflammatory" in component:
        return ("4×/wk oily fish, 7×/wk berries, eliminate ultra-processed foods, "
                "turmeric + black pepper daily, 8 h sleep target. Re-check hs-CRP in 8 weeks.")
    if "APOE" in component:
        return ("ε4 carrier — strict saturated-fat cap (<7%), MIND diet adherence, "
                "150+ min weekly aerobic, sleep 8 h to support glymphatic Aβ clearance, "
                "weekly oily fish for DHA.")
    if "Athletic" in component:
        return ("Single biggest lever — couch-to-active alone shifts all-cause mortality "
                "30-45%. Even 7000 daily steps + 2 strength sessions/wk yield massive gains.")
    if "Musculoskeletal" in component:
        return ("2-3 strength sessions/wk (especially heavy hinge and squat), "
                "daily 10-min mobility, eccentric tendon work for high-risk regions, "
                "vitamin D >40 ng/mL.")
    if "Recovery" in component:
        return ("Sleep ≥8 h, every-6-week deload, omega-3 2 g/day, anti-inflammatory diet, "
                "HRV-guided training intensity.")
    if "Sleep" in component or "circadian" in component:
        return ("Strict sleep/wake times (±30 min), morning sunlight 10-15 min, "
                "caffeine cutoff 8 h pre-bed, dark cool bedroom (16-19 °C), "
                "no screens 30 min pre-sleep.")
    return "Generic lifestyle: sleep, food quality, movement, social connection."


# ── Year-Long Integrated Plan ──────────────────────────────────────────────

def year_long_plan(nutrition: dict, exercise: dict) -> list[dict]:
    bias = exercise.get("power_endurance", {}).get("bias", "Balanced")
    sat = nutrition.get("saturated_fat", {}).get("apoe_genotype", "")
    inflam_tier = nutrition.get("inflammation", {}).get("tier", "")
    if bias.startswith("Power"):
        annual = [
            "Hypertrophy", "Strength", "Power/peak", "Active recovery + sport"
        ]
    elif bias.startswith("Endurance"):
        annual = [
            "Aerobic base", "Threshold build", "Race-prep peak", "Off-season + strength"
        ]
    else:
        annual = [
            "GPP / hypertrophy", "Strength + endurance build", "Specific peak", "Active recovery"
        ]
    blocks = []
    for i, focus in enumerate(annual):
        block_n = i + 1
        weeks = f"Weeks {(i*13)+1}-{(i+1)*13}"
        nut_focus = _quarterly_nutrition_focus(block_n, sat, inflam_tier)
        ex_focus = _quarterly_exercise_focus(block_n, bias)
        labs = _quarterly_labs(block_n)
        blocks.append({
            "mesocycle": f"Q{block_n} — {focus}",
            "weeks": weeks,
            "exercise_focus": ex_focus,
            "nutrition_focus": nut_focus,
            "labs_to_recheck": labs,
            "goal_milestone": _quarterly_milestone(block_n, bias),
        })
    return blocks


def _quarterly_nutrition_focus(q: int, apoe: str, inflam: str) -> str:
    base = {
        1: "Establish 30-day meal-plan adherence; lock breakfast/lunch routine; macros at target.",
        2: "Reintroduce dietary diversity; test 16:8 TRF if appropriate; track sleep + AM glucose.",
        3: "Race-prep or strength-peak fueling — carb periodisation around training days.",
        4: "Diet-break / maintenance: maintenance calories, broader food variety, fewer rules.",
    }[q]
    if apoe.startswith("ε4"):
        base += " MIND-diet items every week. Sat-fat <7% strict."
    if "Elevated" in inflam:
        base += " Anti-inflammatory emphasis: oily fish 4×/wk, polyphenols 1500 mg/day."
    return base


def _quarterly_exercise_focus(q: int, bias: str) -> str:
    if bias.startswith("Power"):
        return [
            "Volume accumulation 8-12 reps; movement quality; build base aerobic 2×/wk.",
            "Heavy strength 3-6 reps; reduce metabolic conditioning to maintenance.",
            "Power/peaking 1-3 reps + plyo; sport-specific; competition simulation.",
            "Active rest — alternate sport, mobility focus, restore tendons.",
        ][q-1]
    if bias.startswith("Endurance"):
        return [
            "Aerobic base 70-80% Z2; introduce 1×/wk strides + 1 strength day.",
            "Threshold work 2×/wk; weekly long progression; concurrent strength maintenance.",
            "VO2max + race-specific intensity; taper into goal event.",
            "Off-season cross-training; rebuild strength base; reduce running volume 50%.",
        ][q-1]
    return [
        "Balanced GPP: 3 lifts + 3 aerobic. Build movement quality and base fitness.",
        "Strength block dominant (4 lift / 2 aerobic) with 1 long Z2.",
        "Mixed-modal peak with specific competition prep.",
        "Active recovery: alternate sport, yoga, mobility focus.",
    ][q-1]


def _quarterly_labs(q: int) -> list[str]:
    schedule = {
        1: ["Full lipid panel (Total/LDL/HDL/Triglycerides + ApoB)",
            "hs-CRP", "Fasting glucose + HbA1c + insulin", "Vitamin D 25(OH)",
            "Ferritin + transferrin saturation", "CBC + comprehensive metabolic panel",
            "Magnesium, B12, folate, homocysteine"],
        2: ["Skip — establish behavioural changes; no remeasure until 12 weeks in."],
        3: ["hs-CRP (track inflammation response)", "HbA1c", "Vitamin D",
            "Lipid panel if making major dietary shifts"],
        4: ["Full panel repeat — track 12-month change vs Q1 baseline.",
            "Optional: ApoB, Lp(a) — best lipid markers if not done initially.",
            "DEXA scan (body comp + bone density)", "VO2max test if endurance-focused."],
    }
    return schedule[q]


def _quarterly_milestone(q: int, bias: str) -> str:
    if bias.startswith("Power"):
        return ["Movement quality: clean reps at moderate loads.",
                "Strength PRs: +5-10% on squat/dead/bench from Q1.",
                "Peak: hit competition-day numbers.",
                "Maintain 90% of peak with reduced volume."][q-1]
    if bias.startswith("Endurance"):
        return ["Aerobic base: Z2 pace ~10s/km faster at same HR.",
                "Threshold pace improved 3-5%.",
                "Goal race performance.",
                "Maintain 80% of fitness with reduced volume."][q-1]
    return ["Body comp shift (5% leaner or +2 kg muscle).",
            "Strength + endurance benchmarks improved.",
            "Goal event or test outcome.",
            "Maintain across active rest."][q-1]


# ── Body Composition Trajectory ────────────────────────────────────────────

def body_composition_trajectory(
    body_weight_kg: float, body_fat_pct: float | None,
    goal: str, nutrition: dict, exercise: dict
) -> dict:
    """
    goal: 'fat_loss' | 'lean_gain' | 'maintenance' | 'recomposition'
    Returns 12-week realistic projection.
    """
    caloric = nutrition.get("caloric", {})
    strength_tier = exercise.get("strength_trainability", {}).get("tier", "")
    fto_risk = nutrition.get("satiety", {}).get("appetite_phenotype", "Standard") != "Standard"

    if goal == "fat_loss":
        rate_kg_per_wk = 0.5 if fto_risk else 0.7
        wks = 12
        expected_loss_kg = round(rate_kg_per_wk * wks, 1)
        guidance = (
            f"Realistic 12-week fat loss: {expected_loss_kg} kg (~{expected_loss_kg/wks:.1f} kg/wk). "
            f"This is sustainable; faster rates accelerate muscle loss and metabolic adaptation. "
            f"With strength training, expect to PRESERVE or even gain a small amount of lean mass."
        )
        if fto_risk:
            guidance += " Elevated-appetite genotype: slower deficit is essential for adherence."
        return {
            "goal": "Fat loss",
            "rate": f"{rate_kg_per_wk} kg/wk",
            "expected_12wk_change": f"-{expected_loss_kg} kg total",
            "caloric_advice": (f"Subtract {caloric.get('loss_deficit_kcal', 500)} kcal/day from TDEE. "
                               f"Protein floor {caloric.get('protein_g_per_kg', 1.8)} g/kg body weight."),
            "training": "Maintain strength training 2-3×/wk to preserve lean mass; add aerobic for adherence + cardio health.",
            "expected_body_fat_drop": "≈3-5 percentage points achievable over 12 weeks at this rate.",
            "guidance": guidance,
        }
    if goal == "lean_gain":
        rate = 0.4 if "High" in strength_tier else 0.25
        gain_total = round(rate * 12, 1)
        return {
            "goal": "Lean gain",
            "rate": f"{rate} kg/wk",
            "expected_12wk_change": f"+{gain_total} kg total (mostly lean for high-hypertrophy responders)",
            "caloric_advice": f"+{caloric.get('gain_surplus_kcal', 250)} kcal/day surplus.",
            "training": "Strength training 4×/wk; progressive overload; ≥1.8 g/kg protein; ≥7 h sleep.",
            "expected_body_fat_change": "+1-2 percentage points if calorie surplus is moderate.",
            "guidance": "Expect 60-80% of gain to be lean mass with proper training; the rest is fat + water.",
        }
    if goal == "recomposition":
        return {
            "goal": "Recomposition (gain muscle while losing fat)",
            "feasibility": "Best for novice trainees, those returning from a layoff, or anyone with high body fat. "
                           "Advanced lean lifters can't easily recomp — choose a phase instead.",
            "expected_12wk_change": "-1 to -3 kg fat, +0.5-2 kg lean mass.",
            "caloric_advice": "Maintenance or very slight deficit (-100 to -200 kcal/day) + high protein (2.0 g/kg).",
            "training": "Strength training 3-4×/wk with progressive overload is non-negotiable.",
        }
    return {
        "goal": "Maintenance",
        "expected_12wk_change": "Body composition stable. Performance gains independent of body comp.",
        "caloric_advice": "Match TDEE; don't fixate on numbers — energy + performance + recovery are the markers.",
    }


# ── Executive One-Page Summary ─────────────────────────────────────────────

def executive_summary(nutrition: dict, exercise: dict, longevity: dict) -> dict:
    return {
        "headline_score": longevity["composite_score"],
        "headline_tier": longevity["tier"],
        "key_facts": [
            f"Macro target: {nutrition['macros']['pct_carbs']}C / "
            f"{nutrition['macros']['pct_fat']}F / {nutrition['macros']['pct_protein']}P",
            f"Caffeine: {nutrition['caffeine'].get('metabolism','—')} — "
            f"cap {nutrition['caffeine'].get('limit_mg','—')} mg/day, cutoff {nutrition['caffeine'].get('cutoff_time','—')}",
            f"Alcohol risk: {nutrition['alcohol'].get('risk','—')}",
            f"APOE: {nutrition.get('saturated_fat',{}).get('apoe_genotype','—')} — "
            f"sat-fat cap {nutrition.get('saturated_fat',{}).get('saturated_fat_cap_g','—')} g/day",
            f"Glycaemic ceiling: {nutrition.get('glycemic_threshold',{}).get('max_carbs_per_meal_g','—')} g/meal",
            f"Power/Endurance bias: {exercise['power_endurance']['bias']} "
            f"({exercise['power_endurance']['ratio_pct_power']}/{exercise['power_endurance']['ratio_pct_endurance']})",
            f"VO2max trainability: {exercise.get('vo2max',{}).get('tier','—')}",
            f"Recovery: {exercise['recovery']['speed']}",
            f"Chronotype: {exercise['chronotype']['chronotype']} — train {exercise['chronotype']['optimal_window']}",
            f"Top sport match: {exercise.get('composite_profile',{}).get('ranked_sports',[{}])[0].get('sport','—')}",
        ],
        "top_three_actions_this_week": [
            longevity["biggest_levers"][0]["improvement_action"],
            longevity["biggest_levers"][1]["improvement_action"],
            longevity["biggest_levers"][2]["improvement_action"],
        ],
        "non_negotiable_daily": [
            f"≥{nutrition.get('caloric',{}).get('protein_g_per_kg',1.6)} g/kg protein",
            f"≥{nutrition.get('fiber',{}).get('target_g',38)} g fibre",
            "Sleep 7.5-9 h",
            "8-10k steps minimum",
            "10 min mobility",
            f"Water ≈{nutrition.get('hydration',{}).get('target_ml_per_kg',33)} mL/kg",
        ],
        "weekly_floors": [
            "≥150 min moderate aerobic OR 75 min vigorous",
            "≥2 strength sessions (compound lifts)",
            "≥4 servings oily fish OR 1 tbsp ground flax/chia daily",
            "≥1 cup leafy greens daily, ≥3 servings cruciferous/wk",
            "≥1 deload week every 4-6 weeks",
        ],
    }


# ── Public synthesis API ───────────────────────────────────────────────────

def integrated_longevity_plan(
    nutrition_result: dict, exercise_result: dict,
    body_weight_kg: float | None = None,
    body_fat_pct: float | None = None,
    body_comp_goal: str = "maintenance",
) -> dict:
    longevity = longevity_composite(nutrition_result, exercise_result)
    annual_plan = year_long_plan(nutrition_result, exercise_result)
    summary = executive_summary(nutrition_result, exercise_result, longevity)
    trajectory = None
    if body_weight_kg is not None:
        trajectory = body_composition_trajectory(
            body_weight_kg, body_fat_pct, body_comp_goal,
            nutrition_result, exercise_result,
        )
    return {
        "executive_summary": summary,
        "longevity_composite": longevity,
        "year_long_plan": annual_plan,
        "body_composition_trajectory": trajectory,
    }
