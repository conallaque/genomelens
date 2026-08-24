"""
Advanced Nutrition Analytics
============================

Polygenic-score synthesis, cardiometabolic dashboard, inflammation index,
histamine intolerance, phase I/II detoxification, personal glycemic threshold,
training-day vs rest-day macro periodisation, a quantitative weekly food matrix,
a shopping list, and three generated recipes that match the genotype constraints.

These are *evidence-aligned* simplifications of published GWAS effect sizes —
the scores are clinically suggestive, not diagnostic. They expose the variance
that is reasonably attributable to common SNPs (the "tip of the iceberg"
typically explaining 5-25% of heritable risk).
"""

from __future__ import annotations

import pandas as pd


def _gt(snps_df: pd.DataFrame | None, rsid: str) -> str | None:
    if snps_df is None or rsid not in snps_df.index:
        return None
    raw = snps_df.loc[rsid].get("genotype")
    if raw is None:
        return None
    s = str(raw).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN"):
        return None
    return s


def _dose(gt: str | None, risk: str) -> int:
    """Return 0, 1, or 2 — number of risk alleles."""
    if not gt:
        return 0
    return gt.count(risk)


# ── Polygenic Score Engine ──────────────────────────────────────────────────
#
# Each PGS sums weighted log-odds (or β) for the risk allele across loci. We
# present a *relative score* (vs typed-loci maximum) and a population-percentile
# estimate using a Gaussian approximation. Weights are from large-cohort meta-
# analyses (DIAGRAM, GIANT, GLGC, etc.) — magnitudes simplified for tractability.

_PGS_DEFS = {
    "T2D": [
        ("rs7903146",  "T", 0.34, "TCF7L2"),
        ("rs10830963", "G", 0.10, "MTNR1B"),
        ("rs5219",     "T", 0.08, "KCNJ11"),
        ("rs13266634", "C", 0.10, "SLC30A8"),
        ("rs1801282",  "C", 0.14, "PPARG"),
        ("rs7754840",  "C", 0.08, "CDKAL1"),
        ("rs9939609",  "A", 0.10, "FTO"),
    ],
    "Obesity_BMI": [
        ("rs9939609",  "A", 0.39, "FTO"),
        ("rs17782313", "C", 0.23, "MC4R"),
        ("rs1137101",  "G", 0.10, "LEPR"),
        ("rs6548238",  "C", 0.13, "TMEM18"),
        ("rs7359397",  "T", 0.07, "SH2B1"),
        ("rs10938397", "G", 0.19, "GNPDA2"),
    ],
    "LDL_cholesterol": [
        ("rs429358",   "C", 0.34, "APOE ε4"),
        ("rs6511720",  "G", 0.18, "LDLR"),
        ("rs646776",   "T", 0.10, "SORT1"),
        ("rs12740374", "G", 0.10, "CELSR2"),
        ("rs11206510", "T", 0.07, "PCSK9"),
        ("rs5082",     "C", 0.08, "APOA2"),
    ],
    "Triglycerides": [
        ("rs964184",   "G", 0.28, "ZPR1/APOA5"),
        ("rs662799",   "G", 0.16, "APOA5"),
        ("rs1260326",  "T", 0.10, "GCKR"),
        ("rs1748195",  "C", 0.07, "ANGPTL3"),
        ("rs328",      "G", 0.18, "LPL"),
        ("rs174547",   "T", 0.07, "FADS1"),
    ],
    "HDL_cholesterol": [
        ("rs3764261",  "C", 0.30, "CETP"),  # C = HDL-lowering
        ("rs1800775",  "C", 0.10, "CETP -629"),
        ("rs328",      "C", 0.18, "LPL"),   # C = lower HDL
        ("rs1800588",  "C", 0.10, "LIPC"),
    ],
    "Caffeine_consumption": [
        ("rs2472297",  "T", 0.18, "CYP1A2/15q24"),
        ("rs4410790",  "C", 0.13, "AHR"),
    ],
    "Alcohol_dependence_protection": [
        ("rs1229984",  "A", 0.85, "ADH1B*2"),   # A = strongly protective
        ("rs671",      "A", 1.20, "ALDH2*2"),
    ],
    "Coeliac_risk": [
        ("rs2187668",  "T", 0.85, "HLA-DQ2.5"),
        ("rs7454108",  "C", 0.55, "HLA-DQ8"),
    ],
    "Vitamin_D_deficiency": [
        ("rs10741657", "G", 0.13, "CYP2R1"),
        ("rs2282679",  "C", 0.16, "GC/VDBP"),
        ("rs12785878", "G", 0.10, "DHCR7"),
    ],
    "BP_systolic": [
        ("rs699",      "G", 0.10, "AGT"),
        ("rs4341",     "T", 0.07, "ACE"),
        ("rs5068",     "A", 0.05, "NPPA"),
        ("rs17367504", "G", 0.10, "MTHFR-NPPA"),
    ],
}


def _score_pgs(snps_df, defn: list[tuple]) -> dict:
    raw = 0.0
    max_possible = 0.0
    loci_typed: list[str] = []
    loci_missing: list[str] = []
    for rsid, risk_allele, weight, label in defn:
        gt = _gt(snps_df, rsid)
        max_possible += 2 * weight
        if gt is None:
            loci_missing.append(label)
            continue
        loci_typed.append(f"{label} ({rsid})={gt}")
        raw += weight * _dose(gt, risk_allele)
    coverage = (len(loci_typed) / len(defn)) if defn else 0
    if max_possible == 0 or coverage == 0:
        return {
            "score": 0.0, "max_possible": max_possible, "coverage": 0,
            "z": None, "percentile": None,
            "typed": loci_typed, "missing": loci_missing,
        }
    # Normalise: assume population mean ≈ max/2, SD ≈ max/6 (Gaussian approx)
    pop_mean = max_possible / 2
    pop_sd = max_possible / 6 if max_possible > 0 else 1
    z = (raw - pop_mean) / pop_sd
    # Approximate Gaussian CDF (Abramowitz formula)
    pct = _phi(z) * 100
    return {
        "score": round(raw, 3),
        "max_possible": round(max_possible, 3),
        "coverage": round(coverage * 100, 1),
        "z": round(z, 2),
        "percentile": round(pct, 1),
        "typed": loci_typed,
        "missing": loci_missing,
    }


def _phi(z: float) -> float:
    """Standard normal CDF (Abramowitz & Stegun 26.2.17)."""
    import math
    t = 1 / (1 + 0.2316419 * abs(z))
    pdf = math.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    poly = (((1.330274429 * t - 1.821255978) * t + 1.781477937) * t
            - 0.356563782) * t + 0.319381530
    cdf = 1 - pdf * poly * t
    return cdf if z >= 0 else 1 - cdf


def _tier_from_percentile(p: float | None) -> str:
    if p is None:
        return "Insufficient data"
    if p >= 90:
        return "Very high"
    if p >= 75:
        return "High"
    if p >= 60:
        return "Above average"
    if p >= 40:
        return "Average"
    if p >= 25:
        return "Below average"
    if p >= 10:
        return "Low"
    return "Very low"


def analyze_polygenic_scores(snps_df) -> dict:
    out: dict[str, dict] = {}
    for trait, defn in _PGS_DEFS.items():
        s = _score_pgs(snps_df, defn)
        s["tier"] = _tier_from_percentile(s.get("percentile"))
        out[trait] = s
    return out


# ── Cardiometabolic Dashboard ───────────────────────────────────────────────
#
# Synthesises five lipid/glycaemic axes into a quantitative profile with
# specific dietary leverage points.

def cardiometabolic_dashboard(pgs: dict) -> dict:
    axes = []
    for key, label, low_advice, high_advice in [
        ("LDL_cholesterol", "LDL response to saturated fat",
         "LDL stays low even on higher-sat-fat diets — Mediterranean fine, dairy OK.",
         "LDL spikes sharply with saturated fat — cap sat fat ≤7%, prefer MUFA/PUFA, "
         "add 2 g plant sterols/day (fortified margarine), 25 g oats β-glucan."),
        ("Triglycerides", "Postprandial triglyceride response",
         "TG handling robust — moderate refined carbs tolerated.",
         "TG spike postprandially — minimise added sugar (<25 g/day), alcohol "
         "(major TG driver), refined grains. Push omega-3 to 2 g EPA+DHA/day."),
        ("HDL_cholesterol", "HDL maintenance",
         "HDL maintenance favourable — standard exercise/diet sustains it.",
         "Genetic tendency toward low HDL — high-intensity exercise (HIIT, lifting), "
         "MUFA emphasis (olive, almonds, avocado), moderate alcohol (if no ALDH2 issue)."),
        ("T2D", "Glycaemic vulnerability",
         "Robust glucose tolerance — standard carb intake fine.",
         "Reduced β-cell capacity / insulin sensitivity — protein-first meal sequencing "
         "(eat protein+veg before carbs cuts postprandial glucose 25%), 10-min "
         "post-meal walks, vinegar 1 tbsp before high-carb meals, chromium-rich foods."),
        ("BP_systolic", "Blood-pressure response to sodium",
         "BP stays low on standard sodium intake.",
         "Sodium-sensitive BP — DASH pattern (≤2.3 g Na, 4.7 g K from leafy greens, "
         "potato, beans, banana), nitrate-rich beets/leafy greens (lower SBP 4-5 mmHg), "
         "magnesium 400 mg from nuts/seeds."),
    ]:
        s = pgs.get(key, {})
        pct = s.get("percentile")
        if pct is None:
            advice = "Insufficient typed variants — assume average."
        elif pct >= 70:
            advice = high_advice
        elif pct <= 30:
            advice = low_advice
        else:
            advice = "Average — standard guidance applies."
        axes.append({
            "axis": label,
            "percentile": pct,
            "tier": s.get("tier", "Insufficient data"),
            "leverage": advice,
        })
    return {"axes": axes}


# ── Inflammation Index ──────────────────────────────────────────────────────

def inflammation_index(snps_df) -> dict:
    il6 = _gt(snps_df, "rs1800795")
    crp = _gt(snps_df, "rs2794520") or _gt(snps_df, "rs1205")
    tnf = _gt(snps_df, "rs1800629")
    score = 0
    factors: list[str] = []
    if il6:
        factors.append(f"rs1800795 (IL6 -174) {il6}")
        if "G" in il6:
            score += il6.count("G")
    if crp:
        factors.append(f"CRP {crp}")
        if "C" in crp:
            score += crp.count("C") * 0.5
    if tnf:
        factors.append(f"rs1800629 (TNF-α -308) {tnf}")
        if "A" in tnf:
            score += tnf.count("A")
    if score >= 3:
        tier = "Elevated inflammatory tone"
        diet = (
            "Adopt a strict anti-inflammatory diet: oily fish 4×/week, 1 tbsp ground "
            "flaxseed daily, 1-2 cups berries, 2-3 cups cruciferous veg/week, turmeric "
            "(1 tsp + black pepper + fat), green tea 3 cups/day, walnuts 1 oz/day, "
            "extra-virgin olive oil 2-3 tbsp/day. AVOID: sugar-sweetened drinks, refined "
            "carbs, industrial seed oils (corn, soybean as primary), processed meats, "
            "fried foods. Target Mediterranean-DASH-MIND hybrid pattern. Specific aim: "
            "reduce hs-CRP by 30-50% over 12 weeks."
        )
    elif score >= 1:
        tier = "Mild inflammatory tendency"
        diet = "Mediterranean pattern with 2× weekly fish, daily berries + nuts, 1 tsp turmeric, minimise processed foods."
    else:
        tier = "Low baseline inflammation"
        diet = "Standard balanced diet; no special anti-inflammatory emphasis needed."
    return {"tier": tier, "score": score, "factors": factors or ["Inflammation SNPs not typed"],
            "guidance": diet}


# ── Histamine intolerance ───────────────────────────────────────────────────

def histamine_intolerance(snps_df) -> dict:
    dao1 = _gt(snps_df, "rs10156191")
    dao2 = _gt(snps_df, "rs1049742")
    hnmt = _gt(snps_df, "rs11558538")
    score = 0
    factors: list[str] = []
    if dao1:
        factors.append(f"rs10156191 (DAO/AOC1) {dao1}")
        if "T" in dao1:
            score += dao1.count("T")
    if dao2:
        factors.append(f"rs1049742 (DAO) {dao2}")
        if "T" in dao2:
            score += dao2.count("T")
    if hnmt:
        factors.append(f"rs11558538 (HNMT) {hnmt}")
        if "T" in hnmt:
            score += hnmt.count("T")
    risk = score >= 2
    return {
        "elevated_risk": risk, "score": score,
        "factors": factors or ["Histamine SNPs not typed"],
        "guidance": (
            "Reduced histamine clearance — if you experience flushing, headaches, "
            "hives, congestion after fermented/aged foods, try a 4-week low-histamine "
            "trial. AVOID/LIMIT: aged cheeses, cured meats, fermented foods (sauerkraut, "
            "kombucha, soy sauce), wine, vinegar, spinach, tomatoes, eggplant, avocado "
            "(unripe), shellfish, chocolate, citrus. SAFER: fresh-cooked meats, white "
            "rice, most cooked vegetables, apples, pears, fresh dairy. Vitamin C and "
            "quercetin (in onions, apples) support DAO activity."
            if risk else
            "Normal histamine clearance — fermented and aged foods well tolerated; "
            "they contribute meaningfully to gut-microbiome diversity."
        ),
    }


# ── Detoxification (Phase I + Phase II) ─────────────────────────────────────

def detoxification_profile(snps_df) -> dict:
    # Phase I — CYP enzymes
    cyp1a2 = _gt(snps_df, "rs762551")           # caffeine, polycyclic aromatic hydrocarbons
    cyp1a1 = _gt(snps_df, "rs1048943")          # estrogen 2-hydroxylation
    cyp2d6_proxy = _gt(snps_df, "rs1065852")
    cyp3a4 = _gt(snps_df, "rs2740574")
    # Phase II — conjugation
    gstm1_proxy = _gt(snps_df, "rs366631")
    gstp1 = _gt(snps_df, "rs1695")              # I105V; G = lower activity
    nat2 = _gt(snps_df, "rs1799929") or _gt(snps_df, "rs1799930")
    sult1a1 = _gt(snps_df, "rs9282861")
    comt = _gt(snps_df, "rs4680")
    phase1: list[str] = []
    phase2: list[str] = []
    for label, gt in [("CYP1A2", cyp1a2), ("CYP1A1", cyp1a1), ("CYP2D6", cyp2d6_proxy),
                      ("CYP3A4", cyp3a4)]:
        if gt:
            phase1.append(f"{label} {gt}")
    for label, gt in [("GSTM1", gstm1_proxy), ("GSTP1", gstp1), ("NAT2", nat2),
                      ("SULT1A1", sult1a1), ("COMT", comt)]:
        if gt:
            phase2.append(f"{label} {gt}")
    # Slow phase II flag — GSTP1 GG or COMT AA
    slow_phase2 = (gstp1 and "G" in gstp1 and gstp1.count("G") == 2) or \
                  (comt and "A" in comt and comt.count("A") == 2)
    return {
        "phase1_typed": phase1 or ["Phase I SNPs not typed"],
        "phase2_typed": phase2 or ["Phase II SNPs not typed"],
        "slow_phase2_clearance": bool(slow_phase2),
        "guidance": (
            "Slow Phase II conjugation flagged. Phase I activates carcinogens — Phase II "
            "neutralises them; imbalance increases reactive-intermediate exposure. "
            "Daily support: 1 cup cruciferous veg (broccoli sprouts ideal — sulforaphane "
            "induces NQO1/GST), 2-3 alliums (garlic, onion — sulfur for sulfation), "
            "1 cup berries (anthocyanins), 1 tsp turmeric, glycine-rich foods (bone broth, "
            "gelatin) for glycine conjugation, NAC- or whey-derived cysteine for "
            "glutathione synthesis. Minimise alcohol, char-grilled meats, smoking, and "
            "high-dose acetaminophen."
            if slow_phase2 else
            "Phase I/II balance acceptable. Support with cruciferous veg 3×/week + "
            "alliums + adequate protein for amino-acid conjugation pool."
        ),
        "cruciferous_target_servings_per_week": 5 if slow_phase2 else 3,
    }


# ── Personal Glycemic Threshold ─────────────────────────────────────────────
#
# Estimates the maximum carb load (g) per single meal at which postprandial
# glucose is likely to stay <140 mg/dL. Drawn from TCF7L2/MTNR1B/GCKR effect
# sizes (each risk allele ~5-10 g/meal tolerance reduction).

def glycemic_threshold(snps_df) -> dict:
    base = 75  # g/meal — population average for a healthy adult
    factors: list[str] = []
    adj = 0
    tcf = _gt(snps_df, "rs7903146")
    mtnr = _gt(snps_df, "rs10830963")
    gckr = _gt(snps_df, "rs1260326")
    pparg = _gt(snps_df, "rs1801282")
    if tcf:
        factors.append(f"TCF7L2 {tcf}")
        adj -= 10 * tcf.count("T")
    if mtnr:
        factors.append(f"MTNR1B {mtnr}")
        adj -= 6 * mtnr.count("G")
    if gckr:
        factors.append(f"GCKR {gckr}")
        adj -= 4 * gckr.count("T")
    if pparg:
        factors.append(f"PPARG Pro12Ala {pparg}")
        # Ala12 (G allele) protective
        if "G" in pparg:
            adj += 5 * pparg.count("G")
    threshold = max(30, base + adj)
    dinner_threshold = max(25, threshold - 15)  # evening glucose tolerance lower
    return {
        "max_carbs_per_meal_g": threshold,
        "max_carbs_dinner_g": dinner_threshold,
        "factors": factors or ["Glycaemic SNPs not typed"],
        "guidance": (
            f"Personal carb ceiling: ~{threshold} g/meal at breakfast/lunch, "
            f"~{dinner_threshold} g at dinner (circadian glucose tolerance drops PM). "
            f"Reference portions: 1 cup cooked rice ≈ 45 g, 1 medium banana ≈ 27 g, "
            f"1 slice bread ≈ 15 g, 1 cup oats ≈ 27 g, 1 medium potato ≈ 37 g, "
            f"1 cup berries ≈ 14 g. Pairing with 25-30 g protein + 10 g fibre + "
            f"healthy fat blunts glucose excursion 30-50%."
        ),
    }


# ── Training-day vs Rest-day Macro Periodisation ────────────────────────────

def macro_periodisation(macros: dict, glycemic: dict, satiety: dict) -> dict:
    base_carb = macros["pct_carbs"]
    train_day_carb = min(60, base_carb + 15)
    rest_day_carb = max(15, base_carb - 10)
    train_fat = max(20, 100 - train_day_carb - 25)
    rest_fat = max(30, 100 - rest_day_carb - 30)
    return {
        "training_day": {
            "pct_carbs": train_day_carb,
            "pct_fat": train_fat,
            "pct_protein": 100 - train_day_carb - train_fat,
            "carb_timing": (
                f"60% of carbs in pre+post-workout window ({glycemic['max_carbs_per_meal_g']} g "
                "pre-session 60 min before, refill within 2 h post). Remaining 40% spread "
                "across other meals."
            ),
        },
        "rest_day": {
            "pct_carbs": rest_day_carb,
            "pct_fat": rest_fat,
            "pct_protein": 100 - rest_day_carb - rest_fat,
            "carb_timing": "Lower-carb day — emphasise protein, vegetables, healthy fats. Useful for insulin sensitivity restoration.",
        },
        "guidance": (
            f"Cycling carbs by training load improves body composition and insulin "
            f"sensitivity vs flat macros. {('Elevated-appetite genotype: limit carb cycling depth — large '+'-'+'-day swings can trigger evening cravings.' if satiety['appetite_phenotype'] != 'Standard' else 'Standard appetite — carb cycling well-tolerated.')}"
        ),
    }


# ── Quantitative Weekly Food Matrix ─────────────────────────────────────────

def weekly_food_matrix(result_partial: dict) -> dict:
    """
    Convert genotype findings into a concrete weekly serving target per food
    group. Result_partial is the in-progress nutrition dict (must contain
    omega3, iron, choline, antioxidants, satiety, saturated_fat, fiber).
    """
    omega3 = result_partial.get("omega3", {})
    iron = result_partial.get("iron", {})
    choline = result_partial.get("choline", {})
    antiox = result_partial.get("antioxidants", {})
    satiety = result_partial.get("satiety", {})
    sat_fat = result_partial.get("saturated_fat", {})
    lact = result_partial.get("lactose", {})
    gluten = result_partial.get("gluten", {})
    detox = result_partial.get("detoxification", {})
    inflam = result_partial.get("inflammation", {})
    cof = result_partial.get("caffeine", {})
    alc = result_partial.get("alcohol", {})

    # Defaults (servings/week)
    matrix = {
        "Oily fish (salmon/sardines/mackerel)": 2,
        "Lean protein (chicken/turkey/lean beef)": 3,
        "Whole eggs": 4,
        "Legumes (beans/lentils/chickpeas)": 4,
        "Tofu/tempeh": 1,
        "Greek yoghurt / kefir": 3,
        "Aged cheese (matchbox-sized)": 2,
        "Leafy greens (1 cup cooked or 2 raw)": 7,
        "Cruciferous veg (broccoli/cauliflower/Brussels)": 3,
        "Allium veg (garlic/onion/leek)": 5,
        "Colorful veg (bell pepper/carrot/beet)": 5,
        "Berries (1/2 cup)": 5,
        "Whole fruit (apple/pear/citrus)": 7,
        "Nuts (1 oz)": 5,
        "Seeds (chia/flax/pumpkin, 1 tbsp)": 7,
        "Whole grains (oats/quinoa/brown rice)": 5,
        "Extra-virgin olive oil (tbsp)": 14,
        "Red meat (lean cut, palm-sized)": 1,
        "Processed/cured meat": 0,
        "Sugary drinks": 0,
        "Coffee (8 oz)": 7,
        "Green tea (cups)": 3,
        "Dark chocolate (≥75%, 1 oz)": 3,
        "Alcohol (standard drink)": 3,
    }

    if omega3.get("ala_conversion") == "Poor":
        matrix["Oily fish (salmon/sardines/mackerel)"] = 4
        matrix["Seeds (chia/flax/pumpkin, 1 tbsp)"] = 10
    if iron.get("overload_risk", "Low").startswith("High") or "Moderate" in iron.get("overload_risk", ""):
        matrix["Lean protein (chicken/turkey/lean beef)"] = 2
        matrix["Red meat (lean cut, palm-sized)"] = 0
    if choline.get("increased_need"):
        matrix["Whole eggs"] = 7
    if antiox.get("reduced_capacity") or (detox and detox.get("slow_phase2_clearance")):
        matrix["Cruciferous veg (broccoli/cauliflower/Brussels)"] = 5
        matrix["Berries (1/2 cup)"] = 7
        matrix["Green tea (cups)"] = 5
    if satiety.get("appetite_phenotype") != "Standard":
        matrix["Legumes (beans/lentils/chickpeas)"] = 6
        matrix["Leafy greens (1 cup cooked or 2 raw)"] = 10
    if sat_fat.get("apoe_genotype", "").startswith("ε4"):
        matrix["Red meat (lean cut, palm-sized)"] = 1
        matrix["Aged cheese (matchbox-sized)"] = 1
        matrix["Extra-virgin olive oil (tbsp)"] = 21
        matrix["Nuts (1 oz)"] = 7
    if lact.get("tolerance", "").startswith("Intolerant"):
        matrix["Greek yoghurt / kefir"] = 5  # cultured = lower lactose
        matrix["Aged cheese (matchbox-sized)"] = 3
    if gluten.get("celiac_risk_haplotype"):
        matrix["Whole grains (oats/quinoa/brown rice)"] = 3  # gluten-free emphasis
    if inflam.get("score", 0) >= 2:
        matrix["Oily fish (salmon/sardines/mackerel)"] = max(matrix["Oily fish (salmon/sardines/mackerel)"], 4)
        matrix["Berries (1/2 cup)"] = max(matrix["Berries (1/2 cup)"], 7)
        matrix["Processed/cured meat"] = 0
    if "Slow" in (cof.get("metabolism") or ""):
        matrix["Coffee (8 oz)"] = 4
    if alc.get("risk") in ("Avoid entirely", "Strongly limit"):
        matrix["Alcohol (standard drink)"] = 0

    return {"servings_per_week": matrix}


# ── Concrete Shopping List ──────────────────────────────────────────────────

def shopping_list(matrix: dict) -> list[dict]:
    sv = matrix["servings_per_week"]
    items = [
        {"category": "Proteins",
         "items": [
             f"{sv['Oily fish (salmon/sardines/mackerel)'] * 5} oz wild salmon/sardines/mackerel",
             f"{sv['Lean protein (chicken/turkey/lean beef)'] * 5} oz chicken/turkey breast",
             f"{sv['Whole eggs']} large eggs",
             f"{sv['Legumes (beans/lentils/chickpeas)']} cans (15 oz) of mixed beans/lentils/chickpeas",
             ("Tofu/tempeh block (14 oz)" if sv["Tofu/tempeh"] else "—"),
         ]},
        {"category": "Dairy / alternatives",
         "items": [
             f"{sv['Greek yoghurt / kefir']} × 6 oz Greek yoghurt or kefir",
             f"{sv['Aged cheese (matchbox-sized)'] * 1} oz aged cheese (parmesan/cheddar/manchego)",
         ]},
        {"category": "Produce — vegetables",
         "items": [
             f"{sv['Leafy greens (1 cup cooked or 2 raw)']} cups spinach/kale/arugula",
             f"{sv['Cruciferous veg (broccoli/cauliflower/Brussels)']} crowns broccoli OR Brussels sprouts bag",
             f"{sv['Allium veg (garlic/onion/leek)']} bulbs garlic / large onions",
             f"{sv['Colorful veg (bell pepper/carrot/beet)']} mixed bell peppers / carrots / beets",
         ]},
        {"category": "Produce — fruit",
         "items": [
             f"{sv['Berries (1/2 cup)']} cups frozen mixed berries (cheapest year-round)",
             f"{sv['Whole fruit (apple/pear/citrus)']} pieces of fruit (apples, pears, citrus)",
         ]},
        {"category": "Nuts / seeds / fats",
         "items": [
             f"{sv['Nuts (1 oz)']} oz mixed walnuts/almonds (raw, unsalted)",
             f"{sv['Seeds (chia/flax/pumpkin, 1 tbsp)']} tbsp ground flaxseed + chia jar",
             f"{sv['Extra-virgin olive oil (tbsp)']} tbsp extra-virgin olive oil (~{sv['Extra-virgin olive oil (tbsp)']//4} oz)",
         ]},
        {"category": "Grains / staples",
         "items": [
             f"{sv['Whole grains (oats/quinoa/brown rice)']} servings: rolled oats, quinoa, brown rice",
         ]},
        {"category": "Beverages / extras",
         "items": [
             f"Coffee — {sv['Coffee (8 oz)']} cups/week",
             f"Green tea — {sv['Green tea (cups)']} cups/week",
             (f"Dark chocolate (≥75% cacao) — {sv['Dark chocolate (≥75%, 1 oz)']} oz" if sv['Dark chocolate (≥75%, 1 oz)'] else "—"),
             "Turmeric (ground), black pepper, ginger, garlic powder, oregano, cinnamon",
         ]},
    ]
    return items


# ── Generated Recipes ───────────────────────────────────────────────────────

def generate_recipes(result_partial: dict) -> list[dict]:
    low_carb = result_partial["macros"]["pct_carbs"] <= 35
    e4 = result_partial.get("saturated_fat", {}).get("apoe_genotype", "").startswith("ε4")
    high_omega3 = result_partial.get("omega3", {}).get("ala_conversion") == "Poor"
    iron_overload = result_partial.get("iron", {}).get("overload_risk", "").startswith("High") or \
                    "Moderate" in result_partial.get("iron", {}).get("overload_risk", "")
    lactose = result_partial.get("lactose", {}).get("tolerance", "").startswith("Intolerant")
    gluten = result_partial.get("gluten", {}).get("celiac_risk_haplotype", False)
    bitter = result_partial.get("taste", {}).get("bitter", "") == "Super-taster"
    appetite = result_partial.get("satiety", {}).get("appetite_phenotype", "Standard") != "Standard"

    # Recipe 1 — Breakfast bowl (adjusts for ApoE, lactose, gluten, choline)
    breakfast = {
        "name": "High-Protein Breakfast Bowl",
        "macros_est": "~30 g protein · 35 g carbs · 18 g fat · 12 g fibre" if not low_carb
                      else "~35 g protein · 12 g carbs · 28 g fat · 8 g fibre",
        "ingredients": [
            ("3 large eggs" if not e4 else "1 whole egg + 3 egg whites"),
            ("1 cup steel-cut oats (cooked)" if not low_carb and not gluten
             else ("1 cup gluten-free oats (cooked)" if not low_carb else "½ cup riced cauliflower sautéed")),
            ("1 cup Greek yoghurt" if not lactose else "1 cup coconut/almond yoghurt"),
            "½ cup mixed berries (blueberries + raspberries)",
            "1 tbsp ground flaxseed" + (" + 1 tbsp chia" if high_omega3 else ""),
            "1 tbsp chopped walnuts",
            ("1 tbsp extra-virgin olive oil drizzle" if e4 else "¼ avocado"),
            "Pinch cinnamon + dash of vanilla",
        ],
        "why_for_you": [
            "Eggs supply choline (PEMT)" if not e4 else "Reduced whole-egg dose limits saturated fat (APOE ε4)",
            "Flax+chia give ALA for FADS poor converters" if high_omega3 else "Flax adds fibre + ALA",
            "Greek yoghurt = aged-fermented, lower-lactose protein" if lactose else "Greek yoghurt = high-protein satiety anchor",
            "Berries deliver anthocyanins (Nrf2 antioxidant pathway)",
        ],
        "method": (
            "1) Cook oats with water + cinnamon. 2) Scramble eggs in olive oil. "
            "3) Assemble bowl: oats base, top with yoghurt, eggs alongside, scatter "
            "berries, flax, chia, walnuts. 4) Drizzle olive oil. Eat within 30 min."
        ),
    }

    # Recipe 2 — Anti-inflammatory dinner sheet-pan
    protein = "1 lb wild salmon fillets" if high_omega3 else "1 lb chicken thighs (boneless)"
    if iron_overload:
        protein = "1 lb wild salmon fillets" if not high_omega3 else "1 lb wild cod or salmon"
    crucifer_treatment = ("roasted hard at 425 °F with olive oil, balsamic + honey glaze "
                         "(super-taster bitter mitigation)" if bitter else "roasted at 400 °F with olive oil")
    dinner = {
        "name": "Mediterranean Sheet-Pan Dinner",
        "macros_est": "~40 g protein · 40 g carbs · 22 g fat · 14 g fibre"
                      if not low_carb else "~45 g protein · 18 g carbs · 30 g fat · 11 g fibre",
        "ingredients": [
            protein,
            "1 lb broccoli + Brussels sprouts (Phase II / cruciferous)",
            ("1 large sweet potato cubed" if not low_carb else "1 cup cauliflower rice"),
            "1 red onion + 1 lemon (sliced)",
            "3 tbsp extra-virgin olive oil",
            "4 garlic cloves minced + 1 tsp turmeric + ½ tsp black pepper + 1 tsp oregano",
            "1 tbsp tahini + 2 tbsp lemon juice + water (drizzle sauce)",
        ],
        "why_for_you": [
            ("Salmon = direct EPA+DHA (FADS poor converter)" if high_omega3 else "Lean protein"),
            "Cruciferous + " + crucifer_treatment + " (sulforaphane induces Phase II detox enzymes)",
            "Turmeric + pepper + olive oil: 2000% bioavailability boost for curcumin",
            ("ApoE ε4 friendly — fish + olive oil, no butter" if e4 else "Olive-oil base supports HDL"),
        ],
        "method": (
            "1) Preheat oven 425 °F. 2) Toss cruciferous + sweet potato + onion in olive oil "
            "+ spices on sheet pan; roast 18 min. 3) Add protein + lemon slices; roast 12 "
            "more min (salmon) or 18 (chicken). 4) Whisk tahini + lemon + water. 5) Plate, "
            "drizzle tahini sauce, finish with fresh parsley."
        ),
    }

    # Recipe 3 — High-satiety lunch (appetite-elevated genotype)
    lunch = {
        "name": ("Volume-Density Power Salad" if appetite else "Mediterranean Bowl"),
        "macros_est": "~35 g protein · 30 g carbs · 22 g fat · 16 g fibre",
        "ingredients": [
            "2 cups mixed greens (arugula + spinach + romaine)",
            ("1 cup chickpeas (rinsed)" if not iron_overload else "1 cup white beans"),
            "5 oz grilled chicken or canned tuna (water-packed)",
            "½ cucumber + 1 cup cherry tomatoes + ¼ red onion",
            ("¼ avocado" if e4 else "¼ avocado + 1 oz feta"),
            "2 tbsp pumpkin seeds",
            "Dressing: 2 tbsp olive oil + 1 tbsp lemon juice + Dijon + garlic",
        ],
        "why_for_you": [
            ("FTO/MC4R appetite genotype — volume + fibre + protein density" if appetite
             else "Balanced satiating lunch"),
            ("Beans replace chickpeas to keep iron load down" if iron_overload else "Chickpeas: fibre + protein + folate"),
            "Avocado MUFA + olive oil for HDL maintenance",
        ],
        "method": (
            "1) Toss greens with dressing in large bowl. 2) Layer beans/chickpeas, "
            "cucumber, tomatoes, onion, avocado. 3) Top with protein + seeds. 4) Mix "
            "thoroughly so flavour reaches every bite (eating speed matters — aim 20 min)."
        ),
    }

    return [breakfast, lunch, dinner]


# ── Public synthesis API ────────────────────────────────────────────────────

def analyze_advanced_nutrition(snps_df, base_result: dict) -> dict:
    """Run all advanced analyzers and return a unified extension dict."""
    pgs = analyze_polygenic_scores(snps_df)
    dashboard = cardiometabolic_dashboard(pgs)
    inflam = inflammation_index(snps_df)
    histamine = histamine_intolerance(snps_df)
    detox = detoxification_profile(snps_df)
    glycemic = glycemic_threshold(snps_df)
    macros = base_result["macros"]
    satiety = base_result.get("satiety", {"appetite_phenotype": "Standard"})
    periodisation = macro_periodisation(macros, glycemic, satiety)

    partial_for_matrix = dict(base_result)
    partial_for_matrix["inflammation"] = inflam
    partial_for_matrix["detoxification"] = detox
    matrix = weekly_food_matrix(partial_for_matrix)
    shop = shopping_list(matrix)
    recipes = generate_recipes(partial_for_matrix)

    return {
        "polygenic_scores": pgs,
        "cardiometabolic_dashboard": dashboard,
        "inflammation": inflam,
        "histamine": histamine,
        "detoxification": detox,
        "glycemic_threshold": glycemic,
        "macro_periodisation": periodisation,
        "weekly_food_matrix": matrix,
        "shopping_list": shop,
        "recipes": recipes,
    }
