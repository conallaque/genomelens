"""
Personalized Nutrition Plan
===========================

Genotype-driven dietary prescription:

  • macro ratios (carb / fat / protein) — FTO, TCF7L2, PPARG, APOE, FADS
  • foods to emphasize / avoid
  • caffeine guidance — CYP1A2, ADORA2A
  • alcohol guidance — ALDH2, ADH1B
  • salt sensitivity — ACE, AGT
  • lactose tolerance — LCT (rs4988235)
  • gluten / celiac risk — HLA-DQ2/8 tag SNPs
  • methylation diet — MTHFR
  • vitamin D from food — VDR / CYP2R1

Output dict:
  {
    status,
    macros: {pct_carbs, pct_fat, pct_protein, rationale},
    emphasize: [...], avoid: [...],
    caffeine: {...}, alcohol: {...}, salt: {...},
    lactose: {...}, gluten: {...},
    methylation: {...}, vitamin_d_food: {...},
    daily_template: [...],
  }
"""

from __future__ import annotations

import pandas as pd

try:
    from .nutrition_advanced import analyze_advanced_nutrition
except ImportError:
    analyze_advanced_nutrition = None
try:
    from .nutrition_protocols import analyze_nutrition_protocols
except ImportError:
    analyze_nutrition_protocols = None


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


# ── Macronutrient ratios ────────────────────────────────────────────────────

def _analyze_macros(snps_df) -> dict:
    fto = _gt(snps_df, "rs9939609")       # A = obesity risk
    tcf7l2 = _gt(snps_df, "rs7903146")     # T = T2D risk + worse carb response
    ppara = _gt(snps_df, "rs1800206")      # PPARA — fat metabolism
    apoe1 = _gt(snps_df, "rs429358")       # APOE e4 risk SNP
    apoe2 = _gt(snps_df, "rs7412")
    fads = _gt(snps_df, "rs174547")

    carb_pressure = 0
    fat_pressure = 0
    factors: list[str] = []

    if fto:
        if "A" in fto:
            carb_pressure -= 1
            factors.append(f"rs9939609 (FTO) {fto} — better response to lower-carb diet")
        else:
            factors.append(f"rs9939609 (FTO) {fto} — favorable FTO; carb-tolerant")

    if tcf7l2 and "T" in tcf7l2:
        carb_pressure -= 2
        factors.append(f"rs7903146 (TCF7L2) {tcf7l2} — reduced glucose tolerance; minimize refined carbs")

    # rs429358 C + rs7412 C = ε4 carrier — saturated-fat sensitive
    if apoe1 and apoe2 and "C" in apoe1 and "C" in apoe2:
        fat_pressure -= 1
        factors.append("APOE ε4 carrier — lower saturated fat, prefer mono/poly")

    if fads and "T" in fads:
        factors.append(f"rs174547 (FADS1) {fads} — favor direct EPA/DHA from fish, not ALA")

    if ppara and "G" in ppara:
        factors.append(f"rs1800206 (PPARA) {ppara} — wild-type fat metabolism")

    # Base 45/30/25 (C/F/P), then shift
    pct_carbs = 45 + (carb_pressure * 5)
    pct_fat = 30 + (fat_pressure * -3)   # if fat pressure negative, drop fat slightly
    pct_protein = 100 - pct_carbs - pct_fat

    # Clamp
    pct_carbs = max(20, min(55, pct_carbs))
    pct_fat = max(20, min(45, pct_fat))
    pct_protein = max(15, min(35, 100 - pct_carbs - pct_fat))
    # Re-balance to 100
    leftover = 100 - (pct_carbs + pct_fat + pct_protein)
    pct_carbs += leftover

    if pct_carbs <= 35:
        rationale = "Lower-carb skew (~30%) — your variants reduce glucose tolerance."
    elif pct_carbs >= 50:
        rationale = "Carb-tolerant — moderate-to-high quality complex carbs work well."
    else:
        rationale = "Balanced macros — no strong genetic pressure toward keto or high-carb."

    return {
        "pct_carbs": pct_carbs,
        "pct_fat": pct_fat,
        "pct_protein": pct_protein,
        "rationale": rationale,
        "factors": factors or ["No macro-relevant SNPs typed"],
    }


# ── Caffeine ────────────────────────────────────────────────────────────────

def _analyze_caffeine(snps_df) -> dict:
    cyp1a2 = _gt(snps_df, "rs762551")     # A allele = fast metaboliser
    adora2a = _gt(snps_df, "rs5751876")    # T allele = anxiety-prone with caffeine

    fast = cyp1a2 and "A" in cyp1a2 and cyp1a2.count("A") == 2
    slow = cyp1a2 and "C" in cyp1a2

    anxiety = adora2a and "T" in adora2a

    if not cyp1a2:
        return {
            "metabolism": "Unknown",
            "limit_mg": 400,
            "cutoff_time": "14:00",
            "factors": ["CYP1A2 not typed"],
            "guidance": "Use standard 400 mg/day adult limit; stop caffeine after 14:00.",
        }

    if fast:
        return {
            "metabolism": "Fast (*1A/*1A)",
            "limit_mg": 400,
            "cutoff_time": "16:00",
            "factors": [f"rs762551 (CYP1A2) {cyp1a2}"] +
                       ([f"rs5751876 (ADORA2A) {adora2a} — anxiety-prone"] if anxiety else []),
            "guidance": (
                "Fast metaboliser — caffeine cleared in ~4-5 h. Up to 400 mg/day safe; "
                "coffee/CV-disease risk is neutral or protective for you."
                + (" ADORA2A T-allele present — drop dose if jittery." if anxiety else "")
            ),
        }
    if slow:
        return {
            "metabolism": "Slow (*1F carrier)",
            "limit_mg": 200,
            "cutoff_time": "12:00",
            "factors": [f"rs762551 (CYP1A2) {cyp1a2} — slow"] +
                       ([f"rs5751876 (ADORA2A) {adora2a} — anxiety-prone"] if anxiety else []),
            "guidance": (
                "Slow metaboliser — caffeine half-life prolonged (8+ h). Cap at 200 mg/day "
                "(≈ 1 large coffee) and avoid after lunch. Higher dose links to elevated "
                "CV risk in slow metabolisers."
                + (" L-Theanine 100-200 mg paired with caffeine smooths anxiety response."
                   if anxiety else "")
            ),
        }
    return {
        "metabolism": "Intermediate",
        "limit_mg": 300,
        "cutoff_time": "13:00",
        "factors": [f"rs762551 (CYP1A2) {cyp1a2}"],
        "guidance": "Intermediate metaboliser — cap at ~300 mg/day; stop by early afternoon.",
    }


# ── Alcohol ─────────────────────────────────────────────────────────────────

def _analyze_alcohol(snps_df) -> dict:
    aldh2 = _gt(snps_df, "rs671")          # A allele (East Asian variant) → flushing
    adh1b = _gt(snps_df, "rs1229984")      # A allele → fast ethanol→acetaldehyde
    factors: list[str] = []
    risk = "Standard"

    if aldh2 and "A" in aldh2:
        if aldh2.count("A") == 2:
            risk = "Avoid entirely"
            factors.append("rs671 (ALDH2*2/*2) — non-functional ALDH2; acetaldehyde toxic accumulation")
        else:
            risk = "Strongly limit"
            factors.append("rs671 (ALDH2*1/*2) — partial deficiency; flushing, cancer-risk elevated")

    if adh1b and "A" in adh1b:
        factors.append("rs1229984 (ADH1B*2) — fast acetaldehyde production; aldehyde load higher")
        if risk == "Standard":
            risk = "Reduce"

    if risk == "Avoid entirely":
        guidance = (
            "ALDH2 homozygous variant — alcohol intake elevates oesophageal cancer "
            "risk dramatically. Strongly recommend abstinence."
        )
    elif risk == "Strongly limit":
        guidance = (
            "ALDH2 heterozygote — flushing reaction is a warning sign of acetaldehyde "
            "accumulation. WHO classifies this with significantly elevated cancer risk. "
            "Maximum 1-2 standard drinks per week, ideally none."
        )
    elif risk == "Reduce":
        guidance = "ADH1B fast-conversion genotype — keep intake light, hydrate well."
    else:
        guidance = "No alcohol-risk variants typed — standard moderation applies (≤1-2 drinks/day for men, ≤1 for women)."

    return {"risk": risk, "factors": factors or ["No alcohol-related SNPs typed"],
            "guidance": guidance}


# ── Salt sensitivity ────────────────────────────────────────────────────────

def _analyze_salt(snps_df) -> dict:
    ace = _gt(snps_df, "rs4341") or _gt(snps_df, "rs4646994")
    agt = _gt(snps_df, "rs699")         # AGT M235T; T = salt-sensitive
    factors: list[str] = []
    sensitive = False
    if agt and "G" in agt:
        sensitive = True
        factors.append(f"rs699 (AGT M235T) {agt} — salt-sensitive hypertension allele")
    if ace and "T" in ace and ace.count("T") == 2:
        sensitive = True
        factors.append(f"ACE D/D ({ace}) — sodium-sensitive BP response")
    return {
        "sensitive": sensitive,
        "factors": factors or ["No salt-sensitivity SNPs typed"],
        "guidance": (
            "Reduce sodium to ≤2.3 g/day; emphasize potassium-rich foods (leafy greens, "
            "beans, potato). DASH-style pattern especially beneficial."
            if sensitive else
            "No genetic salt-sensitivity flagged — keep sodium within general "
            "guidelines (≤2.3 g/day adult)."
        ),
    }


# ── Lactose ─────────────────────────────────────────────────────────────────

def _analyze_lactose(snps_df) -> dict:
    lct = _gt(snps_df, "rs4988235")        # T = persistence; CC = intolerant
    if not lct:
        return {"tolerance": "Unknown", "factors": ["LCT rs4988235 not typed"],
                "guidance": "Lactase-persistence variant not typed."}
    if "T" in lct:
        return {
            "tolerance": "Persistent",
            "factors": [f"rs4988235 (LCT) {lct}"],
            "guidance": (
                "Lactase persistence — dairy is well tolerated lifelong. Greek yoghurt, "
                "kefir, and aged cheeses are excellent protein/probiotic sources."
            ),
        }
    return {
        "tolerance": "Intolerant (non-persistent)",
        "factors": [f"rs4988235 (LCT) {lct}"],
        "guidance": (
            "Lactase non-persistence — symptomatic lactose intolerance likely. Choose "
            "lactose-free milk, aged cheeses (parmesan, cheddar), fermented dairy "
            "(yoghurt, kefir) which are partially pre-digested. Plant milks for liquid dairy."
        ),
    }


# ── Gluten / DQ2-DQ8 risk ───────────────────────────────────────────────────

def _analyze_gluten(snps_df) -> dict:
    dq2 = _gt(snps_df, "rs2187668")        # HLA-DQ2.5 tag
    dq8 = _gt(snps_df, "rs7454108")        # HLA-DQ8 tag
    carrier = (dq2 and "T" in dq2) or (dq8 and "C" in dq8)
    factors: list[str] = []
    if dq2:
        factors.append(f"rs2187668 (HLA-DQ2 tag) {dq2}")
    if dq8:
        factors.append(f"rs7454108 (HLA-DQ8 tag) {dq8}")
    return {
        "celiac_risk_haplotype": bool(carrier),
        "factors": factors or ["DQ2/DQ8 tag SNPs not typed"],
        "guidance": (
            "DQ2/DQ8 carrier — ~3% lifetime risk of developing celiac disease (vs <0.1% "
            "without). If you have GI symptoms, fatigue, or family history, ask your "
            "physician for serology (tTG-IgA) before starting gluten-free diet."
            if carrier else
            "No DQ2/DQ8 risk haplotype typed/detected — celiac disease is essentially "
            "ruled out. Gluten avoidance unnecessary for autoimmune reasons."
        ),
    }


# ── Methylation diet ────────────────────────────────────────────────────────

def _analyze_methylation_diet(snps_df) -> dict:
    mthfr_c677t = _gt(snps_df, "rs1801133")
    mthfr_a1298c = _gt(snps_df, "rs1801131")
    factors: list[str] = []
    needs_extra = False
    if mthfr_c677t:
        factors.append(f"rs1801133 (MTHFR C677T) {mthfr_c677t}")
        if "T" in mthfr_c677t or "A" in mthfr_c677t:
            needs_extra = True
    if mthfr_a1298c:
        factors.append(f"rs1801131 (MTHFR A1298C) {mthfr_a1298c}")
        if mthfr_a1298c.count("G") >= 1 or mthfr_a1298c.count("C") >= 1:
            needs_extra = True
    return {
        "needs_methylation_support": needs_extra,
        "factors": factors or ["MTHFR not typed"],
        "guidance": (
            "Emphasize folate-rich foods: leafy greens (spinach, kale, romaine), liver, "
            "lentils, asparagus, broccoli. Avoid synthetic folic-acid fortified products "
            "(many cereals, breads) — they compete with active folate. Choose "
            "methylfolate-supplemented options if available."
            if needs_extra else
            "Standard varied diet supplies sufficient folate; aim for 1-2 cups leafy "
            "greens daily."
        ),
    }


# ── Vitamin D from food ─────────────────────────────────────────────────────

def _analyze_vitamin_d_food(snps_df) -> dict:
    cyp = _gt(snps_df, "rs10741657")
    gc = _gt(snps_df, "rs2282679")
    vdr_fok = _gt(snps_df, "rs2228570")
    factors = []
    needs_more = False
    if cyp and "G" in cyp:
        factors.append(f"rs10741657 (CYP2R1) {cyp}")
        needs_more = True
    if gc and "C" in gc:
        factors.append(f"rs2282679 (GC/VDBP) {gc}")
        needs_more = True
    if vdr_fok and "T" in vdr_fok:
        factors.append(f"rs2228570 (VDR FokI) {vdr_fok}")
    return {
        "needs_more_intake": needs_more,
        "factors": factors or ["Vitamin D SNPs not typed"],
        "guidance": (
            "Prioritize high-D foods: wild salmon, sardines, herring, egg yolks, "
            "UV-treated mushrooms, fortified dairy. Combine with daytime sun exposure "
            "where possible (10-30 min, depending on skin tone, latitude, season)."
            if needs_more else
            "Modest D intake from fatty fish 2×/week + occasional sun is sufficient."
        ),
    }


# ── Foods to emphasize / avoid (synthesised) ───────────────────────────────

def _build_food_lists(macros: dict, alcohol: dict, lactose: dict, salt: dict, gluten: dict) -> dict:
    emphasize: list[str] = []
    avoid: list[str] = []

    if macros["pct_carbs"] <= 35:
        avoid.extend(["Sugary drinks", "White bread / pastries", "Sweetened breakfast cereals"])
        emphasize.extend(["Non-starchy vegetables", "Berries", "Quinoa / oats (modest portions)"])
    else:
        emphasize.extend(["Whole grains (oats, brown rice, barley)", "Legumes", "Sweet potato"])

    emphasize.extend(["Fatty fish (salmon, sardines) 2-3×/week",
                      "Olive oil (extra-virgin) — primary cooking fat",
                      "Leafy greens daily (folate + nitrates)",
                      "Nuts & seeds (small handful daily)",
                      "Berries (antioxidant load)"])

    if alcohol["risk"] in ("Avoid entirely", "Strongly limit"):
        avoid.append("Alcohol (genetic ALDH2 / ADH1B contraindication)")
    if salt["sensitive"]:
        avoid.extend(["Processed meats / deli", "Salted snacks", "High-sodium soups & sauces"])
    if lactose["tolerance"] == "Intolerant (non-persistent)":
        avoid.append("Fresh milk / ice cream (lactose intolerance)")
        emphasize.append("Lactose-free dairy, aged cheeses, kefir/yoghurt")
    if gluten["celiac_risk_haplotype"]:
        avoid.append("Standard gluten products if symptomatic (DQ2/DQ8 carrier)")

    return {"emphasize": emphasize, "avoid": avoid}


# ── Daily template ──────────────────────────────────────────────────────────

def _build_daily_template(macros: dict, caffeine: dict, alcohol: dict) -> list[dict]:
    low_carb = macros["pct_carbs"] <= 35
    return [
        {"meal": "Breakfast", "example": (
            "3 eggs + sautéed spinach + ½ avocado + black coffee"
            if low_carb else
            "Steel-cut oats with berries + walnuts + Greek yoghurt + coffee"
        )},
        {"meal": "Mid-morning", "example": "Almonds + apple OR 1 boiled egg"},
        {"meal": "Lunch", "example": (
            "Grilled salmon + large mixed-greens salad + olive oil + 1/2 sweet potato"
        )},
        {"meal": "Afternoon", "example": (
            f"Light snack — {'Greek yoghurt or hummus + veggies' if not low_carb else 'cheese + olives'}; "
            f"caffeine cutoff {caffeine['cutoff_time']}"
        )},
        {"meal": "Dinner", "example": (
            "Lean protein (chicken, turkey, or tofu) + cruciferous vegetables (broccoli, "
            "Brussels) + " + ("modest quinoa or wild rice" if not low_carb else "extra olive oil")
        )},
        {"meal": "Evening", "example": (
            "Herbal tea (chamomile/ginger). Avoid alcohol per genotype recommendation."
            if alcohol["risk"] in ("Avoid entirely", "Strongly limit") else
            "Optional: small glass of red wine; finish ≥2 h before sleep."
        )},
    ]


# ── Omega-3 (FADS1/FADS2 ALA→EPA/DHA conversion) ────────────────────────────

def _analyze_omega3(snps_df) -> dict:
    fads1 = _gt(snps_df, "rs174547")
    fads2 = _gt(snps_df, "rs174537") or _gt(snps_df, "rs174575")
    elovl2 = _gt(snps_df, "rs953413")
    factors: list[str] = []
    poor = False
    if fads1:
        factors.append(f"rs174547 (FADS1) {fads1}")
        if "T" in fads1:
            poor = True
    if fads2:
        factors.append(f"rs174537/rs174575 (FADS2) {fads2}")
        if "T" in fads2 or "G" in fads2:
            poor = True
    if elovl2:
        factors.append(f"rs953413 (ELOVL2) {elovl2}")
        if "G" in elovl2:
            poor = True
    if poor:
        return {
            "ala_conversion": "Poor",
            "epa_dha_target_mg": 1500,
            "factors": factors,
            "guidance": (
                "Reduced-activity FADS desaturase variant — ALA from flax/chia/walnut "
                "converts inefficiently (<5%) to EPA/DHA. Hit 1.5–2 g combined EPA+DHA "
                "daily from oily fish (salmon, sardines, mackerel, herring) 3–4×/week, "
                "or algal-/fish-oil supplementation. Vegetarians: algae-oil 500–1000 mg "
                "EPA+DHA daily."
            ),
        }
    return {
        "ala_conversion": "Normal" if factors else "Unknown",
        "epa_dha_target_mg": 1000,
        "factors": factors or ["FADS SNPs not typed"],
        "guidance": (
            "Normal desaturase activity — plant ALA sources (flax, chia, walnut) "
            "contribute meaningfully. Still aim for ~1 g EPA+DHA from fish 2×/week."
        ),
    }


# ── Iron overload risk (hemochromatosis) ────────────────────────────────────

def _analyze_iron(snps_df) -> dict:
    c282y = _gt(snps_df, "rs1800562")
    h63d = _gt(snps_df, "rs1799945")
    tmprss6 = _gt(snps_df, "rs855791")
    factors: list[str] = []
    overload = "Low"
    if c282y:
        factors.append(f"rs1800562 (HFE C282Y) {c282y}")
        if "A" in c282y and c282y.count("A") == 2:
            overload = "High (homozygous)"
        elif "A" in c282y:
            overload = "Moderate (heterozygous)"
    if h63d:
        factors.append(f"rs1799945 (HFE H63D) {h63d}")
        if "G" in h63d:
            if overload == "Low":
                overload = "Mild"
            elif "Moderate" in overload:
                overload = "Compound heterozygous"
    if tmprss6:
        factors.append(f"rs855791 (TMPRSS6) {tmprss6}")
    if overload.startswith("High"):
        guidance = (
            "HFE C282Y homozygous — hereditary hemochromatosis genotype. ~28% penetrance "
            "for iron overload in men, lower in women. Action: ask physician for serum "
            "ferritin + transferrin saturation. Limit red meat to 1×/week, avoid iron-"
            "fortified cereals, do NOT take iron-containing multivitamins, avoid "
            "vitamin-C megadoses with meals (enhances iron absorption), moderate alcohol. "
            "Tea/coffee with meals reduces iron absorption — useful."
        )
    elif "Moderate" in overload or "Compound" in overload:
        guidance = (
            "Heterozygous/compound HFE — modest overload risk. Avoid iron-fortified "
            "supplements unless ferritin documented low. Monitor ferritin every 2–3 yr."
        )
    elif "Mild" in overload:
        guidance = "Single H63D copy — minimal clinical risk; no action needed."
    else:
        guidance = (
            "No hemochromatosis risk variants typed/detected. Standard iron intake; "
            "vegetarians should pair plant iron with vitamin C; menstruating women "
            "may need extra (18 mg/day RDA)."
        )
    return {"overload_risk": overload, "factors": factors or ["HFE not typed"],
            "guidance": guidance}


# ── Choline (PEMT) ──────────────────────────────────────────────────────────

def _analyze_choline(snps_df) -> dict:
    pemt = _gt(snps_df, "rs7946")
    mthfd1 = _gt(snps_df, "rs2236225")
    factors: list[str] = []
    needs = False
    if pemt:
        factors.append(f"rs7946 (PEMT) {pemt}")
        if "T" in pemt:
            needs = True
    if mthfd1:
        factors.append(f"rs2236225 (MTHFD1) {mthfd1}")
        if "A" in mthfd1:
            needs = True
    return {
        "increased_need": needs,
        "target_mg": 550 if needs else 425,
        "factors": factors or ["PEMT/MTHFD1 not typed"],
        "guidance": (
            "Reduced endogenous phosphatidylcholine synthesis — dietary choline becomes "
            "essential. Target 550 mg/day from: 2 whole eggs (≈250 mg), beef liver (1 oz "
            "≈100 mg), soybeans, salmon, chicken, cruciferous veg. Inadequate choline "
            "→ NAFLD risk."
            if needs else
            "Adequate intake target 425 mg/day (women) / 550 mg/day (men). 1–2 whole "
            "eggs daily covers most."
        ),
    }


# ── Vitamin B12 (FUT2 secretor status, TCN2) ────────────────────────────────

def _analyze_b12(snps_df) -> dict:
    fut2 = _gt(snps_df, "rs601338")
    tcn2 = _gt(snps_df, "rs1801198")
    factors: list[str] = []
    lower_status = False
    if fut2:
        factors.append(f"rs601338 (FUT2 secretor) {fut2}")
        # AA = non-secretor — higher serum B12 paradoxically but altered gut absorption
        if "A" in fut2 and fut2.count("A") == 2:
            lower_status = True
    if tcn2:
        factors.append(f"rs1801198 (TCN2) {tcn2}")
        if "G" in tcn2:
            lower_status = True
    return {
        "lower_functional_b12": lower_status,
        "factors": factors or ["B12 transport SNPs not typed"],
        "guidance": (
            "Variants in B12 transport/secretor status — prioritize high-bioavailable "
            "B12 sources: clams, beef liver, sardines, eggs, dairy. Consider 500 µg "
            "methylcobalamin sublingual 2–3×/week if vegetarian/vegan."
            if lower_status else
            "Standard B12 intake from animal foods (meat, fish, dairy, eggs) suffices. "
            "Strict vegans should supplement 250–500 µg cyanocobalamin daily."
        ),
    }


# ── Vitamin A conversion (BCO1) ─────────────────────────────────────────────

def _analyze_vitamin_a(snps_df) -> dict:
    bco1a = _gt(snps_df, "rs7501331")
    bco1b = _gt(snps_df, "rs12934922")
    factors: list[str] = []
    poor = False
    if bco1a:
        factors.append(f"rs7501331 (BCO1) {bco1a}")
        if "T" in bco1a:
            poor = True
    if bco1b:
        factors.append(f"rs12934922 (BCO1) {bco1b}")
        if "A" in bco1b:
            poor = True
    return {
        "beta_carotene_converter": "Poor" if poor else "Normal",
        "factors": factors or ["BCO1 not typed"],
        "guidance": (
            "BCO1 reduced-activity — β-carotene from carrots/sweet potato converts "
            "poorly to retinol (~30% of normal). Include preformed retinol weekly: "
            "egg yolks, dairy, liver (1 oz/week), oily fish. Cooking + fat with "
            "carotenoid veg improves absorption."
            if poor else
            "Normal β-carotene conversion — colorful plant sources (sweet potato, "
            "carrot, kale, pumpkin) cover retinol needs."
        ),
    }


# ── Vitamin C, E ────────────────────────────────────────────────────────────

def _analyze_vitamin_c(snps_df) -> dict:
    slc = _gt(snps_df, "rs33972313") or _gt(snps_df, "rs6596473")
    factors: list[str] = []
    higher_need = False
    if slc:
        factors.append(f"SLC23A1/A2 {slc}")
        if "A" in slc or "T" in slc:
            higher_need = True
    return {
        "higher_need": higher_need,
        "target_mg": 200 if higher_need else 90,
        "factors": factors or ["SLC23A SNPs not typed"],
        "guidance": (
            "Higher vitamin-C requirement — aim 200 mg/day from food: 1 red bell "
            "pepper (~150 mg), kiwi (~70 mg), broccoli, citrus, strawberries. Split "
            "across the day (absorption saturates ~200 mg per dose)."
            if higher_need else
            "Standard 90 mg/day (men)/75 mg (women) from varied produce."
        ),
    }


def _analyze_vitamin_e(snps_df) -> dict:
    cyp4f2 = _gt(snps_df, "rs2108622")
    factors: list[str] = []
    higher = False
    if cyp4f2:
        factors.append(f"rs2108622 (CYP4F2) {cyp4f2}")
        if "T" in cyp4f2:
            higher = True
    return {
        "higher_retention": higher,
        "factors": factors or ["CYP4F2 not typed"],
        "guidance": (
            "CYP4F2 reduced metabolism — vitamin E (α-tocopherol) accumulates more. "
            "AVOID high-dose vitamin-E supplements (>200 IU/day). Get from food: "
            "almonds, sunflower seeds, avocado, olive oil."
            if higher else
            "Standard vitamin-E from nuts/seeds/oils. Supplementation usually unnecessary."
        ),
    }


# ── Taste perception (TAS2R38 bitter, CD36 fat, sweet preference) ───────────

def _analyze_taste(snps_df) -> dict:
    tas2r38 = _gt(snps_df, "rs713598")          # G = PAV taster, C = AVI non-taster
    cd36 = _gt(snps_df, "rs1761667")            # A = poor fat-taste sensitivity
    factors: list[str] = []
    bitter = "Unknown"
    fat_taste = "Unknown"
    if tas2r38:
        factors.append(f"rs713598 (TAS2R38) {tas2r38}")
        if "G" in tas2r38 and tas2r38.count("G") == 2:
            bitter = "Super-taster"
        elif "G" in tas2r38:
            bitter = "Taster"
        else:
            bitter = "Non-taster"
    if cd36:
        factors.append(f"rs1761667 (CD36) {cd36}")
        if "A" in cd36 and cd36.count("A") == 2:
            fat_taste = "Reduced fat perception (overconsumption risk)"
        else:
            fat_taste = "Normal fat perception"
    bitter_advice = {
        "Super-taster": (
            "PAV/PAV super-taster — cruciferous veg (Brussels, kale, broccoli) and "
            "coffee taste markedly bitter. Mitigate with: roasting (caramelises), "
            "olive-oil/lemon dressings, blanching, adding modest sweetness (balsamic, "
            "honey-glaze). Don't skip these vegetables — they remain critical for "
            "phytonutrient intake."
        ),
        "Taster": "Mild bitter sensitivity — most vegetables palatable with light seasoning.",
        "Non-taster": "Bitter compounds barely register — bonus: black coffee, dark chocolate, kale taste fine.",
        "Unknown": "Bitter-taste variant not typed.",
    }
    fat_advice = (
        "Reduced fatty-acid oral detection — risk of unconsciously overeating high-fat "
        "foods. Pre-portion calorie-dense items (nuts, oils, cheese) by weight instead "
        "of estimating."
        if "Reduced" in fat_taste else
        "Normal fat-taste — internal satiety cues for high-fat foods are reliable."
    )
    return {
        "bitter": bitter, "fat_taste": fat_taste,
        "factors": factors or ["Taste-perception SNPs not typed"],
        "bitter_guidance": bitter_advice[bitter],
        "fat_guidance": fat_advice,
    }


# ── Satiety / appetite (FTO, MC4R, LEPR) ────────────────────────────────────

def _analyze_satiety(snps_df) -> dict:
    fto = _gt(snps_df, "rs9939609")
    mc4r = _gt(snps_df, "rs17782313")
    lepr = _gt(snps_df, "rs1137101")
    factors: list[str] = []
    appetite = "Standard"
    score = 0
    if fto:
        factors.append(f"rs9939609 (FTO) {fto}")
        if "A" in fto:
            score += fto.count("A")
    if mc4r:
        factors.append(f"rs17782313 (MC4R) {mc4r}")
        if "C" in mc4r:
            score += mc4r.count("C")
    if lepr:
        factors.append(f"rs1137101 (LEPR) {lepr}")
        if "G" in lepr:
            score += 1
    if score >= 3:
        appetite = "Elevated hunger drive"
    elif score >= 1:
        appetite = "Mildly elevated hunger"
    strategies = (
        "High-satiety eating pattern is critical. Anchor every meal with: (1) 30–40 g "
        "protein (eggs, Greek yoghurt, fish, lean meat, tofu), (2) 8–10 g fiber "
        "(legumes, berries, vegetables), (3) volume from non-starchy veg. Front-load "
        "calories to breakfast/lunch; keep dinner lighter. Avoid liquid calories "
        "(juice, smoothies as meals — they bypass satiety). Pre-meal water 500 mL + "
        "salad reduces total intake. 12-hour overnight fast helps; longer (16:8) only "
        "if it doesn't trigger evening binges."
        if appetite != "Standard" else
        "Normal appetite regulation — hunger/fullness cues are reliable. Standard "
        "guidelines apply."
    )
    return {
        "appetite_phenotype": appetite,
        "satiety_score": score,
        "factors": factors or ["Appetite SNPs not typed"],
        "guidance": strategies,
    }


# ── Saturated-fat sub-typing (APOE detailed + APOA2) ────────────────────────

def _analyze_saturated_fat(snps_df) -> dict:
    apoe1 = _gt(snps_df, "rs429358")
    apoe2 = _gt(snps_df, "rs7412")
    apoa2 = _gt(snps_df, "rs5082")
    factors: list[str] = []
    apoe_geno = "Unknown"
    if apoe1 and apoe2:
        # ε2: rs429358 T + rs7412 T; ε3: T+C; ε4: C+C
        e4 = "C" in apoe1
        e2 = "T" in apoe2
        if e4 and not e2:
            apoe_geno = "ε4 carrier" if apoe1.count("C") == 1 else "ε4/ε4"
        elif e2 and not e4:
            apoe_geno = "ε2 carrier"
        else:
            apoe_geno = "ε3/ε3 (typical)"
        factors.append(f"APOE: rs429358={apoe1}, rs7412={apoe2} ({apoe_geno})")
    high_sat_risk = "C" in (apoe1 or "") and "T" not in (apoe2 or "")
    if apoa2:
        factors.append(f"rs5082 (APOA2) {apoa2}")
    sat_cap_g = 15 if high_sat_risk else 22  # ~7% vs ~10% of 2000 kcal
    if high_sat_risk:
        guidance = (
            f"APOE ε4 carrier — markedly elevated LDL response to saturated fat. Cap "
            f"saturated fat at ≤{sat_cap_g} g/day (~7% of calories). Replace butter/"
            f"coconut oil with extra-virgin olive oil and avocado; limit red meat to "
            f"1–2×/week, choose poultry/fish. Prioritize MUFA (olive, avocado, almonds) "
            f"and PUFA (walnuts, fatty fish). Mediterranean pattern is best-evidence."
        )
    else:
        guidance = (
            f"Standard saturated-fat handling — keep ≤{sat_cap_g} g/day (~10% of "
            f"calories). Full-fat dairy and moderate red meat acceptable; emphasize "
            f"MUFA/PUFA still."
        )
    return {
        "apoe_genotype": apoe_geno,
        "saturated_fat_cap_g": sat_cap_g,
        "factors": factors or ["APOE/APOA2 not typed"],
        "guidance": guidance,
    }


# ── Meal-timing / chrononutrition ───────────────────────────────────────────

def _analyze_meal_timing(snps_df) -> dict:
    clock = _gt(snps_df, "rs1801260")
    melatonin = _gt(snps_df, "rs10830963")  # MTNR1B G allele = impaired glucose tolerance evening
    factors: list[str] = []
    eating_window = "08:00–20:00"
    note = ""
    if clock:
        factors.append(f"rs1801260 (CLOCK) {clock}")
        if "T" in clock:  # evening chronotype
            eating_window = "10:00–20:00"
            note += "Evening chronotype — push first meal later, but cap last meal ≤20:00 to preserve overnight fast. "
    if melatonin:
        factors.append(f"rs10830963 (MTNR1B) {melatonin}")
        if "G" in melatonin:
            note += (
                "MTNR1B G-allele — glucose tolerance drops sharply in evening. Eat largest "
                "carb load at breakfast/lunch; keep dinner protein+veg-heavy with minimal "
                "starch after 19:00. "
            )
    if not note:
        note = "Standard 10–12 h eating window; finish dinner ≥2 h before sleep."
    return {
        "eating_window": eating_window,
        "factors": factors or ["Chrono-nutrition SNPs not typed"],
        "guidance": note.strip(),
    }


# ── Antioxidant capacity ────────────────────────────────────────────────────

def _analyze_antioxidants(snps_df) -> dict:
    sod2 = _gt(snps_df, "rs4880")
    gpx1 = _gt(snps_df, "rs1050450")
    nqo1 = _gt(snps_df, "rs1800566")
    factors: list[str] = []
    low = 0
    if sod2:
        factors.append(f"rs4880 (SOD2 Ala16Val) {sod2}")
        if "G" in sod2:
            low += sod2.count("G") * 0.5
    if gpx1:
        factors.append(f"rs1050450 (GPX1) {gpx1}")
        if "T" in gpx1:
            low += 1
    if nqo1:
        factors.append(f"rs1800566 (NQO1) {nqo1}")
        if "T" in nqo1:
            low += 1
    needs = low >= 1.5
    return {
        "reduced_capacity": needs,
        "factors": factors or ["Antioxidant SNPs not typed"],
        "guidance": (
            "Reduced endogenous antioxidant enzymes — supply abundant dietary "
            "antioxidants. Daily targets: 2 cups dark leafy greens, 1 cup berries, "
            "cruciferous veg 3–4×/week (sulforaphane induces Nrf2), green tea 2–3 "
            "cups, herbs/spices generously (turmeric, rosemary, oregano). Selenium: "
            "2 Brazil nuts/week. Glutathione precursors: whey, eggs, alliums. Avoid "
            "antioxidant megadose pills — paradoxically blunt training adaptations."
            if needs else
            "Standard endogenous antioxidant capacity — a normal varied diet suffices. "
            "Aim for diverse plant colors weekly."
        ),
    }


# ── Fiber target ────────────────────────────────────────────────────────────

def _analyze_fiber(snps_df, macros: dict, satiety: dict, salt: dict) -> dict:
    base = 38  # men default
    extra: list[str] = []
    if macros["pct_carbs"] >= 50:
        base = 45
        extra.append("Higher-carb prescription — push fiber to anchor glycemic response.")
    if satiety["appetite_phenotype"] != "Standard":
        base = max(base, 40)
        extra.append("Elevated appetite — fiber is the most powerful free satiety lever.")
    if salt["sensitive"]:
        extra.append("Salt-sensitive — potassium-rich high-fiber foods (legumes, potato) double-up benefit.")
    return {
        "target_g": base,
        "factors": extra or ["Standard fiber target"],
        "guidance": (
            f"Target {base} g fiber/day — most adults eat 15. Practical: 1 cup oats "
            "(8 g), 1 cup berries (8 g), 1 cup beans/lentils (15 g), 2 cups veg (10 g), "
            "1 oz chia/flax (10 g). Ramp up 5 g/week to avoid GI distress; pair with "
            "water increase."
        ),
    }


# ── Hydration ───────────────────────────────────────────────────────────────

def _analyze_hydration(snps_df, salt: dict) -> dict:
    avp = _gt(snps_df, "rs1042615")
    note = "Baseline: 30–35 mL/kg body weight daily, +500–750 mL per hour exercise."
    if salt["sensitive"]:
        note += " Salt-sensitive: pre-load 500 mL water on waking; emphasize potassium-rich fluids (coconut water, broth) over electrolyte powders."
    return {
        "target_ml_per_kg": 33,
        "guidance": note,
        "factors": ([f"rs1042615 (AVPR1A) {avp}"] if avp else ["Hydration SNPs not typed"]),
    }


# ── Caloric framework (no body data, give worksheet) ────────────────────────

def _caloric_framework(macros: dict, satiety: dict) -> dict:
    return {
        "tdee_formula": "Mifflin-St Jeor: BMR(♂) = 10·kg + 6.25·cm − 5·age + 5; (♀) −161. TDEE = BMR × activity (1.4 sedentary, 1.55 moderate, 1.75 active, 1.9 athlete).",
        "loss_deficit_kcal": 400 if satiety["appetite_phenotype"] != "Standard" else 500,
        "gain_surplus_kcal": 250,
        "protein_g_per_kg": 1.8 if macros["pct_carbs"] <= 35 else 1.6,
        "guidance": (
            "Compute TDEE from formula. For fat loss: subtract 400–500 kcal (smaller "
            "deficit if appetite-elevated genotype — sustainability beats speed). For "
            "lean gain: +250 kcal with strength training. Protein floor regardless of "
            "goal: 1.6–2.0 g/kg body weight, distributed across 3–4 meals."
        ),
    }


# ── Public API ──────────────────────────────────────────────────────────────

def analyze_nutrition(snps_df: pd.DataFrame | None) -> dict:
    if snps_df is None:
        return {"status": "no_data"}

    macros = _analyze_macros(snps_df)
    caffeine = _analyze_caffeine(snps_df)
    alcohol = _analyze_alcohol(snps_df)
    salt = _analyze_salt(snps_df)
    lactose = _analyze_lactose(snps_df)
    gluten = _analyze_gluten(snps_df)
    methyl = _analyze_methylation_diet(snps_df)
    vit_d = _analyze_vitamin_d_food(snps_df)
    omega3 = _analyze_omega3(snps_df)
    iron = _analyze_iron(snps_df)
    choline = _analyze_choline(snps_df)
    b12 = _analyze_b12(snps_df)
    vit_a = _analyze_vitamin_a(snps_df)
    vit_c = _analyze_vitamin_c(snps_df)
    vit_e = _analyze_vitamin_e(snps_df)
    taste = _analyze_taste(snps_df)
    satiety = _analyze_satiety(snps_df)
    sat_fat = _analyze_saturated_fat(snps_df)
    meal_timing = _analyze_meal_timing(snps_df)
    antiox = _analyze_antioxidants(snps_df)
    fiber = _analyze_fiber(snps_df, macros, satiety, salt)
    hydration = _analyze_hydration(snps_df, salt)
    caloric = _caloric_framework(macros, satiety)
    foods = _build_food_lists(macros, alcohol, lactose, salt, gluten)
    # Augment food lists with new insights
    if iron["overload_risk"].startswith("High") or "Moderate" in iron["overload_risk"]:
        foods["avoid"].extend(["Iron-fortified cereals/breads", "Red meat >1×/week",
                               "Iron-containing multivitamins"])
    if omega3["ala_conversion"] == "Poor":
        foods["emphasize"].append("Oily fish 3–4×/week (salmon, sardines, mackerel) — FADS poor converter")
    if sat_fat["apoe_genotype"].startswith("ε4"):
        foods["avoid"].extend(["Butter/coconut oil as primary fats", "Processed red meat"])
        foods["emphasize"].append("Extra-virgin olive oil (primary fat), Mediterranean pattern")
    if choline["increased_need"]:
        foods["emphasize"].append("Whole eggs 1–2/day, liver 1 oz weekly (PEMT choline need)")
    if antiox["reduced_capacity"]:
        foods["emphasize"].append("Cruciferous vegetables 3–4×/week (sulforaphane → Nrf2)")
    if taste["bitter"] == "Super-taster":
        foods["emphasize"].append("Roasted/glazed cruciferous veg (super-taster mitigation)")
    template = _build_daily_template(macros, caffeine, alcohol)

    result = {
        "status": "ok",
        "macros": macros,
        "caffeine": caffeine,
        "alcohol": alcohol,
        "salt": salt,
        "lactose": lactose,
        "gluten": gluten,
        "methylation": methyl,
        "vitamin_d_food": vit_d,
        "omega3": omega3,
        "iron": iron,
        "choline": choline,
        "b12": b12,
        "vitamin_a": vit_a,
        "vitamin_c": vit_c,
        "vitamin_e": vit_e,
        "taste": taste,
        "satiety": satiety,
        "saturated_fat": sat_fat,
        "meal_timing": meal_timing,
        "antioxidants": antiox,
        "fiber": fiber,
        "hydration": hydration,
        "caloric": caloric,
        "emphasize": foods["emphasize"],
        "avoid": foods["avoid"],
        "daily_template": template,
    }

    if analyze_advanced_nutrition is not None:
        try:
            advanced = analyze_advanced_nutrition(snps_df, result)
            result.update(advanced)
        except Exception as exc:
            result["advanced_error"] = str(exc)

    if analyze_nutrition_protocols is not None:
        try:
            protocols = analyze_nutrition_protocols(result)
            result["protocols"] = protocols
        except Exception as exc:
            result["protocols_error"] = str(exc)

    return result


# ── HTML rendering ──────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_NU_CSS = """
<style>
.nu-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           color:#222; max-width: 1100px; margin: 24px auto; padding: 0 16px; }
.nu-wrap h1 { font-size: 1.6em; border-bottom: 2px solid #333; padding-bottom: 6px; }
.nu-wrap h2 { font-size: 1.2em; margin-top: 28px; padding-bottom:4px;
              border-bottom: 1px solid #eee; }
.nu-card { background:#fcfcfd; border:1px solid #e2e2e6; border-radius:10px;
           padding:14px 16px; margin:10px 0; }
.nu-macro-bar { display:flex; height:24px; border-radius:12px; overflow:hidden; margin: 8px 0; }
.nu-macro-bar .c { background:#5a8f3a; }
.nu-macro-bar .f { background:#c08327; }
.nu-macro-bar .p { background:#3a5a8f; }
.nu-macro-bar span { color:white; padding:0 10px; line-height:24px; font-size:0.9em; }
.nu-factors { font-family: Menlo, monospace; font-size:0.85em; color:#555;
              background:#f6f6f7; padding:6px 10px; border-radius:6px; margin-top:6px; }
.nu-two { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
@media (max-width:700px){ .nu-two { grid-template-columns: 1fr; } }
.nu-em { color:#2c7a30; }
.nu-av { color:#a32a2a; }
ul.nu-list { margin: 4px 0 0 18px; padding: 0; }
table.nu { width:100%; border-collapse: collapse; }
table.nu th, table.nu td { padding:8px 10px; border-bottom:1px solid #eee; text-align:left; }
table.nu th { background:#f9f9f9; }
.pgs-bar { background:#eee; border-radius:6px; height:14px; position:relative; overflow:hidden; }
.pgs-bar > div { height:100%; background:linear-gradient(90deg,#3a5a8f,#a32a2a); }
.pgs-row { display:grid; grid-template-columns: 200px 80px 1fr; gap:8px; align-items:center;
           padding:4px 0; border-bottom:1px solid #f0f0f0; font-size:0.9em; }
.dash-axis { display:grid; grid-template-columns: 240px 90px 1fr; gap:10px; align-items:start;
             padding:8px 0; border-bottom:1px solid #f0f0f0; }
.dash-tier { font-weight:600; }
.recipe-card { background:#fdfcf7; border:1px solid #e9e3c8; border-radius:10px;
               padding:14px 16px; margin:10px 0; }
.recipe-card h3 { margin: 0 0 4px; font-size: 1.05em; }
.matrix-table td:last-child { text-align:right; font-variant-numeric: tabular-nums; }
</style>
"""


def _render_advanced_sections(result: dict) -> str:
    pgs = result.get("polygenic_scores")
    dash = result.get("cardiometabolic_dashboard")
    inflam = result.get("inflammation")
    histamine = result.get("histamine")
    detox = result.get("detoxification")
    glycemic = result.get("glycemic_threshold")
    period = result.get("macro_periodisation")
    matrix = result.get("weekly_food_matrix")
    shop = result.get("shopping_list")
    recipes = result.get("recipes")

    # NOTE: every section below (dashboard, inflammation, histamine, detox,
    # glycemic, periodisation, matrix, shopping, recipes) and _render_protocols
    # guards on its OWN data. Do NOT gate the whole function on `pgs` — that
    # previously dropped all of them, plus the entire protocols layer (30-day
    # plan, glucose sim, minerals, cycle phase), whenever no polygenic scores
    # were computed, even though those sections were all computed and present.
    # Polygenic scores table (rendered only when polygenic scores exist).
    pgs_html = ""
    pgs_rows = []
    for trait, s in (pgs or {}).items():
        pct = s.get("percentile")
        bar_width = pct if pct is not None else 0
        pct_text = f"{pct}%" if pct is not None else "—"
        cov_text = f"<small style='color:#888'>coverage {s.get('coverage','—')}%</small>"
        pgs_rows.append(
            f'<div class="pgs-row">'
            f'<div>{_esc(trait.replace("_"," "))} <br>{cov_text}</div>'
            f'<div><strong>{pct_text}</strong><br><small>{_esc(s.get("tier",""))}</small></div>'
            f'<div><div class="pgs-bar"><div style="width:{bar_width}%"></div></div></div>'
            f'</div>'
        )
    if pgs:
        pgs_html = f"""
<h2>Polygenic Scores (percentile vs general population)</h2>
<div class="nu-card">
  <p style="font-size:0.88em;color:#666">Composite gene scores from GWAS-weighted SNPs.
  Higher percentile = stronger genetic predisposition (good or bad depending on trait).
  Confidence depends on coverage — most full-genome panels cover &gt;80%.</p>
  {"".join(pgs_rows)}
</div>
"""

    # Cardiometabolic dashboard
    dash_html = ""
    if dash:
        rows = "".join(
            f'<div class="dash-axis">'
            f'<div><strong>{_esc(a["axis"])}</strong></div>'
            f'<div class="dash-tier">{_esc(a.get("tier",""))}<br><small>{a.get("percentile","—")}%</small></div>'
            f'<div>{_esc(a.get("leverage",""))}</div>'
            f'</div>'
            for a in dash["axes"]
        )
        dash_html = f"""
<h2>Cardiometabolic Dashboard</h2>
<div class="nu-card">{rows}</div>
"""

    # Inflammation
    inflam_html = ""
    if inflam:
        f = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in inflam.get("factors", []))
        inflam_html = f"""
<h2>Inflammation Index</h2>
<div class="nu-card">
  <div><strong>{_esc(inflam.get("tier",""))}</strong> (score {inflam.get("score",0)})</div>
  {f}
  <p>{_esc(inflam.get("guidance",""))}</p>
</div>
"""

    # Histamine
    hist_html = ""
    if histamine:
        f = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in histamine.get("factors", []))
        hist_html = f"""
<h2>Histamine Tolerance</h2>
<div class="nu-card">
  <div><strong>{"Elevated risk" if histamine.get("elevated_risk") else "Normal clearance"}</strong></div>
  {f}<p>{_esc(histamine.get("guidance",""))}</p>
</div>
"""

    # Detox
    detox_html = ""
    if detox:
        p1 = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in detox.get("phase1_typed", []))
        p2 = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in detox.get("phase2_typed", []))
        detox_html = f"""
<h2>Detoxification (Phase I / II)</h2>
<div class="nu-card">
  <strong>Phase I (CYP enzymes)</strong>{p1}
  <strong>Phase II (conjugation)</strong>{p2}
  <p>Cruciferous target: <strong>{detox.get("cruciferous_target_servings_per_week","—")} servings/week</strong></p>
  <p>{_esc(detox.get("guidance",""))}</p>
</div>
"""

    # Glycemic threshold
    gly_html = ""
    if glycemic:
        f = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in glycemic.get("factors", []))
        gly_html = f"""
<h2>Personal Glycemic Threshold</h2>
<div class="nu-card">
  <div>Carb ceiling per meal: <strong>{glycemic.get("max_carbs_per_meal_g","—")} g</strong>
       (dinner: <strong>{glycemic.get("max_carbs_dinner_g","—")} g</strong>)</div>
  {f}<p>{_esc(glycemic.get("guidance",""))}</p>
</div>
"""

    # Macro periodisation
    period_html = ""
    if period:
        td = period["training_day"]
        rd = period["rest_day"]
        period_html = f"""
<h2>Macro Periodisation (Training vs Rest Day)</h2>
<div class="nu-card">
  <table class="nu">
    <tr><th>Day type</th><th>Carbs</th><th>Fat</th><th>Protein</th><th>Carb timing</th></tr>
    <tr><td>Training</td><td>{td['pct_carbs']}%</td><td>{td['pct_fat']}%</td>
        <td>{td['pct_protein']}%</td><td>{_esc(td['carb_timing'])}</td></tr>
    <tr><td>Rest</td><td>{rd['pct_carbs']}%</td><td>{rd['pct_fat']}%</td>
        <td>{rd['pct_protein']}%</td><td>{_esc(rd['carb_timing'])}</td></tr>
  </table>
  <p>{_esc(period.get("guidance",""))}</p>
</div>
"""

    # Weekly food matrix
    matrix_html = ""
    if matrix:
        sv = matrix["servings_per_week"]
        rows = "".join(f"<tr><td>{_esc(k)}</td><td>{v}</td></tr>" for k, v in sv.items())
        matrix_html = f"""
<h2>Quantitative Weekly Food Matrix (servings per week)</h2>
<div class="nu-card">
  <table class="nu matrix-table"><tr><th>Food group</th><th>Servings/wk</th></tr>{rows}</table>
</div>
"""

    # Shopping list
    shop_html = ""
    if shop:
        groups = "".join(
            f"<div class='nu-card'><strong>{_esc(g['category'])}</strong><ul class='nu-list'>"
            + "".join(f"<li>{_esc(i)}</li>" for i in g['items'] if i and i != "—")
            + "</ul></div>"
            for g in shop
        )
        shop_html = f"<h2>Weekly Shopping List (1 person)</h2>{groups}"

    # Recipes
    recipes_html = ""
    if recipes:
        cards = ""
        for r in recipes:
            ing = "".join(f"<li>{_esc(x)}</li>" for x in r["ingredients"])
            why = "".join(f"<li>{_esc(x)}</li>" for x in r["why_for_you"])
            cards += f"""
<div class="recipe-card">
  <h3>{_esc(r["name"])}</h3>
  <div style="color:#666;font-size:0.88em">{_esc(r["macros_est"])}</div>
  <strong>Ingredients</strong><ul>{ing}</ul>
  <strong>Why this matches your genotype</strong><ul>{why}</ul>
  <strong>Method</strong><p>{_esc(r["method"])}</p>
</div>"""
        recipes_html = f"<h2>Generated Recipes Tailored to Your Genotype</h2>{cards}"

    advanced_blob = (pgs_html + dash_html + inflam_html + hist_html + detox_html
                     + gly_html + period_html + matrix_html + shop_html + recipes_html)
    advanced_blob += _render_protocols(result)
    return advanced_blob


def _render_protocols(result: dict) -> str:
    p = result.get("protocols")
    if not p:
        return ""
    out = ""

    # Glucose simulator table
    gs = p.get("glucose_simulator")
    if gs:
        rows = "".join(
            f"<tr><td>{_esc(m['meal'])}</td><td>{m['available_carbs_g']} g</td>"
            f"<td>{m['estimated_peak_mg_dL']}</td><td>{m['estimated_iAUC']}</td>"
            f"<td>{_esc(m['verdict'])}</td></tr>" for m in gs
        )
        out += f"""
<h2>Postprandial Glucose Simulator</h2>
<div class="nu-card">
<p style="color:#666;font-size:0.88em">Estimated glucose response to common meals at your personal threshold.
Peak ≥140 mg/dL = caution; ≥180 = restructure.</p>
<table class="nu"><tr><th>Meal</th><th>Avail. carbs</th><th>Est. peak</th><th>iAUC</th><th>Verdict</th></tr>{rows}</table>
</div>
"""

    # 30-day meal plan
    mp = p.get("meal_plan_30d")
    if mp:
        rows = "".join(
            f"<tr><td>{d['day']}</td><td>{_esc(d['breakfast'])}</td>"
            f"<td>{_esc(d['snack_am'])}</td><td>{_esc(d['lunch'])}</td>"
            f"<td>{_esc(d['snack_pm'])}</td><td>{_esc(d['dinner'])}</td></tr>"
            for d in mp
        )
        out += f"""
<h2>30-Day Meal Plan</h2>
<div class="nu-card">
<table class="nu" style="font-size:0.85em">
<tr><th>Day</th><th>Breakfast</th><th>AM snack</th><th>Lunch</th><th>PM snack</th><th>Dinner</th></tr>
{rows}
</table>
</div>
"""

    # Cooking
    ck = p.get("cooking_methods")
    if ck:
        rule_rows = "".join(
            f"<tr><td>{_esc(r['method'])}</td><td><strong>{_esc(r['AGE_load'])}</strong></td>"
            f"<td>{_esc(r['use_for'])}</td></tr>" for r in ck["rules"]
        )
        oils = "".join(f"<li>{_esc(k)}: {v} °C</li>" for k, v in ck["oil_smoke_points_C"].items())
        practical = "".join(f"<li>{_esc(x)}</li>" for x in ck["practical"])
        out += f"""
<h2>Cooking-Method Optimiser (AGE minimisation)</h2>
<div class="nu-card">
<p>{_esc(ck['rationale'])}</p>
<table class="nu"><tr><th>Method</th><th>AGE load</th><th>Use for</th></tr>{rule_rows}</table>
<strong>Oil smoke points</strong><ul>{oils}</ul>
<strong>Practical</strong><ul>{practical}</ul>
</div>
"""

    # Restaurant guides
    rg = p.get("restaurant_guides")
    if rg:
        cards = ""
        for g in rg:
            orders = "".join(f"<li>{_esc(o)}</li>" for o in g["order"])
            cards += (f'<div class="nu-card"><strong>{_esc(g["cuisine"])}</strong>'
                      f'<ul>{orders}</ul><p><em>Your tweak: {_esc(g["your_tweak"])}</em></p></div>')
        out += f"<h2>Restaurant Ordering Guides</h2>{cards}"

    # Fasting
    fast = p.get("fasting")
    if fast:
        caut = "".join(f"<li>{_esc(c)}</li>" for c in fast["cautions"])
        out += f"""
<h2>Intermittent-Fasting Matchmaker</h2>
<div class="nu-card">
<p><strong>Recommended:</strong> {_esc(fast['recommended_protocol'])}</p>
<p>{_esc(fast['rationale'])}</p>
<strong>Cautions</strong><ul>{caut}</ul>
<p>{_esc(fast['break_fast_meal'])}</p>
</div>
"""

    # Polyphenols
    pp = p.get("polyphenols")
    if pp:
        food_rows = "".join(
            f"<tr><td>{_esc(f['food'])}</td><td>{f['mg']}</td>"
            f"<td>{_esc(f['key_polyphenols'])}</td></tr>" for f in pp["foods"]
        )
        out += f"""
<h2>Polyphenol Prescription</h2>
<div class="nu-card">
<p>Target: <strong>{pp['target_mg_per_day']} mg/day</strong></p>
<table class="nu"><tr><th>Food</th><th>mg per serving</th><th>Key polyphenols</th></tr>{food_rows}</table>
<p>{_esc(pp['guidance'])}</p>
</div>
"""

    # Minerals
    mn = p.get("minerals")
    if mn:
        rows = "".join(
            f"<tr><td>{_esc(m['mineral'])}</td>"
            f"<td>{m.get('rda_mg', m.get('rda_mcg','—'))} {'mg' if 'rda_mg' in m else 'µg'}</td>"
            f"<td>{_esc(m['sources'])}</td><td>{_esc(m['note'])}</td></tr>"
            for m in mn
        )
        out += f"""
<h2>Mineral Panel (quantitative)</h2>
<div class="nu-card">
<table class="nu"><tr><th>Mineral</th><th>Target</th><th>Sources (mg/serving)</th><th>Note</th></tr>{rows}</table>
</div>
"""

    # MIND
    md = p.get("mind_diet")
    if md and md.get("applicable"):
        inc = "".join(
            f"<tr><td>{_esc(i['food'])}</td><td>{_esc(i['servings'])}</td><td>{_esc(i['why'])}</td></tr>"
            for i in md["include_weekly"]
        )
        lim = "".join(
            f"<tr><td>{_esc(row['food'])}</td><td>{_esc(row['limit'])}</td></tr>"
            for row in md["limit_strictly"]
        )
        leverage = "".join(f"<li>{_esc(x)}</li>" for x in md["high_leverage_actions"])
        e4spec = "".join(f"<li>{_esc(x)}</li>" for x in md.get("additional_ε4_specific", []))
        out += f"""
<h2>MIND Diet Protocol (APOE ε4 cognitive nutrition)</h2>
<div class="nu-card">
<p>{_esc(md['header'])}</p>
<strong>Include weekly</strong>
<table class="nu"><tr><th>Food</th><th>Servings</th><th>Why</th></tr>{inc}</table>
<strong>Limit</strong>
<table class="nu"><tr><th>Food</th><th>Limit</th></tr>{lim}</table>
<strong>Highest-leverage daily actions</strong><ul>{leverage}</ul>
<strong>ε4-specific additions</strong><ul>{e4spec}</ul>
</div>
"""

    # Cycle phase
    cp = p.get("cycle_phase")
    if cp:
        out += f"""
<h2>Female Cycle-Phase Nutrition (reference)</h2>
<div class="nu-card">
<p><em>{_esc(cp.get('note',''))}</em></p>
<ul>
  <li><strong>Follicular (d1-14):</strong> {_esc(cp['follicular_d1_14']['macro_emphasis'])}; training: {_esc(cp['follicular_d1_14']['training'])}</li>
  <li><strong>Ovulation (d14-16):</strong> {_esc(cp['ovulation_d14_16']['macro_emphasis'])}; training: {_esc(cp['ovulation_d14_16']['training'])}</li>
  <li><strong>Luteal (d17-28):</strong> {_esc(cp['luteal_d17_28']['macro_emphasis'])}; training: {_esc(cp['luteal_d17_28']['training'])}</li>
  <li><strong>Menstrual (d1-5):</strong> {_esc(cp['menstrual_d1_5']['macro_emphasis'])}; training: {_esc(cp['menstrual_d1_5']['training'])}</li>
</ul>
<p><strong>RED-S:</strong> {_esc(cp['RED_S_screening'])}</p>
</div>
"""

    # Travel
    tj = p.get("travel_jetlag")
    if tj:
        steps = "".join(
            f"<li><strong>{_esc(s['step'])}:</strong> {_esc(s['action'])}</li>"
            for s in tj["steps"]
        )
        sup = "".join(f"<li>{_esc(x)}</li>" for x in tj["supplements_considered"])
        out += f"""
<h2>Travel / Jet-Lag Nutrition Protocol</h2>
<div class="nu-card">
<p>{_esc(tj['header'])}</p>
<ol>{steps}</ol>
<strong>Supplements considered</strong><ul>{sup}</ul>
</div>
"""

    # Pre-bed
    pb = p.get("pre_bed_mps")
    if pb:
        opts = "".join(f"<li>{_esc(x)}</li>" for x in pb["options"])
        avoid = "".join(f"<li>{_esc(x)}</li>" for x in pb["avoid_pre_bed"])
        out += f"""
<h2>Pre-Bed Muscle Protein Synthesis</h2>
<div class="nu-card">
<p>{_esc(pb['rationale'])}</p>
<strong>Options</strong><ul>{opts}</ul>
<strong>Avoid pre-bed</strong><ul>{avoid}</ul>
</div>
"""

    return out


def render_nutrition_html(result: dict, file_label: str = "") -> str:
    if not result or result.get("status") != "ok":
        return f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Nutrition</title>{_NU_CSS}</head>
<body><div class="nu-wrap"><h1>Personalized Nutrition Plan</h1>
<p>Insufficient genetic data for nutrition recommendations.</p></div></body></html>"""

    m = result["macros"]
    macro_bar = (
        f'<div class="nu-macro-bar">'
        f'<span class="c" style="width:{m["pct_carbs"]}%">{m["pct_carbs"]}% Carbs</span>'
        f'<span class="f" style="width:{m["pct_fat"]}%">{m["pct_fat"]}% Fat</span>'
        f'<span class="p" style="width:{m["pct_protein"]}%">{m["pct_protein"]}% Protein</span>'
        f'</div>'
    )
    macro_factors = "".join(f'<div class="nu-factors">{_esc(f)}</div>' for f in m["factors"])

    em_html = "".join(f"<li>{_esc(x)}</li>" for x in result["emphasize"])
    av_html = "".join(f"<li>{_esc(x)}</li>" for x in result["avoid"]) or "<li>—</li>"

    def section(title, d, key="guidance"):
        factors = "".join(f'<div class="nu-factors">{_esc(f)}</div>' for f in d.get("factors", []))
        return f"""
<div class="nu-card">
  <strong>{_esc(title)}</strong>
  {factors}
  <p>{_esc(d.get(key,''))}</p>
</div>"""

    cells = [
        section("☕ Caffeine", result["caffeine"]),
        section("🍷 Alcohol", result["alcohol"]),
        section("🧂 Salt sensitivity", result["salt"]),
        section("🥛 Lactose", result["lactose"]),
        section("🌾 Gluten / celiac risk", result["gluten"]),
        section("🧬 Methylation (folate)", result["methylation"]),
        section("☀ Vitamin D from food", result["vitamin_d_food"]),
        section("🐟 Omega-3 (FADS conversion)", result.get("omega3", {})),
        section("🩸 Iron / hemochromatosis risk", result.get("iron", {})),
        section("🥚 Choline (PEMT)", result.get("choline", {})),
        section("💊 Vitamin B12", result.get("b12", {})),
        section("🥕 Vitamin A (β-carotene conversion)", result.get("vitamin_a", {})),
        section("🍊 Vitamin C", result.get("vitamin_c", {})),
        section("🌰 Vitamin E", result.get("vitamin_e", {})),
        section("🥦 Antioxidant capacity", result.get("antioxidants", {})),
    ]

    # Build extended cards
    sat = result.get("saturated_fat", {})
    sat_card = ""
    if sat:
        sf_factors = "".join(f'<div class="nu-factors">{_esc(f)}</div>' for f in sat.get("factors", []))
        sat_card = f"""
<div class="nu-card">
  <strong>🥩 Saturated-fat handling (APOE)</strong>
  <div>Genotype: <strong>{_esc(sat.get('apoe_genotype','—'))}</strong> ·
       cap: <strong>≤ {sat.get('saturated_fat_cap_g','—')} g/day</strong></div>
  {sf_factors}
  <p>{_esc(sat.get('guidance',''))}</p>
</div>"""

    tas = result.get("taste", {})
    taste_card = ""
    if tas:
        t_factors = "".join(f'<div class="nu-factors">{_esc(f)}</div>' for f in tas.get("factors", []))
        taste_card = f"""
<div class="nu-card">
  <strong>👅 Taste perception</strong>
  <div>Bitter: <strong>{_esc(tas.get('bitter','—'))}</strong> ·
       Fat-taste: <strong>{_esc(tas.get('fat_taste','—'))}</strong></div>
  {t_factors}
  <p>{_esc(tas.get('bitter_guidance',''))}</p>
  <p>{_esc(tas.get('fat_guidance',''))}</p>
</div>"""

    sat_y = result.get("satiety", {})
    satiety_card = ""
    if sat_y:
        f = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in sat_y.get("factors", []))
        satiety_card = f"""
<div class="nu-card">
  <strong>🍽 Appetite & satiety (FTO / MC4R / LEPR)</strong>
  <div>Phenotype: <strong>{_esc(sat_y.get('appetite_phenotype','—'))}</strong>
       (score {sat_y.get('satiety_score','—')})</div>
  {f}
  <p>{_esc(sat_y.get('guidance',''))}</p>
</div>"""

    mt = result.get("meal_timing", {})
    timing_card = ""
    if mt:
        f = "".join(f'<div class="nu-factors">{_esc(x)}</div>' for x in mt.get("factors", []))
        timing_card = f"""
<div class="nu-card">
  <strong>🕒 Meal timing / chrononutrition</strong>
  <div>Eating window: <strong>{_esc(mt.get('eating_window','—'))}</strong></div>
  {f}
  <p>{_esc(mt.get('guidance',''))}</p>
</div>"""

    fib = result.get("fiber", {})
    hyd = result.get("hydration", {})
    cal = result.get("caloric", {})
    targets_card = f"""
<div class="nu-card">
  <strong>📏 Daily targets</strong>
  <ul class="nu-list">
    <li>Fiber target: <strong>{fib.get('target_g','—')} g/day</strong> — {_esc(fib.get('guidance',''))}</li>
    <li>Hydration: ~{hyd.get('target_ml_per_kg','—')} mL/kg/day — {_esc(hyd.get('guidance',''))}</li>
    <li>Protein floor: <strong>{cal.get('protein_g_per_kg','—')} g/kg body weight</strong></li>
    <li>Saturated fat cap: <strong>≤ {sat.get('saturated_fat_cap_g','—')} g/day</strong></li>
  </ul>
  <p><em>Calorie framework:</em> {_esc(cal.get('tdee_formula',''))}</p>
  <p>{_esc(cal.get('guidance',''))}</p>
</div>"""

    template_rows = "".join(
        f"<tr><td>{_esc(d['meal'])}</td><td>{_esc(d['example'])}</td></tr>"
        for d in result["daily_template"]
    )

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Personalized Nutrition Plan{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_NU_CSS}</head><body><div class="nu-wrap">
<h1>Personalized Nutrition Plan</h1>

<h2>Macronutrient Ratio</h2>
<div class="nu-card">
  {macro_bar}
  <p>{_esc(m["rationale"])}</p>
  {macro_factors}
</div>

<h2>Foods</h2>
<div class="nu-two">
  <div class="nu-card">
    <div class="nu-em"><strong>Emphasize</strong></div>
    <ul class="nu-list">{em_html}</ul>
  </div>
  <div class="nu-card">
    <div class="nu-av"><strong>Avoid / limit</strong></div>
    <ul class="nu-list">{av_html}</ul>
  </div>
</div>

<h2>Daily Targets & Caloric Framework</h2>
{targets_card}

<h2>Saturated Fat & Cardiovascular Lipid Response</h2>
{sat_card}

<h2>Appetite, Satiety & Meal Timing</h2>
{satiety_card}
{timing_card}

<h2>Taste Perception & Adherence Strategy</h2>
{taste_card}

<h2>Stimulants, Sensitivities & Micronutrients</h2>
{"".join(cells)}

<h2>Example Daily Pattern</h2>
<div class="nu-card">
<table class="nu">
  <tr><th>Meal</th><th>Example</th></tr>
  {template_rows}
</table>
</div>

{_render_advanced_sections(result)}

<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. These are evidence-aligned starting points; refine with a
registered dietitian if you have diabetes, kidney disease, an eating disorder,
or are pregnant/breastfeeding. Caloric needs not specified here — adjust to
your goals (loss / maintenance / gain).
</p>
</div></body></html>"""
