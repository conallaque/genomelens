"""
Nutrition Protocols & Operational Tools
=======================================

Concrete, deployable layers on top of the analytic core:

  • postprandial glucose simulator (iAUC estimation for arbitrary meals)
  • 30-day meal plan generator
  • cooking-method optimiser (AGE minimisation)
  • restaurant cuisine ordering guides
  • intermittent-fasting protocol matchmaker
  • polyphenol & mineral quantitative panels
  • MIND-diet protocol (APOE ε4 detailed cognitive nutrition)
  • female cycle-phase nutrition (when reproductive state supplied)
  • travel / jet-lag nutrition protocol
  • pre-bed muscle-protein-synthesis prescription
"""

from __future__ import annotations

# ── Postprandial Glucose Simulator ─────────────────────────────────────────
#
# Predicts incremental area-under-curve (iAUC) above fasting for a meal based
# on its macro composition, the user's personal glycaemic threshold, and the
# presence of glycaemic-attenuating factors (vinegar, protein-first sequencing,
# walking after meal). Calibration: ~5 mg/dL·min per gram of available carbs
# at the population average, scaled inversely with personal threshold.

def simulate_postprandial_glucose(
    meal: dict, glycemic: dict, time_of_day: str = "lunch"
) -> dict:
    """
    meal: {
       'name': str, 'carbs_g': float, 'protein_g': float, 'fat_g': float,
       'fibre_g': float, 'vinegar': bool, 'protein_first': bool,
       'walking_after_min': int, 'gi': int (0-100)
    }
    """
    base_threshold = glycemic.get("max_carbs_per_meal_g", 75)
    ceiling = glycemic.get("max_carbs_per_meal_g", 75) if time_of_day != "dinner" \
        else glycemic.get("max_carbs_dinner_g", base_threshold - 15)

    carbs = max(0, meal.get("carbs_g", 0))
    fibre = max(0, meal.get("fibre_g", 0))
    protein = max(0, meal.get("protein_g", 0))
    fat = max(0, meal.get("fat_g", 0))
    gi = meal.get("gi", 55)  # default mixed meal

    available_carbs = max(0, carbs - 0.5 * fibre)
    gi_factor = gi / 55.0

    # Base iAUC (mg/dL·min) — simplified
    iauc = 5 * available_carbs * gi_factor
    # Personal threshold attenuation
    iauc *= (75 / base_threshold) ** 1.2

    # Modifiers
    if meal.get("protein_first") and protein >= 20:
        iauc *= 0.75            # protein-first sequencing -25%
    if meal.get("vinegar"):
        iauc *= 0.80
    if fibre >= 10:
        iauc *= 0.85
    if fat >= 15:
        iauc *= 0.92            # slows absorption
    walk = meal.get("walking_after_min", 0)
    if walk >= 10:
        iauc *= max(0.55, 1 - 0.04 * walk)  # up to ~−40% with 10 min walk
    if time_of_day == "dinner":
        iauc *= 1.18            # circadian glucose intolerance

    peak = 80 + 0.8 * available_carbs * gi_factor * (75 / base_threshold)
    if meal.get("protein_first"):
        peak *= 0.85
    if walk >= 10:
        peak *= 0.88
    if time_of_day == "dinner":
        peak *= 1.08

    # Verdict
    if peak < 140:
        verdict = "Safe — peak <140 mg/dL"
    elif peak < 180:
        verdict = "Caution — modify (split carbs, walk, protein-first)"
    else:
        verdict = "Excessive — restructure meal"

    return {
        "meal": meal.get("name", "(unnamed)"),
        "available_carbs_g": round(available_carbs, 1),
        "estimated_peak_mg_dL": round(peak),
        "estimated_iAUC": round(iauc),
        "verdict": verdict,
        "modifiers_applied": [
            ("protein-first sequencing" if meal.get("protein_first") else None),
            ("vinegar 1 tbsp pre-meal" if meal.get("vinegar") else None),
            (f"fibre {fibre} g" if fibre >= 10 else None),
            (f"fat {fat} g slows absorption" if fat >= 15 else None),
            (f"{walk} min post-meal walk" if walk >= 10 else None),
            ("evening circadian penalty" if time_of_day == "dinner" else None),
        ],
        "carb_ceiling_for_this_meal_g": ceiling,
    }


def simulate_common_meals(glycemic: dict) -> list[dict]:
    meals = [
        # Naïve carb-heavy
        {"name": "Bagel + orange juice (naïve)", "carbs_g": 75, "protein_g": 10,
         "fat_g": 3, "fibre_g": 3, "gi": 75},
        # Same carbs but engineered
        {"name": "Bagel + OJ — engineered (eggs first, vinegar, walk)",
         "carbs_g": 75, "protein_g": 30, "fat_g": 12, "fibre_g": 5, "gi": 75,
         "protein_first": True, "vinegar": True, "walking_after_min": 15},
        # White rice + chicken
        {"name": "White rice + grilled chicken + broccoli",
         "carbs_g": 55, "protein_g": 35, "fat_g": 10, "fibre_g": 6, "gi": 70,
         "protein_first": True},
        # Steel-cut oats
        {"name": "Steel-cut oats + berries + Greek yoghurt",
         "carbs_g": 45, "protein_g": 25, "fat_g": 10, "fibre_g": 9, "gi": 42},
        # Power salad
        {"name": "Power salad (greens + chickpeas + chicken + olive oil)",
         "carbs_g": 30, "protein_g": 35, "fat_g": 22, "fibre_g": 14, "gi": 35,
         "vinegar": True},
        # Pasta dinner
        {"name": "Spaghetti bolognese (dinner)",
         "carbs_g": 85, "protein_g": 30, "fat_g": 18, "fibre_g": 6, "gi": 55},
    ]
    return [
        simulate_postprandial_glucose(m, glycemic,
                                      "dinner" if "dinner" in m["name"].lower() else "lunch")
        for m in meals
    ]


# ── 30-day Meal Plan ───────────────────────────────────────────────────────
#
# Generates 30 unique day plans across 4-week rotation. Three breakfasts × 4
# lunches × 4 dinners × variations = ample variety with controlled constraints.

_BREAKFASTS_HC = [
    "Steel-cut oats + blueberries + walnuts + Greek yoghurt + cinnamon",
    "Sweet-potato hash + 2 eggs + spinach + avocado",
    "Overnight oats + chia + raspberries + almond butter",
    "Whole-grain toast + avocado + 2 eggs + smoked salmon + tomato",
    "Buckwheat pancakes + Greek yoghurt + mixed berries",
]
_BREAKFASTS_LC = [
    "3-egg veggie omelette + ½ avocado + side spinach",
    "Smoked-salmon + cream-cheese roll-ups + cucumber + capers",
    "Greek yoghurt + walnuts + raspberries + chia seeds",
    "Cottage cheese bowl + flax + walnuts + cinnamon",
    "Tofu scramble + mushrooms + spinach + nutritional yeast",
]
_LUNCHES = [
    "Mediterranean bowl: quinoa, chickpeas, cucumber, tomato, feta, olive oil, lemon",
    "Power salad: arugula, grilled chicken, white beans, avocado, pumpkin seeds",
    "Tuna Niçoise: greens, hard-boiled eggs, olives, green beans, potato, anchovy",
    "Lentil soup + side salad + olive oil + sourdough (or GF crackers)",
    "Salmon poke bowl: brown rice, edamame, cucumber, avocado, sesame, nori",
    "Turkey + hummus wrap (or lettuce wraps) + carrot sticks + bell pepper",
    "Cold soba + grilled tofu + sesame + scallions + miso dressing",
    "Stuffed bell peppers: ground turkey, quinoa, tomato, herbs",
]
_DINNERS = [
    "Sheet-pan salmon + broccoli + sweet potato + olive oil + lemon",
    "Slow-cooker chili: beans, ground turkey, tomato, peppers, cumin",
    "Grilled chicken thighs + roasted Brussels sprouts + farro",
    "Cod en papillote + cherry tomatoes + olives + capers + spinach",
    "Tofu stir-fry + bok choy + shiitake + brown rice + ginger-tamari",
    "Lamb kebabs + tzatziki (or coconut yoghurt) + cucumber salad + bulgur",
    "Shrimp scampi + zucchini noodles + garlic + parsley",
    "Bean & vegetable curry + cauliflower rice + lime + cilantro",
]
_SNACKS = [
    "Apple + 1 oz almonds",
    "Greek yoghurt + berries",
    "Carrots + hummus",
    "Hard-boiled egg + 1 oz cheese",
    "Edamame + sea salt",
    "Cottage cheese + cinnamon + walnuts",
    "Dark chocolate (≥75%) + few walnuts",
    "Beef/turkey jerky + olives",
]


def thirty_day_meal_plan(macros: dict, gluten: bool, lactose: bool,
                         iron_overload: bool, e4: bool,
                         high_omega3: bool) -> list[dict]:
    low_carb = macros["pct_carbs"] <= 35
    breakfasts = _BREAKFASTS_LC if low_carb else _BREAKFASTS_HC
    lunches = list(_LUNCHES)
    dinners = list(_DINNERS)

    def adapt(meal: str) -> str:
        m = meal
        if gluten:
            m = m.replace("Whole-grain toast", "Gluten-free toast")
            m = m.replace("sourdough", "gluten-free crackers")
            m = m.replace("Buckwheat pancakes", "Buckwheat pancakes (naturally GF)")
            m = m.replace("farro", "quinoa")
            m = m.replace("bulgur", "quinoa")
            m = m.replace("wrap", "lettuce wrap")
            m = m.replace("soba", "rice noodles")
        if lactose:
            m = m.replace("Greek yoghurt", "lactose-free Greek yoghurt")
            m = m.replace("feta", "aged pecorino")
            m = m.replace("cream-cheese", "lactose-free cream cheese")
            m = m.replace("Cottage cheese", "Skyr (lactose-free)")
            m = m.replace("tzatziki", "coconut-yoghurt tzatziki")
            m = m.replace("cheese", "aged cheese")
        if iron_overload:
            m = m.replace("ground turkey", "white beans + chicken")
            m = m.replace("lamb", "chicken")
            m = m.replace("Beef", "Turkey")
        if e4:
            m = m.replace("cream-cheese", "ricotta-light")
        if high_omega3 and "chicken" in m.lower() and "salmon" not in m.lower():
            # bias one chicken meal toward salmon weekly — leave as is, balance via plan
            pass
        return m

    plan: list[dict] = []
    for day in range(1, 31):
        b = adapt(breakfasts[(day - 1) % len(breakfasts)])
        lunch = adapt(lunches[(day - 1) % len(lunches)])
        d = adapt(dinners[(day - 1) % len(dinners)])
        s1 = _SNACKS[(day - 1) % len(_SNACKS)]
        s2 = _SNACKS[(day + 3) % len(_SNACKS)]
        # Ensure salmon ≥ 3×/week (for FADS poor converters)
        if high_omega3 and day % 2 == 0 and "salmon" not in d.lower() and "cod" not in d.lower() and "shrimp" not in d.lower():
            d = "Sheet-pan salmon + " + d.split(":", 1)[-1] if ":" in d else "Sheet-pan salmon + roasted vegetables"
        plan.append({
            "day": day,
            "breakfast": b,
            "snack_am": s1,
            "lunch": lunch,
            "snack_pm": s2,
            "dinner": d,
        })
    return plan


# ── Cooking-method Optimizer (AGE / oxidation) ─────────────────────────────

def cooking_method_optimizer() -> dict:
    return {
        "rationale": (
            "Advanced glycation end-products (AGEs) form when proteins/fats are "
            "exposed to dry, high heat. They accelerate vascular ageing and oxidative "
            "stress. Same food, different method, can differ 10-100× in AGE load."
        ),
        "rules": [
            {"method": "Steam, boil, poach, sous-vide", "AGE_load": "Lowest",
             "use_for": "Fish, vegetables, eggs (poached), chicken, beans"},
            {"method": "Stew, braise, slow-cook (wet, <100 °C)", "AGE_load": "Low",
             "use_for": "Tough cuts (chuck, shank), legumes, root vegetables"},
            {"method": "Stir-fry (high oil, short time)", "AGE_load": "Moderate",
             "use_for": "Vegetables + tofu; use a smoke-point-stable oil (avocado, refined olive)"},
            {"method": "Bake/roast at ≤180 °C with marinade", "AGE_load": "Moderate",
             "use_for": "Marinade (acid + olive oil) cuts AGE formation 30-50%"},
            {"method": "Grill, broil, pan-sear, deep-fry", "AGE_load": "High",
             "use_for": "Use sparingly; never char meat. Pair with antioxidant veg/herbs."},
            {"method": "Char-grill / barbecue", "AGE_load": "Very high",
             "use_for": "Avoid. If used: trim charred bits, marinate in lemon/vinegar/garlic ≥30 min first."},
        ],
        "oil_smoke_points_C": {
            "Extra-virgin olive oil": 190,
            "Refined olive oil": 240,
            "Avocado oil": 271,
            "Ghee": 250,
            "Coconut oil": 175,
            "Butter": 150,
            "Industrial seed oils (corn, soybean)": "Avoid as primary",
        },
        "practical": [
            "Default to steam, sous-vide, slow-cook, and lower-temp roasting.",
            "Marinade meats 30+ min in acid (lemon, vinegar, yoghurt) before any dry heat.",
            "Pair high-AGE meals with high-polyphenol foods (berries, herbs, tea).",
            "Use a thermometer — most home ovens overshoot by 15-20 °C.",
        ],
    }


# ── Restaurant / Cuisine Ordering Guides ───────────────────────────────────

def restaurant_guides(result: dict) -> list[dict]:
    e4 = result.get("saturated_fat", {}).get("apoe_genotype", "").startswith("ε4")
    appetite = result.get("satiety", {}).get("appetite_phenotype", "Standard") != "Standard"
    lactose = result.get("lactose", {}).get("tolerance", "").startswith("Intolerant")
    gluten = result.get("gluten", {}).get("celiac_risk_haplotype", False)
    iron_overload = "High" in result.get("iron", {}).get("overload_risk", "") or \
                    "Moderate" in result.get("iron", {}).get("overload_risk", "")
    histamine = result.get("histamine", {}).get("elevated_risk", False)

    guides = [
        {"cuisine": "Italian",
         "order": [
             "Antipasto: grilled vegetables, olives, prosciutto (skip if iron overload)",
             "Primi: minestrone, ribollita, or simple tomato pasta (½ portion if appetite-elevated)",
             "Secondi: grilled fish (branzino, salmon), chicken piccata, osso buco",
             "Contorni: sautéed greens, roasted vegetables, white-bean salad",
             "Skip: heavy cream sauces, deep-fried calamari, charred meats",
         ],
         "your_tweak": ("ApoE ε4 — choose tomato/olive-oil sauces over cream; skip butter." if e4 else "Standard ordering OK."),
         },
        {"cuisine": "Japanese / sushi",
         "order": [
             "Edamame, miso soup, seaweed salad to start",
             "Sashimi (salmon, tuna, mackerel) or chirashi over white rice",
             "Brown-rice rolls if available; skip tempura batter",
             "Grilled fish + steamed vegetables main",
             "Skip: spicy mayo rolls, cream-cheese rolls",
         ],
         "your_tweak": (
             "Histamine intolerance — choose freshly prepared sushi (white fish), avoid tuna/mackerel/soy sauce/wasabi."
             if histamine else
             "Sushi/sashimi excellent for FADS poor-converters (direct EPA/DHA)."
         ),
         },
        {"cuisine": "Mexican",
         "order": [
             "Guacamole + jicama sticks instead of chips",
             "Fish tacos on corn tortillas (gluten-free)",
             "Bowls: black beans + chicken + salsa + avocado + greens (skip rice if low-carb)",
             "Side: salsa verde, pico de gallo",
             "Skip: deep-fried chimichangas, queso, sour-cream lakes",
         ],
         "your_tweak": ("Lactose: ask for no cheese/sour cream" if lactose else "Standard ordering OK.")
         + (" Gluten: corn tortillas only." if gluten else ""),
         },
        {"cuisine": "Indian",
         "order": [
             "Tandoori chicken or fish tikka (marinated, lower AGEs than fried)",
             "Daal (lentil) + chana masala (chickpea)",
             "Side: cucumber raita (skip if lactose), kachumber salad",
             "Naan ½ portion or skip; brown basmati if available",
             "Skip: korma/butter chicken (heavy cream), samosas, pakoras",
         ],
         "your_tweak": ("ApoE ε4 — avoid butter chicken/korma; pick tikka, daal, vegetable curries." if e4 else "Standard.")
         + (" Gluten: skip naan; rice instead." if gluten else ""),
         },
        {"cuisine": "American / Steakhouse",
         "order": [
             "Wedge or house salad (dressing on side)",
             "Filet mignon or grilled fish — ≤6 oz portion",
             "Sides: roasted Brussels, asparagus, steamed broccoli, sweet potato (no marshmallow)",
             "Skip: fried appetisers, creamed spinach, mac & cheese, dessert bread",
         ],
         "your_tweak": ("Iron overload — choose fish over red meat; limit beef to 1×/month." if iron_overload else "Choose leaner cuts, smaller portions.")
         + (" Appetite-elevated — order an appetiser salad first." if appetite else ""),
         },
        {"cuisine": "Mediterranean / Greek",
         "order": [
             "Greek salad, tabbouleh, hummus, baba ghanoush",
             "Grilled lamb or chicken souvlaki, grilled fish",
             "Side: rice pilaf (½ portion), grilled vegetables",
             "Skip: spanakopita (phyllo/cheese), gyro (high-AGE)",
         ],
         "your_tweak": "Generally well-suited — Mediterranean pattern is best-evidence for most genotypes.",
         },
        {"cuisine": "Thai / Vietnamese",
         "order": [
             "Tom yum, pho (broth) to start",
             "Curries with vegetables + protein, brown rice if available",
             "Grilled fish, larb, summer rolls (rice paper = GF)",
             "Skip: pad thai (sugar-heavy), deep-fried spring rolls",
         ],
         "your_tweak": ("Gluten: confirm fish sauce / soy substitution to tamari." if gluten else "Standard."),
         },
        {"cuisine": "Fast / casual (when stuck)",
         "order": [
             "Chipotle: bowl, double protein, fajita veg, black beans, fresh salsa, guac",
             "Sweetgreen / Cava: kale + greens base, protein + chickpeas + veg, olive-oil-based dressing",
             "McDonald's last resort: grilled chicken salad, apple slices, water",
             "Starbucks: egg-white bites + Americano (skip pastries)",
         ],
         "your_tweak": "Pre-portion mentally: half the rice, skip the bread, double the vegetables.",
         },
    ]
    return guides


# ── Fasting Protocol Matchmaker ────────────────────────────────────────────

def fasting_matchmaker(result: dict) -> dict:
    appetite = result.get("satiety", {}).get("appetite_phenotype", "Standard")
    chrono = (result.get("meal_timing", {}).get("eating_window", "")
              if result.get("meal_timing") else "")
    t2d_pct = result.get("polygenic_scores", {}).get("T2D", {}).get("percentile", 50)
    iron_overload = "High" in (result.get("iron", {}).get("overload_risk", "") or "")
    # Female state would be required to fully assess — keep general

    if appetite != "Standard":
        recommended = "12:12 (gentle) or 14:10"
        rationale = (
            "Elevated-appetite genotype — long fasts can trigger evening binges. "
            "Start with 12-h overnight fast; progress to 14:10 only if it doesn't drive "
            "compensatory overeating. AVOID OMAD/20:4 — these undermine adherence."
        )
    elif t2d_pct is not None and t2d_pct >= 70:
        recommended = "16:8 (TRF) — early eating window"
        rationale = (
            "Elevated T2D genotype — time-restricted eating (08:00-16:00 or 09:00-17:00) "
            "consistently shown to improve insulin sensitivity, blood pressure, and "
            "weight. Earlier eating window beats later (eTRF > lTRF)."
        )
    elif chrono and "10:00" in chrono:
        recommended = "16:8 (mid-day window)"
        rationale = "Evening chronotype — 10:00-18:00 or 11:00-19:00 fits circadian cortisol/insulin profile."
    else:
        recommended = "14:10 baseline; optional 16:8 1-2×/week"
        rationale = "Moderate TRF works without compromising training fuelling."

    cautions = [
        "Pregnant/breastfeeding → no fasting beyond ≤12 h overnight.",
        "History of eating disorder → no restrictive windows.",
        "On glucose-lowering meds → coordinate with physician.",
        "Heavy training day → no fasted workouts >75 min unless adapted.",
    ]
    if iron_overload:
        cautions.append("Long fasts paradoxically raise hepcidin variability — monitor ferritin if cycling fasts >18 h.")
    return {
        "recommended_protocol": recommended,
        "rationale": rationale,
        "cautions": cautions,
        "break_fast_meal": (
            "Break-fast meal: 30-40 g protein + 8-10 g fibre + healthy fat + low-GI "
            "carbs ≤30 g. e.g. 3-egg omelette + spinach + avocado + ½ cup berries."
        ),
    }


# ── Polyphenol Quantitative Prescription ───────────────────────────────────

def polyphenol_panel(antioxidants: dict, inflam: dict) -> dict:
    boost = antioxidants.get("reduced_capacity") or (inflam.get("score", 0) >= 2)
    target_mg_day = 1500 if boost else 800
    foods = [
        {"food": "Berries (½ cup)", "mg": 380, "key_polyphenols": "Anthocyanins, ellagitannins"},
        {"food": "Coffee (8 oz)", "mg": 350, "key_polyphenols": "Chlorogenic acid"},
        {"food": "Dark chocolate (1 oz, ≥75%)", "mg": 250, "key_polyphenols": "Catechins, procyanidins"},
        {"food": "Green tea (8 oz)", "mg": 200, "key_polyphenols": "EGCG"},
        {"food": "Red onion (½ medium)", "mg": 130, "key_polyphenols": "Quercetin"},
        {"food": "Extra-virgin olive oil (1 tbsp)", "mg": 90, "key_polyphenols": "Oleocanthal, hydroxytyrosol"},
        {"food": "Apples (1 medium, with skin)", "mg": 110, "key_polyphenols": "Quercetin, catechins"},
        {"food": "Cherries (½ cup)", "mg": 130, "key_polyphenols": "Anthocyanins"},
        {"food": "Walnuts (1 oz)", "mg": 90, "key_polyphenols": "Ellagic acid, urolithins (precursor)"},
        {"food": "Pomegranate (½ fruit)", "mg": 600, "key_polyphenols": "Punicalagins"},
        {"food": "Turmeric (1 tsp + pepper + fat)", "mg": 200, "key_polyphenols": "Curcumin"},
        {"food": "Red wine (5 oz — if no alcohol contraindication)", "mg": 200, "key_polyphenols": "Resveratrol, anthocyanins"},
    ]
    return {
        "target_mg_per_day": target_mg_day,
        "foods": foods,
        "guidance": (
            f"Polyphenol target: {target_mg_day} mg/day. Mixing diverse sources beats "
            "single-source mega-dosing. Daily template: 1 cup berries (380), 8 oz coffee "
            "(350), 2 tbsp olive oil (180), 1 cup green tea (200), 1 apple (110), 1 oz "
            "dark chocolate (250). Total ≈1500 mg. Polyphenols cluster: anthocyanins "
            "(berries) + flavanols (cocoa, tea) + hydroxytyrosols (olive) + lignans "
            "(flax, sesame) + curcuminoids (turmeric) = broad signalling."
            + (" Reduced antioxidant capacity flagged — double down on diversity." if boost else "")
        ),
    }


# ── Mineral Quantitative Panel ─────────────────────────────────────────────

def mineral_panel(result: dict) -> list[dict]:
    salt_sensitive = result.get("salt", {}).get("sensitive", False)
    appetite = result.get("satiety", {}).get("appetite_phenotype", "Standard") != "Standard"
    return [
        {"mineral": "Magnesium", "rda_mg": 420 if not appetite else 500,
         "sources": "Pumpkin seeds (1 oz=170), dark chocolate (1 oz=65), almonds (1 oz=80), "
                   "spinach (cup=150), avocado (1=60), black beans (cup=120)",
         "note": "Most adults eat <300 mg. Glycinate form best-absorbed if supplementing 200-400 mg/night."},
        {"mineral": "Potassium", "rda_mg": 4700 if salt_sensitive else 3400,
         "sources": "Potato (1 medium=900), spinach (cup cooked=840), banana (450), "
                   "salmon (4 oz=500), lentils (cup=730), avocado (1=975), yoghurt (1 cup=380)",
         "note": "Salt-sensitive: lean hard on potassium — lowers SBP 4-5 mmHg as much as a drug."},
        {"mineral": "Zinc", "rda_mg": 11,
         "sources": "Oysters (3 oz=74!), beef (3 oz=7), pumpkin seeds (1 oz=2.2), "
                   "chickpeas (cup=2.5), cashews (1 oz=1.6)",
         "note": "Vegetarians need +50% (phytate). Don't exceed 40 mg supplemented — depletes copper."},
        {"mineral": "Selenium", "rda_mcg": 55,
         "sources": "Brazil nuts (1 nut=68-90 — TWO nuts = full week's need!), tuna (3 oz=92), "
                   "sardines (3 oz=45), eggs (15)",
         "note": "2 Brazil nuts/week is the simplest selenium strategy. Toxic above 400 µg/day."},
        {"mineral": "Copper", "rda_mg": 0.9,
         "sources": "Liver (3 oz=12!!), oysters (3 oz=4.8), cashews (1 oz=0.6), sunflower seeds (1 oz=0.5)",
         "note": "Excess zinc (>40 mg/day chronic) depletes copper — watch this if supplementing zinc."},
        {"mineral": "Iodine", "rda_mcg": 150,
         "sources": "Seaweed (1 sheet nori=20), cod (3 oz=99), iodised salt (1 g=77), yoghurt (cup=75), eggs (24)",
         "note": "Vegans avoiding seaweed and iodised salt commonly deficient — track carefully."},
        {"mineral": "Chromium", "rda_mcg": 35,
         "sources": "Broccoli (½ cup=11), grape juice (cup=8), whole grains, beef",
         "note": "Glucose tolerance support — minor effect; food sources sufficient."},
        {"mineral": "Calcium", "rda_mg": 1000,
         "sources": "Yoghurt (cup=300), milk (cup=300), sardines w/bones (3 oz=325), "
                   "kale (cup cooked=180), tofu set with Ca (½ cup=250-400)",
         "note": "If avoiding dairy: fortified plant milks (250-450 mg/cup), canned fish with bones, leafy greens."},
        {"mineral": "Sodium", "rda_mg": 2300 if not salt_sensitive else 1500,
         "sources": "Processed foods dominate — bread, soup, sauce. Be aware not avoid.",
         "note": ("Salt-sensitive — strict ≤1500 mg/day" if salt_sensitive else "Within 2300 mg easily achievable.")},
    ]


# ── MIND Diet (APOE ε4 detailed cognitive protocol) ────────────────────────

def mind_diet_protocol(e4: bool) -> dict:
    if not e4:
        return {"applicable": False,
                "note": "MIND-diet protocol is highest-leverage for ApoE ε4 carriers. "
                        "You're not ε4 — the standard Mediterranean pattern is sufficient."}
    return {
        "applicable": True,
        "header": (
            "MIND diet (Mediterranean-DASH-Intervention for Neurodegenerative Delay). "
            "Strongest evidence base for slowing cognitive decline; effect size in MIND "
            "trial ≈ 53% lower Alzheimer's risk with high adherence vs low."
        ),
        "include_weekly": [
            {"food": "Leafy greens", "servings": "≥6/week",
             "why": "Folate, lutein, vitamin K1 — slows cognitive decline ~11 yrs equivalent"},
            {"food": "Other vegetables", "servings": "≥1/day", "why": "Carotenoids, vitamin C"},
            {"food": "Berries", "servings": "≥2/week",
             "why": "Anthocyanins cross blood-brain barrier; flavonoids improve memory"},
            {"food": "Nuts (especially walnuts)", "servings": "≥5/week",
             "why": "α-linolenic acid, vitamin E, polyphenols"},
            {"food": "Olive oil", "servings": "Primary cooking fat",
             "why": "Oleocanthal — natural ibuprofen-like anti-inflammatory; MUFA"},
            {"food": "Whole grains", "servings": "≥3/day", "why": "B vitamins, fibre"},
            {"food": "Fish (esp oily)", "servings": "≥1/week (aim 2-3)",
             "why": "DHA — neuronal membrane; ε4 carriers especially benefit"},
            {"food": "Beans/legumes", "servings": "≥3/week", "why": "Folate, fibre, B vitamins"},
            {"food": "Poultry", "servings": "≥2/week", "why": "Lean protein"},
            {"food": "Wine (red, if no contraindication)", "servings": "≤1 glass/day",
             "why": "Resveratrol — discontinue if any ALDH2 variant"},
        ],
        "limit_strictly": [
            {"food": "Red meat (including beef, pork)", "limit": "≤3 servings/week"},
            {"food": "Butter / margarine", "limit": "<1 tbsp/day — olive oil instead"},
            {"food": "Cheese", "limit": "≤1 serving/week"},
            {"food": "Pastries / sweets", "limit": "≤4 servings/week"},
            {"food": "Fried / fast food", "limit": "≤1 serving/week — ideally 0"},
        ],
        "high_leverage_actions": [
            "Daily 1 cup leafy greens — biggest single win",
            "Walnuts daily (1 oz) — only nut with meaningful ALA",
            "Salmon or sardines 2-3×/week — direct DHA",
            "Eliminate butter; cook in EVOO",
            "Berries 2× weekly minimum (frozen counts)",
        ],
        "additional_ε4_specific": [
            "Choline 550 mg/day (eggs, liver) — precursor for acetylcholine",
            "Sleep 7-9 h — glymphatic Aβ clearance happens in deep sleep",
            "Aerobic exercise 150+ min/week — single most-evidenced cognitive intervention",
            "Avoid coconut oil and saturated tropical fats despite hype — ε4 LDL response is harsh",
            "Consider DHA-rich algae oil 1 g/day if fish intake limited",
        ],
    }


# ── Female cycle-phase nutrition (placeholder; activates if user supplied) ──

def cycle_phase_nutrition() -> dict:
    return {
        "note": "Activate by passing reproductive_state. Tailors macros, iron, and "
                "training-fuelling to follicular/ovulatory/luteal/menstrual phases.",
        "follicular_d1_14": {
            "macro_emphasis": "Higher carb tolerance (insulin sensitivity peak)",
            "iron": "Replace menstrual losses — vitamin C with iron-rich meals",
            "training": "Best window for strength PRs and high-intensity work",
        },
        "ovulation_d14_16": {
            "macro_emphasis": "Antioxidants for ovulatory inflammation — berries, leafy greens",
            "training": "Peak power output — schedule key sessions",
        },
        "luteal_d17_28": {
            "macro_emphasis": "+5-10% calories (BMR rises 100-300 kcal), more complex carbs and magnesium",
            "cravings": "Magnesium 400 mg/day from food (pumpkin seeds, dark chocolate, leafy greens)",
            "training": "Endurance feels harder (higher core temp); shorten intervals",
        },
        "menstrual_d1_5": {
            "macro_emphasis": "Iron-rich foods + vitamin C cofactor",
            "training": "Reduce volume 20-30% if symptoms; gentle movement helps cramps",
        },
        "RED_S_screening": (
            "Relative Energy Deficiency in Sport — if cycle becomes irregular or absent "
            "while training hard, ENERGY AVAILABILITY is likely too low. Target >40 "
            "kcal/kg lean mass/day. Long-term bone density, immunity, performance all "
            "suffer. This is a red-flag medical issue, not a body-comp 'win'."
        ),
    }


# ── Travel / Jet-lag nutrition ─────────────────────────────────────────────

def travel_jetlag_protocol() -> dict:
    return {
        "header": "5-step protocol to compress jet-lag adjustment from 1 day/timezone to ~half.",
        "steps": [
            {
                "step": "Pre-flight (24-48 h before)",
                "action": (
                    "Begin shifting meal times 1 h/day toward destination. East travel: "
                    "earlier; west: later. Hydrate aggressively. Stop caffeine 8 h before "
                    "intended sleep window."
                ),
            },
            {
                "step": "On the plane",
                "action": (
                    "Set watch to destination time immediately. Eat ONLY at destination "
                    "meal times — even if not hungry, force a small protein-fat meal. "
                    "Skip the in-flight cocktail (alcohol worsens circadian disruption). "
                    "240-360 mL water per hour."
                ),
            },
            {
                "step": "Arrival breakfast",
                "action": (
                    "First destination-morning meal: heavy protein (30 g) + caffeine + "
                    "10-15 min sunlight exposure. Anchors cortisol awakening response to "
                    "new time zone."
                ),
            },
            {
                "step": "First evening",
                "action": (
                    "Carb-leaning dinner finished by 19:00 local. Tryptophan-rich (turkey, "
                    "salmon, pumpkin seeds) + complex carbs to spike serotonin → melatonin "
                    "conversion. Melatonin 0.3-0.5 mg 4 h before intended sleep — micro-dose, "
                    "not the 3-5 mg sold OTC (too much sedates without phase-shifting)."
                ),
            },
            {
                "step": "Days 2-3",
                "action": (
                    "Maintain destination meal timing strictly. Don't nap >20 min before "
                    "14:00 local. Daily morning sunlight 15 min. Caffeine cutoff 8 h before "
                    "intended bedtime."
                ),
            },
        ],
        "supplements_considered": [
            "Melatonin 0.3 mg micro-dose 4 h before bed at destination",
            "L-theanine 200 mg with caffeine for smoother wake-up",
            "Electrolyte mix on flight — dehydration amplifies jet-lag fatigue",
        ],
    }


# ── Pre-bed muscle protein synthesis ───────────────────────────────────────

def pre_bed_mps_protocol() -> dict:
    return {
        "rationale": (
            "Overnight is the longest fasting window. 30-40 g slow-digesting protein "
            "1-2 h before sleep elevates overnight MPS ~22% (Trommelen 2016) and "
            "doesn't impair sleep quality. Especially valuable for those training PM "
            "or with strength/hypertrophy goals."
        ),
        "options": [
            "1 cup Greek yoghurt (or skyr) + 1 tbsp almond butter + cinnamon (~30 g protein)",
            "Casein shake 30 g + water",
            "Cottage cheese ½ cup + walnuts (~28 g protein)",
            "If lactose-intolerant: pea-isolate or rice-isolate 40 g + 1 tbsp tahini",
        ],
        "avoid_pre_bed": [
            "High-fat very-late meals (>20 g fat within 2 h of sleep can fragment sleep)",
            "Alcohol (suppresses REM, blunts MPS)",
            "Large carb loads if glycaemic-vulnerable — slow-release fine, sugar-spike not",
        ],
    }


# ── Public synthesis API ───────────────────────────────────────────────────

def analyze_nutrition_protocols(result: dict) -> dict:
    glycemic = result.get("glycemic_threshold", {})
    macros = result.get("macros", {"pct_carbs": 45})
    antiox = result.get("antioxidants", {})
    inflam = result.get("inflammation", {})

    lactose = result.get("lactose", {}).get("tolerance", "").startswith("Intolerant")
    gluten = result.get("gluten", {}).get("celiac_risk_haplotype", False)
    iron_overload = "High" in (result.get("iron", {}).get("overload_risk") or "") or \
                    "Moderate" in (result.get("iron", {}).get("overload_risk") or "")
    e4 = result.get("saturated_fat", {}).get("apoe_genotype", "").startswith("ε4")
    high_omega3 = result.get("omega3", {}).get("ala_conversion") == "Poor"

    return {
        "glucose_simulator": simulate_common_meals(glycemic) if glycemic else None,
        "meal_plan_30d": thirty_day_meal_plan(macros, gluten, lactose, iron_overload, e4, high_omega3),
        "cooking_methods": cooking_method_optimizer(),
        "restaurant_guides": restaurant_guides(result),
        "fasting": fasting_matchmaker(result),
        "polyphenols": polyphenol_panel(antiox, inflam),
        "minerals": mineral_panel(result),
        "mind_diet": mind_diet_protocol(e4),
        "cycle_phase": cycle_phase_nutrition(),
        "travel_jetlag": travel_jetlag_protocol(),
        "pre_bed_mps": pre_bed_mps_protocol(),
    }
