"""
Personalised Nutrition Plan
===========================

Genotype-driven dietary prescription:

  • macro ratios (carb / fat / protein) — FTO, TCF7L2, PPARG, APOE, FADS
  • foods to emphasise / avoid
  • caffeine guidance — CYP1A2, ADORA2A
  • alcohol guidance — ALDH2, ADH1B
  • salt sensitivity — ACE, AGT
  • lactose tolerance — LCT (rs4988235)
  • gluten / coeliac risk — HLA-DQ2/8 tag SNPs
  • methylation diet — MTHFR
  • vitamin D from food — VDR / CYP2R1

Output dict:
  {
    status,
    macros: {pct_carbs, pct_fat, pct_protein, rationale},
    emphasise: [...], avoid: [...],
    caffeine: {...}, alcohol: {...}, salt: {...},
    lactose: {...}, gluten: {...},
    methylation: {...}, vitamin_d_food: {...},
    daily_template: [...],
  }
"""

from __future__ import annotations

from typing import Dict, List, Optional
import pandas as pd


def _gt(snps_df: Optional[pd.DataFrame], rsid: str) -> Optional[str]:
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

def _analyze_macros(snps_df) -> Dict:
    fto = _gt(snps_df, "rs9939609")       # A = obesity risk
    tcf7l2 = _gt(snps_df, "rs7903146")     # T = T2D risk + worse carb response
    ppara = _gt(snps_df, "rs1800206")      # PPARA — fat metabolism
    apoe1 = _gt(snps_df, "rs429358")       # APOE e4 risk SNP
    apoe2 = _gt(snps_df, "rs7412")
    fads = _gt(snps_df, "rs174547")

    carb_pressure = 0
    fat_pressure = 0
    factors: List[str] = []

    if fto:
        if "A" in fto:
            carb_pressure -= 1
            factors.append(f"rs9939609 (FTO) {fto} — better response to lower-carb diet")
        else:
            factors.append(f"rs9939609 (FTO) {fto} — favourable FTO; carb-tolerant")

    if tcf7l2:
        if "T" in tcf7l2:
            carb_pressure -= 2
            factors.append(f"rs7903146 (TCF7L2) {tcf7l2} — reduced glucose tolerance; minimise refined carbs")

    if apoe1 and apoe2:
        # rs429358 C + rs7412 C = ε4 carrier — saturated-fat sensitive
        if "C" in apoe1 and "C" in apoe2:
            fat_pressure -= 1
            factors.append(f"APOE ε4 carrier — lower saturated fat, prefer mono/poly")

    if fads and "T" in fads:
        factors.append(f"rs174547 (FADS1) {fads} — favour direct EPA/DHA from fish, not ALA")

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

def _analyze_caffeine(snps_df) -> Dict:
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

def _analyze_alcohol(snps_df) -> Dict:
    aldh2 = _gt(snps_df, "rs671")          # A allele (East Asian variant) → flushing
    adh1b = _gt(snps_df, "rs1229984")      # A allele → fast ethanol→acetaldehyde
    factors: List[str] = []
    risk = "Standard"

    if aldh2 and "A" in aldh2:
        if aldh2.count("A") == 2:
            risk = "Avoid entirely"
            factors.append(f"rs671 (ALDH2*2/*2) — non-functional ALDH2; acetaldehyde toxic accumulation")
        else:
            risk = "Strongly limit"
            factors.append(f"rs671 (ALDH2*1/*2) — partial deficiency; flushing, cancer-risk elevated")

    if adh1b and "A" in adh1b:
        factors.append(f"rs1229984 (ADH1B*2) — fast acetaldehyde production; aldehyde load higher")
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

def _analyze_salt(snps_df) -> Dict:
    ace = _gt(snps_df, "rs4341") or _gt(snps_df, "rs4646994")
    agt = _gt(snps_df, "rs699")         # AGT M235T; T = salt-sensitive
    factors: List[str] = []
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
            "Reduce sodium to ≤2.3 g/day; emphasise potassium-rich foods (leafy greens, "
            "beans, potato). DASH-style pattern especially beneficial."
            if sensitive else
            "No genetic salt-sensitivity flagged — keep sodium within general "
            "guidelines (≤2.3 g/day adult)."
        ),
    }


# ── Lactose ─────────────────────────────────────────────────────────────────

def _analyze_lactose(snps_df) -> Dict:
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

def _analyze_gluten(snps_df) -> Dict:
    dq2 = _gt(snps_df, "rs2187668")        # HLA-DQ2.5 tag
    dq8 = _gt(snps_df, "rs7454108")        # HLA-DQ8 tag
    carrier = (dq2 and "T" in dq2) or (dq8 and "C" in dq8)
    factors: List[str] = []
    if dq2:
        factors.append(f"rs2187668 (HLA-DQ2 tag) {dq2}")
    if dq8:
        factors.append(f"rs7454108 (HLA-DQ8 tag) {dq8}")
    return {
        "celiac_risk_haplotype": bool(carrier),
        "factors": factors or ["DQ2/DQ8 tag SNPs not typed"],
        "guidance": (
            "DQ2/DQ8 carrier — ~3% lifetime risk of developing coeliac disease (vs <0.1% "
            "without). If you have GI symptoms, fatigue, or family history, ask your "
            "physician for serology (tTG-IgA) before starting gluten-free diet."
            if carrier else
            "No DQ2/DQ8 risk haplotype typed/detected — coeliac disease is essentially "
            "ruled out. Gluten avoidance unnecessary for autoimmune reasons."
        ),
    }


# ── Methylation diet ────────────────────────────────────────────────────────

def _analyze_methylation_diet(snps_df) -> Dict:
    mthfr_c677t = _gt(snps_df, "rs1801133")
    mthfr_a1298c = _gt(snps_df, "rs1801131")
    factors: List[str] = []
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
            "Emphasise folate-rich foods: leafy greens (spinach, kale, romaine), liver, "
            "lentils, asparagus, broccoli. Avoid synthetic folic-acid fortified products "
            "(many cereals, breads) — they compete with active folate. Choose "
            "methylfolate-supplemented options if available."
            if needs_extra else
            "Standard varied diet supplies sufficient folate; aim for 1-2 cups leafy "
            "greens daily."
        ),
    }


# ── Vitamin D from food ─────────────────────────────────────────────────────

def _analyze_vitamin_d_food(snps_df) -> Dict:
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
            "Prioritise high-D foods: wild salmon, sardines, herring, egg yolks, "
            "UV-treated mushrooms, fortified dairy. Combine with daytime sun exposure "
            "where possible (10-30 min, depending on skin tone, latitude, season)."
            if needs_more else
            "Modest D intake from fatty fish 2×/week + occasional sun is sufficient."
        ),
    }


# ── Foods to emphasise / avoid (synthesised) ───────────────────────────────

def _build_food_lists(macros: Dict, alcohol: Dict, lactose: Dict, salt: Dict, gluten: Dict) -> Dict:
    emphasise: List[str] = []
    avoid: List[str] = []

    if macros["pct_carbs"] <= 35:
        avoid.extend(["Sugary drinks", "White bread / pastries", "Sweetened breakfast cereals"])
        emphasise.extend(["Non-starchy vegetables", "Berries", "Quinoa / oats (modest portions)"])
    else:
        emphasise.extend(["Whole grains (oats, brown rice, barley)", "Legumes", "Sweet potato"])

    emphasise.extend(["Fatty fish (salmon, sardines) 2-3×/week",
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
        emphasise.append("Lactose-free dairy, aged cheeses, kefir/yoghurt")
    if gluten["celiac_risk_haplotype"]:
        avoid.append("Standard gluten products if symptomatic (DQ2/DQ8 carrier)")

    return {"emphasise": emphasise, "avoid": avoid}


# ── Daily template ──────────────────────────────────────────────────────────

def _build_daily_template(macros: Dict, caffeine: Dict, alcohol: Dict) -> List[Dict]:
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


# ── Public API ──────────────────────────────────────────────────────────────

def analyze_nutrition(snps_df: Optional[pd.DataFrame]) -> Dict:
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
    foods = _build_food_lists(macros, alcohol, lactose, salt, gluten)
    template = _build_daily_template(macros, caffeine, alcohol)

    return {
        "status": "ok",
        "macros": macros,
        "caffeine": caffeine,
        "alcohol": alcohol,
        "salt": salt,
        "lactose": lactose,
        "gluten": gluten,
        "methylation": methyl,
        "vitamin_d_food": vit_d,
        "emphasise": foods["emphasise"],
        "avoid": foods["avoid"],
        "daily_template": template,
    }


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
</style>
"""


def render_nutrition_html(result: Dict, file_label: str = "") -> str:
    if not result or result.get("status") != "ok":
        return f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Nutrition</title>{_NU_CSS}</head>
<body><div class="nu-wrap"><h1>Personalised Nutrition Plan</h1>
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

    em_html = "".join(f"<li>{_esc(x)}</li>" for x in result["emphasise"])
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
        section("🌾 Gluten / coeliac risk", result["gluten"]),
        section("🧬 Methylation", result["methylation"]),
        section("☀ Vitamin D from food", result["vitamin_d_food"]),
    ]

    template_rows = "".join(
        f"<tr><td>{_esc(d['meal'])}</td><td>{_esc(d['example'])}</td></tr>"
        for d in result["daily_template"]
    )

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Personalised Nutrition Plan{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_NU_CSS}</head><body><div class="nu-wrap">
<h1>Personalised Nutrition Plan</h1>

<h2>Macronutrient Ratio</h2>
<div class="nu-card">
  {macro_bar}
  <p>{_esc(m["rationale"])}</p>
  {macro_factors}
</div>

<h2>Foods</h2>
<div class="nu-two">
  <div class="nu-card">
    <div class="nu-em"><strong>Emphasise</strong></div>
    <ul class="nu-list">{em_html}</ul>
  </div>
  <div class="nu-card">
    <div class="nu-av"><strong>Avoid / limit</strong></div>
    <ul class="nu-list">{av_html}</ul>
  </div>
</div>

<h2>Stimulants, Sensitivities & Special Considerations</h2>
{"".join(cells)}

<h2>Example Daily Pattern</h2>
<div class="nu-card">
<table class="nu">
  <tr><th>Meal</th><th>Example</th></tr>
  {template_rows}
</table>
</div>

<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. These are evidence-aligned starting points; refine with a
registered dietitian if you have diabetes, kidney disease, an eating disorder,
or are pregnant/breastfeeding. Caloric needs not specified here — adjust to
your goals (loss / maintenance / gain).
</p>
</div></body></html>"""
