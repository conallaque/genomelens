"""
Environmental Optimization — behavioural protocols the report should generate
=============================================================================

Cross-panel, *actionable* behavioural protocols derived from genotype — the
"what should I actually do differently" layer that most reports leave implicit.
Three domains:

  1. **Circadian light timing** — from chronotype variants (CLOCK), a concrete
     morning-light / evening-dark / melatonin-timing protocol tuned to whether
     the person leans morning or evening.
  2. **Exercise-modality fit** — from ACTN3 R577X and an ACE I/D proxy, a
     power-vs-endurance lean with training-emphasis recommendations. Framed as
     a *tendency that shifts what responds fastest*, not a ceiling — everyone
     should do both strength and cardio for health.
  3. **Vitamin-D seasonality** — vitamin-D-pathway variants (GC, CYP2R1,
     DHCR7, VDR) combined with **latitude** to produce a seasonal
     supplementation protocol. Latitude drives the recommendation (winter
     cutaneous synthesis effectively stops above ~37°N from ~Oct-Mar
     regardless of genotype); genotype modulates how strongly.

Honesty notes
-------------
* Latitude is an explicit parameter (default 40.0 °N, "temperate northern
  hemisphere"). The assumption is always stated in the output.
* Vitamin-D effect alleles are taken from the SUNLIGHT consortium
  (Wang 2010; Jiang 2018), but consumer-chip strand conventions can flip these
  calls — so the genotype contribution is reported as a modest *tendency tier*
  and the actionable recommendation leans on latitude, which is unambiguous.

References
----------
Katzenberg 1998; Roenneberg 2007 (chronotype & light). Yang 2009; MacArthur
2007 (ACTN3 R577X & muscle performance). Montgomery 1998; Puthucheary 2011
(ACE I/D & endurance). Wang 2010 & Jiang 2018 SUNLIGHT (vitamin-D loci).
Holick 2007; Engelsen 2010 (latitude & cutaneous vitamin-D synthesis).
"""

from __future__ import annotations

import pandas as pd


def _gt(df, rsid) -> str | None:
    if rsid not in df.index:
        return None
    row = df.loc[rsid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    g = row.get("genotype")
    if g is None:
        return None
    s = str(g).upper().replace(" ", "").replace("-", "")
    return s or None


# ─── 1. Circadian light timing ────────────────────────────────────────────────

def analyze_circadian_light(df) -> dict | None:
    gt = _gt(df, "rs1801260")   # CLOCK 3111T/C; C → evening lean
    if gt is None:
        return None
    n_c = gt.count("C")
    if n_c == 2:
        lean = "Evening-leaning"
        protocol = [
            "Get bright light within 30-60 min of waking — 20-40 min outdoors, "
            "or a 10,000-lux light box — to phase-advance your clock earlier.",
            "Dim household lights and cut blue light (screens on night mode / "
            "glasses) from ~2 h before target bedtime.",
            "If you need to shift earlier, 0.3-0.5 mg melatonin ~5-6 h before "
            "desired sleep acts as a chronobiotic (timing matters more than dose; "
            "this is NOT a bedtime sleeping-pill dose).",
            "Anchor wake time 7 days/week — evening types drift latest on free days.",
            "Schedule demanding cognitive work for late morning/early afternoon, "
            "not first thing.",
        ]
        detail = ("CLOCK 3111 C/C associates on average with an evening "
                  "chronotype. Your circadian system resists early schedules; "
                  "deliberate morning light + evening dark is the highest-"
                  "leverage lever to align with a conventional day.")
    elif n_c == 0:
        lean = "Morning-leaning"
        protocol = [
            "Protect your morning alertness peak for your most important work.",
            "Get outdoor light around midday to prevent an over-early evening "
            "melatonin rise and the associated early-evening energy crash.",
            "Be cautious with very early obligations that further advance you — "
            "morning types can drift into waking uncomfortably early.",
            "Keep evening light moderate but don't force late nights; your "
            "system is primed to sleep earlier.",
        ]
        detail = ("CLOCK 3111 T/T associates on average with a morning "
                  "chronotype. Your alertness front-loads; build your schedule "
                  "around an early peak rather than fighting it.")
    else:
        lean = "Intermediate chronotype"
        protocol = [
            "Consistent morning outdoor light (15-20 min) stabilises your rhythm.",
            "Consistent wake time is the single biggest lever for an "
            "intermediate type — pick one and hold it across the week.",
            "Cut bright/blue light in the last hour before bed.",
        ]
        detail = ("CLOCK 3111 T/C — an intermediate chronotype. You have "
                  "flexibility; consistency (not timing tricks) is your lever.")
    return {"gene": "CLOCK", "rsid": "rs1801260", "genotype": gt,
            "lean": lean, "detail": detail, "protocol": protocol,
            "citation": "Katzenberg 1998; Roenneberg 2007"}


# ─── 2. Exercise-modality fit ─────────────────────────────────────────────────

def analyze_exercise_modality(df) -> dict | None:
    actn3 = _gt(df, "rs1815739")   # C=R577 (power), T=X577 (endurance-lean)
    ace = _gt(df, "rs4343")        # ACE I/D proxy; G≈D (power), A≈I (endurance)
    if actn3 is None and ace is None:
        return None

    power_score = 0.0
    endurance_score = 0.0
    basis: list[dict] = []

    if actn3 is not None:
        n_r = actn3.count("C")  # R allele
        if n_r == 2:
            power_score += 1.0
            a_call = "ACTN3 R/R — functional α-actinin-3; fast-twitch/power favoured"
        elif n_r == 1:
            power_score += 0.5; endurance_score += 0.5
            a_call = "ACTN3 R/X — mixed fast/slow-twitch profile"
        else:
            endurance_score += 1.0
            a_call = "ACTN3 X/X — no α-actinin-3; endurance-leaning muscle profile"
        basis.append({"gene": "ACTN3", "rsid": "rs1815739", "genotype": actn3,
                      "call": a_call})
    if ace is not None:
        n_d = ace.count("G")   # G ≈ D allele (power/strength)
        if n_d == 2:
            power_score += 0.6
            e_call = "ACE D/D (proxy) — associated with power/strength & hypertrophy response"
        elif n_d == 1:
            power_score += 0.3; endurance_score += 0.3
            e_call = "ACE I/D (proxy) — mixed"
        else:
            endurance_score += 0.6
            e_call = "ACE I/I (proxy) — associated with endurance & aerobic efficiency"
        basis.append({"gene": "ACE", "rsid": "rs4343", "genotype": ace,
                      "call": e_call})

    if power_score > endurance_score + 0.25:
        lean = "Power / strength-leaning"
        emphasis = [
            "You likely see fast gains from strength, power and sprint work — "
            "explosive compound lifts, plyometrics, short high-intensity sprints.",
            "Keep 2-3 zone-2 cardio sessions/week for cardiovascular health — "
            "your genotype favours power but heart health needs aerobic base.",
            "Higher relative recovery need after heavy eccentric work — program "
            "adequate rest between power sessions.",
        ]
    elif endurance_score > power_score + 0.25:
        lean = "Endurance-leaning"
        emphasis = [
            "You likely adapt well to and recover quickly from endurance work — "
            "distance running/cycling/rowing, tempo and threshold sessions.",
            "Add 2× resistance training/week to preserve power, bone density and "
            "muscle mass — endurance genotypes still need strength for longevity.",
            "You can typically tolerate higher aerobic volume than power types.",
        ]
    else:
        lean = "Mixed / versatile"
        emphasis = [
            "Balanced fast/slow-twitch profile — you can succeed at both power "
            "and endurance; periodise between strength blocks and aerobic blocks.",
            "Concurrent training (mixing both) works well for you; use it.",
        ]
    return {"lean": lean, "power_score": round(power_score, 2),
            "endurance_score": round(endurance_score, 2), "basis": basis,
            "emphasis": emphasis,
            "caveat": ("Genetics shifts what responds *fastest*, not your ceiling. "
                       "Training, consistency and recovery dominate outcomes; "
                       "everyone needs both strength and cardio for health."),
            "citation": "Yang 2009 (ACTN3); Montgomery 1998 (ACE)"}


# ─── 3. Vitamin-D seasonality (genotype × latitude) ───────────────────────────

# SUNLIGHT-consortium low-25(OH)D effect alleles. Consumer-chip strand can flip
# these — treated as a modest tendency, not a hard call.
_VITD_LOW_ALLELE = {
    "rs2282679":  ("GC (DBP)", "G", "Wang 2010"),
    "rs10741657": ("CYP2R1", "A", "Jiang 2018"),
    "rs12785878": ("DHCR7/NADSYN1", "G", "Wang 2010"),
    "rs2228570":  ("VDR (FokI)", "A", "Uitterlinden 2004"),
}


def analyze_vitamin_d_seasonal(df, latitude: float = 40.0) -> dict | None:
    typed = []
    n_low = 0
    n_typed = 0
    for rsid, (gene, low_allele, cite) in _VITD_LOW_ALLELE.items():
        gt = _gt(df, rsid)
        if gt is None:
            continue
        n_typed += 1
        n = gt.count(low_allele)
        n_low += n
        typed.append({"gene": gene, "rsid": rsid, "genotype": gt,
                      "low_d_alleles": n, "citation": cite})
    if n_typed == 0:
        # Still give latitude guidance even with no genotype data.
        tendency = "unknown (no vitamin-D variants typed)"
    else:
        frac = n_low / (2 * n_typed)
        if frac >= 0.5:
            tendency = "Higher genetic tendency to lower baseline 25(OH)D"
        elif frac >= 0.25:
            tendency = "Moderate genetic tendency to lower baseline 25(OH)D"
        else:
            tendency = "Lower genetic tendency (baseline 25(OH)D less affected)"

    abs_lat = abs(latitude)
    # Cutaneous synthesis effectively ceases above ~37° from roughly Oct-Mar
    # (N hemisphere) due to solar zenith angle (Engelsen 2010; Holick 2007).
    if abs_lat >= 50:
        winter = ("At ≥50° latitude, cutaneous vitamin-D synthesis is negligible "
                  "for ~5-6 months (roughly October-March in the N hemisphere). "
                  "Year-round attention with definite winter supplementation.")
        months = "October–March (and often into April)"
    elif abs_lat >= 37:
        winter = ("At ~37-50° latitude, cutaneous synthesis effectively stops "
                  "for ~4-5 winter months. Winter supplementation is broadly "
                  "recommended regardless of genotype.")
        months = "November–February (extend if genetically lower)"
    elif abs_lat >= 23:
        winter = ("At subtropical latitude, some winter synthesis persists but "
                  "is reduced; supplementation still commonly warranted for "
                  "indoor lifestyles or darker skin.")
        months = "December–January (individual)"
    else:
        winter = ("Near the tropics, year-round synthesis is possible with "
                  "regular midday sun exposure; genotype/skin/lifestyle drive need.")
        months = "usually not needed from latitude alone"

    protocol = [
        f"Assumed latitude: {abs_lat:.0f}° — {winter}",
        "Get 10-30 min midday sun on arms/legs when the UV index ≥3 (spring-"
        "autumn) as your primary source; supplement in the low-UV months.",
        "Typical maintenance supplementation in deficit months is 1000-2000 "
        "IU/day D3 for most adults; confirm with a serum 25(OH)D test rather "
        "than guessing (target ~30-50 ng/mL).",
        "Take D3 with the largest fat-containing meal for absorption; pair with "
        "vitamin K2 and adequate magnesium if supplementing at the higher end.",
    ]
    if "Higher" in tendency:
        protocol.append(
            "Your genotype leans toward lower baseline 25(OH)D — bias toward "
            "the upper end of the maintenance range in winter and definitely "
            "test serum levels rather than assuming sun exposure suffices.")

    return {
        "latitude_assumed": abs_lat,
        "tendency": tendency,
        "n_typed": n_typed,
        "n_low_alleles": n_low,
        "variants": typed,
        "supplement_months": months,
        "protocol": protocol,
        "caveat": ("Consumer-chip strand conventions can flip vitamin-D effect "
                   "alleles, so the genetic tendency is modest evidence; the "
                   "latitude-based recommendation is the reliable part. A serum "
                   "25(OH)D test settles it definitively."),
        "citation": "Wang 2010; Jiang 2018 (SUNLIGHT); Holick 2007 (latitude)",
    }


# ─── Master ───────────────────────────────────────────────────────────────────

def analyze_environmental_optimization(df: pd.DataFrame,
                                       latitude: float = 40.0) -> dict:
    circadian = None
    exercise = None
    vitamin_d = None
    try:
        circadian = analyze_circadian_light(df)
    except Exception:
        pass
    try:
        exercise = analyze_exercise_modality(df)
    except Exception:
        pass
    try:
        vitamin_d = analyze_vitamin_d_seasonal(df, latitude=latitude)
    except Exception:
        pass

    return {
        "available": any([circadian, exercise, vitamin_d]),
        "latitude_assumed": abs(latitude),
        "circadian": circadian,
        "exercise": exercise,
        "vitamin_d": vitamin_d,
        "note": ("Actionable behavioural protocols from genotype. Latitude is an "
                 "explicit assumption (default 40°N); pass your real latitude for "
                 "precise vitamin-D seasonality."),
    }
