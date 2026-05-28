"""
Exercise Protocols & Operational Tools
======================================

Concrete operational layers on top of the analytic core:

  • VO2max prediction & fitness-test calculators (Cooper, Rockport, beep)
  • HR-zone calculator (formula explained)
  • Running pace zones + race-time predictor (Riegel)
  • Cycling FTP zones
  • Sport-specific 12-week plans for top-matched sports
  • Lifting cue library for the main barbell lifts
  • Recovery modality matrix (quantified)
  • Master-athlete adjustments by decade
  • Pre-workout / post-workout supplement stacks
  • Movement-screen self-test
  • Detraining & retraining timeline
  • HRV-guided week structure
"""

from __future__ import annotations

from typing import Dict, List


# ── Fitness test calculators ───────────────────────────────────────────────

def fitness_test_calculators() -> Dict:
    return {
        "cooper_12min_run": {
            "formula": "VO2max (ml/kg/min) = (distance_m − 504.9) / 44.73",
            "interpretation_men_30y": {
                "Excellent (>2800m)": ">52", "Good (2400-2800m)": "44-52",
                "Average (2200-2400m)": "39-44", "Below (<2200m)": "<39",
            },
            "interpretation_women_30y": {
                "Excellent (>2300m)": ">45", "Good (2000-2300m)": "37-45",
                "Average (1800-2000m)": "31-37", "Below (<1800m)": "<31",
            },
            "how_to": "Warm up 10 min. Run as far as possible in exactly 12 min on flat track.",
        },
        "rockport_1mile_walk": {
            "formula": "VO2max = 132.853 − 0.0769·weight_lb − 0.3877·age + 6.315·(1 if male else 0) − 3.2649·time_min − 0.1565·HR_at_finish",
            "how_to": "Walk 1 mile (1.6 km) as briskly as possible. Record time + HR at finish.",
            "best_for": "De-conditioned, older athletes, or anyone unable to run.",
        },
        "beep_test_multistage": {
            "formula": "Level/shuttle reached → table lookup → VO2max",
            "how_to": "20 m shuttle runs, pace increases every minute. Continue to failure.",
            "best_for": "Team-sport athletes (soccer, basketball, rugby).",
        },
        "30min_time_trial_threshold": {
            "formula": "Lactate threshold HR ≈ average HR over final 20 min of 30-min all-out time trial (Friel).",
            "how_to": "30 min all-out solo TT on familiar route (running, rowing, or cycling). Wear HR monitor.",
        },
        "two_km_row": {
            "formula": "Erg score (time) → table for VO2max. World class men: <6:00 / women: <7:00.",
            "how_to": "All-out 2 km on Concept2. Damper 4-6, drag factor 120-130 men / 110-120 women.",
        },
        "ymca_step_test": {
            "formula": "3 min stepping on 12-in box at 24 steps/min; recovery HR for 1 min predicts VO2max.",
            "best_for": "Quick, equipment-light screen.",
        },
    }


# ── HR-Zone & Pace calculator (function user can call with age) ────────────

def hr_zones_for_age(age: int, resting_hr: int = 60) -> Dict:
    """Return BPM ranges per zone using Tanaka + Karvonen reserve method."""
    hrmax = round(208 - 0.7 * age)
    hrr = hrmax - resting_hr
    zones = []
    for label, lo_pct, hi_pct in [
        ("Z1 Recovery", 0.50, 0.60),
        ("Z2 Aerobic base", 0.60, 0.70),
        ("Z3 Tempo", 0.70, 0.80),
        ("Z4 Threshold", 0.80, 0.90),
        ("Z5 VO2max", 0.90, 1.00),
    ]:
        karv_lo = round(resting_hr + lo_pct * hrr)
        karv_hi = round(resting_hr + hi_pct * hrr)
        zones.append({"zone": label,
                      "pct_hrmax": f"{int(lo_pct*100)}-{int(hi_pct*100)}%",
                      "bpm_range_karvonen": f"{karv_lo}-{karv_hi}"})
    return {"age": age, "resting_hr": resting_hr, "HRmax_tanaka": hrmax,
            "HR_reserve": hrr, "zones": zones,
            "note": "Karvonen (using HR reserve) is more individually accurate than flat %HRmax."}


def running_pace_zones(threshold_pace_min_per_km: float) -> Dict:
    """Given lactate-threshold pace (min/km), return pace targets per zone."""
    t = threshold_pace_min_per_km  # in minutes per km
    return {
        "threshold_pace_min_km": round(t, 2),
        "zones": [
            {"zone": "Z1 Recovery", "pace_min_km": f"{t*1.30:.2f}-{t*1.40:.2f}"},
            {"zone": "Z2 Aerobic base", "pace_min_km": f"{t*1.15:.2f}-{t*1.25:.2f}"},
            {"zone": "Z3 Tempo (mara/half pace)", "pace_min_km": f"{t*1.05:.2f}-{t*1.10:.2f}"},
            {"zone": "Z4 Threshold (10k pace)", "pace_min_km": f"{t*0.97:.2f}-{t*1.03:.2f}"},
            {"zone": "Z5 VO2max (5k/3k pace)", "pace_min_km": f"{t*0.88:.2f}-{t*0.95:.2f}"},
            {"zone": "Z6 Anaerobic (1k/400m)", "pace_min_km": f"{t*0.78:.2f}-{t*0.85:.2f}"},
        ],
    }


def race_time_predictor(known_distance_km: float, known_time_min: float,
                        target_distance_km: float) -> Dict:
    """Riegel formula: T2 = T1 × (D2/D1)^1.06"""
    predicted = known_time_min * ((target_distance_km / known_distance_km) ** 1.06)
    return {
        "known": f"{known_distance_km} km in {known_time_min:.1f} min",
        "predicted_target": f"{target_distance_km} km in {predicted:.1f} min "
                            f"({predicted/target_distance_km:.2f} min/km pace)",
        "caveat": "Riegel assumes adequate endurance training for the longer distance. Less accurate beyond ~3× extrapolation.",
    }


def cycling_ftp_zones(ftp_watts: int) -> Dict:
    """Coggan power zones from FTP."""
    z = lambda lo, hi: f"{int(ftp_watts*lo)}-{int(ftp_watts*hi)} W"
    return {
        "ftp_watts": ftp_watts,
        "zones": [
            {"zone": "Z1 Active recovery", "watts": f"<{int(ftp_watts*0.55)}"},
            {"zone": "Z2 Endurance", "watts": z(0.55, 0.75)},
            {"zone": "Z3 Tempo", "watts": z(0.76, 0.90)},
            {"zone": "Z4 Lactate threshold", "watts": z(0.91, 1.05)},
            {"zone": "Z5 VO2max", "watts": z(1.06, 1.20)},
            {"zone": "Z6 Anaerobic", "watts": z(1.21, 1.50)},
            {"zone": "Z7 Neuromuscular", "watts": f">{int(ftp_watts*1.5)}"},
        ],
        "ftp_test_protocols": [
            "20-min all-out × 0.95 = FTP estimate (Allen & Coggan)",
            "Ramp test (Zwift/TrainerRoad): 1-min steps + 75% of final 1-min power",
            "8-min × 2 with 10-min recovery × 0.90",
        ],
    }


# ── Sport-Specific 12-Week Plans ───────────────────────────────────────────

def _plan_powerlifting() -> Dict:
    return {
        "sport": "Powerlifting",
        "weeks": [
            {"phase": "Hypertrophy (wk 1-4)",
             "schedule": "Squat 2×/wk, Bench 2×/wk, Deadlift 1×/wk. 4×8 @ 70-75%.",
             "accessory": "RDL, pause squats, close-grip bench, rows, walking lunges. 3×10-12."},
            {"phase": "Strength (wk 5-8)",
             "schedule": "Same lifts. 5×5 @ 80%, +2.5-5 kg/wk on top set.",
             "accessory": "Reduce volume; same exercises 3×6-8."},
            {"phase": "Peak (wk 9-11)",
             "schedule": "Top set 1-3 reps @ 87-93%, back-off 3×5 @ 80%.",
             "accessory": "Minimal — pause work, lockouts."},
            {"phase": "Taper (wk 12)",
             "schedule": "−40% volume, single @ 90% midweek, full rest 48 h pre-meet."},
        ],
        "meet_day": "1 light opener at 92% of PR, 2nd at PR, 3rd a +2.5 kg PR attempt.",
    }


def _plan_marathon() -> Dict:
    return {
        "sport": "Marathon (sub-3:30 target as example)",
        "weeks": [
            {"phase": "Base (wk 1-6)",
             "weekly_mileage_km": "40 → 60",
             "key_workouts": "Long run +1 km/wk to 25 km; 1 tempo 6-8 km; 1 easy 10-12 km. Strength 2×."},
            {"phase": "Build (wk 7-10)",
             "weekly_mileage_km": "60 → 80",
             "key_workouts": "Long 28-32 km with last 8 km @ race pace; threshold 4×8 min; "
                             "easy aerobic 14 km. Strength 1× (maintenance)."},
            {"phase": "Peak (wk 11)",
             "weekly_mileage_km": "65",
             "key_workouts": "Last 30 km long run with mid-section MP; VO2max 5×3 min."},
            {"phase": "Taper (wk 12-13)",
             "weekly_mileage_km": "45 then 30",
             "key_workouts": "Short race-pace efforts, 4-6 strides daily, full rest 48 h pre-race."},
        ],
        "race_day_fueling": (
            "Carb-load 8-10 g/kg in last 36 h. Race-morning breakfast 3 h pre: 1.5-2 g/kg "
            "carbs (oats + banana + honey, low-fat). Race fuelling: 60-90 g carbs/h "
            "(gels every 30-35 min) + water + electrolytes 500 mg Na/L."
        ),
    }


def _plan_crossfit_hybrid() -> Dict:
    return {
        "sport": "CrossFit / Hybrid",
        "weeks": [
            {"phase": "GPP base (wk 1-4)",
             "weekly_template": "M strength, T conditioning, W skill+aerobic, Th strength, "
                                "F conditioning, Sat long aerobic or sport, Sun off.",
             "focus": "Build aerobic base + perfect Olympic-lift technique at light loads."},
            {"phase": "Strength block (wk 5-8)",
             "weekly_template": "Compound strength 4×/wk, 1 metcon, 1 long Z2.",
             "focus": "Squat/deadlift/bench/press progression. Cap metcons at 2/wk."},
            {"phase": "Mixed modal (wk 9-11)",
             "weekly_template": "3 strength + 3 metcons + 1 long aerobic.",
             "focus": "Combine: e.g. heavy snatch + 1 km row + thrusters."},
            {"phase": "Competition prep (wk 12)",
             "weekly_template": "Deload + 1 simulated comp workout.",
             "focus": "Recover, sharpen mental game, scout judging standards."},
        ],
    }


def _plan_olympic_weightlifting() -> Dict:
    return {
        "sport": "Olympic Weightlifting",
        "weeks": [
            {"phase": "Technique block (wk 1-4)",
             "schedule": "Snatch + variants 3×/wk, C&J + variants 3×/wk, Squat 3×/wk. 70-80%.",
             "focus": "Bar path, receiving positions, foot work. Video every session."},
            {"phase": "Strength block (wk 5-8)",
             "schedule": "Same lift split. Squat 5×5 @ 80-85%. Olympic lifts 70-85%."},
            {"phase": "Intensity (wk 9-11)",
             "schedule": "Snatch + C&J at 85-95%, singles + doubles."},
            {"phase": "Peak/competition (wk 12)",
             "schedule": "Heavy single mid-week, deload, comp day."},
        ],
    }


def _plan_5k_runner() -> Dict:
    return {
        "sport": "5k Runner",
        "weeks": [
            {"phase": "Base (wk 1-4)",
             "schedule": "M off, T 6×400 @ 5k pace, W 8 km easy, Th tempo 4 km, F off/easy, "
                         "Sat 12 km long, Sun easy 6 km."},
            {"phase": "Build (wk 5-8)",
             "schedule": "Intervals 5×800 @ 5k, tempo 6 km, long 14 km, threshold 4×8 min."},
            {"phase": "VO2max (wk 9-11)",
             "schedule": "Norwegian 4×4, 12×400 @ 3k pace, tempo with surge intervals."},
            {"phase": "Sharpen + race (wk 12)",
             "schedule": "Short race-pace, 6×200, race day."},
        ],
    }


def sport_specific_plans(ranked_sports: List[Dict]) -> List[Dict]:
    table = {
        "Powerlifting": _plan_powerlifting,
        "Olympic weightlifting": _plan_olympic_weightlifting,
        "Marathon / ultra": _plan_marathon,
        "Distance running (5-10k)": _plan_5k_runner,
        "CrossFit / Hybrid": _plan_crossfit_hybrid,
        "Sprinting (100-400m)": _plan_powerlifting,
        "Cycling (road)": _plan_marathon,
        "Triathlon": _plan_marathon,
    }
    top = [s["sport"] for s in ranked_sports[:3]]
    plans = []
    for sport in top:
        gen = table.get(sport)
        if gen:
            plans.append(gen())
        else:
            plans.append({
                "sport": sport,
                "weeks": [
                    {"phase": "Base (wk 1-4)", "focus": "Build general fitness in sport-specific patterns."},
                    {"phase": "Build (wk 5-8)", "focus": "Increase intensity and sport-specific skill."},
                    {"phase": "Peak (wk 9-11)", "focus": "Competition-specific intensity, full speed/scenarios."},
                    {"phase": "Taper (wk 12)", "focus": "−40% volume, sharpen, recover."},
                ],
            })
    return plans


# ── Lifting Cue Library ────────────────────────────────────────────────────

def lifting_cues() -> List[Dict]:
    return [
        {"lift": "Back squat",
         "setup": "Bar across mid-traps. Brace 360°. Feet shoulder-width, slight toe-out 15-20°.",
         "cues": ["Spread the floor with feet (external rotation creates hip stability)",
                  "Knees track over middle toes — drive them out at the bottom",
                  "Sit between your hips, not back onto your heels",
                  "Maintain a tall chest — eyes forward, not up"],
         "common_faults": ["Knee valgus (caving in) — strengthen glute med",
                           "Butt wink (lumbar flexion at bottom) — improve hip mobility, reduce depth temporarily",
                           "Forward chest drop — strengthen upper back, work paused squats"]},
        {"lift": "Deadlift (conventional)",
         "setup": "Mid-foot under bar. Hinge to grip; shins touch bar. Squeeze armpits.",
         "cues": ["'Bend the bar around your shins' — lat engagement",
                  "Take the slack out of the bar before pulling (creates total-body tension)",
                  "Push the floor away — don't yank up",
                  "Hips and shoulders rise together — no early hip shoot"],
         "common_faults": ["Round upper back (acceptable if controlled, dangerous if uncontrolled)",
                           "Hyperextending at top — neutral spine, glute squeeze, no leaning back",
                           "Bar drifting forward — keep lats tight, bar against legs"]},
        {"lift": "Bench press",
         "setup": "Eyes under bar. Shoulder blades retracted + depressed. Arch lower back, feet planted.",
         "cues": ["'Bend the bar' — internal rotation cue for chest engagement",
                  "Elbows ~45-70° from torso (not 90°)",
                  "Touch bar at lower-sternum/upper-abdomen — not nipples",
                  "Drive feet through floor on press (leg drive)"],
         "common_faults": ["Flared elbows — shoulder impingement risk",
                           "Bouncing off chest — lose tension, eccentric value, and shoulder stability",
                           "Shoulders rolling forward at lockout — keep blades pinned"]},
        {"lift": "Overhead press",
         "setup": "Bar in front rack. Glutes + abs braced. Feet hip-width.",
         "cues": ["'Push your head through the window' once bar passes face",
                  "Squeeze glutes hard to prevent lumbar extension",
                  "Bar path straight up — bar finishes over middle of foot"],
         "common_faults": ["Excessive lower-back arch — soft brace, weaker core",
                           "Pressing forward — bar drifts away from base of support"]},
        {"lift": "Romanian deadlift (RDL)",
         "setup": "Stand tall with bar at hip. Soft-knee.",
         "cues": ["Push hips BACK (not down) — hinge, don't squat",
                  "Feel hamstrings stretch; bar slides along thighs",
                  "Stop when you feel a strong stretch — not at the floor"],
         "common_faults": ["Rounding lower back — reduce ROM, strengthen lats",
                           "Squatting instead of hinging — restart pattern with PVC pipe"]},
        {"lift": "Pull-up",
         "setup": "Hang from bar, shoulder-width grip. Shoulder blades pulled down/back BEFORE pulling.",
         "cues": ["'Pull elbows down toward back pockets'",
                  "Chest to bar, not chin over",
                  "Full hang at bottom (active hang, not dead hang) for full ROM"],
         "common_faults": ["Kipping when working strict — eliminate momentum",
                           "Shrugging into ears — re-set scaps each rep"]},
    ]


# ── Recovery Modality Matrix (quantified) ──────────────────────────────────

def recovery_modality_matrix() -> List[Dict]:
    return [
        {"modality": "Sleep 8-9 h",
         "effect_size": "★★★★★ — single highest-leverage intervention",
         "evidence": "Sleep <7 h drops max strength 4-7%, REM testosterone rises ~15% per extra hour to 8 h",
         "protocol": "Hard floor: never below 7 h on training weeks. Cool room (16-19 °C), dark, devices off 30 min pre-bed."},
        {"modality": "Carbs + protein post-workout",
         "effect_size": "★★★★",
         "evidence": "Glycogen resynthesis 50% faster with 1 g/kg carbs + 0.3 g/kg protein within 30 min vs 2 h delay",
         "protocol": "0.3 g/kg protein + 0.5-1 g/kg carbs within 1-2 h of finishing."},
        {"modality": "Hydration ≥1.5× sweat lost",
         "effect_size": "★★★★",
         "evidence": "2% dehydration = ~10% endurance performance drop; full rehydration requires sodium",
         "protocol": "Weigh before/after; replace 150% with electrolyte-containing fluid in 4 h."},
        {"modality": "Active recovery (Z1 20-30 min)",
         "effect_size": "★★★",
         "evidence": "Improves lactate clearance, reduces DOMS perception, no fitness cost",
         "protocol": "Easy bike/walk/swim 20-30 min day after hard session — RPE 3-4."},
        {"modality": "Sauna (Finnish, dry)",
         "effect_size": "★★★",
         "evidence": "Plasma volume +7% in 10 sessions (Scoon), endurance perf +7% trained athletes, cardiovascular adaptation",
         "protocol": "80-90 °C, 20-30 min, 3-4×/wk POST-exercise. Rehydrate."},
        {"modality": "Cold-water immersion",
         "effect_size": "★★ (paradoxical for hypertrophy)",
         "evidence": "Reduces inflammation/DOMS but blunts strength/hypertrophy adaptations 12+ h post-lift",
         "protocol": "Use only on heavy endurance/competition days. AVOID within 6 h of strength training."},
        {"modality": "Massage / soft-tissue work",
         "effect_size": "★★",
         "evidence": "Modest DOMS reduction; perceptual rather than mechanical benefit",
         "protocol": "1-2×/wk during heavy blocks. Foam roll daily 10 min self-massage."},
        {"modality": "Compression garments",
         "effect_size": "★",
         "evidence": "Small reduction in DOMS, no clear performance benefit",
         "protocol": "Optional — useful during long travel or back-to-back competition."},
        {"modality": "Stretching (static, post-session)",
         "effect_size": "★★",
         "evidence": "ROM gains, no acute DOMS reduction. Performed pre-strength acutely DECREASES force.",
         "protocol": "Static stretching post-workout only or in separate mobility session."},
        {"modality": "Tart cherry juice 12 oz",
         "effect_size": "★★",
         "evidence": "Reduces DOMS markers (CK, IL-6) and improves sleep onset",
         "protocol": "Pre-bed during heavy training weeks (consider sugar load)."},
        {"modality": "Magnesium 200-400 mg pre-bed",
         "effect_size": "★★",
         "evidence": "Improves sleep onset and quality; muscle relaxation",
         "protocol": "Glycinate or threonate form, 200-400 mg 30 min before bed."},
        {"modality": "Caffeine cycling (off 5-7 d pre-event)",
         "effect_size": "★★★ (for competition only)",
         "evidence": "Restores adenosine-receptor sensitivity; greater ergogenic effect on race day",
         "protocol": "Caffeine taper 7 days pre-event, then 3-6 mg/kg 45 min pre-race."},
        {"modality": "HRV-guided autoregulation",
         "effect_size": "★★★★",
         "evidence": "16-22% greater fitness gains in HRV-guided vs fixed-plan groups (Vesterinen 2016)",
         "protocol": "Track AM HRV; hard sessions only when HRV ≥ rolling mean; easy or rest below."},
    ]


# ── Master Athlete Adjustments by Decade ───────────────────────────────────

def master_athlete_adjustments() -> List[Dict]:
    return [
        {"decade": "30s",
         "key_changes": "Recovery slightly slower vs 20s. Mostly indistinguishable from younger.",
         "adjustments": ["Add 1 mandatory deload every 5-6 weeks",
                         "Prioritise sleep extension on hard days",
                         "Mobility 10 min/day non-negotiable"]},
        {"decade": "40s",
         "key_changes": "Testosterone -1% per year. Mitochondrial decline begins. Connective tissue stiffens.",
         "adjustments": ["Heavy strength 2×/wk minimum — counters sarcopenia",
                         "Joint-friendly modalities (rowing, cycling) substitute high-impact",
                         "Pre-workout mobility extended to 12-15 min",
                         "Protein target 1.8-2.0 g/kg (vs 1.6 for younger adults)"]},
        {"decade": "50s",
         "key_changes": "VO2max declines ~1% per year if untrained, ~0.5% if trained. Bone density drops post-menopause.",
         "adjustments": ["VO2max work 2×/wk preserves cardiac function (Joyner)",
                         "Heavy lifting (≥80% 1RM) 1-2×/wk — bone-loading",
                         "Recovery: 72 h between hard same-system sessions",
                         "Power output (jumps, throws, sprints) preserved with explicit speed work weekly"]},
        {"decade": "60s+",
         "key_changes": "Sarcopenia accelerates without strength training. Balance becomes critical.",
         "adjustments": ["Strength training 3×/wk (compound + accessories)",
                         "Balance work daily (single-leg work, eyes-closed stands)",
                         "Protein 2.0-2.4 g/kg to overcome anabolic resistance",
                         "Vitamin D ≥40 ng/mL, calcium 1200 mg — bone health",
                         "Cardio remains important: 150 min moderate + 2 vigorous sessions/wk"]},
    ]


# ── Pre / Post workout supplement stacks ───────────────────────────────────

def workout_supplement_stacks(caffeine_ergogenic: Dict, strength_tier: str,
                              vo2_tier: str) -> Dict:
    caff = caffeine_ergogenic.get("responder", "")
    return {
        "pre_workout": [
            (f"Caffeine 3-6 mg/kg 45 min pre" if caff.startswith("Strong")
             else f"Caffeine ≤1 mg/kg or skip" if caff.startswith("Null")
             else "Caffeine 2-3 mg/kg 45 min pre (test response)"),
            "Beta-alanine 3-5 g daily (loading 2-4 wk; benefits 1-4 min efforts)",
            "Creatine monohydrate 5 g daily (timing doesn't matter; total stack)",
            "Beetroot juice 500 ml 2-3 h pre — endurance + repeated sprints (nitric oxide)",
            "Citrulline malate 6-8 g 60 min pre — strength volume, blood flow",
        ],
        "intra_workout": [
            "Sessions <60 min: water only",
            "60-90 min: 30 g carbs/h (sports drink or gel)",
            "90+ min endurance: 60-90 g carbs/h + 500-750 ml fluid + 300-500 mg Na",
            "Hot conditions: + 200-300 mg sodium per L fluid above standard",
        ],
        "post_workout": [
            "Protein 0.3 g/kg within 1-2 h (whey/leucine ideal for MPS spike)",
            "Carbs 0.5-1 g/kg if glycogen-depleted",
            "Tart cherry juice 12 oz on heavy days (DOMS + sleep)",
            "Vitamin D 1000-2000 IU with the meal (if deficient)",
            "Magnesium 200-400 mg if going to bed within 4 h",
        ],
        "evidence_tiers": {
            "A_proven": ["Creatine", "Caffeine (responder genotype)", "Carbs intra-workout >90 min",
                         "Beta-alanine", "Beetroot/nitrate", "Protein post-workout"],
            "B_some_evidence": ["Citrulline", "L-theanine", "Tart cherry", "Sodium bicarbonate"],
            "C_avoid_or_unproven": ["BCAAs (whole protein superior)", "Glutamine for performance",
                                     "Most 'fat-burner' stacks", "Excess vitamin C/E (blunts adaptation)",
                                     "Testosterone-boosters (almost all useless)"],
        },
        "personalised_priority": (
            "Strength-focused: creatine + caffeine + protein. "
            "Endurance-focused: beetroot + carbs + caffeine + sodium. "
            f"Strength tier: {strength_tier}. VO2 tier: {vo2_tier}."
        ),
    }


# ── Movement Screen (FMS-style self test) ──────────────────────────────────

def movement_screen() -> List[Dict]:
    return [
        {"test": "Overhead squat (PVC overhead, feet shoulder-width)",
         "checks": "Heels stay down, knees track over feet, chest tall, bar over feet",
         "common_fail": "Heels rise / knees cave → ankle + hip mobility limits",
         "remedy": "Calf stretches, banded ankle mob, goblet squats with elevated heels"},
        {"test": "In-line lunge (heel-to-toe, bar vertical along spine)",
         "checks": "Vertical torso, knee tracks straight, no lateral shift",
         "common_fail": "Knee caves / lateral shift → glute med weak / hip mobility",
         "remedy": "Single-leg work, lateral band walks, hip airplane"},
        {"test": "Active straight-leg raise (lying, leg straight up)",
         "checks": "Raised leg passes hip line; opposite leg stays flat",
         "common_fail": "<70° → hamstring tightness / posterior chain restriction",
         "remedy": "Banded supine hamstring stretch, 90/90 work"},
        {"test": "Shoulder mobility (reach behind back, fingertip distance)",
         "checks": "Fingertips touch or within 1 fist-width",
         "common_fail": "Asymmetry > 1 hand-width → thoracic + shoulder mobility deficit",
         "remedy": "Wall slides, T-spine extensions, sleeper stretch"},
        {"test": "Trunk stability push-up (push-up with feet together, hands under shoulders)",
         "checks": "Body rises as one unit — no sag, no buckling",
         "common_fail": "Hip sag → core weakness",
         "remedy": "Plank progressions, dead bugs, bird-dogs (McGill big 3)"},
        {"test": "Rotary stability (quadruped: same-side hand+knee out, then under)",
         "checks": "Maintain spine alignment, controlled tempo",
         "common_fail": "Hip shift / rotation → anti-rotation core weak",
         "remedy": "Pallof presses, side planks with rotation"},
        {"test": "Single-leg balance (eyes open then closed, 30 sec each)",
         "checks": "30 sec eyes open easy; 20+ sec eyes closed",
         "common_fail": "<10 sec eyes closed → proprioceptive deficit, fall risk",
         "remedy": "Daily 1-min single-leg stand each leg, bosu work, tai chi"},
    ]


# ── Detraining / Retraining Timeline ───────────────────────────────────────

def detraining_timeline() -> List[Dict]:
    return [
        {"period": "Day 1-3 off", "what_happens": "Glycogen stores stay full; soreness clears. ZERO detraining."},
        {"period": "Week 1", "what_happens": "Plasma volume drops 5-10%; HR rises 5-10 bpm at given pace. Strength UNCHANGED."},
        {"period": "Week 2", "what_happens": "VO2max -4-6%. Strength still preserved (neural).  "
                                              "Hypertrophy maintained 2-3 weeks even without training."},
        {"period": "Weeks 3-4", "what_happens": "VO2max -10%. Strength begins to drop. Re-introduce 2× lifting/wk to preserve."},
        {"period": "Weeks 5-8", "what_happens": "Substantial detraining. VO2max -15-20%. Hypertrophy losses begin (but slow)."},
        {"period": "Beyond 12 weeks", "what_happens": "Detraining largely complete for endurance. Strength can hold longer (esp. with maintenance volume)."},
        {"period": "RETRAINING — muscle memory", "what_happens": "Myonuclei from prior hypertrophy persist 10+ years — return to previous size 2-3× faster than initial gains. Strength can return in 4-8 weeks vs months."},
    ]


# ── HRV-Guided Week Structure ──────────────────────────────────────────────

def hrv_guided_week() -> Dict:
    return {
        "tracking": (
            "AM HRV (rMSSD or proprietary) immediately on waking, daily. Track 14-day "
            "rolling mean and 7-day rolling mean. Daily reading interpreted relative to "
            "baseline."
        ),
        "decision_rules": [
            {"signal": "Today ≥ 7-day mean AND trending up",
             "action": "Green light — schedule HARDEST session this cycle."},
            {"signal": "Today within ±5% of 7-day mean",
             "action": "Normal — execute as planned."},
            {"signal": "Today < 7-day mean by 5-10%",
             "action": "Yellow — reduce volume 25%, intensity cap RPE 7."},
            {"signal": "Today < 7-day mean by >10% (1 day)",
             "action": "Orange — easy aerobic only (Z1 20-30 min)."},
            {"signal": "3 consecutive days < baseline",
             "action": "RED — mandatory full deload or rest week. Don't push through."},
        ],
        "context_modifiers": [
            "Alcohol, illness, poor sleep all crash HRV — diagnose cause before reducing training.",
            "Caffeine before measurement skews readings — measure fasted.",
            "Hormonal cycle (female): HRV dips luteal phase — normal, not a red flag.",
            "Travel + altitude crash HRV 24-72 h — give yourself a buffer."
        ],
        "implementation": (
            "Whoop, Oura, Garmin, HRV4Training, Elite HRV — all valid. Consistency matters "
            "more than device. Take 4 weeks to build baseline before acting on day-to-day "
            "deviations."
        ),
    }


# ── Public synthesis API ───────────────────────────────────────────────────

def analyze_exercise_protocols(result: Dict) -> Dict:
    profile = result.get("composite_profile", {})
    ranked = profile.get("ranked_sports", [])
    caff_erg = result.get("caffeine_ergogenic", {})
    strength_tier = result.get("strength_trainability", {}).get("tier", "")
    vo2_tier = result.get("vo2max", {}).get("tier", "")
    return {
        "fitness_tests": fitness_test_calculators(),
        "hr_zone_demo_age_35": hr_zones_for_age(35, 60),
        "running_pace_demo_threshold_5min_km": running_pace_zones(5.0),
        "race_time_demo": race_time_predictor(5, 24.0, 21.1),  # 5k 24:00 → half
        "cycling_ftp_demo": cycling_ftp_zones(250),
        "sport_specific_plans": sport_specific_plans(ranked) if ranked else [],
        "lifting_cues": lifting_cues(),
        "recovery_matrix": recovery_modality_matrix(),
        "master_athlete": master_athlete_adjustments(),
        "supplement_stacks": workout_supplement_stacks(caff_erg, strength_tier, vo2_tier),
        "movement_screen": movement_screen(),
        "detraining_timeline": detraining_timeline(),
        "hrv_guided_week": hrv_guided_week(),
    }
