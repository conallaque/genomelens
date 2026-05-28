"""
Advanced Exercise Analytics
===========================

Composite athletic profile, ranked sport-matching, quantitative per-region
injury risk map, daily readiness formula, fully spelled-out sample workouts,
concurrent-training interference model, lactate-threshold estimate, tapering
and deload protocols, mental-skills profile (COMT-driven), cold/heat
adaptation, mobility prescription, and a 12-week plyometric progression.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple
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


# ── Composite Athletic Profile + Sport Match ────────────────────────────────

_SPORTS = [
    # (name, weights: power, endurance, recovery, hypertrophy, vo2max, pain, low-injury)
    ("Powerlifting",                {"power": 0.45, "hypertrophy": 0.30, "pain": 0.15, "low_injury": 0.10}),
    ("Olympic weightlifting",       {"power": 0.55, "hypertrophy": 0.15, "pain": 0.10, "low_injury": 0.20}),
    ("Sprinting (100-400m)",        {"power": 0.55, "vo2max": 0.10, "low_injury": 0.20, "pain": 0.15}),
    ("CrossFit / Hybrid",           {"power": 0.30, "endurance": 0.25, "hypertrophy": 0.20, "vo2max": 0.15, "pain": 0.10}),
    ("Distance running (5-10k)",    {"endurance": 0.45, "vo2max": 0.30, "low_injury": 0.15, "pain": 0.10}),
    ("Marathon / ultra",            {"endurance": 0.50, "vo2max": 0.25, "recovery": 0.10, "pain": 0.15}),
    ("Cycling (road)",              {"endurance": 0.45, "vo2max": 0.30, "power": 0.15, "low_injury": 0.10}),
    ("Triathlon",                   {"endurance": 0.40, "vo2max": 0.25, "recovery": 0.15, "low_injury": 0.10, "pain": 0.10}),
    ("Swimming",                    {"endurance": 0.30, "vo2max": 0.25, "power": 0.20, "low_injury": 0.25}),
    ("Rowing",                      {"endurance": 0.35, "power": 0.25, "vo2max": 0.20, "hypertrophy": 0.20}),
    ("Football (American)",         {"power": 0.40, "hypertrophy": 0.25, "vo2max": 0.10, "pain": 0.15, "low_injury": 0.10}),
    ("Soccer / football",           {"endurance": 0.30, "power": 0.25, "vo2max": 0.25, "low_injury": 0.20}),
    ("Basketball",                  {"power": 0.35, "vo2max": 0.25, "endurance": 0.20, "low_injury": 0.20}),
    ("Martial arts / boxing",       {"power": 0.30, "endurance": 0.25, "vo2max": 0.20, "pain": 0.25}),
    ("Rock climbing",               {"power": 0.30, "endurance": 0.25, "pain": 0.20, "low_injury": 0.25}),
    ("Yoga / flexibility-led",      {"recovery": 0.40, "low_injury": 0.40, "endurance": 0.20}),
]


def _normalise(v: float, lo: float = 0, hi: float = 100) -> float:
    return max(0, min(100, (v - lo) * 100 / (hi - lo)))


def composite_athletic_profile(base: Dict) -> Dict:
    """Synthesise base+advanced exercise data into normalised 0-100 attributes."""
    pe = base["power_endurance"]
    recovery_speed = base["recovery"]["speed"]
    injury_risks = base["injury_risk"]["risks"]
    vo2 = base.get("vo2max", {})
    strength = base.get("strength_trainability", {})
    pain = base.get("pain_tolerance", {})

    power = pe.get("ratio_pct_power", 50)
    endurance = pe.get("ratio_pct_endurance", 50)
    recovery_score = {"Fast": 85, "Moderate": 60, "Slow": 35}.get(recovery_speed, 50)
    hypertrophy_map = {
        "High hypertrophy responder": 85,
        "Average responder": 60,
        "Low hypertrophy / strength-leaning responder": 35,
    }
    hypertrophy = hypertrophy_map.get(strength.get("tier", ""), 50)
    vo2_map = {
        "High responder": 85, "Average responder": 60,
        "Low responder": 35, "Unknown": 50,
    }
    vo2max_score = vo2_map.get(vo2.get("tier", "Unknown"), 50)
    pain_score = 50
    if "Higher" in (pain.get("tolerance", "") or ""):
        pain_score = 80
    elif "Lower" in (pain.get("tolerance", "") or ""):
        pain_score = 30
    # Low-injury = inverse of how many risks are "Elevated"
    elev = sum(1 for r in injury_risks if "elev" in r["level"].lower())
    low_injury = max(20, 90 - 25 * elev)

    attrs = {
        "power": power,
        "endurance": endurance,
        "recovery": recovery_score,
        "hypertrophy": hypertrophy,
        "vo2max": vo2max_score,
        "pain": pain_score,
        "low_injury": low_injury,
    }

    # Score each sport
    sport_scores: List[Dict] = []
    for name, weights in _SPORTS:
        score = sum(weights.get(k, 0) * attrs.get(k, 50) for k in attrs)
        sport_scores.append({"sport": name, "fit_score": round(score, 1)})
    sport_scores.sort(key=lambda x: x["fit_score"], reverse=True)

    return {
        "attributes": attrs,
        "ranked_sports": sport_scores[:8],
        "bottom_sports": sport_scores[-3:],
        "overall_score": round(sum(attrs.values()) / len(attrs), 1),
        "summary": (
            f"Composite athletic index: {round(sum(attrs.values()) / len(attrs), 1)}/100. "
            f"Top-matched sport: {sport_scores[0]['sport']} (fit {sport_scores[0]['fit_score']}). "
            f"Strongest attribute: {max(attrs, key=attrs.get)} ({max(attrs.values()):.0f}). "
            f"Weakest: {min(attrs, key=attrs.get)} ({min(attrs.values()):.0f})."
        ),
    }


# ── Quantitative Body-Region Injury Risk Map ────────────────────────────────

def injury_risk_map(snps_df) -> Dict:
    col1a1 = _gt(snps_df, "rs1800012")
    col5a1 = _gt(snps_df, "rs12722")
    mmp3 = _gt(snps_df, "rs679620")
    vdr = _gt(snps_df, "rs2228570")
    gdf5 = _gt(snps_df, "rs143383")

    def region_score(base: int, *adds) -> int:
        return max(5, min(95, base + sum(adds)))

    knee_acl = region_score(
        20,
        20 if col1a1 and "T" in col1a1 else 0,
        10 if col1a1 and col1a1.count("T") == 2 else 0,
        5 if gdf5 and "T" in gdf5 else 0,
    )
    achilles = region_score(
        20,
        25 if col5a1 and col5a1.count("C") == 2 else (-10 if col5a1 and col5a1.count("T") == 2 else 0),
        10 if mmp3 and "A" in mmp3 else 0,
    )
    patellar = region_score(
        18,
        20 if col5a1 and col5a1.count("C") == 2 else 0,
        10 if col1a1 and "T" in col1a1 else 0,
    )
    rotator_cuff = region_score(
        20,
        12 if col1a1 and "T" in col1a1 else 0,
        8 if mmp3 and "A" in mmp3 else 0,
    )
    stress_fx = region_score(
        15,
        18 if vdr and "T" in vdr else 0,
        12 if col1a1 and "T" in col1a1 else 0,
    )
    lumbar_disc = region_score(
        20,
        15 if col1a1 and "T" in col1a1 else 0,
        10 if gdf5 and "T" in gdf5 else 0,
    )

    rows = [
        {"region": "Knee / ACL", "risk_pct": knee_acl,
         "interventions": "Nordic curls 2×/wk, single-leg RDLs, hop-stick drills, FMS screening."},
        {"region": "Achilles tendon", "risk_pct": achilles,
         "interventions": "Daily heel drops (eccentric), slow gradual mileage build, calf raises 3×/wk."},
        {"region": "Patellar tendon", "risk_pct": patellar,
         "interventions": "Spanish squats / decline squats 2×/wk, isometric wall sits, manage running volume jumps."},
        {"region": "Rotator cuff / shoulder", "risk_pct": rotator_cuff,
         "interventions": "Band external rotations, face pulls 100/week, scapular control, avoid kipping pullups early."},
        {"region": "Bone / stress fracture", "risk_pct": stress_fx,
         "interventions": "Vitamin D >40 ng/mL, calcium 1000 mg, heavy strength 2×/wk, mileage <10%/wk increase."},
        {"region": "Lumbar disc", "risk_pct": lumbar_disc,
         "interventions": "McGill big-3 daily (curl-up, side bridge, bird-dog), hip-hinge mastery before loading, no flexed-spine loading."},
    ]
    rows.sort(key=lambda r: r["risk_pct"], reverse=True)
    return {"regions": rows,
            "overall_index": round(sum(r["risk_pct"] for r in rows) / len(rows), 1)}


# ── Daily Readiness Formula ─────────────────────────────────────────────────

def readiness_formula(base: Dict) -> Dict:
    recovery_speed = base["recovery"]["speed"]
    slow = recovery_speed == "Slow"
    return {
        "formula": (
            "Daily readiness (0-100) = 0.30·HRV_norm + 0.25·Sleep_norm + "
            "0.15·RestingHR_norm + 0.15·SubjectiveEnergy + 0.10·Soreness_inverse + "
            "0.05·Mood. Each component normalised vs your 14-day rolling baseline."
        ),
        "thresholds": [
            {"range": "85-100", "action": "Green — push hard. PRs encouraged. Heavy day or VO2max."},
            {"range": "65-84", "action": "Yellow — normal session, hold RPE ≤8. Skip max efforts."},
            {"range": "50-64", "action": "Orange — reduce volume 30%. Easy aerobic or technique work only."},
            {"range": "<50", "action": "Red — full rest or 20 min walk. Forcing it now buys nothing."},
        ],
        "tracking": (
            "Daily HRV (rMSSD on waking via Whoop/Oura/HRV4Training), resting HR, "
            "subjective 1-5 sliders for energy/soreness/mood. Take a 14-day baseline "
            "before computing deviations. 3 consecutive days <baseline = mandatory deload."
        ),
        "personal_note": (
            "Slow-recovery genotype — be MORE conservative with thresholds. Push red "
            "to 55, treat 50-64 as orange." if slow else
            "Standard thresholds apply."
        ),
    }


# ── Detailed Sample Workouts ────────────────────────────────────────────────

def sample_workouts(base: Dict) -> List[Dict]:
    bias = base["power_endurance"]["bias"]
    fast_recovery = base["recovery"]["speed"] == "Fast"
    window = base["chronotype"].get("optimal_window", "anytime")

    if bias.startswith("Power"):
        return [
            {
                "name": "A — Lower Body Strength (Heavy)",
                "duration": "60-70 min",
                "best_window": window,
                "blocks": [
                    "Warm-up (10 min): bike 5 min, dynamic mobility, glute activation, 3 ramping squat sets",
                    "Back squat — 5×3 @ 85% 1RM, 3-min rest",
                    "Romanian deadlift — 4×6 @ 75%, 2-min rest",
                    "Bulgarian split squat — 3×8/leg @ moderate DBs, 90 sec",
                    "Nordic curl — 3×5 (eccentric 4 sec)",
                    "Standing calf raise — 3×12",
                    "Finisher: weighted plank 3×45 sec",
                    "Cool-down: 5 min walk + hip flexor + quad stretch",
                ],
                "rpe_target": "Top set RPE 8 (2 reps in reserve)",
            },
            {
                "name": "B — Power / Plyometric Day",
                "duration": "45-55 min",
                "best_window": window,
                "blocks": [
                    "Warm-up (10 min): pogo hops, A-skips, band activation",
                    "Box jump — 5×3 (max height attainable cleanly), 90 sec rest",
                    "Hang power clean — 5×3 @ 70%, 2-min rest",
                    "Push press — 4×4 @ 75%",
                    "Broad jump — 4×3, 90 sec",
                    "Med-ball rotational throw — 3×6/side",
                    "Cool-down: 5 min easy + ankle/hip mobility",
                ],
                "rpe_target": "Bar-speed focused — no grinding reps. RPE 7-8.",
            },
            {
                "name": "C — Conditioning (Power-Friendly)",
                "duration": "25-30 min",
                "best_window": "Anytime",
                "blocks": [
                    "Warm-up: 5 min easy bike",
                    "8 × 30-sec sprint @ 90% effort / 90-sec walk recovery",
                    "Cool-down: 5 min easy + breathing reset",
                ],
                "rpe_target": "Each sprint RPE 9; full rest between."
            },
        ]

    if bias.startswith("Endurance"):
        return [
            {
                "name": "A — Zone 2 Long Aerobic",
                "duration": "75-90 min",
                "best_window": "Anytime (morning common)",
                "blocks": [
                    "Warm-up: 10 min easy progression to Z2",
                    "Steady Z2 (60-70% HRmax) — 60-75 min nasal breathing or conversational pace",
                    "Cool-down: 5 min easy + 5 min mobility",
                ],
                "rpe_target": "RPE 4-5. If breathing through the nose is impossible, slow down.",
            },
            {
                "name": "B — Threshold Intervals",
                "duration": "60 min",
                "best_window": window,
                "blocks": [
                    "Warm-up: 15 min easy + 4×30 sec strides",
                    "Main: 4 × 8 min @ lactate threshold (RPE 7-8), 3-min easy between",
                    "Cool-down: 10 min easy",
                ],
                "rpe_target": "Threshold = 'comfortably hard' — sustainable for ~60 min if it were a race.",
            },
            {
                "name": "C — VO2max / Norwegian 4×4",
                "duration": "45 min",
                "best_window": window,
                "blocks": [
                    "Warm-up: 12 min easy + 4 ramping efforts",
                    "Main: 4 × 4 min @ 90-95% HRmax, 3-min easy jog between",
                    "Cool-down: 8 min easy",
                ],
                "rpe_target": "RPE 9 in last minute of each interval.",
            },
            {
                "name": "D — Strength Maintenance",
                "duration": "40 min",
                "best_window": "After easy days only",
                "blocks": [
                    "Warm-up: 10 min",
                    "Trap bar deadlift — 3×5 @ moderate",
                    "Walking lunge — 3×8/leg",
                    "Push-up or DB press — 3×8",
                    "Pull-up or row — 3×8",
                    "Plank + side plank — 3 rounds",
                ],
                "rpe_target": "RPE 7. Bone & connective-tissue insurance — don't fry yourself."
            },
        ]

    # Balanced/mixed
    return [
        {
            "name": "A — Full-Body Strength",
            "duration": "55 min",
            "best_window": window,
            "blocks": [
                "Warm-up: 8 min + activation",
                "Back squat — 4×5",
                "Bench press or DB press — 4×6",
                "Pendlay row — 4×6",
                "Trap-bar deadlift — 3×5",
                "Plank + hanging leg raise — 3 rounds",
            ],
            "rpe_target": "RPE 7-8.",
        },
        {
            "name": "B — Mixed Conditioning",
            "duration": "35 min",
            "best_window": "Anytime",
            "blocks": [
                "Warm-up: 5 min easy",
                "EMOM 20 min: odd min — 10 KB swings @ moderate; even min — 8 burpees",
                "Cool-down: 5 min walk + mobility",
            ],
            "rpe_target": "Steady; should finish strong, not staggering.",
        },
        {
            "name": "C — Long Aerobic + Strides",
            "duration": "55 min",
            "best_window": "Anytime",
            "blocks": [
                "45 min Z2 (hike, ride, or run)",
                "6 × 20-sec strides at end (relaxed fast running)",
                "Cool-down: 5 min mobility",
            ],
            "rpe_target": "Z2 RPE 4. Strides crisp but submaximal.",
        },
    ]


# ── Concurrent Training Interference Model ──────────────────────────────────

def concurrent_training_model(base: Dict) -> Dict:
    bias = base["power_endurance"]["bias"]
    recovery = base["recovery"]["speed"]
    if bias.startswith("Power") and recovery == "Slow":
        interference = "High"
        rule = ("Separate strength and endurance by ≥6 h on same day, or alternate days. "
                "Keep weekly endurance ≤ 90 min total when in a strength block.")
    elif bias.startswith("Power"):
        interference = "Moderate"
        rule = ("Do strength FIRST when concurrent. Cap endurance at ≤ 150 min/week during "
                "hypertrophy/strength blocks.")
    elif bias.startswith("Endurance"):
        interference = "Low"
        rule = ("Concurrent strength minimally hurts endurance and adds bone/injury "
                "resilience. 2 short strength sessions/wk after key endurance days is "
                "best.")
    else:
        interference = "Moderate"
        rule = ("Hybrid training works — sequence by daily goal. Hard endurance and hard "
                "strength on same day only if ≥4 h separation; otherwise alternate.")
    return {"interference_level": interference, "rule": rule}


# ── Lactate Threshold & VO2max Rough Estimates ──────────────────────────────

def aerobic_estimates(base: Dict) -> Dict:
    vo2_tier = base.get("vo2max", {}).get("tier", "Unknown")
    estimate = {
        "High responder": ("Strong genetic ceiling — well-trained endurance athletes from this profile reach VO2max 60+ ml/kg/min (men), 55+ (women).",
                           "Lactate threshold typically falls at ~85% of HRmax with training."),
        "Average responder": ("Trained ceiling typically VO2max 50-58 (men), 45-52 (women).",
                              "Lactate threshold ~80% HRmax with consistent threshold work."),
        "Low responder": ("Lower trained ceiling (~45-52 men, 40-48 women) — but threshold trainability is independent and often very high.",
                          "Lactate threshold improvements may outpace VO2max gains 2:1."),
    }
    vo2_text, lt_text = estimate.get(vo2_tier, ("Trainability unknown — assume average.",
                                                "Lactate threshold typically 75-85% HRmax."))
    return {
        "tier": vo2_tier,
        "vo2max_estimate": vo2_text,
        "lactate_threshold_estimate": lt_text,
        "test_protocols": [
            "30-30 sec ramp test or 2400m time trial → estimate VO2max",
            "30-min all-out time trial; average HR ≈ lactate threshold HR (Joe Friel)",
            "Lab gold standard: graded incremental treadmill / cycle ergometer + gas analysis",
        ],
    }


# ── Tapering Protocol ───────────────────────────────────────────────────────

def tapering_protocol(base: Dict) -> Dict:
    bias = base["power_endurance"]["bias"]
    if bias.startswith("Power"):
        return {
            "duration_days": 7,
            "volume_change": "−40% sets, intensity maintained",
            "structure": [
                "D-7: heavy single + light volume",
                "D-6: rest",
                "D-5: speed work, low volume",
                "D-4: rest",
                "D-3: opener — top set @ 90%, 2-3 sets only",
                "D-2: full rest",
                "D-1: short activation, no fatigue",
                "Competition day",
            ],
        }
    if bias.startswith("Endurance"):
        return {
            "duration_days": 14,
            "volume_change": "Week-1 −30%, week-2 −60%; intensity preserved",
            "structure": [
                "2 weeks out: cut volume 30%, keep one threshold + one VO2 session",
                "1 week out: cut volume 60%, replace with short race-pace efforts",
                "3 days out: easy 20-min jog + 4 strides",
                "2 days out: rest or 15-min walk",
                "1 day out: shake-out 10 min + 4 strides",
                "Race day",
            ],
        }
    return {
        "duration_days": 10,
        "volume_change": "−40% volume across the board",
        "structure": ["Cut volume 40% 10 days out; keep one short high-quality session per modality."],
    }


# ── Deload Protocol ─────────────────────────────────────────────────────────

def deload_protocol(base: Dict) -> Dict:
    return {
        "frequency": "Every 4-6 weeks (slow recovery: every 4; fast: every 6).",
        "structure": [
            "Volume −50% (cut sets in half, keep load).",
            "Intensity −20% on top sets (RPE 6-7 instead of 8-9).",
            "Cut highest-impact work (sprints, plyo, max-effort lifts) first.",
            "Sleep extension +30 min; emphasise pre-bed nutrition (40 g protein + carbs).",
            "Mobility/aerobic flush 2-3 sessions; full rest 1-2 days.",
        ],
        "trigger_signals": (
            "If you see: HRV ↓ >7% for 3+ days, RPE drift (same load feels harder), "
            "sleep quality dropping, joint pain accumulating, plateau in PRs ≥ 2 weeks → "
            "deload immediately, don't wait for the scheduled one."
        ),
    }


# ── Mental Skills Profile (COMT-driven) ─────────────────────────────────────

def mental_skills_profile(snps_df) -> Dict:
    comt = _gt(snps_df, "rs4680")
    drd2 = _gt(snps_df, "rs1800497")
    bdnf = _gt(snps_df, "rs6265")
    factors: List[str] = []
    profile = "Balanced"
    strategies: List[str] = []
    if comt:
        factors.append(f"rs4680 (COMT) {comt}")
        if "A" in comt and comt.count("A") == 2:
            profile = "Worrier (Met/Met)"
            strategies += [
                "Sharp focus under low/moderate stress, but performance drops sharply under high-pressure stimulus.",
                "Pre-event protocol: long warm-up to reduce novelty stress, breathing 4-7-8 box ×5 min, familiar cues only.",
                "Avoid caffeine spikes pre-competition — they amplify anxiety here. L-theanine 200 mg may help.",
                "Train under simulated pressure regularly so 'high stakes' becomes familiar territory.",
            ]
        elif "G" in comt and comt.count("G") == 2:
            profile = "Warrior (Val/Val)"
            strategies += [
                "Performance HOLDS UP or improves under high-pressure conditions — naturally clutch.",
                "May underperform in monotonous training — keep sessions varied or competitive.",
                "Caffeine helps focus; arousal upregulation strategies (loud music, aggressive cues) work.",
                "Risk: confidence can drift into recklessness — externally regulate workload.",
            ]
        else:
            profile = "Val/Met balanced"
            strategies += ["Flexible stress-performance curve — both calm and aroused states work."]
    if drd2 and "A" in drd2:
        strategies.append("DRD2 A1 — long-term motivation needs external structure (training partner, scheduled times, public log).")
    if bdnf and "A" in bdnf:
        strategies.append("BDNF Met — skill-acquisition is slower; use spaced practice (3×15 min beats 1×45 min for skill retention).")
    return {"profile": profile, "factors": factors or ["Mental-skills SNPs not typed"],
            "strategies": strategies or ["Standard mental-skills approach."]}


# ── Cold / Heat Adaptation ──────────────────────────────────────────────────

def thermal_adaptation(snps_df) -> Dict:
    ucp1 = _gt(snps_df, "rs1800592")
    ucp3 = _gt(snps_df, "rs1800849")
    adrb3 = _gt(snps_df, "rs4994")
    factors: List[str] = []
    cold = "Average"
    if ucp1:
        factors.append(f"rs1800592 (UCP1) {ucp1}")
        if "G" in ucp1: cold = "Reduced cold tolerance (lower brown-fat thermogenesis)"
    if ucp3:
        factors.append(f"rs1800849 (UCP3) {ucp3}")
    if adrb3:
        factors.append(f"rs4994 (ADRB3) {adrb3}")
        if "G" in adrb3: cold = "Reduced lipolytic cold response"
    return {
        "cold_tolerance": cold,
        "factors": factors or ["Thermal SNPs not typed"],
        "cold_protocol": (
            "Cold exposure 2-3 min @ 50-59 °F (10-15 °C) 3-5×/wk pre-breakfast — "
            "increases brown-fat activity over 4-6 weeks (Søberg protocol: 11 min/wk "
            "minimum across multiple sessions). Avoid IMMEDIATELY post-strength training "
            "(blunts hypertrophy adaptation 12+ h)."
        ),
        "heat_protocol": (
            "Sauna 80-90 °C, 20-30 min, 3-4×/wk post-exercise — heat-shock protein "
            "induction, plasma volume +7%, cardiac output adaptation. Endurance "
            "performance improves ~7% in trained athletes after 10 sauna sessions "
            "(Scoon et al). Rehydrate 1.5× sweat lost."
        ),
    }


# ── Mobility / Flexibility Protocol ─────────────────────────────────────────

def mobility_protocol(base: Dict) -> Dict:
    bias = base["power_endurance"]["bias"]
    elev_injury = any("elev" in r["level"].lower() for r in base["injury_risk"]["risks"])
    items = [
        "Hip flexor — couch stretch 2 min/side daily",
        "T-spine — open-book / quadruped rotation 2 min daily",
        "Ankle dorsiflexion — wall test + banded mobilisation 2 min/side, 3×/wk",
        "Shoulder — wall slides + 90/90 ER, 5 min, 3×/wk",
    ]
    if bias.startswith("Power"):
        items += ["Hip mobility (90/90 transitions, frog stretch) — 5 min daily — heavy hip work compresses ROM"]
    if bias.startswith("Endurance"):
        items += ["Calf/Achilles soft tissue (foam roll + lacrosse ball) 5 min, 3×/wk", "Hip flexor + glute med activation post-run"]
    if elev_injury:
        items += ["Eccentric tendon protocol (heel drops, Nordics, decline squat) — non-negotiable 3×/wk"]
    return {"daily_minutes": 10 if elev_injury else 7, "items": items,
            "guidance": "Mobility is a SKILL — daily small doses outperform weekly 30-min sessions. Pair to existing cues (post-shower, post-workout)."}


# ── Plyometric Progression ──────────────────────────────────────────────────

def plyometric_progression(base: Dict, injury_map: Dict) -> Dict:
    elev = any(r["risk_pct"] >= 50 for r in injury_map.get("regions", [])
               if r["region"] in ("Knee / ACL", "Patellar tendon", "Achilles tendon"))
    if elev:
        return {
            "ready": False,
            "preparation_weeks": "4-6 weeks of basic strength + isometrics before initiating jumps.",
            "phases": [
                "Weeks 1-2: Isometric loading — wall sit 5×45 sec, isometric squat hold, calf isometric.",
                "Weeks 3-4: Slow eccentric — tempo squat 4-1-0, Nordic curl 3×5, heel drop 3×12.",
                "Weeks 5-6: Reintroduce basics — pogo hops 3×10, low box step-off 3×5.",
                "Then re-assess; continue to phase 1 below.",
            ],
        }
    return {
        "ready": True,
        "phases": [
            {"weeks": "1-3", "level": "Foundational",
             "drills": "Pogo hops 3×10, line hops 3×20, ankle bounce 3×15, low box step-down 3×5/side"},
            {"weeks": "4-6", "level": "Bilateral jumps",
             "drills": "Box jump 4×3 (focus on soft landing), broad jump 4×3, squat jump 3×5"},
            {"weeks": "7-9", "level": "Unilateral / reactive",
             "drills": "Single-leg box jump 3×3/leg, lateral bound 3×4/side, depth jump (12 in box) 3×3"},
            {"weeks": "10-12", "level": "Advanced reactive",
             "drills": "Depth jump (18-24 in) 3×3, drop-to-broad-jump 3×3, hurdle hops 3×5"},
        ],
        "rules": [
            "Plyometrics first in session, before fatigue.",
            "Quality over quantity — stop a set the moment landing mechanics degrade.",
            "≥72 h between max-effort plyometric sessions.",
            "Total foot contacts/week: beginner <80, intermediate 80-150, advanced 150-250.",
        ],
    }


# ── Public synthesis API ────────────────────────────────────────────────────

def analyze_advanced_exercise(snps_df, base_result: Dict) -> Dict:
    profile = composite_athletic_profile(base_result)
    injury_map = injury_risk_map(snps_df)
    readiness = readiness_formula(base_result)
    workouts = sample_workouts(base_result)
    concurrent = concurrent_training_model(base_result)
    aerobic = aerobic_estimates(base_result)
    taper = tapering_protocol(base_result)
    deload = deload_protocol(base_result)
    mental = mental_skills_profile(snps_df)
    thermal = thermal_adaptation(snps_df)
    mobility = mobility_protocol(base_result)
    plyo = plyometric_progression(base_result, injury_map)
    return {
        "composite_profile": profile,
        "injury_risk_map": injury_map,
        "daily_readiness": readiness,
        "sample_workouts": workouts,
        "concurrent_training": concurrent,
        "aerobic_estimates": aerobic,
        "tapering": taper,
        "deload": deload,
        "mental_skills": mental,
        "thermal_adaptation": thermal,
        "mobility": mobility,
        "plyometric_progression": plyo,
    }
