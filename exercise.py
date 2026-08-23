"""
Personalised Exercise Programming
=================================

Reads ACTN3, ACE, COL1A1, COL5A1, PPARGC1A, IL6, CRP, BDNF, and CLOCK
genotypes and returns a training-prescription dict:

  • power_vs_endurance      — bias score & narrative
  • injury_risk             — tendon / ligament / bone with mitigation
  • recovery_speed          — fast / moderate / slow + protocol notes
  • optimal_training_window — chronotype-derived clock-time window
  • cognitive_exercise      — BDNF-specific cognitive-benefit notes
  • weekly_template         — example 7-day split tailored to the profile

The module degrades cleanly: each SNP is independently called, and missing
genotypes are skipped with explicit "not tested" notes rather than fabricated
defaults.
"""

from __future__ import annotations

import pandas as pd

try:
    from exercise_advanced import analyze_advanced_exercise
except ImportError:
    analyze_advanced_exercise = None
try:
    from exercise_protocols import analyze_exercise_protocols
except ImportError:
    analyze_exercise_protocols = None


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


# ── Power vs endurance ──────────────────────────────────────────────────────

def _analyze_power_endurance(snps_df) -> dict:
    actn3 = _gt(snps_df, "rs1815739")        # C=R (power), T=X (stop codon)
    ace = _gt(snps_df, "rs4646994") or _gt(snps_df, "rs4341")   # I/D proxy
    ppargc1a = _gt(snps_df, "rs8192678")     # G=Gly482, A=Ser482 (lower endurance trainability)

    power_score = 0
    endurance_score = 0
    factors: list[str] = []

    if actn3:
        if actn3.count("C") == 2:
            power_score += 2
            factors.append("ACTN3 CC (R/R) — sprint/power dominant fast-twitch fibres")
        elif "C" in actn3 and "T" in actn3:
            power_score += 1
            factors.append("ACTN3 CT (R/X) — mixed fibre type, balanced")
        elif actn3.count("T") == 2:
            endurance_score += 2
            factors.append("ACTN3 TT (X/X) — α-actinin-3 deficiency; endurance-biased")

    if ace and "G" in ace:                  # G ~ I allele on - strand reporting
        endurance_score += 1
        factors.append("ACE I-allele present — endurance-favouring vascular response")
    elif ace and "T" in ace:
        power_score += 1
        factors.append("ACE D-allele present — power/strength bias")

    if ppargc1a:
        if "G" in ppargc1a and ppargc1a.count("G") == 2:
            endurance_score += 1
            factors.append("PPARGC1A Gly482 — strong aerobic trainability")
        elif "A" in ppargc1a and ppargc1a.count("A") == 2:
            power_score += 0.5
            factors.append("PPARGC1A Ser482 — reduced mitochondrial biogenesis with training")

    if not factors:
        return {
            "bias": "Unknown",
            "ratio_pct_power": 50, "ratio_pct_endurance": 50,
            "factors": [],
            "recommendation": "No bias SNPs available — train balanced 50/50 power/endurance.",
        }

    total = power_score + endurance_score
    ratio_power = 50 if total == 0 else round(100 * power_score / total)
    ratio_endurance = 100 - ratio_power

    if ratio_power >= 70:
        bias = "Power-dominant"
        rec = (
            "Prioritise heavy strength (3-6 reps) and explosive work (jumps, sprints, "
            "Olympic lifts). 2-3 short conditioning sessions/week is sufficient. "
            "Avoid chronic high-volume aerobic which blunts your strength response."
        )
    elif ratio_power >= 55:
        bias = "Mixed, power-leaning"
        rec = (
            "Balanced split favouring strength (~4 sessions) over endurance "
            "(~2 sessions). Hybrid workouts (e.g. CrossFit, athletic prep) suit you."
        )
    elif ratio_power >= 45:
        bias = "Balanced"
        rec = (
            "Genuinely mixed phenotype — concurrent training works well. "
            "Programme polarised: 80% easy aerobic + 20% intense, with 2-3 lifts/week."
        )
    elif ratio_power >= 30:
        bias = "Mixed, endurance-leaning"
        rec = (
            "Favour aerobic base building (Zone 2 4-5×/week) with 2 strength sessions "
            "as injury-prevention insurance. Tempo work and threshold intervals respond well."
        )
    else:
        bias = "Endurance-dominant"
        rec = (
            "High aerobic-volume tolerance — long Zone 2 sessions, threshold work, "
            "and progressive interval training will yield large adaptations. "
            "Keep at least 1-2 strength sessions to preserve bone & power."
        )

    return {
        "bias": bias,
        "ratio_pct_power": ratio_power,
        "ratio_pct_endurance": ratio_endurance,
        "factors": factors,
        "recommendation": rec,
    }


# ── Injury risk ─────────────────────────────────────────────────────────────

def _analyze_injury_risk(snps_df) -> dict:
    col1a1 = _gt(snps_df, "rs1800012")       # Sp1 site; T allele = soft-tissue risk
    col5a1 = _gt(snps_df, "rs12722")         # CC vs TT; T = better tendon, C = stiffer
    mmp3 = _gt(snps_df, "rs679620")
    risks: list[dict] = []

    if col1a1 and "T" in col1a1:                    # one or two T copies
        risks.append({
            "tissue": "Ligaments / ACL",
            "level": "Elevated" if col1a1.count("T") == 2 else "Moderate",
            "marker": f"rs1800012 (COL1A1 Sp1) {col1a1}",
            "mitigation": (
                "Prioritise neuromuscular knee-control drills (Nordics, single-leg "
                "RDL, balance work). Collagen + vitamin C 1 hour pre-training "
                "(15 g + 50 mg) improves tendon collagen synthesis."
            ),
        })

    if col5a1:
        if "C" in col5a1 and col5a1.count("C") == 2:
            risks.append({
                "tissue": "Achilles / patellar tendon",
                "level": "Elevated",
                "marker": "rs12722 (COL5A1) CC",
                "mitigation": (
                    "Higher tendinopathy risk profile. Build eccentric calf/quad capacity "
                    "(slow heel drops, decline squats). Avoid sudden volume jumps in running."
                ),
            })
        elif "T" in col5a1 and col5a1.count("T") == 2:
            risks.append({
                "tissue": "Achilles / patellar tendon",
                "level": "Protective",
                "marker": "rs12722 (COL5A1) TT",
                "mitigation": "Favourable tendon-stiffness profile — running economy advantage.",
            })

    if mmp3 and "A" in mmp3:
        risks.append({
            "tissue": "Tendon (MMP3 remodelling)",
            "level": "Mildly elevated",
            "marker": f"rs679620 (MMP3) {mmp3}",
            "mitigation": "Avoid sudden training-load increases; deload weeks every 4-6 weeks.",
        })

    if not risks:
        risks.append({
            "tissue": "Soft tissue",
            "level": "Baseline",
            "marker": "—",
            "mitigation": "No risk-elevating tissue-genetics variants typed.",
        })

    return {"risks": risks}


# ── Recovery / inflammation ─────────────────────────────────────────────────

def _analyze_recovery(snps_df) -> dict:
    il6 = _gt(snps_df, "rs1800795")
    crp = _gt(snps_df, "rs2794520") or _gt(snps_df, "rs1205")
    sod2 = _gt(snps_df, "rs4880")

    inflam_score = 0
    notes: list[str] = []

    if il6:
        if "G" in il6 and il6.count("G") == 2:
            inflam_score += 2
            notes.append(f"rs1800795 (IL6 -174) {il6} — high baseline inflammatory tone")
        elif "G" in il6:
            inflam_score += 1
            notes.append(f"rs1800795 (IL6 -174) {il6} — moderate inflammatory tendency")

    if crp and "C" in crp and crp.count("C") == 2:
        inflam_score += 1
        notes.append(f"CRP genotype {crp} — higher baseline CRP")

    if sod2 and "G" in sod2:
        notes.append(f"rs4880 (SOD2 Ala16Val) {sod2} — reduced mitochondrial antioxidant capacity")

    if inflam_score >= 2:
        speed = "Slow"
        protocol = (
            "Expect prolonged DOMS and slower clearance of training stress. "
            "48-72 h between hard sessions on the same muscle group. Aggressive "
            "recovery: sleep ≥8 h, omega-3, tart cherry, contrast showers."
        )
    elif inflam_score == 1:
        speed = "Moderate"
        protocol = (
            "Typical recovery curves. 24-48 h between hard same-muscle sessions. "
            "Standard sleep, nutrition, and active-recovery sufficient."
        )
    else:
        speed = "Fast"
        protocol = (
            "Below-average baseline inflammation — tolerate higher training "
            "frequency. Hard sessions back-to-back are feasible with adequate "
            "fuelling and sleep."
        )

    return {
        "speed": speed,
        "inflammation_score": inflam_score,
        "factors": notes or ["No inflammation SNPs typed"],
        "protocol": protocol,
    }


# ── Chronotype / training window ────────────────────────────────────────────

def _analyze_chronotype(snps_df) -> dict:
    clock = _gt(snps_df, "rs1801260")       # T=evening, C=morning (3111C)
    per3 = _gt(snps_df, "rs228697")

    factors: list[str] = []
    score = 0  # negative = morning, positive = evening
    if clock:
        if "T" in clock and clock.count("T") == 2:
            score += 2
            factors.append(f"rs1801260 (CLOCK 3111) {clock} — strong evening chronotype")
        elif "T" in clock:
            score += 1
            factors.append(f"rs1801260 (CLOCK 3111) {clock} — mild evening preference")
        else:
            score -= 1
            factors.append(f"rs1801260 (CLOCK 3111) {clock} — morning preference")
    if per3 and "C" in per3:
        score += 1
        factors.append(f"rs228697 (PER3) {per3} — evening-shifted circadian")

    if score >= 2:
        chronotype = "Evening"
        window = "16:00 – 20:00"
        rationale = (
            "Late-chronotype genotype — peak strength, power, and reaction-time "
            "performance shifts later in the day. Avoid 06:00 sessions if performance "
            "matters; reserve early hours for low-intensity work."
        )
    elif score >= 1:
        chronotype = "Slight evening"
        window = "15:00 – 19:00"
        rationale = "Mild evening bias. Best output mid-afternoon to early evening."
    elif score <= -1:
        chronotype = "Morning"
        window = "07:00 – 11:00"
        rationale = (
            "Morning chronotype — strength and cognitive output peak in the first half "
            "of the day. Long-duration aerobic also tolerated AM."
        )
    else:
        chronotype = "Neutral / unknown"
        window = "Anytime — train at convenient time"
        rationale = "No clear chronotype signal from typed SNPs."

    return {
        "chronotype": chronotype,
        "optimal_window": window,
        "factors": factors,
        "rationale": rationale,
    }


# ── BDNF / cognitive exercise benefits ──────────────────────────────────────

def _analyze_cognitive_exercise(snps_df) -> dict:
    bdnf = _gt(snps_df, "rs6265")           # Val66Met; A = Met (reduced activity-dependent secretion)
    if not bdnf:
        return {
            "genotype": None,
            "summary": "BDNF Val66Met (rs6265) not typed.",
            "training_note": (
                "Aerobic exercise still produces ~2-fold BDNF elevations in most "
                "people — keep 3+ aerobic sessions/week for cognitive benefit."
            ),
        }
    if "A" in bdnf and bdnf.count("A") == 2:
        return {
            "genotype": "Met/Met",
            "summary": (
                "BDNF Met/Met — significantly reduced activity-dependent BDNF release. "
                "Exercise is especially important; cognitive-enhancement effects of "
                "aerobic training are real but slightly muted vs Val/Val."
            ),
            "training_note": (
                "5+ aerobic sessions/week (Zone 2, 45+ min) maximise BDNF compensation. "
                "Skill-based / coordinative training (dance, martial arts, climbing) "
                "produces additional BDNF response. Avoid prolonged inactivity."
            ),
        }
    if "A" in bdnf:
        return {
            "genotype": "Val/Met",
            "summary": (
                "BDNF Val/Met — intermediate activity-dependent secretion. "
                "Cognitive benefits of exercise are robust."
            ),
            "training_note": (
                "Mix steady-state aerobic with 1-2 weekly high-intensity sessions "
                "(HIIT). Both formats elevate BDNF; combination beats either alone."
            ),
        }
    return {
        "genotype": "Val/Val",
        "summary": (
            "BDNF Val/Val — efficient activity-dependent BDNF secretion. Strong "
            "cognitive-benefit response to aerobic exercise."
        ),
        "training_note": (
            "Conventional aerobic training (3-4 ×/week, 40-60 min Zone 2) is highly "
            "effective for mood, memory, and neuroplasticity outcomes."
        ),
    }


# ── Weekly template synthesis ───────────────────────────────────────────────

def _build_weekly_template(power_endurance: dict, recovery: dict, chronotype: dict) -> list[dict]:
    bias = power_endurance["bias"]
    fast_recovery = recovery["speed"] == "Fast"
    days: list[dict] = []
    window = chronotype.get("optimal_window", "")

    if bias.startswith("Power"):
        days = [
            {"day": "Mon", "session": "Lower body strength (squat, deadlift, jumps)"},
            {"day": "Tue", "session": "Conditioning — 20 min Zone 2 + sprint intervals"},
            {"day": "Wed", "session": "Upper body strength (press, pull, accessories)"},
            {"day": "Thu", "session": "Recovery — mobility, walk"},
            {"day": "Fri", "session": "Full-body power (cleans, push press, plyo)"},
            {"day": "Sat", "session": "Athletic / sport-specific play"},
            {"day": "Sun", "session": "Off / mobility"},
        ]
    elif bias.startswith("Endurance"):
        days = [
            {"day": "Mon", "session": "Zone 2 endurance 60 min"},
            {"day": "Tue", "session": "Threshold intervals (4×8 min)"},
            {"day": "Wed", "session": "Zone 2 endurance 45 min + strength accessories"},
            {"day": "Thu", "session": "Recovery / easy spin"},
            {"day": "Fri", "session": "VO2max intervals (5×3 min)"},
            {"day": "Sat", "session": "Long Zone 2 90+ min"},
            {"day": "Sun", "session": "Off"},
        ]
    else:
        days = [
            {"day": "Mon", "session": "Full-body strength A"},
            {"day": "Tue", "session": "Zone 2 endurance 45 min"},
            {"day": "Wed", "session": "Full-body strength B"},
            {"day": "Thu", "session": "HIIT (e.g. 10×1 min hard / 1 min easy)"},
            {"day": "Fri", "session": "Mobility / yoga"},
            {"day": "Sat", "session": "Long mixed session — hike, ride, or sport"},
            {"day": "Sun", "session": "Off"},
        ]

    if not fast_recovery and len(days) > 0:
        # Insert extra recovery if slow inflammation profile
        days[3]["session"] = "Full recovery day (slow inflammation profile)"

    # Annotate timing window
    for d in days:
        if "Off" not in d["session"] and "Recovery" not in d["session"] and "mobility" not in d["session"].lower():
            d["time"] = window
        else:
            d["time"] = "Anytime"

    return days


# ── VO2max trainability ─────────────────────────────────────────────────────

def _analyze_vo2max_trainability(snps_df) -> dict:
    vegfa = _gt(snps_df, "rs2010963")
    hif1a = _gt(snps_df, "rs11549465")
    adrb2 = _gt(snps_df, "rs1042713")
    nrf2 = _gt(snps_df, "rs7181866")
    factors: list[str] = []
    score = 0
    if vegfa:
        factors.append(f"rs2010963 (VEGFA) {vegfa}")
        if "G" in vegfa:
            score += vegfa.count("G")
    if hif1a:
        factors.append(f"rs11549465 (HIF1A) {hif1a}")
        if "T" in hif1a:
            score += 1
    if adrb2:
        factors.append(f"rs1042713 (ADRB2) {adrb2}")
        if "A" in adrb2:
            score += 1
    if nrf2:
        factors.append(f"rs7181866 (NRF2) {nrf2}")
        if "G" in nrf2:
            score += 1
    if score >= 3:
        tier = "High responder"
        guidance = (
            "Strong VO2max trainability — expect +15–20% with 8 weeks of polarised "
            "training. Push interval volume: 4×4 min @ 90–95% HRmax, 1–2×/week."
        )
    elif score >= 1:
        tier = "Average responder"
        guidance = (
            "Typical VO2max response (~+10%) with structured aerobic+interval work. "
            "Combine 2 Zone-2 + 1 VO2max + 1 threshold session/week."
        )
    elif factors:
        tier = "Low responder"
        guidance = (
            "Genetic VO2max-response cluster on the low end. Don't be discouraged — "
            "shift focus to lactate-threshold/economy gains (which improve regardless "
            "of VO2max) and strength-based aerobic transfer (heavy carries, hill work)."
        )
    else:
        tier = "Unknown"
        guidance = "Trainability SNPs not typed — assume average response."
    return {"tier": tier, "score": score, "factors": factors or ["VO2max SNPs not typed"],
            "guidance": guidance}


# ── Strength trainability (MSTN / IGF1 / ACVR1B) ────────────────────────────

def _analyze_strength_trainability(snps_df) -> dict:
    mstn = _gt(snps_df, "rs1805086")
    igf1 = _gt(snps_df, "rs35767")
    actn3 = _gt(snps_df, "rs1815739")
    factors: list[str] = []
    score = 0
    if mstn:
        factors.append(f"rs1805086 (MSTN) {mstn}")
        if "A" in mstn:
            score += 2
    if igf1:
        factors.append(f"rs35767 (IGF1) {igf1}")
        if "A" in igf1:
            score += 1
    if actn3 and "C" in actn3:
        score += actn3.count("C") * 0.5
    if score >= 2:
        tier = "High hypertrophy responder"
        rec = (
            "Strong genetic hypertrophy response — push moderate-rep volume (8–12 reps, "
            "12–20 sets/muscle group/week). Expect visible gains within 8–12 weeks."
        )
    elif score >= 1:
        tier = "Average responder"
        rec = (
            "Typical hypertrophy response — 10–15 sets/muscle/week at 6–12 reps, "
            "progressive overload weekly."
        )
    else:
        tier = "Low hypertrophy / strength-leaning responder"
        rec = (
            "Hypertrophy genes lean low — focus on strength (3–6 reps, heavy) where "
            "neural adaptations drive bigger early gains, plus high-volume blocks "
            "(8–12 weeks) for body composition. Consistency beats genes."
        )
    return {"tier": tier, "score": score, "factors": factors or ["Strength SNPs not typed"],
            "recommendation": rec}


# ── Fat-loss response to exercise (FTO + ADRB2/ADRB3) ───────────────────────

def _analyze_fat_loss_response(snps_df) -> dict:
    fto = _gt(snps_df, "rs9939609")
    adrb2_27 = _gt(snps_df, "rs1042714")
    adrb3 = _gt(snps_df, "rs4994")
    factors: list[str] = []
    fto_risk = False
    if fto:
        factors.append(f"rs9939609 (FTO) {fto}")
        if "A" in fto:
            fto_risk = True
    if adrb2_27:
        factors.append(f"rs1042714 (ADRB2 Q27E) {adrb2_27}")
    if adrb3:
        factors.append(f"rs4994 (ADRB3 Trp64Arg) {adrb3}")
    if fto_risk:
        guidance = (
            "FTO obesity-risk allele — GOOD NEWS: large studies show ≥150 min/week of "
            "moderate aerobic exercise NEUTRALISES the FTO weight-gain effect. Aerobic "
            "exercise is non-negotiable. Strength training adds compounding benefit by "
            "raising RMR. Combine with high-protein/high-fibre eating for best fat-loss "
            "outcomes."
        )
    else:
        guidance = (
            "No FTO risk allele — standard exercise dose for fat loss applies "
            "(150–300 min/week moderate aerobic + 2× strength)."
        )
    return {"fto_risk_allele": fto_risk, "factors": factors or ["Fat-loss SNPs not typed"],
            "guidance": guidance}


# ── Pain tolerance (COMT, OPRM1) ────────────────────────────────────────────

def _analyze_pain_tolerance(snps_df) -> dict:
    comt = _gt(snps_df, "rs4680")
    oprm1 = _gt(snps_df, "rs1799971")
    factors: list[str] = []
    tolerance = "Average"
    if comt:
        factors.append(f"rs4680 (COMT Val158Met) {comt}")
        if "A" in comt and comt.count("A") == 2:
            tolerance = "Lower (Met/Met — 'worrier')"
        elif "G" in comt and comt.count("G") == 2:
            tolerance = "Higher (Val/Val — 'warrior')"
    if oprm1:
        factors.append(f"rs1799971 (OPRM1) {oprm1}")
        if "G" in oprm1:
            tolerance += " · OPRM1 G — higher pain sensitivity"
    if "Lower" in tolerance:
        guidance = (
            "Lower pain threshold — RPE will feel harder than HR/power data suggests. "
            "Use objective metrics (HR, pace, bar speed) not 'how hard it feels'. Build "
            "RPE tolerance gradually with shorter intervals; mental-skill work "
            "(breathing, self-talk) yields disproportionate gains."
        )
    elif "Higher" in tolerance:
        guidance = (
            "Higher pain tolerance — risk of overtraining because you can push through "
            "warning signs. USE objective recovery metrics (HRV, resting HR, sleep) and "
            "schedule mandatory deload weeks."
        )
    else:
        guidance = "Typical pain perception — RPE scales reliably."
    return {"tolerance": tolerance, "factors": factors or ["Pain SNPs not typed"],
            "guidance": guidance}


# ── Motivation / adherence (DRD2, COMT, dopamine) ───────────────────────────

def _analyze_motivation(snps_df) -> dict:
    drd2 = _gt(snps_df, "rs1800497")  # A1 allele = lower D2 receptor density
    comt = _gt(snps_df, "rs4680")
    factors: list[str] = []
    low_drive = False
    if drd2:
        factors.append(f"rs1800497 (DRD2 Taq1A) {drd2}")
        if "A" in drd2:
            low_drive = True
    if comt and "G" in comt and comt.count("G") == 2:
        factors.append("rs4680 (COMT Val/Val) — fast dopamine clearance")
    if low_drive:
        guidance = (
            "Reduced D2 receptor density — exercise-induced 'reward' feels muted and "
            "habit formation takes longer. Strategies: (1) train at fixed time/cue daily "
            "(habit stack), (2) use external accountability (training partner, coach, "
            "public log), (3) chase variety (novel sports/classes hit dopamine harder "
            "than steady-state cardio), (4) explicit short-term rewards post-session "
            "(post-workout coffee, podcast). Expect 90+ days to lock in vs 60 for "
            "average."
        )
    else:
        guidance = (
            "Typical dopaminergic reward — exercise habit usually forms within 6–8 "
            "weeks of consistency. Use that 'runner's high' window."
        )
    return {"low_dopamine_reward": low_drive, "factors": factors or ["DRD2/COMT not typed"],
            "guidance": guidance}


# ── Caffeine ergogenic response (CYP1A2 for performance) ────────────────────

def _analyze_caffeine_ergogenic(snps_df) -> dict:
    cyp1a2 = _gt(snps_df, "rs762551")
    factors: list[str] = []
    if not cyp1a2:
        return {"responder": "Unknown",
                "factors": ["CYP1A2 not typed"],
                "guidance": "Trial 3 mg/kg 45 min pre-exercise to test response."}
    factors.append(f"rs762551 (CYP1A2) {cyp1a2}")
    if "A" in cyp1a2 and cyp1a2.count("A") == 2:
        return {
            "responder": "Strong ergogenic responder",
            "factors": factors,
            "guidance": (
                "Fast metaboliser — caffeine reliably improves endurance/power. "
                "Pre-event: 3–6 mg/kg body weight 45–60 min before competition. "
                "Cycle caffeine off 5–7 days before key events to restore sensitivity."
            ),
        }
    if "C" in cyp1a2:
        return {
            "responder": "Null / negative responder",
            "factors": factors,
            "guidance": (
                "Slow metaboliser — high-dose pre-workout caffeine likely impairs "
                "rather than enhances endurance and elevates BP. Stick to ≤1 mg/kg "
                "if any; consider beta-alanine, beetroot juice, or no stimulant."
            ),
        }
    return {"responder": "Intermediate", "factors": factors,
            "guidance": "Modest pre-exercise dose (2–3 mg/kg) — test individual response."}


# ── Sleep need / recovery quality ───────────────────────────────────────────

def _analyze_sleep(snps_df) -> dict:
    ada = _gt(snps_df, "rs73598374")
    per3_vntr = _gt(snps_df, "rs57875989")
    factors: list[str] = []
    deep_sleep = "Average"
    if ada:
        factors.append(f"rs73598374 (ADA) {ada}")
        if "T" in ada:
            deep_sleep = "Deeper slow-wave sleep (efficient recoverer)"
    if per3_vntr:
        factors.append(f"rs57875989 (PER3 VNTR) {per3_vntr}")
    guidance = (
        "Standard 7.5–9 h target. Sleep-extension protocol on hard training weeks: "
        "+30–60 min vs baseline. Track HRV trend; <baseline 3 days = mandatory deload."
    )
    if "Deeper" in deep_sleep:
        guidance += " ADA T allele — slow-wave sleep is unusually efficient; 7 h may feel adequate but still target 8 h for athletic recovery."
    return {"sleep_phenotype": deep_sleep, "factors": factors or ["Sleep SNPs not typed"],
            "guidance": guidance}


# ── Iron-endurance interaction (HFE) ────────────────────────────────────────

def _analyze_iron_endurance(snps_df) -> dict:
    c282y = _gt(snps_df, "rs1800562")
    h63d = _gt(snps_df, "rs1799945")
    factors: list[str] = []
    risk = False
    if c282y and "A" in c282y:
        factors.append(f"rs1800562 (HFE C282Y) {c282y}")
        risk = True
    if h63d and "G" in h63d:
        factors.append(f"rs1799945 (HFE H63D) {h63d}")
        risk = True
    if risk:
        guidance = (
            "Endurance athletes with HFE variants have higher iron stores than peers "
            "— DO NOT supplement iron without lab confirmation. Annual ferritin + "
            "transferrin saturation check. Excess iron blunts mitochondrial function "
            "and increases oxidative damage in athletes."
        )
    else:
        guidance = (
            "Endurance athletes commonly run low ferritin (foot-strike haemolysis, "
            "sweat losses). Check ferritin annually if you do >5 h/week endurance. "
            "Target ferritin >40 ng/mL for performance."
        )
    return {"hfe_carrier": risk, "factors": factors or ["HFE not typed"],
            "guidance": guidance}


# ── Stress fracture / bone health ───────────────────────────────────────────

def _analyze_stress_fracture(snps_df) -> dict:
    vdr_fok = _gt(snps_df, "rs2228570")
    col1a1 = _gt(snps_df, "rs1800012")
    factors: list[str] = []
    risk = False
    if vdr_fok:
        factors.append(f"rs2228570 (VDR FokI) {vdr_fok}")
        if "T" in vdr_fok:
            risk = True
    if col1a1:
        factors.append(f"rs1800012 (COL1A1) {col1a1}")
        if "T" in col1a1:
            risk = True
    if risk:
        guidance = (
            "Elevated stress-fracture risk markers — for runners/jumpers: ramp running "
            "mileage ≤10%/week, alternate impact days with cycling/swim, prioritise "
            "vitamin D ≥40 ng/mL, calcium 1000 mg/day, and bone-loading strength work "
            "(heavy squats, jumps) 2×/week — these strengthen bone more than running."
        )
    else:
        guidance = "No stress-fracture risk variants typed — standard load progression."
    return {"elevated_risk": risk, "factors": factors or ["Bone SNPs not typed"],
            "guidance": guidance}


# ── HR zones (formula-based since age not provided) ─────────────────────────

def _hr_zones() -> dict:
    return {
        "formula": "HRmax ≈ 208 − 0.7 × age (Tanaka). Zones as % HRmax:",
        "zones": [
            {"zone": "Z1 Recovery", "pct": "50–60%", "use": "Active recovery, warm-up/cool-down"},
            {"zone": "Z2 Aerobic base", "pct": "60–70%", "use": "Mitochondrial / fat oxidation — most weekly volume"},
            {"zone": "Z3 Tempo", "pct": "70–80%", "use": "Lactate-clearance work, 'comfortably hard'"},
            {"zone": "Z4 Threshold", "pct": "80–90%", "use": "Lactate threshold intervals (8–20 min)"},
            {"zone": "Z5 VO2max", "pct": "90–100%", "use": "3–5 min hard intervals — high-leverage capacity work"},
        ],
        "polarised": "80/20 rule: 80% of weekly training in Z1–Z2, 20% in Z3–Z5. Avoid the 'grey zone' (Z3 every day) — neither builds base nor builds peak.",
    }


# ── Training nutrition (pre/intra/post) ─────────────────────────────────────

def _training_nutrition(power_endurance: dict, caffeine_erg: dict) -> dict:
    is_endurance = power_endurance["bias"].startswith("Endurance") or "endurance-leaning" in power_endurance["bias"]
    return {
        "pre_workout": (
            "60–90 min before: 0.5–1 g/kg carbs + 20 g protein (e.g. oats+whey, "
            "rice cake+jam+yoghurt). " +
            (caffeine_erg.get("guidance", "") if caffeine_erg.get("responder", "").startswith("Strong") else "")
        ),
        "intra_workout": (
            "Sessions <60 min: water only. 60–90 min: 30 g carbs/h (sports drink or gel). "
            "90+ min: 60–90 g carbs/h + 500–750 mL fluid + 300–500 mg sodium."
            if is_endurance else
            "Strength sessions: water only is fine. >90 min lift: 20 g intra-EAAs or carbs optional."
        ),
        "post_workout": (
            "Within 1–2 h: 0.3 g/kg protein (20–40 g) + carbs (1 g/kg if depleted). "
            "Examples: Greek yoghurt + berries + granola; chicken + rice + veg. Hydrate "
            "to 150% of fluid lost (weigh in/out)."
        ),
    }


# ── Periodisation overview ──────────────────────────────────────────────────

def _periodisation(power_endurance: dict) -> list[dict]:
    if power_endurance["bias"].startswith("Power"):
        return [
            {"phase": "Weeks 1–4 — Hypertrophy", "focus": "8–12 reps, 12–16 sets/muscle/wk, RPE 7–8"},
            {"phase": "Weeks 5–8 — Strength", "focus": "3–6 reps, heavy compounds, RPE 7–9"},
            {"phase": "Weeks 9–11 — Power/peak", "focus": "1–3 reps + plyometrics + Olympic lifts"},
            {"phase": "Week 12 — Deload", "focus": "−40% volume, technique focus"},
        ]
    if power_endurance["bias"].startswith("Endurance"):
        return [
            {"phase": "Weeks 1–6 — Base", "focus": "Zone 2 volume +10%/wk, 1 tempo session"},
            {"phase": "Weeks 7–9 — Build", "focus": "Threshold intervals 2×/wk + long Z2"},
            {"phase": "Weeks 10–11 — Peak", "focus": "VO2max work (4×4 min), race-pace efforts"},
            {"phase": "Week 12 — Taper/deload", "focus": "−50% volume, maintain intensity"},
        ]
    return [
        {"phase": "Weeks 1–4 — General prep", "focus": "Balanced strength + Z2 aerobic base"},
        {"phase": "Weeks 5–8 — Specific build", "focus": "Strength bias OR endurance bias chosen by goal"},
        {"phase": "Weeks 9–11 — Intensity", "focus": "HIIT + heavy strength, sport-specific"},
        {"phase": "Week 12 — Deload", "focus": "Mobility, lighter loads, recover for next block"},
    ]


# ── Warm-up / mobility prescription ─────────────────────────────────────────

def _warmup(injury_risk: dict, power_endurance: dict) -> dict:
    elevated = any("elev" in r["level"].lower() for r in injury_risk["risks"])
    base = [
        "5 min easy cardio (bike/row/jog) — raise core temp",
        "Dynamic mobility: leg swings, hip openers, T-spine rotations — 3 min",
        "Activation: glute bridges 2×10, band pull-aparts 2×15, scap pushups 2×10",
    ]
    if power_endurance["bias"].startswith("Power"):
        base.append("Specific ramp: 3 ascending sets to working weight (50/70/90%)")
    if elevated:
        base.insert(0, "EXTRA collagen + vit C 60 min pre-session (15 g + 50 mg) — tendon synthesis")
        base.append("Eccentric prep: Nordics 1×5 (ham) / decline squat 1×10 (patellar) / heel drops (Achilles)")
    return {
        "duration_min": 15 if elevated else 10,
        "protocol": base,
        "cooldown": "5 min easy cardio + 5 min static stretching (focus on tight areas — hips, calves, T-spine).",
    }


# ── Public API ──────────────────────────────────────────────────────────────

def analyze_exercise(snps_df: pd.DataFrame | None) -> dict:
    if snps_df is None:
        return {"status": "no_data"}

    pe = _analyze_power_endurance(snps_df)
    injury = _analyze_injury_risk(snps_df)
    recovery = _analyze_recovery(snps_df)
    chronotype = _analyze_chronotype(snps_df)
    cognitive = _analyze_cognitive_exercise(snps_df)
    vo2 = _analyze_vo2max_trainability(snps_df)
    strength = _analyze_strength_trainability(snps_df)
    fat_loss = _analyze_fat_loss_response(snps_df)
    pain = _analyze_pain_tolerance(snps_df)
    motivation = _analyze_motivation(snps_df)
    caff_erg = _analyze_caffeine_ergogenic(snps_df)
    sleep = _analyze_sleep(snps_df)
    iron_end = _analyze_iron_endurance(snps_df)
    stress_fx = _analyze_stress_fracture(snps_df)
    weekly = _build_weekly_template(pe, recovery, chronotype)
    hr_zones = _hr_zones()
    nutrition_timing = _training_nutrition(pe, caff_erg)
    periodisation = _periodisation(pe)
    warmup = _warmup(injury, pe)

    result = {
        "status": "ok",
        "power_endurance": pe,
        "injury_risk": injury,
        "recovery": recovery,
        "chronotype": chronotype,
        "cognitive": cognitive,
        "vo2max": vo2,
        "strength_trainability": strength,
        "fat_loss": fat_loss,
        "pain_tolerance": pain,
        "motivation": motivation,
        "caffeine_ergogenic": caff_erg,
        "sleep": sleep,
        "iron_endurance": iron_end,
        "stress_fracture": stress_fx,
        "hr_zones": hr_zones,
        "training_nutrition": nutrition_timing,
        "periodisation": periodisation,
        "warmup": warmup,
        "weekly_template": weekly,
    }

    if analyze_advanced_exercise is not None:
        try:
            advanced = analyze_advanced_exercise(snps_df, result)
            result.update(advanced)
        except Exception as exc:
            result["advanced_error"] = str(exc)

    if analyze_exercise_protocols is not None:
        try:
            protocols = analyze_exercise_protocols(result)
            result["protocols"] = protocols
        except Exception as exc:
            result["protocols_error"] = str(exc)

    return result


# ── HTML rendering ──────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_EX_CSS = """
<style>
.ex-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           color:#222; max-width: 1100px; margin: 24px auto; padding: 0 16px; }
.ex-wrap h1 { font-size: 1.6em; border-bottom: 2px solid #333; padding-bottom: 6px; }
.ex-wrap h2 { font-size: 1.2em; margin-top: 28px; padding-bottom:4px;
              border-bottom: 1px solid #eee; }
.ex-card { background:#fcfcfd; border:1px solid #e2e2e6; border-radius:10px;
           padding:14px 16px; margin:10px 0; }
.ex-bias-bar { width:100%; height:18px; border-radius:9px; overflow:hidden;
               display:flex; margin: 8px 0 14px 0; }
.ex-bias-bar .pow { background:#a32a2a; }
.ex-bias-bar .end { background:#1e6091; }
.ex-bias-bar span { color:white; padding:0 8px; font-size:0.8em; line-height:18px; }
.ex-factors { font-family: Menlo, monospace; font-size:0.85em; color:#555;
              background:#f6f6f7; padding:6px 10px; border-radius:6px; margin-top:6px; }
.ex-risk-elev { color:#a32a2a; font-weight:600; }
.ex-risk-prot { color:#2c7a30; font-weight:600; }
.ex-risk-base { color:#5a6772; }
table.ex { width:100%; border-collapse: collapse; margin-top:10px; }
table.ex th, table.ex td { padding:8px 10px; border-bottom:1px solid #eee; text-align:left; }
table.ex th { background:#f9f9f9; }
.attr-bar { background:#eee; border-radius:6px; height:14px; overflow:hidden; }
.attr-bar > div { height:100%; background:linear-gradient(90deg,#1e6091,#a32a2a); }
.attr-row { display:grid; grid-template-columns: 160px 60px 1fr; gap:8px; align-items:center;
            padding:4px 0; font-size:0.92em; }
.sport-rank { display:flex; justify-content:space-between; padding:6px 0; border-bottom:1px solid #f0f0f0; }
.workout-card { background:#fdfcf7; border:1px solid #e9e3c8; border-radius:10px;
                padding:14px 16px; margin:10px 0; }
.workout-card h3 { margin: 0 0 4px; font-size: 1.05em; }
.region-risk { display:grid; grid-template-columns: 180px 60px 1fr; gap:8px; padding:6px 0;
               border-bottom:1px solid #f0f0f0; align-items:center; }
.region-risk .bar { background:#eee; border-radius:4px; height:10px; overflow:hidden; }
.region-risk .bar > div { height:100%; background:linear-gradient(90deg,#2c7a30,#c08327,#a32a2a); }
.readiness-table td { font-size:0.92em; }
</style>
"""


def _risk_class(level: str) -> str:
    if "elev" in level.lower():
        return "ex-risk-elev"
    if "prot" in level.lower():
        return "ex-risk-prot"
    return "ex-risk-base"


def _render_advanced_exercise(result: dict) -> str:
    profile = result.get("composite_profile")
    injury_map = result.get("injury_risk_map")
    readiness = result.get("daily_readiness")
    workouts = result.get("sample_workouts")
    concurrent = result.get("concurrent_training")
    aerobic = result.get("aerobic_estimates")
    taper = result.get("tapering")
    deload = result.get("deload")
    mental = result.get("mental_skills")
    thermal = result.get("thermal_adaptation")
    mobility = result.get("mobility")
    plyo = result.get("plyometric_progression")

    if not profile:
        return ""

    # Composite + sport ranking
    attr_rows = "".join(
        f'<div class="attr-row"><div>{_esc(k.replace("_"," ").title())}</div>'
        f'<div><strong>{v:.0f}</strong></div>'
        f'<div><div class="attr-bar"><div style="width:{v:.0f}%"></div></div></div></div>'
        for k, v in profile["attributes"].items()
    )
    sport_rows = "".join(
        f'<div class="sport-rank"><span>{i+1}. {_esc(s["sport"])}</span>'
        f'<span><strong>{s["fit_score"]:.1f}</strong></span></div>'
        for i, s in enumerate(profile["ranked_sports"])
    )
    profile_html = f"""
<h2>Composite Athletic Profile</h2>
<div class="ex-card">
  <p><strong>Overall index: {profile['overall_score']:.1f}/100</strong></p>
  <p>{_esc(profile['summary'])}</p>
  <h3 style="font-size:1em;margin-top:14px">Attribute scores</h3>
  {attr_rows}
  <h3 style="font-size:1em;margin-top:14px">Top sport matches</h3>
  {sport_rows}
</div>
"""

    # Injury map
    injury_html = ""
    if injury_map:
        rows = "".join(
            f'<div class="region-risk">'
            f'<div>{_esc(r["region"])}</div>'
            f'<div><strong>{r["risk_pct"]}%</strong></div>'
            f'<div><div class="bar"><div style="width:{r["risk_pct"]}%"></div></div>'
            f'<small>{_esc(r["interventions"])}</small></div>'
            f'</div>'
            for r in injury_map["regions"]
        )
        injury_html = f"""
<h2>Body-Region Injury Risk Map</h2>
<div class="ex-card">
  <p>Composite index: <strong>{injury_map['overall_index']}</strong> · per-region quantitative scoring + mitigations</p>
  {rows}
</div>
"""

    # Readiness
    readiness_html = ""
    if readiness:
        rows = "".join(
            f'<tr><td><strong>{_esc(t["range"])}</strong></td><td>{_esc(t["action"])}</td></tr>'
            for t in readiness["thresholds"]
        )
        readiness_html = f"""
<h2>Daily Readiness Formula</h2>
<div class="ex-card">
  <p><em>{_esc(readiness["formula"])}</em></p>
  <table class="ex readiness-table">
    <tr><th>Score</th><th>Action</th></tr>
    {rows}
  </table>
  <p>{_esc(readiness.get("tracking",""))}</p>
  <p><strong>For you:</strong> {_esc(readiness.get("personal_note",""))}</p>
</div>
"""

    # Sample workouts
    workouts_html = ""
    if workouts:
        cards = ""
        for w in workouts:
            blocks = "".join(f"<li>{_esc(b)}</li>" for b in w["blocks"])
            cards += f"""
<div class="workout-card">
  <h3>{_esc(w["name"])}</h3>
  <div style="color:#666;font-size:0.88em">{_esc(w["duration"])} · best at {_esc(w["best_window"])}</div>
  <ul>{blocks}</ul>
  <p><strong>Target intensity:</strong> {_esc(w["rpe_target"])}</p>
</div>"""
        workouts_html = f"<h2>Fully Spelled-Out Sample Workouts</h2>{cards}"

    # Concurrent
    concur_html = ""
    if concurrent:
        concur_html = f"""
<h2>Concurrent Training Interference</h2>
<div class="ex-card">
  <div><strong>Interference level: {_esc(concurrent['interference_level'])}</strong></div>
  <p>{_esc(concurrent['rule'])}</p>
</div>
"""

    # Aerobic estimates
    aero_html = ""
    if aerobic:
        tests = "".join(f"<li>{_esc(t)}</li>" for t in aerobic.get("test_protocols", []))
        aero_html = f"""
<h2>Aerobic Capacity Projections</h2>
<div class="ex-card">
  <p><strong>VO2max:</strong> {_esc(aerobic.get('vo2max_estimate',''))}</p>
  <p><strong>Lactate threshold:</strong> {_esc(aerobic.get('lactate_threshold_estimate',''))}</p>
  <strong>Self-test protocols:</strong><ul>{tests}</ul>
</div>
"""

    # Taper
    taper_html = ""
    if taper:
        s = "".join(f"<li>{_esc(x)}</li>" for x in taper.get("structure", []))
        taper_html = f"""
<h2>Competition Taper Protocol</h2>
<div class="ex-card">
  <p><strong>{taper['duration_days']} days · {_esc(taper['volume_change'])}</strong></p>
  <ol>{s}</ol>
</div>
"""

    # Deload
    deload_html = ""
    if deload:
        s = "".join(f"<li>{_esc(x)}</li>" for x in deload.get("structure", []))
        deload_html = f"""
<h2>Deload Protocol</h2>
<div class="ex-card">
  <p><strong>Frequency:</strong> {_esc(deload['frequency'])}</p>
  <ul>{s}</ul>
  <p><strong>Unscheduled-deload triggers:</strong> {_esc(deload.get('trigger_signals',''))}</p>
</div>
"""

    # Mental skills
    mental_html = ""
    if mental:
        f = "".join(f'<div class="ex-factors">{_esc(x)}</div>' for x in mental.get("factors", []))
        s = "".join(f"<li>{_esc(x)}</li>" for x in mental.get("strategies", []))
        mental_html = f"""
<h2>Mental Skills Profile</h2>
<div class="ex-card">
  <div><strong>{_esc(mental.get('profile',''))}</strong></div>
  {f}
  <ul>{s}</ul>
</div>
"""

    # Thermal
    thermal_html = ""
    if thermal:
        f = "".join(f'<div class="ex-factors">{_esc(x)}</div>' for x in thermal.get("factors", []))
        thermal_html = f"""
<h2>Cold &amp; Heat Adaptation Protocols</h2>
<div class="ex-card">
  <div><strong>{_esc(thermal.get('cold_tolerance',''))}</strong></div>
  {f}
  <p><strong>Cold protocol:</strong> {_esc(thermal.get('cold_protocol',''))}</p>
  <p><strong>Heat (sauna) protocol:</strong> {_esc(thermal.get('heat_protocol',''))}</p>
</div>
"""

    # Mobility
    mobility_html = ""
    if mobility:
        items = "".join(f"<li>{_esc(x)}</li>" for x in mobility.get("items", []))
        mobility_html = f"""
<h2>Mobility / Flexibility Prescription</h2>
<div class="ex-card">
  <p><strong>{mobility.get('daily_minutes','—')} min daily</strong></p>
  <ul>{items}</ul>
  <p><em>{_esc(mobility.get('guidance',''))}</em></p>
</div>
"""

    # Plyometric
    plyo_html = ""
    if plyo:
        if not plyo.get("ready", True):
            phases_str = "".join(f"<li>{_esc(p)}</li>" for p in plyo.get("phases", []))
            plyo_html = f"""
<h2>Plyometric Progression (Preparation Required)</h2>
<div class="ex-card">
  <p><strong>Not yet ready for full plyometrics</strong> — elevated joint/tendon risk.</p>
  <p>{_esc(plyo.get('preparation_weeks',''))}</p>
  <ol>{phases_str}</ol>
</div>
"""
        else:
            rows = "".join(
                f'<tr><td>{_esc(p["weeks"])}</td><td>{_esc(p["level"])}</td><td>{_esc(p["drills"])}</td></tr>'
                for p in plyo.get("phases", [])
            )
            rules = "".join(f"<li>{_esc(r)}</li>" for r in plyo.get("rules", []))
            plyo_html = f"""
<h2>12-Week Plyometric Progression</h2>
<div class="ex-card">
  <table class="ex"><tr><th>Weeks</th><th>Level</th><th>Drills</th></tr>{rows}</table>
  <ul>{rules}</ul>
</div>
"""

    blob = (profile_html + injury_html + readiness_html + workouts_html
            + concur_html + aero_html + taper_html + deload_html
            + mental_html + thermal_html + mobility_html + plyo_html)
    blob += _render_exercise_protocols(result)
    return blob


def _render_exercise_protocols(result: dict) -> str:
    p = result.get("protocols")
    if not p:
        return ""
    out = ""

    # Fitness tests
    ft = p.get("fitness_tests")
    if ft:
        cards = ""
        for name, test in ft.items():
            cards += f"<div class='ex-card'><strong>{_esc(name.replace('_',' ').title())}</strong>"
            for k, v in test.items():
                if isinstance(v, dict):
                    sub = "".join(f"<li>{_esc(kk)}: {_esc(vv)}</li>" for kk, vv in v.items())
                    cards += f"<p><em>{_esc(k.replace('_',' '))}:</em><ul>{sub}</ul></p>"
                else:
                    cards += f"<p><em>{_esc(k.replace('_',' '))}:</em> {_esc(v)}</p>"
            cards += "</div>"
        out += f"<h2>Fitness-Test Calculators</h2>{cards}"

    # HR zones for age 35 demo
    hz = p.get("hr_zone_demo_age_35")
    if hz:
        rows = "".join(
            f"<tr><td>{_esc(z['zone'])}</td><td>{_esc(z['pct_hrmax'])}</td>"
            f"<td>{_esc(z['bpm_range_karvonen'])} bpm</td></tr>" for z in hz["zones"]
        )
        out += f"""
<h2>HR Zones (example: age {hz['age']}, RHR {hz['resting_hr']})</h2>
<div class="ex-card">
<p>HRmax (Tanaka): <strong>{hz['HRmax_tanaka']} bpm</strong> · HR reserve: {hz['HR_reserve']}</p>
<table class="ex"><tr><th>Zone</th><th>% HRmax</th><th>BPM (Karvonen)</th></tr>{rows}</table>
<p><em>{_esc(hz['note'])}</em></p>
</div>
"""

    # Running paces
    rp = p.get("running_pace_demo_threshold_5min_km")
    if rp:
        rows = "".join(
            f"<tr><td>{_esc(z['zone'])}</td><td>{_esc(z['pace_min_km'])} min/km</td></tr>"
            for z in rp["zones"]
        )
        out += f"""
<h2>Running Pace Zones (example: 5:00 min/km threshold)</h2>
<div class="ex-card">
<table class="ex"><tr><th>Zone</th><th>Pace</th></tr>{rows}</table>
</div>
"""

    # FTP
    ftp = p.get("cycling_ftp_demo")
    if ftp:
        rows = "".join(
            f"<tr><td>{_esc(z['zone'])}</td><td>{_esc(z['watts'])}</td></tr>"
            for z in ftp["zones"]
        )
        tests = "".join(f"<li>{_esc(t)}</li>" for t in ftp["ftp_test_protocols"])
        out += f"""
<h2>Cycling FTP Zones (example: FTP {ftp['ftp_watts']} W)</h2>
<div class="ex-card">
<table class="ex"><tr><th>Zone</th><th>Watts</th></tr>{rows}</table>
<strong>FTP test protocols</strong><ul>{tests}</ul>
</div>
"""

    # Race time predictor
    rt = p.get("race_time_demo")
    if rt:
        out += f"""
<h2>Race-Time Predictor (Riegel)</h2>
<div class="ex-card">
<p>Known: {_esc(rt['known'])}</p>
<p><strong>Predicted:</strong> {_esc(rt['predicted_target'])}</p>
<p><em>{_esc(rt['caveat'])}</em></p>
</div>
"""

    # Sport-specific 12-wk plans
    plans = p.get("sport_specific_plans")
    if plans:
        cards = ""

        def _week_body(w: dict) -> str:
            # Week dicts are heterogeneous across sports: strength plans use
            # `schedule`(+`accessory`); endurance plans (marathon / cycling /
            # triathlon) use `weekly_mileage_km`(+`key_workouts`). Render every
            # descriptive field that is present rather than a single fallback,
            # so no plan's weeks render blank.
            parts = []
            if w.get("weekly_mileage_km"):
                parts.append(
                    f"<em>Weekly volume (km):</em> {_esc(w['weekly_mileage_km'])}"
                )
            if w.get("schedule"):
                parts.append(_esc(w["schedule"]))
            if w.get("key_workouts"):
                parts.append(f"<em>Key workouts:</em> {_esc(w['key_workouts'])}")
            if w.get("accessory"):
                parts.append(f"<em>Accessory:</em> {_esc(w['accessory'])}")
            if not parts:  # legacy fallback
                body = w.get("weekly_template", w.get("focus", ""))
                if body:
                    parts.append(_esc(body))
            return "<br>".join(parts)

        for pl in plans:
            weeks = "".join(
                f"<li><strong>{_esc(w.get('phase',''))}:</strong> {_week_body(w)}</li>"
                for w in pl["weeks"]
            )
            extra = ""
            if "meet_day" in pl:
                extra = f"<p><strong>Meet day:</strong> {_esc(pl['meet_day'])}</p>"
            if "race_day_fueling" in pl:
                extra = f"<p><strong>Race-day fueling:</strong> {_esc(pl['race_day_fueling'])}</p>"
            cards += (f"<div class='ex-card'><strong>{_esc(pl['sport'])} — 12-week plan</strong>"
                      f"<ul>{weeks}</ul>{extra}</div>")
        out += f"<h2>Sport-Specific 12-Week Plans (top matches)</h2>{cards}"

    # Lifting cues
    lc = p.get("lifting_cues")
    if lc:
        cards = ""
        for card in lc:
            cues = "".join(f"<li>{_esc(c)}</li>" for c in card["cues"])
            faults = "".join(f"<li>{_esc(c)}</li>" for c in card["common_faults"])
            cards += (f"<div class='ex-card'><strong>{_esc(card['lift'])}</strong>"
                      f"<p><em>Setup:</em> {_esc(card['setup'])}</p>"
                      f"<strong>Cues</strong><ul>{cues}</ul>"
                      f"<strong>Common faults</strong><ul>{faults}</ul></div>")
        out += f"<h2>Lifting Cue Library</h2>{cards}"

    # Recovery modality matrix
    rm = p.get("recovery_matrix")
    if rm:
        rows = "".join(
            f"<tr><td>{_esc(m['modality'])}</td><td>{_esc(m['effect_size'])}</td>"
            f"<td>{_esc(m['evidence'])}</td><td>{_esc(m['protocol'])}</td></tr>"
            for m in rm
        )
        out += f"""
<h2>Recovery Modality Matrix</h2>
<div class="ex-card">
<table class="ex"><tr><th>Modality</th><th>Effect</th><th>Evidence</th><th>Protocol</th></tr>{rows}</table>
</div>
"""

    # Master athlete
    ma = p.get("master_athlete")
    if ma:
        cards = ""
        for d in ma:
            adj = "".join(f"<li>{_esc(x)}</li>" for x in d["adjustments"])
            cards += (f"<div class='ex-card'><strong>{_esc(d['decade'])}</strong>"
                      f"<p><em>{_esc(d['key_changes'])}</em></p><ul>{adj}</ul></div>")
        out += f"<h2>Master Athlete Adjustments by Decade</h2>{cards}"

    # Supplement stacks
    ss = p.get("supplement_stacks")
    if ss:
        pre = "".join(f"<li>{_esc(x)}</li>" for x in ss["pre_workout"])
        intra = "".join(f"<li>{_esc(x)}</li>" for x in ss["intra_workout"])
        post = "".join(f"<li>{_esc(x)}</li>" for x in ss["post_workout"])
        ev = ss["evidence_tiers"]
        a = "".join(f"<li>{_esc(x)}</li>" for x in ev["A_proven"])
        b = "".join(f"<li>{_esc(x)}</li>" for x in ev["B_some_evidence"])
        c = "".join(f"<li>{_esc(x)}</li>" for x in ev["C_avoid_or_unproven"])
        out += f"""
<h2>Workout Supplement Stacks</h2>
<div class="ex-card">
  <strong>Pre-workout</strong><ul>{pre}</ul>
  <strong>Intra-workout</strong><ul>{intra}</ul>
  <strong>Post-workout</strong><ul>{post}</ul>
  <strong>Tier A — proven</strong><ul>{a}</ul>
  <strong>Tier B — some evidence</strong><ul>{b}</ul>
  <strong>Tier C — avoid / unproven</strong><ul>{c}</ul>
  <p><em>{_esc(ss['personalised_priority'])}</em></p>
</div>
"""

    # Movement screen
    ms = p.get("movement_screen")
    if ms:
        rows = "".join(
            f"<tr><td>{_esc(t['test'])}</td><td>{_esc(t['checks'])}</td>"
            f"<td>{_esc(t['common_fail'])}</td><td>{_esc(t['remedy'])}</td></tr>"
            for t in ms
        )
        out += f"""
<h2>Movement Screen (self-test)</h2>
<div class="ex-card">
<table class="ex"><tr><th>Test</th><th>Pass criteria</th><th>Common fail</th><th>Remedy</th></tr>{rows}</table>
</div>
"""

    # Detraining timeline
    dt = p.get("detraining_timeline")
    if dt:
        rows = "".join(
            f"<tr><td>{_esc(d['period'])}</td><td>{_esc(d['what_happens'])}</td></tr>"
            for d in dt
        )
        out += f"""
<h2>Detraining / Retraining Timeline</h2>
<div class="ex-card">
<table class="ex"><tr><th>Period</th><th>Physiology</th></tr>{rows}</table>
</div>
"""

    # HRV-guided week
    hg = p.get("hrv_guided_week")
    if hg:
        rules = "".join(
            f"<tr><td>{_esc(r['signal'])}</td><td>{_esc(r['action'])}</td></tr>"
            for r in hg["decision_rules"]
        )
        ctx = "".join(f"<li>{_esc(x)}</li>" for x in hg["context_modifiers"])
        out += f"""
<h2>HRV-Guided Week Structure</h2>
<div class="ex-card">
<p>{_esc(hg['tracking'])}</p>
<table class="ex"><tr><th>Signal</th><th>Action</th></tr>{rules}</table>
<strong>Context modifiers</strong><ul>{ctx}</ul>
<p><em>{_esc(hg['implementation'])}</em></p>
</div>
"""

    return out


def render_exercise_html(result: dict, file_label: str = "") -> str:
    if not result or result.get("status") != "ok":
        return f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Exercise</title>{_EX_CSS}</head><body>
<div class="ex-wrap"><h1>Personalised Exercise Programming</h1>
<p>Insufficient genetic data for exercise programming.</p></div></body></html>"""

    pe = result["power_endurance"]
    pow_pct = pe["ratio_pct_power"]
    end_pct = pe["ratio_pct_endurance"]
    bias_bar = (
        f'<div class="ex-bias-bar">'
        f'<span class="pow" style="width:{pow_pct}%">Power {pow_pct}%</span>'
        f'<span class="end" style="width:{end_pct}%">Endurance {end_pct}%</span>'
        f'</div>'
    )

    risk_rows = "".join(
        f'<tr><td>{_esc(r["tissue"])}</td>'
        f'<td class="{_risk_class(r["level"])}">{_esc(r["level"])}</td>'
        f'<td>{_esc(r["marker"])}</td>'
        f'<td>{_esc(r["mitigation"])}</td></tr>'
        for r in result["injury_risk"]["risks"]
    )

    week_rows = "".join(
        f'<tr><td>{_esc(d["day"])}</td><td>{_esc(d["session"])}</td>'
        f'<td>{_esc(d.get("time",""))}</td></tr>'
        for d in result["weekly_template"]
    )

    factors_html = "".join(
        f'<div class="ex-factors">{_esc(f)}</div>' for f in pe["factors"]
    ) or '<div class="ex-factors">No bias SNPs typed</div>'

    cog = result["cognitive"]

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Personalised Exercise{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_EX_CSS}</head><body><div class="ex-wrap">
<h1>Personalised Exercise Programming</h1>

<h2>Power vs Endurance Bias</h2>
<div class="ex-card">
  <div><strong>{_esc(pe['bias'])}</strong></div>
  {bias_bar}
  {factors_html}
  <p>{_esc(pe['recommendation'])}</p>
</div>

<h2>Injury Risk</h2>
<div class="ex-card">
<table class="ex">
  <tr><th>Tissue</th><th>Risk</th><th>Marker</th><th>Mitigation</th></tr>
  {risk_rows}
</table>
</div>

<h2>Recovery Speed</h2>
<div class="ex-card">
  <div><strong>{_esc(result['recovery']['speed'])}</strong> recovery profile
       (inflammation score {result['recovery']['inflammation_score']})</div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result['recovery']['factors'])}
  <p>{_esc(result['recovery']['protocol'])}</p>
</div>

<h2>Optimal Training Window</h2>
<div class="ex-card">
  <div><strong>{_esc(result['chronotype']['chronotype'])} chronotype</strong> —
       suggested window <strong>{_esc(result['chronotype']['optimal_window'])}</strong></div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result['chronotype']['factors'])}
  <p>{_esc(result['chronotype']['rationale'])}</p>
</div>

<h2>Cognitive Benefit Profile (BDNF)</h2>
<div class="ex-card">
  <div><strong>{_esc(cog.get('genotype') or 'Not tested')}</strong></div>
  <p>{_esc(cog['summary'])}</p>
  <p><em>{_esc(cog['training_note'])}</em></p>
</div>

<h2>VO2max Trainability</h2>
<div class="ex-card">
  <div><strong>{_esc(result.get('vo2max',{}).get('tier','—'))}</strong>
       (score {result.get('vo2max',{}).get('score','—')})</div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('vo2max',{}).get('factors',[]))}
  <p>{_esc(result.get('vo2max',{}).get('guidance',''))}</p>
</div>

<h2>Strength / Hypertrophy Trainability</h2>
<div class="ex-card">
  <div><strong>{_esc(result.get('strength_trainability',{}).get('tier','—'))}</strong></div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('strength_trainability',{}).get('factors',[]))}
  <p>{_esc(result.get('strength_trainability',{}).get('recommendation',''))}</p>
</div>

<h2>Fat-Loss Response to Exercise</h2>
<div class="ex-card">
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('fat_loss',{}).get('factors',[]))}
  <p>{_esc(result.get('fat_loss',{}).get('guidance',''))}</p>
</div>

<h2>Pain Tolerance & RPE Calibration</h2>
<div class="ex-card">
  <div><strong>{_esc(result.get('pain_tolerance',{}).get('tolerance','—'))}</strong></div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('pain_tolerance',{}).get('factors',[]))}
  <p>{_esc(result.get('pain_tolerance',{}).get('guidance',''))}</p>
</div>

<h2>Motivation & Habit Formation (Dopamine)</h2>
<div class="ex-card">
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('motivation',{}).get('factors',[]))}
  <p>{_esc(result.get('motivation',{}).get('guidance',''))}</p>
</div>

<h2>Caffeine as Ergogenic Aid</h2>
<div class="ex-card">
  <div><strong>{_esc(result.get('caffeine_ergogenic',{}).get('responder','—'))}</strong></div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('caffeine_ergogenic',{}).get('factors',[]))}
  <p>{_esc(result.get('caffeine_ergogenic',{}).get('guidance',''))}</p>
</div>

<h2>Sleep & Recovery</h2>
<div class="ex-card">
  <div><strong>{_esc(result.get('sleep',{}).get('sleep_phenotype','—'))}</strong></div>
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('sleep',{}).get('factors',[]))}
  <p>{_esc(result.get('sleep',{}).get('guidance',''))}</p>
</div>

<h2>Iron & Endurance Athletes</h2>
<div class="ex-card">
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('iron_endurance',{}).get('factors',[]))}
  <p>{_esc(result.get('iron_endurance',{}).get('guidance',''))}</p>
</div>

<h2>Stress-Fracture / Bone Health</h2>
<div class="ex-card">
  {"".join(f'<div class="ex-factors">{_esc(f)}</div>' for f in result.get('stress_fracture',{}).get('factors',[]))}
  <p>{_esc(result.get('stress_fracture',{}).get('guidance',''))}</p>
</div>

<h2>Heart-Rate Zones</h2>
<div class="ex-card">
  <p><em>{_esc(result.get('hr_zones',{}).get('formula',''))}</em></p>
  <table class="ex">
    <tr><th>Zone</th><th>% HRmax</th><th>Use</th></tr>
    {"".join(f'<tr><td>{_esc(z["zone"])}</td><td>{_esc(z["pct"])}</td><td>{_esc(z["use"])}</td></tr>' for z in result.get('hr_zones',{}).get('zones',[]))}
  </table>
  <p>{_esc(result.get('hr_zones',{}).get('polarised',''))}</p>
</div>

<h2>Warm-Up Prescription ({result.get('warmup',{}).get('duration_min','—')} min)</h2>
<div class="ex-card">
  <ul>{"".join(f"<li>{_esc(x)}</li>" for x in result.get('warmup',{}).get('protocol',[]))}</ul>
  <p><strong>Cool-down:</strong> {_esc(result.get('warmup',{}).get('cooldown',''))}</p>
</div>

<h2>Training Nutrition (Pre / Intra / Post)</h2>
<div class="ex-card">
  <p><strong>Pre-workout:</strong> {_esc(result.get('training_nutrition',{}).get('pre_workout',''))}</p>
  <p><strong>Intra-workout:</strong> {_esc(result.get('training_nutrition',{}).get('intra_workout',''))}</p>
  <p><strong>Post-workout:</strong> {_esc(result.get('training_nutrition',{}).get('post_workout',''))}</p>
</div>

<h2>12-Week Periodisation</h2>
<div class="ex-card">
  <table class="ex">
    <tr><th>Phase</th><th>Focus</th></tr>
    {"".join(f'<tr><td>{_esc(p["phase"])}</td><td>{_esc(p["focus"])}</td></tr>' for p in result.get('periodisation',[]))}
  </table>
</div>

<h2>Example 7-Day Template</h2>
<div class="ex-card">
<table class="ex">
  <tr><th>Day</th><th>Session</th><th>Time</th></tr>
  {week_rows}
</table>
</div>

{_render_advanced_exercise(result)}

<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. Recommendations are starting points derived from published
sport-genetics literature. Adjust based on training history, injury status,
and lab-confirmed inflammation/recovery markers. Consult a physician before
starting a new exercise programme if you have a cardiovascular or
musculoskeletal condition.
</p>
</div></body></html>"""
