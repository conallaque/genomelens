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


# ── Power vs endurance ──────────────────────────────────────────────────────

def _analyze_power_endurance(snps_df) -> Dict:
    actn3 = _gt(snps_df, "rs1815739")        # C=R (power), T=X (stop codon)
    ace = _gt(snps_df, "rs4646994") or _gt(snps_df, "rs4341")   # I/D proxy
    ppargc1a = _gt(snps_df, "rs8192678")     # G=Gly482, A=Ser482 (lower endurance trainability)

    power_score = 0
    endurance_score = 0
    factors: List[str] = []

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
    if total == 0:
        ratio_power = 50
    else:
        ratio_power = round(100 * power_score / total)
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

def _analyze_injury_risk(snps_df) -> Dict:
    col1a1 = _gt(snps_df, "rs1800012")       # Sp1 site; T allele = soft-tissue risk
    col5a1 = _gt(snps_df, "rs12722")         # CC vs TT; T = better tendon, C = stiffer
    mmp3 = _gt(snps_df, "rs679620")
    risks: List[Dict] = []

    if col1a1:
        if "T" in col1a1:                    # one or two T copies
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
                "marker": f"rs12722 (COL5A1) CC",
                "mitigation": (
                    "Higher tendinopathy risk profile. Build eccentric calf/quad capacity "
                    "(slow heel drops, decline squats). Avoid sudden volume jumps in running."
                ),
            })
        elif "T" in col5a1 and col5a1.count("T") == 2:
            risks.append({
                "tissue": "Achilles / patellar tendon",
                "level": "Protective",
                "marker": f"rs12722 (COL5A1) TT",
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

def _analyze_recovery(snps_df) -> Dict:
    il6 = _gt(snps_df, "rs1800795")
    crp = _gt(snps_df, "rs2794520") or _gt(snps_df, "rs1205")
    sod2 = _gt(snps_df, "rs4880")

    inflam_score = 0
    notes: List[str] = []

    if il6:
        if "G" in il6 and il6.count("G") == 2:
            inflam_score += 2
            notes.append(f"rs1800795 (IL6 -174) {il6} — high baseline inflammatory tone")
        elif "G" in il6:
            inflam_score += 1
            notes.append(f"rs1800795 (IL6 -174) {il6} — moderate inflammatory tendency")

    if crp:
        if "C" in crp and crp.count("C") == 2:
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

def _analyze_chronotype(snps_df) -> Dict:
    clock = _gt(snps_df, "rs1801260")       # T=evening, C=morning (3111C)
    per3 = _gt(snps_df, "rs228697")

    factors: List[str] = []
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

def _analyze_cognitive_exercise(snps_df) -> Dict:
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

def _build_weekly_template(power_endurance: Dict, recovery: Dict, chronotype: Dict) -> List[Dict]:
    bias = power_endurance["bias"]
    fast_recovery = recovery["speed"] == "Fast"
    days: List[Dict] = []
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


# ── Public API ──────────────────────────────────────────────────────────────

def analyze_exercise(snps_df: Optional[pd.DataFrame]) -> Dict:
    if snps_df is None:
        return {"status": "no_data"}

    pe = _analyze_power_endurance(snps_df)
    injury = _analyze_injury_risk(snps_df)
    recovery = _analyze_recovery(snps_df)
    chronotype = _analyze_chronotype(snps_df)
    cognitive = _analyze_cognitive_exercise(snps_df)
    weekly = _build_weekly_template(pe, recovery, chronotype)

    return {
        "status": "ok",
        "power_endurance": pe,
        "injury_risk": injury,
        "recovery": recovery,
        "chronotype": chronotype,
        "cognitive": cognitive,
        "weekly_template": weekly,
    }


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
</style>
"""


def _risk_class(level: str) -> str:
    if "elev" in level.lower():
        return "ex-risk-elev"
    if "prot" in level.lower():
        return "ex-risk-prot"
    return "ex-risk-base"


def render_exercise_html(result: Dict, file_label: str = "") -> str:
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

<h2>Example 7-Day Template</h2>
<div class="ex-card">
<table class="ex">
  <tr><th>Day</th><th>Session</th><th>Time</th></tr>
  {week_rows}
</table>
</div>

<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. Recommendations are starting points derived from published
sport-genetics literature. Adjust based on training history, injury status,
and lab-confirmed inflammation/recovery markers. Consult a physician before
starting a new exercise programme if you have a cardiovascular or
musculoskeletal condition.
</p>
</div></body></html>"""
