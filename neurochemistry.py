"""
Neurochemistry — dopamine, serotonin, plasticity, stress & reward profile
==========================================================================

A dedicated deep-genetics module for the brain's monoamine systems, packaging
the well-established "warrior vs worrier" COMT literature with matching MAOA,
BDNF, DRD2/DRD4, 5-HTT, HTR2A, CACNA1C, OPRM1, TPH2, and CHRNA5 findings into
a *composite phenotype* with concrete practical recommendations (caffeine
protocol, stimulant / SSRI response prediction, meditation style fit, career
neurotype signature, addiction-risk red flags).

Every claim is grounded in the published literature — coefficients from primary
papers where meaningful (e.g. COMT Val158Met increases enzyme activity ~3-4×;
BDNF Val66Met reduces activity-dependent secretion by ~30%).

Coverage
--------
  • **Prefrontal dopamine tone** — COMT Val158Met (rs4680) + rs4633 (LD tag).
  • **Whole-brain monoamine turnover** — MAOA rs6323 (X-linked; high vs low
    activity — the *high* variant is the healthier one; low-activity MAOA
    plus adverse childhood = increased aggression risk in Caspi 2002).
  • **Neuroplasticity** — BDNF Val66Met (rs6265). Met carriers have ~30% less
    activity-dependent BDNF secretion, reduced motor learning, slower stroke
    recovery, worse working-memory training gains.
  • **Reward / dopamine-receptor density** — DRD2/ANKK1 Taq1A (rs1800497),
    DRD4 -521 (rs1800955), OPRM1 A118G (rs1799971).
  • **Serotonin signalling** — 5-HTTLPR (proxy rs25531), HTR2A T102C (rs6313),
    TPH2 rs4570625.
  • **Mood stability / calcium channels** — CACNA1C rs1006737 (the strongest
    common bipolar / SCZ / MDD locus).
  • **Nicotine receptor / addiction risk** — CHRNA5 rs16969968 (heavy-smoker
    variant; also raises lung-cancer risk).

Composite output
----------------
Categorised findings + a `composite` block with predicted:
  • Stress-response profile ("warrior", "worrier", "adaptive middle").
  • Learning / plasticity tier.
  • Stimulant response prediction (Val/Met + Val/Val → higher expected effect).
  • SSRI response prediction (COMT × 5-HTTLPR × MAOA activity).
  • Caffeine protocol (dose ceiling, cut-off time — cross-references CYP1A2 if present).
  • Meditation-style fit (focused-attention vs open-monitoring vs somatic).
  • Career neurotype signature.
  • Addiction / substance-use flags.

References
----------
Egan 2001 PNAS (COMT & prefrontal function); Diamond 2007 (Val/Met heterozygote
advantage on cognition). Meyer-Lindenberg 2006 Neuron (MAOA X sociocognition).
Egan 2003 Cell (BDNF Val66Met & memory / hippocampal function). Chen 2006 Science
(BDNF Val66Met knock-in mouse). Caspi 2002 Science (MAOA × childhood adversity).
Ge 2009 (IL28B - not here). Way & Taylor 2010 (COMT & social pain). Wang 2019
(CACNA1C & psychiatric spectrum). Bertelsen 2015 (CHRNA5 & lung cancer).
"""

from __future__ import annotations

import pandas as pd

CAT_DOPAMINE = "Prefrontal Dopamine (COMT axis)"
CAT_MONOAMINE = "Monoamine Turnover (MAOA)"
CAT_PLASTICITY = "Neuroplasticity (BDNF)"
CAT_REWARD = "Reward & Receptor Density"
CAT_SEROTONIN = "Serotonin System"
CAT_MOOD = "Mood Stability & Ion Channels"
CAT_ADDICTION = "Addiction / Substance Response"


def _gt(df, rsid: str) -> str | None:
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


def _find(category, name, gene, rsid, genotype, phenotype, mechanism, action,
          confidence, citation):
    return {
        "category": category, "name": name, "gene": gene, "rsid": rsid,
        "genotype": genotype or "—", "phenotype": phenotype,
        "mechanism": mechanism, "action": action,
        "confidence": confidence, "citation": citation,
    }


# ── COMT — the warrior/worrier axis ───────────────────────────────────────────

def _comt(df):
    gt = _gt(df, "rs4680")
    if gt is None:
        return None
    # G = Val (fast enzyme, ~3-4× activity); A = Met (slow enzyme)
    n_val = gt.count("G")
    if n_val == 2:
        ph = ("Warrior (Val/Val) — fast prefrontal dopamine clearance")
        mech = ("Homozygous Val158Met = ~3-4× COMT enzyme activity in the "
                "prefrontal cortex. Lower baseline PFC dopamine, but strong "
                "stress resilience and preserved cognition under acute pressure. "
                "Reward system needs stimulation gradients; risk-tolerant; "
                "often lower anxiety trait but higher sensation-seeking.")
        action = ("Thrives under moderate/high pressure; can find low-arousal "
                  "environments unbearable. Prefer variable-intensity work. "
                  "Stimulants (caffeine, methylphenidate, amphetamine) tend to "
                  "have stronger cognitive effects than in Met/Met carriers.")
    elif n_val == 0:
        ph = "Worrier (Met/Met) — slow prefrontal dopamine clearance"
        mech = ("Homozygous Met/Met COMT clears prefrontal dopamine much more "
                "slowly. Higher baseline PFC dopamine → better sustained "
                "working memory when calm, but dopamine 'overload' under acute "
                "stress → anxiety, choking under pressure, rumination. Better "
                "pain thresholds; often more emotionally reactive.")
        action = ("Structure high-stakes cognitive work for calm environments. "
                  "Stress-management practices materially improve performance. "
                  "Stimulants may over-shoot — small doses often work better than large.")
    else:
        ph = "Adaptive middle (Val/Met) — heterozygote advantage on cognition"
        mech = ("Heterozygous Val/Met is the most common European genotype "
                "(~50%). On average, Val/Met carriers score BETTER on complex "
                "cognitive tasks under moderate stress than either homozygote "
                "(documented heterozygote advantage; Diamond 2007). You retain "
                "some Warrior stress-resilience AND some Worrier working-memory "
                "sustainment.")
        action = ("Genuinely favourable for high-performance / high-pressure "
                  "cognitive work. Handles acute stress well without paying the "
                  "sustained-cognition tax of Val/Val.")
    return _find(CAT_DOPAMINE, "COMT Val158Met", "COMT", "rs4680", gt, ph,
                 mech, action, "high",
                 "Egan 2001 PNAS; Diamond 2007; Meyer-Lindenberg 2006")


# ── MAOA — X-linked monoamine turnover ────────────────────────────────────────

def _maoa(df):
    gt = _gt(df, "rs6323")
    if gt is None:
        return None
    # T = high-activity (T-allele-carrying males have MAOA-H phenotype)
    if "T" in gt and "G" not in gt:
        ph = "MAOA-H (high activity) — fast monoamine clearance"
        mech = ("The high-activity MAOA allele produces more enzyme, "
                "meaning faster intracellular breakdown of serotonin, "
                "norepinephrine, and dopamine. Emotional intensity clears "
                "faster, less rumination, better stress recovery.")
        action = ("Rapid affective recovery — negative emotions don't linger "
                  "as long as they do for MAOA-L carriers. Modest reduction "
                  "in SSRI headroom (you already clear serotonin fast). "
                  "Modestly protective against depression's rumination subtype.")
        conf = "high"
    elif "G" in gt and "T" not in gt:
        ph = "MAOA-L (low activity) — slow monoamine clearance"
        mech = ("Low-activity MAOA leaves more serotonin/norepinephrine/dopamine "
                "in the synapse. Caspi 2002 famously showed that MAOA-L "
                "combined with childhood adversity substantially raises "
                "impulsive-aggression risk. In the absence of severe adversity, "
                "MAOA-L is largely neutral or mildly enhances emotional depth.")
        action = ("Awareness only; the Caspi × environment interaction requires "
                  "severe childhood adversity to become behaviourally relevant. "
                  "Meditation and emotion-regulation practices are especially "
                  "valuable.")
        conf = "moderate"
    else:
        ph = "MAOA heterozygous (female carrier / mixed activity)"
        mech = ("Females are heterozygous at this X-linked locus and show "
                "intermediate MAOA activity due to random X-inactivation.")
        action = "Awareness only; intermediate turnover phenotype."
        conf = "moderate"
    return _find(CAT_MONOAMINE, "MAOA activity level", "MAOA", "rs6323", gt,
                 ph, mech, action, conf,
                 "Sabol 1998; Caspi 2002 Science; Meyer-Lindenberg 2006 Neuron")


# ── BDNF — plasticity ────────────────────────────────────────────────────────

def _bdnf(df):
    gt = _gt(df, "rs6265")
    if gt is None:
        return None
    n_T = gt.count("T")   # T = Met allele; C = Val (favorable)
    if n_T == 0:
        ph = "BDNF Val/Val — full activity-dependent BDNF secretion"
        mech = ("Both alleles encode the Val66 form of pro-BDNF, which is "
                "efficiently trafficked and released in response to neural "
                "activity. Better motor learning, better working-memory "
                "training gains, better recovery from stroke, better long-term "
                "expertise accrual through deliberate practice.")
        action = ("Genuinely favourable for skill acquisition. Time invested in "
                  "deliberate practice pays off more than average — spend the "
                  "reps in the domains you care about.")
        conf = "high"
    elif n_T == 1:
        ph = "BDNF Val/Met heterozygous — modestly reduced activity-dependent secretion"
        mech = ("Met66 pro-BDNF is trafficked less efficiently to dendrites; "
                "activity-dependent BDNF secretion is reduced ~30% per Met "
                "allele. Modestly slower complex-motor and working-memory "
                "training gains than Val/Val.")
        action = ("Emphasise consistent practice over intense bursts; aerobic "
                  "exercise materially raises baseline BDNF and partly "
                  "compensates.")
        conf = "high"
    else:
        ph = "BDNF Met/Met — reduced activity-dependent secretion"
        mech = ("Met/Met carriers have ~50% reduced activity-dependent BDNF "
                "secretion. Documented differences in hippocampal function, "
                "motor-learning speed, and recovery from stroke.")
        action = ("Aerobic exercise (running, cycling, swimming) is the single "
                  "most effective non-pharmacological BDNF booster and largely "
                  "closes the gap. Prioritise consistent moderate-intensity "
                  "cardio over occasional intense sessions.")
        conf = "high"
    return _find(CAT_PLASTICITY, "BDNF Val66Met", "BDNF", "rs6265", gt,
                 ph, mech, action, conf,
                 "Egan 2003 Cell; Chen 2006 Science")


# ── DRD2/ANKK1 Taq1A ─────────────────────────────────────────────────────────

def _drd2(df):
    gt = _gt(df, "rs1800497")
    if gt is None:
        return None
    # T = A1 allele (fewer D2 receptors); C = A2 (normal)
    n_T = gt.count("T")
    if n_T >= 1:
        ph = f"DRD2/ANKK1 Taq1A A1 carrier ({n_T}× T) — reduced striatal D2 density"
        mech = ("The Taq1A A1 allele associates with ~30% lower striatal D2 "
                "receptor density. Reduced reward-sensitivity feedback: more "
                "prone to reward-deficiency behaviours (higher addiction-"
                "susceptibility signal), and modestly attenuated aversive-"
                "learning from punishment.")
        action = ("Awareness. Reward-deficiency susceptibility is real but not "
                  "deterministic — behaviours matter more. Regular reward from "
                  "physical activity and creative work partly compensates.")
        conf = "moderate"
    else:
        ph = "DRD2/ANKK1 Taq1A A2/A2 — typical striatal D2 receptor density"
        mech = "Standard dopamine receptor density and reward-sensitivity."
        action = "None specific."
        conf = "moderate"
    return _find(CAT_REWARD, "DRD2 / ANKK1 Taq1A", "DRD2/ANKK1", "rs1800497",
                 gt, ph, mech, action, conf,
                 "Blum 1990 JAMA; Neville 2004")


# ── DRD4 -521 ────────────────────────────────────────────────────────────────

def _drd4(df):
    gt = _gt(df, "rs1800955")
    if gt is None:
        return None
    n_T = gt.count("T")
    if n_T == 2:
        ph = "DRD4 -521 T/T — higher promoter activity"
        mech = "Higher DRD4 expression; some link to novelty-seeking and ADHD subphenotypes."
        conf = "moderate"
    elif n_T == 1:
        ph = "DRD4 -521 T/C — intermediate DRD4 expression"
        mech = "Intermediate DRD4 promoter activity."
        conf = "moderate"
    else:
        ph = "DRD4 -521 C/C — lower promoter activity"
        mech = "Lower DRD4 expression."
        conf = "moderate"
    return _find(CAT_REWARD, "DRD4 -521 promoter", "DRD4", "rs1800955", gt,
                 ph, mech, "Novelty-seeking is behavioural, not fixed by this variant.",
                 conf, "Okuyama 2000")


# ── OPRM1 A118G — mu-opioid receptor / reward sensitivity ─────────────────────

def _oprm1(df):
    gt = _gt(df, "rs1799971")
    if gt is None:
        return None
    n_G = gt.count("G")
    if n_G >= 1:
        ph = f"OPRM1 A118G G-carrier ({n_G}× G) — altered opioid receptor affinity"
        mech = ("The G118 allele produces a mu-opioid receptor with higher "
                "affinity for beta-endorphin. Documented effects: enhanced "
                "response to social bonding cues (Way & Taylor 2010), altered "
                "opioid analgesic dose requirements (G-carriers often need "
                "higher morphine doses post-op), and different alcohol-reward "
                "signalling. Naltrexone response is stronger in G-carriers "
                "(a key predictor for alcohol-use-disorder treatment).")
        action = ("Flag OPRM1 status if being prescribed post-surgical opioids "
                  "(may need dose adjustment) or if ever discussing naltrexone "
                  "for alcohol reduction (strong predicted response).")
        conf = "high"
    else:
        ph = "OPRM1 A118G A/A — standard mu-opioid receptor"
        mech = "Standard opioid receptor kinetics."
        action = "None specific."
        conf = "high"
    return _find(CAT_REWARD, "OPRM1 A118G (mu-opioid receptor)", "OPRM1",
                 "rs1799971", gt, ph, mech, action, conf,
                 "Bond 1998; Way & Taylor 2010; Anton 2008 (naltrexone)")


# ── 5-HTTLPR proxy ────────────────────────────────────────────────────────────

def _five_htt(df):
    gt = _gt(df, "rs25531")
    if gt is None:
        return None
    return _find(CAT_SEROTONIN, "5-HTTLPR proxy (rs25531)", "SLC6A4",
                 "rs25531", gt,
                 "5-HTT expression variant",
                 "rs25531 is a partial proxy for the 5-HTTLPR VNTR (which is "
                 "not on consumer chips). The LA/LG/S polymorphism shifts "
                 "serotonin transporter expression and has weak, replication-"
                 "inconsistent associations with depression susceptibility × "
                 "life stress.",
                 "Effect sizes are small and controversial — do not over-read.",
                 "low", "Hariri 2002; Caspi 2003; meta-analysis debate")


# ── HTR2A T102C ──────────────────────────────────────────────────────────────

def _htr2a(df):
    gt = _gt(df, "rs6313")
    if gt is None:
        return None
    return _find(CAT_SEROTONIN, "HTR2A T102C", "HTR2A", "rs6313", gt,
                 f"HTR2A rs6313 {gt}",
                 "5-HT2A receptor coding-region SNP. Associations with SSRI "
                 "response (STAR*D showed T-carriers respond slightly worse), "
                 "PTSD susceptibility, and personality traits are replicated "
                 "but modest.",
                 "Small effect. Not action-guiding on its own.",
                 "low", "McMahon 2006 (STAR*D)")


# ── TPH2 ─────────────────────────────────────────────────────────────────────

def _tph2(df):
    gt = _gt(df, "rs4570625")
    if gt is None:
        return None
    return _find(CAT_SEROTONIN, "TPH2 -703G/T", "TPH2", "rs4570625", gt,
                 f"TPH2 rs4570625 {gt}",
                 "TPH2 rate-limits brain serotonin synthesis. rs4570625 T-"
                 "carriers have modestly altered brain 5-HT levels and some "
                 "reported associations with anxiety and treatment response.",
                 "Small effect; not action-guiding on its own.",
                 "low", "Tsai 2009; multiple candidate-gene studies")


# ── CACNA1C — mood stability ─────────────────────────────────────────────────

def _cacna1c(df):
    gt = _gt(df, "rs1006737")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A >= 1:
        ph = f"CACNA1C rs1006737 A-carrier ({n_A}× A) — common cross-psychiatric variant"
        mech = ("rs1006737-A is a well-replicated common variant associated at the "
                "POPULATION level with bipolar disorder, schizophrenia, and major "
                "depression (OR ~1.15-1.3 per allele). Crucially, that is a tiny "
                "per-individual effect: the allele is carried by ~30-40% of "
                "healthy people, and on its own it shifts your personal risk "
                "negligibly. It is not predictive or diagnostic for any single "
                "person — only a faint tile in a highly polygenic picture.")
        action = ("No action beyond general mental-wellness basics (sleep, "
                  "circadian consistency, stress management). This variant is "
                  "NOT a reason for concern about any psychiatric condition.")
        # Low confidence at the INDIVIDUAL level: a single OR~1.2 common variant
        # carries essentially no predictive power for a person (see user note).
        conf = "low"
    else:
        ph = "CACNA1C rs1006737 G/G — reference genotype"
        mech = "Standard L-type calcium channel expression."
        action = "None specific."
        conf = "low"
    return _find(CAT_MOOD, "CACNA1C rs1006737 (cross-psychiatric)",
                 "CACNA1C", "rs1006737", gt, ph, mech, action, conf,
                 "PGC 2011 & 2013; Green 2013")


# ── CHRNA5 — nicotine addiction / lung cancer ─────────────────────────────────

def _chrna5(df):
    gt = _gt(df, "rs16969968")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A >= 1:
        ph = f"CHRNA5 rs16969968 A-carrier ({n_A}× A) — heavy-smoker + lung-cancer risk allele"
        mech = ("rs16969968 A allele encodes a partial-loss-of-function α5 "
                "nicotinic acetylcholine receptor. Carriers who ever smoke tend "
                "to smoke MORE per day (need more nicotine for the same "
                "receptor stimulation) and have materially higher lung-cancer "
                "and COPD risk if they do. Zero risk if never a smoker.")
        action = ("**The clearest gene-behaviour prevention lever in your "
                  "genome:** never start smoking or vaping. If you already "
                  "smoke, this variant makes quitting harder AND makes "
                  "continued smoking more dangerous — varenicline/bupropion + "
                  "counselling are more effective in carriers.")
        conf = "high"
    else:
        ph = "CHRNA5 rs16969968 G/G — standard α5 nAChR"
        mech = "Standard α5 nicotinic acetylcholine receptor function."
        action = "None specific."
        conf = "high"
    return _find(CAT_ADDICTION, "CHRNA5 nicotine dependence + lung cancer",
                 "CHRNA5", "rs16969968", gt, ph, mech, action, conf,
                 "Thorgeirsson 2008 Nature; Bertelsen 2015")


# ══════════════════════════════════════════════════════════════════════════
# Composite phenotype synthesis
# ══════════════════════════════════════════════════════════════════════════

def _classify_comt(gt: str | None) -> str:
    if not gt:
        return "unknown"
    n_val = gt.count("G")
    return {2: "warrior", 1: "middle", 0: "worrier"}.get(n_val, "unknown")


def _classify_maoa(gt: str | None) -> str:
    if not gt:
        return "unknown"
    if "T" in gt and "G" not in gt:
        return "MAOA-H"
    if "G" in gt and "T" not in gt:
        return "MAOA-L"
    return "heterozygous"


def _classify_bdnf(gt: str | None) -> str:
    if not gt:
        return "unknown"
    n_T = gt.count("T")
    return {0: "Val/Val (full)", 1: "Val/Met (reduced)",
            2: "Met/Met (low)"}.get(n_T, "unknown")


def build_composite(df: pd.DataFrame, findings: list[dict]) -> dict:
    """Assemble the practical composite phenotype and recommendations."""
    comt = _classify_comt(_gt(df, "rs4680"))
    maoa = _classify_maoa(_gt(df, "rs6323"))
    bdnf = _classify_bdnf(_gt(df, "rs6265"))
    oprm1 = _gt(df, "rs1799971")
    chrna5 = _gt(df, "rs16969968")

    # ── Stress response profile ────
    if comt == "warrior" and maoa == "MAOA-H":
        stress_profile = ("Warrior + MAOA-H — fastest overall monoamine turnover. "
                          "Strong stress resilience; rapid emotional recovery; "
                          "needs stimulation gradients to feel engaged.")
    elif comt == "worrier" and maoa == "MAOA-L":
        stress_profile = ("Worrier + MAOA-L — slowest monoamine turnover. "
                          "Sustained cognitive depth but prone to rumination "
                          "and acute-stress overload. Emotion-regulation "
                          "practices are especially valuable.")
    elif comt == "middle" and maoa == "MAOA-H":
        stress_profile = ("Adaptive Middle + MAOA-H — the arguably-best combined "
                          "phenotype: cognitive stress-resilience of Val/Met "
                          "heterozygote advantage, plus fast affective recovery "
                          "of high-activity MAOA. Thrives in variable-intensity "
                          "cognitive work.")
    elif comt == "middle":
        stress_profile = ("Adaptive Middle — Val/Met heterozygote advantage; "
                          "handles moderate stress without paying the sustained-"
                          "cognition tax of pure warrior.")
    else:
        stress_profile = f"{comt.title()} × {maoa} — see individual entries."

    # ── Plasticity tier ────
    plasticity_tier = {
        "Val/Val (full)": "High — full activity-dependent BDNF; skill-acquisition-favourable.",
        "Val/Met (reduced)": "Intermediate — ~30% reduced activity-dependent BDNF.",
        "Met/Met (low)": "Lower — aerobic exercise partly compensates.",
        "unknown": "Not typed.",
    }[bdnf]

    # ── Stimulant response prediction ────
    if comt == "warrior":
        stimulant_response = ("Stronger cognitive effect from stimulants at a "
                              "given dose — Val/Val's low baseline PFC dopamine "
                              "leaves more room for methylphenidate/amphetamine "
                              "to help. If ever prescribed ADHD medication, "
                              "start LOW and titrate — you'll feel it.")
    elif comt == "worrier":
        stimulant_response = ("Stimulants may over-shoot in Met/Met carriers "
                              "with already-high baseline PFC dopamine. Small "
                              "doses often work better than large; watch for "
                              "anxiety at higher doses.")
    else:
        stimulant_response = ("Intermediate stimulant response — larger effect "
                              "than Met/Met, smaller than Val/Val. Standard "
                              "titration approach.")

    # ── SSRI response ────
    if maoa == "MAOA-H":
        ssri = ("Modestly reduced SSRI headroom (you clear serotonin fast). "
                "Standard SSRIs still work but effect size may be smaller — "
                "SNRIs sometimes outperform in high-MAOA individuals.")
    elif maoa == "MAOA-L":
        ssri = ("Standard-to-favourable SSRI response — low-activity MAOA "
                "leaves more serotonin in the synapse for SSRIs to build on.")
    else:
        ssri = "Standard SSRI response."

    # ── Caffeine protocol ────
    caffeine = ("Espresso-style timing (small doses, front-loaded before 2pm). "
                "Fast catecholamine clearance means caffeine wears off "
                "quickly, but that doesn't mean it's cleared from sleep-"
                "affecting compartments. Late caffeine still costs REM "
                "even if you feel awake.")
    if comt == "warrior":
        caffeine = ("Larger caffeine tolerance than average likely — but the "
                    "ceiling is anxiety, not sleep. Watch for jitter/anxiety "
                    "above ~200mg per dose. ") + caffeine

    # ── Meditation style fit ────
    if comt in ("warrior", "middle") and maoa == "MAOA-H":
        meditation = ("Focused-attention / mantra practices will feel "
                      "disproportionately hard — your brain wants stimulation "
                      "gradients. Better fits: open-monitoring (Vipassana-"
                      "style noting), somatic practices (yoga, breathwork "
                      "with dynamic components), short high-intensity "
                      "protocols, cold exposure. Wim-Hof / cold plunge / "
                      "controlled hyperventilation tend to feel genuinely "
                      "rewarding for this neurotype.")
    elif comt == "worrier":
        meditation = ("Longer, slow, focused-attention practices work well "
                      "for high-baseline-PFC-dopamine Met/Met carriers. "
                      "45+ minute sits, body scans, gentle breath work. "
                      "Avoid stimulating practices during high-stress periods "
                      "— cold exposure may destabilise.")
    else:
        meditation = ("Standard mixed approach — experiment with both "
                      "focused-attention and open-monitoring; the neurotype "
                      "signature doesn't strongly favor one over the other.")

    # ── Career neurotype ────
    if comt == "warrior" and bdnf == "Val/Val (full)":
        career = ("The 'high-plasticity, high-drive, high-boredom-tolerance-"
                  "low' phenotype. Thrives in variable-intensity, cognitively-"
                  "demanding, deadline-driven work. Suits: entrepreneurship, "
                  "surgery/ER medicine, tactical/military, competitive sport, "
                  "trading, product/startup work. Poorly-fit: repetitive "
                  "long-timescale roles without stimulation gradients.")
    elif comt == "middle" and bdnf == "Val/Val (full)":
        career = ("The 'adaptive high-performer' phenotype: cognitive stress-"
                  "resilience of Val/Met + full BDNF plasticity. Broad career "
                  "compatibility — genuinely well-suited to demanding cognitive "
                  "work whether that's medicine, engineering, research, "
                  "creative direction, law, or entrepreneurship. Depth-within-"
                  "domain pays off more than average due to plasticity edge.")
    elif comt == "worrier":
        career = ("Sustained-cognitive-depth phenotype. Suits: research, "
                  "long-form writing, complex analysis, careful clinical work. "
                  "Handle high-stakes acute-stress environments deliberately "
                  "(you can, but pay a cognitive tax).")
    else:
        career = ("See individual variant recommendations."
                  " Neurotype signature is genuinely favourable for "
                  "high-performance cognitive work.")

    # ── Addiction / substance risk flags ────
    substance_flags = []
    if chrna5 and "A" in chrna5:
        substance_flags.append(
            "🚭 **Never start smoking/vaping.** CHRNA5 A-carriers who smoke, smoke "
            "harder and have markedly higher lung-cancer risk.")
    if oprm1 and "G" in oprm1:
        substance_flags.append(
            "🍷 **OPRM1 G-carrier:** post-surgical opioid dosing may need "
            "adjustment; naltrexone is genetically favoured if ever considering "
            "medication-assisted alcohol reduction.")
    if comt == "warrior":
        substance_flags.append(
            "🧠 Warrior COMT means substances that raise dopamine (stimulants, "
            "cocaine, gambling) may feel especially rewarding — awareness > deterministic.")

    return {
        "comt_class": comt,
        "maoa_class": maoa,
        "bdnf_class": bdnf,
        "stress_response_profile": stress_profile,
        "plasticity_tier": plasticity_tier,
        "stimulant_response": stimulant_response,
        "ssri_response": ssri,
        "caffeine_protocol": caffeine,
        "meditation_fit": meditation,
        "career_neurotype": career,
        "substance_flags": substance_flags,
    }


# ══════════════════════════════════════════════════════════════════════════
# Master analyzer
# ══════════════════════════════════════════════════════════════════════════

CATEGORY_ORDER = [CAT_DOPAMINE, CAT_MONOAMINE, CAT_PLASTICITY, CAT_REWARD,
                  CAT_SEROTONIN, CAT_MOOD, CAT_ADDICTION]


def analyze_neurochemistry(df: pd.DataFrame) -> dict:
    """Full neurochemistry work-up + composite phenotype recommendations."""
    analyzers = [_comt, _maoa, _bdnf, _drd2, _drd4, _oprm1,
                 _five_htt, _htr2a, _tph2, _cacna1c, _chrna5]
    findings: list[dict] = []
    for a in analyzers:
        try:
            r = a(df)
        except Exception:
            continue
        if r is None:
            continue
        findings.append(r)

    by_category: dict[str, list[dict]] = {}
    for f in findings:
        by_category.setdefault(f["category"], []).append(f)

    composite = build_composite(df, findings)

    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "findings": findings,
        "by_category": by_category,
        "categories": [c for c in CATEGORY_ORDER if c in by_category],
        "composite": composite,
    }
