"""
Addiction Genetics — alcohol, opioid, nicotine, cannabis susceptibility
=======================================================================

Dedicated module for well-established addiction-related loci. Reports each
substance's per-locus genotype, aggregates into a substance-specific
susceptibility profile, and produces practical clinical flags:

  • **Alcohol** — ADH1B rs1229984, ALDH2 rs671, ADH1C rs1693482, CYP2E1
    rs2031920, GABRA2 rs279858. Metabolism-side protection (East Asian
    flush) plus dependence-side risk (GABRA2 haplotype). Cross-refs OPRM1
    for reward signalling.
  • **Opioid** — OPRM1 A118G rs1799971 (G-carriers: altered mu-opioid
    receptor kinetics, higher post-op morphine dose requirement, stronger
    naltrexone response — a *clinically useful* flag).
  • **Nicotine** — CHRNA5 D398N rs16969968 (heavy-smoker + lung-cancer risk
    if you ever start), CYP2A6 rs1801272 (metabolism rate).
  • **Cannabis** — CNR1 rs2023239 (dependence signal), FAAH rs324420
    (anandamide clearance; A-carriers have altered response and lower
    anxiety at baseline).
  • **General reward / impulse control** — DAT1 rs27072, DRD4 rs1800955,
    COMT rs4680, MAOA rs6323 (cross-referenced from neurochemistry).
  • **Stress × substance-use** — CRHR1 rs110402, FKBP5 rs1360780 (gene ×
    childhood-adversity interaction, per Binder 2008 & Klengel 2013).

Composite output
----------------
Per-substance susceptibility tier, overall risk profile, and a list of
*clinically-useful* flags: naltrexone-response prediction, opioid-dosing
adjustment, never-smoke warning, stress-substance interaction alerts.

Educational, not diagnostic. Genetics contributes ~50% of variance in
substance-use-disorder risk in population studies, but that's variance
across the population, not deterministic within an individual. Behaviour,
environment, and life circumstances dominate individual outcomes.

References
----------
Edenberg 2004 (ADH1B/ALDH2 & alcoholism); Bierut 2010 (CHRNA5); Bond 1998
(OPRM1 A118G); Anton 2008 (OPRM1 × naltrexone response); Chiara 2003 (DRD2 /
reward deficiency); Sabol 1998, Caspi 2002 (MAOA × environment); Binder 2008,
Klengel 2013 (FKBP5 × trauma); Edenberg 2004, Covault 2004 (GABRA2 &
alcoholism); Sipe 2002 (FAAH C385A).
"""

from __future__ import annotations

import pandas as pd

CAT_ALCOHOL = "Alcohol — Metabolism & Dependence"
CAT_OPIOID = "Opioid & Endogenous-reward"
CAT_NICOTINE = "Nicotine & Stimulant"
CAT_CANNABIS = "Cannabis & Endocannabinoid"
CAT_STRESS = "Stress × Substance-use"


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


def _f(cat, name, gene, rsid, gt, impact, verdict, mechanism, action,
       confidence, citation):
    return {"category": cat, "name": name, "gene": gene, "rsid": rsid,
            "genotype": gt or "—", "impact": impact, "verdict": verdict,
            "mechanism": mechanism, "action": action,
            "confidence": confidence, "citation": citation}


# ─── ALCOHOL ────────────────────────────────────────────────────────────────

def _adh1b(df):
    gt = _gt(df, "rs1229984")
    if gt is None:
        return None
    n_A = gt.count("A")   # A = His47 (fast enzyme → aversive flush precursor)
    if n_A >= 2:
        return _f(CAT_ALCOHOL, "ADH1B*2 (Arg47His) — strong protection",
                  "ADH1B", "rs1229984", gt, "protective",
                  "STRONG protection against alcohol-use disorder",
                  "Homozygous ADH1B*2 (fast enzyme). Rapid ethanol→acetaldehyde "
                  "conversion → flushing, nausea and aversive drinking response. "
                  "One of the strongest single-variant protective effects known "
                  "in humans against alcoholism (OR ~0.3-0.4 for AUD). Near-fixed "
                  "in East Asians; rare in Europeans.",
                  "Drinking is unlikely to feel rewarding; this is protective.",
                  "high", "Edenberg 2004; Bierut 2012")
    if n_A == 1:
        return _f(CAT_ALCOHOL, "ADH1B*2 heterozygous — modest alcohol protection",
                  "ADH1B", "rs1229984", gt, "protective",
                  "Modest per-allele protection against AUD",
                  "One copy of the fast-enzyme ADH1B*2 allele. Partial "
                  "acetaldehyde-accumulation phenotype; drinking may feel less "
                  "purely rewarding.",
                  "Standard alcohol-use awareness applies; genetics tilts slightly favorably.",
                  "high", "Edenberg 2004")
    return _f(CAT_ALCOHOL, "ADH1B — standard European allele", "ADH1B",
              "rs1229984", gt, "neutral",
              "No genetic protection from fast alcohol → acetaldehyde",
              "Standard ADH1B (Arg47/Arg47). Ethanol metabolism proceeds "
              "at typical European rate; no built-in flush deterrent.",
              "Behavioural moderation is not offset by biology here.",
              "high", "Edenberg 2004")


def _aldh2(df):
    gt = _gt(df, "rs671")
    if gt is None:
        return None
    n_A = gt.count("A")   # A = Glu487Lys = LOF
    if n_A >= 1:
        return _f(CAT_ALCOHOL, "ALDH2*2 (Asian flush) — very strong protection",
                  "ALDH2", "rs671", gt, "protective",
                  "Very strong protection against alcohol-use disorder",
                  f"{'Homozygous' if n_A == 2 else 'Heterozygous'} for the "
                  "ALDH2*2 loss-of-function allele. Even one copy severely "
                  "impairs acetaldehyde clearance → dramatic flushing, "
                  "tachycardia, nausea after minimal alcohol. Homozygotes "
                  "essentially cannot drink. Near-fixed as a protective allele "
                  "in East Asia; rare in Europeans.",
                  "Alcohol will be aversive; also flags MASSIVELY higher esophageal-"
                  "cancer risk with any regular drinking (acetaldehyde is a Group 1 "
                  "carcinogen).",
                  "high", "Higuchi 1994; Brooks 2009 (esophageal cancer)")
    return _f(CAT_ALCOHOL, "ALDH2 wild-type — standard acetaldehyde clearance",
              "ALDH2", "rs671", gt, "neutral",
              "Standard alcohol tolerance",
              "Wild-type ALDH2; normal acetaldehyde clearance. No 'Asian flush'.",
              "None specific.", "high", "Higuchi 1994")


def _cyp2e1(df):
    gt = _gt(df, "rs2031920")
    if gt is None:
        return None
    return _f(CAT_ALCOHOL, "CYP2E1 c1/c2", "CYP2E1", "rs2031920", gt,
              "informational",
              f"Microsomal ethanol oxidation genotype {gt}",
              "CYP2E1 handles ethanol at higher blood-alcohol levels and "
              "generates oxidative stress. Variants modestly alter "
              "susceptibility to alcohol-related liver injury.",
              "Small effect; do not over-read.", "moderate", "Yin 2007")


def _gabra2(df):
    gt = _gt(df, "rs279858")
    if gt is None:
        return None
    n_C = gt.count("C")   # C is the risk allele in most meta-analyses
    if n_C >= 1:
        return _f(CAT_ALCOHOL, "GABRA2 — alcohol-dependence risk haplotype",
                  "GABRA2", "rs279858", gt,
                  "susceptible",
                  f"{n_C}× C — carries alcohol-dependence risk allele (OR ~1.3-1.5)",
                  "The GABRA2 risk haplotype has been one of the most replicated "
                  "loci for adult alcohol dependence (Edenberg 2004; COGA cohort). "
                  "It's a *susceptibility* signal, not deterministic — modifies "
                  "GABA-A receptor subunit assembly and alcohol's rewarding "
                  "effect. Interacts with childhood conduct problems and adverse "
                  "environments.",
                  "Not a reason to avoid alcohol; a reason to be aware of "
                  "personal-relationship-with-alcohol patterns as they develop. "
                  "If binge/craving patterns ever emerge, this variant is one "
                  "reason to take them seriously earlier than the population "
                  "average.",
                  "moderate", "Edenberg 2004; Covault 2004")
    return _f(CAT_ALCOHOL, "GABRA2 wild-type", "GABRA2", "rs279858", gt,
              "neutral", "No GABRA2 dependence-risk allele",
              "Common protective allele at this locus.",
              "None specific.", "moderate", "Edenberg 2004")


# ─── OPIOID / REWARD ────────────────────────────────────────────────────────

def _oprm1(df):
    gt = _gt(df, "rs1799971")
    if gt is None:
        return None
    n_G = gt.count("G")   # G = N40D, altered receptor
    if n_G >= 1:
        return _f(CAT_OPIOID, "OPRM1 A118G — altered mu-opioid receptor",
                  "OPRM1", "rs1799971", gt,
                  "clinically-relevant",
                  f"{n_G}× G-allele — clinically useful for pain management + naltrexone",
                  "The G118 (N40D) allele produces a mu-opioid receptor with "
                  "higher affinity for beta-endorphin. Three clinical implications: "
                  "(1) G-carriers often need HIGHER morphine doses for equivalent "
                  "post-surgical analgesia; (2) naltrexone works BETTER in "
                  "G-carriers (Anton 2008 — one of the few positive pharmacogenetic "
                  "predictors in addiction medicine); (3) altered alcohol reward "
                  "signalling — heterogeneous behavioural effect but real.",
                  "**Flag this to any anesthesiologist before surgery** (may need "
                  "higher opioid dose). If alcohol-use ever becomes a concern, "
                  "naltrexone is a genetically-favoured therapy.",
                  "high", "Bond 1998; Way & Taylor 2010; Anton 2008")
    return _f(CAT_OPIOID, "OPRM1 wild-type (A/A)", "OPRM1", "rs1799971", gt,
              "neutral", "Standard mu-opioid receptor kinetics",
              "Standard opioid receptor function.",
              "None specific.", "high", "Bond 1998")


def _dat1(df):
    gt = _gt(df, "rs27072")
    if gt is None:
        return None
    n_T = gt.count("T")
    return _f(CAT_OPIOID, "DAT1 (SLC6A3) dopamine transporter", "SLC6A3",
              "rs27072", gt, "informational",
              f"DAT1 genotype {gt}",
              "The DAT1 dopamine transporter clears dopamine from the synapse. "
              "rs27072 variants have been associated with alcohol withdrawal "
              "severity and stimulant response, but effect sizes are small.",
              "Not action-guiding on its own.", "low", "Ueno 1999")


# ─── NICOTINE ───────────────────────────────────────────────────────────────

def _chrna5(df):
    """Duplicates neurochemistry's CHRNA5 finding but re-framed for the
    addiction lens (and to make this module standalone)."""
    gt = _gt(df, "rs16969968")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A >= 1:
        return _f(CAT_NICOTINE, "CHRNA5 D398N — heavy-smoker + lung-cancer risk",
                  "CHRNA5", "rs16969968", gt,
                  "susceptible",
                  f"{n_A}× A — 'never start' is the clearest genetic prevention lever",
                  "The A allele of CHRNA5 rs16969968 encodes a partial-loss-of-"
                  "function α5 nicotinic receptor. Carriers who smoke smoke MORE "
                  "per day (need more nicotine for equivalent receptor stimulation) "
                  "AND have materially higher lung-cancer + COPD risk. Zero effect "
                  "if never a smoker.",
                  "**Never smoke or vape, including socially.** If you already "
                  "smoke, this makes quitting harder AND makes continued smoking "
                  "more dangerous — varenicline / bupropion + counseling are "
                  "genetically-favored.",
                  "high", "Thorgeirsson 2008 Nature; Bierut 2010")
    return None


def _cyp2a6(df):
    gt = _gt(df, "rs1801272")
    if gt is None:
        return None
    # rs1801272 is a proxy for CYP2A6*2 (LOF, slow nicotine metabolism)
    n_A = gt.count("A")
    if n_A >= 1:
        return _f(CAT_NICOTINE, "CYP2A6 slow-nicotine metabolism",
                  "CYP2A6", "rs1801272", gt,
                  "protective",
                  "Slower nicotine metabolism — reduced smoking-initiation risk",
                  "Carriers of the slow-metabolism CYP2A6 variant get more "
                  "nicotine effect per cigarette, tend to smoke fewer cigarettes, "
                  "and have lower dependence risk if they ever start.",
                  "Reinforces the 'never smoke' rationale — but doesn't remove "
                  "CHRNA5 lung-cancer risk if you do.",
                  "moderate", "Pianezza 1998; Malaiyandi 2005")
    return None


# ─── CANNABIS ───────────────────────────────────────────────────────────────

def _cnr1(df):
    gt = _gt(df, "rs2023239")
    if gt is None:
        return None
    n_C = gt.count("C")
    if n_C >= 1:
        return _f(CAT_CANNABIS, "CNR1 — cannabis-response variant", "CNR1",
                  "rs2023239", gt, "informational",
                  f"{n_C}× C — modest cannabis-dependence signal",
                  "CNR1 encodes the CB1 cannabinoid receptor. rs2023239 C-carriers "
                  "have modestly elevated cannabis-dependence risk in some cohorts "
                  "but effects are small and inconsistent.",
                  "Not action-guiding. Standard cannabis-use awareness applies.",
                  "low", "Zhang 2004; Hopfer 2006")
    return None


def _faah(df):
    gt = _gt(df, "rs324420")
    if gt is None:
        return None
    n_A = gt.count("A")   # A = C385A missense; reduced FAAH activity → higher anandamide
    if n_A >= 1:
        return _f(CAT_CANNABIS, "FAAH C385A — higher endogenous anandamide",
                  "FAAH", "rs324420", gt,
                  "informational",
                  f"{n_A}× A — reduced FAAH activity → higher anandamide baseline",
                  "The A allele produces a less-stable FAAH enzyme, leaving more "
                  "anandamide (the brain's endogenous cannabinoid) in circulation. "
                  "A-carriers have been shown in neuroimaging studies to have "
                  "modestly LOWER anxiety, altered fear-extinction, and altered "
                  "response to exogenous cannabinoids.",
                  "Awareness only; small effects.", "moderate",
                  "Sipe 2002; Dincheva 2015")
    return _f(CAT_CANNABIS, "FAAH wild-type", "FAAH", "rs324420", gt,
              "neutral", "Standard anandamide clearance",
              "Standard FAAH enzyme activity.", "None specific.",
              "moderate", "Sipe 2002")


# ─── STRESS × SUBSTANCE-USE ─────────────────────────────────────────────────

def _crhr1(df):
    gt = _gt(df, "rs110402")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A >= 1:
        return _f(CAT_STRESS, "CRHR1 — stress-drinking interaction locus",
                  "CRHR1", "rs110402", gt, "susceptible",
                  f"{n_A}× A — stress-induced-drinking interaction allele",
                  "CRHR1 rs110402 A-carriers show elevated stress-induced "
                  "alcohol consumption in cohort studies. Effect emerges under "
                  "stress; largely neutral otherwise.",
                  "Actively manage acute-stress → alcohol coupling if a pattern "
                  "emerges. Not a general drinking-risk flag.",
                  "moderate", "Blomeyer 2008; Treutlein 2006")
    return _f(CAT_STRESS, "CRHR1 wild-type", "CRHR1", "rs110402", gt,
              "neutral", "No CRHR1 stress-drinking risk allele",
              "Common protective allele.", "None specific.",
              "moderate", "Blomeyer 2008")


def _fkbp5(df):
    gt = _gt(df, "rs1360780")
    if gt is None:
        return None
    n_T = gt.count("T")
    if n_T >= 1:
        return _f(CAT_STRESS, "FKBP5 — trauma × substance-use interaction locus",
                  "FKBP5", "rs1360780", gt, "conditional",
                  f"{n_T}× T — interaction allele; matters only with childhood adversity",
                  "FKBP5 rs1360780 T-carriers have impaired glucocorticoid-"
                  "receptor feedback; combined with significant childhood "
                  "adversity, this locus substantially raises substance-use "
                  "disorder, PTSD, and MDD risk (Binder 2008; Klengel 2013). "
                  "In the absence of severe childhood adversity, essentially "
                  "neutral. The gene × environment interaction is the key.",
                  "If a personal history of significant childhood adversity "
                  "exists, standard trauma-informed care is materially more "
                  "important given this locus. Otherwise no specific action.",
                  "moderate", "Binder 2008 JAMA; Klengel 2013 Nat Neurosci")
    return _f(CAT_STRESS, "FKBP5 wild-type", "FKBP5", "rs1360780", gt,
              "neutral", "No FKBP5 trauma-interaction allele",
              "Common protective allele; normal glucocorticoid feedback.",
              "None specific.", "moderate", "Binder 2008")


# ─── Composite ──────────────────────────────────────────────────────────────

def _tally(findings, impact):
    return sum(1 for f in findings if f.get("impact") == impact)


def build_composite(findings: list[dict]) -> dict:
    protective = _tally(findings, "protective")
    susceptible = _tally(findings, "susceptible")
    conditional = _tally(findings, "conditional")

    # Alcohol-specific tier
    alc_prot = sum(1 for f in findings if f["category"] == CAT_ALCOHOL
                    and f["impact"] == "protective")
    alc_susc = sum(1 for f in findings if f["category"] == CAT_ALCOHOL
                    and f["impact"] == "susceptible")
    if alc_prot >= 2:
        alcohol_tier = ("Strongly protected",
                        "Multiple genetic protections against alcohol-use disorder. "
                        "Drinking is likely aversive or of low reward.")
    elif alc_prot == 1 and alc_susc == 0:
        alcohol_tier = ("Modestly protected",
                        "One protective variant, no risk variants. Slight genetic tilt "
                        "away from problem drinking.")
    elif alc_susc >= 1 and alc_prot == 0:
        alcohol_tier = ("Modestly susceptible",
                        f"{alc_susc} dependence-risk variant(s), no metabolism protection. "
                        "Genetics tilts modestly toward alcohol-use disorder — behaviour "
                        "and context dominate.")
    elif alc_prot == 0 and alc_susc == 0:
        alcohol_tier = ("Baseline (typical European)",
                        "Neither strong protection nor strong risk. Standard European "
                        "adult risk (~10-15% lifetime AUD). Behaviour dominates.")
    else:
        alcohol_tier = ("Mixed",
                        f"Both protective ({alc_prot}) and risk ({alc_susc}) variants — "
                        "net effect is likely close to baseline.")

    # Clinical flags — the useful ones
    flags = []
    for f in findings:
        if f["gene"] == "OPRM1" and "clinically" in f["impact"]:
            flags.append(("🩹 Opioid dosing",
                          "Flag OPRM1 status to any anesthesiologist before surgery — "
                          "you may need higher opioid doses for equivalent analgesia."))
            flags.append(("💊 Naltrexone response",
                          "OPRM1 G-carriers respond better to naltrexone for alcohol-"
                          "reduction. If problem drinking ever emerges, this is a "
                          "genetically-favored intervention."))
        if f["gene"] == "CHRNA5" and f["impact"] == "susceptible":
            flags.append(("🚭 Never start smoking/vaping",
                          "CHRNA5 A-carrier + smoker = much higher lung-cancer and "
                          "COPD risk. If never a smoker, zero effect — currently your "
                          "most-realised genetic prevention win."))
        if f["gene"] == "ALDH2" and f["impact"] == "protective":
            flags.append(("⚠️ Esophageal cancer",
                          "ALDH2*2 carriers who drink regularly have much higher "
                          "esophageal-cancer risk — treat any drinking as high-cost."))
        if f["gene"] == "FKBP5" and f["impact"] == "conditional":
            flags.append(("🧠 Trauma-informed care",
                          "FKBP5 T-carrier: substance-use / PTSD / depression risk "
                          "materially elevated if paired with childhood adversity. "
                          "Standard trauma-informed care is more valuable given this locus."))

    # Overall
    if protective >= 2 and susceptible == 0:
        overall = ("Low overall susceptibility",
                   "Multiple genetic protections, no strong risk variants. Below-average "
                   "genetic susceptibility to substance-use disorders.")
    elif susceptible >= 2 and protective == 0:
        overall = ("Elevated susceptibility",
                   "Multiple substance-related risk variants. Behavioural moderation "
                   "and mindful use patterns are especially important — genes tilt "
                   "the table.")
    elif protective >= 1 and susceptible >= 1:
        overall = ("Mixed",
                   "Some protective and some susceptible variants. Net effect near "
                   "baseline European risk (~10-15% lifetime substance-use disorder).")
    else:
        overall = ("Baseline",
                   "No strong genetic tilts detected. Standard European adult risk. "
                   "Behaviour dominates.")

    return {
        "alcohol_tier": alcohol_tier[0], "alcohol_narrative": alcohol_tier[1],
        "overall_tier": overall[0], "overall_narrative": overall[1],
        "n_protective": protective, "n_susceptible": susceptible,
        "n_conditional": conditional,
        "clinical_flags": [{"title": t, "text": s} for t, s in flags],
    }


CATEGORY_ORDER = [CAT_ALCOHOL, CAT_OPIOID, CAT_NICOTINE, CAT_CANNABIS, CAT_STRESS]


def analyze_addiction_genetics(df: pd.DataFrame) -> dict:
    analyzers = [
        _adh1b, _aldh2, _cyp2e1, _gabra2,          # alcohol
        _oprm1, _dat1,                             # opioid / reward
        _chrna5, _cyp2a6,                          # nicotine
        _cnr1, _faah,                              # cannabis
        _crhr1, _fkbp5,                            # stress × substance
    ]
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

    composite = build_composite(findings)

    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "findings": findings,
        "by_category": by_category,
        "categories": [c for c in CATEGORY_ORDER if c in by_category],
        "composite": composite,
    }
