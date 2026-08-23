"""
Detoxification & Environmental Resilience Panel
===============================================

A dedicated module for how *your* genome handles environmental insults, with a
deliberate focus on the two exposures you asked about:

  1. Airborne combustion products — **wildfire smoke**, wood smoke and tobacco
     smoke, which all deliver the same core toxicants: fine particulate matter
     (PM2.5), polycyclic aromatic hydrocarbons (PAHs), aromatic amines,
     benzene/quinones and a large oxidative-stress load.
  2. **Heavy metals** — lead, cadmium, arsenic, mercury, copper and iron.

Biology, in one paragraph: a toxicant is first *activated* by Phase I enzymes
(the CYP1 family, switched on by the AhR receptor) into reactive intermediates,
then *neutralised* by Phase II conjugation (glutathione-S-transferases, NAT2,
NQO1, epoxide hydrolase) and mopped up by the antioxidant system (the NRF2
master switch driving SOD2, GPX1, catalase, heme-oxygenase). The genotypes that
matter most are the ones that make Phase I *fast* while Phase II or the
antioxidant response is *slow* — you generate reactive, carcinogenic
intermediates but clear them poorly. This module scores that balance and turns
it into a personalised action plan.

This complements — it does not replace — ``metal_oxidative.py`` (Parkinson's /
metal-handling / oxidative loci) and ``wellness.py`` (SOD2/GPX1). Genes reused
from those modules are re-interpreted here through the specific lens of smoke and
metal exposure, never re-declared.

Nothing here is a diagnostic test. Consumer-chip calls are hints; a true
GSTM1/GSTT1 null needs a PCR/CNV assay, and blood-lead / heavy-metal burden is a
lab measurement, not a genotype. Every finding is labelled with a confidence
level and read conservatively.
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path

import pandas as pd

import snp_registry

# ── genotype helpers (strand-aware via the registry) ──────────────────────────

def _gt(snps_df: pd.DataFrame, rsid: str) -> str | None:
    if rsid not in snps_df.index:
        return None
    row = snps_df.loc[rsid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    gt = row.get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    return s or None


def _dose(snps_df: pd.DataFrame, rsid: str, risk_allele: str,
          ref_allele: str) -> int | None:
    """Strand-aware risk-allele dosage using the shared registry helper."""
    return snp_registry.risk_dose_from_df(
        snps_df, rsid, risk_allele=risk_allele, ref_allele=ref_allele)


# ── domain / category labels ──────────────────────────────────────────────────

CAT_ACTIVATION = "Phase I — Toxicant Activation (smoke / PAH)"
CAT_CONJUGATION = "Phase II — Conjugation & Clearance"
CAT_ANTIOX = "Antioxidant & Oxidative-Stress Response"
CAT_METAL = "Heavy-Metal Handling"

_IMPACT_ORDER = {"higher-load": 3, "reduced-clearance": 3, "reduced": 2,
                 "intermediate": 1, "typical": 0, "protective": -1}


def _finding(category: str, trait: str, gene: str, rsid: str, genotype: str,
             result: str, action: str, confidence: str, impact: str,
             evidence: str = "") -> dict:
    return {
        "category": category, "trait": trait, "gene": gene, "rsid": rsid,
        "genotype": genotype or "—", "result": result, "action": action,
        "confidence": confidence, "impact": impact,
        "evidence": evidence or f"{rsid} genotype {genotype or '—'}",
    }


# ─── Phase I — bioactivation of smoke/PAHs ─────────────────────────────────────

def _cyp1a1(snps):
    # rs1048943 (Ile462Val, m2) and rs4646903 (m1, 3'UTR). Val462 / m1 =
    # higher inducibility → more PAH bioactivation, especially in smokers.
    m2 = _dose(snps, "rs1048943", "G", "A")   # G = Val462 (fast)
    m1 = _dose(snps, "rs4646903", "C", "T")   # C = m1 variant (inducible)
    if m2 is None and m1 is None:
        return None
    hi = (m2 or 0) + (m1 or 0)
    parts = []
    if m2 is not None:
        parts.append(f"rs1048943 Val462 dose {m2}")
    if m1 is not None:
        parts.append(f"rs4646903 m1 dose {m1}")
    if hi >= 1:
        result = ("Carries a high-inducibility CYP1A1 allele. CYP1A1 is the "
                  "front-line enzyme that converts smoke- and PAH-derived "
                  "compounds into reactive intermediates. A more inducible "
                  "enzyme means smoke exposure generates more of these reactive "
                  "species — which is only a problem if Phase II clearance can't "
                  "keep up (see below).")
        impact = "higher-load"
    else:
        result = ("Typical-inducibility CYP1A1 genotype — no evidence of "
                  "unusually fast bioactivation of PAHs at these markers.")
        impact = "typical"
    action = ("The activation step isn't something you slow down directly; the "
              "lever is reducing exposure (smoke, char-grilled meat) and keeping "
              "Phase II / antioxidant capacity high.")
    return _finding(CAT_ACTIVATION, "CYP1A1 PAH Bioactivation", "CYP1A1",
                    "rs1048943", _gt(snps, "rs1048943") or _gt(snps, "rs4646903"),
                    result, action, "moderate", impact, "; ".join(parts))


def _cyp1a2(snps):
    # rs762551 *1F. A allele = rapid inducibility (esp. in smokers/coffee).
    a = _dose(snps, "rs762551", "A", "C")
    if a is None:
        return None
    if a >= 1:
        result = ("Carries the CYP1A2*1F 'rapid' allele. In smokers this enzyme "
                  "is strongly induced, accelerating both caffeine metabolism "
                  "and the activation of aromatic amines / PAHs from smoke.")
        impact = "higher-load"
    else:
        result = ("Slower CYP1A2 inducibility genotype (*1F absent) — less "
                  "smoke-driven induction of this enzyme.")
        impact = "typical"
    action = ("Relevant mainly under exposure — the practical takeaway is the "
              "same: minimise smoke, and support Phase II clearance.")
    return _finding(CAT_ACTIVATION, "CYP1A2 Inducibility (*1F)", "CYP1A2",
                    "rs762551", _gt(snps, "rs762551"), result, action,
                    "moderate", impact, f"rs762551 A(*1F) dose {a}")


def _cyp1b1(snps):
    # rs1056836 Leu432Val. Val = higher activity toward PAHs/estrogens.
    v = _dose(snps, "rs1056836", "G", "C")   # G = Val432
    if v is None:
        return None
    if v >= 1:
        result = ("Carries the CYP1B1 Val432 allele, associated with higher "
                  "activity toward PAHs and other combustion by-products — "
                  "another Phase I bioactivation route for smoke toxicants.")
        impact = "higher-load"
    else:
        result = ("CYP1B1 Leu432 genotype — typical Phase I activity at this "
                  "locus.")
        impact = "typical"
    action = ("No direct action; contributes to the overall activation side of "
              "your smoke-resilience balance.")
    return _finding(CAT_ACTIVATION, "CYP1B1 PAH Activation (L432V)", "CYP1B1",
                    "rs1056836", _gt(snps, "rs1056836"), result, action,
                    "low", impact, f"rs1056836 Val432 dose {v}")


def _ahr(snps):
    # rs2066853 Arg554Lys — modifies AhR signalling / CYP1 induction.
    k = _dose(snps, "rs2066853", "A", "G")   # A = Lys554
    if k is None:
        return None
    result = ("The aryl-hydrocarbon receptor (AhR) is the sensor that switches "
              "on the CYP1 enzymes when it detects smoke/dioxin-like compounds. "
              + ("This variant (Lys554) subtly alters that induction response."
                 if k >= 1 else
                 "Your genotype is the common Arg554 form."))
    action = ("Cruciferous vegetables (indole-3-carbinol / sulforaphane) also "
              "signal through AhR/NRF2 — a food-based lever on this pathway.")
    return _finding(CAT_ACTIVATION, "AhR Smoke-Sensor Signalling", "AHR",
                    "rs2066853", _gt(snps, "rs2066853"), result, action,
                    "low", "higher-load" if k >= 1 else "typical",
                    f"rs2066853 Lys554 dose {k}")


# ─── Phase II — conjugation / clearance ────────────────────────────────────────

def _ephx1(snps):
    # rs1051740 Y113H (exon3): His113 (C) = SLOW; rs2234922 H139R (exon4):
    # Arg139 (G) = FAST. Combined "slow" activity leaves PAH epoxides lingering.
    slow = _dose(snps, "rs1051740", "C", "T")   # C = His113 (slow)
    fast = _dose(snps, "rs2234922", "G", "A")   # G = Arg139 (fast)
    if slow is None and fast is None:
        return None
    # crude activity index: slow alleles reduce, fast alleles raise
    idx = (fast or 0) - (slow or 0)
    parts = []
    if slow is not None:
        parts.append(f"rs1051740 His113(slow) dose {slow}")
    if fast is not None:
        parts.append(f"rs2234922 Arg139(fast) dose {fast}")
    if idx <= -1:
        result = ("Predicted SLOW microsomal epoxide hydrolase (EPHX1). This "
                  "enzyme detoxifies the reactive epoxides that Phase I makes "
                  "from PAHs in wood/tobacco smoke. Slower activity means those "
                  "epoxides hang around longer before being cleared.")
        impact = "reduced-clearance"
    elif idx >= 1:
        result = ("Predicted FAST EPHX1 activity — efficient clearance of PAH "
                  "epoxides at this locus.")
        impact = "protective"
    else:
        result = ("Intermediate EPHX1 activity — typical epoxide-clearance "
                  "capacity.")
        impact = "intermediate"
    action = ("Supports the cruciferous-vegetable / sulforaphane rationale "
              "(induces the broader epoxide/GST detox network) and, above all, "
              "reducing smoke exposure so less epoxide is generated.")
    return _finding(CAT_CONJUGATION, "EPHX1 Epoxide Clearance", "EPHX1",
                    "rs1051740", _gt(snps, "rs1051740"), result, action,
                    "moderate", impact, "; ".join(parts))


def _nat2(snps):
    # Slow-acetylator alleles: *5 (rs1801280 C), *6 (rs1799930 A), *7
    # (rs1799931 A). Any two slow alleles ⇒ slow acetylator phenotype, which
    # handles aromatic amines from smoke more slowly (bladder-cancer relevant).
    a5 = _dose(snps, "rs1801280", "C", "T")
    a6 = _dose(snps, "rs1799930", "A", "G")
    a7 = _dose(snps, "rs1799931", "A", "G")
    if a5 is None and a6 is None and a7 is None:
        return None
    slow_alleles = (a5 or 0) + (a6 or 0) + (a7 or 0)
    parts = [p for p in (
        f"*5 dose {a5}" if a5 is not None else "",
        f"*6 dose {a6}" if a6 is not None else "",
        f"*7 dose {a7}" if a7 is not None else "") if p]
    if slow_alleles >= 2:
        result = ("Predicted SLOW NAT2 acetylator. NAT2 conjugates aromatic "
                  "amines — a class of bladder carcinogens abundant in cigarette "
                  "and, to a lesser degree, wildfire smoke. Slow acetylators "
                  "clear them less efficiently; the classic epidemiological "
                  "finding is higher smoking-related bladder-cancer risk.")
        impact = "reduced-clearance"
        conf = "moderate"
    elif slow_alleles == 1:
        result = ("Intermediate NAT2 acetylator (one slow allele) — between the "
                  "fast and slow phenotypes.")
        impact = "intermediate"
        conf = "moderate"
    else:
        result = ("Predicted FAST/normal NAT2 acetylator at the *5/*6/*7 "
                  "markers — efficient handling of aromatic amines.")
        impact = "typical"
        conf = "moderate"
    action = ("The strongest lever here is avoiding tobacco smoke entirely and "
              "limiting exposure to well-done/char-grilled meats (also a source "
              "of aromatic amines). NAT2 phenotype is a smoking-risk multiplier, "
              "not a stand-alone risk.")
    return _finding(CAT_CONJUGATION, "NAT2 Acetylator (aromatic amines)",
                    "NAT2", "rs1801280", _gt(snps, "rs1801280"), result, action,
                    conf, impact, "; ".join(parts))


def _nqo1(snps):
    # rs1800566 P187S: Ser187 (T on + strand) = reduced/absent NQO1 → poorer
    # detox of quinones (from benzene/smoke) and weaker antioxidant recycling.
    s = _dose(snps, "rs1800566", "T", "C")   # T = Ser187 (reduced)
    if s is None:
        return None
    if s >= 2:
        result = ("NQO1 187 Ser/Ser — greatly reduced NQO1 activity. NQO1 "
                  "detoxifies quinones (from benzene and smoke) and regenerates "
                  "antioxidant CoQ10/vitamin E. Low activity is a meaningful gap "
                  "in defence against combustion by-products.")
        impact = "reduced"
        conf = "moderate"
    elif s == 1:
        result = ("NQO1 187 Pro/Ser — intermediate NQO1 activity (roughly half "
                  "of normal).")
        impact = "intermediate"
        conf = "moderate"
    else:
        result = ("NQO1 187 Pro/Pro — full NQO1 quinone-detox activity.")
        impact = "typical"
        conf = "moderate"
    action = ("NQO1 is a core NRF2 target: sulforaphane (broccoli sprouts) and "
              "other cruciferous compounds up-regulate it. Reduced-activity "
              "carriers benefit most from a steady cruciferous intake and from "
              "minimising benzene/smoke exposure.")
    return _finding(CAT_CONJUGATION, "NQO1 Quinone Detox (P187S)", "NQO1",
                    "rs1800566", _gt(snps, "rs1800566"), result, action,
                    conf, impact, f"rs1800566 Ser187 dose {s}")


def _gst_conjugation(snps):
    # GSTM1/GSTT1 proxy SNPs (true null = whole-gene deletion, not callable
    # here) + GSTP1 rs1695 Ile105Val. GSTs conjugate glutathione onto activated
    # PAHs — the main Phase II route for smoke carcinogens.
    m1 = _dose(snps, "rs4147565", "A", "G")   # GSTM1 proxy
    t1 = _dose(snps, "rs4630", "A", "G")       # GSTT1 proxy
    p1 = _dose(snps, "rs1695", "G", "A")       # GSTP1 Val105 (G)
    if m1 is None and t1 is None and p1 is None:
        return None
    parts = [p for p in (
        f"GSTM1 proxy rs4147565 dose {m1}" if m1 is not None else "",
        f"GSTT1 proxy rs4630 dose {t1}" if t1 is not None else "",
        f"GSTP1 Val105 rs1695 dose {p1}" if p1 is not None else "") if p]
    reduced = (p1 or 0) >= 1
    result = ("Glutathione-S-transferases (GSTM1/GSTT1/GSTP1) are the workhorse "
              "Phase II enzymes that glue glutathione onto activated PAHs so they "
              "can be excreted. "
              + ("Your GSTP1 Val105 allele is associated with altered activity "
                 "toward PAH substrates. "
                 if reduced else
                 "Your GSTP1 genotype is the common Ile105 form. ")
              + "The clinically important GSTM1/GSTT1 'null' (enzyme entirely "
                "absent) is a whole-gene deletion these chip SNPs cannot confirm "
                "— treat the proxies as a weak hint only.")
    action = ("Glutathione is the shared currency of this whole system. Support "
              "it with cruciferous vegetables (sulforaphane induces GSTs), "
              "adequate protein (cysteine/glycine/glutamate) and, if exposure is "
              "high, discussing NAC with a clinician. A true GSTM1/GSTT1-null "
              "determination needs a PCR/CNV assay.")
    return _finding(CAT_CONJUGATION, "Glutathione-S-Transferase Capacity "
                    "(GSTM1/T1/P1)", "GSTM1/GSTT1/GSTP1", "rs1695",
                    _gt(snps, "rs1695"), result, action, "low",
                    "reduced" if reduced else "typical", "; ".join(parts))


# ─── Antioxidant / NRF2 response ───────────────────────────────────────────────

def _nrf2(snps):
    # rs6721961 (NFE2L2 promoter -617): T allele = lower NRF2 expression →
    # blunted induction of the entire Phase II / antioxidant battery.
    t = _dose(snps, "rs6721961", "T", "G")
    if t is None:
        return None
    if t >= 1:
        result = ("NFE2L2/NRF2 promoter variant (-617 T) associated with lower "
                  "baseline expression of NRF2 — the master switch that turns on "
                  "your entire antioxidant and Phase II defence when it senses "
                  "oxidative stress from smoke/particulates. A blunted switch "
                  "means a weaker built-in response.")
        impact = "reduced"
    else:
        result = ("Common NFE2L2/NRF2 genotype — typical inducibility of the "
                  "antioxidant master-switch.")
        impact = "typical"
    action = ("NRF2 is directly activated by sulforaphane (broccoli sprouts are "
              "the richest source), and by exercise and polyphenols. Lower-"
              "expression carriers get the most benefit from a consistent "
              "sulforaphane/cruciferous habit — it partly compensates for the "
              "weaker baseline.")
    return _finding(CAT_ANTIOX, "NRF2 Antioxidant Master-Switch", "NFE2L2",
                    "rs6721961", _gt(snps, "rs6721961"), result, action,
                    "moderate", impact, f"rs6721961 T dose {t}")


def _sod2(snps):
    # rs4880 A16V: Ala16 (C on + strand for the Ala-coding) targets SOD2 to
    # mitochondria more efficiently. Val16 (T) = less efficient import.
    v = _dose(snps, "rs4880", "T", "C")   # T = Val16 (reported orientation)
    if v is None:
        return None
    if v >= 1:
        result = ("SOD2 Val16 allele — mitochondrial superoxide dismutase is "
                  "imported into mitochondria slightly less efficiently, a small "
                  "reduction in front-line defence against the superoxide burst "
                  "that particulate/oxidative exposure triggers.")
        impact = "reduced"
    else:
        result = ("SOD2 Ala16 genotype — efficient mitochondrial targeting of "
                  "this key antioxidant enzyme.")
        impact = "typical"
    action = ("Works downstream with GPX1 and catalase; supported by manganese "
              "(the SOD2 cofactor) sufficiency and an antioxidant-rich diet. See "
              "also the Wellness oxidative-stress trait.")
    return _finding(CAT_ANTIOX, "SOD2 Mitochondrial Superoxide Defense",
                    "SOD2", "rs4880", _gt(snps, "rs4880"), result, action,
                    "moderate", impact, f"rs4880 Val16 dose {v}")


def _gpx1(snps):
    # rs1050450 P198L: Leu198 (T) = reduced glutathione-peroxidase activity;
    # selenium-dependent enzyme that clears peroxides.
    gpx1 = _dose(snps, "rs1050450", "T", "C")   # T = Leu198 (reduced)
    if gpx1 is None:
        return None
    if gpx1 >= 1:
        result = ("GPX1 Leu198 allele — reduced glutathione-peroxidase activity. "
                  "GPX1 uses selenium to neutralise the hydrogen peroxide and "
                  "lipid peroxides generated by smoke/particulate exposure.")
        impact = "reduced"
    else:
        result = ("GPX1 Pro198 genotype — typical glutathione-peroxidase "
                  "activity.")
        impact = "typical"
    action = ("This enzyme is selenium-dependent: ensure adequate selenium "
              "(1–2 Brazil nuts/day is plenty — do not megadose). Pairs with the "
              "glutathione-support measures above.")
    return _finding(CAT_ANTIOX, "GPX1 Peroxide Clearance (selenium)", "GPX1",
                    "rs1050450", _gt(snps, "rs1050450"), result, action,
                    "moderate", impact, f"rs1050450 Leu198 dose {gpx1}")


def _catalase(snps):
    # rs1001179 -262 C>T: T = lower catalase expression / H2O2 clearance.
    t = _dose(snps, "rs1001179", "T", "C")
    if t is None:
        return None
    if t >= 1:
        result = ("CAT -262 T allele — lower catalase expression, so hydrogen "
                  "peroxide (a major reactive species from particulate exposure) "
                  "is cleared a little more slowly.")
        impact = "reduced"
    else:
        result = ("CAT -262 CC — typical catalase expression and H2O2 "
                  "clearance.")
        impact = "typical"
    action = ("Reinforces the antioxidant-rich-diet rationale; catalase works "
              "alongside GPX1 and SOD2 as the peroxide-clearing team.")
    return _finding(CAT_ANTIOX, "Catalase H2O2 Clearance", "CAT",
                    "rs1001179", _gt(snps, "rs1001179"), result, action,
                    "moderate", impact, f"rs1001179 T dose {t}")


def _hmox1(snps):
    # rs2071746 (-413 A>T promoter): modulates HMOX1 (heme-oxygenase-1), a
    # cytoprotective enzyme strongly induced by particulate/oxidative stress.
    t = _dose(snps, "rs2071746", "T", "A")
    if t is None:
        return None
    result = ("Heme-oxygenase-1 (HMOX1) is one of the most protective genes "
              "against particulate-matter and oxidative lung injury. This "
              "promoter variant modestly influences how strongly it is induced.")
    action = ("General antioxidant support applies; HMOX1 is also induced by "
              "exercise and by NRF2 activators such as sulforaphane.")
    return _finding(CAT_ANTIOX, "HMOX1 Cytoprotective Response", "HMOX1",
                    "rs2071746", _gt(snps, "rs2071746"), result, action,
                    "low", "intermediate", f"rs2071746 T dose {t}")


# ─── Heavy-metal handling ──────────────────────────────────────────────────────

def _alad_lead(snps):
    # rs1800435 K59N (ALAD2): alters lead binding in blood → higher blood-lead
    # retention for a given exposure in some studies.
    d = _dose(snps, "rs1800435", "G", "C")   # G = Asn59 (ALAD2)
    if d is None:
        return None
    if d >= 1:
        result = ("Carries the ALAD2 (Asn59) allele. ALAD is the enzyme lead "
                  "binds to in blood; the ALAD2 variant changes lead binding and "
                  "in several cohorts is linked to higher blood-lead levels for a "
                  "given exposure. Research-grade signal, not a lead test.")
        impact = "reduced"
    else:
        result = ("Common ALAD1 genotype — typical lead-binding kinetics.")
        impact = "typical"
    action = ("Minimise lead exposure (old paint/pipes, some imported ceramics, "
              "shooting ranges). Adequate iron, calcium and zinc reduce "
              "gut absorption of lead. If exposure is plausible, a blood-lead "
              "test is the actual measurement.")
    return _finding(CAT_METAL, "ALAD Lead Handling (ALAD2)", "ALAD",
                    "rs1800435", _gt(snps, "rs1800435"), result, action,
                    "low", impact, f"rs1800435 Asn59 dose {d}")


def _as3mt_arsenic(snps):
    # rs11191439 (AS3MT Met287Thr): influences arsenic methylation efficiency
    # and the mono-/di-methyl-arsenic ratio (excretion).
    d = _dose(snps, "rs11191439", "C", "T")   # C = Thr287
    if d is None:
        return None
    if d >= 1:
        result = ("Carries an AS3MT variant (Thr287) that shifts arsenic "
                  "methylation. AS3MT methylates inorganic arsenic so it can be "
                  "excreted; some AS3MT genotypes are associated with a less "
                  "efficient methylation profile and higher retained arsenic.")
        impact = "reduced"
    else:
        result = ("Common AS3MT genotype at this marker — typical arsenic-"
                  "methylation profile.")
        impact = "typical"
    action = ("Arsenic exposure is mostly dietary/water: filter private-well "
              "water (test it), and vary grain sources (rice concentrates "
              "arsenic). Folate/B-vitamin sufficiency supports the methylation "
              "that excretes arsenic.")
    return _finding(CAT_METAL, "AS3MT Arsenic Methylation", "AS3MT",
                    "rs11191439", _gt(snps, "rs11191439"), result, action,
                    "low", impact, f"rs11191439 Thr287 dose {d}")


def _pon1(snps):
    # rs662 Q192R: paraoxonase-1, detoxifies organophosphates and hydrolyses
    # oxidised lipids; Arg192 (G) vs Gln192 (A) shift substrate specificity.
    g = _dose(snps, "rs662", "G", "A")   # G = Arg192
    if g is None:
        return None
    result = ("Paraoxonase-1 (PON1) both detoxifies organophosphate pesticides "
              "and clears oxidised lipids that accumulate under smoke/oxidative "
              "stress. "
              + ("Your Arg192 allele hydrolyses some organophosphates faster but "
                 "certain oxidised lipids more slowly — a trade-off rather than "
                 "'better' or 'worse'."
                 if g >= 1 else
                 "Your Gln192 genotype is the common form."))
    action = ("Minimise organophosphate-pesticide exposure (wash produce, "
              "caution with garden/agricultural chemicals). PON1 activity is "
              "supported by a Mediterranean-style, polyphenol-rich diet.")
    return _finding(CAT_METAL, "PON1 Organophosphate / Oxidised-Lipid Defense",
                    "PON1", "rs662", _gt(snps, "rs662"), result, action,
                    "low", "intermediate" if g >= 1 else "typical",
                    f"rs662 Arg192 dose {g}")


def _mt_metals(snps):
    # Metallothionein MT1A/MT2A — cadmium/lead buffering. Reused from
    # metal_oxidative through the exposure lens.
    mt1 = _dose(snps, "rs8052394", "G", "A")
    mt2 = _dose(snps, "rs28366003", "G", "A")
    if mt1 is None and mt2 is None:
        return None
    score = (mt1 or 0) + (mt2 or 0)
    parts = [p for p in (
        f"MT1A rs8052394 dose {mt1}" if mt1 is not None else "",
        f"MT2A rs28366003 dose {mt2}" if mt2 is not None else "") if p]
    if score >= 2:
        result = ("Metallothionein (MT1A/MT2A) variants that some cohorts link "
                  "to lower cadmium/lead buffering capacity. Cadmium is notably "
                  "concentrated in tobacco smoke. Research-grade only.")
        impact = "reduced"
    else:
        result = ("Metallothionein genotype not associated with reduced metal "
                  "buffering at the variants tested.")
        impact = "typical"
    action = ("Zinc induces metallothionein — ensure adequate (not excessive) "
              "dietary zinc. Avoiding tobacco smoke is the single biggest lever "
              "on cadmium body-burden.")
    return _finding(CAT_METAL, "Metallothionein Cadmium/Lead Buffering",
                    "MT1A/MT2A", "rs8052394", _gt(snps, "rs8052394"), result,
                    action, "low", impact, "; ".join(parts))


# ─── Composite smoke-resilience score ─────────────────────────────────────────

def _smoke_resilience(findings: list[dict]) -> dict:
    """Aggregate the Phase I / Phase II / antioxidant findings into a single
    interpretable susceptibility index. The genotype that matters is FAST Phase I
    (more activation) combined with SLOW Phase II / weak antioxidant defence
    (poor clearance) — that mismatch, not any single gene, is what raises
    smoke-related risk."""
    activation = sum(1 for f in findings
                     if f["category"] == CAT_ACTIVATION and f["impact"] == "higher-load")
    clearance_deficit = sum(1 for f in findings
                            if f["category"] == CAT_CONJUGATION
                            and f["impact"] in ("reduced", "reduced-clearance"))
    antiox_deficit = sum(1 for f in findings
                         if f["category"] == CAT_ANTIOX and f["impact"] == "reduced")

    # Weighted: a clearance/antioxidant deficit counts double when Phase I is
    # fast (the dangerous "activate-but-don't-clear" combination).
    mismatch_bonus = 1 if (activation >= 1 and (clearance_deficit + antiox_deficit) >= 2) else 0
    raw = activation + clearance_deficit + antiox_deficit + mismatch_bonus

    if raw >= 5:
        tier, color, headline = ("Higher susceptibility", "#f85149",
            "Your genotype leans toward activating combustion toxicants faster "
            "than you clear them — the profile that benefits most from serious "
            "smoke-exposure precautions.")
    elif raw >= 2:
        tier, color, headline = ("Typical / mixed", "#d29922",
            "A mixed profile — some steps favourable, some less so. Standard "
            "smoke-exposure precautions apply, with the personalised emphases "
            "below.")
    else:
        tier, color, headline = ("Lower susceptibility", "#3fb950",
            "Your genotype shows no strong activate-but-don't-clear mismatch at "
            "the markers tested. Good baseline resilience — the behavioural "
            "precautions still matter during heavy smoke.")

    return {
        "tier": tier, "color": color, "headline": headline, "raw_score": raw,
        "activation_hits": activation,
        "clearance_deficit_hits": clearance_deficit,
        "antioxidant_deficit_hits": antiox_deficit,
        "activate_but_dont_clear": bool(mismatch_bonus),
    }


# ─── Personalised protocol ─────────────────────────────────────────────────────

def _build_protocol(findings: list[dict], score: dict) -> dict:
    """Turn the genotype picture into a concrete, tiered action plan. Behavioural
    measures apply to everyone; the nutrient emphases scale with the specific
    deficits found."""
    by_gene = {f["gene"]: f for f in findings}

    def _has_reduced(*genes):
        return any(by_gene.get(g, {}).get("impact") in
                   ("reduced", "reduced-clearance") for g in genes)

    nutrition: list[dict] = []

    # Sulforaphane / NRF2 — the single highest-leverage food lever, scaled up
    # when the NRF2/Phase-II axis is weak.
    nrf2_weak = _has_reduced("NFE2L2", "NQO1", "GSTM1/GSTT1/GSTP1")
    nutrition.append({
        "item": "Sulforaphane (broccoli sprouts / cruciferous vegetables)",
        "emphasis": "high" if nrf2_weak else "standard",
        "detail": (
            "Broccoli sprouts are the richest dietary source; mature broccoli, "
            "kale, cabbage and Brussels sprouts also count. Sulforaphane is the "
            "most potent dietary activator of NRF2, which up-regulates the whole "
            "Phase II / antioxidant battery (NQO1, GSTs, HMOX1).") +
            (" Because your NRF2 / Phase-II genotype runs on the weaker side, "
             "this is your top nutritional priority — aim for a daily serving, "
             "and note that a brief hold of chopped raw cruciferous before "
             "cooking preserves the active enzyme." if nrf2_weak else
             " A few servings per week is a sensible baseline."),
    })

    # Glutathione support
    gst_weak = _has_reduced("GSTM1/GSTT1/GSTP1", "EPHX1")
    nutrition.append({
        "item": "Glutathione support (NAC, adequate protein, allium vegetables)",
        "emphasis": "high" if gst_weak else "standard",
        "detail": (
            "Glutathione is the currency Phase II spends to conjugate smoke "
            "toxicants. Support it with sulphur-rich foods (garlic, onions, "
            "eggs) and adequate protein; N-acetylcysteine (NAC) is a direct "
            "cysteine donor worth discussing with a clinician if smoke exposure "
            "is sustained.") +
            (" Your glutathione-transferase / epoxide-clearance genotype makes "
             "this especially relevant." if gst_weak else ""),
    })

    # Selenium for GPX1
    if _has_reduced("GPX1"):
        nutrition.append({
            "item": "Selenium (1–2 Brazil nuts/day)",
            "emphasis": "high",
            "detail": ("Your GPX1 genotype reduces glutathione-peroxidase "
                       "activity, and that enzyme is selenium-dependent. One to "
                       "two Brazil nuts daily is sufficient — selenium is toxic "
                       "in excess, so do not megadose or stack supplements."),
        })

    # Manganese cofactor for SOD2
    if _has_reduced("SOD2"):
        nutrition.append({
            "item": "Manganese-containing whole foods (nuts, wholegrains, tea)",
            "emphasis": "standard",
            "detail": ("SOD2 uses manganese as its cofactor; your Val16 genotype "
                       "makes efficient SOD2 function slightly more valuable. "
                       "Food sources are ample — supplementation is not needed."),
        })

    nutrition.append({
        "item": "Antioxidant-rich, polyphenol-heavy diet + hydration",
        "emphasis": "standard",
        "detail": ("Colourful vegetables/fruit, green tea, olive oil and berries "
                   "supply the broad antioxidant network (vitamins C/E, "
                   "polyphenols) that backs up your enzymatic defences. Stay well "
                   "hydrated during smoke events to support clearance."),
    })

    behavioural = [
        {"item": "Track the air quality index (AQI) daily during smoke events",
         "detail": ("Use AirNow.gov or the AirNow app (EPA) for real-time PM2.5 "
                    "and AQI where you are. Wildfire smoke drifting down from "
                    "Canadian fires has repeatedly pushed Michigan — including "
                    "the Upper Peninsula — into the Unhealthy (AQI >150) range in "
                    "recent summers. Let the number drive your day.")},
        {"item": "Seal and filter your indoor air when PM2.5 is high",
         "detail": ("Close windows, run a HEPA purifier sized to the room, and "
                    "use MERV-13 filters in HVAC/furnace systems. A DIY box-fan + "
                    "MERV-13 (\"Corsi-Rosenthal\") filter is a cheap, effective "
                    "backup. Indoor PM2.5 can be a small fraction of outdoor "
                    "levels with good filtration.")},
        {"item": "Wear a well-fitted N95/KN95 outdoors during heavy smoke",
         "detail": ("Cloth and surgical masks do little against PM2.5. A snug "
                    "N95/KN95 does. Keep a supply for AQI spikes.")},
        {"item": "Shift or scale back outdoor exertion by AQI",
         "detail": ("Exercise multiplies the volume of air — and particulates — "
                    "you inhale. When AQI >100, move hard workouts indoors "
                    "(filtered air); above ~150, keep outdoor exertion brief. "
                    "Early morning is often (not always) cleaner — check the "
                    "number rather than assuming.")},
        {"item": "Zero tobacco / vape smoke",
         "detail": ("This is non-negotiable given your profile: tobacco smoke "
                    "delivers the same PAHs and aromatic amines as wildfire "
                    "smoke plus cadmium, and stacks directly with your Phase I / "
                    "NAT2 genotype. It is the single largest controllable source "
                    "of every toxicant this section covers.")},
    ]

    metal = [
        {"item": "Reduce heavy-metal exposure routes",
         "detail": ("Lead: old paint/pipes, some imported ceramics/spices, "
                    "shooting ranges. Cadmium: tobacco smoke above all. Arsenic: "
                    "private-well water (test it) and rice — vary grains. Mercury: "
                    "limit high-mercury fish (shark, swordfish, king mackerel).")},
        {"item": "Nutrients that blunt metal uptake",
         "detail": ("Adequate iron, calcium and zinc compete with lead and "
                    "cadmium at the gut, lowering absorption. This matters more "
                    "given your ALAD/metallothionein genotype. Folate/B-vitamins "
                    "support the methylation that excretes arsenic.")},
        {"item": "Measure, don't guess",
         "detail": ("Genotype is not body-burden. If a real exposure is "
                    "plausible, a blood-lead level or a urine heavy-metal panel "
                    "ordered by a clinician is the actual test.")},
    ]

    return {"nutrition": nutrition, "behavioural": behavioural, "metal": metal}


MICHIGAN_CONTEXT = (
    "Why this section is front-and-centre for you: the Upper Peninsula and the "
    "wider Great Lakes region have taken repeated hits of drifting Canadian "
    "wildfire smoke over the last few summers, with PM2.5 spikes that pushed "
    "local air into the Unhealthy range for days at a time. Wildfire smoke is "
    "not generic 'bad air' — it is a concentrated dose of the exact toxicants "
    "this panel scores: fine particulate matter that lodges deep in the lungs, "
    "PAHs and aromatic amines that your Phase I/II enzymes must activate and "
    "clear, and a heavy oxidative-stress load. Your genotype doesn't change "
    "whether you should protect yourself during a smoke event — everyone should "
    "— but it does tell you which levers (NRF2/sulforaphane, glutathione, "
    "selenium) give you the most personal return."
)


# ─── Master analyzer ───────────────────────────────────────────────────────────

def analyze_detox(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        # Phase I
        _cyp1a1, _cyp1a2, _cyp1b1, _ahr,
        # Phase II
        _ephx1, _nat2, _nqo1, _gst_conjugation,
        # Antioxidant
        _nrf2, _sod2, _gpx1, _catalase, _hmox1,
        # Heavy metals
        _alad_lead, _as3mt_arsenic, _pon1, _mt_metals,
    ]
    findings: list[dict] = []
    for a in analyzers:
        try:
            r = a(snps_df)
            if r is not None:
                findings.append(r)
        except Exception:
            continue

    by_category: dict[str, list[dict]] = {}
    for f in findings:
        by_category.setdefault(f["category"], []).append(f)

    score = _smoke_resilience(findings)
    protocol = _build_protocol(findings, score)

    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "findings": findings,
        "by_category": by_category,
        "category_order": [CAT_ACTIVATION, CAT_CONJUGATION, CAT_ANTIOX, CAT_METAL],
        "smoke_resilience": score,
        "protocol": protocol,
        "michigan_context": MICHIGAN_CONTEXT,
    }


# ── Cross-check against the unified SNP registry ──────────────────────────────

def _scan_rsids_referenced() -> list[str]:
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> dict[str, list[str]]:
    """Presence-only audit: every rsID referenced here must be registered."""
    registered, missing = [], []
    for rsid in _scan_rsids_referenced():
        (registered if snp_registry.get(rsid) is not None else missing).append(rsid)
    return {"registered": registered, "missing": missing}
