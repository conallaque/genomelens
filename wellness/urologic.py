"""
Urologic & Genitourinary Panel
==============================

A dedicated genotype-based screen for urologic conditions that were missing from
the report's other panels — the specific ask from the user's chat session
("your report doesn't cover OAB"). Coverage:

  • Overactive bladder (OAB) & detrusor function — ADRB3, CHRM3
  • Benign prostatic hyperplasia (BPH) & 5-alpha-reductase — SRD5A2
  • Prostate cancer — HOXB13 G84E (gold-standard hereditary), 8q24 (rs1447295,
    rs6983267), MSMB (rs10993994)
  • Kidney stones (nephrolithiasis) — CLDN14, CASR, SLC34A1
  • Testicular germ-cell cancer — KITLG (rs995030), SPRY4 (rs4324715)
  • Hypogonadism / erectile function — SHBG (rs1799941), androgen SRD5A2
  • Bladder-cancer smoking amplifier — NAT2 slow acetylator (re-interpreted
    from detox.py for the urologic lens)

The output follows the same shape as ``metal_oxidative.py`` /
``gut_health.py`` (``{category, trait, gene, rsid, genotype, result, action,
confidence, impact}``), so it plugs into the existing renderer / cross-check
tooling. Every SNP is registered in the unified ``snp_registry``.

Educational only. Consumer-chip calls are hints; a positive HOXB13 G84E on a
consumer chip is worth confirming with clinical-grade sequencing (or a
family-history-driven prostate work-up) before acting on it clinically.
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path

import pandas as pd

import snp_registry

CAT_BLADDER = "Bladder & Detrusor (OAB / continence)"
CAT_PROSTATE = "Prostate — BPH & Cancer Risk"
CAT_STONES = "Kidney Stones (Nephrolithiasis)"
CAT_TESTIS = "Testicular / Reproductive"
CAT_HORMONES = "Androgen & DHT Metabolism"


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


def _dose(snps_df, rsid, risk_allele, ref_allele):
    return snp_registry.risk_dose_from_df(
        snps_df, rsid, risk_allele=risk_allele, ref_allele=ref_allele)


def _finding(category, trait, gene, rsid, genotype, result, action,
             confidence, impact, evidence=""):
    return {
        "category": category, "trait": trait, "gene": gene, "rsid": rsid,
        "genotype": genotype or "—", "result": result, "action": action,
        "confidence": confidence, "impact": impact,
        "evidence": evidence or f"{rsid} genotype {genotype or '—'}",
    }


# ─── Bladder / OAB ────────────────────────────────────────────────────────────

def _adrb3_oab(snps):
    """ADRB3 Trp64Arg (rs4994) — the beta-3 adrenergic receptor is the target of
    mirabegron (an OAB drug). Arg64 carriers have altered receptor signalling
    associated with obesity and, in some cohorts, higher OAB symptom risk and
    modified response to beta-3 agonists."""
    d = _dose(snps, "rs4994", "G", "A")   # G = Arg64
    if d is None:
        return None
    gt = _gt(snps, "rs4994")
    if d >= 1:
        result = ("Carries the ADRB3 Arg64 allele (rs4994 G). The beta-3 "
                  "adrenergic receptor is the drug target of mirabegron for "
                  "overactive bladder; carriers show altered receptor signalling "
                  "and, in some cohorts, higher susceptibility to OAB symptoms "
                  "and modified response to beta-3-agonist therapy.")
        impact = "reduced"
    else:
        result = ("Common ADRB3 Trp64/Trp64 genotype — typical beta-3 "
                  "adrenergic-receptor function at this locus.")
        impact = "typical"
    action = ("Awareness only. If OAB symptoms (urgency, frequency, nocturia) "
              "develop, this variant is worth mentioning to a urologist — "
              "mirabegron response and dosing may differ.")
    return _finding(CAT_BLADDER, "ADRB3 Trp64Arg (OAB / β3 receptor)",
                    "ADRB3", "rs4994", gt, result, action,
                    "moderate", impact)


def _chrm3_oab(snps):
    """CHRM3 (M3 muscarinic receptor) rs2229870 — CHRM3 is the target of
    antimuscarinics (oxybutynin, tolterodine, solifenacin) used for OAB.
    Variation has been associated with urinary urgency in some cohorts."""
    d = _dose(snps, "rs2229870", "T", "C")
    if d is None:
        return None
    gt = _gt(snps, "rs2229870")
    if d >= 1:
        result = ("Carries a CHRM3 (M3 muscarinic receptor) coding variant. "
                  "CHRM3 mediates detrusor contraction and is the target of "
                  "antimuscarinic OAB drugs (oxybutynin, solifenacin, "
                  "tolterodine). Effect on drug response is not clinically "
                  "established.")
        impact = "intermediate"
    else:
        result = ("Common CHRM3 genotype at this locus.")
        impact = "typical"
    action = ("If OAB is a concern, mention this variant to a urologist — "
              "antimuscarinic response and side-effect profile (dry mouth, "
              "cognitive effects) are worth tracking.")
    return _finding(CAT_BLADDER, "CHRM3 M3 Muscarinic Receptor (OAB drugs)",
                    "CHRM3", "rs2229870", gt, result, action, "low", impact)


# ─── Prostate — BPH ───────────────────────────────────────────────────────────

def _srd5a2_bph(snps):
    """SRD5A2 V89L (rs523349) — 5α-reductase type-2 converts testosterone to
    DHT. Leu89 (G) has lower enzyme activity → lower prostate-DHT exposure →
    lower BPH progression in some cohorts, and modified finasteride/dutasteride
    response. Also relevant for androgenetic alopecia."""
    d = _dose(snps, "rs523349", "G", "C")
    if d is None:
        return None
    gt = _gt(snps, "rs523349")
    if d >= 1:
        result = ("Carries the SRD5A2 Leu89 allele (rs523349 G). 5α-reductase "
                  "type-2 is the enzyme that converts testosterone to DHT — the "
                  "androgen that drives BPH progression and pattern hair loss. "
                  "Leu89 lowers activity, and cohort studies link it to slower "
                  "BPH progression and modified response to 5-ARI drugs "
                  "(finasteride, dutasteride).")
        impact = "protective"
    else:
        result = ("Common SRD5A2 Val89 genotype — typical 5α-reductase-2 "
                  "activity and DHT synthesis.")
        impact = "typical"
    action = ("Relevant for BPH, prostate cancer trajectory, and hair-loss drug "
              "response (finasteride, dutasteride). Worth flagging if a urologist "
              "considers 5-ARI therapy.")
    return _finding(CAT_PROSTATE, "SRD5A2 V89L (5α-reductase-2 activity)",
                    "SRD5A2", "rs523349", gt, result, action,
                    "moderate", impact)


def _srd5a2_a49t(snps):
    """SRD5A2 A49T (rs9282858) — a gain-of-function variant that in some studies
    correlates with more aggressive prostate cancer. Rare."""
    d = _dose(snps, "rs9282858", "A", "G")
    if d is None:
        return None
    gt = _gt(snps, "rs9282858")
    if d >= 1:
        result = ("Carries the SRD5A2 Thr49 allele (rs9282858 A) — a gain-of-"
                  "function variant that raises 5α-reductase-2 activity. "
                  "Associated in some cohorts with more aggressive prostate "
                  "cancer. Rare in most populations.")
        impact = "reduced"
    else:
        result = ("Common SRD5A2 Ala49 genotype — typical enzyme activity at "
                  "this residue.")
        impact = "typical"
    action = ("No independent action, but flag to a urologist if prostate "
              "cancer is diagnosed — Grade/aggressiveness monitoring matters.")
    return _finding(CAT_PROSTATE, "SRD5A2 A49T (prostate-cancer aggressiveness)",
                    "SRD5A2", "rs9282858", gt, result, action,
                    "low", impact)


# ─── Prostate — Cancer ────────────────────────────────────────────────────────

def _hoxb13_g84e(snps):
    """HOXB13 G84E (rs138213197 T) — the single most important common hereditary
    prostate-cancer variant known. Even heterozygous carriers have ~3-4× lifetime
    prostate-cancer risk with earlier onset. Rare (<1% of Europeans)."""
    d = _dose(snps, "rs138213197", "T", "C")
    if d is None:
        return None
    gt = _gt(snps, "rs138213197")
    if d >= 1:
        result = ("Carries HOXB13 G84E (rs138213197 T) — the highest-penetrance "
                  "COMMON hereditary prostate-cancer variant. Even a single copy "
                  "raises lifetime prostate-cancer risk ~3–4×, with earlier "
                  "onset and higher chance of aggressive disease. Rare in the "
                  "general population (~1% of Europeans, higher in some "
                  "families).")
        impact = "higher-load"
        conf = "high"
    else:
        result = ("HOXB13 G84E not detected (common genotype). This is the "
                  "single most important hereditary prostate-cancer variant on "
                  "consumer chips; a negative result is genuinely reassuring "
                  "for that specific risk.")
        impact = "typical"
        conf = "high"
    action = ("If POSITIVE: confirm with clinical-grade sequencing, consider "
              "PSA screening starting at 40 (rather than 50-55), and share the "
              "result with first-degree male relatives. Discuss with a urologist "
              "and consider genetic counseling.")
    return _finding(CAT_PROSTATE, "HOXB13 G84E (hereditary prostate cancer)",
                    "HOXB13", "rs138213197", gt, result, action,
                    conf, impact)


def _prostate_8q24(snps):
    """8q24 prostate cancer variants — rs1447295 (A) and rs6983267 (G) are two
    independent risk loci in the 8q24 gene desert. Each contributes a modest
    OR (~1.15–1.25 per allele) but they are among the earliest and best-
    replicated common prostate-cancer variants."""
    d1 = _dose(snps, "rs1447295", "A", "C")   # A = risk
    d2 = _dose(snps, "rs6983267", "G", "T")   # G = risk
    if d1 is None and d2 is None:
        return None
    gt1 = _gt(snps, "rs1447295") or "—"
    gt2 = _gt(snps, "rs6983267") or "—"
    risk = (d1 or 0) + (d2 or 0)
    if risk >= 2:
        result = ("Carries multiple 8q24 prostate-cancer risk alleles "
                  "(rs1447295 A and/or rs6983267 G). This gene desert harbours "
                  "some of the earliest replicated common prostate-cancer risk "
                  "variants; each modest on its own, cumulatively contributing "
                  "to overall risk.")
        impact = "reduced"
    elif risk == 1:
        result = ("Carries one 8q24 prostate-cancer risk allele "
                  f"(rs1447295 {gt1}, rs6983267 {gt2}). Modest per-allele "
                  "effect; interpret alongside family history and PRS.")
        impact = "intermediate"
    else:
        result = ("No 8q24 prostate-cancer risk alleles detected at these "
                  "two markers.")
        impact = "typical"
    action = ("These add small increments to prostate-cancer risk on top of "
              "family history and PSA screening. No independent action, but "
              "combine with HOXB13 and PRS for a fuller risk picture.")
    return _finding(CAT_PROSTATE, "8q24 Prostate-Cancer Risk Loci (rs1447295, rs6983267)",
                    "8q24 region", "rs1447295", f"{gt1} / {gt2}",
                    result, action, "moderate", impact,
                    evidence=f"rs1447295 A-dose {d1}; rs6983267 G-dose {d2}")


def _msmb_prostate(snps):
    """MSMB rs10993994 — β-microseminoprotein promoter variant; T allele reduces
    MSMB expression and raises prostate-cancer risk."""
    d = _dose(snps, "rs10993994", "T", "C")
    if d is None:
        return None
    gt = _gt(snps, "rs10993994")
    if d >= 1:
        result = ("Carries the MSMB rs10993994 T allele. This β-microseminoprotein "
                  "promoter variant reduces MSMB expression and is one of the "
                  "most-replicated common prostate-cancer risk variants "
                  "(per-allele OR ~1.2–1.3).")
        impact = "reduced"
    else:
        result = ("Common MSMB genotype — typical β-microseminoprotein expression.")
        impact = "typical"
    action = ("Small independent risk contributor; combine with HOXB13, 8q24 "
              "and family history for a full picture.")
    return _finding(CAT_PROSTATE, "MSMB rs10993994 (β-microseminoprotein / PC risk)",
                    "MSMB", "rs10993994", gt, result, action,
                    "moderate", impact)


# ─── Kidney stones ────────────────────────────────────────────────────────────

def _cldn14_stones(snps):
    """CLDN14 rs219780 — claudin-14 tight-junction gene in the thick ascending
    limb of the loop of Henle. C allele is the risk allele for calcium
    nephrolithiasis, especially with dietary calcium/sodium load."""
    d = _dose(snps, "rs219780", "C", "T")
    if d is None:
        return None
    gt = _gt(snps, "rs219780")
    if d >= 1:
        result = ("Carries the CLDN14 rs219780 C allele. Claudin-14 controls "
                  "paracellular calcium reabsorption in the thick ascending "
                  "limb; risk-allele carriers have higher urinary calcium "
                  "excretion and 40–60% higher risk of calcium kidney stones, "
                  "particularly under high-sodium/high-calcium diets.")
        impact = "reduced"
    else:
        result = ("Common CLDN14 genotype — typical renal calcium handling at "
                  "this locus.")
        impact = "typical"
    action = ("Prudent lifestyle if C-carrier: adequate hydration (>2.5 L/day), "
              "moderate sodium (<2300 mg), dietary calcium (paradoxically "
              "protective when consumed WITH oxalate-containing foods), and "
              "limit oxalate-heavy foods (spinach, chocolate, nuts) if a stone "
              "has been passed.")
    return _finding(CAT_STONES, "CLDN14 rs219780 (calcium-stone risk)",
                    "CLDN14", "rs219780", gt, result, action,
                    "moderate", impact)


def _slc34a1_stones(snps):
    """SLC34A1 rs4074995 — sodium-phosphate cotransporter in the proximal
    tubule. G allele associates with higher kidney-stone risk and reduced eGFR
    in large GWAS."""
    d = _dose(snps, "rs4074995", "G", "A")
    if d is None:
        return None
    gt = _gt(snps, "rs4074995")
    if d >= 1:
        result = ("Carries the SLC34A1 rs4074995 G allele. This proximal-"
                  "tubule sodium-phosphate cotransporter variant is associated "
                  "with kidney-stone risk and slightly lower eGFR in large "
                  "population studies.")
        impact = "reduced"
    else:
        result = ("Common SLC34A1 genotype at this locus.")
        impact = "typical"
    action = ("Hydration and standard stone-prevention diet if there is a "
              "personal or family history of stones.")
    return _finding(CAT_STONES, "SLC34A1 rs4074995 (phosphate transport, stones)",
                    "SLC34A1", "rs4074995", gt, result, action, "low", impact)


def _casr_stones(snps):
    """CASR rs1042636 (Arg990Gly) — calcium-sensing receptor. G (Gly990) is a
    gain-of-function allele associated with higher urinary calcium and stone
    risk."""
    d = _dose(snps, "rs1042636", "G", "A")
    if d is None:
        return None
    gt = _gt(snps, "rs1042636")
    if d >= 1:
        result = ("Carries the CASR Gly990 allele (rs1042636 G). This gain-of-"
                  "function variant in the calcium-sensing receptor raises "
                  "urinary calcium excretion and is associated with kidney-"
                  "stone risk and, at higher copy numbers, mild hypercalciuria.")
        impact = "reduced"
    else:
        result = ("Common CASR Arg990 genotype — typical calcium-sensing "
                  "receptor activity.")
        impact = "typical"
    action = ("Combine with CLDN14 and dietary calcium/sodium context. Adequate "
              "hydration and moderating sodium/animal-protein intake is the "
              "general lever.")
    return _finding(CAT_STONES, "CASR R990G (calcium-sensing receptor)",
                    "CASR", "rs1042636", gt, result, action,
                    "moderate", impact)


# ─── Testicular / reproductive ────────────────────────────────────────────────

def _kitlg_testicular(snps):
    """KITLG rs995030 — a top testicular germ-cell tumor (TGCT) risk locus.
    G allele confers modestly higher risk (per-allele OR ~1.3-1.4)."""
    d = _dose(snps, "rs995030", "G", "A")
    if d is None:
        return None
    gt = _gt(snps, "rs995030")
    if d >= 1:
        result = ("Carries the KITLG rs995030 G allele — one of the strongest "
                  "common testicular germ-cell cancer risk variants known "
                  "(per-allele OR ~1.3–1.4). Absolute risk still low.")
        impact = "reduced"
    else:
        result = ("Common KITLG rs995030 genotype at this locus.")
        impact = "typical"
    action = ("Practical action: know your baseline — monthly testicular self-"
              "exam and prompt evaluation of any painless lump. Family history "
              "of testicular cancer materially raises risk regardless of this "
              "variant.")
    return _finding(CAT_TESTIS, "KITLG rs995030 (testicular cancer risk)",
                    "KITLG", "rs995030", gt, result, action,
                    "moderate", impact)


def _spry4_testicular(snps):
    """SPRY4 rs4324715 — another TGCT locus."""
    d = _dose(snps, "rs4324715", "C", "T")
    if d is None:
        return None
    gt = _gt(snps, "rs4324715")
    if d >= 1:
        result = ("Carries the SPRY4 rs4324715 C allele — a second replicated "
                  "testicular germ-cell tumor risk variant with a modest per-"
                  "allele effect.")
        impact = "reduced"
    else:
        result = ("Common SPRY4 rs4324715 genotype at this locus.")
        impact = "typical"
    action = ("Same as KITLG — self-exam and prompt lump evaluation.")
    return _finding(CAT_TESTIS, "SPRY4 rs4324715 (testicular cancer risk)",
                    "SPRY4", "rs4324715", gt, result, action,
                    "low", impact)


# ─── Androgen / hormone metabolism ────────────────────────────────────────────

def _shbg_variant(snps):
    """SHBG rs1799941 A allele — associates with higher SHBG levels, hence lower
    bioavailable (free) testosterone at any given total-T."""
    d = _dose(snps, "rs1799941", "A", "G")
    if d is None:
        return None
    gt = _gt(snps, "rs1799941")
    if d >= 1:
        result = ("Carries the SHBG rs1799941 A allele. Associated with higher "
                  "SHBG (sex-hormone-binding globulin) levels — meaning at any "
                  "given TOTAL testosterone, the bioavailable FREE testosterone "
                  "is lower. Relevant if hypogonadism or ED symptoms arise: "
                  "total-T alone can be misleading; free-T or a calculated free-T "
                  "from total-T + SHBG is a better read.")
        impact = "reduced"
    else:
        result = ("Common SHBG rs1799941 genotype — typical SHBG regulation.")
        impact = "typical"
    action = ("If symptoms of low testosterone (fatigue, low libido, ED) "
              "prompt lab testing, always pair a total-T with an SHBG so free-T "
              "can be calculated. This variant makes total-T alone especially "
              "unreliable.")
    return _finding(CAT_HORMONES, "SHBG rs1799941 (free-testosterone bioavailability)",
                    "SHBG", "rs1799941", gt, result, action,
                    "moderate", impact)


def _pkd1_ident(snps):
    """PKD1 rs2072499 — a common polymorphism in the polycystic kidney disease
    1 region. NOT a pathogenic PKD1 mutation (those are rare LOFs) — this common
    SNP is a modest kidney-function modifier in GWAS."""
    d = _dose(snps, "rs2072499", "A", "G")
    if d is None:
        return None
    gt = _gt(snps, "rs2072499")
    if d >= 1:
        result = ("Common PKD1 polymorphism (rs2072499). NOT a pathogenic PKD1 "
                  "mutation — those are rare loss-of-function variants that "
                  "cause autosomal-dominant polycystic kidney disease. This is "
                  "a common SNP with modest effect on kidney-function traits "
                  "in GWAS.")
        impact = "typical"
    else:
        result = ("Common PKD1 rs2072499 genotype — no independent implication.")
        impact = "typical"
    action = ("No independent action. A family history of PKD warrants proper "
              "clinical evaluation (renal ultrasound, PKD1/PKD2 sequencing), "
              "not this SNP.")
    return _finding(CAT_STONES, "PKD1 rs2072499 (common polymorphism)",
                    "PKD1", "rs2072499", gt, result, action, "low", impact)


# ─── Bladder-cancer smoking amplifier (cross-referenced from detox) ───────────

def _nat2_bladder(snps):
    """NAT2 slow acetylator (*5/*6/*7) — the classic gene×environment risk for
    smoking-related bladder cancer. Reused from detox.py for the urologic lens."""
    a5 = _dose(snps, "rs1801280", "C", "T")
    a6 = _dose(snps, "rs1799930", "A", "G")
    a7 = _dose(snps, "rs1799931", "A", "G")
    if a5 is None and a6 is None and a7 is None:
        return None
    slow = (a5 or 0) + (a6 or 0) + (a7 or 0)
    if slow >= 2:
        result = ("Predicted SLOW NAT2 acetylator. NAT2 conjugates the aromatic "
                  "amines in cigarette smoke, which are the class of carcinogens "
                  "most linked to bladder cancer. Slow acetylators exposed to "
                  "tobacco smoke have materially higher bladder-cancer risk — "
                  "this is one of the best-established gene×environment "
                  "interactions in cancer epidemiology.")
        impact = "reduced"
        conf = "moderate"
    elif slow == 1:
        result = ("Intermediate NAT2 acetylator — between fast and slow.")
        impact = "intermediate"
        conf = "moderate"
    else:
        result = ("Predicted FAST/normal NAT2 acetylator — efficient aromatic-"
                  "amine handling.")
        impact = "typical"
        conf = "moderate"
    action = ("If slow: zero tobacco (single largest lever), limit exposure to "
              "hair-dye and rubber-industry aromatic amines, and mention to a "
              "urologist that bladder-cancer screening symptoms (painless "
              "hematuria) deserve a lower threshold to investigate.")
    return _finding(CAT_BLADDER, "NAT2 Acetylator × Bladder-Cancer Risk",
                    "NAT2", "rs1801280", _gt(snps, "rs1801280") or "—",
                    result, action, conf, impact,
                    evidence=f"*5:{a5}, *6:{a6}, *7:{a7}")


# ─── Master analyzer ──────────────────────────────────────────────────────────

CATEGORY_ORDER = [CAT_BLADDER, CAT_PROSTATE, CAT_STONES, CAT_TESTIS, CAT_HORMONES]


def analyze_urologic(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        _adrb3_oab, _chrm3_oab, _nat2_bladder,
        _srd5a2_bph, _srd5a2_a49t, _hoxb13_g84e, _prostate_8q24, _msmb_prostate,
        _cldn14_stones, _slc34a1_stones, _casr_stones, _pkd1_ident,
        _kitlg_testicular, _spry4_testicular,
        _shbg_variant,
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

    n_flagged = sum(1 for f in findings if f["impact"] in
                    ("reduced", "reduced-clearance", "higher-load"))
    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "n_flagged": n_flagged,
        "findings": findings,
        "by_category": by_category,
        "categories": [c for c in CATEGORY_ORDER if c in by_category],
    }


# ── Registry cross-check ─────────────────────────────────────────────────────

def _scan_rsids_referenced() -> list[str]:
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> dict[str, list[str]]:
    registered, missing = [], []
    for rsid in _scan_rsids_referenced():
        (registered if snp_registry.get(rsid) is not None else missing).append(rsid)
    return {"registered": registered, "missing": missing}
