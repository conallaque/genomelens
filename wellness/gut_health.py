"""
Gut Health & Microbiome Genetics Panel
======================================

Host-genetic variants that influence digestion and *shape* the gut microbiome.

IMPORTANT — what this is NOT
----------------------------
A consumer SNP chip genotypes **your** (host) germline DNA. It cannot sequence
the bacteria living in your gut — that is metagenomic (16S / shotgun) data from
a stool sample, which this pipeline never sees. So this panel does **not**
report "your microbiome." It reports human germline variants with replicated
associations to:

  * carbohydrate digestion  — LCT/MCM6 lactase persistence
  * microbiome shaping      — FUT2 secretor status → mucosal fucosylation, a
                              well-replicated *host* determinant of
                              *Bifidobacterium* abundance
  * gluten / coeliac risk   — HLA-DQ2.5 / DQ8 tag SNPs
  * food intolerance        — AOC1 (DAO) histamine-degradation capacity
  * inflammatory-bowel risk — NOD2 (risk) and IL23R (protective) modifiers

These are **predispositions**, not diagnoses or measurements; each listed trait
explains only a fraction of phenotype variance and is heavily modified by diet,
fibre intake, antibiotics, and environment. Microbiome-shaping findings are
tiered low-moderate confidence accordingly.

Several SNPs are reused from existing registry records rather than re-declared:
FUT2 (rs602662 — wellness.py reads it for B12) and the coeliac HLA-DQ tags
(rs2187668 / rs7454108 — nutrition.py reads them for gluten). This panel reads
them for the gut angle and **complements — does not replace** those sections.

Strand handling
---------------
Unlike wellness.py / metal_oxidative.py (which count a literal allele), this
module defers dosage to ``snp_registry.risk_dose_from_df`` so chips that report
a SNP on the - strand are counted correctly. rs4988235 (lactase persistence) is
the classic strand-flip trap — commonly reported as C/T (- strand) rather than
the registry's canonical + strand G/A — and a literal ``count("A")`` would read
it backwards. The registry's ``derived`` allele is the risk/effect allele.

Output shape matches wellness.py / metal_oxidative.py so the renderer can group
by category:  {category, trait, result, action, evidence, confidence}
"""

from __future__ import annotations

import re as _re
from pathlib import Path as _Path

import pandas as pd

from core import snp_registry  # strand-aware dose + import-time audit (see bottom)


def _gt(snps_df: pd.DataFrame, rsid: str) -> str | None:
    """Raw genotype string for chip-gap detection (presence vs absence)."""
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _dose(snps_df: pd.DataFrame, rsid: str) -> int | None:
    """Strand-aware derived(=effect)-allele dose via the unified registry."""
    return snp_registry.risk_dose_from_df(snps_df, rsid)


CAT_CARB = "Carbohydrate Digestion"
CAT_MICRO = "Microbiome Shaping"
CAT_GLUTEN = "Gluten & Coeliac Risk"
CAT_INTOL = "Food Intolerance"
CAT_IBD = "Inflammatory-Bowel Predisposition"


def _chip_gap(category: str, trait: str, rsid: str, gene: str) -> dict:
    """Surface a variant the chip didn't type — never silently return nothing.

    A silent ``return None`` is indistinguishable to the user from "checked and
    typical". This makes the gap explicit (see README contributing notes)."""
    return {
        "category": category, "trait": trait,
        "result": f"Not evaluable — {gene} {rsid} is not genotyped on this file.",
        "action": ("No call possible from this chip. A panel that types this "
                   "position would be needed to evaluate it."),
        "evidence": f"{rsid} absent from input", "confidence": "none",
    }


# ─── Carbohydrate digestion ───────────────────────────────────────────────

def _lactase(snps):
    rsid = "rs4988235"  # derived A = lactase-persistence allele (LCT/MCM6)
    if _gt(snps, rsid) is None:
        return _chip_gap(CAT_CARB, "Lactase Persistence (LCT/MCM6)", rsid, "MCM6")
    dose = _dose(snps, rsid)
    if dose is None:
        return None
    if dose >= 1:
        result = ("Lactase-persistence allele present — lactase production "
                  "typically continues into adulthood, so fresh dairy is usually "
                  "well tolerated.")
        action = ("No restriction indicated on genetic grounds. Secondary "
                  "intolerance can still occur (illness, ageing, gut injury).")
    else:
        result = ("Lactase non-persistence genotype (no persistence allele) — "
                  "lactase activity declines after weaning in most carriers, so "
                  "fresh milk / ice-cream may cause bloating, gas, or diarrhoea.")
        action = ("Favour lactose-free milk, hard aged cheeses (parmesan, "
                  "cheddar) and fermented dairy (yoghurt, kefir); lactase enzyme "
                  "tablets help. Penetrance varies — many non-persisters tolerate "
                  "small amounts.")
    return {"category": CAT_CARB, "trait": "Lactase Persistence (LCT/MCM6)",
            "result": result, "action": action,
            "evidence": f"{rsid} persistence-allele dose: {dose}", "confidence": "high"}


# ─── Microbiome shaping ────────────────────────────────────────────────────

def _fut2_secretor(snps):
    rsid = "rs602662"  # derived A = non-secretor-associated allele
    if _gt(snps, rsid) is None:
        return _chip_gap(CAT_MICRO, "FUT2 Secretor Status", rsid, "FUT2")
    dose = _dose(snps, rsid)
    if dose is None:
        return None
    if dose >= 2:
        result = ("Likely non-secretor (FUT2). Without secreted ABO/H antigens "
                  "in gut mucus, cohorts consistently show lower "
                  "Bifidobacterium abundance and altered mucosal microbiota. "
                  "Non-secretors also resist most norovirus strains and may "
                  "absorb B12 less efficiently (see Wellness).")
        action = ("Consider Bifidobacterium-containing fermented foods / "
                  "probiotics and prebiotic fibre to support bifidobacteria; "
                  "monitor B12 status. This is a host-genetic tendency, not a "
                  "microbiome measurement.")
    else:
        result = ("Likely secretor (FUT2) — secretes ABO/H antigens into gut "
                  "mucus (the more common state), associated with typical "
                  "Bifidobacterium-bearing mucosal microbiota.")
        action = ("No specific action. Secretor status is one of many host "
                  "factors shaping the microbiome; daily diet and fibre intake "
                  "dominate.")
    return {"category": CAT_MICRO, "trait": "FUT2 Secretor Status",
            "result": result, "action": action,
            "evidence": f"{rsid} non-secretor-allele dose: {dose}", "confidence": "moderate"}


# ─── Gluten / coeliac ──────────────────────────────────────────────────────

def _coeliac(snps):
    have_25 = _gt(snps, "rs2187668") is not None
    have_8 = _gt(snps, "rs7454108") is not None
    if not have_25 and not have_8:
        return _chip_gap(CAT_GLUTEN, "Coeliac HLA-DQ2.5 / DQ8 Risk",
                         "rs2187668 / rs7454108", "HLA-DQA1/DQB1")
    dq25 = _dose(snps, "rs2187668")  # derived T tags DQ2.5
    dq8 = _dose(snps, "rs7454108")   # derived C tags DQ8
    carries = (dq25 or 0) >= 1 or (dq8 or 0) >= 1
    tags = []
    if dq25 is not None:
        tags.append(f"DQ2.5(rs2187668) dose {dq25}")
    if dq8 is not None:
        tags.append(f"DQ8(rs7454108) dose {dq8}")
    if carries:
        result = ("Carries a coeliac-permissive HLA haplotype tag (DQ2.5 and/or "
                  "DQ8). ~30-40% of the general population carry one and most "
                  "never develop coeliac disease — the haplotype is necessary "
                  "but far from sufficient. Tag-SNP method, not full HLA typing.")
        action = ("If you have GI / iron-deficiency / skin symptoms, ask a "
                  "physician for tTG-IgA serology BEFORE removing gluten "
                  "(going gluten-free first invalidates the test). Do not "
                  "self-diagnose from this tag.")
    else:
        result = ("Neither the DQ2.5 nor DQ8 tag SNP detected. >99% of coeliac "
                  "patients carry one of these haplotypes, so classic coeliac "
                  "disease is genetically very unlikely.")
        action = ("Coeliac disease essentially ruled out on genetic grounds; "
                  "gluten avoidance is not warranted for autoimmune reasons. "
                  "Non-coeliac gluten/wheat sensitivity is a separate, non-HLA "
                  "entity.")
    return {"category": CAT_GLUTEN, "trait": "Coeliac HLA-DQ2.5 / DQ8 Risk",
            "result": result, "action": action,
            "evidence": "; ".join(tags), "confidence": "moderate"}


# ─── Food intolerance ──────────────────────────────────────────────────────

def _histamine_dao(snps):
    rsid = "rs10156191"  # derived T (Thr16Met) = lower DAO activity
    if _gt(snps, rsid) is None:
        return _chip_gap(CAT_INTOL, "Histamine Degradation (AOC1/DAO)", rsid, "AOC1")
    dose = _dose(snps, rsid)
    if dose is None:
        return None
    if dose >= 1:
        result = ("Carries the AOC1 (DAO) Thr16Met allele associated in some "
                  "cohorts with lower serum diamine-oxidase activity — the "
                  "enzyme that clears dietary histamine. May predispose to "
                  "histamine-intolerance symptoms (flushing, headache, hives, "
                  "loose stool) from histamine-rich foods. Research-grade.")
        action = ("If symptoms fit, trial a temporary low-histamine diet (limit "
                  "aged cheese, cured meat, wine, fermented foods) with a "
                  "clinician/dietitian. DAO genetics explain only part of "
                  "histamine tolerance — gut health and some medications also "
                  "lower DAO.")
    else:
        result = ("No reduced-DAO allele at this SNP — typical genetic "
                  "histamine-degradation capacity (this single variant only).")
        action = "No histamine-specific dietary restriction indicated genetically."
    return {"category": CAT_INTOL, "trait": "Histamine Degradation (AOC1/DAO)",
            "result": result, "action": action,
            "evidence": f"{rsid} reduced-DAO-allele dose: {dose}", "confidence": "low"}


# ─── Inflammatory-bowel predisposition ─────────────────────────────────────

def _nod2_crohn(snps):
    rsid = "rs2066844"  # derived T = R702W Crohn's risk allele
    if _gt(snps, rsid) is None:
        return _chip_gap(CAT_IBD, "NOD2 R702W (Crohn's risk)", rsid, "NOD2")
    dose = _dose(snps, rsid)
    if dose is None:
        return None
    if dose >= 1:
        zyg = "homozygous" if dose >= 2 else "heterozygous"
        result = (f"Carries NOD2 R702W ({zyg}), one of three classic "
                  "Crohn's-disease risk alleles. Heterozygotes have ~2-4x "
                  "relative risk and homozygotes/compound heterozygotes more — "
                  "but absolute lifetime risk stays low; most carriers never "
                  "develop Crohn's.")
        action = ("Awareness only — no validated gene-based prevention. Not "
                  "smoking is the strongest modifiable Crohn's factor. Report "
                  "persistent diarrhoea, abdominal pain, or rectal bleeding to a "
                  "physician. The other two NOD2 variants (G908R, 1007fs) are "
                  "not evaluated here.")
    else:
        result = ("No NOD2 R702W risk allele at this SNP — typical (the G908R "
                  "and 1007fs NOD2 variants are not tested here).")
        action = "No NOD2-specific action indicated from this variant."
    return {"category": CAT_IBD, "trait": "NOD2 R702W (Crohn's risk)",
            "result": result, "action": action,
            "evidence": f"{rsid} R702W-allele dose: {dose}", "confidence": "moderate"}


def _il23r_protective(snps):
    rsid = "rs11209026"  # derived A = R381Q protective allele
    if _gt(snps, rsid) is None:
        return _chip_gap(CAT_IBD, "IL23R R381Q (protective)", rsid, "IL23R")
    dose = _dose(snps, rsid)
    if dose is None:
        return None
    if dose >= 1:
        result = ("Carries IL23R R381Q (Gln381) — a well-replicated *protective* "
                  "variant that lowers risk of Crohn's disease, ulcerative "
                  "colitis, ankylosing spondylitis and psoriasis by dampening "
                  "IL-23 signalling. Present in only ~5-7% of Europeans.")
        action = ("Favourable finding; no action needed. Does not guarantee "
                  "protection — IBD is highly polygenic.")
    else:
        result = ("Does not carry the IL23R R381Q protective allele — the common "
                  "state (~93% of people); confers no protection or extra risk by "
                  "itself.")
        action = "No action; this is the common genotype."
    return {"category": CAT_IBD, "trait": "IL23R R381Q (protective)",
            "result": result, "action": action,
            "evidence": f"{rsid} protective-allele dose: {dose}", "confidence": "moderate"}


# ─── Master analyzer ──────────────────────────────────────────────────────

def analyze_gut_health(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        # Carbohydrate digestion
        _lactase,
        # Microbiome shaping
        _fut2_secretor,
        # Gluten / coeliac
        _coeliac,
        # Food intolerance
        _histamine_dao,
        # Inflammatory-bowel predisposition
        _nod2_crohn, _il23r_protective,
    ]
    predictions: list[dict] = []
    for a in analyzers:
        try:
            r = a(snps_df)
            if r is not None:
                predictions.append(r)
        except Exception:
            continue

    by_category: dict[str, list[dict]] = {}
    for p in predictions:
        by_category.setdefault(p["category"], []).append(p)

    return {
        "predictions": predictions,
        "by_category": by_category,
        "n_predictions": len(predictions),
        "categories": list(by_category.keys()),
    }


# ── Cross-check against the unified SNP registry ──────────────────────────

def _scan_rsids_referenced() -> list[str]:
    """Every standalone rsID literal referenced in this module must be
    registered (combined-string literals like ``"rsA / rsB"`` are intentionally
    not matched by this regex; their components appear separately)."""
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> dict[str, list[str]]:
    """Presence-only audit: every rsID referenced here must be in the registry.
    Returns ``{"registered": [...], "missing": [...]}``."""
    registered, missing = [], []
    for rsid in _scan_rsids_referenced():
        (registered if snp_registry.get(rsid) is not None else missing).append(rsid)
    return {"registered": registered, "missing": missing}
