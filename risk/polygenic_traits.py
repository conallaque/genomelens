"""
Trait Genetics — honest genotype-level trait reporting
======================================================

**A deliberate design decision, stated up front:** this module does NOT invent
polygenic percentiles from a handful of hand-picked SNPs. A real polygenic
score requires thousands of variants with per-variant GWAS effect weights,
applied to a matched ancestry, standardized against a reference-population
distribution. The project's `pgs_catalog.py` does exactly that — but only for
the ~15 *disease* scores whose scoring files are downloaded locally. No
published trait PGS scoring files (chronotype, sleep duration, personality,
educational attainment) are present, so computing a trait percentile here would
be fabrication, and no caveat paragraph rescues a fabricated number.

Instead this module reports **genotype-level trait calls** for traits where a
*single variant* genuinely carries interpretable signal — taste, smell, hair,
eye color, chronotype, earwax, a few sleep loci — with the direction and the
primary citation, and nothing more precise than the evidence supports.

For traits that are irreducibly **polygenic** (height) or **socially fraught
and low-signal** (cognition, personality), it does something different: it
explains *why* no number is given, states the ceiling of what genetics can
honestly say, and points to the single-variant loci that are already reported
in the neurochemistry module. There is a hard structural cap here — no score,
no percentile, no ranking — not a caveat bolted onto a number.

References
----------
Kim 2003 (TAS2R38 PTC bitter); Eriksson 2012 23andMe (OR6A2 cilantro,
asparagus anosmia); Yoshiura 2006 (ABCC11 earwax/body odour); Eriksson 2010
(photic sneeze); Sturm 2008 (HERC2 eye color); Han 2008 / Sulem 2007 (MC1R,
IRF4 pigmentation); Kimura 2009 (EDAR hair); Katzenberg 1998 & Garcia-Rios 2014
(CLOCK chronotype); Allebrandt 2013 (ABCC9 sleep duration); Okbay 2016 &
Lee 2018 (educational-attainment PGS — cited to explain the structural cap,
not to score); Lo 2017 (personality GWAS effect sizes).
"""

from __future__ import annotations

import pandas as pd

CAT_CHRONO = "Chronotype & Sleep"
CAT_TASTE = "Taste, Smell & Diet Perception"
CAT_APPEAR = "Appearance & Physical Traits"


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


def _t(cat, trait, gene, rsid, gt, call, detail, citation, confidence="moderate"):
    return {"category": cat, "trait": trait, "gene": gene, "rsid": rsid,
            "genotype": gt or "—", "call": call, "detail": detail,
            "citation": citation, "confidence": confidence}


# ─── CHRONOTYPE & SLEEP ───────────────────────────────────────────────────────

def _clock_chronotype(df):
    gt = _gt(df, "rs1801260")
    if gt is None:
        return None
    # C allele (minor) → evening preference; T/T → morning-leaning
    n_c = gt.count("C")
    if n_c == 0:
        call = "Morning-leaning tendency (T/T)"
        detail = ("The CLOCK 3111 T/T genotype is associated on average with a "
                  "slightly earlier chronotype and morning preference. Effect is "
                  "small — sleep timing is mostly behavioral and light-driven.")
    elif n_c == 1:
        call = "Intermediate chronotype (T/C)"
        detail = "One CLOCK 3111-C allele; intermediate morning/evening tendency."
    else:
        call = "Evening-leaning tendency (C/C)"
        detail = ("The CLOCK 3111 C/C genotype is associated on average with a "
                  "later chronotype / evening preference and slightly higher "
                  "eveningness in questionnaire studies.")
    return _t(CAT_CHRONO, "Chronotype (morning vs evening)", "CLOCK", "rs1801260",
              gt, call, detail, "Katzenberg 1998; Garcia-Rios 2014")


def _per2_morningness(df):
    gt = _gt(df, "rs934945")
    if gt is None:
        return None
    n_a = gt.count("A")
    call = ("PER2 variant carrier" if n_a else "PER2 common genotype")
    detail = ("PER2 rs934945 has weak associations with morningness and "
              "diurnal-preference measures; effect sizes are small.")
    return _t(CAT_CHRONO, "Circadian PER2 variant", "PER2", "rs934945", gt,
              call, detail, "Carpen 2005", confidence="low")


def _abcc9_sleep_duration(df):
    gt = _gt(df, "rs11046205")
    if gt is None:
        return None
    detail = ("ABCC9 has been linked to habitual sleep *duration* in a "
              "pan-European GWAS — carriers of one allele slept ~30 min less "
              "per night on average. This is a population tendency, not a "
              "prescription; your actual sleep need is best found empirically.")
    return _t(CAT_CHRONO, "Sleep-duration tendency", "ABCC9", "rs11046205", gt,
              f"ABCC9 rs11046205 {gt}", detail, "Allebrandt 2013",
              confidence="low")


def _dec2_short_sleeper(df):
    gt = _gt(df, "rs121912617")
    if gt is None:
        return None
    # The famous natural-short-sleeper mutation — almost always absent.
    if any(a in gt for a in ("A", "T")) and gt not in ("GG", "CC"):
        return _t(CAT_CHRONO, "Natural short-sleeper variant (RARE)", "BHLHE41 (DEC2)",
                  "rs121912617", gt, "Possible short-sleeper allele present",
                  "The DEC2/BHLHE41 P384R mutation lets rare carriers thrive on "
                  "~6 h sleep. Extremely rare — confirm before acting on it.",
                  "He 2009", confidence="low")
    return None


# ─── TASTE, SMELL & DIET PERCEPTION ───────────────────────────────────────────

def _tas2r38_bitter(df):
    gt = _gt(df, "rs713598")
    if gt is None:
        return None
    # G (Pro) = taster; C (Ala) = non-taster
    n_taster = gt.count("G")
    if n_taster == 2:
        call = "Strong bitter taster (PAV/PAV)"
        detail = ("You likely perceive PTC/PROP and cruciferous bitterness "
                  "(raw broccoli, Brussels sprouts, coffee, dark beer) strongly. "
                  "Super-tasters often prefer milder foods.")
    elif n_taster == 1:
        call = "Intermediate bitter taster"
        detail = "Intermediate sensitivity to bitter PTC-class compounds."
    else:
        call = "Bitter non-taster (AVI/AVI)"
        detail = ("You likely perceive PTC-class bitterness weakly — bitter "
                  "vegetables and coffee taste milder to you than to tasters.")
    return _t(CAT_TASTE, "Bitter (PTC) tasting", "TAS2R38", "rs713598", gt,
              call, detail, "Kim 2003", confidence="high")


def _or6a2_cilantro(df):
    gt = _gt(df, "rs72921001")
    if gt is None:
        return None
    n_a = gt.count("A")
    if n_a >= 1:
        call = "Cilantro may taste soapy"
        detail = ("Near the OR6A2 olfactory-receptor gene; A-allele carriers are "
                  "more likely to perceive cilantro/coriander as soapy (it detects "
                  "the aldehydes responsible).")
    else:
        call = "Cilantro likely tastes normal"
        detail = "You likely do not have the OR6A2 soapy-cilantro perception."
    return _t(CAT_TASTE, "Cilantro soapy-taste", "OR6A2 region", "rs72921001",
              gt, call, detail, "Eriksson 2012 (23andMe)", confidence="moderate")


def _abcc11_earwax_odor(df):
    gt = _gt(df, "rs17822931")
    if gt is None:
        return None
    n_t = gt.count("T")   # T (538A) = dry earwax + reduced body odour
    if n_t == 2:
        call = "Dry earwax + low body odour"
        detail = ("ABCC11 T/T → dry, flaky earwax and markedly reduced axillary "
                  "body odour (this variant is near-fixed in East Asians). Also "
                  "predicts low colostrum-type apocrine secretion.")
    elif n_t == 1:
        call = "Wet earwax (carrier of dry allele)"
        detail = "One dry-earwax allele; the wet-earwax allele is dominant, so wet earwax."
    else:
        call = "Wet earwax + typical body odour"
        detail = ("ABCC11 C/C → wet (sticky) earwax and typical apocrine body "
                  "odour — the common European/African genotype.")
    return _t(CAT_TASTE, "Earwax type & body-odour", "ABCC11", "rs17822931", gt,
              call, detail, "Yoshiura 2006", confidence="high")


def _asparagus_anosmia(df):
    gt = _gt(df, "rs4481887")
    if gt is None:
        return None
    n_a = gt.count("A")
    if n_a >= 1:
        call = "Reduced ability to smell asparagus metabolites"
        detail = ("Near OR2M7; A-allele carriers are less able to detect the "
                  "sulphurous asparagus-urine metabolites — you may be an "
                  "'asparagus anosmic'.")
    else:
        call = "Can likely smell asparagus metabolites"
        detail = "You likely detect the characteristic asparagus-urine odour."
    return _t(CAT_TASTE, "Asparagus-odour detection", "OR2M7 region", "rs4481887",
              gt, call, detail, "Eriksson 2012 (23andMe)", confidence="moderate")


def _photic_sneeze(df):
    gt = _gt(df, "rs10427255")
    if gt is None:
        return None
    n_c = gt.count("C")
    if n_c >= 1:
        call = "Photic sneeze reflex likely (ACHOO)"
        detail = ("C-allele carriers are more likely to sneeze on sudden bright-"
                  "light exposure (the 'ACHOO' photic sneeze reflex).")
    else:
        call = "Photic sneeze reflex less likely"
        detail = "You are less likely to have the bright-light sneeze reflex."
    return _t(CAT_TASTE, "Photic sneeze reflex", "intergenic 2q22", "rs10427255",
              gt, call, detail, "Eriksson 2010 (23andMe)", confidence="moderate")


# ─── APPEARANCE & PHYSICAL ────────────────────────────────────────────────────

def _herc2_eye(df):
    gt = _gt(df, "rs12913832")
    if gt is None:
        return None
    n_g = gt.count("G")   # G = blue-eye allele (recessive)
    if n_g == 2:
        call = "Blue / light eyes likely"
        detail = ("HERC2 rs12913832 G/G strongly predicts blue or light eye "
                  "color (down-regulates OCA2 → less iris melanin). ~99% of "
                  "blue-eyed Europeans are G/G.")
    elif n_g == 1:
        call = "Green / hazel / intermediate likely"
        detail = "One blue-eye allele — often green, hazel or intermediate iris color."
    else:
        call = "Brown eyes likely"
        detail = "HERC2 A/A predicts brown eyes (functional OCA2, more iris melanin)."
    return _t(CAT_APPEAR, "Eye color", "HERC2 / OCA2", "rs12913832", gt,
              call, detail, "Sturm 2008", confidence="high")


def _mc1r_red_hair(df):
    g7 = _gt(df, "rs1805007")   # R151C
    g8 = _gt(df, "rs1805008")   # R160W
    if g7 is None and g8 is None:
        return None
    n_r = (g7.count("T") if g7 else 0) + (g8.count("T") if g8 else 0)
    if n_r >= 2:
        call = "Red-hair / very-fair-skin alleles (2+ copies)"
        detail = ("Two or more MC1R red-hair variants → high chance of red or "
                  "strawberry-blond hair, very fair freckled skin, higher UV "
                  "sensitivity and lidocaine/anaesthetic considerations.")
    elif n_r == 1:
        call = "One MC1R red-hair allele (carrier)"
        detail = ("One MC1R 'R' variant — often not red-haired but fair-skinned, "
                  "freckle-prone, and a red-hair carrier. Modestly higher UV "
                  "sensitivity.")
    else:
        call = "No MC1R red-hair alleles tested"
        detail = "No R151C/R160W red-hair variants detected."
    return _t(CAT_APPEAR, "Hair color / skin (MC1R)", "MC1R",
              "rs1805007+rs1805008", f"{g7 or '-'}/{g8 or '-'}", call, detail,
              "Sulem 2007; Han 2008", confidence="moderate")


def _irf4_pigment(df):
    gt = _gt(df, "rs12203592")
    if gt is None:
        return None
    n_t = gt.count("T")
    if n_t >= 1:
        call = "Lighter hair / freckling / photoaging tendency"
        detail = ("IRF4 rs12203592 T-allele is associated with lighter hair in "
                  "childhood, more freckling, sun-sensitivity and photoaging.")
    else:
        call = "Common IRF4 pigment genotype"
        detail = "Common genotype; no strong IRF4 lightening/freckling signal."
    return _t(CAT_APPEAR, "Freckling & photoaging", "IRF4", "rs12203592", gt,
              call, detail, "Han 2008", confidence="moderate")


def _edar_hair(df):
    gt = _gt(df, "rs3827760")
    if gt is None:
        return None
    n_g = gt.count("G")   # G (derived, East-Asian) → thick straight hair
    if n_g >= 1:
        call = "Thicker / straighter hair tendency (EDAR derived allele)"
        detail = ("EDAR 370A carriers tend to have thicker hair shafts, straighter "
                  "hair, shovel-shaped incisors and more eccrine sweat glands "
                  "(near-fixed in East Asians and Native Americans).")
    else:
        call = "Ancestral EDAR (typical European/African hair morphology)"
        detail = "Ancestral EDAR allele — the common European/African hair morphology."
    return _t(CAT_APPEAR, "Hair thickness / morphology", "EDAR", "rs3827760", gt,
              call, detail, "Kimura 2009", confidence="moderate")


# ─── Polygenic / fraught traits — STRUCTURAL CAP (no scores) ──────────────────

def _polygenic_notes() -> list[dict]:
    """These traits are either irreducibly polygenic (height) or low-signal and
    socially fraught (cognition, personality). By design this returns
    *explanations*, never scores, percentiles, or rankings."""
    return [
        {
            "trait": "Height",
            "why_no_number": (
                "Height is ~80% heritable but spread across 10,000+ variants of "
                "tiny individual effect. No consumer-chip subset predicts your "
                "adult height with usable accuracy — a genome-wide height PGS "
                "explains ~40% of variance in ideal research conditions and far "
                "less on a genotyping array. Any single-SNP 'tall/short' call "
                "would be meaningless."),
            "honest_statement": (
                "The honest answer for height is: measure it. Genetics confirms "
                "it's highly heritable, but can't give you a trustworthy number "
                "from this data."),
        },
        {
            "trait": "Cognitive ability / working memory",
            "why_no_number": (
                "Educational-attainment and cognitive PGS are real research tools "
                "(Lee 2018 EA4), but even genome-wide they explain only ~10-12% "
                "of variance in educational attainment and less in measured "
                "cognition — and they are heavily confounded by environment and "
                "population structure. A handful of SNPs explains essentially "
                "nothing. **This module will not output a cognition score, "
                "percentile, or ranking — that is a hard design limit, not a "
                "caveat.**"),
            "honest_statement": (
                "The dopamine/plasticity variants that DO carry individual signal "
                "for cognitive *style* (COMT, BDNF) are reported in the "
                "Neurochemistry section — as style, not as an intelligence "
                "measure. Nothing here is an IQ test."),
        },
        {
            "trait": "Personality (extraversion, neuroticism, openness…)",
            "why_no_number": (
                "Personality GWAS exist, but genome-wide personality PGS explain "
                "only ~2-4% of trait variance — near-useless for an individual "
                "prediction. Single variants (e.g. near CADM2 for extraversion) "
                "carry effectively zero individual signal."),
            "honest_statement": (
                "Genetics cannot tell you your personality. Your own experience, "
                "and a validated questionnaire, vastly out-predict any DNA-based "
                "estimate. No score is given here on purpose."),
        },
    ]


# ─── Master ───────────────────────────────────────────────────────────────────

CATEGORY_ORDER = [CAT_CHRONO, CAT_TASTE, CAT_APPEAR]


def analyze_polygenic_traits(df: pd.DataFrame) -> dict:
    analyzers = [
        _clock_chronotype, _per2_morningness, _abcc9_sleep_duration, _dec2_short_sleeper,
        _tas2r38_bitter, _or6a2_cilantro, _abcc11_earwax_odor, _asparagus_anosmia,
        _photic_sneeze,
        _herc2_eye, _mc1r_red_hair, _irf4_pigment, _edar_hair,
    ]
    findings: list[dict] = []
    for a in analyzers:
        try:
            r = a(df)
        except Exception:
            continue
        if r is not None:
            findings.append(r)

    by_category: dict[str, list[dict]] = {}
    for f in findings:
        by_category.setdefault(f["category"], []).append(f)

    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "findings": findings,
        "by_category": by_category,
        "categories": [c for c in CATEGORY_ORDER if c in by_category],
        "polygenic_notes": _polygenic_notes(),
        "methodology_note": (
            "Genotype-level single-variant trait calls only. No polygenic "
            "percentiles are computed — no trait PGS scoring files are available "
            "locally, and fabricating a percentile from a few SNPs would be "
            "dishonest. Height, cognition and personality are handled with an "
            "explicit no-score explanation rather than a number."),
    }
