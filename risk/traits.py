"""
Trait Predictions
-----------------

Concrete phenotype calls from genotype data — not disease risk, but
observable traits. Each prediction includes the genotype evidence, the
prediction, and a confidence indicator.

Traits covered:
  * Lactose persistence (LCT)
  * Alcohol flush (ALDH2)
  * Caffeine metabolism speed (CYP1A2)
  * Bitter taste perception (TAS2R38 haplotype)
  * Earwax type (ABCC11)
  * Eye color (multi-SNP — HERC2, OCA2, TYR, SLC24A4)
  * Hair color tendency (MC1R + others)
  * Chronotype (PER3 / CLOCK)
  * Short-sleeper allele (BHLHE41)
  * Muscle fiber composition (ACTN3)
  * Vitamin D synthesis efficiency (DBP / VDR)
  * Caffeine-induced anxiety (ADORA2A)
"""

import re as _re
from pathlib import Path as _Path

import pandas as pd

from core import snp_registry  # V8 cross-check; see audit_against_registry below


def _gt(snps_df: pd.DataFrame, rsid: str) -> str | None:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> int | None:
    gt = _gt(snps_df, rsid)
    if gt is None or len(gt) != 2:
        return None
    return gt.count(allele.upper())


# ─── Individual trait analyzers ───────────────────────────────────────────────

def _trait_lactose(snps_df) -> dict:
    # T allele = lactase persistence (continued lactase production)
    g = _gt(snps_df, "rs4988235") or _gt(snps_df, "rs182549")
    if g is None:
        return {"trait": "Lactose Persistence (LCT)", "result": "Not tested",
                "evidence": "rs4988235 / rs182549 not called", "confidence": "n/a"}
    if "T" in g:
        return {"trait": "Lactose Persistence (LCT)",
                "result": "Likely lactose tolerant (lactase persists into adulthood)",
                "evidence": f"rs4988235 genotype: {g} (T allele present)",
                "confidence": "high"}
    return {"trait": "Lactose Persistence (LCT)",
            "result": "Likely lactose intolerant (lactase non-persistence after weaning)",
            "evidence": f"rs4988235 genotype: {g} (no T allele)",
            "confidence": "high",
            "detail": (
                "Symptoms (bloating, diarrhea) often emerge in late childhood/teenage years. "
                "Cultured dairy (yogurt, kefir) and aged cheeses are usually well tolerated. "
                "Lactase enzyme supplements or lactose-free milk as needed. Most of East Asia, "
                "Sub-Saharan Africa, and many indigenous populations are non-persistent."
            )}


def _trait_alcohol_flush(snps_df) -> dict:
    g = _gt(snps_df, "rs671")  # A = Lys487 = deficient
    if g is None:
        return {"trait": "Alcohol Flush (ALDH2)", "result": "Not tested",
                "evidence": "rs671 not called", "confidence": "n/a"}
    if g == "AA":
        return {"trait": "Alcohol Flush (ALDH2)",
                "result": "Severe alcohol intolerance — homozygous ALDH2 deficiency",
                "evidence": "rs671 genotype: AA (Lys/Lys)",
                "confidence": "high",
                "detail": (
                    "Essentially cannot metabolize acetaldehyde. Alcohol consumption "
                    "causes intense flushing, tachycardia, nausea. Strongly elevated "
                    "esophageal cancer risk with any drinking. Best abstained."
                )}
    if "A" in g:
        return {"trait": "Alcohol Flush (ALDH2)",
                "result": "Asian flush — partial ALDH2 deficiency",
                "evidence": f"rs671 genotype: {g} (heterozygous)",
                "confidence": "high",
                "detail": (
                    "Flushing with alcohol; reduced tolerance. Drinkers have ~5× "
                    "esophageal squamous cell carcinoma risk. Minimize alcohol."
                )}
    return {"trait": "Alcohol Flush (ALDH2)",
            "result": "Normal ALDH2 activity — no flushing predicted",
            "evidence": f"rs671 genotype: {g}",
            "confidence": "high"}


def _trait_caffeine_speed(snps_df) -> dict:
    g = _gt(snps_df, "rs762551")  # CYP1A2 *1F — C is slow
    if g is None:
        return {"trait": "Caffeine Metabolism Speed", "result": "Not tested",
                "evidence": "rs762551 not called", "confidence": "n/a"}
    if g == "CC":
        return {"trait": "Caffeine Metabolism Speed",
                "result": "Slow caffeine metabolizer (CYP1A2 *1F/*1F)",
                "evidence": "rs762551 genotype: CC",
                "confidence": "high",
                "detail": (
                    "Caffeine half-life is ~50% longer. Coffee in the afternoon will "
                    "disrupt sleep. Some studies link slow metabolism + high coffee "
                    "intake to elevated MI risk. Recommend <200 mg/day, none after noon."
                )}
    if g in ("AC", "CA"):
        return {"trait": "Caffeine Metabolism Speed",
                "result": "Intermediate caffeine metabolizer",
                "evidence": "rs762551 genotype: AC",
                "confidence": "high"}
    return {"trait": "Caffeine Metabolism Speed",
            "result": "Fast caffeine metabolizer (CYP1A2 *1F/*1A or *1A/*1A)",
            "evidence": f"rs762551 genotype: {g}",
            "confidence": "high",
            "detail": "Caffeine cleared quickly. Higher intake tolerated (still cap ≤400 mg/day)."}


def _trait_bitter_taste(snps_df) -> dict:
    # TAS2R38: PAV haplotype (rs713598 C, rs1726866 G, rs10246939 T) = strong taster
    # AVI haplotype (G, A, C) = non-taster
    a = _gt(snps_df, "rs713598")
    b = _gt(snps_df, "rs1726866")
    c = _gt(snps_df, "rs10246939")
    if None in (a, b, c):
        return {"trait": "Bitter Taste (PTC / PROP)", "result": "Not fully tested",
                "evidence": "Need all three of rs713598, rs1726866, rs10246939",
                "confidence": "n/a"}
    pav_count = (a.count("C") + b.count("G") + c.count("T")) // 3
    if pav_count >= 2:
        return {"trait": "Bitter Taste (PTC / PROP)",
                "result": "Strong taster (PAV/PAV)",
                "evidence": "TAS2R38 haplotype: PAV/PAV",
                "confidence": "high",
                "detail": (
                    "Bitter compounds (PTC, PROP, glucosinolates in cruciferous "
                    "vegetables) taste intensely bitter. Often associated with vegetable "
                    "aversion in childhood, lower alcohol/coffee preference. Light "
                    "steaming + olive oil/lemon, roasting (caramelizes bitterness), "
                    "or blending into smoothies makes greens more palatable."
                )}
    if pav_count == 1:
        return {"trait": "Bitter Taste (PTC / PROP)",
                "result": "Intermediate taster (PAV/AVI)",
                "evidence": "TAS2R38 haplotype: heterozygous PAV/AVI",
                "confidence": "moderate"}
    return {"trait": "Bitter Taste (PTC / PROP)",
            "result": "Non-taster (AVI/AVI)",
            "evidence": "TAS2R38 haplotype: AVI/AVI",
            "confidence": "high",
            "detail": "Cannot taste PTC/PROP; cruciferous vegetables and bitter foods are tolerated easily."}


def _trait_earwax(snps_df) -> dict:
    g = _gt(snps_df, "rs17822931")  # T/T = dry earwax
    if g is None:
        return {"trait": "Earwax Type (ABCC11)", "result": "Not tested",
                "evidence": "rs17822931 not called", "confidence": "n/a"}
    if g == "TT":
        return {"trait": "Earwax Type (ABCC11)",
                "result": "Dry earwax + minimal axillary odor",
                "evidence": "rs17822931 genotype: TT",
                "confidence": "high",
                "detail": (
                    "Almost universal in East Asians and Native Americans (>95%). "
                    "Apocrine glands produce minimal odor — many T/T individuals "
                    "don't need deodorant. Strong ancestry marker."
                )}
    if "T" in g:
        return {"trait": "Earwax Type (ABCC11)",
                "result": "Wet earwax (heterozygous CT — wet phenotype is dominant)",
                "evidence": f"rs17822931 genotype: {g}",
                "confidence": "high"}
    return {"trait": "Earwax Type (ABCC11)",
            "result": "Wet earwax (CC)",
            "evidence": "rs17822931 genotype: CC",
            "confidence": "high"}


def _trait_eye_color(snps_df) -> dict:
    # Primary determinant: HERC2 rs12913832 G allele = blue eyes
    herc2 = _gt(snps_df, "rs12913832")
    tyr = _dose(snps_df, "rs1042602", "A")
    if herc2 is None:
        return {"trait": "Eye Color", "result": "Not tested",
                "evidence": "rs12913832 (HERC2) not called", "confidence": "n/a"}
    if herc2 == "GG":
        # Add gradient via TYR/SLC24A4
        modifier = ""
        if tyr is not None and tyr > 0:
            modifier = " (with possible lighter/gray gradient via TYR modifier)"
        return {"trait": "Eye Color",
                "result": f"Blue eyes likely{modifier}",
                "evidence": "rs12913832: GG (HERC2 blue-eye haplotype)",
                "confidence": "high"}
    if herc2 == "AG" or herc2 == "GA":
        return {"trait": "Eye Color",
                "result": "Green, hazel, or light brown — intermediate phenotype",
                "evidence": "rs12913832: AG",
                "confidence": "moderate",
                "detail": "Heterozygous HERC2 typically yields green or hazel eyes."}
    return {"trait": "Eye Color",
            "result": "Brown eyes likely",
            "evidence": f"rs12913832: {herc2}",
            "confidence": "high"}


def _trait_hair_color(snps_df) -> dict:
    # MC1R RHC variants increase red/blonde tendency
    r151c = _dose(snps_df, "rs1805007", "T")
    r160w = _dose(snps_df, "rs1805008", "T")
    d294h = _dose(snps_df, "rs1805009", "C")
    rhc_total = sum(x for x in [r151c, r160w, d294h] if x is not None)
    irf4 = _dose(snps_df, "rs12203592", "T")  # lighter pigmentation

    if rhc_total is None and irf4 is None:
        return {"trait": "Hair Color Tendency", "result": "Not tested", "evidence": "",
                "confidence": "n/a"}

    if rhc_total >= 2:
        return {"trait": "Hair Color Tendency",
                "result": "Red hair likely (≥2 MC1R RHC alleles)",
                "evidence": f"MC1R RHC alleles: {rhc_total}",
                "confidence": "high",
                "detail": (
                    "Red hair carriers have higher melanoma risk per UV dose. "
                    "Use rigorous UV protection: daily broad-spectrum SPF 30+, "
                    "hats, sunglasses. Annual skin checks if multiple atypical moles."
                )}
    if rhc_total == 1:
        return {"trait": "Hair Color Tendency",
                "result": "Red/blonde tinge possible (single MC1R RHC carrier)",
                "evidence": "MC1R RHC alleles: 1",
                "confidence": "moderate"}
    if irf4 and irf4 >= 1:
        return {"trait": "Hair Color Tendency",
                "result": "Light/blonde hair tendency",
                "evidence": "IRF4 rs12203592 T allele present",
                "confidence": "moderate"}
    return {"trait": "Hair Color Tendency",
            "result": "Dark hair likely",
            "evidence": "No MC1R red-hair alleles or IRF4 lightening allele detected",
            "confidence": "moderate"}


def _trait_chronotype(snps_df) -> dict:
    clock = _gt(snps_df, "rs1801260")  # C = evening preference
    per3 = _gt(snps_df, "rs2230912")
    if clock is None and per3 is None:
        return {"trait": "Chronotype", "result": "Not tested", "evidence": "", "confidence": "n/a"}
    if clock and "C" in clock:
        return {"trait": "Chronotype",
                "result": "Tendency toward evening chronotype (night owl)",
                "evidence": f"rs1801260 (CLOCK) genotype: {clock}",
                "confidence": "moderate",
                "detail": (
                    "Functions better with later sleep/wake times. Forced early "
                    "schedules cause social jet lag, mood and metabolic issues. "
                    "Advocate for later work start when possible. Bright morning "
                    "light + evening dimming helps any chronotype."
                )}
    return {"trait": "Chronotype",
            "result": "Likely morning or neutral chronotype",
            "evidence": f"rs1801260 (CLOCK) genotype: {clock or 'not called'}",
            "confidence": "low"}


def _trait_short_sleeper(snps_df) -> dict:
    g = _gt(snps_df, "rs77086077")  # BHLHE41 / DEC2 short-sleeper variant
    if g is None:
        return {"trait": "Short-Sleeper Allele (BHLHE41/DEC2)", "result": "Not tested",
                "evidence": "rs77086077 not called", "confidence": "n/a"}
    if "T" in g:  # rare allele
        return {"trait": "Short-Sleeper Allele (BHLHE41/DEC2)",
                "result": "Carrier of short-sleeper variant",
                "evidence": f"rs77086077 genotype: {g}",
                "confidence": "high",
                "detail": (
                    "May function on 5–6 h sleep without obvious deficits. Rare. "
                    "Most people reporting 'short sleep' are actually chronically "
                    "sleep-deprived — test by sleeping 7–9 h consistently for "
                    "4 weeks and assess."
                )}
    return {"trait": "Short-Sleeper Allele (BHLHE41/DEC2)",
            "result": "Non-carrier — typical 7–9 h sleep need",
            "evidence": f"rs77086077 genotype: {g}",
            "confidence": "high"}


def _trait_muscle_fiber(snps_df) -> dict:
    g = _gt(snps_df, "rs1815739")  # ACTN3 R577X — T allele = stop codon
    if g is None:
        return {"trait": "Muscle Fiber Composition (ACTN3)", "result": "Not tested",
                "evidence": "rs1815739 not called", "confidence": "n/a"}
    if g == "TT":
        return {"trait": "Muscle Fiber Composition (ACTN3)",
                "result": "α-actinin-3 deficient (XX) — slight endurance bias",
                "evidence": "rs1815739 genotype: TT",
                "confidence": "high",
                "detail": (
                    "No α-actinin-3 in fast-twitch fibers (~18% of people). Sprint/"
                    "power performance is slightly attenuated at elite level; "
                    "recreational training response is fully normal."
                )}
    if g == "CC":
        return {"trait": "Muscle Fiber Composition (ACTN3)",
                "result": "Full α-actinin-3 (RR) — slight power/sprint bias",
                "evidence": "rs1815739 genotype: CC",
                "confidence": "high",
                "detail": "Elite sprint and power athletes are enriched for RR."}
    return {"trait": "Muscle Fiber Composition (ACTN3)",
            "result": "Heterozygous (RX) — typical composition",
            "evidence": f"rs1815739 genotype: {g}",
            "confidence": "high"}


def _trait_vitd_efficiency(snps_df) -> dict:
    dbp = _dose(snps_df, "rs2282679", "C")
    cyp2r1 = _dose(snps_df, "rs10741657", "G")
    if dbp is None and cyp2r1 is None:
        return {"trait": "Vitamin D Synthesis Efficiency", "result": "Not tested",
                "evidence": "", "confidence": "n/a"}
    score = (dbp or 0) + (cyp2r1 or 0)
    if score >= 3:
        return {"trait": "Vitamin D Synthesis Efficiency",
                "result": "Lower-than-average vitamin D synthesis / serum 25-OH-D",
                "evidence": f"DBP + CYP2R1 risk-allele dosage: {score}",
                "confidence": "moderate",
                "detail": (
                    "May require higher dietary or supplemental vitamin D to reach "
                    "the same 25-OH-D level. Test 25-OH-D; aim 30–50 ng/mL."
                )}
    return {"trait": "Vitamin D Synthesis Efficiency",
            "result": "Typical vitamin D synthesis",
            "evidence": f"DBP + CYP2R1 risk-allele dosage: {score}",
            "confidence": "moderate"}


def _trait_caffeine_anxiety(snps_df) -> dict:
    g = _gt(snps_df, "rs5751876")
    if g is None:
        return {"trait": "Caffeine-Induced Anxiety Susceptibility", "result": "Not tested",
                "evidence": "rs5751876 not called", "confidence": "n/a"}
    if g == "TT":
        return {"trait": "Caffeine-Induced Anxiety Susceptibility",
                "result": "Increased susceptibility to caffeine-induced anxiety (ADORA2A T/T)",
                "evidence": "rs5751876 genotype: TT",
                "confidence": "moderate",
                "detail": "Strong anxiety response to caffeine. Limit to <100 mg/day or switch to decaf. L-theanine 200 mg can blunt caffeine-related anxiety."}
    return {"trait": "Caffeine-Induced Anxiety Susceptibility",
            "result": "Typical caffeine response",
            "evidence": f"rs5751876 genotype: {g}",
            "confidence": "moderate"}


def _trait_smoking_dependence(snps_df) -> dict:
    chrna3 = _dose(snps_df, "rs1051730", "T")
    chrna5 = _dose(snps_df, "rs16969968", "A")
    score = (chrna3 or 0) + (chrna5 or 0)
    if chrna3 is None and chrna5 is None:
        return {"trait": "Nicotine Dependence Susceptibility", "result": "Not tested",
                "evidence": "", "confidence": "n/a"}
    if score >= 2:
        return {"trait": "Nicotine Dependence Susceptibility",
                "result": "Elevated nicotine-dependence and lung-cancer-in-smokers risk",
                "evidence": f"CHRNA3/5 risk-allele dosage: {score}",
                "confidence": "high",
                "detail": (
                    "If smoking: vigorous multimodal cessation (varenicline + NRT + "
                    "behavioral support). If non-smoker: maintain — initiation risk "
                    "is elevated."
                )}
    return {"trait": "Nicotine Dependence Susceptibility",
            "result": "Typical risk",
            "evidence": f"CHRNA3/5 risk-allele dosage: {score}",
            "confidence": "moderate"}


# ─── V4 trait analyzers (additional ~30) ──────────────────────────────────────

def _trait_male_pattern_baldness(snps_df) -> dict:
    g = _gt(snps_df, "rs6152")
    if g is None:
        return {"trait": "Male-Pattern Baldness Tendency", "result": "Not tested",
                "evidence": "rs6152 not called", "confidence": "n/a"}
    # AR locus — risk allele varies but A (vs G) usually associates with higher AR activity
    if "A" in g:
        return {"trait": "Male-Pattern Baldness Tendency",
                "result": "Elevated androgenetic alopecia susceptibility",
                "evidence": f"AR rs6152 genotype: {g}", "confidence": "moderate",
                "detail": "Topical minoxidil, oral/topical finasteride options. Consult dermatology if concerned."}
    return {"trait": "Male-Pattern Baldness Tendency",
            "result": "Lower androgenetic alopecia tendency",
            "evidence": f"AR rs6152 genotype: {g}", "confidence": "moderate"}


def _trait_photic_sneeze(snps_df) -> dict:
    g = _gt(snps_df, "rs10427255")
    if g is None:
        return {"trait": "Photic Sneeze Reflex (ACHOO)", "result": "Not tested",
                "evidence": "rs10427255 not called", "confidence": "n/a"}
    if "T" in g:
        return {"trait": "Photic Sneeze Reflex (ACHOO)",
                "result": "Likely photic sneezer",
                "evidence": f"rs10427255: {g}", "confidence": "moderate",
                "detail": "Sneezing when stepping into bright sunlight is harmless; sunglasses mitigate."}
    return {"trait": "Photic Sneeze Reflex (ACHOO)",
            "result": "Unlikely to be a photic sneezer",
            "evidence": f"rs10427255: {g}", "confidence": "moderate"}


def _trait_cilantro_aversion(snps_df) -> dict:
    g = _gt(snps_df, "rs72921001")
    if g is None:
        return {"trait": "Cilantro Aversion (OR6A2)", "result": "Not tested",
                "evidence": "rs72921001 not called", "confidence": "n/a"}
    if "C" in g:
        return {"trait": "Cilantro Aversion (OR6A2)",
                "result": "Likely to perceive cilantro as soapy",
                "evidence": f"OR6A2 rs72921001: {g}", "confidence": "high",
                "detail": "Substitute parsley, basil, or mint where cilantro is called for."}
    return {"trait": "Cilantro Aversion (OR6A2)",
            "result": "Likely to enjoy cilantro",
            "evidence": f"OR6A2 rs72921001: {g}", "confidence": "high"}


def _trait_asparagus_smell(snps_df) -> dict:
    g = _gt(snps_df, "rs4481887")
    if g is None:
        return {"trait": "Asparagus Urine Smell Detection",
                "result": "Not tested", "evidence": "rs4481887 not called", "confidence": "n/a"}
    if "G" in g:
        return {"trait": "Asparagus Urine Smell Detection",
                "result": "Likely to detect the characteristic asparagus urine smell",
                "evidence": f"rs4481887: {g}", "confidence": "moderate"}
    return {"trait": "Asparagus Urine Smell Detection",
            "result": "Likely anosmic to asparagus urine metabolites",
            "evidence": f"rs4481887: {g}", "confidence": "moderate",
            "detail": "Phenotypic curiosity — only some people smell asparagus metabolites in urine."}


def _trait_body_odor(snps_df) -> dict:
    g = _gt(snps_df, "rs17822931")
    if g is None:
        return {"trait": "Body Odor Type (ABCC11)", "result": "Not tested",
                "evidence": "rs17822931 not called", "confidence": "n/a"}
    if g == "TT":
        return {"trait": "Body Odor Type (ABCC11)",
                "result": "Minimal axillary odor (dry earwax phenotype)",
                "evidence": f"ABCC11: {g}", "confidence": "high",
                "detail": "Common in East Asian and Native American populations."}
    return {"trait": "Body Odor Type (ABCC11)",
            "result": "Wet earwax / typical apocrine secretion",
            "evidence": f"ABCC11: {g}", "confidence": "high"}


def _trait_freckling(snps_df) -> dict:
    irf4 = _dose(snps_df, "rs12203592", "T")
    mc1r_total = sum(d for d in [
        _dose(snps_df, "rs1805007", "T"),
        _dose(snps_df, "rs1805008", "T"),
        _dose(snps_df, "rs1805009", "C"),
    ] if d is not None)
    if irf4 is None and mc1r_total == 0:
        return {"trait": "Freckling Tendency", "result": "Not tested",
                "evidence": "Key freckling SNPs not called", "confidence": "n/a"}
    if (irf4 and irf4 >= 1) or mc1r_total >= 1:
        return {"trait": "Freckling Tendency",
                "result": "Higher freckling tendency",
                "evidence": f"IRF4 dose={irf4}, MC1R RHC alleles={mc1r_total}",
                "confidence": "moderate",
                "detail": "UV protection is particularly important; freckling correlates with melanoma risk in fair-skinned individuals."}
    return {"trait": "Freckling Tendency", "result": "Lower freckling tendency",
            "evidence": f"IRF4 dose={irf4}, MC1R RHC alleles={mc1r_total}",
            "confidence": "moderate"}


def _trait_empathy(snps_df) -> dict:
    g = _gt(snps_df, "rs53576")
    if g is None:
        return {"trait": "Empathy / Social Sensitivity (OXTR)", "result": "Not tested",
                "evidence": "rs53576 not called", "confidence": "n/a"}
    if g == "GG":
        return {"trait": "Empathy / Social Sensitivity (OXTR)",
                "result": "Higher empathic accuracy tendency",
                "evidence": "OXTR: GG", "confidence": "low",
                "detail": "Modest effect; lived experience dominates."}
    if "G" in g:
        return {"trait": "Empathy / Social Sensitivity (OXTR)",
                "result": "Intermediate empathy tendency",
                "evidence": f"OXTR: {g}", "confidence": "low"}
    return {"trait": "Empathy / Social Sensitivity (OXTR)",
            "result": "Lower empathy tendency (modest effect)",
            "evidence": f"OXTR: {g}", "confidence": "low"}


def _trait_salt_sensitivity(snps_df) -> dict:
    g = _gt(snps_df, "rs5443")
    if g is None:
        return {"trait": "Salt-Sensitivity of Blood Pressure", "result": "Not tested",
                "evidence": "GNB3 rs5443 not called", "confidence": "n/a"}
    if "T" in g:
        return {"trait": "Salt-Sensitivity of Blood Pressure",
                "result": "Higher salt-sensitivity of BP",
                "evidence": f"GNB3 rs5443: {g}", "confidence": "moderate",
                "detail": "DASH diet; sodium <2300 mg/day (1500 mg/day if hypertensive); potassium-rich foods."}
    return {"trait": "Salt-Sensitivity of Blood Pressure",
            "result": "Lower salt-sensitivity",
            "evidence": f"GNB3 rs5443: {g}", "confidence": "moderate"}


def _trait_endurance_vs_power(snps_df) -> dict:
    actn3 = _gt(snps_df, "rs1815739")
    if actn3 is None:
        return {"trait": "Endurance vs Power Bias", "result": "Not tested",
                "evidence": "ACTN3 not called", "confidence": "n/a"}
    if actn3 == "TT":
        return {"trait": "Endurance vs Power Bias",
                "result": "Endurance-typical (ACTN3 null)",
                "evidence": "ACTN3: TT", "confidence": "high",
                "detail": "Recreational training response is normal; effects matter most at elite level."}
    if actn3 == "CC":
        return {"trait": "Endurance vs Power Bias",
                "result": "Power/sprint-typical (ACTN3 R/R)",
                "evidence": "ACTN3: CC", "confidence": "high"}
    return {"trait": "Endurance vs Power Bias",
            "result": "Mixed type (heterozygous ACTN3)",
            "evidence": f"ACTN3: {actn3}", "confidence": "high"}


def _trait_injury_susceptibility(snps_df) -> dict:
    col1 = _dose(snps_df, "rs1800012", "T")
    col5 = _dose(snps_df, "rs12722", "T")
    total = (col1 or 0) + (col5 or 0)
    if col1 is None and col5 is None:
        return {"trait": "Tendon/Ligament Injury Susceptibility",
                "result": "Not tested", "evidence": "", "confidence": "n/a"}
    if total >= 2:
        return {"trait": "Tendon/Ligament Injury Susceptibility",
                "result": "Elevated soft-tissue injury risk",
                "evidence": f"COL1A1 + COL5A1 risk-allele dose: {total}",
                "confidence": "moderate",
                "detail": "Thorough warmup, progressive loading, eccentric strength work, adequate protein/vitamin C for collagen synthesis."}
    return {"trait": "Tendon/Ligament Injury Susceptibility",
            "result": "Typical injury risk",
            "evidence": f"COL1A1 + COL5A1 dose: {total}",
            "confidence": "moderate"}


def _trait_bone_density(snps_df) -> dict:
    lrp5 = _dose(snps_df, "rs3736228", "T")
    vdr = _dose(snps_df, "rs2228570", "T")
    score = (lrp5 or 0) + (vdr or 0)
    if lrp5 is None and vdr is None:
        return {"trait": "Bone Density Tendency", "result": "Not tested",
                "evidence": "", "confidence": "n/a"}
    if score >= 2:
        return {"trait": "Bone Density Tendency",
                "result": "Lower bone density tendency",
                "evidence": f"LRP5+VDR risk dose: {score}",
                "confidence": "moderate",
                "detail": "Weight-bearing/resistance exercise; calcium 1000-1200 mg/day; vitamin D ≥30 ng/mL; vitamin K2; DEXA from menopause or age 65."}
    return {"trait": "Bone Density Tendency", "result": "Typical bone density tendency",
            "evidence": f"LRP5+VDR risk dose: {score}", "confidence": "moderate"}


def _trait_pain_sensitivity(snps_df) -> dict:
    comt = _gt(snps_df, "rs4680")
    if comt is None:
        return {"trait": "Pain Sensitivity (COMT)", "result": "Not tested",
                "evidence": "rs4680 not called", "confidence": "n/a"}
    if comt == "AA":
        # Met/Met — slow COMT, generally higher pain sensitivity
        return {"trait": "Pain Sensitivity (COMT)",
                "result": "Higher pain sensitivity (Met/Met)",
                "evidence": "COMT: AA", "confidence": "moderate",
                "detail": "May benefit from multimodal pain management; address anxiety component."}
    if comt == "GG":
        return {"trait": "Pain Sensitivity (COMT)",
                "result": "Lower pain sensitivity (Val/Val)",
                "evidence": "COMT: GG", "confidence": "moderate"}
    return {"trait": "Pain Sensitivity (COMT)",
            "result": "Intermediate pain sensitivity",
            "evidence": f"COMT: {comt}", "confidence": "moderate"}


def _trait_novelty_seeking(snps_df) -> dict:
    drd4 = _gt(snps_df, "rs1800955")
    if drd4 is None:
        return {"trait": "Novelty Seeking Tendency", "result": "Not tested",
                "evidence": "rs1800955 not called", "confidence": "n/a"}
    if "T" in drd4:
        return {"trait": "Novelty Seeking Tendency",
                "result": "Higher novelty-seeking trait",
                "evidence": f"DRD4 rs1800955: {drd4}", "confidence": "low",
                "detail": "Modest effect; behavioral patterns and context dominate."}
    return {"trait": "Novelty Seeking Tendency",
            "result": "Lower novelty-seeking tendency",
            "evidence": f"DRD4 rs1800955: {drd4}", "confidence": "low"}


def _trait_height_polygenic(snps_df) -> dict:
    # Quick polygenic-ish estimate from a handful of well-known height SNPs
    height_snps = [
        ("rs6060369", "T"), ("rs2562784", "G"), ("rs143384", "T"),
        ("rs2275035", "C"),
    ]
    doses = []
    for rsid, ea in height_snps:
        d = _dose(snps_df, rsid, ea)
        if d is not None:
            doses.append(d)
    if not doses:
        return {"trait": "Height Polygenic Tendency", "result": "Not tested",
                "evidence": "Height SNPs not called", "confidence": "n/a"}
    avg = sum(doses) / len(doses)
    if avg > 1.2:
        result = "Above-average height tendency"
    elif avg < 0.8:
        result = "Below-average height tendency"
    else:
        result = "Average height tendency"
    return {"trait": "Height Polygenic Tendency",
            "result": result,
            "evidence": f"Avg risk-allele dose across {len(doses)} height SNPs: {avg:.2f}",
            "confidence": "low",
            "detail": "Polygenic — many SNPs contribute. Strong genetic component but childhood nutrition matters too."}


def _trait_memory_bdnf(snps_df) -> dict:
    g = _gt(snps_df, "rs6265")
    if g is None:
        return {"trait": "Memory / Neuroplasticity (BDNF)", "result": "Not tested",
                "evidence": "rs6265 not called", "confidence": "n/a"}
    if "A" in g:
        return {"trait": "Memory / Neuroplasticity (BDNF)",
                "result": "BDNF Val66Met — reduced activity-dependent BDNF secretion",
                "evidence": f"BDNF: {g}", "confidence": "moderate",
                "detail": "Aerobic exercise robustly elevates BDNF and counteracts the effect. Sleep, mental stimulation, omega-3."}
    return {"trait": "Memory / Neuroplasticity (BDNF)",
            "result": "Val/Val — typical BDNF secretion",
            "evidence": f"BDNF: {g}", "confidence": "moderate"}


def _trait_warfarin_sensitivity(snps_df) -> dict:
    vkorc = _gt(snps_df, "rs9923231")
    if vkorc is None:
        return {"trait": "Warfarin Sensitivity", "result": "Not tested",
                "evidence": "VKORC1 not called", "confidence": "n/a"}
    if vkorc == "AA":
        return {"trait": "Warfarin Sensitivity",
                "result": "High warfarin sensitivity — low starting dose if ever prescribed",
                "evidence": "VKORC1 -1639G>A: AA", "confidence": "high"}
    if "A" in vkorc:
        return {"trait": "Warfarin Sensitivity",
                "result": "Intermediate warfarin sensitivity",
                "evidence": f"VKORC1: {vkorc}", "confidence": "high"}
    return {"trait": "Warfarin Sensitivity",
            "result": "Typical warfarin sensitivity",
            "evidence": f"VKORC1: {vkorc}", "confidence": "high"}


def _trait_cortisol_response(snps_df) -> dict:
    g = _gt(snps_df, "rs41423247")
    if g is None:
        return {"trait": "Cortisol / Stress Response (NR3C1)", "result": "Not tested",
                "evidence": "rs41423247 not called", "confidence": "n/a"}
    if g == "GG":
        return {"trait": "Cortisol / Stress Response (NR3C1)",
                "result": "More sensitive cortisol response",
                "evidence": "NR3C1 BclI: GG", "confidence": "moderate",
                "detail": "Stress-management practices matter more; avoid chronic overtraining."}
    return {"trait": "Cortisol / Stress Response (NR3C1)",
            "result": "Typical cortisol response",
            "evidence": f"NR3C1 BclI: {g}", "confidence": "moderate"}


def _trait_thyroid_t3_response(snps_df) -> dict:
    g = _gt(snps_df, "rs225014")
    if g is None:
        return {"trait": "Thyroid T4→T3 Conversion (DIO2)", "result": "Not tested",
                "evidence": "rs225014 not called", "confidence": "n/a"}
    if g == "CC":
        return {"trait": "Thyroid T4→T3 Conversion (DIO2)",
                "result": "Potentially suboptimal T4→T3 conversion",
                "evidence": "DIO2 T92A: CC", "confidence": "low",
                "detail": "If hypothyroid and feeling suboptimal on levothyroxine, discuss combination T4/T3 with endocrinologist."}
    return {"trait": "Thyroid T4→T3 Conversion (DIO2)",
            "result": "Typical T4→T3 conversion",
            "evidence": f"DIO2 T92A: {g}", "confidence": "low"}


def _trait_endurance_vo2(snps_df) -> dict:
    g = _gt(snps_df, "rs8192678")
    if g is None:
        return {"trait": "Aerobic Trainability (PPARGC1A)", "result": "Not tested",
                "evidence": "rs8192678 not called", "confidence": "n/a"}
    if "A" in g:
        return {"trait": "Aerobic Trainability (PPARGC1A)",
                "result": "Less-responsive PGC-1α — emphasises training consistency",
                "evidence": f"PPARGC1A Gly482Ser: {g}", "confidence": "moderate",
                "detail": "Need more total training volume / consistency to maximize adaptation. HIIT effective."}
    return {"trait": "Aerobic Trainability (PPARGC1A)",
            "result": "Responsive PGC-1α — typical aerobic trainability",
            "evidence": f"PPARGC1A: {g}", "confidence": "moderate"}


def _trait_lifespan_longevity_snps(snps_df) -> dict:
    foxo3 = _dose(snps_df, "rs2802292", "G")
    apoe_e2 = _dose(snps_df, "rs7412", "T")
    if foxo3 is None and apoe_e2 is None:
        return {"trait": "Longevity-Associated Variants", "result": "Not tested",
                "evidence": "", "confidence": "n/a"}
    score = (foxo3 or 0) + (apoe_e2 or 0)
    if score >= 2:
        result = "Multiple longevity-associated alleles present"
    elif score == 1:
        result = "Some longevity-associated alleles present"
    else:
        result = "Standard genetic longevity profile"
    return {"trait": "Longevity-Associated Variants",
            "result": result,
            "evidence": f"FOXO3 G dose: {foxo3}; APOE ε2 dose: {apoe_e2}",
            "confidence": "low",
            "detail": "Lifestyle dominates; Mediterranean diet, exercise, sleep, social engagement are the strongest interventions."}


# ─── Master analyzer ──────────────────────────────────────────────────────────

def predict_traits(snps_df: pd.DataFrame) -> dict:
    analyzers = [
        _trait_lactose, _trait_alcohol_flush, _trait_caffeine_speed,
        _trait_bitter_taste, _trait_earwax, _trait_eye_color,
        _trait_hair_color, _trait_chronotype, _trait_short_sleeper,
        _trait_muscle_fiber, _trait_vitd_efficiency, _trait_caffeine_anxiety,
        _trait_smoking_dependence,
        # V4 additions
        _trait_male_pattern_baldness, _trait_photic_sneeze,
        _trait_cilantro_aversion, _trait_asparagus_smell, _trait_body_odor,
        _trait_freckling, _trait_empathy, _trait_salt_sensitivity,
        _trait_endurance_vs_power, _trait_injury_susceptibility,
        _trait_bone_density, _trait_pain_sensitivity, _trait_novelty_seeking,
        _trait_height_polygenic, _trait_memory_bdnf, _trait_warfarin_sensitivity,
        _trait_cortisol_response, _trait_thyroid_t3_response,
        _trait_endurance_vo2, _trait_lifespan_longevity_snps,
    ]
    results = []
    for a in analyzers:
        try:
            results.append(a(snps_df))
        except Exception as e:
            results.append({"trait": a.__name__, "result": f"Analysis failed: {e}",
                            "evidence": "", "confidence": "error"})
    return {
        "predictions": results,
        "n_predictions": sum(1 for r in results if r.get("confidence") not in ("n/a", "error")),
        "n_not_tested": sum(1 for r in results if r.get("confidence") == "n/a"),
    }


# ── V8: cross-check against the unified SNP registry ──────────────────────
# Same shape as wellness.audit_against_registry — traits.py uses inline
# rsID literals in rule conditionals, so the audit is presence-only.

def _scan_rsids_referenced() -> list[str]:
    src = _Path(__file__).read_text()
    return sorted(set(_re.findall(r'"(rs\d+)"', src)))


def audit_against_registry() -> dict[str, list[str]]:
    """Return ``{"registered": [...], "missing": [...]}`` for every rsID
    referenced anywhere in this module.

    traits.py covers many phenotype SNPs (eye/hair color, taste, earwax,
    skin pigmentation, MC1R variants) that are well-characterized in the
    literature but had not been added to the unified registry in V8.
    The 32 deferred rsIDs are documented in ``CHANGELOG.md`` under
    "V8.1 follow-ups" — each needs literature-cited ancestral/derived +
    GRCh37/38 coordinates before joining the registry.
    """
    referenced = _scan_rsids_referenced()
    registered: list[str] = []
    missing: list[str] = []
    for r in referenced:
        if snp_registry.get(r) is not None:
            registered.append(r)
        else:
            missing.append(r)
    return {"registered": registered, "missing": missing}


_AUDIT = audit_against_registry()
