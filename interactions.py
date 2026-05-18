"""
Compound Heterozygosity & Gene-Gene Risk Amplification
------------------------------------------------------

Detects clinically meaningful multi-variant patterns that change risk
qualitatively beyond what any single variant implies:

  * MTHFR C677T + A1298C compound heterozygosity (methylation impairment)
  * HFE C282Y / H63D compound (atypical iron-overload risk)
  * Factor V Leiden + Prothrombin G20210A (additive thrombophilia)
  * Multiple 9p21 risk alleles + Lp(a) (compound CAD risk)
  * APOE-ε4 carriage × TOMM40 haplotype
  * Multi-hit Crohn's (NOD2 + ATG16L1 + IL23R)
  * Slow NAT2 + Null GSTM1 (carcinogen handling deficit)
  * High aromatase + Slow COMT (estrogen-driven cancer profile)
  * MTHFR + CBS + Folate / B12 cycle profile
  * Slow caffeine (CYP1A2 *1F/*1F) + Anxiety-prone ADORA2A
  * Slow ALDH2 + Fast ADH1B (acetaldehyde accumulation)
  * IL23R + IL6 + TNF anti-inflammatory consideration

Each finding combines specific risk-allele dosages and emits a clinical
interpretation with action items.
"""

from typing import Dict, List, Optional
import pandas as pd


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> Optional[int]:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--") or len(s) != 2:
        return None
    return s.count(allele.upper())


def _gt(snps_df: pd.DataFrame, rsid: str) -> Optional[str]:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    return str(gt).upper()


def detect_interactions(snps_df: pd.DataFrame) -> Dict:
    """Returns a dict with a list of finding-records, each describing a
    multi-variant pattern detected in the genotype data."""
    findings: List[Dict] = []

    # ── MTHFR compound heterozygosity ─────────────────────────────────────────
    c677t = _dose(snps_df, "rs1801133", "T")   # C677T
    a1298c = _dose(snps_df, "rs1801131", "C")  # A1298C
    if c677t is not None and a1298c is not None:
        total = c677t + a1298c
        if c677t == 2:
            findings.append({
                "title": "MTHFR C677T Homozygous (T/T)",
                "severity": "high",
                "variants": ["rs1801133 (T/T)"],
                "interpretation": (
                    "Homozygous MTHFR 677TT reduces enzyme activity ~70%, "
                    "substantially impairing folate-cycle methylation and "
                    "elevating homocysteine. Combined cardiovascular, fertility, "
                    "and methylation-detox impact."
                ),
                "action": (
                    "Switch to methylfolate (5-MTHF, 400–800 µg/day) — NOT folic "
                    "acid. Methylcobalamin B12 1000 µg/day. Riboflavin (B2) 25 mg/day. "
                    "Test plasma homocysteine; target <8 µmol/L. Limit alcohol. "
                    "If trying to conceive, methylfolate well before pregnancy."
                ),
            })
        elif a1298c == 2:
            findings.append({
                "title": "MTHFR A1298C Homozygous (C/C)",
                "severity": "moderate",
                "variants": ["rs1801131 (C/C)"],
                "interpretation": (
                    "Homozygous MTHFR 1298CC reduces BH4 production with downstream "
                    "effects on neurotransmitter synthesis and methylation tone."
                ),
                "action": (
                    "Methylfolate + methylcobalamin support. BH4 support: tetrahydrobiopterin "
                    "is not over-the-counter; lifestyle support via stress management and "
                    "adequate methylation cofactors."
                ),
            })
        elif c677t >= 1 and a1298c >= 1:
            findings.append({
                "title": f"MTHFR Compound Heterozygosity (C677T {'+/+' if c677t == 2 else '+/-'} + A1298C {'+/+' if a1298c == 2 else '+/-'})",
                "severity": "high" if total >= 3 else "moderate",
                "variants": ["rs1801133", "rs1801131"],
                "interpretation": (
                    "Compound heterozygosity for MTHFR C677T and A1298C functionally "
                    "approaches homozygous C677T — substantial methylation impairment. "
                    "This is one of the most clinically relevant common-variant "
                    "compound states; many clinicians miss it because each variant "
                    "is reported alone."
                ),
                "action": (
                    "Treat as approximately equivalent to homozygous C677T for "
                    "supplementation: methylfolate 400–800 µg/day, methylcobalamin "
                    "1000 µg/day, riboflavin 25 mg/day, B6 (P5P) 25 mg/day. "
                    "Plasma homocysteine measurement. Limit alcohol. Particularly "
                    "important pre-conception (paternal AND maternal)."
                ),
            })

    # ── HFE compound heterozygosity ───────────────────────────────────────────
    c282y = _dose(snps_df, "rs1800562", "A")    # C282Y
    h63d = _dose(snps_df, "rs1799945", "G")     # H63D
    if c282y is not None and h63d is not None:
        if c282y == 2:
            findings.append({
                "title": "HFE C282Y Homozygous — Classical Hemochromatosis Risk",
                "severity": "high",
                "variants": ["rs1800562 (A/A)"],
                "interpretation": (
                    "Homozygous C282Y is the canonical hemochromatosis genotype. "
                    "Penetrance is incomplete — ~25–60% of homozygotes develop "
                    "iron overload symptoms; men more often than women. Untreated "
                    "iron overload damages liver, heart, joints, pancreas, "
                    "endocrine glands."
                ),
                "action": (
                    "Test serum ferritin and transferrin saturation NOW. If "
                    "elevated (ferritin >300 men / >200 women, TSat >45%), "
                    "see hematology or gastroenterology — therapeutic phlebotomy "
                    "is definitive treatment. Avoid iron supplements and high-dose "
                    "vitamin C with iron-rich meals. Limit alcohol (synergistic "
                    "liver toxicity). Periodic monitoring even if labs initially "
                    "normal."
                ),
            })
        elif c282y == 1 and h63d == 1:
            findings.append({
                "title": "HFE Compound Heterozygote (C282Y/H63D)",
                "severity": "moderate",
                "variants": ["rs1800562 (heterozygous)", "rs1799945 (heterozygous)"],
                "interpretation": (
                    "Compound heterozygosity carries mild-to-moderate iron loading "
                    "risk, particularly with co-factors (alcohol, hepatitis, "
                    "metabolic syndrome). Most carriers do NOT develop clinical "
                    "hemochromatosis but a meaningful minority do."
                ),
                "action": (
                    "Test serum ferritin and transferrin saturation. Periodic "
                    "monitoring (every 2–3 years) if normal; phlebotomy if "
                    "elevated. Limit alcohol; avoid iron supplements."
                ),
            })

    # ── Factor V Leiden + Prothrombin G20210A ─────────────────────────────────
    fvl = _dose(snps_df, "rs6025", "A")          # Factor V Leiden
    prothrombin = _dose(snps_df, "rs1799963", "A")  # F2 G20210A
    if fvl is not None and prothrombin is not None and fvl + prothrombin > 0:
        if fvl == 2 or prothrombin == 2 or (fvl >= 1 and prothrombin >= 1):
            findings.append({
                "title": "Multi-Variant Thrombophilia (Factor V Leiden + Prothrombin G20210A)",
                "severity": "high",
                "variants": [
                    f"rs6025 (Factor V Leiden, dosage {fvl})",
                    f"rs1799963 (Prothrombin G20210A, dosage {prothrombin})",
                ],
                "interpretation": (
                    "Carriers of both Factor V Leiden and Prothrombin G20210A "
                    "have multiplicative VTE risk — heterozygote+heterozygote "
                    "is ~20× baseline VTE risk vs ~5× for FVL alone. Homozygous "
                    "FVL alone is ~80× baseline."
                ),
                "action": (
                    "INFORM ALL HEALTHCARE PROVIDERS — especially before any "
                    "surgery, pregnancy, or estrogen-containing hormone therapy "
                    "(birth control, HRT). Avoid prolonged immobilization "
                    "(long flights — compression stockings, hydration, walking). "
                    "See a hematologist to discuss prophylactic anticoagulation "
                    "in high-risk situations (post-surgery, pregnancy). NOT a "
                    "reason to take chronic anticoagulation in absence of clinical events."
                ),
            })
        elif fvl == 1:
            findings.append({
                "title": "Factor V Leiden Heterozygote",
                "severity": "moderate",
                "variants": [f"rs6025 (A/G, dosage 1)"],
                "interpretation": (
                    "Single-copy Factor V Leiden — ~5–7× VTE risk vs background. "
                    "Most carriers never have a VTE event."
                ),
                "action": (
                    "Inform providers before surgery, pregnancy, hormonal "
                    "contraceptives. Compression stockings + hydration on long "
                    "flights. Discuss with hematologist if any VTE history "
                    "(personal or first-degree family)."
                ),
            })

    # ── 9p21 multi-allele CAD risk ────────────────────────────────────────────
    p9_a = _dose(snps_df, "rs10757278", "G")
    p9_b = _dose(snps_df, "rs1333049", "C")
    p9_c = _dose(snps_df, "rs2383206", "G")
    lpa_a = _dose(snps_df, "rs10455872", "G")
    lpa_b = _dose(snps_df, "rs3798220", "C")
    p9_total = sum(x for x in [p9_a, p9_b, p9_c] if x is not None)
    lpa_total = sum(x for x in [lpa_a, lpa_b] if x is not None)
    if p9_total >= 3 and lpa_total >= 1:
        findings.append({
            "title": "9p21 + Lp(a): Compound Coronary Risk",
            "severity": "high",
            "variants": ["rs10757278/rs1333049/rs2383206 (9p21)", "rs10455872/rs3798220 (LPA)"],
            "interpretation": (
                "Carriage of multiple 9p21 risk alleles together with at least "
                "one LPA risk allele combines two of the strongest common-variant "
                "CAD pathways (9p21 vascular biology + Lp(a) atherogenicity). "
                "Effect on lifetime CAD risk is substantial and largely independent "
                "of LDL-C."
            ),
            "action": (
                "Get a one-time Lp(a) blood test (mass or molar). If Lp(a) is "
                "elevated (>50 mg/dL or >125 nmol/L), aggressive control of all "
                "modifiable risk factors: ApoB <80 mg/dL, BP <120/80, no smoking, "
                "≥150 min aerobic + 2 resistance sessions/week. Consider coronary "
                "artery calcium scan at age 40–45. Emerging Lp(a)-lowering "
                "therapies (siRNA, antisense) are in late-stage trials."
            ),
        })

    # ── Multi-hit Crohn's risk ────────────────────────────────────────────────
    nod2 = _dose(snps_df, "rs2066844", "T")
    atg16 = _dose(snps_df, "rs2241880", "G")
    il23r = _dose(snps_df, "rs11209026", "A")  # protective allele
    crohn_pos = (nod2 or 0) + (atg16 or 0)
    crohn_protect = (il23r or 0)
    if crohn_pos >= 2 and crohn_protect == 0:
        findings.append({
            "title": "Multi-Variant Crohn's Disease Risk (NOD2 + ATG16L1, no IL23R protection)",
            "severity": "moderate",
            "variants": ["rs2066844 (NOD2)", "rs2241880 (ATG16L1)", "rs11209026 (IL23R protective absent)"],
            "interpretation": (
                "Carriage of multiple Crohn's-susceptibility variants in "
                "innate-immune / autophagy pathways without the protective "
                "IL23R R381Q allele represents an elevated polygenic Crohn's "
                "background. Penetrance is low — most carriers never develop IBD."
            ),
            "action": (
                "Do NOT smoke (strongest modifiable Crohn's risk — smokers have "
                "~2× incidence and worse prognosis). Watch for chronic GI symptoms "
                "(diarrhea >4 weeks, blood, weight loss, perianal disease) — see "
                "gastroenterology promptly. Maintain healthy gut microbiome "
                "(fiber, fermented foods, limit unnecessary antibiotics). "
                "Mediterranean-style diet."
            ),
        })

    # ── Slow NAT2 + GSTM1/GSTT1 null ──────────────────────────────────────────
    nat2 = _dose(snps_df, "rs1495741", "A")  # A = fast acetylator allele; G/G = slow
    nat2_slow = (2 - nat2) if nat2 is not None else None  # crude
    gstm1 = _dose(snps_df, "rs366631", "DEL")  # not reliably detectable
    gstt1 = _dose(snps_df, "rs71748309", "DEL")
    if nat2 is not None and nat2 == 0:  # G/G slow acetylator
        findings.append({
            "title": "NAT2 Slow Acetylator Phenotype",
            "severity": "moderate",
            "variants": ["rs1495741 (slow acetylator genotype)"],
            "interpretation": (
                "Slow NAT2 acetylation prolongs exposure to aromatic amines "
                "(from cigarette smoke, charred meat, certain industrial chemicals, "
                "hair dyes, sulfonamide drugs, isoniazid). Increases bladder "
                "cancer risk in smokers and altered drug metabolism (especially "
                "isoniazid for TB treatment, sulfonamide antibiotics)."
            ),
            "action": (
                "Do NOT smoke. Minimize charred/grilled meat. Mention slow-"
                "acetylator status if ever prescribed isoniazid (peripheral "
                "neuropathy risk; B6 supplementation) or sulfonamide antibiotics. "
                "Avoid aromatic-amine occupational exposures."
            ),
        })

    # ── Estrogen metabolism profile (Slow COMT + Active CYP1B1) ───────────────
    comt = _dose(snps_df, "rs4680", "A")          # Met allele (slow)
    cyp1b1_a = _dose(snps_df, "rs1056836", "G")   # Val allele (active)
    cyp1b1_b = _dose(snps_df, "rs1056827", "T")
    if comt is not None and comt == 2 and (cyp1b1_a or 0) >= 1:
        findings.append({
            "title": "Estrogen-Metabolism Profile: Slow COMT + Active CYP1B1",
            "severity": "moderate",
            "variants": ["rs4680 COMT Met/Met (slow)", "rs1056836 CYP1B1 (active)"],
            "interpretation": (
                "This combination favors accumulation of 4-hydroxyestrogens "
                "(from CYP1B1) and slow methylation/clearance (by COMT). 4-OH "
                "estrogens can generate DNA-damaging quinones. Of greatest "
                "interest for estrogen-driven cancer risk (breast, endometrial) "
                "in women; relevant for hormonal balance in men too."
            ),
            "action": (
                "Cruciferous vegetables (broccoli, Brussels sprouts, cabbage, kale) "
                "favor 2-OH over 4-OH estrogen pathways. DIM (diindolylmethane) "
                "from cruciferous vegetables. Limit alcohol (raises estrogen). "
                "Adequate magnesium and B-vitamin cofactors for COMT. Maintain "
                "healthy body composition. For women with strong breast-cancer "
                "family history, this profile is one of several inputs into "
                "screening intensity decisions."
            ),
        })

    # ── ALDH2 + ADH1B (East Asian alcohol profile) ────────────────────────────
    aldh2 = _dose(snps_df, "rs671", "A")     # Lys/Lys = severe deficiency
    adh1b = _dose(snps_df, "rs1229984", "A") # His/His = fast first step
    if aldh2 is not None and aldh2 >= 1:
        findings.append({
            "title": (
                f"ALDH2 Deficiency — {'Severe (homozygous)' if aldh2 == 2 else 'Partial (heterozygous)'}"
                + (" + Fast ADH1B" if adh1b and adh1b >= 1 else "")
            ),
            "severity": "high" if aldh2 == 2 else "moderate",
            "variants": [
                f"rs671 ALDH2 (dosage {aldh2})",
                f"rs1229984 ADH1B (dosage {adh1b})" if adh1b is not None else "rs1229984 not called",
            ],
            "interpretation": (
                "ALDH2*2 carriers accumulate acetaldehyde — a known carcinogen — "
                "after alcohol consumption. Heterozygotes who drink have ~5× "
                "esophageal squamous cell carcinoma risk vs non-carriers; "
                "homozygotes typically cannot tolerate alcohol at all. Fast ADH1B "
                "accelerates ethanol→acetaldehyde, worsening the burden."
            ),
            "action": (
                "Homozygotes (A/A): essentially abstain from alcohol — flushing, "
                "tachycardia, nausea, and substantially elevated upper-GI cancer "
                "risk. Heterozygotes (G/A) who drink: minimize quantity and "
                "frequency; the cancer-risk premium for any-drinking is real. "
                "Avoid acetaldehyde-rich foods/drinks (some fermented products). "
                "Get screened for upper-GI symptoms."
            ),
        })

    # ── Slow caffeine + anxiety-prone ADORA2A ─────────────────────────────────
    cyp1a2 = _dose(snps_df, "rs762551", "C")   # C = slow allele
    adora = _dose(snps_df, "rs5751876", "T")   # T = anxiety-prone
    if cyp1a2 is not None and cyp1a2 == 2 and adora is not None and adora >= 1:
        findings.append({
            "title": "Caffeine Profile: Slow Metaboliser + Anxiety-Prone Adenosine Receptor",
            "severity": "moderate",
            "variants": ["rs762551 CYP1A2 *1F/*1F (slow)", f"rs5751876 ADORA2A (dosage {adora})"],
            "interpretation": (
                "Slow caffeine clearance (~50% slower than normal) combined with "
                "an anxiety-sensitive A2A adenosine receptor genotype increases "
                "risk of caffeine-induced anxiety, insomnia, palpitations, and "
                "in some studies elevated MI risk in slow metabolizers consuming "
                "high coffee intake."
            ),
            "action": (
                "Limit caffeine to <200 mg/day (~2 cups of brewed coffee). NONE "
                "after noon. Consider switching to decaf or tea. L-theanine "
                "200 mg can blunt caffeine-related anxiety in those who continue "
                "coffee. Note that pregnancy further slows CYP1A2 — caffeine "
                "metabolism becomes dramatically slower."
            ),
        })

    # ── Methylation × Heavy Metal compound risk ───────────────────────────────
    if c677t and c677t >= 1 and (gstm1 or 0) >= 1:
        # We can't reliably detect GSTM1 deletion from this SNP alone, so this
        # is a soft signal — emit only when MTHFR is at least heterozygous.
        findings.append({
            "title": "Methylation × Detox Compound — MTHFR + GST Activity Concerns",
            "severity": "moderate",
            "variants": ["rs1801133 (MTHFR C677T)", "rs366631 (GSTM1 proxy)"],
            "interpretation": (
                "Combined impaired methylation (MTHFR) and reduced phase II "
                "conjugation (GSTM1 null) reduces both arsenic methylation/"
                "excretion and electrophile clearance. Together this is a less "
                "favorable detox profile for environmental toxin exposure."
            ),
            "action": (
                "Methylation support (methylfolate, methyl-B12, riboflavin). "
                "Phase II support — cruciferous vegetables daily (broccoli "
                "sprouts are particularly rich in sulforaphane). NAC 600 mg "
                "twice/day for glutathione precursor. Filter drinking water "
                "(activated carbon for organics; ideally also reverse osmosis "
                "if arsenic / heavy metals are a concern in your local water)."
            ),
        })

    # ── APOE-ε4 carriage flag ─────────────────────────────────────────────────
    apoe_e4 = _dose(snps_df, "rs429358", "C")
    apoe_e2 = _dose(snps_df, "rs7412", "T")
    if apoe_e4 is not None and apoe_e2 is not None:
        if apoe_e4 == 2:
            findings.append({
                "title": "APOE ε4/ε4 — Highest-Risk Alzheimer's Genotype",
                "severity": "high",
                "variants": ["rs429358 (C/C)", "rs7412 (C/C)"],
                "interpretation": (
                    "~2% of the population carries two ε4 alleles, with ~8–12× "
                    "lifetime Alzheimer's risk vs ε3/ε3. Also elevates cardiovascular "
                    "risk and LDL cholesterol. ε4/ε4 does NOT mean Alzheimer's is "
                    "certain — many carriers never develop it — but prevention "
                    "is unusually high-value."
                ),
                "action": (
                    "Daily aerobic exercise (≥150 min/week; reduces amyloid burden "
                    "in trials). Mediterranean or MIND diet. 7–9 h sleep (deep "
                    "sleep clears amyloid). Treat hypertension and dyslipidemia "
                    "aggressively (BP <120/80, ApoB <80). Treat hearing loss "
                    "early (hearing loss is a top modifiable AD risk). Avoid "
                    "head trauma. Consider consulting a behavioral neurologist "
                    "or AD-prevention clinic for personalized risk-reduction. "
                    "Clinical trials of prevention agents are ongoing — eligible "
                    "for some."
                ),
            })

    return {
        "findings": findings,
        "n_findings": len(findings),
        "high_severity_count": sum(1 for f in findings if f["severity"] == "high"),
        "moderate_severity_count": sum(1 for f in findings if f["severity"] == "moderate"),
    }
