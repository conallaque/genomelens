"""
Personalised Supplement Stack
=============================

Generates a ranked, evidence-tagged supplement plan from the genetic profile.

Inputs (any may be None — the module gracefully degrades):
  • snps_df          — raw genotype DataFrame indexed by rsID
  • pgx_result       — output of pgx.analyze_pgx
  • wellness_result  — output of wellness.analyze_wellness
  • carrier_result   — output of carrier.analyze_carriers
  • phewas_result    — output of phewas.analyze_phewas
  • bloodwork_result — output of bloodwork.compare_bloodwork

Output:
  {
    status,
    tiers: {essential: [...], recommended: [...], optional: [...]},
    avoid: [...],
    total_estimated_cost_usd_monthly,
    n_supplements,
  }

Each supplement entry:
  {name, dose, timing, form, reasoning, snps, tier,
   monthly_cost_usd, interactions, category}
"""

from __future__ import annotations

import pandas as pd

from core import snp_registry  # V7: single source of truth for rsID metadata + strand-aware dose

# ── Low-level genotype helpers ──────────────────────────────────────────────
#
# `_risk_dose` previously lived here as a strand-aware dose helper. It now
# delegates to `snp_registry.risk_dose_from_df`, which is the canonical
# implementation used across all migrated modules. The local function is
# kept as a thin shim so external callers (and existing internal call-sites)
# don't break — but the rsID must now be in the registry, which forces every
# downstream rule to declare its variant metadata in one place.

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


def _dose(snps_df: pd.DataFrame | None, rsid: str, allele: str) -> int | None:
    gt = _gt(snps_df, rsid)
    if gt is None or len(gt) != 2:
        return None
    return gt.count(allele.upper())


def _risk_dose(
    snps_df: pd.DataFrame | None,
    rsid: str,
    risk_allele: str | None = None,
    ref_allele: str | None = None,
) -> int | None:
    """
    Strand-aware risk-allele dosage. **Delegates to the unified SNP registry.**

    Behaviour invariants (covered by tests/unit/test_supplements.py and
    tests/registry/test_snp_registry.py):

      * Returns 0/1/2 regardless of whether the chip reports on + or − strand.
      * Returns None for missing rsID, no-call ("--"), or single-char hemizygous
        Y-chromosome calls.
      * If ``risk_allele`` / ``ref_allele`` are omitted, the registry's
        canonical assignments are used (preferred — single source of truth).
    """
    return snp_registry.risk_dose_from_df(
        snps_df, rsid,
        risk_allele=risk_allele, ref_allele=ref_allele,
    )


def _chip_gap(rsid: str, gene_label: str, what_it_would_say: str) -> dict:
    """
    Returns a placeholder entry surfaced when a rule's defining SNP is not on
    the chip. Lets the user see *why* a recommendation was missing rather than
    a silent skip.
    """
    return {
        "name": f"[Not tested] {gene_label}",
        "dose": "—", "timing": "—", "form": "—",
        "category": "Chip gap",
        "tier": "chip_gap",
        "reasoning": (
            f"This rule was skipped because {rsid} ({gene_label}) is not on your "
            f"chip — chip can neither rule it in nor out. Imputation or a "
            f"different chip might cover it. {what_it_would_say}"
        ),
        "snps": [f"{rsid} ({gene_label}) — not typed"],
        "monthly_cost_usd": 0,
        "interactions": "",
    }


# ── Supplement catalogue ─────────────────────────────────────────────────────
#
# Each `_rule_*` function inspects a slice of the genetic data and returns one
# or more supplement dicts (or None / []). Rules are kept short and explicit.

CATEGORY_METHYL = "Methylation & B-vitamins"
CATEGORY_VITAMIN = "Vitamins & Minerals"
CATEGORY_HORMONE = "Hormones & Endocrine"
CATEGORY_INFLAM = "Inflammation & Recovery"
CATEGORY_DETOX = "Detox & Antioxidants"
CATEGORY_NEURO = "Neuro & Mood"
CATEGORY_FITNESS = "Performance & Recovery"
CATEGORY_GUT = "Gut & Digestion"

TIER_ESSENTIAL = "essential"
TIER_RECOMMENDED = "recommended"
TIER_OPTIONAL = "optional"


# ── Methylation / B-vitamin rules ────────────────────────────────────────────

def _rule_mthfr(snps_df) -> list[dict]:
    # MTHFR C677T (rs1801133) — T allele reduces enzyme activity 30-70% per copy.
    # Risk = T (+ strand) / A (− strand); ancestral = C (+) / G (−).
    t_dose = _risk_dose(snps_df, "rs1801133", "T", "C")
    if t_dose is None:
        return []
    if t_dose >= 2:
        return [{
            "name": "L-Methylfolate (5-MTHF)",
            "dose": "800-1000 mcg",
            "timing": "Morning with food",
            "form": "5-MTHF (e.g. Quatrefolic, Metafolin) — NOT folic acid",
            "category": CATEGORY_METHYL,
            "tier": TIER_ESSENTIAL,
            "reasoning": (
                "Homozygous MTHFR 677TT reduces methylenetetrahydrofolate reductase "
                "activity ~70%. Pre-methylated 5-MTHF bypasses the enzymatic step."
            ),
            "snps": ["rs1801133 (MTHFR C677T) TT"],
            "monthly_cost_usd": 18,
            "interactions": "Avoid pairing with folic acid (>400 mcg); some get over-methylation "
                            "irritability — pair with B6 (P5P) and methylcobalamin.",
        }]
    if t_dose == 1:
        return [{
            "name": "L-Methylfolate (5-MTHF)",
            "dose": "400 mcg",
            "timing": "Morning with food",
            "form": "5-MTHF",
            "category": CATEGORY_METHYL,
            "tier": TIER_RECOMMENDED,
            "reasoning": (
                "Heterozygous MTHFR 677CT — ~30-40% reduction in enzyme activity. "
                "Pre-methylated folate ensures consistent methylation cycle flux."
            ),
            "snps": ["rs1801133 (MTHFR C677T) CT"],
            "monthly_cost_usd": 12,
            "interactions": "Replace any folic-acid containing multivitamin.",
        }]
    return []


def _rule_mthfr_a1298c(snps_df) -> list[dict]:
    # MTHFR A1298C (rs1801131) — A→C; risk = C (+) / G (−); ancestral = A (+) / T (−).
    risk = _risk_dose(snps_df, "rs1801131", "C", "A")
    if risk is None:
        return []
    if risk >= 1:
        return [{
            "name": "Methyl-B-Complex (B6 P5P + Methylcobalamin)",
            "dose": "B6 25 mg + B12 1000 mcg",
            "timing": "Morning",
            "form": "Pyridoxal-5-phosphate (P5P) + methylcobalamin (not cyanocobalamin)",
            "category": CATEGORY_METHYL,
            "tier": TIER_RECOMMENDED,
            "reasoning": (
                "MTHFR A1298C reduces BH4 regeneration, affecting neurotransmitter "
                "and methylation pathways. Active B-vitamin forms support the cycle."
            ),
            "snps": ["rs1801131 (MTHFR A1298C)"],
            "monthly_cost_usd": 15,
            "interactions": "Methyl-donors can over-stimulate in COMT slow homozygotes — start low.",
        }]
    return []


def _rule_comt(snps_df) -> list[dict]:
    # COMT Val158Met (rs4680) — G→A; risk = A (Met, slow) on + strand or T on − strand.
    met_dose = _risk_dose(snps_df, "rs4680", "A", "G")
    val_dose = (2 - met_dose) if met_dose is not None else None
    gt = _gt(snps_df, "rs4680")
    if met_dose is None or gt is None:
        return []
    if met_dose == 2:
        return [{
            "name": "Magnesium Glycinate",
            "dose": "300-400 mg elemental",
            "timing": "Evening (calming)",
            "form": "Glycinate or threonate",
            "category": CATEGORY_NEURO,
            "tier": TIER_ESSENTIAL,
            "reasoning": (
                "COMT Met/Met homozygotes break down dopamine and catecholamines "
                "slowly — buildup amplifies stress sensitivity. Magnesium is a "
                "calming NMDA modulator and supports sleep."
            ),
            "snps": [f"rs4680 (COMT Val158Met) {gt} — Met/Met"],
            "monthly_cost_usd": 14,
            "interactions": "Avoid high-dose methyl-donors (TMG, SAM-e) in evening — can "
                            "drive anxiety in slow-COMT carriers.",
        }]
    if val_dose == 2:
        return [{
            "name": "SAM-e (optional)",
            "dose": "200-400 mg",
            "timing": "Morning, empty stomach",
            "form": "Enteric-coated S-adenosyl-methionine",
            "category": CATEGORY_NEURO,
            "tier": TIER_OPTIONAL,
            "reasoning": (
                "COMT Val/Val (fast metaboliser) clears catecholamines rapidly. "
                "SAM-e tends to be well tolerated and supports mood, methylation."
            ),
            "snps": [f"rs4680 (COMT Val158Met) {gt} — Val/Val"],
            "monthly_cost_usd": 28,
            "interactions": "Do not combine with SSRIs/MAOIs without physician.",
        }]
    return []


# ── Vitamin metabolism rules ────────────────────────────────────────────────

def _rule_vitamin_d(snps_df, wellness, phewas) -> list[dict]:
    cyp2r1 = _dose(snps_df, "rs10741657", "G")
    gc = _dose(snps_df, "rs2282679", "C")
    vdr_taq = _gt(snps_df, "rs731236")
    vdr_fok = _gt(snps_df, "rs2228570")

    score = (cyp2r1 or 0) + (gc or 0)
    vdr_active = vdr_fok and "T" in vdr_fok
    snps_used: list[str] = []
    if cyp2r1 is not None:
        snps_used.append(f"rs10741657 (CYP2R1) dose {cyp2r1}")
    if gc is not None:
        snps_used.append(f"rs2282679 (GC/VDBP) dose {gc}")
    if vdr_fok:
        snps_used.append(f"rs2228570 (VDR FokI) {vdr_fok}")
    if vdr_taq:
        snps_used.append(f"rs731236 (VDR TaqI) {vdr_taq}")

    if score >= 3 or vdr_active:
        dose = "4000 IU"
        tier = TIER_ESSENTIAL
        reasoning = (
            "Multiple low-D genotypes (CYP2R1 + GC/VDBP) plus/or active VDR FokI — "
            "expect persistently low serum 25(OH)D unless dosing aggressively."
        )
    elif score >= 2:
        dose = "2000-3000 IU"
        tier = TIER_RECOMMENDED
        reasoning = "Moderately reduced vitamin D status genotype."
    elif score >= 1:
        dose = "1000-2000 IU"
        tier = TIER_RECOMMENDED
        reasoning = "Mildly reduced vitamin D synthesis or transport."
    else:
        dose = "1000 IU"
        tier = TIER_OPTIONAL
        reasoning = "Standard maintenance; favourable VDR/CYP2R1/GC genotype."

    return [{
        "name": "Vitamin D3 + K2 (MK-7)",
        "dose": f"{dose} D3 / 100-200 mcg K2",
        "timing": "With fattiest meal of the day",
        "form": "Cholecalciferol (D3) with menaquinone-7",
        "category": CATEGORY_VITAMIN,
        "tier": tier,
        "reasoning": reasoning,
        "snps": snps_used,
        "monthly_cost_usd": 12,
        "interactions": "Co-supplement magnesium (D activation requires Mg); avoid high-dose D "
                        "without K2 (calcium routing).",
    }]


def _rule_b12(snps_df) -> list[dict]:
    fut2 = _gt(snps_df, "rs602662")
    tcn2 = _gt(snps_df, "rs1801198")
    needs = []
    if fut2 and "A" in fut2:
        needs.append(f"rs602662 (FUT2) {fut2} → reduced absorption")
    if tcn2 and "G" in tcn2 and tcn2.count("G") >= 1:
        needs.append(f"rs1801198 (TCN2) {tcn2} → reduced transport")
    if not needs:
        return []
    return [{
        "name": "Methylcobalamin B12",
        "dose": "1000 mcg sublingual",
        "timing": "Morning",
        "form": "Methylcobalamin (avoid cyanocobalamin)",
        "category": CATEGORY_METHYL,
        "tier": TIER_RECOMMENDED,
        "reasoning": (
            "Genotype predicts reduced B12 status due to absorption (FUT2) and/or "
            "transport (TCN2) variants. Sublingual methyl-form bypasses both."
        ),
        "snps": needs,
        "monthly_cost_usd": 10,
        "interactions": "Vegans and PPI/metformin users need higher doses; recheck serum B12 yearly.",
    }]


def _rule_iron(snps_df, carrier) -> list[dict]:
    # HFE C282Y / H63D carriers should NOT routinely supplement iron
    if carrier:
        for entry in carrier.get("affected", []) + carrier.get("carriers", []):
            if entry.get("gene") == "HFE":
                return [{
                    "name": "AVOID iron supplements",
                    "dose": "—",
                    "timing": "—",
                    "form": "—",
                    "category": CATEGORY_VITAMIN,
                    "tier": "avoid",
                    "reasoning": (
                        f"HFE variant ({entry.get('variant', 'C282Y/H63D')}) carrier — "
                        "iron supplementation can drive overload. Verify ferritin first "
                        "and avoid multivitamins containing iron."
                    ),
                    "snps": [f"{entry['rsid']} ({entry.get('variant','HFE')})"],
                    "monthly_cost_usd": 0,
                    "interactions": "Skip standard 'iron' multis; choose iron-free formulations.",
                }]
    # TMPRSS6 / TFR2 low-iron genotype → opportunity supplement
    tmprss6 = _dose(snps_df, "rs855791", "A")
    if tmprss6 is not None and tmprss6 >= 1:
        return [{
            "name": "Gentle Iron (bisglycinate)",
            "dose": "18-25 mg elemental",
            "timing": "Empty stomach with vitamin C, away from coffee/tea",
            "form": "Ferrous bisglycinate (lower GI side-effects)",
            "category": CATEGORY_VITAMIN,
            "tier": TIER_OPTIONAL,
            "reasoning": (
                "TMPRSS6 variant associated with mildly lower hemoglobin / MCV. Only "
                "supplement if ferritin <50 ng/mL or symptomatic — test first."
            ),
            "snps": [f"rs855791 (TMPRSS6) dose {tmprss6}"],
            "monthly_cost_usd": 8,
            "interactions": "Calcium, magnesium, zinc compete for absorption — separate by 2h.",
        }]
    return []


def _rule_omega3(snps_df) -> list[dict]:
    # FADS1 rs174547 — T allele = slower conversion of ALA → EPA/DHA.
    # Risk = T (+) / A (−); ancestral = C (+) / G (−).
    t_dose = _risk_dose(snps_df, "rs174547", "T", "C")
    gt = _gt(snps_df, "rs174547")
    if t_dose is None or gt is None:
        return []
    if t_dose >= 1:
        dose = "2 g EPA+DHA" if t_dose == 2 else "1.5 g EPA+DHA"
        return [{
            "name": "EPA/DHA (fish or algal oil)",
            "dose": dose,
            "timing": "With largest meal",
            "form": "Triglyceride-form fish oil or algal oil (vegetarian)",
            "category": CATEGORY_INFLAM,
            "tier": TIER_ESSENTIAL if t_dose == 2 else TIER_RECOMMENDED,
            "reasoning": (
                "FADS1 reduced-conversion genotype — endogenous ALA→EPA/DHA pathway "
                "is impaired. Direct preformed EPA/DHA supplementation needed."
            ),
            "snps": [f"rs174547 (FADS1) {gt} — T-allele dose {t_dose}"],
            "monthly_cost_usd": 22,
            "interactions": "Mild anticoagulant; check with physician if on warfarin/DOACs.",
        }]
    return []


# ── Hormones / inflammation rules ───────────────────────────────────────────

def _rule_zinc_dhea(snps_df) -> list[dict]:
    # SRD5A2 / AR / SHBG context — pragmatic: low SHBG genotypes often benefit
    # from zinc + boron support
    out = []
    shbg = _gt(snps_df, "rs1799941")
    if shbg and shbg.count("A") >= 1:
        out.append({
            "name": "Zinc Picolinate",
            "dose": "15 mg",
            "timing": "Evening with food",
            "form": "Picolinate or bisglycinate",
            "category": CATEGORY_HORMONE,
            "tier": TIER_OPTIONAL,
            "reasoning": (
                "SHBG genotype predicts higher SHBG → lower free testosterone. "
                "Zinc supports endogenous T production and aromatase modulation."
            ),
            "snps": [f"rs1799941 (SHBG) {shbg}"],
            "monthly_cost_usd": 7,
            "interactions": "Long-term high zinc depletes copper — keep dose ≤25 mg/day.",
        })
    return out


def _rule_curcumin(snps_df, phewas) -> list[dict]:
    # IL6 -174 G/G (rs1800795 CC on - strand reported as GG) → high inflammation
    il6 = _gt(snps_df, "rs1800795")
    crp_pred = None
    if phewas and "C-Reactive Protein" in phewas.get("traits", {}):
        crp_pred = phewas["traits"]["C-Reactive Protein"]["result"].get("tier")
    high_inflam = (il6 and "G" in il6) or crp_pred in ("Very high", "Above average")
    if not high_inflam:
        return []
    return [{
        "name": "Curcumin (with piperine or phytosome)",
        "dose": "500-1000 mg curcuminoids",
        "timing": "With meals 1-2× daily",
        "form": "Meriva / Longvida / BCM-95 or curcumin + 5-10 mg piperine",
        "category": CATEGORY_INFLAM,
        "tier": TIER_RECOMMENDED,
        "reasoning": (
            "IL6 promoter variant and/or elevated genetic CRP prediction — chronic "
            "low-grade inflammation tendency. Curcumin inhibits NF-κB and lowers CRP."
        ),
        "snps": [f"rs1800795 (IL6) {il6}" if il6 else "phewas-CRP elevated"],
        "monthly_cost_usd": 20,
        "interactions": "Mild antiplatelet; pause 1 week pre-surgery; can interact with tacrolimus.",
    }]


def _rule_glutathione(snps_df, tier1) -> list[dict]:
    # GSTM1/GSTT1 null and GSTP1 rs1695 — reduced detox capacity
    gstp = _gt(snps_df, "rs1695")
    nrf2 = _gt(snps_df, "rs6721961")
    if not (gstp and "G" in gstp) and not (nrf2 and "T" in nrf2):
        return []
    return [{
        "name": "N-Acetyl Cysteine (NAC)",
        "dose": "600-1200 mg",
        "timing": "Empty stomach 1-2× daily",
        "form": "NAC capsule",
        "category": CATEGORY_DETOX,
        "tier": TIER_RECOMMENDED,
        "reasoning": (
            "Reduced glutathione-system capacity (GSTP1 / NRF2 variant). NAC is a "
            "direct cysteine donor for glutathione synthesis — supports detox and "
            "antioxidant defense."
        ),
        "snps": [f"rs1695 (GSTP1) {gstp}" if gstp else "",
                 f"rs6721961 (NRF2) {nrf2}" if nrf2 else ""],
        "monthly_cost_usd": 15,
        "interactions": "Smells of sulfur; can blunt nitroglycerin response.",
    }]


# ── Drug-metabolism awareness ───────────────────────────────────────────────

def _rule_caffeine(snps_df) -> list[dict]:
    # CYP1A2 rs762551 — *1F (C) reduces inducibility → slow metaboliser.
    # Risk = C (+) / G (−); ancestral = A (+) / T (−).
    c_dose = _risk_dose(snps_df, "rs762551", "C", "A")
    gt = _gt(snps_df, "rs762551")
    if c_dose is None or gt is None or c_dose == 0:
        return []
    return [{
        "name": "L-Theanine",
        "dose": "100-200 mg with caffeine",
        "timing": "Paired with morning coffee",
        "form": "Suntheanine or generic L-theanine",
        "category": CATEGORY_NEURO,
        "tier": TIER_OPTIONAL,
        "reasoning": (
            "CYP1A2 slow-metaboliser genotype (rs762551 *1F) — caffeine half-life "
            "prolonged; L-theanine blunts jitteriness without reducing alertness."
        ),
        "snps": [f"rs762551 (CYP1A2) {gt}"],
        "monthly_cost_usd": 9,
        "interactions": "Generally well tolerated; mildly hypotensive in large doses.",
    }]


def _rule_creatine(snps_df, wellness) -> list[dict]:
    # ACTN3 R577X (rs1815739) — C (+) / G (−) is R-allele (power); T (+) / A (−) is X.
    r_dose = _risk_dose(snps_df, "rs1815739", "C", "T")
    actn3 = _gt(snps_df, "rs1815739")
    if r_dose is None or actn3 is None:
        return []
    return [{
        "name": "Creatine Monohydrate",
        "dose": "3-5 g daily",
        "timing": "Any time, consistency > timing",
        "form": "Micronised monohydrate (skip exotic forms)",
        "category": CATEGORY_FITNESS,
        "tier": TIER_RECOMMENDED if r_dose >= 1 else TIER_OPTIONAL,
        "reasoning": (
            "ACTN3 R-allele carriers respond particularly well to creatine for "
            "power output and muscle volume. Universal cognitive benefits apply "
            "across genotypes."
        ),
        "snps": [f"rs1815739 (ACTN3) {actn3}"],
        "monthly_cost_usd": 8,
        "interactions": "Very safe; transient water-weight gain; hydrate well.",
    }]


# ── Aggregator ──────────────────────────────────────────────────────────────

_TIER_ORDER = {
    TIER_ESSENTIAL: 0, TIER_RECOMMENDED: 1, TIER_OPTIONAL: 2,
    "avoid": 3, "chip_gap": 4,
}

# Rules whose absence from a chip should be surfaced explicitly to the user.
# Each entry: (rsid, gene_label, what_it_would_inform).
_KEY_RULE_SNPS: list[tuple] = [
    ("rs1801133", "MTHFR C677T",       "Methylfolate vs folic-acid recommendation."),
    ("rs1801131", "MTHFR A1298C",      "Methyl-B-complex recommendation."),
    ("rs4680",    "COMT Val158Met",    "Magnesium / SAM-e mood support."),
    ("rs10741657","CYP2R1 (vit D)",    "Vitamin D dose tier."),
    ("rs2282679", "GC/VDBP (vit D)",   "Vitamin D dose tier."),
    ("rs2228570", "VDR FokI",          "Vitamin D dose tier."),
    ("rs602662",  "FUT2 (B12 absorp)", "Methylcobalamin recommendation."),
    ("rs1801198", "TCN2 (B12 transp)", "Methylcobalamin recommendation."),
    ("rs855791",  "TMPRSS6 (iron)",    "Gentle-iron recommendation."),
    ("rs174547",  "FADS1 (omega-3)",   "EPA/DHA dose."),
    ("rs1799941", "SHBG",              "Zinc/free-T recommendation."),
    ("rs1800795", "IL6 promoter",      "Curcumin / anti-inflammation."),
    ("rs1695",    "GSTP1 (detox)",     "NAC recommendation."),
    ("rs6721961", "NRF2",              "NAC recommendation."),
    ("rs762551",  "CYP1A2 (caffeine)", "L-theanine pairing."),
    ("rs1815739", "ACTN3 (power)",     "Creatine tier."),
]


def build_supplement_stack(
    snps_df: pd.DataFrame | None = None,
    pgx_result: dict | None = None,
    wellness_result: dict | None = None,
    carrier_result: dict | None = None,
    phewas_result: dict | None = None,
    tier1_results: list[dict] | None = None,
    bloodwork_result: dict | None = None,
) -> dict:
    """Run every rule and consolidate into a tiered, deduplicated stack."""
    all_recs: list[dict] = []
    chip_gaps: list[dict] = []
    if snps_df is not None:
        all_recs.extend(_rule_mthfr(snps_df))
        all_recs.extend(_rule_mthfr_a1298c(snps_df))
        all_recs.extend(_rule_comt(snps_df))
        all_recs.extend(_rule_vitamin_d(snps_df, wellness_result, phewas_result))
        all_recs.extend(_rule_b12(snps_df))
        all_recs.extend(_rule_iron(snps_df, carrier_result))
        all_recs.extend(_rule_omega3(snps_df))
        all_recs.extend(_rule_zinc_dhea(snps_df))
        all_recs.extend(_rule_curcumin(snps_df, phewas_result))
        all_recs.extend(_rule_glutathione(snps_df, tier1_results))
        all_recs.extend(_rule_caffeine(snps_df))
        all_recs.extend(_rule_creatine(snps_df, wellness_result))

        # Surface chip gaps for rules whose defining SNP is missing — gives
        # the user an auditable "we checked but couldn't see this" rather
        # than a silent skip indistinguishable from "checked and ancestral".
        for rsid, label, hint in _KEY_RULE_SNPS:
            if _gt(snps_df, rsid) is None:
                chip_gaps.append(_chip_gap(rsid, label, hint))

    # Deduplicate by name, keeping the highest-priority entry
    by_name: dict[str, dict] = {}
    for rec in all_recs:
        key = rec["name"]
        if key not in by_name or _TIER_ORDER[rec["tier"]] < _TIER_ORDER[by_name[key]["tier"]]:
            by_name[key] = rec
    final = list(by_name.values())

    # Sort within each tier by category
    final.sort(key=lambda r: (_TIER_ORDER[r["tier"]], r["category"], r["name"]))

    tiers: dict[str, list[dict]] = {
        TIER_ESSENTIAL: [], TIER_RECOMMENDED: [], TIER_OPTIONAL: [], "avoid": [],
        "chip_gap": chip_gaps,
    }
    for r in final:
        tiers[r["tier"]].append(r)

    # Bloodwork-driven refinements: if user actually has high CRP, bump curcumin up
    if bloodwork_result and bloodwork_result.get("status") == "ok":
        for row in bloodwork_result.get("rows", []):
            if row["trait"] == "C-Reactive Protein" and row["actual_tier"] in ("Very high", "Above average"):
                for r in final:
                    if r["name"].startswith("Curcumin") and r["tier"] != TIER_ESSENTIAL:
                        tiers[r["tier"]].remove(r)
                        r["tier"] = TIER_ESSENTIAL
                        r["reasoning"] += " [Actual CRP confirms genetic prediction.]"
                        tiers[TIER_ESSENTIAL].append(r)
            if row["trait"] == "25-OH Vitamin D" and row["actual_tier"] in ("Very low", "Below average"):
                for r in final:
                    if r["name"].startswith("Vitamin D3"):
                        r["reasoning"] += " [Actual 25(OH)D confirms low.]"

    monthly_cost = sum(r.get("monthly_cost_usd", 0) for r in final if r["tier"] != "avoid")

    return {
        "status": "ok" if final or chip_gaps else "no_data",
        "tiers": tiers,
        "avoid": tiers["avoid"],
        "chip_gaps": chip_gaps,
        "n_supplements": sum(1 for r in final if r["tier"] != "avoid"),
        "n_essential": len(tiers[TIER_ESSENTIAL]),
        "n_recommended": len(tiers[TIER_RECOMMENDED]),
        "n_optional": len(tiers[TIER_OPTIONAL]),
        "n_chip_gaps": len(chip_gaps),
        "total_estimated_cost_usd_monthly": monthly_cost,
        "notes": (
            "Genetic recommendations are starting points — confirm with labs and a "
            "physician. Doses are typical adult ranges; pregnant/breastfeeding/under-18 "
            "users should not self-supplement without clinical guidance."
        ),
    }


# ── HTML rendering ───────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_SUPP_CSS = """
<style>
.sp-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           color:#222; max-width: 1100px; margin: 24px auto; padding: 0 16px; }
.sp-wrap h1 { font-size: 1.6em; border-bottom: 2px solid #333; padding-bottom: 6px; }
.sp-wrap h2 { font-size: 1.2em; margin-top: 28px; padding-bottom:4px;
              border-bottom: 1px solid #eee; }
.sp-summary { display:flex; gap:14px; flex-wrap:wrap; margin: 18px 0; }
.sp-stat { background:#f6f6f7; border:1px solid #ddd; border-radius:8px;
           padding:10px 14px; min-width:120px; }
.sp-stat .v { font-size:1.5em; font-weight:600; }
.sp-stat .l { font-size:0.8em; color:#666; text-transform:uppercase; letter-spacing:0.05em; }
.sp-card { background:#fcfcfd; border:1px solid #e2e2e6; border-radius:10px;
           padding:14px 16px; margin:10px 0; }
.sp-card.sp-essential { border-left: 5px solid #2c7a30; }
.sp-card.sp-recommended { border-left: 5px solid #b48a00; }
.sp-card.sp-optional { border-left: 5px solid #5a6772; }
.sp-card.sp-avoid { border-left: 5px solid #a32a2a; background:#fdf2f2; }
.sp-card.sp-chip_gap { border-left: 5px solid #888; background:#f4f4f6;
                       color:#555; font-style:italic; }
.sp-name { font-size:1.1em; font-weight:600; }
.sp-meta { color:#555; font-size:0.9em; margin-top:4px; }
.sp-snps { font-size:0.85em; color:#666; font-family:Menlo, monospace; margin-top:6px; }
.sp-reason { margin-top:6px; }
.sp-inter { font-size:0.85em; color:#883333; margin-top:6px; }
</style>
"""


def _render_card(r: dict) -> str:
    snps = [s for s in r.get("snps", []) if s]
    return f"""
<div class="sp-card sp-{r['tier']}">
  <div class="sp-name">{_esc(r['name'])}
    <span style="color:#888;font-weight:400;font-size:0.85em">— {_esc(r['category'])}</span></div>
  <div class="sp-meta"><strong>{_esc(r['dose'])}</strong> · {_esc(r['timing'])} · {_esc(r['form'])}
    · <strong>~${r.get('monthly_cost_usd', 0)}/mo</strong></div>
  <div class="sp-reason">{_esc(r['reasoning'])}</div>
  {f'<div class="sp-snps">SNPs: {"; ".join(_esc(s) for s in snps)}</div>' if snps else ''}
  <div class="sp-inter">⚠ {_esc(r['interactions'])}</div>
</div>
""".strip()


def render_supplements_html(result: dict, file_label: str = "") -> str:
    if not result or result.get("status") != "ok":
        body = "<p>No supplement recommendations could be generated (insufficient genetic data).</p>"
    else:
        sections = []
        for tier_key, tier_label in [
            (TIER_ESSENTIAL, "Essential"),
            (TIER_RECOMMENDED, "Recommended"),
            (TIER_OPTIONAL, "Optional"),
        ]:
            items = result["tiers"].get(tier_key, [])
            if not items:
                continue
            sections.append(f"<h2>{tier_label} ({len(items)})</h2>")
            sections.extend(_render_card(r) for r in items)
        if result.get("avoid"):
            sections.append(f"<h2>Avoid ({len(result['avoid'])})</h2>")
            sections.extend(_render_card(r) for r in result["avoid"])
        gaps = result.get("chip_gaps") or result["tiers"].get("chip_gap", [])
        if gaps:
            sections.append(
                f"<h2 style='color:#666'>Chip gaps — couldn't evaluate ({len(gaps)})</h2>"
                f"<p style='color:#666;font-size:0.9em'>These SNPs are not on your "
                f"chip, so the corresponding recommendation could not be checked. "
                f"Imputation (<code>--impute</code>) or a different chip may cover them.</p>"
            )
            sections.extend(_render_card(r) for r in gaps)
        body = f"""
          <div class="sp-summary">
            <div class="sp-stat"><div class="v">{result['n_supplements']}</div>
              <div class="l">Total</div></div>
            <div class="sp-stat"><div class="v" style="color:#2c7a30">{result['n_essential']}</div>
              <div class="l">Essential</div></div>
            <div class="sp-stat"><div class="v" style="color:#b48a00">{result['n_recommended']}</div>
              <div class="l">Recommended</div></div>
            <div class="sp-stat"><div class="v" style="color:#5a6772">{result['n_optional']}</div>
              <div class="l">Optional</div></div>
            <div class="sp-stat"><div class="v">~${result['total_estimated_cost_usd_monthly']}</div>
              <div class="l">Est. monthly</div></div>
            <div class="sp-stat"><div class="v" style="color:#888">{result.get('n_chip_gaps', 0)}</div>
              <div class="l">Chip gaps</div></div>
          </div>
          {''.join(sections)}
        """
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Personalised Supplement Stack{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_SUPP_CSS}</head><body><div class="sp-wrap">
<h1>Personalised Supplement Stack</h1>
<p style="color:#666">{_esc(result.get('notes',''))}</p>
{body}
<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. Doses and timing are general adult ranges from peer-reviewed
genetic-nutrition literature; individual response varies. Confirm with serum
labs before sustained high-dose supplementation, and discuss with a physician
if pregnant, breastfeeding, under 18, or taking prescription medications.
</p>
</div></body></html>"""
