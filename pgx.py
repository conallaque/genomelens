"""
Pharmacogenomic Phenotyping (CPIC-style)
----------------------------------------

Genotype-to-phenotype translation using the CPIC activity-score model:

    activity = baseline + Σ (function_impact_i × dosage_i)

Each variant contributes a function impact (e.g. LOF = -1.0 per allele,
reduced-function = -0.5, increased-function = +0.5). Summed across both
chromosomes, the resulting activity score is classified into a metabolizer
phenotype, which in turn drives drug-specific recommendations grounded in
published CPIC guidelines.

Caveats (called out in the report):
  * Star-allele calling from SNP arrays alone is approximate — true
    pharmacogenomic testing also detects gene duplications/deletions
    (most relevant for CYP2D6 ultra-rapid metabolizers) and rare alleles.
  * CYP2D6 ultra-rapid status cannot be detected without CNV data.
  * Phasing (which star allele is on which chromosome) is inferred from
    activity sums and is correct for clinical purposes in nearly all cases
    but is not guaranteed.
"""

from pathlib import Path
from typing import Dict, List, Optional
import json
import pandas as pd

_DRUG_DB_PATH = Path(__file__).resolve().parent / "drug_database.json"
_DRUG_DB_CACHE: Optional[List[Dict]] = None


def _load_drug_database() -> List[Dict]:
    global _DRUG_DB_CACHE
    if _DRUG_DB_CACHE is None:
        try:
            with open(_DRUG_DB_PATH) as fh:
                data = json.load(fh)
            _DRUG_DB_CACHE = data.get("drugs", []) if isinstance(data, dict) else list(data)
        except Exception:
            _DRUG_DB_CACHE = []
    return _DRUG_DB_CACHE


def analyze_drug_database(snps_df: pd.DataFrame) -> List[Dict]:
    """For every drug entry in drug_database.json, report a finding if the
    user has a called genotype for any of the drug's snp_markers.

    Returns one finding per drug entry that has at least one matched SNP.
    """
    findings: List[Dict] = []
    drugs = _load_drug_database()
    if not drugs:
        return findings

    if snps_df is None or snps_df.empty:
        return findings
    if 'rsid' in snps_df.columns and snps_df.index.name != 'rsid':
        snps_df = snps_df.set_index('rsid')
    snps_index = snps_df.index

    for entry in drugs:
        markers = entry.get("snp_markers") or []
        matched: List[Dict] = []
        for rsid in markers:
            if rsid not in snps_index:
                continue
            gt = snps_df.loc[rsid].get("genotype")
            gt_str = str(gt).upper() if gt is not None else ""
            gt_clean = gt_str.replace(" ", "").replace("-", "")
            if gt_clean in ("", "NAN", "--"):
                continue
            matched.append({"rsid": rsid, "genotype": gt_str})

        if not matched:
            continue

        findings.append({
            "drug": entry.get("drug_name", ""),
            "genes": entry.get("genes", []) or [],
            "snp_markers": markers,
            "matched_snps": matched,
            "n_matched": len(matched),
            "n_markers": len(markers),
            "phenotypes": entry.get("phenotypes", []) or [],
            "recommendation": entry.get("dosing_recommendation", ""),
        })

    return findings


# ─── Helper ───────────────────────────────────────────────────────────────────
def _dosage(genotype: object, allele: str) -> Optional[int]:
    if genotype is None:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--"):
        return None
    if len(gt) != 2:
        return None
    return gt.count(allele.upper())


# ─── Per-gene activity-score definitions ──────────────────────────────────────
#
# Each gene config:
#   baseline_activity     — full-function diplotype score (usually 2.0)
#   variants              — list of {rsid, effect_allele, star_allele, function_impact,
#                                    name, notes}
#   phenotype_bins        — ordered list of (max_activity, phenotype, code, class)
#   drug_recs             — list of {drug, phenotype_code → recommendation}
#   cpic_guideline        — reference URL/citation string
#

GENES: Dict[str, Dict] = {

    # ── CYP2D6 ────────────────────────────────────────────────────────────────
    "CYP2D6": {
        "long_name": "Cytochrome P450 2D6",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs3892097",  "effect_allele": "A", "star_allele": "*4",  "function_impact": -1.0,
             "name": "1846G>A (splice defect)", "notes": "Most common LOF in Europeans (~20%)."},
            {"rsid": "rs1065852",  "effect_allele": "A", "star_allele": "*10", "function_impact": -0.5,
             "name": "100C>T (P34S)", "notes": "Common reduced-function allele in East Asians."},
            {"rsid": "rs28371725", "effect_allele": "T", "star_allele": "*41", "function_impact": -0.5,
             "name": "2988G>A (intron)", "notes": "Reduced-function allele present at ~9% in Europeans."},
            {"rsid": "rs16947",    "effect_allele": "A", "star_allele": "*2",  "function_impact": 0.0,
             "name": "2850C>T (R296C)", "notes": "Generally normal-function tag SNP."},
            {"rsid": "rs1135840",  "effect_allele": "G", "star_allele": "*2/*4 tag", "function_impact": 0.0,
             "name": "4180G>C (S486T)", "notes": "Haplotype-defining; not independently scored."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM)", "PM", "pheno-pm"),
            (1.0,  "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
            (10.0, "Rapid / Ultra-rapid Metabolizer", "UM", "pheno-um"),
        ],
        "um_caveat": (
            "Note: true ultra-rapid metabolizer (UM) status requires detection of "
            "CYP2D6 gene duplications, which a SNP array cannot resolve. A 'high "
            "activity' score from SNPs alone usually indicates normal metabolism."
        ),
        "drug_recs": [
            {
                "drug": "Atomoxetine (Strattera)",
                "PM": "Use 50% of standard dose. Slow titration. Monitor for cardiac/sleep side effects.",
                "IM": "Start at lower end of dose range; titrate to response.",
                "NM": "Standard dosing per label.",
                "UM": "May have reduced efficacy at standard dose; clinical monitoring matters more than dose escalation.",
            },
            {
                "drug": "Codeine, Tramadol",
                "PM": "AVOID. Pro-drugs require CYP2D6 for activation; PMs get minimal analgesia. Use morphine, hydromorphone, or non-opioid analgesics.",
                "IM": "Reduced analgesia. Consider alternatives.",
                "NM": "Standard dosing.",
                "UM": "AVOID — UMs produce dangerously high morphine levels and risk respiratory depression. Lethal in children post-tonsillectomy.",
            },
            {
                "drug": "Paroxetine, Fluoxetine, Venlafaxine (SSRIs/SNRIs)",
                "PM": "Lower starting dose or alternative (citalopram, sertraline less 2D6-dependent).",
                "IM": "Standard or modestly reduced dose; monitor side effects.",
                "NM": "Standard dosing.",
                "UM": "May need higher dose; consider non-2D6 alternative if inadequate response.",
            },
            {
                "drug": "Tamoxifen (endocrine breast cancer therapy)",
                "PM": "Major concern — PMs activate tamoxifen to endoxifen poorly. Discuss aromatase inhibitor alternative with oncologist.",
                "IM": "Possible reduced efficacy. Discuss with oncologist.",
                "NM": "Standard dosing.",
                "UM": "Standard dosing.",
            },
            {
                "drug": "Metoprolol, Carvedilol (beta-blockers)",
                "PM": "Start at lower dose; titrate cautiously (higher exposure).",
                "IM": "Start at slightly reduced dose.",
                "NM": "Standard dosing.",
                "UM": "May need higher dose for clinical response.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for CYP2D6 (Crews 2021, Bell 2017, Goetz 2018).",
    },

    # ── CYP2C9 ────────────────────────────────────────────────────────────────
    "CYP2C9": {
        "long_name": "Cytochrome P450 2C9",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs1799853", "effect_allele": "T", "star_allele": "*2", "function_impact": -0.5,
             "name": "R144C", "notes": "Reduced-function. ~12% European allele frequency."},
            {"rsid": "rs1057910", "effect_allele": "C", "star_allele": "*3", "function_impact": -1.0,
             "name": "I359L", "notes": "Severely reduced function. ~7% European frequency."},
            {"rsid": "rs9332377", "effect_allele": "G", "star_allele": "*5", "function_impact": -1.0,
             "name": "*5 variant", "notes": "Less common LOF allele."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM)", "PM", "pheno-pm"),
            (1.0,  "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Warfarin",
                "PM": "Substantial dose reduction (often 30–50% of typical). Genotype-guided dosing recommended.",
                "IM": "Moderate dose reduction. Genotype-guided dosing recommended.",
                "NM": "Standard initial dose, INR-guided.",
            },
            {
                "drug": "Phenytoin",
                "PM": "Reduce maintenance dose; toxicity risk at standard doses.",
                "IM": "Modest dose reduction; monitor levels.",
                "NM": "Standard dosing.",
            },
            {
                "drug": "NSAIDs (celecoxib, ibuprofen, naproxen)",
                "PM": "Lower starting dose for CYP2C9 substrates with narrow therapeutic index (celecoxib).",
                "IM": "Modest caution.",
                "NM": "Standard dosing.",
            },
            {
                "drug": "Sulfonylureas (glipizide, glyburide)",
                "PM": "Start at lowest dose; titrate carefully — hypoglycemia risk.",
                "IM": "Use lower starting dose.",
                "NM": "Standard dosing.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for CYP2C9 (Karnes 2021).",
    },

    # ── CYP2C19 ───────────────────────────────────────────────────────────────
    "CYP2C19": {
        "long_name": "Cytochrome P450 2C19",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs4244285",   "effect_allele": "A", "star_allele": "*2",  "function_impact": -1.0,
             "name": "681G>A (splice defect)", "notes": "Most common LOF. ~15% European frequency."},
            {"rsid": "rs12248560",  "effect_allele": "T", "star_allele": "*17", "function_impact": +0.5,
             "name": "-806C>T (promoter)", "notes": "Increased-function. ~21% European frequency."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM)", "PM", "pheno-pm"),
            (1.25, "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
            (2.75, "Rapid Metabolizer (RM)", "RM", "pheno-rm"),
            (10.0, "Ultra-rapid Metabolizer (UM)", "UM", "pheno-um"),
        ],
        "drug_recs": [
            {
                "drug": "Clopidogrel (Plavix)",
                "PM": "AVOID. Clopidogrel requires CYP2C19 activation; PMs get minimal antiplatelet effect → high stent thrombosis risk. Use prasugrel or ticagrelor.",
                "IM": "AVOID standard dose. Use prasugrel or ticagrelor unless contraindicated.",
                "NM": "Standard dosing.",
                "RM": "Standard dosing (slightly enhanced antiplatelet effect).",
                "UM": "Standard dosing; effective antiplatelet response.",
            },
            {
                "drug": "PPIs (omeprazole, esomeprazole, lansoprazole, pantoprazole)",
                "PM": "Standard or somewhat reduced dose — higher plasma levels and better acid suppression.",
                "IM": "Standard dose; good response.",
                "NM": "Standard dose.",
                "RM": "Standard dose may have reduced efficacy. Consider higher dose, twice-daily, or alternative.",
                "UM": "Standard dose often inadequate. Use higher-dose / twice-daily PPI or alternative class (e.g., H2 blocker, vonoprazan).",
            },
            {
                "drug": "Citalopram, Escitalopram (SSRIs)",
                "PM": "Reduce dose by ~50%; risk of QT prolongation at higher levels.",
                "IM": "Standard or modestly reduced dose; monitor.",
                "NM": "Standard dosing.",
                "RM": "May need higher dose or alternative SSRI.",
                "UM": "Standard dose often subtherapeutic; consider sertraline or non-2C19 antidepressant.",
            },
            {
                "drug": "Voriconazole (antifungal)",
                "PM": "Higher exposure; toxicity risk. Reduce dose or use alternative.",
                "IM": "Standard or modestly reduced dose.",
                "NM": "Standard.",
                "RM": "May need higher dose for efficacy.",
                "UM": "Standard dose often subtherapeutic; therapeutic drug monitoring recommended.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for CYP2C19 (Lee 2022).",
    },

    # ── CYP3A5 ────────────────────────────────────────────────────────────────
    "CYP3A5": {
        "long_name": "Cytochrome P450 3A5",
        "baseline_activity": 2.0,  # interpret as "expression units" rather than activity
        "variants": [
            {"rsid": "rs776746", "effect_allele": "C", "star_allele": "*3", "function_impact": -1.0,
             "name": "6986A>G (splice defect)", "notes": "Loss-of-function. ~85% European frequency."},
        ],
        "phenotype_bins": [
            (0.25, "Non-expressor (PM)", "PM", "pheno-pm"),
            (1.25, "Intermediate expressor (IM)", "IM", "pheno-im"),
            (2.25, "Expressor (NM)", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Tacrolimus (immunosuppressant)",
                "PM": "Start at standard dose (most Europeans are PMs / non-expressors).",
                "IM": "Increase initial dose ~1.5×; titrate to trough.",
                "NM": "Increase initial dose ~2×; expressors metabolize tacrolimus rapidly.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for CYP3A5 / Tacrolimus (Birdwell 2015).",
    },

    # ── TPMT ──────────────────────────────────────────────────────────────────
    "TPMT": {
        "long_name": "Thiopurine S-methyltransferase",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs1142345", "effect_allele": "C", "star_allele": "*3C", "function_impact": -1.0,
             "name": "Y240C", "notes": "Most common European LOF (~4–5% allele frequency)."},
            {"rsid": "rs1800460", "effect_allele": "T", "star_allele": "*3A/*3B", "function_impact": -1.0,
             "name": "A154T", "notes": "Co-occurs with *3C on the *3A haplotype."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM)", "PM", "pheno-pm"),
            (1.25, "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Azathioprine, 6-Mercaptopurine, Thioguanine",
                "PM": "Use ALTERNATIVE drug or 10× dose reduction with intensive CBC monitoring. Risk of life-threatening neutropenia at standard doses.",
                "IM": "Start at 30–80% of standard dose. Slow titration. Monitor CBC weekly initially.",
                "NM": "Standard starting dose. Monitor CBC.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for TPMT and NUDT15 (Relling 2019).",
    },

    # ── NUDT15 ────────────────────────────────────────────────────────────────
    "NUDT15": {
        "long_name": "Nudix hydrolase 15 (thiopurine sensitivity)",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs116855232", "effect_allele": "T", "star_allele": "*3", "function_impact": -1.0,
             "name": "R139C", "notes": "Critical LOF variant in East/South Asian populations (~10%)."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM)", "PM", "pheno-pm"),
            (1.25, "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Azathioprine, 6-Mercaptopurine, Thioguanine",
                "PM": "Use ALTERNATIVE drug or ≥10× dose reduction. Severe myelosuppression risk.",
                "IM": "Start at 30–50% of standard dose. Monitor CBC weekly.",
                "NM": "Standard dosing (in the absence of TPMT variants).",
            },
        ],
        "cpic_guideline": "CPIC Guideline for TPMT and NUDT15 (Relling 2019).",
    },

    # ── SLCO1B1 ───────────────────────────────────────────────────────────────
    "SLCO1B1": {
        "long_name": "Solute carrier organic anion transporter 1B1",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs4149056", "effect_allele": "C", "star_allele": "*5/*15", "function_impact": -1.0,
             "name": "V174A", "notes": "Decreased transporter function. ~15% European allele frequency."},
            {"rsid": "rs2306283", "effect_allele": "G", "star_allele": "*1B", "function_impact": +0.0,
             "name": "N130D", "notes": "Largely normal-function on its own; defines haplotypes with *5."},
        ],
        "phenotype_bins": [
            (0.25, "Decreased Function (PM-like)", "PM", "pheno-pm"),
            (1.25, "Intermediate Function (IM-like)", "IM", "pheno-im"),
            (2.25, "Normal Function", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Simvastatin",
                "PM": "AVOID 80 mg dose. Limit to 20 mg or switch to rosuvastatin/pravastatin/atorvastatin. Severe myopathy risk.",
                "IM": "Limit simvastatin to 40 mg/day max. Consider alternative.",
                "NM": "Standard dosing per cardiovascular indication.",
            },
            {
                "drug": "Atorvastatin, Pitavastatin",
                "PM": "Use lower dose or switch to rosuvastatin/pravastatin (less SLCO1B1-dependent).",
                "IM": "Standard or modestly reduced dose; monitor for muscle symptoms.",
                "NM": "Standard dosing.",
            },
            {
                "drug": "Rosuvastatin, Pravastatin",
                "PM": "Less SLCO1B1-dependent; usually safer choice. Still monitor for muscle symptoms.",
                "IM": "Generally well-tolerated.",
                "NM": "Standard dosing.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for Statins and SLCO1B1 (Cooper-DeHoff 2022).",
    },

    # ── VKORC1 ────────────────────────────────────────────────────────────────
    "VKORC1": {
        "long_name": "Vitamin K Epoxide Reductase Complex 1",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs9923231", "effect_allele": "A", "star_allele": "Hap A (-1639 G>A)", "function_impact": -1.0,
             "name": "-1639 G>A promoter", "notes": "A allele lowers VKORC1 expression → increased warfarin sensitivity."},
        ],
        "phenotype_bins": [
            (0.25, "High Warfarin Sensitivity", "PM", "pheno-pm"),
            (1.25, "Intermediate Warfarin Sensitivity", "IM", "pheno-im"),
            (2.25, "Normal Warfarin Sensitivity", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Warfarin",
                "PM": "Substantial dose reduction (often 50–70% of standard). Combined with CYP2C9 variants requires careful genotype-guided dosing.",
                "IM": "Modest dose reduction.",
                "NM": "Standard starting dose, INR-guided.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for Warfarin (Johnson 2017).",
    },

    # ── UGT1A1 ────────────────────────────────────────────────────────────────
    "UGT1A1": {
        "long_name": "UDP-glucuronosyltransferase 1A1",
        "baseline_activity": 2.0,
        "variants": [
            {"rsid": "rs887829", "effect_allele": "T", "star_allele": "*28 tag", "function_impact": -1.0,
             "name": "Gilbert syndrome tag", "notes": "Tags the *28 TA7 promoter repeat reducing UGT1A1 expression."},
        ],
        "phenotype_bins": [
            (0.25, "Poor Metabolizer (PM, *28/*28 — Gilbert syndrome phenotype)", "PM", "pheno-pm"),
            (1.25, "Intermediate Metabolizer (IM)", "IM", "pheno-im"),
            (2.25, "Normal Metabolizer (NM)", "NM", "pheno-nm"),
        ],
        "drug_recs": [
            {
                "drug": "Irinotecan (chemotherapy)",
                "PM": "Reduce dose 30–50%. Severe neutropenia and diarrhea risk at standard doses.",
                "IM": "Consider modest dose reduction. Monitor for toxicity.",
                "NM": "Standard dosing.",
            },
            {
                "drug": "Atazanavir (HIV protease inhibitor)",
                "PM": "Higher rates of clinically detectable jaundice (cosmetic). No dose change needed.",
                "IM": "Mild jaundice possible.",
                "NM": "Standard.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for Irinotecan and UGT1A1 (Bachtiar 2015 region).",
    },

    # ── HLA-B*57:01 (binary, not activity-based) ──────────────────────────────
    "HLA-B*57:01": {
        "long_name": "HLA-B*57:01 (abacavir hypersensitivity allele)",
        "baseline_activity": 2.0,  # unused — we override classification below
        "binary_test": True,
        "variants": [
            {"rsid": "rs763035", "effect_allele": "T", "star_allele": "B*57:01 tag", "function_impact": 0,
             "name": "HLA-B*57:01 proxy", "notes": "Tag SNP for HLA-B*57:01. Direct HLA typing is the gold standard."},
        ],
        "phenotype_bins": [],
        "drug_recs": [
            {
                "drug": "Abacavir (HIV)",
                "POS": "CONTRAINDICATED. Severe hypersensitivity risk (5–8% of carriers; can be fatal). HLA-B*57:01 status MUST be confirmed by direct HLA typing before abacavir prescription — do not rely on SNP-array proxy alone.",
                "NEG": "Standard abacavir use; routine hypersensitivity monitoring.",
            },
        ],
        "cpic_guideline": "CPIC Guideline for HLA-B and Abacavir (Martin 2014).",
    },
}


# ─── Core analysis ────────────────────────────────────────────────────────────
def _pgx_confidence(callable_variants: int, total_variants: int) -> tuple:
    """Map defining-variant coverage to an explicit confidence level.

    PGx phenotyping from a SNP array is inherently approximate (it cannot see
    gene duplications/deletions or rare alleles), so coverage maps to at most
    "high" only when every defining variant for the gene was typed.
    Returns (confidence, note).
    """
    if callable_variants == 0:
        return "none", (
            f"0 of {total_variants} defining variant(s) typed — phenotype cannot "
            "be determined and is reported as indeterminate."
        )
    pct = 100 * callable_variants / max(total_variants, 1)
    base = f"{callable_variants} of {total_variants} defining variant(s) typed"
    if pct < 50:
        return "low", (
            f"{base} (<50%). Untyped variants may carry star alleles that would "
            "change this call; treat as provisional."
        )
    if pct < 100:
        return "moderate", (
            f"{base}. Some defining variants were not typed; rare alleles on "
            "untyped sites cannot be excluded."
        )
    return "high", (
        f"{base} (full panel). Note: SNP arrays still cannot resolve gene "
        "duplications/deletions or rare alleles."
    )


def _classify(activity: float, bins: List) -> tuple:
    for max_act, label, code, cls in bins:
        if activity <= max_act:
            return label, code, cls
    return bins[-1][1], bins[-1][2], bins[-1][3]


def _analyze_gene(gene_name: str, gene_def: Dict, snps_df: pd.DataFrame) -> Dict:
    activity = gene_def.get("baseline_activity", 2.0)
    variant_calls: List[Dict] = []
    callable_variants = 0
    total_variants = len(gene_def["variants"])

    for v in gene_def["variants"]:
        rsid = v["rsid"]
        if rsid not in snps_df.index:
            variant_calls.append({**v, "genotype": None, "dosage": None,
                                  "called": False, "applied_impact": 0.0})
            continue
        gt = snps_df.loc[rsid].get("genotype")
        dose = _dosage(gt, v["effect_allele"])
        if dose is None:
            variant_calls.append({**v, "genotype": str(gt), "dosage": None,
                                  "called": False, "applied_impact": 0.0})
            continue
        applied = v["function_impact"] * dose
        activity += applied
        callable_variants += 1
        variant_calls.append({
            **v,
            "genotype": str(gt).upper(),
            "dosage": dose,
            "called": True,
            "applied_impact": applied,
        })

    activity = max(activity, 0.0)
    callability_pct = round(100 * callable_variants / max(total_variants, 1), 1)
    confidence, conf_note = _pgx_confidence(callable_variants, total_variants)

    # Binary tests (HLA-B*57:01) use a different classification.
    if gene_def.get("binary_test"):
        # A negative call requires the proxy SNP to actually be typed — absence
        # of data is NOT a negative result and must never be reported as one.
        if callable_variants == 0:
            return {
                "gene": gene_name,
                "long_name": gene_def["long_name"],
                "callable_variants": 0,
                "total_variants": total_variants,
                "callability_pct": 0.0,
                "confidence": "none",
                "confidence_note": conf_note,
                "indeterminate": True,
                "binary_result": "Indeterminate (proxy SNP not typed)",
                "phenotype": "Indeterminate — tag SNP not on this chip",
                "phenotype_code": "IND",
                "phenotype_class": "pheno-indeterminate",
                "activity_score": None,
                "variant_calls": variant_calls,
                "drug_recs": gene_def["drug_recs"],
                "cpic_guideline": gene_def["cpic_guideline"],
                "is_binary": True,
            }
        positive = any(c.get("dosage", 0) and c["dosage"] > 0 for c in variant_calls)
        return {
            "gene": gene_name,
            "long_name": gene_def["long_name"],
            "callable_variants": callable_variants,
            "total_variants": total_variants,
            "callability_pct": callability_pct,
            "confidence": confidence,
            "confidence_note": conf_note,
            "indeterminate": False,
            "binary_result": "POSITIVE (proxy)" if positive else "Negative (proxy)",
            "phenotype": "POSITIVE for HLA-B*57:01 tag" if positive else "Negative for HLA-B*57:01 tag",
            "phenotype_code": "POS" if positive else "NEG",
            "phenotype_class": "pheno-pm" if positive else "pheno-nm",
            "activity_score": None,
            "variant_calls": variant_calls,
            "drug_recs": gene_def["drug_recs"],
            "cpic_guideline": gene_def["cpic_guideline"],
            "is_binary": True,
        }

    # Activity-score genes: with NO defining variant typed the activity is just
    # the untouched baseline, which would masquerade as a confident "Normal
    # Metabolizer". Report it honestly as indeterminate instead.
    if callable_variants == 0:
        return {
            "gene": gene_name,
            "long_name": gene_def["long_name"],
            "callable_variants": 0,
            "total_variants": total_variants,
            "callability_pct": 0.0,
            "confidence": "none",
            "confidence_note": conf_note,
            "indeterminate": True,
            "activity_score": None,
            "baseline_activity": gene_def.get("baseline_activity", 2.0),
            "phenotype": "Indeterminate — no defining variants typed",
            "phenotype_code": "IND",
            "phenotype_class": "pheno-indeterminate",
            "variant_calls": variant_calls,
            "drug_recs": gene_def["drug_recs"],
            "cpic_guideline": gene_def["cpic_guideline"],
            "um_caveat": gene_def.get("um_caveat", ""),
            "is_binary": False,
        }

    label, code, cls = _classify(activity, gene_def["phenotype_bins"])
    return {
        "gene": gene_name,
        "long_name": gene_def["long_name"],
        "callable_variants": callable_variants,
        "total_variants": total_variants,
        "callability_pct": callability_pct,
        "confidence": confidence,
        "confidence_note": conf_note,
        "indeterminate": False,
        "activity_score": round(activity, 2),
        "baseline_activity": gene_def.get("baseline_activity", 2.0),
        "phenotype": label,
        "phenotype_code": code,
        "phenotype_class": cls,
        "variant_calls": variant_calls,
        "drug_recs": gene_def["drug_recs"],
        "cpic_guideline": gene_def["cpic_guideline"],
        "um_caveat": gene_def.get("um_caveat", ""),
        "is_binary": False,
    }


# ─── Novelty / exploratory gene panels (non-clinical) ────────────────────────
#
# These panels are presented in the report for curiosity only. They mix
# real HGNC-approved genes with informal/non-standard symbols that do not
# correspond to validated clinical findings. They MUST be rendered with a
# clear "non-clinical / for fun" disclaimer and never feed into actionable
# recommendations.
#
NOVELTY_PANELS: List[Dict] = [
    {
        "section": "Electro-Biological Sensitivity",
        "description": "Genetic factors affecting electrical/neurological sensitivity.",
        "genes": [
            {"gene": "ASMT",  "rsid": "rs4446909",
             "note": "Melatonin biosynthesis; promoter variation has been studied in circadian phenotypes."},
            {"gene": "NOS3",  "rsid": "rs1799983",
             "note": "Endothelial nitric oxide synthase; G894T (Glu298Asp) tweaks NO signaling."},
            {"gene": "COMT",  "rsid": "rs4680",
             "note": "Val158Met; intermediate dopamine clearance in heterozygotes (one Val158 \"warrior\" copy, one Met158 \"worrier\" copy)."},
            {"gene": "BDNF",  "rsid": "rs6265",
             "note": "Val66Met; major-allele homozygotes show standard activity-dependent BDNF secretion and neuroplasticity."},
            {"gene": "CACNA1C", "rsid": "rs1006737",
             "note": "Voltage-gated L-type calcium channel (Cav1.2); most-studied neuropsychiatric variant."},
            {"gene": "HTR2A", "rsid": "rs6313",
             "note": "Serotonin 2A receptor; T102C variant, major-allele homozygotes have baseline signaling."},
            {"gene": "CRHR2", "rsid": None,
             "note": "Corticotropin-releasing hormone receptor 2; stress-axis receptor involved in HPA-axis regulation and autonomic response. Not tested on this chip."},
            {"gene": "ISCA1", "rsid": None,
             "note": "Iron-sulfur cluster assembly 1; mitochondrial Fe-S cluster biogenesis and electron transport chain function. Not tested on this chip."},
        ],
    },
    {
        "section": "Photic / Circadian Sensitivity",
        "description": "Serotonin-melatonin biosynthesis and retinal transcription factors influencing light/circadian response.",
        "genes": [
            {"gene": "AANAT", "rsid": None,
             "note": "Aralkylamine N-acetyltransferase; rate-limiting enzyme in melatonin biosynthesis. Not tested on this chip."},
            {"gene": "ASMT",  "rsid": None,
             "note": "Acetylserotonin O-methyltransferase; final step of melatonin biosynthesis. Not tested on this chip."},
            {"gene": "TPH1",  "rsid": "rs1800532",
             "note": "Tryptophan hydroxylase 1; rate-limiting enzyme for peripheral serotonin synthesis. Major-allele homozygotes have baseline activity."},
            {"gene": "DDC",   "rsid": "rs921451",
             "note": "DOPA decarboxylase (AADC); catalyzes 5-HTP→serotonin and L-DOPA→dopamine."},
            {"gene": "CRX",   "rsid": None,
             "note": "Cone-rod homeobox transcription factor; photoreceptor differentiation and maintenance. Not tested on this chip."},
            {"gene": "OTX2",  "rsid": None,
             "note": "Orthodenticle homeobox 2; transcription factor required for retinal and brain development. Not tested on this chip."},
            {"gene": "PAX6",  "rsid": None,
             "note": "Paired box 6; master regulator transcription factor for eye and neural development. Not tested on this chip."},
        ],
    },
    {
        "section": "Pressure Sensitivity",
        "description": "Genetic factors affecting blood pressure and vascular response.",
        "genes": [
            {"gene": "MTHFR",  "rsid": "rs1801133",
             "note": "C677T; T allele reduces enzyme activity → mild homocysteine elevation."},
            {"gene": "PPARA",  "rsid": "rs1800206",
             "note": "L162V; V allele linked to altered lipid handling and vascular response."},
        ],
    },
]


def analyze_novelty_panels(snps_df: pd.DataFrame) -> List[Dict]:
    """Light-touch lookup for novelty panels. Returns one dict per section
    with per-gene genotype calls (or 'Not tested on this chip')."""
    sections = []
    for panel in NOVELTY_PANELS:
        gene_results = []
        any_called = False
        for g in panel["genes"]:
            rsid = g.get("rsid")
            if rsid and rsid in snps_df.index:
                gt = snps_df.loc[rsid].get("genotype")
                gt_str = str(gt).upper() if gt is not None and str(gt).strip() not in ("", "nan", "--") else None
                if gt_str:
                    any_called = True
                    gene_results.append({
                        "gene": g["gene"], "rsid": rsid, "tested": True,
                        "genotype": gt_str, "note": g["note"],
                    })
                    continue
            gene_results.append({
                "gene": g["gene"], "rsid": rsid, "tested": False,
                "genotype": None, "note": g["note"],
            })
        sections.append({
            "section": panel["section"],
            "description": panel["description"],
            "any_called": any_called,
            "genes": gene_results,
        })
    return sections


def analyze_pgx(snps_df: pd.DataFrame) -> Dict:
    """Run all PGx genes. Returns a dict with per-gene phenotypes and a
    consolidated list of clinically-actionable drug findings."""
    per_gene: Dict[str, Dict] = {}
    for gene_name, gene_def in GENES.items():
        per_gene[gene_name] = _analyze_gene(gene_name, gene_def, snps_df)

    # Aggregate actionable findings — phenotypes that change standard dosing
    actionable: List[Dict] = []
    actionable_codes = {"PM", "IM", "UM", "RM", "POS"}
    for gene_name, result in per_gene.items():
        if result["phenotype_code"] not in actionable_codes:
            continue
        for drug_rec in result["drug_recs"]:
            rec_for_pheno = drug_rec.get(result["phenotype_code"])
            if not rec_for_pheno:
                continue
            actionable.append({
                "gene": gene_name,
                "phenotype": result["phenotype"],
                "phenotype_class": result["phenotype_class"],
                "drug": drug_rec["drug"],
                "recommendation": rec_for_pheno,
                "guideline": result["cpic_guideline"],
            })

    try:
        novelty_panels = analyze_novelty_panels(snps_df)
    except Exception:
        novelty_panels = []

    try:
        database_findings = analyze_drug_database(snps_df)
    except Exception:
        database_findings = []

    return {
        "per_gene": per_gene,
        "actionable_findings": actionable,
        "n_genes_tested": len(per_gene),
        "n_actionable_findings": len(actionable),
        "novelty_panels": novelty_panels,
        "database_findings": database_findings,
        "n_database_findings": len(database_findings),
    }
