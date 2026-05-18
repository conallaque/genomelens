"""
HLA Imputation via Tag SNPs
===========================

HLA genes (chromosome 6p21) are the most polymorphic in the human genome.
Full-resolution HLA typing requires a reference panel (T1DGC, ~5000 phased
reference samples with both SNP and HLA typing) plus SNP2HLA/HIBAG software.

This module implements a simplified but useful approximation: for each
*clinically actionable* HLA allele, we use well-established tag SNPs from
the literature to predict carrier status. Tag SNPs are surrounding common
SNPs in tight linkage disequilibrium with the HLA allele of interest.

What we can detect well (single high-LD tag SNP):
  * HLA-B*57:01 — abacavir hypersensitivity / flucloxacillin hepatotoxicity
  * HLA-B*27 — ankylosing spondylitis, uveitis
  * HLA-B*58:01 — allopurinol SCAR
  * HLA-DRB1*15:01 — multiple sclerosis
  * HLA-DQ2.5 / DQ8 — celiac disease
  * HLA-DQB1*06:02 — narcolepsy type 1
  * HLA-C*06:02 — psoriasis
  * HLA-A*31:01 — carbamazepine SCAR

What we CANNOT reliably detect from common SNPs:
  * HLA-B*15:02 — only an Asian-specific tag exists; warn the user explicitly
  * Full 4-digit typing — needs T1DGC-quality reference panel + sequencing

This is enough for the vast majority of clinically actionable HLA findings.
"""

from __future__ import annotations

from typing import Dict, List, Optional
import pandas as pd


def _gt(snps_df: pd.DataFrame, rsid: str) -> Optional[str]:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> Optional[int]:
    gt = _gt(snps_df, rsid)
    if gt is None or len(gt) != 2:
        return None
    return gt.count(allele.upper())


# ─── Tag-SNP definitions for clinically actionable HLA alleles ───────────────
# Each allele lists one or more tag SNPs. The "primary" tag is the highest-LD
# marker. Confidence reflects literature D' or r² values.
#
# Sources: Karnes 2017, de Bakker 2006, Hetherington 2002, Karlin-Neumann
# HLA-tag literature, and CPIC documentation.

HLA_TAGS: List[Dict] = [
    {
        "allele": "HLA-B*57:01",
        "tags": [
            {"rsid": "rs2395029", "carrier_allele": "G", "ld_quality": "very high (r²>0.95)",
             "note": "HCP5 G2735C — perfect tag for B*57:01 in Europeans"},
        ],
        "frequency": "European ~6-8%, African ~3-4%, East Asian <1%",
        "clinical": [
            ("Drug — Abacavir (HIV)",
             "CONTRAINDICATED in B*57:01 carriers — severe (5-8%) and potentially fatal "
             "hypersensitivity reaction. CPIC Level A. Direct HLA typing required for confirmation "
             "before any abacavir prescription."),
            ("Drug — Flucloxacillin (antibiotic)",
             "B*57:01 carriers have ~80× risk of drug-induced liver injury. Limit use."),
            ("HIV control",
             "B*57:01 carriers are over-represented among HIV elite controllers — "
             "associated with slower progression to AIDS."),
        ],
    },
    {
        "allele": "HLA-B*27",
        "tags": [
            {"rsid": "rs4349859", "carrier_allele": "A", "ld_quality": "very high (r²>0.9)",
             "note": "Strongest tag for B*27 in Europeans"},
            {"rsid": "rs13202464", "carrier_allele": "G", "ld_quality": "high",
             "note": "Secondary B*27 tag"},
        ],
        "frequency": "European ~6-8%, Nordic peoples ~10-15%, African <1%, East Asian ~2-6%",
        "clinical": [
            ("Disease — Ankylosing spondylitis",
             "~90% of AS patients are B27+, but only ~6-8% of B27 carriers develop AS. "
             "Be aware of chronic inflammatory back pain (>3 months, age <45, morning stiffness >30 min, "
             "alternating buttock pain, improves with exercise) — see rheumatology if present."),
            ("Disease — Acute anterior uveitis",
             "Carriers have markedly elevated uveitis risk. Sudden eye pain/redness → urgent ophthalmology."),
            ("Disease — Reactive arthritis, psoriatic arthritis",
             "Modestly elevated risk; B27+ AS patients respond well to IL-17/TNF biologics."),
        ],
    },
    {
        "allele": "HLA-B*58:01",
        "tags": [
            {"rsid": "rs9263726", "carrier_allele": "T", "ld_quality": "very high (perfect tag in Asians)",
             "note": "Strong tag for B*58:01"},
        ],
        "frequency": "Han Chinese ~10-15%, Thai ~8-10%, European ~1-2%",
        "clinical": [
            ("Drug — Allopurinol (gout)",
             "Severe cutaneous adverse reactions (SJS/TEN) markedly elevated in B*58:01 carriers — "
             "Asian populations especially. CPIC Level A. If gout in B*58:01 carrier, use febuxostat or "
             "uricosurics instead."),
        ],
    },
    {
        "allele": "HLA-A*31:01",
        "tags": [
            {"rsid": "rs1633021", "carrier_allele": "T", "ld_quality": "moderate (r²~0.7)",
             "note": "European/Japanese tag for A*31:01"},
        ],
        "frequency": "European ~2-5%, Japanese ~10%, Native American ~10-15%",
        "clinical": [
            ("Drug — Carbamazepine",
             "A*31:01 carriers have elevated risk of carbamazepine-induced hypersensitivity reactions "
             "(DRESS, maculopapular eruption, SJS). CPIC Level A. Consider alternative anticonvulsants "
             "or test HLA before carbamazepine initiation."),
        ],
    },
    {
        "allele": "HLA-B*15:02",
        "tags": [],  # No reliable common-SNP tag in non-Asian populations
        "frequency": "Han Chinese ~5-10%, Thai/Malay ~10-20%, European <0.1%",
        "clinical": [
            ("Drug — Carbamazepine, oxcarbazepine, phenytoin, lamotrigine",
             "B*15:02 carriers have >1000× risk of SJS/TEN with these aromatic anticonvulsants. "
             "FDA black-box warning for Han Chinese / SE Asian patients. CPIC Level A. "
             "DIRECT HLA TYPING is required — this SNP-array proxy cannot rule out B*15:02 "
             "reliably outside of Asian populations."),
        ],
        "chip_caveat": True,
    },
    {
        "allele": "HLA-DRB1*15:01",
        "tags": [
            {"rsid": "rs3135388", "carrier_allele": "A", "ld_quality": "very high (r²>0.95)",
             "note": "Best tag for DRB1*15:01 in Europeans"},
            {"rsid": "rs9271100", "carrier_allele": "A", "ld_quality": "high"},
        ],
        "frequency": "European ~10-15%, African ~10%, East Asian <5%",
        "clinical": [
            ("Disease — Multiple sclerosis",
             "DRB1*15:01 is the strongest single MS risk factor (per-allele OR ~3, accounting "
             "for ~10% of MS heritability). Most carriers never develop MS. Vitamin D sufficiency "
             "(≥40 ng/mL) is associated with reduced MS risk. Don't smoke. Watch for episodic "
             "neurological symptoms (vision, sensation, balance)."),
            ("Disease — Narcolepsy modifier",
             "Tagged secondarily for some narcolepsy haplotypes."),
        ],
    },
    {
        "allele": "HLA-DQ2.5 (DQA1*05:01 / DQB1*02:01)",
        "tags": [
            {"rsid": "rs2187668", "carrier_allele": "T", "ld_quality": "very high (perfect tag)",
             "note": "Definitive DQ2.5 tag"},
        ],
        "frequency": "European ~25-30%, North African ~15-20%, East Asian <5%",
        "clinical": [
            ("Disease — Celiac disease",
             "~95% of celiac patients carry DQ2 or DQ8. Carriage is necessary but not sufficient — "
             "only ~1% of carriers develop celiac. If GI symptoms, fatigue, anemia, "
             "or first-degree celiac family history, test tTG-IgA antibodies on a gluten-containing diet. "
             "Being negative for both DQ2 and DQ8 essentially RULES OUT celiac."),
            ("Disease — Dermatitis herpetiformis",
             "Same celiac-spectrum association."),
        ],
    },
    {
        "allele": "HLA-DQ8 (DQA1*03:01 / DQB1*03:02)",
        "tags": [
            {"rsid": "rs7454108", "carrier_allele": "C", "ld_quality": "very high",
             "note": "Strong DQ8 tag"},
        ],
        "frequency": "European ~10%, Native American higher, African lower",
        "clinical": [
            ("Disease — Celiac disease (second major risk haplotype)",
             "~5% of celiac patients carry only DQ8 (not DQ2). Together DQ2+DQ8 cover ~99% of celiac."),
            ("Disease — Type 1 diabetes",
             "DQ8 + DQ2 (DR3/DR4 heterozygous) is the strongest T1D HLA risk profile."),
        ],
    },
    {
        "allele": "HLA-DQB1*06:02",
        "tags": [
            {"rsid": "rs2858884", "carrier_allele": "G", "ld_quality": "high",
             "note": "DQB1*06:02 tag (narcolepsy)"},
            {"rsid": "rs3104373", "carrier_allele": "C", "ld_quality": "moderate"},
        ],
        "frequency": "European ~25%, African ~30%, East Asian ~10%",
        "clinical": [
            ("Disease — Narcolepsy type 1 (with cataplexy)",
             ">98% of narcolepsy-with-cataplexy patients are DQB1*06:02 positive. However, only "
             "~0.05% of carriers develop narcolepsy — the variant is necessary but far from sufficient. "
             "Triggered by viral infections (H1N1, post-vaccination cases reported with one specific "
             "adjuvanted H1N1 vaccine). If unexplained daytime sleepiness, see sleep medicine."),
        ],
    },
    {
        "allele": "HLA-C*06:02",
        "tags": [
            {"rsid": "rs10484554", "carrier_allele": "T", "ld_quality": "very high (r²>0.9)",
             "note": "Cw6 tag for psoriasis"},
        ],
        "frequency": "European ~10-15%, lower in African/East Asian",
        "clinical": [
            ("Disease — Psoriasis vulgaris",
             "Cw6+ individuals have ~10× psoriasis risk; ~60% of psoriasis patients carry Cw6. "
             "Cw6+ patients tend to have earlier onset and respond particularly well to ustekinumab "
             "(IL-12/23 inhibitor)."),
            ("Disease — Psoriatic arthritis",
             "Modestly elevated risk in Cw6 carriers."),
        ],
    },
    {
        "allele": "HLA-B*51",
        "tags": [
            {"rsid": "rs116799036", "carrier_allele": "T", "ld_quality": "moderate",
             "note": "B*51 tag (Behcet's)"},
        ],
        "frequency": "Mediterranean / Silk Road populations highest; lower elsewhere",
        "clinical": [
            ("Disease — Behçet's syndrome",
             "B*51 is the strongest single-gene risk factor for Behçet's, a multi-system inflammatory "
             "disease with oral and genital ulcers, uveitis, skin involvement. Most B*51 carriers do not "
             "develop Behçet's."),
        ],
    },
]


# ─── Transplant compatibility context ────────────────────────────────────────
TRANSPLANT_CONTEXT = (
    "HLA matching is critical for organ transplantation (kidney, heart, liver, lung) and "
    "particularly for hematopoietic stem cell (bone marrow) transplantation. Donors and "
    "recipients are matched on HLA-A, -B, -C (class I) and HLA-DR, -DQ, -DP (class II). "
    "Closer HLA matching reduces graft rejection and graft-versus-host disease risk. "
    "Rare HLA types make finding matched unrelated donors more difficult, especially "
    "for individuals of mixed or non-European ancestry where bone marrow donor registries "
    "have lower representation. If transplantation is ever needed, full clinical HLA typing "
    "(not chip-based prediction) is required."
)


# ─── HIV / infection context ─────────────────────────────────────────────────
INFECTION_CONTEXT = [
    {"hla": "HLA-B*27, B*57:01, B*51",
     "association": "Slower HIV progression / elite control. CD8+ T-cell responses against "
                    "conserved HIV epitopes are particularly effective from these alleles."},
    {"hla": "HLA-B*35:01",
     "association": "Faster HIV progression."},
    {"hla": "HLA-B*46:01",
     "association": "Associated with worse SARS / SARS-CoV-2 (COVID-19) outcomes in some studies."},
    {"hla": "HLA-DRB1*15:01",
     "association": "Influences EBV response (relevant for MS risk pathway)."},
    {"hla": "HLA-B*15:01",
     "association": "Recently identified as the strongest known genetic factor for asymptomatic "
                    "SARS-CoV-2 infection."},
]


# ─── Core analysis ───────────────────────────────────────────────────────────

def _impute_allele(snps_df: pd.DataFrame, allele_def: Dict) -> Dict:
    """Estimate carrier status for one HLA allele using its tag SNP(s).

    Returns dict with status, dosage, confidence, called_tags.
    """
    if not allele_def["tags"]:
        return {
            "allele": allele_def["allele"],
            "status": "no_tag_available",
            "dosage": None,
            "confidence": "unable",
            "called_tags": [],
            "untested_tags": [],
        }

    called: List[Dict] = []
    untested: List[Dict] = []
    dosages: List[int] = []
    qualities: List[str] = []

    for tag in allele_def["tags"]:
        rsid = tag["rsid"]
        dose = _dose(snps_df, rsid, tag["carrier_allele"])
        if dose is None:
            untested.append(tag)
            continue
        called.append({**tag, "dosage": dose,
                       "genotype": _gt(snps_df, rsid)})
        dosages.append(dose)
        qualities.append(tag["ld_quality"])

    if not called:
        return {
            "allele": allele_def["allele"],
            "status": "untested",
            "dosage": None,
            "confidence": "unable",
            "called_tags": [],
            "untested_tags": untested,
        }

    # If multiple tag SNPs called, take maximum dosage (most informative for carrier status).
    # In reality the proper approach is consensus / haplotype phasing, but for clinically
    # actionable carrier status this approximation works.
    dosage = max(dosages)
    if dosage == 0:
        status = "negative"
    elif dosage == 1:
        status = "carrier (heterozygous)"
    else:
        status = "homozygous"

    # Confidence based on best tag quality + number of tags called
    if any("very high" in q for q in qualities):
        confidence = "high"
    elif any("high" in q for q in qualities):
        confidence = "moderate"
    else:
        confidence = "low"

    return {
        "allele": allele_def["allele"],
        "status": status,
        "dosage": dosage,
        "confidence": confidence,
        "called_tags": called,
        "untested_tags": untested,
        "frequency": allele_def.get("frequency", ""),
        "clinical": allele_def.get("clinical", []),
        "chip_caveat": allele_def.get("chip_caveat", False),
    }


def analyze_hla(snps_df: pd.DataFrame) -> Dict:
    """Impute clinically actionable HLA alleles. Returns a structured dict for
    the HTML report and downstream consumers (PGx, counseling, emergency card).
    """
    results: List[Dict] = []
    for allele_def in HLA_TAGS:
        r = _impute_allele(snps_df, allele_def)
        results.append(r)

    # Summary statistics
    n_total = len(results)
    n_called = sum(1 for r in results if r["status"] not in ("untested", "no_tag_available"))
    n_carrier = sum(1 for r in results if r["status"] in ("carrier (heterozygous)", "homozygous"))

    # Carrier alleles list for downstream consumers
    carrier_alleles = [r["allele"] for r in results
                       if r["status"] in ("carrier (heterozygous)", "homozygous")]

    return {
        "alleles": results,
        "carrier_alleles": carrier_alleles,
        "n_alleles_tested": n_total,
        "n_alleles_called": n_called,
        "n_carrier_alleles": n_carrier,
        "transplant_context": TRANSPLANT_CONTEXT,
        "infection_context": INFECTION_CONTEXT,
    }
