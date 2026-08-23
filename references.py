"""
References & Evidence Catalog
-----------------------------

Maps key variants to:
  * PubMed IDs (PMIDs) for primary GWAS / clinical studies
  * Evidence levels (CPIC A/B/C/D for PGx; GWAS / replicated meta-analysis;
    ClinGen-style assertion strength)
  * Source consortium / guideline body

The report renders both:
  * Per-variant evidence chips (inline)
  * A consolidated References appendix grouped by category
"""


# CPIC evidence levels (A = strong, B = moderate, C = optional, D = informational)
# GWAS levels: GWAS-A = genome-wide significant + replicated; GWAS-B = GWS + meta-analysis
# Clinical levels: ClinVar Path/Likely Path; ACMG; OMIM

REFERENCES: dict[str, dict] = {
    # ── APOE / Alzheimer's ───────────────────────────────────────────────────
    "rs429358": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "APOE ε4 — the strongest common-variant Alzheimer's risk factor.",
        "pmids": ["8346443", "23300843", "30820047"],
        "guidelines": ["NIA-AA recommendations", "IGAP 2019"],
        "clinvar": "drug response, risk factor",
    },
    "rs7412": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "APOE ε2 — protective for Alzheimer's, hyperlipidemia risk for E2/E2.",
        "pmids": ["8346443", "30820047"],
        "guidelines": ["NIA-AA recommendations"],
    },
    "rs75932628": {
        "evidence_level": "Clinical-Validated (rare variant)",
        "evidence_summary": "TREM2 R47H — ~3× Alzheimer's risk per copy, similar magnitude to APOE-ε4.",
        "pmids": ["23150934", "23150908"],
    },

    # ── Hemochromatosis ──────────────────────────────────────────────────────
    "rs1800562": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "HFE C282Y — defining variant for hereditary hemochromatosis type 1.",
        "pmids": ["8696333", "29754290"],
        "guidelines": ["AASLD Hereditary Hemochromatosis guidelines"],
        "clinvar": "Pathogenic (recessive)",
    },
    "rs1799945": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "HFE H63D — mild iron-loading risk; relevant compound het with C282Y.",
        "pmids": ["8696333"],
        "guidelines": ["AASLD"],
    },

    # ── Thrombophilia ────────────────────────────────────────────────────────
    "rs6025": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "Factor V Leiden — primary inherited thrombophilia.",
        "pmids": ["8164741", "20308722"],
        "guidelines": ["ACMG / ACOG"],
        "clinvar": "Risk factor",
    },
    "rs1799963": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "Prothrombin G20210A — second most common inherited thrombophilia.",
        "pmids": ["8916933"],
        "guidelines": ["ACMG / ACOG"],
    },

    # ── Methylation ──────────────────────────────────────────────────────────
    "rs1801133": {
        "evidence_level": "Population-Validated",
        "evidence_summary": "MTHFR C677T — 70% activity loss in T/T; elevated homocysteine.",
        "pmids": ["7647779", "12068374"],
    },
    "rs1801131": {
        "evidence_level": "Population-Validated",
        "evidence_summary": "MTHFR A1298C — moderate enzyme impairment; relevant in compound state.",
        "pmids": ["10416283"],
    },

    # ── PGx CYP2D6 ────────────────────────────────────────────────────────────
    "rs3892097": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2D6 *4 (1846G>A) — most common European LOF; affects codeine, tramadol, atomoxetine, tamoxifen, antidepressants.",
        "pmids": ["32189290", "29385227"],
        "guidelines": ["CPIC Guideline for CYP2D6 (Crews 2021)"],
    },
    "rs1065852": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2D6 *10 (100C>T) — reduced function, common in East Asians.",
        "pmids": ["32189290"],
        "guidelines": ["CPIC Guideline for CYP2D6"],
    },
    "rs28371725": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2D6 *41 — reduced-function allele.",
        "pmids": ["32189290"],
        "guidelines": ["CPIC"],
    },
    "rs16947": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2D6 *2 family tag.",
        "pmids": ["32189290"],
    },

    # ── CYP2C9 ────────────────────────────────────────────────────────────────
    "rs1799853": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2C9 *2 — reduced function; warfarin, NSAIDs, phenytoin.",
        "pmids": ["29193033"],
        "guidelines": ["CPIC Guideline for CYP2C9 / Warfarin"],
    },
    "rs1057910": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2C9 *3 — severely reduced function; warfarin sensitivity.",
        "pmids": ["29193033"],
        "guidelines": ["CPIC"],
    },

    # ── CYP2C19 ───────────────────────────────────────────────────────────────
    "rs4244285": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2C19 *2 — primary LOF, critical for clopidogrel, PPIs, citalopram.",
        "pmids": ["23695185", "34159645"],
        "guidelines": ["CPIC Guideline for CYP2C19 / Clopidogrel"],
    },
    "rs12248560": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP2C19 *17 — increased-function allele.",
        "pmids": ["23695185"],
        "guidelines": ["CPIC"],
    },

    # ── TPMT / NUDT15 ─────────────────────────────────────────────────────────
    "rs1142345": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "TPMT *3C — LOF; severe thiopurine toxicity if homozygous.",
        "pmids": ["30447069"],
        "guidelines": ["CPIC Guideline for TPMT/NUDT15"],
    },
    "rs1800460": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "TPMT *3A/*3B variant.",
        "pmids": ["30447069"],
        "guidelines": ["CPIC"],
    },
    "rs116855232": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "NUDT15 *3 — critical thiopurine sensitivity variant in East Asians.",
        "pmids": ["30447069"],
        "guidelines": ["CPIC"],
    },

    # ── SLCO1B1 ───────────────────────────────────────────────────────────────
    "rs4149056": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "SLCO1B1 *5 — simvastatin-induced myopathy risk.",
        "pmids": ["35034351"],
        "guidelines": ["CPIC Guideline for Statins"],
    },

    # ── VKORC1 ────────────────────────────────────────────────────────────────
    "rs9923231": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "VKORC1 -1639G>A — major warfarin dose determinant.",
        "pmids": ["28198005"],
        "guidelines": ["CPIC Guideline for Warfarin"],
    },

    # ── CYP3A5 ────────────────────────────────────────────────────────────────
    "rs776746": {
        "evidence_level": "CPIC Level A",
        "evidence_summary": "CYP3A5 *3 — non-expressor; tacrolimus dosing.",
        "pmids": ["25801146"],
        "guidelines": ["CPIC Guideline for Tacrolimus"],
    },

    # ── HLA-B*57:01 ───────────────────────────────────────────────────────────
    "rs763035": {
        "evidence_level": "CPIC Level A (verify by HLA typing)",
        "evidence_summary": "HLA-B*57:01 proxy — abacavir hypersensitivity contraindication.",
        "pmids": ["24447389"],
        "guidelines": ["CPIC Guideline for Abacavir"],
    },

    # ── CAD ──────────────────────────────────────────────────────────────────
    "rs10455872": {
        "evidence_level": "GWAS-A (replicated meta-analysis)",
        "evidence_summary": "LPA — strongest common variant for elevated Lp(a) and CAD.",
        "pmids": ["20032323", "33892491"],
    },
    "rs3798220": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "LPA I4399M — large effect on Lp(a) levels.",
        "pmids": ["20032323"],
    },
    "rs10757278": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "9p21 — strongest non-lipid CAD locus.",
        "pmids": ["17478679", "26343387"],
    },
    "rs1333049": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "9p21 — replicated CAD locus.",
        "pmids": ["17478679"],
    },
    "rs11591147": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "PCSK9 R46L — LOF lowering LDL; protective.",
        "pmids": ["16554528"],
    },
    "rs6511720": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "LDLR — protective for LDL/CAD.",
        "pmids": ["20660468"],
    },

    # ── T2D ──────────────────────────────────────────────────────────────────
    "rs7903146": {
        "evidence_level": "GWAS-A (strongest T2D variant)",
        "evidence_summary": "TCF7L2 — per-allele OR ~1.4 for T2D.",
        "pmids": ["16415884", "24509480"],
    },
    "rs9939609": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "FTO — strongest common variant for BMI / obesity.",
        "pmids": ["17434869"],
    },
    "rs10830963": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "MTNR1B — fasting glucose, T2D, gestational diabetes.",
        "pmids": ["19151713"],
    },
    "rs13266634": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "SLC30A8 — beta-cell zinc transporter.",
        "pmids": ["17463249"],
    },

    # ── Cancer ───────────────────────────────────────────────────────────────
    "rs2981582": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "FGFR2 — strongest common breast cancer locus.",
        "pmids": ["17529967", "31104728"],
    },
    "rs17879961": {
        "evidence_level": "Clinical-Moderate",
        "evidence_summary": "CHEK2 I157T — moderate-penetrance breast/colon/prostate/kidney cancer.",
        "pmids": ["18781189"],
    },
    "rs1447295": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "8q24 — replicated prostate cancer locus.",
        "pmids": ["17401365"],
    },
    "rs10993994": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "MSMB — prostate cancer susceptibility.",
        "pmids": ["18264097"],
    },

    # ── Atrial Fib ────────────────────────────────────────────────────────────
    "rs2200733": {
        "evidence_level": "GWAS-A",
        "evidence_summary": "PITX2 / 4q25 — strongest common AF locus.",
        "pmids": ["17603472"],
    },

    # ── Celiac / Autoimmune ──────────────────────────────────────────────────
    "rs2187668": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "HLA-DQ2.5 tag — celiac disease primary risk haplotype.",
        "pmids": ["18568571"],
        "guidelines": ["ESPGHAN celiac guidelines"],
    },
    "rs7454108": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "HLA-DQ8 tag — second celiac-risk haplotype.",
        "pmids": ["18568571"],
        "guidelines": ["ESPGHAN"],
    },
    "rs4349859": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "HLA-B*27 tag — ankylosing spondylitis, uveitis.",
        "pmids": ["21743469"],
    },

    # ── Misc ─────────────────────────────────────────────────────────────────
    "rs671": {
        "evidence_level": "Clinical-Validated",
        "evidence_summary": "ALDH2 E487K — acetaldehyde accumulation; cancer risk with alcohol.",
        "pmids": ["19874574"],
    },
    "rs1815739": {
        "evidence_level": "GWAS / Functional",
        "evidence_summary": "ACTN3 R577X — fast-twitch fiber actinin function.",
        "pmids": ["12704393"],
    },
}


def get_reference(rsid: str) -> dict:
    """Return the reference entry for an rsID, or a default 'not catalogued' record."""
    return REFERENCES.get(rsid, {
        "evidence_level": "Catalogued (GWAS/lit)",
        "evidence_summary": "Variant included from peer-reviewed GWAS or clinical literature.",
        "pmids": [],
        "guidelines": [],
    })


def collect_references_used(tier1_results: list[dict]) -> list[dict]:
    """Return the catalogued references for variants that appeared in tier-1
    results, sorted by category and gene."""
    refs = []
    for r in tier1_results:
        ref = REFERENCES.get(r["rsid"])
        if ref:
            refs.append({
                "rsid": r["rsid"],
                "gene": r["gene"],
                "variant_name": r["variant_name"],
                "category": r["category"],
                **ref,
            })
    # Sort by category then gene
    refs.sort(key=lambda x: (x["category"], x["gene"]))
    return refs


# Level → CSS class (so report can color-code chips)
LEVEL_CLASSES = {
    "CPIC Level A": "lvl-cpic-a",
    "CPIC Level A (verify by HLA typing)": "lvl-cpic-a",
    "CPIC Level B": "lvl-cpic-b",
    "CPIC Level C": "lvl-cpic-c",
    "Clinical-Validated": "lvl-clin",
    "Clinical-Validated (rare variant)": "lvl-clin",
    "Clinical-Moderate": "lvl-clin-mod",
    "GWAS-A (replicated meta-analysis)": "lvl-gwas-a",
    "GWAS-A (strongest T2D variant)": "lvl-gwas-a",
    "GWAS-A": "lvl-gwas-a",
    "Population-Validated": "lvl-pop",
    "GWAS / Functional": "lvl-gwas-a",
    "Catalogued (GWAS/lit)": "lvl-other",
}


def level_class(level: str) -> str:
    return LEVEL_CLASSES.get(level, "lvl-other")
