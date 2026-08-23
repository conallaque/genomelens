"""
Health Economics — Clinical & Payer ROI for Genomic Interventions
=================================================================

Translates the genomic findings already produced by the pipeline (PGx
phenotypes, polygenic-risk tiers, APOE genotype) into a health-economic
view: what does acting on each finding cost, what adverse outcome does it
avert, and what is the return on that spend — for an individual, for a
clinic, and for a payer (insurer) covering a large member population.

This module is **decision-support economics, not accounting**. Every
dollar figure is an order-of-magnitude estimate drawn from published
pharmacoeconomic and screening literature. Actual ROI varies enormously by
clinic, payer, patient adherence, drug mix, and local costs. The numbers
here are meant to size opportunity, not to bill against.

Input shape
-----------
``analyze_health_economics(findings, snps_df)`` consumes the same dict the
pipeline already persists to ``tier1_results.json`` — specifically the
``pgx_summary``, ``prs_summary`` and ``apoe_genotype`` keys — so it can be
run either inside the pipeline or, for testing, straight against that JSON.
``snps_df`` is accepted for SNP-level confirmation (APOE / ACTN3 via the
unified ``snp_registry``); it is optional and the module degrades
gracefully when a SNP is absent.

Output shape
------------
``{status, findings_with_economics, clinic_dashboard, payer_impact,
disclaimer, ...}`` — see ``analyze_health_economics`` for the contract.

A finding is suppressed (not an error) when its economics are unknown, when
the genotype is normal/indeterminate, or when the polygenic tier is not
elevated. Missing data never raises; it just yields fewer findings.
"""

from __future__ import annotations

import pandas as pd

try:
    import snp_registry  # optional SNP-level confirmation (APOE / ACTN3)
except Exception:  # pragma: no cover - registry should always import
    snp_registry = None  # type: ignore


# ─── Economic model constants (documented assumptions) ─────────────────────
# All are estimates; surfaced in the disclaimer and editable in one place.

DISCOUNT_RATE = 0.03          # standard health-economic real discount rate
NPV_HORIZON_YEARS = 3         # NPV evaluated over a 3-year window

# Clinic revenue model — a genomics-enabled clinic offering this as a
# subscription service. These are business-model assumptions, not clinical.
CLINIC_REVENUE_MONTHLY = 500.0   # $/patient/month subscription
CLINIC_GROSS_MARGIN = 0.85       # gross margin on that subscription

DEFAULT_CLINIC_PATIENTS = 100
DEFAULT_PAYER_MEMBERS = 100_000


# ─── Pharmacogenomic intervention economics (one-time test cost) ───────────
# Keyed by the CPIC gene as it appears in pgx_summary. Each entry:
#   cost            — $ for the PGx test + consult (one-time, upfront)
#   outcome_value   — $ health-economic value of the averted adverse event
#   clinical_benefit— plain-language description of what is prevented
#   prevalence      — fraction of a general population with an *actionable*
#                     (atypical) phenotype for this gene (payer scaling)
#   qaly_gain       — QALYs gained per averted event (cost-per-QALY)
# Sources: CPIC guidelines + published PGx cost-effectiveness analyses;
# rounded to order-of-magnitude. ROI is derived (outcome_value / cost).
PGX_ECONOMICS: dict[str, dict] = {
    "CYP2C9": {
        "drug": "warfarin",
        "clinical_benefit": "Prevent major bleeding event from warfarin overdosing",
        "cost": 300, "p_rx": 0.08, "p_adr": 0.25, "rrr": 0.60, "adr_cost": 15_000,
        "prevalence": 0.35, "qaly_gain": 0.30,
        "src": "Johnson et al. (2017) Clin Pharmacol Ther — CPIC warfarin guideline",
    },
    "VKORC1": {
        "drug": "warfarin",
        "clinical_benefit": "Improve warfarin dosing accuracy, fewer INR excursions",
        "cost": 150, "p_rx": 0.08, "p_adr": 0.20, "rrr": 0.50, "adr_cost": 8_000,
        "prevalence": 0.37, "qaly_gain": 0.15,
        "src": "Johnson et al. (2017) Clin Pharmacol Ther — CPIC warfarin guideline",
    },
    "CYP2C19": {
        "drug": "clopidogrel",
        "clinical_benefit": "Avoid clopidogrel non-response / stent thrombosis (MACE)",
        "cost": 250, "p_rx": 0.10, "p_adr": 0.20, "rrr": 0.50, "adr_cost": 30_000,
        "prevalence": 0.30, "qaly_gain": 0.35,
        "src": "Kazi et al. (2014) Ann Intern Med — genotype-guided antiplatelet",
    },
    "CYP2D6": {
        "drug": "codeine / antidepressants",
        "clinical_benefit": "Avoid opioid toxicity or ineffective analgesia/SSRI dosing",
        "cost": 250, "p_rx": 0.15, "p_adr": 0.15, "rrr": 0.60, "adr_cost": 12_000,
        "prevalence": 0.30, "qaly_gain": 0.20,
        "src": "Hicks et al. (2015) Clin Pharmacol Ther — CPIC codeine; Bousman (2023) CPIC SSRI",
    },
    "CYP3A5": {
        "drug": "tacrolimus",
        "clinical_benefit": "Hit tacrolimus target faster, avoid rejection/toxicity",
        "cost": 250, "p_rx": 0.05, "p_adr": 0.30, "rrr": 0.55, "adr_cost": 10_000,
        "prevalence": 0.40, "qaly_gain": 0.20,
        "src": "Birdwell et al. (2015) Clin Pharmacol Ther — CPIC tacrolimus",
    },
    "TPMT": {
        "drug": "thiopurines (azathioprine / 6-MP)",
        "clinical_benefit": "Prevent serious myelosuppression from thiopurines",
        "cost": 200, "p_rx": 0.03, "p_adr": 0.35, "rrr": 0.70, "adr_cost": 25_000,
        "prevalence": 0.10, "qaly_gain": 0.40,
        "src": "Relling et al. (2019) Clin Pharmacol Ther — CPIC thiopurines",
    },
    "NUDT15": {
        "drug": "thiopurines (azathioprine / 6-MP)",
        "clinical_benefit": "Prevent thiopurine-induced leukopenia (esp. East Asian ancestry)",
        "cost": 200, "p_rx": 0.03, "p_adr": 0.35, "rrr": 0.70, "adr_cost": 25_000,
        "prevalence": 0.10, "qaly_gain": 0.40,
        "src": "Relling et al. (2019) Clin Pharmacol Ther — CPIC thiopurines",
    },
    "SLCO1B1": {
        "drug": "simvastatin",
        "clinical_benefit": "Prevent statin-induced myopathy / rhabdomyolysis",
        "cost": 150, "p_rx": 0.20, "p_adr": 0.10, "rrr": 0.50, "adr_cost": 6_000,
        "prevalence": 0.25, "qaly_gain": 0.15,
        "src": "Ramsey et al. (2014) Clin Pharmacol Ther — CPIC simvastatin/SLCO1B1",
    },
    "UGT1A1": {
        "drug": "irinotecan",
        "clinical_benefit": "Prevent severe irinotecan neutropenia / diarrhea",
        "cost": 200, "p_rx": 0.02, "p_adr": 0.40, "rrr": 0.60, "adr_cost": 14_000,
        "prevalence": 0.15, "qaly_gain": 0.30,
        "src": "Innocenti et al. (2009) J Clin Oncol — UGT1A1 irinotecan dosing",
    },
    "HLA-B*57:01": {
        "drug": "abacavir",
        "clinical_benefit": "Prevent abacavir hypersensitivity reaction",
        "cost": 150, "p_rx": 0.02, "p_adr": 0.55, "rrr": 0.95, "adr_cost": 20_000,
        "prevalence": 0.06, "qaly_gain": 0.50,
        "src": "Schackman et al. (2008) Ann Intern Med — cost-saving",
    },
}

# Phenotype keywords that make a PGx result *actionable* (i.e. it would
# change prescribing). Normal-function / indeterminate results are not
# counted as realized economic findings.
_ACTIONABLE_PGX_KEYWORDS = (
    "poor", "intermediate", "rapid", "ultrarapid", "decreased",
    "increased", "non-expressor", "high warfarin sensitivity",
    "positive", "deficien",
)


# ─── Polygenic / lifestyle intervention economics (recurring cost) ─────────
# Keyed by the PRS panel name in prs_summary (and APOE handled separately).
# Cost here is an *annual* program cost; outcome_value is the averted
# event's value. Triggered when the polygenic tier is Elevated or High.
PRS_ECONOMICS: dict[str, dict] = {
    "Coronary Artery Disease": {
        "finding": "Elevated coronary-artery-disease polygenic risk",
        "clinical_benefit": "Intensive lipid management (statin) to prevent MI",
        "intervention": "Intensive statin therapy + lipid monitoring",
        "cost": 500, "outcome_value": 250_000,
        "prevalence": 0.20, "qaly_gain": 1.50, "recurring": True,
        "src": "Mega et al. (2015) Lancet — PRS-guided statin benefit",
    },
    "Type 2 Diabetes": {
        "finding": "Elevated type-2-diabetes polygenic risk",
        "clinical_benefit": "Prevent / delay T2D onset",
        "intervention": "CGM + lifestyle coaching",
        "cost": 3_600, "outcome_value": 50_000,
        "prevalence": 0.20, "qaly_gain": 0.80, "recurring": True,
        "src": "Knowler (2002) NEJM — DPP; Khera (2016) NEJM — PRS × lifestyle",
    },
    "BMI / Obesity Tendency": {
        "finding": "Elevated obesity-tendency polygenic risk",
        "clinical_benefit": "Prevent obesity-driven T2D via structured program",
        "intervention": "Structured weight-management program",
        "cost": 1_200, "outcome_value": 30_000,
        "prevalence": 0.20, "qaly_gain": 0.40, "recurring": True,
        "src": "Locke et al. (2015) Nature — BMI GWAS; NICE CG189 weight management",
    },
    "Breast Cancer": {
        "finding": "Elevated breast-cancer polygenic risk",
        "clinical_benefit": "Risk-stratified screening (annual MRI + mammography)",
        "intervention": "Enhanced breast screening program",
        "cost": 1_200, "outcome_value": 120_000,
        "prevalence": 0.12, "qaly_gain": 1.80, "recurring": True,
        "src": "Pashayan et al. (2018) Genet Med — PRS-stratified screening CEA",
    },
    "Prostate Cancer": {
        "finding": "Elevated prostate-cancer polygenic risk",
        "clinical_benefit": "Risk-stratified PSA screening (earlier, more frequent)",
        "intervention": "Annual PSA + MRI surveillance",
        "cost": 800, "outcome_value": 80_000,
        "prevalence": 0.12, "qaly_gain": 1.20, "recurring": True,
        "src": "Callender et al. (2019) Ann Intern Med — PRS-informed PSA screening",
    },
    "Alzheimer's Disease": {
        "finding": "Elevated late-onset Alzheimer's polygenic risk",
        "clinical_benefit": "Intensive vascular risk management + cognitive reserve building",
        "intervention": "Multi-domain dementia prevention (FINGER-style)",
        "cost": 2_400, "outcome_value": 200_000,
        "prevalence": 0.15, "qaly_gain": 2.00, "recurring": True,
        "src": "Ngandu et al. (2015) Lancet — FINGER trial; Livingston (2020) Lancet — dementia prevention",
    },
    "Atrial Fibrillation": {
        "finding": "Elevated atrial-fibrillation polygenic risk",
        "clinical_benefit": "Proactive rhythm monitoring + stroke prevention",
        "intervention": "Annual ECG + pulse oximetry + anticoagulation if detected",
        "cost": 600, "outcome_value": 150_000,
        "prevalence": 0.08, "qaly_gain": 1.50, "recurring": True,
        "src": "Lubitz et al. (2017) Circ — PRS for AF; Hart (2007) Ann Intern Med — anticoagulation NNT",
    },
}

# ─── Expanded PGS Catalog economics (conditions beyond the core PRS panels) ──
# These map to conditions scored by pgs_catalog.py that are NOT already in
# PRS_ECONOMICS. Keyed by the label string from pgs_catalog.CONDITIONS.
EXPANDED_PGS_ECONOMICS: dict[str, dict] = {
    "Hypertension / Systolic BP": {
        "finding": "Elevated hypertension polygenic risk",
        "clinical_benefit": "Early antihypertensive therapy + lifestyle modification",
        "intervention": "Home BP monitoring + medication titration",
        "cost": 400, "outcome_value": 60_000,
        "prevalence": 0.30, "qaly_gain": 0.80, "recurring": True,
        "src": "Sun et al. (2021) JAMA Cardiol — PRS × BP treatment benefit",
    },
    "Ischemic Stroke": {
        "finding": "Elevated ischemic-stroke polygenic risk",
        "clinical_benefit": "Aggressive vascular risk-factor management",
        "intervention": "Antiplatelet + statin intensification + BP control",
        "cost": 600, "outcome_value": 180_000,
        "prevalence": 0.06, "qaly_gain": 2.00, "recurring": True,
        "src": "Abraham et al. (2019) Stroke — PRS for ischemic stroke",
    },
    "Chronic Kidney Disease": {
        "finding": "Elevated chronic-kidney-disease polygenic risk",
        "clinical_benefit": "Early nephroprotection (SGLT2i / ACEi) + monitoring",
        "intervention": "Annual eGFR + uACR screening + SGLT2 inhibitor if indicated",
        "cost": 1_800, "outcome_value": 100_000,
        "prevalence": 0.10, "qaly_gain": 1.20, "recurring": True,
        "src": "Wuttke et al. (2019) Nat Commun — CKD GWAS; Heerspink (2020) NEJM — SGLT2i",
    },
    "Major Depressive Disorder": {
        "finding": "Elevated major-depression polygenic risk",
        "clinical_benefit": "Proactive mental health screening + early CBT referral",
        "intervention": "Annual PHQ-9 screening + accessible therapy",
        "cost": 1_500, "outcome_value": 35_000,
        "prevalence": 0.15, "qaly_gain": 0.60, "recurring": True,
        "src": "Howard et al. (2019) Nat Neurosci — MDD GWAS; Chisholm (2016) Lancet Psych — CBT CEA",
    },
    "Inflammatory Bowel Disease": {
        "finding": "Elevated inflammatory-bowel-disease polygenic risk",
        "clinical_benefit": "Early gastroenterology referral + fecal calprotectin monitoring",
        "intervention": "Annual calprotectin + symptom surveillance",
        "cost": 500, "outcome_value": 45_000,
        "prevalence": 0.02, "qaly_gain": 0.80, "recurring": True,
        "src": "de Lange et al. (2017) Nat Genet — IBD GWAS; van der Valk (2016) IBD — IBD cost burden",
    },
    "Asthma": {
        "finding": "Elevated asthma polygenic risk",
        "clinical_benefit": "Environmental trigger reduction + early controller therapy",
        "intervention": "Allergen avoidance counseling + ICS if symptomatic",
        "cost": 600, "outcome_value": 25_000,
        "prevalence": 0.08, "qaly_gain": 0.40, "recurring": True,
        "src": "Demenais et al. (2018) Nat Genet — asthma GWAS; GINA 2023 guidelines",
    },
    "Rheumatoid Arthritis": {
        "finding": "Elevated rheumatoid-arthritis polygenic risk",
        "clinical_benefit": "Early rheumatology referral if joint symptoms emerge",
        "intervention": "Awareness + annual anti-CCP if symptomatic",
        "cost": 300, "outcome_value": 40_000,
        "prevalence": 0.02, "qaly_gain": 0.60, "recurring": True,
        "src": "Okada et al. (2014) Nature — RA GWAS; Finckh (2006) Arthritis Rheum — early DMARD CEA",
    },
}

# ─── HLA drug-hypersensitivity economics ────────────────────────────────────
# Keyed by allele name as it appears in hla_result["carrier_alleles"].
HLA_ECONOMICS: dict[str, dict] = {
    "HLA-B*58:01": {
        "finding": "HLA-B*58:01 carrier — allopurinol hypersensitivity risk",
        "clinical_benefit": "Avoid allopurinol SJS/TEN (use febuxostat instead)",
        "drug": "allopurinol",
        "cost": 150, "outcome_value": 45_000,
        "prevalence": 0.05, "qaly_gain": 0.40,
        "src": "Stamp et al. (2016) Intern Med J — HLA-B*58:01 CEA; CPIC Level A",
    },
    "HLA-A*31:01": {
        "finding": "HLA-A*31:01 carrier — carbamazepine hypersensitivity risk",
        "clinical_benefit": "Avoid carbamazepine DRESS/SJS (use alternative anticonvulsant)",
        "drug": "carbamazepine",
        "cost": 150, "outcome_value": 35_000,
        "prevalence": 0.04, "qaly_gain": 0.35,
        "src": "Plumpton et al. (2015) Epilepsia — HLA-A*31:01 screening CEA; CPIC Level A",
    },
    "HLA-B*15:02": {
        "finding": "HLA-B*15:02 carrier — carbamazepine/phenytoin SJS/TEN risk",
        "clinical_benefit": "Avoid aromatic anticonvulsants (SJS/TEN risk >1000×)",
        "drug": "carbamazepine/phenytoin",
        "cost": 150, "outcome_value": 60_000,
        "prevalence": 0.01, "qaly_gain": 0.50,
        "src": "Chen et al. (2011) Pharmacogenomics — HLA-B*15:02 screening CEA; FDA black-box",
    },
}

# ─── Carrier-screening economics (per condition affected/carrier) ────────────
CARRIER_ECONOMICS: dict[str, dict] = {
    "Hereditary Hemochromatosis (HH, type 1)": {
        "finding_affected": "HFE C282Y homozygous — hemochromatosis monitoring",
        "finding_carrier": "HFE C282Y carrier — reproductive awareness",
        "clinical_benefit_affected": "Early phlebotomy prevents organ damage",
        "clinical_benefit_carrier": "Partner testing for reproductive risk assessment",
        "cost_affected": 200, "outcome_affected": 60_000, "qaly_affected": 1.50,
        "cost_carrier": 100, "outcome_carrier": 2_000, "qaly_carrier": 0.10,
        "prev_affected": 0.005, "prev_carrier": 0.10,
        "src": "Adams et al. (2005) NEJM — hemochromatosis screening; phlebotomy is curative",
    },
    "Venous Thromboembolism (VTE) Susceptibility": {
        "finding_affected": "Factor V Leiden / Prothrombin homozygous — VTE prevention",
        "finding_carrier": "Factor V Leiden / Prothrombin carrier — situational anticoagulation",
        "clinical_benefit_affected": "Avoid estrogen-containing contraceptives; prophylactic anticoagulation perioperatively",
        "clinical_benefit_carrier": "Awareness for high-risk situations (surgery, pregnancy, long flights)",
        "cost_affected": 300, "outcome_affected": 35_000, "qaly_affected": 0.80,
        "cost_carrier": 100, "outcome_carrier": 10_000, "qaly_carrier": 0.30,
        "prev_affected": 0.002, "prev_carrier": 0.05,
        "src": "Cohn et al. (2013) J Thromb Haemost — FVL screening; Rosendaal (2005) Lancet",
    },
    "Hereditary Breast/Colon/Prostate/Kidney Cancer Susceptibility": {
        "finding_affected": "CHEK2 I157T homozygous — enhanced cancer surveillance",
        "finding_carrier": "CHEK2 I157T carrier — moderate cancer risk awareness",
        "clinical_benefit_affected": "Enhanced screening (mammography, colonoscopy)",
        "clinical_benefit_carrier": "Risk-aware screening schedule",
        "cost_affected": 800, "outcome_affected": 80_000, "qaly_affected": 1.20,
        "cost_carrier": 200, "outcome_carrier": 15_000, "qaly_carrier": 0.30,
        "prev_affected": 0.001, "prev_carrier": 0.05,
        "src": "Cybulski et al. (2011) J Clin Oncol — CHEK2 I157T cancer risks",
    },
    "Cystic Fibrosis": {
        "finding_affected": "CFTR homozygous/compound het — CF management",
        "finding_carrier": "CFTR carrier — reproductive genetic counseling",
        "clinical_benefit_affected": "Early modulator therapy (Trikafta)",
        "clinical_benefit_carrier": "Partner screening to assess 25% offspring risk",
        "cost_affected": 5_000, "outcome_affected": 300_000, "qaly_affected": 5.00,
        "cost_carrier": 250, "outcome_carrier": 5_000, "qaly_carrier": 0.15,
        "prev_affected": 0.0004, "prev_carrier": 0.04,
        "src": "Middleton et al. (2019) NEJM — Trikafta; ACOG — carrier screening CEA",
    },
    "Celiac Disease Susceptibility": {
        "finding_affected": "HLA-DQ2/DQ8 positive — celiac awareness",
        "finding_carrier": "HLA-DQ2/DQ8 carrier — celiac rule-out capability",
        "clinical_benefit_affected": "tTG-IgA testing if GI symptoms; early gluten-free diet",
        "clinical_benefit_carrier": "Negative DQ2/DQ8 rules out celiac",
        "cost_affected": 200, "outcome_affected": 15_000, "qaly_affected": 0.50,
        "cost_carrier": 50, "outcome_carrier": 3_000, "qaly_carrier": 0.10,
        "prev_affected": 0.01, "prev_carrier": 0.30,
        "src": "Green et al. (2015) Ann Intern Med — celiac screening; NICE NG20",
    },
    "Broad Autoimmunity Susceptibility (T1D, RA, lupus, Graves')": {
        "finding_affected": "PTPN22 R620W — broad autoimmunity awareness",
        "finding_carrier": "PTPN22 R620W carrier — autoimmune screening awareness",
        "clinical_benefit_affected": "Anti-inflammatory lifestyle; thyroid + autoimmune panel monitoring",
        "clinical_benefit_carrier": "Awareness for autoimmune symptom clusters",
        "cost_affected": 400, "outcome_affected": 25_000, "qaly_affected": 0.60,
        "cost_carrier": 100, "outcome_carrier": 5_000, "qaly_carrier": 0.10,
        "prev_affected": 0.005, "prev_carrier": 0.08,
        "src": "Bottini et al. (2004) Nat Genet — PTPN22; general autoimmune CEA literature",
    },
}

# ─── Compound-interaction economics ──────────────────────────────────────────
INTERACTION_ECONOMICS: dict[str, dict] = {
    "MTHFR C677T Homozygous (T/T)": {
        "finding": "MTHFR 677TT — methylfolate supplementation",
        "clinical_benefit": "Prevent hyperhomocysteinemia-driven cardiovascular / fertility risk",
        "cost": 300, "outcome_value": 20_000,
        "prevalence": 0.10, "qaly_gain": 0.40,
        "src": "Hickey et al. (2013) J Inherit Metab Dis — MTHFR; Klerk (2002) JAMA — homocysteine × CVD",
    },
    "MTHFR A1298C Homozygous (C/C)": {
        "finding": "MTHFR 1298CC — methylation support",
        "clinical_benefit": "BH4 pathway support; neurotransmitter synthesis optimization",
        "cost": 200, "outcome_value": 8_000,
        "prevalence": 0.08, "qaly_gain": 0.15,
        "src": "Weisberg et al. (2001) Mol Genet Metab — A1298C functional impact",
    },
    "HFE C282Y Homozygous — Classical Hemochromatosis Risk": {
        "finding": "HFE C282Y/C282Y — hemochromatosis iron monitoring",
        "clinical_benefit": "Prevent iron-overload organ damage via therapeutic phlebotomy",
        "cost": 200, "outcome_value": 60_000,
        "prevalence": 0.005, "qaly_gain": 1.50,
        "src": "Adams et al. (2005) NEJM; Allen et al. (2008) NEJM — C282Y penetrance",
    },
}

# APOE ε4 carriers — high CAD/AD risk; same statin-style intervention as the
# CAD PRS path but driven by genotype. Kept separate so it can fire even
# when no CAD PRS panel is present.
APOE_E4_ECONOMICS = {
    "finding": "APOE ε4 carrier (elevated cardiovascular / Alzheimer's risk)",
    "clinical_benefit": "Intensive cardiovascular risk reduction (statin) to prevent MI",
    "intervention": "Intensive statin therapy + lipid monitoring",
    "cost": 500, "outcome_value": 250_000,
    "prevalence": 0.25, "qaly_gain": 1.50, "recurring": True,
}

# Optional exercise / longevity inputs (only used if the caller supplies
# them in `findings`). Kept here so the economic assumptions live together.
VO2MAX_ECONOMICS = {
    "finding": "Low cardiorespiratory fitness (VO2max) capacity",
    "clinical_benefit": "Structured training to raise VO2max (all-cause mortality)",
    "intervention": "Supervised aerobic training program",
    "cost": 2_000, "outcome_value": 40_000,
    "prevalence": 0.30, "qaly_gain": 0.50, "recurring": True,
}
# LONGEVITY_VALUE_PER_PERCENTILE was removed in the 2026-08-22 provenance pass.
# It priced a percentile of a composite longevity score at a flat $10,000 with
# no source, and because the composite re-aggregates variants valued elsewhere
# in this module, every dollar it produced was also a double count. The
# longevity signal is still reported; see _exercise_longevity_findings.


# ─── Addiction genetics economics (keyed by substance category) ────────────
ADDICTION_ECONOMICS: dict[str, dict] = {
    "Alcohol": {
        "finding": "Alcohol-metabolism susceptibility variant",
        "clinical_benefit": "Personalized AUD risk awareness + naltrexone/acamprosate matching",
        "cost": 200, "outcome_value": 25_000,
        "prevalence": 0.10, "qaly_gain": 0.80,
        "src": "Kranzler & Soyka (2018) JAMA — pharmacogenomics of AUD treatment",
    },
    "Opioid": {
        "finding": "Opioid-response pharmacogenomic variant",
        "clinical_benefit": "OPRM1-guided opioid dosing + naltrexone response prediction",
        "cost": 200, "outcome_value": 40_000,
        "prevalence": 0.05, "qaly_gain": 1.20,
        "src": "Crist et al. (2018) Mol Psychiatry — OPRM1 mu-opioid receptor PGx",
    },
    "Nicotine": {
        "finding": "Nicotine-dependence susceptibility variant",
        "clinical_benefit": "CHRNA5-guided smoking cessation (varenicline dosing)",
        "cost": 500, "outcome_value": 35_000,
        "prevalence": 0.15, "qaly_gain": 1.50,
        "src": "Chen et al. (2012) Am J Psychiatry — CHRNA5 × cessation pharmacotherapy",
    },
    "Cannabis": {
        "finding": "Cannabis-sensitivity variant (CNR1/FAAH)",
        "clinical_benefit": "Endocannabinoid system awareness for pain management",
        "cost": 100, "outcome_value": 5_000,
        "prevalence": 0.10, "qaly_gain": 0.10,
        "src": "Hindocha et al. (2019) Transl Psychiatry — FAAH rs324420",
    },
    "Stress": {
        "finding": "Stress-axis variant (CRHR1/FKBP5)",
        "clinical_benefit": "Trauma-informed care + HPA-axis-targeted therapy selection",
        "cost": 200, "outcome_value": 12_000,
        "prevalence": 0.12, "qaly_gain": 0.40,
        "src": "Binder et al. (2008) JAMA Psych — FKBP5 × PTSD treatment response",
    },
}

# ─── Metal/oxidative stress economics (keyed by gene) ──────────────────────
METAL_OXIDATIVE_ECONOMICS: dict[str, dict] = {
    "G6PD": {
        "finding": "G6PD deficiency — oxidative drug avoidance",
        "clinical_benefit": "Avoid primaquine/dapsone/rasburicase hemolytic crisis",
        "cost": 150, "outcome_value": 20_000,
        "prevalence": 0.05, "qaly_gain": 0.40,
        "src": "WHO G6PD testing guideline (2018); CPIC rasburicase",
    },
    "LRRK2": {
        "finding": "LRRK2 G2019S carrier — Parkinson's early intervention",
        "clinical_benefit": "Prodromal PD monitoring + neuroprotective trial eligibility",
        "cost": 500, "outcome_value": 100_000,
        "prevalence": 0.01, "qaly_gain": 2.00,
        "src": "Healy et al. (2008) Lancet Neurol — LRRK2 penetrance + PD trials",
    },
    "ATP7B": {
        "finding": "ATP7B variant — Wilson's disease copper monitoring",
        "clinical_benefit": "Serum copper/ceruloplasmin monitoring + early chelation",
        "cost": 200, "outcome_value": 50_000,
        "prevalence": 0.003, "qaly_gain": 1.50,
        "src": "EASL CPG Wilson's disease (2012); Roberts & Schilsky (2023)",
    },
}

# ─── Mendelian randomization outcome costs (maps MR outcome → healthcare cost) ─
MR_OUTCOME_COSTS: dict[str, dict] = {
    "CAD": {"cost": 85_000, "qaly_loss": 1.50, "baseline_10yr": 0.10,
            "src": "Heidenreich (2011) Circulation"},
    "T2D": {"cost": 85_000, "qaly_loss": 2.00, "baseline_10yr": 0.15,
            "src": "ADA (2018) Diabetes Care"},
    "Stroke": {"cost": 100_000, "qaly_loss": 2.50, "baseline_10yr": 0.05,
               "src": "Xu et al. (2019) JAHA"},
    "AD": {"cost": 200_000, "qaly_loss": 3.00, "baseline_10yr": 0.12,
           "src": "Hurd et al. (2013) NEJM"},
    "CKD": {"cost": 100_000, "qaly_loss": 1.20, "baseline_10yr": 0.10,
            "src": "Honeycutt et al. (2013) JASN"},
    "Breast Cancer": {"cost": 120_000, "qaly_loss": 1.80, "baseline_10yr": 0.12,
                      "src": "Blumen et al. (2016) JNCCN"},
    "Prostate Cancer": {"cost": 80_000, "qaly_loss": 1.20, "baseline_10yr": 0.12,
                        "src": "Mariotto et al. (2011) JNCI"},
    "Atrial Fibrillation": {"cost": 50_000, "qaly_loss": 1.00, "baseline_10yr": 0.08,
                            "src": "Kim et al. (2011) Circ Cardiovasc Qual"},
}

# ─── Neurochemistry economics (keyed by gene) ─────────────────────────────
NEUROCHEMISTRY_ECONOMICS: dict[str, dict] = {
    "COMT": {
        "finding": "COMT Val158Met — psychiatric medication optimization",
        "clinical_benefit": "COMT-guided SSRI/SNRI selection (catecholamine clearance)",
        "cost": 200, "outcome_value": 15_000,
        "prevalence": 0.25, "qaly_gain": 0.30,
        "src": "Bousman (2023) CPIC antidepressant; Lachman et al. (1996)",
    },
    "BDNF": {
        "finding": "BDNF Val66Met — antidepressant response prediction",
        "clinical_benefit": "Exercise-first + BDNF-aware therapy selection",
        "cost": 200, "outcome_value": 10_000,
        "prevalence": 0.20, "qaly_gain": 0.25,
        "src": "Niitsu et al. (2013) J Affect Disord — BDNF × antidepressant",
    },
    "DRD2": {
        "finding": "DRD2/ANKK1 Taq1A — reward-pathway variant",
        "clinical_benefit": "Dopamine-aware addiction treatment + antipsychotic response",
        "cost": 200, "outcome_value": 18_000,
        "prevalence": 0.15, "qaly_gain": 0.40,
        "src": "Blum et al. (2014) J Genet Syndr Gene Ther — DRD2 reward deficiency",
    },
    "OPRM1": {
        "finding": "OPRM1 A118G — opioid/naltrexone response",
        "clinical_benefit": "OPRM1-guided analgesic dosing + AUD treatment matching",
        "cost": 200, "outcome_value": 25_000,
        "prevalence": 0.12, "qaly_gain": 0.60,
        "src": "Kranzler et al. (2013) Neuropsychopharmacology — OPRM1 × naltrexone",
    },
}

# ─── Urologic/genitourinary economics (keyed by category) ─────────────────
UROLOGIC_ECONOMICS: dict[str, dict] = {
    "BPH/Prostate": {
        "finding": "Prostate-risk variant — early screening",
        "clinical_benefit": "PRS-augmented PSA + MRI surveillance",
        "cost": 600, "outcome_value": 60_000,
        "prevalence": 0.12, "qaly_gain": 1.00,
        "src": "Callender et al. (2019) Ann Intern Med — prostate screening",
    },
    "Kidney Stones": {
        "finding": "Kidney-stone susceptibility variant",
        "clinical_benefit": "Dietary calcium/oxalate optimization + hydration protocol",
        "cost": 200, "outcome_value": 15_000,
        "prevalence": 0.10, "qaly_gain": 0.20,
        "src": "Curhan et al. (2004) JASN — nephrolithiasis prevention",
    },
    "OAB/Bladder": {
        "finding": "Overactive-bladder susceptibility variant",
        "clinical_benefit": "Anticholinergic PGx + behavioral therapy",
        "cost": 300, "outcome_value": 8_000,
        "prevalence": 0.15, "qaly_gain": 0.15,
        "src": "Irwin et al. (2009) BJU Int — OAB economic burden",
    },
    "Testicular": {
        "finding": "Testicular-cancer risk variant",
        "clinical_benefit": "Self-exam awareness + early ultrasound if symptomatic",
        "cost": 100, "outcome_value": 25_000,
        "prevalence": 0.01, "qaly_gain": 0.80,
        "src": "McGlynn & Cook (2009) JNCI — testicular cancer screening",
    },
}

# ─── ACMG actionable gene economics (ClinVar P/LP → intervention) ─────────
ACMG_GENE_ECONOMICS: dict[str, dict] = {
    "BRCA1": {"finding": "BRCA1 pathogenic variant", "cost": 2_000, "outcome_value": 150_000,
              "qaly_gain": 3.00, "clinical_benefit": "Enhanced screening + risk-reducing surgery option",
              "src": "Manchanda et al. (2015) J Clin Oncol"},
    "BRCA2": {"finding": "BRCA2 pathogenic variant", "cost": 2_000, "outcome_value": 140_000,
              "qaly_gain": 2.80, "clinical_benefit": "Enhanced screening + risk-reducing surgery",
              "src": "Manchanda et al. (2015) J Clin Oncol"},
    "MLH1": {"finding": "MLH1 Lynch syndrome variant", "cost": 1_500, "outcome_value": 100_000,
             "qaly_gain": 2.50, "clinical_benefit": "Annual colonoscopy from age 25 + aspirin",
             "src": "Ladabaum et al. (2011) Ann Intern Med"},
    "MSH2": {"finding": "MSH2 Lynch syndrome variant", "cost": 1_500, "outcome_value": 100_000,
             "qaly_gain": 2.50, "clinical_benefit": "Annual colonoscopy + gynecologic surveillance",
             "src": "Ladabaum et al. (2011) Ann Intern Med"},
    "MSH6": {"finding": "MSH6 Lynch syndrome variant", "cost": 1_500, "outcome_value": 80_000,
             "qaly_gain": 2.00, "clinical_benefit": "Enhanced colonoscopy + endometrial screening",
             "src": "Ladabaum et al. (2011) Ann Intern Med"},
    "PMS2": {"finding": "PMS2 Lynch syndrome variant", "cost": 1_500, "outcome_value": 60_000,
             "qaly_gain": 1.50, "clinical_benefit": "Colonoscopy surveillance program",
             "src": "Ladabaum et al. (2011) Ann Intern Med"},
    "LDLR": {"finding": "LDLR familial hypercholesterolemia", "cost": 500, "outcome_value": 200_000,
             "qaly_gain": 3.50, "clinical_benefit": "High-intensity statin + cascade screening",
             "src": "Nherera et al. (2011) Heart — FH CEA"},
    "APOB": {"finding": "APOB familial hypercholesterolemia", "cost": 500, "outcome_value": 180_000,
             "qaly_gain": 3.00, "clinical_benefit": "High-intensity statin therapy",
             "src": "Nherera et al. (2011) Heart — FH CEA"},
    "PCSK9": {"finding": "PCSK9 familial hypercholesterolemia", "cost": 3_000, "outcome_value": 200_000,
              "qaly_gain": 3.50, "clinical_benefit": "PCSK9 inhibitor + cascade screening",
              "src": "Kazi et al. (2017) JAMA Cardiol"},
    "SCN5A": {"finding": "SCN5A channelopathy", "cost": 1_000, "outcome_value": 120_000,
              "qaly_gain": 4.00, "clinical_benefit": "Cardiac monitoring + beta-blocker/ICD",
              "src": "Kaufman et al. (2014) Circ Cardiovasc Genet"},
    "KCNQ1": {"finding": "KCNQ1 long-QT syndrome", "cost": 800, "outcome_value": 110_000,
              "qaly_gain": 3.80, "clinical_benefit": "Beta-blocker + activity restriction",
              "src": "Kaufman et al. (2014) Circ Cardiovasc Genet"},
    "KCNH2": {"finding": "KCNH2 long-QT syndrome", "cost": 800, "outcome_value": 110_000,
              "qaly_gain": 3.80, "clinical_benefit": "Beta-blocker + QT-prolonging drug avoidance",
              "src": "Kaufman et al. (2014) Circ Cardiovasc Genet"},
    "RET": {"finding": "RET MEN2 variant", "cost": 5_000, "outcome_value": 250_000,
            "qaly_gain": 5.00, "clinical_benefit": "Prophylactic thyroidectomy + calcitonin",
            "src": "Wells et al. (2015) Thyroid; Brandi (2001) JCEM"},
    "TP53": {"finding": "TP53 Li-Fraumeni variant", "cost": 3_000, "outcome_value": 200_000,
             "qaly_gain": 4.00, "clinical_benefit": "Annual whole-body MRI (Toronto protocol)",
             "src": "Villani et al. (2016) Lancet Oncol"},
    "RB1": {"finding": "RB1 retinoblastoma variant", "cost": 500, "outcome_value": 80_000,
            "qaly_gain": 3.00, "clinical_benefit": "Pediatric eye exams + family screening",
            "src": "Soliman et al. (2016) J AAPOS"},
    "MYH7": {"finding": "MYH7 hypertrophic cardiomyopathy", "cost": 1_200, "outcome_value": 90_000,
             "qaly_gain": 2.50, "clinical_benefit": "Echo surveillance + exercise restriction",
             "src": "Maron et al. (2014) Circulation — HCM guidelines"},
    "MYBPC3": {"finding": "MYBPC3 hypertrophic cardiomyopathy", "cost": 1_200, "outcome_value": 90_000,
               "qaly_gain": 2.50, "clinical_benefit": "Echo surveillance + cascade screening",
               "src": "Maron et al. (2014) Circulation"},
}

# ─── PheWAS extreme-tier category economics ───────────────────────────────
PHEWAS_CATEGORY_ECONOMICS: dict[str, dict] = {
    "Lipids": {"finding": "Extreme-tier lipid prediction", "cost": 500,
               "outcome_value": 30_000, "qaly_gain": 0.50,
               "clinical_benefit": "Lipid-lowering therapy",
               "src": "CTT Collaboration (2010) Lancet"},
    "Glucose/Diabetes": {"finding": "Extreme-tier glucose prediction", "cost": 1_200,
                         "outcome_value": 50_000, "qaly_gain": 0.80,
                         "clinical_benefit": "Diabetes prevention program",
                         "src": "DPP Research Group (2002) NEJM"},
    "Cardiovascular": {"finding": "Extreme-tier cardiovascular prediction", "cost": 800,
                       "outcome_value": 40_000, "qaly_gain": 0.60,
                       "clinical_benefit": "Cardiovascular risk management",
                       "src": "Yusuf et al. (2004) Lancet — INTERHEART"},
    "Kidney": {"finding": "Extreme-tier kidney prediction", "cost": 600,
               "outcome_value": 35_000, "qaly_gain": 0.40,
               "clinical_benefit": "Nephroprotection monitoring",
               "src": "Wuttke et al. (2019) Nat Commun"},
    "Inflammation": {"finding": "Extreme-tier inflammatory prediction", "cost": 400,
                     "outcome_value": 15_000, "qaly_gain": 0.30,
                     "clinical_benefit": "Anti-inflammatory lifestyle + monitoring",
                     "src": "Ridker (2003) Circulation — CRP"},
    "Liver": {"finding": "Extreme-tier liver prediction", "cost": 300,
              "outcome_value": 25_000, "qaly_gain": 0.40,
              "clinical_benefit": "Hepatoprotection + NAFLD monitoring",
              "src": "GBD (2020) J Hepatol"},
    "Hematology": {"finding": "Extreme-tier hematology prediction", "cost": 200,
                   "outcome_value": 8_000, "qaly_gain": 0.15,
                   "clinical_benefit": "Iron/B12/folate monitoring",
                   "src": "Kassebaum et al. (2014) Blood"},
    "Body Composition": {"finding": "Extreme-tier body composition prediction", "cost": 800,
                         "outcome_value": 20_000, "qaly_gain": 0.30,
                         "clinical_benefit": "Weight management program",
                         "src": "Locke et al. (2015) Nature — BMI"},
}

# ─── Immunogenetics economics (keyed by gene) ─────────────────────────────
IMMUNOGENETICS_ECONOMICS: dict[str, dict] = {
    "CCR5": {
        "finding": "CCR5-delta32 — HIV resistance",
        "clinical_benefit": "PrEP decision support (natural HIV resistance factor)",
        "cost": 100, "outcome_value": 30_000,
        "prevalence": 0.10, "qaly_gain": 0.50,
        "src": "Samson et al. (1996) Nature — CCR5-Δ32",
    },
    "IL28B": {
        "finding": "IL28B rs12979860 — HCV treatment response",
        "clinical_benefit": "IL28B-guided HCV DAA treatment selection",
        "cost": 200, "outcome_value": 50_000,
        "prevalence": 0.05, "qaly_gain": 1.00,
        "src": "Ge et al. (2009) Nature — IL28B × HCV treatment",
    },
    "HbS": {
        "finding": "HbS sickle-cell trait — anesthesia/altitude awareness",
        "clinical_benefit": "Sickle-trait perioperative management protocol",
        "cost": 100, "outcome_value": 15_000,
        "prevalence": 0.08, "qaly_gain": 0.20,
        "src": "Tsaras et al. (2009) Am J Med — sickle trait complications",
    },
}

# ─── Wellness prediction economics (keyed by category) ────────────────────
WELLNESS_ECONOMICS: dict[str, dict] = {
    "Nutrition": {
        "finding": "Nutrient-metabolism variant — deficiency prevention",
        "clinical_benefit": "Targeted supplementation (vitamin D/B12/folate/iron)",
        "cost": 300, "outcome_value": 8_000,
        "prevalence": 0.25, "qaly_gain": 0.15,
        "src": "Holick (2007) NEJM — vitamin D; WHO nutrition guidelines",
    },
    "Fitness": {
        "finding": "Exercise-response variant — training optimization",
        "clinical_benefit": "Genotype-matched training program (injury prevention)",
        "cost": 600, "outcome_value": 12_000,
        "prevalence": 0.20, "qaly_gain": 0.20,
        "src": "Bouchard et al. (2011) Med Sci Sports Exerc — HERITAGE",
    },
    "Sleep": {
        "finding": "Sleep-architecture variant — chronotherapy",
        "clinical_benefit": "Chronotype-aligned schedule + sleep hygiene protocol",
        "cost": 200, "outcome_value": 6_000,
        "prevalence": 0.30, "qaly_gain": 0.10,
        "src": "Jones et al. (2019) Nat Commun — chronotype GWAS",
    },
}

# ─── Detoxification economics (keyed by gene) ─────────────────────────────
DETOX_ECONOMICS: dict[str, dict] = {
    "NAT2": {
        "finding": "NAT2 slow acetylator — occupational cancer risk awareness",
        "clinical_benefit": "Aromatic amine exposure monitoring + bladder screening",
        "cost": 400, "outcome_value": 35_000,
        "prevalence": 0.40, "qaly_gain": 0.60,
        "src": "Garcia-Closas et al. (2005) Lancet — NAT2 × bladder cancer",
    },
    "PON1": {
        "finding": "PON1 Q192R — organophosphate susceptibility",
        "clinical_benefit": "Pesticide exposure avoidance + cholinesterase monitoring",
        "cost": 200, "outcome_value": 15_000,
        "prevalence": 0.25, "qaly_gain": 0.30,
        "src": "Costa et al. (2013) Annu Rev Pharmacol Toxicol — PON1",
    },
    "GSTT1": {
        "finding": "GSTT1 null genotype — reduced glutathione conjugation",
        "clinical_benefit": "Cruciferous vegetable intake + antioxidant monitoring",
        "cost": 100, "outcome_value": 6_000,
        "prevalence": 0.20, "qaly_gain": 0.10,
        "src": "Economopoulos & Sergentanis (2010) Cancer Epidemiol — GST × cancer",
    },
}

# ─── Family planning economics (partner testing ROI per gene) ──────────────
PARTNER_TESTING_COST = 300
GENETIC_COUNSELING_COST = 400
_AFFECTED_CHILD_COST = {
    "Cystic Fibrosis": 300_000,
    "Sickle Cell Disease": 200_000,
    "Tay-Sachs Disease": 150_000,
    "Spinal Muscular Atrophy": 250_000,
    "Phenylketonuria": 80_000,
    "Hemochromatosis": 60_000,
    "Beta-Thalassemia": 180_000,
    "Congenital Adrenal Hyperplasia": 120_000,
}


DISCLAIMER = (
    "These are estimates based on published pharmacoeconomic literature, not "
    "guarantees. Actual ROI varies by clinic, payer, patient adherence, drug "
    "mix, and local costs. Intervention costs are one-time for pharmacogenomic "
    "tests and annual for lifestyle programs; outcome values are the modeled "
    "health-economic value of an averted adverse event. Figures size "
    "opportunity for decision support — they are not billing or actuarial values."
)


# ─── Core economic math ────────────────────────────────────────────────────

def calculate_roi(cost: float, outcome_value: float) -> float | None:
    """ROI as a simple ratio (outcome value / cost). None if cost is zero."""
    if not cost or cost <= 0:
        return None
    return round(outcome_value / cost, 1)


def calculate_payback_months(cost: float, outcome_value: float,
                             recurring_cost: bool = False,
                             rate: float = DISCOUNT_RATE) -> float | None:
    """Months until cumulative discounted benefit exceeds cumulative cost.

    Walks month-by-month: benefit accrues at ``outcome_value / 12`` per month
    (discounted), cost is either upfront or spread over months. Returns the
    first month where cumulative benefit ≥ cumulative cost, or None if it
    never breaks even within 10 years.
    """
    if not outcome_value or outcome_value <= 0:
        return None
    monthly_benefit = outcome_value / 12.0
    monthly_cost = cost / 12.0 if recurring_cost else 0.0
    cum_benefit = 0.0
    cum_cost = cost if not recurring_cost else 0.0
    monthly_rate = (1.0 + rate) ** (1.0 / 12.0) - 1.0
    for m in range(1, 121):
        disc = 1.0 / (1.0 + monthly_rate) ** m
        cum_benefit += monthly_benefit * disc
        cum_cost += monthly_cost * disc
        if cum_benefit >= cum_cost:
            return float(m)
    return None


def calculate_npv(
    cost: float,
    outcome_value: float,
    recurring_cost: bool = False,
    horizon: int = NPV_HORIZON_YEARS,
    rate: float = DISCOUNT_RATE,
) -> float:
    """Net present value over ``horizon`` years at ``rate``.

    Benefit (``outcome_value``) accrues each year and is discounted. Cost is
    either a one-time upfront spend at t=0 (pharmacogenomic tests) or a
    recurring annual spend discounted over the horizon (lifestyle programs).
    """
    benefit = sum(outcome_value / (1 + rate) ** t for t in range(1, horizon + 1))
    if recurring_cost:
        spend = sum(cost / (1 + rate) ** t for t in range(1, horizon + 1))
    else:
        spend = cost  # upfront, undiscounted at t=0
    return round(benefit - spend, 2)


def _econ_record(
    finding: str,
    clinical_benefit: str,
    cost: float,
    outcome_value: float,
    confidence: str,
    *,
    recurring: bool = False,
    source: str = "",
    prevalence: float = 0.0,
    qaly_gain: float = 0.0,
    evidence: str = "",
) -> dict:
    """Assemble one finding's economics record with derived metrics."""
    return {
        "finding": finding,
        "clinical_benefit": clinical_benefit,
        "intervention_cost": cost,
        "cost_basis": "annual" if recurring else "one-time",
        "outcome_value": outcome_value,
        "roi": calculate_roi(cost, outcome_value),
        "payback_months": calculate_payback_months(cost, outcome_value,
                                                    recurring_cost=recurring),
        "npv_3year": calculate_npv(cost, outcome_value, recurring_cost=recurring),
        "confidence": confidence,
        "category": source,
        "prevalence": prevalence,
        "qaly_gain": qaly_gain,
        "evidence": evidence,
    }


# ─── Finding extractors ────────────────────────────────────────────────────

def _is_actionable_pgx(phenotype: str) -> bool:
    p = (phenotype or "").lower()
    if not p or "indeterminate" in p or "normal" in p:
        # "Normal Warfarin Sensitivity" etc. — not actionable.
        return False
    return any(k in p for k in _ACTIONABLE_PGX_KEYWORDS)


def _pgx_findings(pgx_summary: dict) -> list[dict]:
    out: list[dict] = []
    for gene, info in (pgx_summary or {}).items():
        econ = PGX_ECONOMICS.get(gene)
        if econ is None:
            continue
        phenotype = (info or {}).get("phenotype", "")
        if not _is_actionable_pgx(phenotype):
            continue
        # Decomposed PGx valuation aligned with PGX_CEA in value_of_information:
        #   expected_benefit = p_rx × p_adr × rrr × adr_cost
        # where p_rx = probability patient is ever prescribed this drug,
        # p_adr = probability of ADR given Rx + actionable genotype,
        # rrr = relative risk reduction from genotype-guided prescribing.
        p_rx = econ.get("p_rx", 0.10)
        p_adr = econ.get("p_adr", 0.15)
        rrr = econ.get("rrr", 0.50)
        adr_cost = econ.get("adr_cost", 10_000)
        conditional_value = round(p_rx * p_adr * rrr * adr_cost)
        out.append(_econ_record(
            finding=f"{gene} {phenotype} ({econ['drug']})",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=conditional_value,
            confidence="high", source="Pharmacogenomics",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=(f"{gene} phenotype: {phenotype} — "
                      f"p_rx={p_rx} × p_adr={p_adr} × rrr={rrr} × ${adr_cost:,}"),
        ))
    return out


def _prs_findings(prs_summary: dict) -> list[dict]:
    out: list[dict] = []
    for panel, info in (prs_summary or {}).items():
        econ = PRS_ECONOMICS.get(panel)
        if econ is None:
            continue
        tier = (info or {}).get("tier")
        if tier not in ("Elevated", "High"):
            continue
        pct = (info or {}).get("percentile")
        ev = f"{panel} polygenic tier: {tier}"
        if pct is not None:
            ev += f" ({pct:g}th percentile)"
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=f"{econ['clinical_benefit']} — {econ['intervention']}",
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate", source="Polygenic Risk", recurring=econ["recurring"],
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=ev,
        ))
    return out


def _apoe_genotype_from_snps(snps_df: pd.DataFrame | None) -> str | None:
    """Best-effort APOE ε4 confirmation via the registry SNP, used only when
    the caller did not supply an apoe_genotype. Returns None if undetermined."""
    if snps_df is None or snp_registry is None:
        return None
    try:
        if "rs429358" not in snps_df.index:
            return None
        gt = str(snps_df.loc["rs429358"].get("genotype", "")).upper()
        gt = gt.replace(" ", "").replace("-", "")
        if "C" in gt:  # C allele at rs429358 marks an ε4 haplotype
            return "e4-carrier (rs429358 C)"
    except Exception:
        return None
    return None


def _apoe_findings(apoe_genotype: str | None, snps_df) -> list[dict]:
    geno = apoe_genotype
    confirm = ""
    if not geno:
        geno = _apoe_genotype_from_snps(snps_df)
        confirm = " (rs429358-derived)"
    if not geno or ("4" not in str(geno) and "e4" not in str(geno).lower()):
        return []
    e = APOE_E4_ECONOMICS
    return [_econ_record(
        finding=e["finding"],
        clinical_benefit=f"{e['clinical_benefit']} — {e['intervention']}",
        cost=e["cost"], outcome_value=e["outcome_value"],
        confidence="moderate", source="Genotype", recurring=e["recurring"],
        prevalence=e["prevalence"], qaly_gain=e["qaly_gain"],
        evidence=f"APOE genotype: {geno}{confirm}",
    )]


def _exercise_longevity_findings(findings: dict) -> list[dict]:
    """Optional findings the caller may supply: low VO2max tier, and a
    longevity percentile to quantify improvement headroom."""
    out: list[dict] = []
    vo2 = str(findings.get("vo2max_tier", "")).lower()
    if vo2 in ("low", "poor", "below average"):
        e = VO2MAX_ECONOMICS
        out.append(_econ_record(
            finding=e["finding"],
            clinical_benefit=f"{e['clinical_benefit']} — {e['intervention']}",
            cost=e["cost"], outcome_value=e["outcome_value"],
            confidence="moderate", source="Exercise / Lifestyle",
            recurring=e["recurring"], prevalence=e["prevalence"],
            qaly_gain=e["qaly_gain"], evidence=f"VO2max tier: {vo2}",
        ))
    pct = findings.get("longevity_percentile")
    if isinstance(pct, (int, float)) and 0 <= pct < 50:
        headroom = 50 - pct  # percentiles of achievable improvement to median
        # NOT MONETISED, deliberately. This used to be valued at a flat
        # $10,000 per percentile of headroom, which made it the single largest
        # line in the report — on one measured genome it produced $431,000 of
        # claimed benefit and 54% of the entire modelled total. Two independent
        # reasons to stop:
        #
        #   1. The longevity percentile is a COMPOSITE of variants that are
        #      already valued individually elsewhere in this module. Valuing it
        #      too counts the same genotypes a second time.
        #   2. The $10,000/percentile rate had no source. A percentile of a
        #      composite score is not a disease, has no cost of illness, and no
        #      published dose-response converts one into dollars.
        #
        # It stays in the report as a signal worth acting on — it just no
        # longer carries a number the model can be wrong about.
        out.append(_econ_record(
            finding="Below-median longevity composite (improvement headroom)",
            clinical_benefit=(
                f"~{headroom:g} percentile points of headroom to the population "
                f"median. Reported as a signal, not monetised — the composite "
                f"aggregates variants already valued individually above."),
            cost=0, outcome_value=0,
            confidence="low", source="Longevity", recurring=False,
            prevalence=0.0, qaly_gain=0.0,
            evidence=f"Longevity percentile: {pct:g} (not monetised)",
        ))
    return out


def _expanded_pgs_findings(expanded_pgs_result: dict | None) -> list[dict]:
    """Findings from expanded PGS catalog conditions not already in PRS_ECONOMICS."""
    out: list[dict] = []
    if not expanded_pgs_result or not expanded_pgs_result.get("panels"):
        return out
    already_keyed = set(PRS_ECONOMICS.keys())
    for panel_name, info in expanded_pgs_result["panels"].items():
        if panel_name in already_keyed:
            continue
        econ = EXPANDED_PGS_ECONOMICS.get(panel_name)
        if econ is None:
            continue
        tier = (info.get("result") or info).get("tier")
        if tier not in ("Elevated", "High"):
            continue
        pct = (info.get("result") or info).get("percentile")
        ev = f"{panel_name} expanded PGS tier: {tier}"
        if pct is not None:
            ev += f" ({pct:g}th percentile)"
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=f"{econ['clinical_benefit']} — {econ['intervention']}",
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="low", source="Expanded Polygenic Score",
            recurring=econ["recurring"], prevalence=econ["prevalence"],
            qaly_gain=econ["qaly_gain"], evidence=ev,
        ))
    return out


def _hla_findings(hla_result: dict | None) -> list[dict]:
    """Drug-hypersensitivity economics from HLA typing (beyond HLA-B*57:01
    which is already in PGX_ECONOMICS via the PGx pathway)."""
    out: list[dict] = []
    if not hla_result:
        return out
    carrier_alleles = hla_result.get("carrier_alleles", [])
    for allele in carrier_alleles:
        econ = HLA_ECONOMICS.get(allele)
        if econ is None:
            continue
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="high", source="HLA Pharmacogenomics",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{allele} carrier — avoid {econ['drug']}",
        ))
    return out


def _carrier_findings(carrier_result: dict | None) -> list[dict]:
    """Economic findings from carrier screening — both affected (homozygous)
    and carrier (heterozygous) states with distinct valuations."""
    out: list[dict] = []
    if not carrier_result:
        return out
    for record in (carrier_result.get("affected") or []):
        disease = record.get("disease", "")
        econ = CARRIER_ECONOMICS.get(disease)
        if econ is None:
            continue
        out.append(_econ_record(
            finding=econ["finding_affected"],
            clinical_benefit=econ["clinical_benefit_affected"],
            cost=econ["cost_affected"], outcome_value=econ["outcome_affected"],
            confidence="high", source="Carrier Screening",
            prevalence=econ.get("prev_affected", 0.005),
            qaly_gain=econ["qaly_affected"],
            evidence=f"{record.get('gene', '?')} {record.get('variant', '')} — homozygous affected",
        ))
    for record in (carrier_result.get("carriers") or []):
        disease = record.get("disease", "")
        econ = CARRIER_ECONOMICS.get(disease)
        if econ is None:
            continue
        out.append(_econ_record(
            finding=econ["finding_carrier"],
            clinical_benefit=econ["clinical_benefit_carrier"],
            cost=econ["cost_carrier"], outcome_value=econ["outcome_carrier"],
            confidence="moderate", source="Carrier Screening",
            prevalence=econ.get("prev_carrier", 0.04),
            qaly_gain=econ["qaly_carrier"],
            evidence=f"{record.get('gene', '?')} {record.get('variant', '')} — heterozygous carrier",
        ))
    return out


def _interaction_findings(interactions_result: dict | None) -> list[dict]:
    """Economic findings from compound multi-variant interactions."""
    out: list[dict] = []
    if not interactions_result:
        return out
    for finding in (interactions_result.get("findings") or []):
        title = finding.get("title", "")
        econ = INTERACTION_ECONOMICS.get(title)
        if econ is None:
            continue
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate" if finding.get("severity") == "high" else "low",
            source="Compound Interaction",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{title} (severity: {finding.get('severity', 'unknown')})",
        ))
    return out


def _addiction_findings(addiction_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not addiction_result or not addiction_result.get("available"):
        return out
    seen_cats: set = set()
    for finding in (addiction_result.get("findings") or []):
        if finding.get("impact") not in ("susceptible", "clinically_useful"):
            continue
        cat = finding.get("category", "")
        if not cat or cat in seen_cats:
            continue
        econ = ADDICTION_ECONOMICS.get(cat)
        if econ is None:
            continue
        seen_cats.add(cat)
        gene = finding.get("gene", "")
        out.append(_econ_record(
            finding=f"{gene} — {econ['finding']}",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate" if finding.get("confidence") == "high" else "low",
            source="Addiction Genetics",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{gene} {finding.get('genotype', '')} — {finding.get('verdict', '')}",
        ))
    return out


def _metal_oxidative_findings(metal_oxidative_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not metal_oxidative_result:
        return out
    seen_genes: set = set()
    for pred in (metal_oxidative_result.get("predictions") or []):
        gene = ""
        cv = pred.get("clinical_variant")
        if cv:
            gene = cv.get("gene", "")
        if not gene:
            trait = pred.get("trait", "")
            for g in METAL_OXIDATIVE_ECONOMICS:
                if g.lower() in trait.lower():
                    gene = g
                    break
        if not gene or gene in seen_genes:
            continue
        econ = METAL_OXIDATIVE_ECONOMICS.get(gene)
        if econ is None:
            continue
        seen_genes.add(gene)
        out.append(_econ_record(
            finding=econ["finding"],
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate" if pred.get("confidence") == "high" else "low",
            source="Metal/Oxidative",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{pred.get('trait', '')} — {pred.get('result', '')}",
        ))
    return out


def _mr_findings(mr_result: dict | None) -> list[dict]:
    """Economic findings from Mendelian randomization: causal risk → cost."""
    out: list[dict] = []
    if not mr_result:
        return out
    seen: set = set()
    for f in (mr_result.get("findings") or []):
        if f.get("status") != "computed":
            continue
        outcome = f.get("outcome", "")
        if outcome in seen:
            continue
        rr = f.get("outcome_relative_risk")
        if rr is None or rr <= 1.0:
            continue
        oc = None
        for key, val in MR_OUTCOME_COSTS.items():
            if key.lower() in outcome.lower() or outcome.lower() in key.lower():
                oc = val
                break
        if oc is None:
            continue
        seen.add(outcome)
        excess_risk = (rr - 1.0) * oc["baseline_10yr"]
        avoidable_cost = round(excess_risk * oc["cost"])
        qaly = round(excess_risk * oc["qaly_loss"], 2)
        if avoidable_cost < 100:
            continue
        exposure = f.get("exposure", "")
        out.append(_econ_record(
            finding=f"MR: elevated {exposure} → {outcome} (RR {rr:.2f})",
            clinical_benefit=f"Reduce {exposure} exposure to lower {outcome} risk",
            cost=500, outcome_value=avoidable_cost,
            confidence="moderate", source="Mendelian Randomization",
            prevalence=oc["baseline_10yr"], qaly_gain=qaly,
            evidence=(f"{exposure} → {outcome}: RR {rr:.2f}, "
                      f"excess 10yr risk {excess_risk:.3f} ({f.get('n_used', 0)} instruments)"),
        ))
    return out


def _neurochemistry_findings(neurochemistry_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not neurochemistry_result or not neurochemistry_result.get("available"):
        return out
    seen: set = set()
    for finding in (neurochemistry_result.get("findings") or []):
        if finding.get("impact") in ("neutral",):
            continue
        gene = finding.get("gene", "")
        if not gene or gene in seen:
            continue
        econ = NEUROCHEMISTRY_ECONOMICS.get(gene)
        if econ is None:
            continue
        seen.add(gene)
        out.append(_econ_record(
            finding=f"{gene} {finding.get('name', '')} — {econ['finding']}",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="low", source="Neurochemistry",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{gene} {finding.get('genotype', '')} — {finding.get('verdict', '')}",
        ))
    return out


def _urologic_findings(urologic_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not urologic_result or not urologic_result.get("available"):
        return out
    seen_cats: set = set()
    for finding in (urologic_result.get("findings") or []):
        if finding.get("impact") in ("normal", "typical"):
            continue
        cat = finding.get("category", "")
        if not cat or cat in seen_cats:
            continue
        econ = UROLOGIC_ECONOMICS.get(cat)
        if econ is None:
            continue
        seen_cats.add(cat)
        out.append(_econ_record(
            finding=f"{finding.get('gene', '')} — {econ['finding']}",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="low" if finding.get("confidence") != "high" else "moderate",
            source="Urologic/GU",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{finding.get('trait', '')} — {finding.get('result', '')}",
        ))
    return out


def _clinical_variant_findings(clinical_variants_result: dict | None) -> list[dict]:
    """Economics for ClinVar P/LP variants in ACMG actionable genes."""
    out: list[dict] = []
    if not clinical_variants_result or not clinical_variants_result.get("available"):
        return out
    seen_genes: set = set()
    for v in (clinical_variants_result.get("findings") or []):
        sig = (v.get("significance") or "").lower()
        if "pathogenic" not in sig:
            continue
        gene = v.get("gene", "")
        if not gene or gene in seen_genes:
            continue
        econ = ACMG_GENE_ECONOMICS.get(gene)
        if econ is None:
            continue
        seen_genes.add(gene)
        zyg = v.get("zygosity", "heterozygous")
        out.append(_econ_record(
            finding=f"{econ['finding']} ({zyg})",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="high", source="Clinical Variant (ClinVar)",
            prevalence=0.005, qaly_gain=econ["qaly_gain"],
            evidence=f"{gene} {v.get('ref', '')}{v.get('pos', '')} → ClinVar {sig}",
        ))
    return out


def _phewas_findings(phewas_result: dict | None) -> list[dict]:
    """Economics for extreme-tier PheWAS biomarker predictions."""
    out: list[dict] = []
    if not phewas_result:
        return out
    seen_cats: set = set()
    for h in (phewas_result.get("headline") or []):
        cat = h.get("category", "")
        if not cat or cat in seen_cats:
            continue
        econ = PHEWAS_CATEGORY_ECONOMICS.get(cat)
        if econ is None:
            continue
        tier = h.get("tier", "")
        if tier not in ("Very High", "Very Low", "High", "Low"):
            continue
        seen_cats.add(cat)
        out.append(_econ_record(
            finding=f"{h.get('trait', cat)} — {econ['finding']}",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="low", source="PheWAS Biomarker",
            prevalence=0.10, qaly_gain=econ["qaly_gain"],
            evidence=(f"{h.get('trait', '')} predicted {h.get('predicted_value', '')} "
                      f"{h.get('unit', '')} (p{h.get('percentile', '')} — {tier})"),
        ))
    return out


def _immunogenetics_findings(immunogenetics_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not immunogenetics_result or not immunogenetics_result.get("available"):
        return out
    seen: set = set()
    for finding in (immunogenetics_result.get("findings") or []):
        if finding.get("impact") in ("intermediate", "neutral"):
            continue
        gene = finding.get("gene", "")
        if not gene or gene in seen:
            continue
        econ = IMMUNOGENETICS_ECONOMICS.get(gene)
        if econ is None:
            continue
        seen.add(gene)
        out.append(_econ_record(
            finding=f"{gene} — {econ['finding']}",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=econ["outcome_value"],
            confidence="moderate" if finding.get("confidence") == "high" else "low",
            source="Immunogenetics",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=f"{gene} {finding.get('genotype', '')} — {finding.get('verdict', '')}",
        ))
    return out


def _wellness_findings(wellness_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not wellness_result:
        return out
    seen_cats: set = set()
    for pred in (wellness_result.get("predictions") or []):
        cat = pred.get("category", "")
        if not cat or cat in seen_cats:
            continue
        econ = WELLNESS_ECONOMICS.get(cat)
        if econ is None:
            continue
        result_text = (pred.get("result") or "").lower()
        if any(k in result_text for k in ("reduced", "impaired", "low", "poor",
                                           "slow", "deficien", "elevated risk")):
            seen_cats.add(cat)
            out.append(_econ_record(
                finding=f"{pred.get('trait', cat)} — {econ['finding']}",
                clinical_benefit=econ["clinical_benefit"],
                cost=econ["cost"], outcome_value=econ["outcome_value"],
                confidence="low", source="Wellness Prediction",
                prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
                evidence=f"{pred.get('trait', '')} — {pred.get('result', '')}",
            ))
    return out


def _detox_findings(detox_result: dict | None) -> list[dict]:
    out: list[dict] = []
    if not detox_result or not detox_result.get("available"):
        return out
    seen: set = set()
    for finding in (detox_result.get("findings") or []):
        gene = finding.get("gene", "")
        if not gene or gene in seen:
            continue
        econ = DETOX_ECONOMICS.get(gene)
        if econ is None:
            continue
        result_text = (finding.get("result") or "").lower()
        if any(k in result_text for k in ("reduced", "slow", "null", "deficien",
                                           "impaired", "absent")):
            seen.add(gene)
            out.append(_econ_record(
                finding=f"{gene} — {econ['finding']}",
                clinical_benefit=econ["clinical_benefit"],
                cost=econ["cost"], outcome_value=econ["outcome_value"],
                confidence="low" if finding.get("confidence") != "high" else "moderate",
                source="Detoxification",
                prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
                evidence=f"{finding.get('trait', '')} — {finding.get('result', '')}",
            ))
    return out


def _family_planning_findings(family_planning_result: dict | None) -> list[dict]:
    """Partner testing ROI from reproductive genetics — cost of test vs
    expected cost of affected offspring."""
    out: list[dict] = []
    if not family_planning_result or not family_planning_result.get("available"):
        return out
    for item in (family_planning_result.get("recessive_items") or []):
        disease = item.get("disease", "")
        p_carrier = item.get("partner_carrier_freq")
        p_child = item.get("child_clinical_risk")
        if not p_carrier or not p_child:
            continue
        try:
            p_carrier = float(str(p_carrier).rstrip("%")) / (
                100.0 if "%" in str(item.get("partner_carrier_freq", "")) else 1.0)
        except (ValueError, TypeError):
            continue
        try:
            p_child = float(str(p_child).rstrip("%")) / (
                100.0 if "%" in str(item.get("child_clinical_risk", "")) else 1.0)
        except (ValueError, TypeError):
            continue
        affected_cost = _AFFECTED_CHILD_COST.get(disease, 150_000)
        expected_avoidable = round(p_carrier * p_child * affected_cost)
        if expected_avoidable < 50:
            continue
        test_cost = PARTNER_TESTING_COST + GENETIC_COUNSELING_COST
        out.append(_econ_record(
            finding=f"{disease} carrier — partner testing value",
            clinical_benefit=f"Partner carrier testing for {disease} reproductive planning",
            cost=test_cost, outcome_value=expected_avoidable,
            confidence="moderate", source="Family Planning",
            prevalence=p_carrier, qaly_gain=round(p_child * 2.0, 2),
            evidence=(f"{item.get('gene', '')} carrier × partner freq {p_carrier:.3f} "
                      f"× child risk {p_child:.3f} × ${affected_cost:,}"),
        ))
    return out


def _top_drugs_findings(top_drugs_result: dict | None) -> list[dict]:
    """High-value PGx findings from the top-prescribed-drugs screen:
    for drugs that are genotype-actionable, the p_rx is effectively 1.0
    because these are commonly prescribed medications the patient is
    likely to encounter."""
    out: list[dict] = []
    if not top_drugs_result or not top_drugs_result.get("available"):
        return out
    seen_genes: set = set()
    for drug_info in (top_drugs_result.get("actionable") or []):
        gene = drug_info.get("gene", "")
        if not gene or gene in seen_genes:
            continue
        econ = PGX_ECONOMICS.get(gene)
        if econ is None:
            continue
        seen_genes.add(gene)
        p_adr = econ.get("p_adr", 0.15)
        rrr = econ.get("rrr", 0.50)
        adr_cost = econ.get("adr_cost", 10_000)
        boosted_value = round(0.50 * p_adr * rrr * adr_cost)
        if boosted_value <= econ.get("outcome_value", 0):
            continue
        out.append(_econ_record(
            finding=f"{gene} — top-drug actionable ({drug_info.get('drug', '')})",
            clinical_benefit=econ["clinical_benefit"],
            cost=econ["cost"], outcome_value=boosted_value,
            confidence="high", source="Top-Drugs PGx Screen",
            prevalence=econ["prevalence"], qaly_gain=econ["qaly_gain"],
            evidence=(f"High-prevalence drug {drug_info.get('drug', '')} × "
                      f"{gene} actionable phenotype (boosted p_rx ≈ 0.50)"),
        ))
    return out


# ─── Scaling ───────────────────────────────────────────────────────────────

def scale_to_clinic(findings_econ: dict, patient_count: int = DEFAULT_CLINIC_PATIENTS) -> dict:
    """Scale per-finding economics to a clinic population, weighting each
    finding by the fraction of patients who actually have it (prevalence)."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return {"patient_count": patient_count, "n_findings": 0,
                "note": "No actionable findings with economics for this profile."}

    # Prevalence-weighted cost and benefit per patient: each finding
    # contributes its cost × prevalence (fraction of patients who have it).
    weighted_cost = sum(f["intervention_cost"] * f.get("prevalence", 0.10)
                        for f in findings)
    weighted_benefit = sum(f["outcome_value"] * f.get("prevalence", 0.10)
                          for f in findings)
    avg_cost = round(weighted_cost, 2)
    avg_benefit = round(weighted_benefit, 2)

    monthly_margin = CLINIC_REVENUE_MONTHLY * CLINIC_GROSS_MARGIN
    payback = round(avg_cost / monthly_margin, 1) if monthly_margin else None

    return {
        "patient_count": patient_count,
        "n_findings": len(findings),
        "avg_cost_per_patient": avg_cost,
        "avg_benefit_per_patient": avg_benefit,
        "avg_roi": calculate_roi(avg_cost, avg_benefit),
        "total_cost": round(avg_cost * patient_count, 2),
        "total_benefit": round(avg_benefit * patient_count, 2),
        "revenue_model_monthly": CLINIC_REVENUE_MONTHLY,
        "gross_margin": CLINIC_GROSS_MARGIN,
        "payback_period_months": payback,
        "summary": (
            f"Applied to {patient_count} patients: "
            f"prevalence-weighted cost ${round(avg_cost * patient_count):,}, "
            f"modeled benefit ${round(avg_benefit * patient_count):,}, "
            f"ROI {calculate_roi(avg_cost, avg_benefit)}:1"
        ),
    }


def scale_to_payer(findings_econ: dict, member_population: int = DEFAULT_PAYER_MEMBERS) -> dict:
    """Scale per-finding economics to a payer's member population using each
    finding's population prevalence."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return {"member_population": member_population, "affected_members": 0,
                "note": "No actionable findings with economics for this profile."}

    affected_total = 0
    total_cost = 0.0
    total_benefit = 0.0
    total_qalys = 0.0
    per_finding: list[dict] = []
    for f in findings:
        affected = round(member_population * f.get("prevalence", 0.0))
        cost = affected * f["intervention_cost"]
        benefit = affected * f["outcome_value"]
        qalys = affected * f.get("qaly_gain", 0.0)
        affected_total += affected
        total_cost += cost
        total_benefit += benefit
        total_qalys += qalys
        per_finding.append({
            "finding": f["finding"], "affected_members": affected,
            "total_cost": round(cost), "total_benefit": round(benefit),
        })

    # Cap affected_total at plan size — independent prevalence sums
    # double-count members who qualify for multiple interventions.
    unique_affected = min(affected_total, member_population)

    cost_per_qaly = round(total_cost / total_qalys) if total_qalys else None
    return {
        "member_population": member_population,
        "affected_members": unique_affected,
        "affected_member_interventions": affected_total,
        "total_cost": round(total_cost),
        "total_benefit": round(total_benefit),
        "roi": calculate_roi(total_cost, total_benefit),
        "cost_per_qaly": cost_per_qaly,
        "net_savings": round(total_benefit - total_cost),
        "per_finding": per_finding,
        "summary": (
            f"Applied to {member_population:,} members: "
            f"{unique_affected:,} unique members affected, "
            f"{affected_total:,} intervention-events, cost ${round(total_cost):,}, "
            f"modeled savings ${round(total_benefit):,}, "
            f"ROI {calculate_roi(total_cost, total_benefit)}:1"
        ),
    }


# ─── Markdown summary ──────────────────────────────────────────────────────

def generate_economics_summary(findings_econ: dict) -> str:
    """Human-readable markdown for the HTML report."""
    findings = findings_econ.get("findings_with_economics", [])
    if not findings:
        return ("## Health Economics\n\nNo actionable genomic findings with "
                "modeled economics for this profile.\n")

    # Main summary shows only high-confidence findings; if none are
    # high-confidence we fall back to the best available with a note so the
    # summary is never silently empty.
    high_conf = findings_econ.get("high_confidence") or [
        f for f in findings if f.get("confidence") == "high"
    ]
    main = high_conf if high_conf else findings
    top = sorted(main, key=lambda f: (f.get("roi") or 0), reverse=True)[:3]
    lines = ["## Health Economics — Clinical & Payer ROI", ""]
    heading = ("### Top high-confidence interventions by ROI" if high_conf
               else "### Top interventions by ROI (no high-confidence findings)")
    lines.append(heading)
    for f in top:
        lines.append(
            f"- **{f['finding']}** — {f['clinical_benefit']}. "
            f"Cost ${f['intervention_cost']:,} ({f['cost_basis']}), "
            f"value ${f['outcome_value']:,}, **ROI {f['roi']}:1**, "
            f"payback {f['payback_months']} mo, NPV(3y) ${f['npv_3year']:,.0f}."
        )

    clinic = findings_econ.get("clinic_dashboard", {})
    if clinic.get("n_findings"):
        lines += ["", "### Clinic dashboard", clinic.get("summary", ""),
                  f"- Avg cost/patient: ${clinic['avg_cost_per_patient']:,} · "
                  f"avg benefit/patient: ${clinic['avg_benefit_per_patient']:,} · "
                  f"avg ROI {clinic['avg_roi']}:1",
                  f"- Subscription revenue ${clinic['revenue_model_monthly']:,.0f}/mo "
                  f"at {int(clinic['gross_margin']*100)}% margin · "
                  f"payback {clinic['payback_period_months']} mo"]

    payer = findings_econ.get("payer_impact", {})
    if payer.get("affected_members"):
        lines += ["", "### Payer impact", payer.get("summary", ""),
                  f"- Cost per QALY: "
                  f"${payer['cost_per_qaly']:,}" if payer.get("cost_per_qaly")
                  else "- Cost per QALY: n/a",
                  f"- Net modeled savings: ${payer['net_savings']:,}"]

    lines += ["", f"_{DISCLAIMER}_"]
    return "\n".join(lines)


# ─── Master analyzer ───────────────────────────────────────────────────────

def analyze_health_economics(findings: dict, snps_df: pd.DataFrame,
                             expanded_pgs_result: dict | None = None,
                             hla_result: dict | None = None,
                             carrier_result: dict | None = None,
                             interactions_result: dict | None = None,
                             addiction_result: dict | None = None,
                             metal_oxidative_result: dict | None = None,
                             mr_result: dict | None = None,
                             neurochemistry_result: dict | None = None,
                             urologic_result: dict | None = None,
                             clinical_variants_result: dict | None = None,
                             phewas_result: dict | None = None,
                             immunogenetics_result: dict | None = None,
                             wellness_result: dict | None = None,
                             detox_result: dict | None = None,
                             family_planning_result: dict | None = None,
                             top_drugs_result: dict | None = None) -> dict:
    """Compute clinical & payer ROI for genomic interventions.

    Accepts results from all pipeline analysis modules and produces
    economic findings for every data source that has a valuation pathway.
    """
    findings = findings or {}
    econ_findings: list[dict] = []
    try:
        econ_findings += _pgx_findings(findings.get("pgx_summary", {}))
        econ_findings += _prs_findings(findings.get("prs_summary", {}))
        econ_findings += _apoe_findings(findings.get("apoe_genotype"), snps_df)
        econ_findings += _exercise_longevity_findings(findings)
        econ_findings += _expanded_pgs_findings(expanded_pgs_result)
        econ_findings += _hla_findings(hla_result)
        econ_findings += _carrier_findings(carrier_result)
        econ_findings += _interaction_findings(interactions_result)
        econ_findings += _addiction_findings(addiction_result)
        econ_findings += _metal_oxidative_findings(metal_oxidative_result)
        econ_findings += _mr_findings(mr_result)
        econ_findings += _neurochemistry_findings(neurochemistry_result)
        econ_findings += _urologic_findings(urologic_result)
        econ_findings += _clinical_variant_findings(clinical_variants_result)
        econ_findings += _phewas_findings(phewas_result)
        econ_findings += _immunogenetics_findings(immunogenetics_result)
        econ_findings += _wellness_findings(wellness_result)
        econ_findings += _detox_findings(detox_result)
        econ_findings += _family_planning_findings(family_planning_result)
        econ_findings += _top_drugs_findings(top_drugs_result)
    except Exception as e:  # never raise from the pipeline
        return {"status": "error", "error": str(e),
                "findings_with_economics": [], "disclaimer": DISCLAIMER}

    # Deduplicate: when the same condition appears from both carrier and
    # interaction extractors (e.g. HFE C282Y), keep the higher-value one.
    seen_conditions: dict[str, int] = {}
    deduped: list[dict] = []
    for f in econ_findings:
        key = f["finding"].lower().split("—")[0].strip()
        if key in seen_conditions:
            existing_idx = seen_conditions[key]
            if f["outcome_value"] > deduped[existing_idx]["outcome_value"]:
                deduped[existing_idx] = f
            continue
        seen_conditions[key] = len(deduped)
        deduped.append(f)
    econ_findings = deduped

    # Rank by ROI (highest first); keep None ROIs last.
    econ_findings.sort(key=lambda f: (f.get("roi") or 0), reverse=True)

    result: dict = {
        "status": "success" if econ_findings else "no_findings",
        "n_findings": len(econ_findings),
        "findings_with_economics": econ_findings,
        "high_confidence": [f for f in econ_findings if f["confidence"] == "high"],
        "disclaimer": DISCLAIMER,
    }
    result["clinic_dashboard"] = scale_to_clinic(result)
    result["payer_impact"] = scale_to_payer(result)
    return result


# ══════════════════════════════════════════════════════════════════════════
# PERSONAL ECONOMIC-IMPACT ANALYSIS (standalone economic_analysis.html)
# ══════════════════════════════════════════════════════════════════════════
#
# Models the individual's 10-year modeled economic impact of ACTING on their
# results — expected medical-cost avoidance + monetised quality-of-life (QALY)
# gains, net of intervention cost — across genomic (PGx / carrier / PRS) and
# blood-work (PREVENT cardiovascular risk, prediabetes, biological age)
# findings. All figures are population-average model estimates, not a promise
# for any individual; the disclaimer is surfaced prominently in the report.

VALUE_PER_QALY = 100_000          # standard US cost-effectiveness threshold
PERSONAL_HORIZON_YEARS = 10
# Expected event costs below are one-time amounts realised at some point inside
# PERSONAL_HORIZON_YEARS, not annual flows, so they are discounted to present
# value at the horizon midpoint — the standard simplification when event timing
# within a window is unmodelled. Without this the report labels its output
# "Over 10 years" while summing undiscounted future dollars as if they were
# present dollars, which overstates the benefit.
_MIDPOINT_DISCOUNT = 1.0 / (1.0 + DISCOUNT_RATE) ** (PERSONAL_HORIZON_YEARS / 2.0)
_ANALYSIS_COST = 700              # ~$200 genome + ~$500 lab panel, one-time

# Cost-of-illness / intervention assumptions (US, order-of-magnitude, 10-yr).
_MACE_COST = 85_000               # acute + 1-yr MI/stroke care
_MACE_QALY = 1.5
_STATIN_RRR = 0.27                # relative risk reduction, primary prevention
_STATIN_COST_10YR = 500
_T2D_COST = 85_000                # lifetime excess medical cost of T2D
_T2D_QALY = 2.0
_DPP_RRR = 0.58                   # lifestyle-program incidence reduction (DPP)
_DPP_COST_10YR = 3_000
_PREDIAB_PROGRESSION_10YR = 0.35  # prediabetes → diabetes over ~10 yr

# Shared with value_of_information: averting a case frees the MARGINAL cost of that
# case, not its AVERAGE lifetime cost (fixed capacity persists whether or not one
# case is prevented). Single source of truth, with a fallback so this module stays
# importable standalone.
try:
    from .value_of_information import MARGINAL_COST_FRACTION as _MARGINAL_COST_FRACTION
except Exception:                 # pragma: no cover — standalone import
    _MARGINAL_COST_FRACTION = 0.60


def _money(x) -> str:
    x = round(x)
    return f"-${abs(x):,}" if x < 0 else f"${x:,}"



# Personal-economics categories mapped to the same three adherence archetypes
# the pooled engine uses. Keyed on what acting on the finding asks of the
# person — a daily tablet, an appointment you attend, or a habit you keep —
# because that, not the disease, is what predicts whether they keep doing it.
# Anything unmapped falls back to ``adherence_default`` rather than to 1.0, so
# a gap in this table shows up as a discount rather than as a flattering
# silence.
_CATEGORY_ADHERENCE: dict[str, str] = {
    "Pharmacogenomic / genomic": "adherence_pharmacological",
    "HLA Pharmacogenomics":      "adherence_pharmacological",
    "Neurochemistry":            "adherence_pharmacological",
    "Cardiovascular":            "adherence_pharmacological",
    "Clinical Variant":          "adherence_screening",
    "Carrier Screening":         "adherence_screening",
    "PheWAS Biomarker":          "adherence_screening",
    "Expanded Polygenic Score":  "adherence_screening",
    "Mendelian Randomization":   "adherence_screening",
    "Compound Interaction":      "adherence_screening",
    "Metabolic":                 "adherence_lifestyle",
    "Addiction Genetics":        "adherence_lifestyle",
    "Wellness Genetics":         "adherence_lifestyle",
    "Biological aging":          "adherence_lifestyle",
}


def _adherence_for_category(category: str) -> float:
    """Real-world adherence multiplier for a personal-economics category."""
    try:
        from . import params as _ep
        return float(_ep.value(
            _CATEGORY_ADHERENCE.get(category, "adherence_default")))
    except Exception:
        # The page must still produce a number if the registry is missing. It
        # then produces the older, undiscounted one — which is the reason the
        # item carries its adherence factor, so the report can say which.
        return 1.0


def analyze_personal_economics(economics_result: dict | None = None,
                               bloodwork_result: dict | None = None,
                               genetic_age_result: dict | None = None,
                               meta: dict | None = None,
                               carrier_result: dict | None = None,
                               hla_result: dict | None = None,
                               interactions_result: dict | None = None,
                               expanded_pgs_result: dict | None = None,
                               addiction_result: dict | None = None,
                               neurochemistry_result: dict | None = None,
                               mr_result: dict | None = None,
                               clinical_variants_result: dict | None = None,
                               family_planning_result: dict | None = None,
                               phewas_result: dict | None = None,
                               wellness_result: dict | None = None) -> dict:
    """Build the personal 10-year economic-impact model from the run's results."""
    items: list[dict] = []
    # Findings that carry a real decision but deliberately no dollar figure.
    # Kept as a first-class output rather than dropped, so "not monetised" is
    # visible in the report as a choice instead of looking like an omission.
    not_monetised: list[dict] = []

    def add(category, finding, avoided, qaly, intervention, confidence, basis):
        # REAL-WORLD ADHERENCE. Everything below is trial efficacy: the benefit
        # if the person does the thing. The pooled payer analysis in
        # value_of_information already charges the efficacy-to-effectiveness
        # gap; this page was still reporting the undiscounted figure, so one
        # report described the same genome two ways. Same three archetypes,
        # keyed on what acting actually asks of the person.
        _adh = _adherence_for_category(category)
        # MARGINAL vs AVERAGE: the per-condition cost constants here (_MACE_COST,
        # _T2D_COST, the lifetime COI figures) are AVERAGE costs. Averting one case
        # frees only the marginal cost — a large share of average cost is fixed
        # capacity that persists. Scale the cash side down accordingly; this only
        # ever reduces the claimed saving. See value_of_information.MARGINAL_COST_FRACTION.
        avoided = avoided * _MARGINAL_COST_FRACTION
        # Discount to present value over the stated horizon (see
        # _MIDPOINT_DISCOUNT) so "over N years" describes the arithmetic.
        avoided = avoided * _MIDPOINT_DISCOUNT
        qv = qaly * VALUE_PER_QALY * _MIDPOINT_DISCOUNT
        # Benefit and ongoing cost are scaled by the same factor: someone who
        # stops taking the statin stops paying for it. Scaling only the benefit
        # would be as wrong in the other direction.
        avoided, qv = avoided * _adh, qv * _adh
        qaly, intervention = qaly * _adh, intervention * _adh
        items.append({
            "category": category, "finding": finding,
            "avoided": round(avoided), "qaly": round(qaly, 2),
            "qaly_value": round(qv), "intervention": round(intervention),
            "net": round(avoided + qv - intervention),
            "confidence": confidence, "basis": basis,
            "adherence": round(_adh, 3),
        })

    # ── Genomic actionable findings (reuse the curated per-condition econ) ──
    # outcome_value from _pgx_findings is already conditional (p_rx × p_adr ×
    # rrr × adr_cost), so we do NOT apply an additional 0.30 discount — that
    # was a bare constant standing in for the decomposition that is now explicit.
    if economics_result:
        for f in economics_result.get("findings_with_economics", []):
            outcome = f.get("outcome_value") or f.get("benefit") or 0
            prev = f.get("prevalence", 0.15)
            qaly = (f.get("qaly_gain") or f.get("qaly") or 0) * prev
            cost = f.get("cost", 200)
            avoided = outcome * prev
            label = (f.get("clinical_benefit") or f.get("finding")
                     or f.get("drug") or "Genomic finding")
            if avoided <= 0 and qaly <= 0:
                continue
            add("Pharmacogenomic / genomic", label, avoided, qaly, cost,
                f.get("confidence", "moderate"),
                "Avoided adverse event × probability of relevant exposure (curated per-condition model).")

    # ── Blood-work derived ──
    adv = ((bloodwork_result or {}).get("clinical") or {}).get("advanced") or {}
    prevent = next((i["value"] for i in adv.get("indices", [])
                    if i.get("id") == "prevent_ascvd"), None)
    if prevent is not None:
        p = prevent / 100.0
        avoided = p * _MACE_COST * _STATIN_RRR
        qaly = p * _MACE_QALY * _STATIN_RRR
        add("Cardiovascular", f"10-yr ASCVD risk {prevent:.1f}% — risk-factor management",
            avoided, qaly, _STATIN_COST_10YR, "high" if prevent >= 7.5 else "moderate",
            f"PREVENT 10-yr ASCVD probability × ${_MACE_COST:,} event cost × {_STATIN_RRR:.0%} "
            "relative risk reduction from statin/lifestyle.")

    # Prediabetes → T2D (from biological-age input glucose or an HbA1c flag)
    bio = adv.get("biological_age") or {}
    glucose = (bio.get("inputs") or {}).get("glucose")
    a1c_flag = any("HbA1c" in (fl.get("name", "")) or "Glucose" in (fl.get("name", ""))
                   for fl in (((bloodwork_result or {}).get("clinical") or {}).get("flags") or []))
    prediabetic = (glucose is not None and 100 <= glucose < 126) or a1c_flag
    if prediabetic:
        avoided = _PREDIAB_PROGRESSION_10YR * _T2D_COST * _DPP_RRR
        qaly = _PREDIAB_PROGRESSION_10YR * _T2D_QALY * _DPP_RRR
        add("Metabolic", "Prediabetes flagged — intensive lifestyle (DPP)",
            avoided, qaly, _DPP_COST_10YR, "high",
            f"~{_PREDIAB_PROGRESSION_10YR:.0%} 10-yr progression risk × ${_T2D_COST:,} lifetime "
            f"T2D cost × {_DPP_RRR:.0%} reduction from a diabetes-prevention program.")

    # ── HLA drug-hypersensitivity (personal) ──
    if hla_result:
        for allele in (hla_result.get("carrier_alleles") or []):
            econ = HLA_ECONOMICS.get(allele)
            if econ is None:
                continue
            avoided = econ["outcome_value"] * 0.15
            qaly = econ["qaly_gain"]
            add("HLA Pharmacogenomics", econ["finding"],
                avoided, qaly, econ["cost"], "high",
                f"{allele} carrier — published cost-effectiveness for pre-prescription testing "
                f"({econ.get('src', 'CPIC')}).")

    # ── Carrier screening (personal) ──
    if carrier_result:
        for record in (carrier_result.get("affected") or []):
            disease = record.get("disease", "")
            econ = CARRIER_ECONOMICS.get(disease)
            if econ is None:
                continue
            add("Carrier Screening", econ["finding_affected"],
                econ["outcome_affected"] * 0.30, econ["qaly_affected"],
                econ["cost_affected"], "high",
                f"{record.get('gene', '?')} homozygous — {econ['clinical_benefit_affected']} "
                f"({econ.get('src', '')}).")
        for record in (carrier_result.get("carriers") or []):
            disease = record.get("disease", "")
            econ = CARRIER_ECONOMICS.get(disease)
            if econ is None:
                continue
            add("Carrier Screening", econ["finding_carrier"],
                econ["outcome_carrier"] * 0.05, econ["qaly_carrier"],
                econ["cost_carrier"], "moderate",
                f"{record.get('gene', '?')} heterozygous carrier — "
                f"{econ['clinical_benefit_carrier']}.")

    # ── Compound interactions (personal) ──
    if interactions_result:
        for finding in (interactions_result.get("findings") or []):
            title = finding.get("title", "")
            econ = INTERACTION_ECONOMICS.get(title)
            if econ is None:
                continue
            avoided = econ["outcome_value"] * 0.25
            add("Compound Interaction", econ["finding"],
                avoided, econ["qaly_gain"], econ["cost"],
                "moderate" if finding.get("severity") == "high" else "low",
                f"{title} — {econ['clinical_benefit']} ({econ.get('src', '')}).")

    # ── Expanded PGS (personal) ──
    if expanded_pgs_result:
        already_seen = set()
        for panel_name, info in (expanded_pgs_result.get("panels") or {}).items():
            if panel_name in PRS_ECONOMICS:
                continue
            econ = EXPANDED_PGS_ECONOMICS.get(panel_name)
            if econ is None or panel_name in already_seen:
                continue
            tier = (info.get("result") or info).get("tier")
            if tier not in ("Elevated", "High"):
                continue
            already_seen.add(panel_name)
            avoided = econ["outcome_value"] * 0.10
            add("Expanded Polygenic Score", econ["finding"],
                avoided, econ["qaly_gain"] * 0.5, econ["cost"], "low",
                f"{panel_name} PGS tier {tier} — {econ['clinical_benefit']} "
                f"({econ.get('src', '')}).")

    # Biological aging — grounded in the PhenoAge clock's own validated 10-year
    # mortality output (mortality-calibrated; each year of acceleration ≈ 9%
    # higher all-cause mortality, HR 1.09/yr, Levine 2018). We compare the
    # person's modeled 10-yr mortality risk against the baseline risk for their
    # chronological age (PhenoAge = age), and value the excess/deficit in QALYs.
    accel = bio.get("accel")
    mort_pct = bio.get("mortality_10yr_pct")
    chrono = bio.get("chronological")
    if accel is not None and mort_pct is not None and chrono is not None:
        import math
        # 10-yr mortality for PhenoAge == chronological age (inverse of the clock)
        L = (chrono - 141.50225) * 0.090165
        base_m = 1.0 - math.exp(-math.exp(L) / 0.00553)
        cur_m = mort_pct / 100.0
        excess = cur_m - base_m                      # >0 if biologically older
        remaining_qalys = 15.0                        # discounted remaining healthspan
        qaly = abs(excess) * remaining_qalys
        avoided = abs(excess) * 40_000               # excess lifetime medical cost fraction
        intervention = 1_500 if excess > 0 else 0    # lifestyle cost to reverse
        if excess > 0:
            label = (f"Biological age +{accel:.1f} yr — {excess*100:.1f} pts excess 10-yr "
                     f"mortality risk, reversible via lifestyle")
        else:
            label = (f"Biological age {accel:.1f} yr — {abs(excess)*100:.1f} pts lower 10-yr "
                     f"mortality risk (value already banked)")
        if qaly > 0.001 or avoided > 1:
            add("Biological aging", label, avoided, qaly, intervention, "high",
                "Derived from your PhenoAge 10-year mortality risk vs the baseline for your "
                "chronological age. PhenoAge is a mortality-calibrated clock — each year of "
                "acceleration ≈ 9% higher all-cause mortality (HR 1.09/yr; Levine, Aging 2018).")

    # ── Addiction genetics (personal) ──
    if addiction_result and addiction_result.get("available"):
        composite = addiction_result.get("composite", {})
        for cat_name, econ in ADDICTION_ECONOMICS.items():
            tier = composite.get(f"{cat_name.lower()}_tier",
                                composite.get("overall_tier", ""))
            if tier in ("elevated", "high", "Elevated", "High"):
                avoided = econ["outcome_value"] * 0.12
                add("Addiction Genetics", f"{cat_name} susceptibility — {econ['finding']}",
                    avoided, econ["qaly_gain"] * 0.3, econ["cost"], "low",
                    f"{cat_name} tier {tier} — {econ['clinical_benefit']} "
                    f"({econ.get('src', '')}).")
                break

    # ── Neurochemistry (personal) ──
    if neurochemistry_result and neurochemistry_result.get("available"):
        composite = neurochemistry_result.get("composite", {})
        comt = composite.get("comt_class", "")
        if comt and comt != "normal":
            econ = NEUROCHEMISTRY_ECONOMICS.get("COMT", {})
            if econ:
                avoided = econ["outcome_value"] * 0.10
                add("Neurochemistry", f"COMT {comt} — {econ['finding']}",
                    avoided, econ["qaly_gain"] * 0.3, econ["cost"], "low",
                    f"COMT class {comt} — {econ['clinical_benefit']} ({econ.get('src', '')}).")

    # ── Mendelian randomization (personal) ──
    if mr_result:
        for f in (mr_result.get("findings") or [])[:3]:
            if f.get("status") != "computed":
                continue
            rr = f.get("outcome_relative_risk")
            if rr is None or rr <= 1.05:
                continue
            outcome = f.get("outcome", "")
            exposure = f.get("exposure", "")
            oc = None
            for key, val in MR_OUTCOME_COSTS.items():
                if key.lower() in outcome.lower():
                    oc = val
                    break
            if oc is None:
                continue
            excess_risk = (rr - 1.0) * oc["baseline_10yr"]
            avoided = round(excess_risk * oc["cost"] * 0.3)
            qaly = round(excess_risk * oc["qaly_loss"] * 0.3, 2)
            if avoided > 100:
                add("Mendelian Randomization",
                    f"MR: {exposure} → {outcome} (RR {rr:.2f})",
                    avoided, qaly, 500, "moderate",
                    f"Causal projection: reducing {exposure} lowers {outcome} risk "
                    f"(MR RR {rr:.2f}, {f.get('n_used', 0)} instruments).")

    # ── Clinical variants (personal) ──
    if clinical_variants_result and clinical_variants_result.get("available"):
        for v in (clinical_variants_result.get("findings") or []):
            sig = (v.get("significance") or "").lower()
            if "pathogenic" not in sig:
                continue
            gene = v.get("gene", "")
            econ = ACMG_GENE_ECONOMICS.get(gene)
            if econ is None:
                continue
            zyg = v.get("zygosity", "het")
            avoided = econ["outcome_value"] * (0.50 if "homo" in zyg else 0.25)
            add("Clinical Variant", f"{econ['finding']} ({zyg})",
                avoided, econ["qaly_gain"], econ["cost"], "high",
                f"ClinVar {sig} in {gene} — {econ['clinical_benefit']} "
                f"({econ.get('src', '')}).")

    # ── Family planning (personal) — NOT MONETISED ──
    # This block used to add $5,000 and 0.1 QALY per carrier condition. That
    # contradicted the policy stated in value_of_information.NOT_VALUED, which
    # says reproductive findings are deliberately never given a dollar value —
    # two parts of the same report were answering the same question
    # differently, and the one with a number was winning.
    #
    # The reasoning in NOT_VALUED holds here: attaching a figure to an affected
    # birth prices a prospective child, and the $5,000 rate embedded an uptake
    # assumption that is really a reproductive preference. The finding is
    # surfaced as a decision to consider, and the carrier-panel recommendation
    # driven by runs of homozygosity is reported categorically alongside it.
    if family_planning_result and family_planning_result.get("available"):
        n_actionable = family_planning_result.get("n_actionable", 0)
        if n_actionable > 0:
            not_monetised.append({
                "category": "Family Planning",
                "finding": f"Reproductive genetics — {n_actionable} carrier "
                           f"condition(s) where partner testing would be "
                           f"informative",
                "decision": "Consider partner carrier testing and genetic "
                            "counselling before conception.",
                "indicative_cost": PARTNER_TESTING_COST + GENETIC_COUNSELING_COST,
                "reason": "Reproductive outcomes are deliberately not "
                          "monetised. The cost of testing is shown because it "
                          "is a real price; no benefit figure is attached "
                          "because valuing one would price a prospective "
                          "child and embed one set of reproductive "
                          "preferences as if it were universal.",
            })

    # ── PheWAS extreme predictions (personal) ──
    if phewas_result:
        n_extreme = len(phewas_result.get("headline") or [])
        if n_extreme > 0:
            for h in (phewas_result.get("headline") or [])[:3]:
                cat = h.get("category", "")
                econ = PHEWAS_CATEGORY_ECONOMICS.get(cat)
                if econ is None:
                    continue
                avoided = econ["outcome_value"] * 0.08
                add("PheWAS Biomarker",
                    f"{h.get('trait', cat)} extreme prediction ({h.get('tier', '')})",
                    avoided, econ["qaly_gain"] * 0.3, econ["cost"], "low",
                    f"Genetically predicted {h.get('trait', '')} at "
                    f"p{h.get('percentile', '')} — {econ['clinical_benefit']} "
                    f"({econ.get('src', '')}).")

    # ── Wellness predictions (personal) ──
    if wellness_result:
        n_actionable_wellness = 0
        for pred in (wellness_result.get("predictions") or []):
            result_text = (pred.get("result") or "").lower()
            if any(k in result_text for k in ("reduced", "impaired", "low",
                                               "poor", "deficien")):
                n_actionable_wellness += 1
        if n_actionable_wellness > 0:
            avoided = n_actionable_wellness * 2_000
            add("Wellness Genetics",
                f"{n_actionable_wellness} actionable wellness variant(s)",
                avoided, 0.05 * n_actionable_wellness, 300, "low",
                "Nutrient metabolism, sleep, and fitness optimization variants "
                "with published intervention economics.")

    # ── Deduplicate personal items: same condition from multiple sections ──
    seen_keys: dict[str, int] = {}
    deduped_items: list[dict] = []
    for item in items:
        key = item["finding"].lower().split("—")[0].strip()
        if key in seen_keys:
            idx = seen_keys[key]
            if item["net"] > deduped_items[idx]["net"]:
                deduped_items[idx] = item
            continue
        seen_keys[key] = len(deduped_items)
        deduped_items.append(item)
    items = deduped_items

    # CORRELATED-TARGET POOLING. The exact-key dedup above only catches items
    # that are literally identical. It does not catch the same genotype
    # surfacing through two different panels — a COMT line from the
    # neurochemistry module and a COMT-guided-prescribing line from the
    # pharmacogenomics module were both counted in full, which on one measured
    # genome valued that single variant at $16,550. Same for a carrier result
    # appearing once as a carrier finding and again as symptom awareness.
    #
    # Rank the lines that share a target, keep the strongest at full value and
    # discount the rest. The discount is recorded on the item rather than
    # applied invisibly, so the report can show which lines were reduced.
    try:
        from . import engine as _ee
        # Match against the gene symbols THIS module actually emits, rather
        # than a regex for gene-shaped words — the latter reads "MI", "MACE"
        # and "B12" as genes and invents targets that do not exist.
        _vocab = frozenset(_ee.DEFAULT_GENE_VOCABULARY) | frozenset(
            k.split("*")[0].split(":")[0].upper()
            for d in (ACMG_GENE_ECONOMICS, PGX_ECONOMICS, HLA_ECONOMICS,
                      NEUROCHEMISTRY_ECONOMICS, METAL_OXIDATIVE_ECONOMICS,
                      IMMUNOGENETICS_ECONOMICS, DETOX_ECONOMICS)
            for k in d
            if k.replace("-", "").replace("*", "").replace(":", "").isalnum())
        items = _ee.deduplicate_by_target(items, value_key="net",
                                          text_key="finding",
                                          fallback_key="category",
                                          vocabulary=_vocab)
        items = [i for i in items if i.get("net", 0) or i.get("qaly", 0)]
        pooled_targets = sorted({i["pool_target"] for i in items
                                 if i.get("pool_rank", 0) > 0})
    except Exception:
        # The economics must still produce a number if the engine is missing;
        # it just produces the older, higher one, and says so via the flag.
        pooled_targets = []
        _pooling_applied = False
    else:
        _pooling_applied = True

    total_avoided = sum(i["avoided"] for i in items)
    total_qaly = sum(i["qaly"] for i in items)
    total_qaly_value = sum(i["qaly_value"] for i in items)
    total_intervention = sum(i["intervention"] for i in items)
    total_net = sum(i["net"] for i in items)
    gross_benefit = total_avoided + total_qaly_value
    roi = round((total_net) / _ANALYSIS_COST, 1) if _ANALYSIS_COST else None

    # COST-SAVING vs COST-EFFECTIVE (honesty split).
    # total_net is NET MONETARY BENEFIT — dominated by monetised QALYs (health value
    # you'd pay for), NOT cash returned. Calling total_net/cost "ROI" implies money
    # back. Separate the two: the CASH side is averted cost minus what you spend;
    # only if THAT is positive is the analysis genuinely cost-saving. Most prevention
    # is cost-effective (great value per QALY) but cost-adding (net cash out), and the
    # report must not blur the two.
    net_cash = total_avoided - total_intervention          # money only, no QALYs
    cash_roi = (round(net_cash / _ANALYSIS_COST, 1)
                if _ANALYSIS_COST else None)
    value_to_cost_ratio = roi                              # NMB per $ — mostly health
    if net_cash > 0:
        verdict = "cost-saving"                            # genuine money back
    elif total_net > 0:
        verdict = "cost-effective (adds cost, worth it per QALY)"
    else:
        verdict = "not cost-effective at this threshold"

    # The counterfactual is "these findings, pooled the same way, with everyone
    # following through" — so it is derived from the post-pooling net rather
    # than from a figure captured before the correlated-target discount, which
    # would fold the pooling correction into the adherence one.
    _efficacy_net = sum(i["net"] / (i.get("adherence") or 1.0) for i in items)

    items.sort(key=lambda i: -i["net"])
    return {
        "available": bool(items),
        "n_items": len(items),
        "horizon_years": PERSONAL_HORIZON_YEARS,
        "items": items,
        "not_monetised": not_monetised,
        "n_not_monetised": len(not_monetised),
        # Real-world adherence, reported rather than applied invisibly. The
        # efficacy figure is what this page used to headline; keeping both
        # means the discount is visible instead of just making the number
        # smaller for no stated reason.
        "adherence_applied": any(i.get("adherence", 1.0) < 1.0 for i in items),
        "efficacy_net": round(_efficacy_net),
        "adherence_drag": round(_efficacy_net - total_net),
        "mean_adherence": (
            round(sum(i.get("adherence", 1.0) for i in items) / len(items), 3)
            if items else 1.0),
        # Correlated-target pooling, reported rather than applied invisibly.
        "pooling_applied": _pooling_applied,
        "pooled_targets": pooled_targets,
        "n_pooled_targets": len(pooled_targets),
        "total_avoided": round(total_avoided),
        "total_qaly": round(total_qaly, 2),
        "total_qaly_value": round(total_qaly_value),
        "total_intervention": round(total_intervention),
        "total_net": round(total_net),
        "net_low": round(total_net * 0.5),          # ±50% sensitivity band
        "net_high": round(total_net * 1.5),
        "gross_benefit": round(gross_benefit),
        # Cash side, kept explicitly separate from monetised health value.
        "net_cash": round(net_cash),
        "cash_roi": cash_roi,
        "value_to_cost_ratio": value_to_cost_ratio,
        "is_cost_saving": net_cash > 0,
        "verdict": verdict,
        "top_preventable": (max(items, key=lambda i: i["avoided"])["finding"]
                            if items else None),
        "analysis_cost": _ANALYSIS_COST,
        "roi": roi,
        "value_per_qaly": VALUE_PER_QALY,
        "disclaimer": DISCLAIMER,
    }


def build_cost_consequence_analysis(econ: dict | None) -> dict:
    """**Cost-consequence analysis (CCA)** — the disaggregated table.

    *In plain English:* a cost-effectiveness analysis crushes everything into one
    number (cost per QALY). A cost-consequence analysis deliberately does not: it
    lays out each cost and each consequence side by side, in its own natural units,
    and lets the reader weigh them. Nothing is hidden inside a ratio.

    NICE prefers this format for digital health precisely because such tools produce
    several different kinds of benefit (some monetary, some health, some neither),
    and collapsing them into a single ICER discards information the decision-maker
    needs. It is also the honest format when the evidence is too thin to justify a
    confident single number — which is the case here.

    Deliberately reports NO summary ratio: that is the point of the format.
    Ref: NICE Evidence Standards Framework for digital health technologies;
    cost-consequence analysis (health-economics-metrics topic).
    """
    if not econ or not econ.get("available"):
        return {"available": False,
                "reason": "No modelled economic items — needs actionable findings "
                          "and/or a blood-work panel."}

    items = econ.get("items", [])
    # Disaggregated consequences, each in its OWN unit — never summed together.
    rows = [
        {"kind": "Cost", "measure": "Testing / analysis cost",
         "value": econ.get("analysis_cost", 0), "unit": "$ one-off",
         "note": "What you pay to obtain the information."},
        {"kind": "Cost", "measure": "Downstream intervention cost",
         "value": econ.get("total_intervention", 0), "unit": f"$ over {econ.get('horizon_years','—')} yr",
         "note": "Screening, prophylaxis, lifestyle programmes triggered by findings."},
        {"kind": "Consequence", "measure": "Medical costs averted (marginal)",
         "value": econ.get("total_avoided", 0), "unit": f"$ over {econ.get('horizon_years','—')} yr",
         "note": "Marginal, not average, cost of the cases avoided."},
        {"kind": "Consequence", "measure": "Quality-adjusted life-years gained",
         "value": econ.get("total_qaly", 0), "unit": "QALYs",
         "note": "Health gain in its own unit — deliberately NOT converted here."},
        {"kind": "Consequence", "measure": "Actionable findings identified",
         "value": len(items), "unit": "findings",
         "note": "Count of distinct modelled findings driving the analysis."},
        {"kind": "Consequence", "measure": "High-confidence findings",
         "value": sum(1 for i in items if i.get("confidence") == "high"),
         "unit": "findings",
         "note": "Subset resting on stronger evidence."},
    ]
    # The cash bottom line is reported as a fact, still not blended with QALYs.
    net_cash = econ.get("net_cash", 0)
    return {
        "available": True,
        "rows": rows,
        "net_cash": net_cash,
        "total_qaly": econ.get("total_qaly", 0),
        "verdict": econ.get("verdict", ""),
        "horizon_years": econ.get("horizon_years"),
        "plain_english": (
            f"Over {econ.get('horizon_years','—')} years this profile models "
            f"{econ.get('total_qaly', 0)} quality-adjusted life-years gained and a net "
            f"cash position of {_money(net_cash)}, from {len(items)} findings. Costs and "
            f"consequences are listed separately, in their own units, rather than "
            f"collapsed into a single cost-per-QALY figure — so you can weigh them "
            f"yourself instead of trusting one ratio."),
        "why_this_format": (
            "Cost-consequence analysis is NICE's preferred presentation for digital "
            "health technologies: it keeps monetary and non-monetary consequences "
            "visible side by side instead of hiding them inside an ICER. It is also "
            "the appropriate format when the underlying parameters are illustrative."),
        "caveat": ("No summary ratio is given by design. Figures are modelled on "
                   "illustrative parameters — not measurements, not medical advice."),
    }


def _render_cca_html(cca: dict | None) -> str:
    """Render the cost-consequence table (NICE's preferred digital-health format)."""
    if not cca or not cca.get("available"):
        return ""
    rows = ""
    for r in cca["rows"]:
        is_cost = r["kind"] == "Cost"
        col = "#8a6100" if is_cost else "#1a7f37"
        val = r["value"]
        shown = _money(val) if str(r["unit"]).startswith("$") else f"{val}"
        rows += (
            f'<tr>'
            f'<td style="padding:6px 10px;color:{col};font-weight:600">{r["kind"]}</td>'
            f'<td style="padding:6px 10px">{_esc_econ(r["measure"])}</td>'
            f'<td style="padding:6px 10px;text-align:right;font-weight:700;'
            f'font-variant-numeric:tabular-nums">{shown}</td>'
            f'<td style="padding:6px 10px;color:#8a94a3">{_esc_econ(r["unit"])}</td>'
            f'<td style="padding:6px 10px;color:#5b6673;font-size:.92em">{_esc_econ(r["note"])}</td>'
            f'</tr>')
    return f"""
    <div style="margin:16px 0;border:1px solid #dbe3ec;border-radius:10px;padding:12px 14px">
      <div style="font-weight:700;color:#12467a">Cost-consequence analysis
        <span style="font-weight:400;color:#8a94a3;font-size:.85em">
        — costs and consequences, disaggregated</span></div>
      <div style="font-size:.86em;color:#5b6673;margin:4px 0 8px">
        {_esc_econ(cca['plain_english'])}</div>
      <table style="width:100%;border-collapse:collapse;font-size:.88em">
        <thead><tr style="text-align:left;color:#8a94a3;border-bottom:1px solid #dbe3ec">
          <th style="padding:6px 10px">Type</th><th style="padding:6px 10px">Measure</th>
          <th style="padding:6px 10px;text-align:right">Value</th>
          <th style="padding:6px 10px">Unit</th><th style="padding:6px 10px">Note</th>
        </tr></thead>
        <tbody>{rows}</tbody>
      </table>
      <div style="margin-top:8px;font-size:.82em;color:#6b5a1e;background:#fff8e6;
                  border:1px solid #f0e0a8;border-radius:6px;padding:8px 10px">
        <strong>No single cost-per-QALY figure is given, on purpose.</strong>
        {_esc_econ(cca['why_this_format'])} {_esc_econ(cca['caveat'])}
      </div>
    </div>"""


def _esc_econ(s) -> str:
    s = "" if s is None else str(s)
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _render_not_monetised_html(rows: list[dict] | None) -> str:
    """Findings that carry a real decision but deliberately no dollar figure.

    Showing these matters. If a reproductive finding simply vanished from the
    economic page, the omission would look like an oversight; stated as a
    choice, with the reason attached, it reads as what it is.
    """
    if not rows:
        return ""
    items = "".join(f"""
      <div style="border-top:1px solid #eef1f4;padding:9px 0">
        <div style="font-weight:600;color:#2b3440">{_esc_econ(r.get('finding',''))}</div>
        <div style="font-size:.88em;color:#48545f;margin-top:3px">
          <strong>Decision:</strong> {_esc_econ(r.get('decision',''))}</div>
        <div style="font-size:.85em;color:#6a7683;margin-top:3px">
          Indicative cost of acting: {_money(r.get('indicative_cost', 0))}
          &middot; no benefit figure attached</div>
        <div style="font-size:.82em;color:#8a94a3;margin-top:4px;font-style:italic">
          {_esc_econ(r.get('reason',''))}</div>
      </div>""" for r in rows)
    return f"""
    <div style="border:1px solid #e3e7ec;border-left:4px solid #8a5cf6;border-radius:10px;
                padding:13px 15px;margin:14px 0;background:#fbfcfe">
      <div style="display:flex;justify-content:space-between;align-items:baseline;
                  gap:12px;flex-wrap:wrap">
        <div style="font-weight:700;color:#8a5cf6">Reported without a dollar figure</div>
        <div style="font-size:.72em;color:#8a94a3;border:1px solid #dfe4ea;
                    border-radius:20px;padding:2px 9px;white-space:nowrap">
          deliberately not monetised</div>
      </div>
      <div style="font-size:.86em;color:#48545f;margin-top:5px">
        These findings are excluded from every total on this page. That is a
        modelling decision, not a gap.</div>
      {items}
    </div>"""


def _render_adherence_basis_html(econ: dict) -> str:
    """State that the totals above are effectiveness, not efficacy.

    Without this the page reports a smaller number than it used to with no
    stated reason, which is its own kind of dishonesty — and the pooled payer
    analysis elsewhere in the report would be describing the same genome on a
    different basis.
    """
    if not econ.get("adherence_applied"):
        return ""
    drag = econ.get("adherence_drag", 0)
    eff = econ.get("efficacy_net", 0)
    pct = round(100.0 * drag / eff, 1) if eff else 0.0
    return f"""
    <div style="margin-top:10px;padding:9px 12px;background:#f7fafd;
                border:1px solid #d8e2ee;border-radius:8px;font-size:.85em;
                color:#42566b">
      <strong>These are real-world figures, not trial figures.</strong>
      Every benefit above has been multiplied by the share of people who
      actually keep doing the thing — roughly half for daily preventive
      medication, less for sustained behaviour change, more for a screening
      appointment you attend once. Ongoing intervention costs are discounted by
      the same factor, because someone who stops taking a statin stops paying
      for it. At full adherence these findings would total
      <strong>{_money(eff)}</strong>; charging realistic adherence
      (mean {econ.get('mean_adherence', 1.0):.0%}) removes
      <strong>{_money(drag)}</strong> ({pct}%).
      <div style="font-size:.92em;color:#7b8794;margin-top:4px">
        WHO (2003), Adherence to Long-Term Therapies. Screening uptake and
        behavioural maintenance are declared assumptions.</div>
    </div>"""


def render_economic_analysis_html(econ: dict, file_label: str = "") -> str:
    """Standalone economic-impact sheet (economic_analysis.html)."""
    if not econ or not econ.get("available"):
        body = ("<p>No modeled economic-impact items for this profile — this needs "
                "actionable genomic findings and/or a supplied blood-work panel "
                "(<code>--bloodwork labs.json</code>).</p>")
    else:
        conf_col = {"high": "#1a7f37", "moderate": "#8a6100", "low": "#8a94a3"}
        rows = ""
        for i in econ["items"]:
            c = conf_col.get(i["confidence"], "#8a94a3")
            net_col = "#1a7f37" if i["net"] >= 0 else "#b3261e"
            rows += f"""
            <tr>
              <td><strong>{_esc_econ(i['finding'])}</strong>
                  <div style="color:#8a94a3;font-size:.82em">{_esc_econ(i['category'])} ·
                  <span style="color:{c}">{_esc_econ(i['confidence'])} confidence</span></div>
                  <div style="color:#9aa4b0;font-size:.78em">{_esc_econ(i['basis'])}</div></td>
              <td style="text-align:right;font-variant-numeric:tabular-nums">{_money(i['avoided'])}</td>
              <td style="text-align:right;font-variant-numeric:tabular-nums">{i['qaly']:g}<br>
                  <span style="color:#8a94a3;font-size:.8em">{_money(i['qaly_value'])}</span></td>
              <td style="text-align:right;font-variant-numeric:tabular-nums;color:#8a6100">{_money(i['intervention'])}</td>
              <td style="text-align:right;font-variant-numeric:tabular-nums;font-weight:700;color:{net_col}">{_money(i['net'])}</td>
            </tr>"""
        roi = econ.get("roi")
        net = econ["total_net"]
        net_col = "#1a7f37" if net >= 0 else "#b3261e"

        # SVG contribution chart — net value per finding (green +, red −)
        chart = ""
        ch_items = [i for i in econ["items"] if i["net"]]
        if ch_items:
            maxabs = max(abs(i["net"]) for i in ch_items) or 1
            barW, rowH, midX = 260, 26, 250
            svg_rows = ""
            for idx, i in enumerate(ch_items):
                y = idx * rowH + 4
                w = barW * abs(i["net"]) / maxabs
                col = "#2fae57" if i["net"] >= 0 else "#e0524a"
                x = midX if i["net"] >= 0 else midX - w
                lab = (i["finding"][:34] + "…") if len(i["finding"]) > 35 else i["finding"]
                svg_rows += (
                    f'<text x="{midX-8}" y="{y+15}" font-size="10.5" fill="#33404d" '
                    f'text-anchor="end">{_esc_econ(lab)}</text>'
                    f'<rect x="{x:.0f}" y="{y}" width="{w:.0f}" height="16" rx="3" fill="{col}"/>'
                    f'<text x="{(x+w+5) if i["net"]>=0 else (x-5):.0f}" y="{y+13}" font-size="10" '
                    f'fill="{col}" text-anchor="{"start" if i["net"]>=0 else "end"}">{_money(i["net"])}</text>')
            H = len(ch_items) * rowH + 10
            chart = (
                f'<div style="margin:14px 0"><div style="font-weight:700;color:#12467a;'
                f'margin-bottom:4px">Net value by finding</div>'
                f'<svg viewBox="0 0 560 {H}" width="100%" style="max-width:640px" '
                f'xmlns="http://www.w3.org/2000/svg">'
                f'<line x1="{midX}" y1="0" x2="{midX}" y2="{H}" stroke="#dbe3ec"/>{svg_rows}</svg></div>')

        rng = (f'<div style="color:#5b6673;font-size:.85em;margin-top:2px">'
               f'plausible range {_money(econ["net_low"])} – {_money(econ["net_high"])} (±50%)</div>')
        top_prev = (f'<div style="color:#8a94a3;font-size:.82em;margin-top:6px">'
                    f'Largest avoidable cost: <strong>{_esc_econ(econ.get("top_preventable",""))}</strong></div>'
                    if econ.get("top_preventable") else "")
        body = f"""
        {chart}
        <div style="display:flex;gap:18px;flex-wrap:wrap;align-items:center;
             background:linear-gradient(135deg,#f2f9f4,#eef4fb);border:1px solid #dbe8dd;
             border-radius:14px;padding:18px 22px;margin:14px 0">
          <div style="text-align:center;min-width:200px">
            <div style="font-size:.8em;color:#5b6673;text-transform:uppercase;letter-spacing:.04em">
              Modeled {econ['horizon_years']}-yr net value</div>
            <div style="font-size:2.8em;font-weight:800;color:{net_col}">{_money(net)}</div>
            <div style="color:#8a94a3;font-size:.82em">of acting on your results</div>
            {rng}
          </div>
          <div style="flex:1;min-width:240px;display:grid;grid-template-columns:1fr 1fr;gap:10px">
            <div><div style="font-size:.78em;color:#8a94a3">Medical cost avoided</div>
              <div style="font-size:1.4em;font-weight:700">{_money(econ['total_avoided'])}</div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Quality-of-life value</div>
              <div style="font-size:1.4em;font-weight:700">{_money(econ['total_qaly_value'])}
              <span style="font-size:.55em;color:#8a94a3">({econ['total_qaly']:g} QALYs)</span></div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Intervention cost</div>
              <div style="font-size:1.4em;font-weight:700;color:#8a6100">{_money(econ['total_intervention'])}</div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Value per $ of analysis
              (mostly health, not cash)</div>
              <div style="font-size:1.4em;font-weight:700;color:{net_col}">{roi}×</div></div>
            <div><div style="font-size:.78em;color:#8a94a3">Net <em>cash</em>
              (averted cost − spend)</div>
              <div style="font-size:1.4em;font-weight:700;color:{"#1a7f37" if econ.get('is_cost_saving') else "#8a6100"}">{_money(econ.get('net_cash', 0))}</div></div>
          </div>
          <div style="margin-top:10px;padding:8px 12px;background:#fff8e6;
                      border:1px solid #f0e0a8;border-radius:8px;font-size:.85em;color:#6b5a1e">
            <strong>Verdict: {_esc_econ(econ.get('verdict', ''))}.</strong>
            Most of the headline figure is <em>monetised health</em> (quality-adjusted
            life-years valued at {_money(econ.get('value_per_qaly', 0))}/QALY), not money
            returned to you. Prevention is typically cost-<em>effective</em> — excellent
            value per healthy year — while still costing more cash than it saves. The
            two are shown separately so neither is mistaken for the other.
          </div>
        </div>
        {_render_adherence_basis_html(econ)}
        {_render_cca_html(build_cost_consequence_analysis(econ))}
        {_render_not_monetised_html(econ.get("not_monetised"))}
        {top_prev}
        <table style="width:100%;border-collapse:collapse;font-size:.9em;margin-top:8px">
          <thead><tr style="background:#f7f9fb;text-align:right">
            <th style="text-align:left;padding:8px 10px">Finding</th>
            <th style="padding:8px 10px">Cost avoided</th>
            <th style="padding:8px 10px">QALY (value)</th>
            <th style="padding:8px 10px">Intervention</th>
            <th style="padding:8px 10px">Net value</th></tr></thead>
          <tbody>{rows}</tbody>
        </table>"""

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Economic-Impact Analysis{(' — ' + _esc_econ(file_label)) if file_label else ''}</title>
<style>
body{{font-family:-apple-system,"Segoe UI",Roboto,sans-serif;color:#1a1f26;max-width:1000px;
 margin:24px auto;padding:0 16px;line-height:1.5}}
h1{{font-size:1.7em;border-bottom:2px solid #12467a;padding-bottom:6px;color:#12467a}}
td,th{{padding:8px 10px;border-bottom:1px solid #eef1f4;vertical-align:top}}
.disc{{background:#fff8e6;border:1px solid #f0e0b0;border-radius:8px;padding:10px 14px;
 font-size:.85em;color:#6a5b2a;margin:12px 0}}
</style></head><body>
<h1>Economic-Impact Analysis</h1>
<p style="color:#667">A model of the 10-year economic impact of acting on your genomic and
blood-work results — expected medical-cost avoidance plus monetised quality-of-life gains,
net of intervention cost.</p>
<div class="disc">⚠️ <strong>Illustrative model, not financial or medical advice.</strong>
Figures use population-average cost-of-illness and risk-reduction estimates applied to your
findings; individual outcomes vary widely. Quality-of-life gains are monetised at the standard
${econ.get('value_per_qaly', 100000):,}/QALY threshold. Not a guarantee of savings.</div>
{body}
<p style="color:#9aa4b0;font-size:.8em;margin-top:16px">{_esc_econ(econ.get('disclaimer',''))}</p>
</body></html>"""
