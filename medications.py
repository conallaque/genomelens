"""
Medication Lookup
=================

Cross-references a user-specified medication list against PGx phenotypes
already calculated in pgx.py. Generates a focused "Medication Review"
section in the report.

Drug-name normalisation: brand → generic. The internal lookup table maps
each drug to (metabolizing_gene, optional_pathway_note). When the
phenotype is found in pgx_result, we surface the recommendation;
otherwise we emit a "not flagged by PGx" note explaining why.
"""

from __future__ import annotations

from typing import Dict, List, Optional


# Brand → generic. Multiple brands may map to the same generic.
BRAND_TO_GENERIC: Dict[str, str] = {
    # SSRIs/SNRIs/antidepressants
    "zoloft": "sertraline",
    "lexapro": "escitalopram",
    "cipralex": "escitalopram",
    "celexa": "citalopram",
    "cipramil": "citalopram",
    "prozac": "fluoxetine",
    "sarafem": "fluoxetine",
    "paxil": "paroxetine",
    "seroxat": "paroxetine",
    "luvox": "fluvoxamine",
    "effexor": "venlafaxine",
    "cymbalta": "duloxetine",
    "wellbutrin": "bupropion",
    "elavil": "amitriptyline",
    "tofranil": "imipramine",
    "norpramin": "desipramine",
    "remeron": "mirtazapine",
    # ADHD stimulants
    "vyvanse": "lisdexamfetamine",
    "elvanse": "lisdexamfetamine",
    "adderall": "amphetamine",
    "mydayis": "amphetamine",
    "dexedrine": "dextroamphetamine",
    "ritalin": "methylphenidate",
    "concerta": "methylphenidate",
    "focalin": "dexmethylphenidate",
    "strattera": "atomoxetine",
    "intuniv": "guanfacine",
    "kapvay": "clonidine",
    # Statins
    "lipitor": "atorvastatin",
    "zocor": "simvastatin",
    "crestor": "rosuvastatin",
    "pravachol": "pravastatin",
    "livalo": "pitavastatin",
    "mevacor": "lovastatin",
    "lescol": "fluvastatin",
    # Anticoagulants / antiplatelets
    "coumadin": "warfarin",
    "jantoven": "warfarin",
    "plavix": "clopidogrel",
    "effient": "prasugrel",
    "brilinta": "ticagrelor",
    "eliquis": "apixaban",
    "xarelto": "rivaroxaban",
    "pradaxa": "dabigatran",
    "savaysa": "edoxaban",
    # PPIs
    "prilosec": "omeprazole",
    "nexium": "esomeprazole",
    "prevacid": "lansoprazole",
    "protonix": "pantoprazole",
    "aciphex": "rabeprazole",
    # Pain / opioids
    "tylenol": "acetaminophen",
    "codeine": "codeine",
    "ultram": "tramadol",
    "vicodin": "hydrocodone",
    "norco": "hydrocodone",
    "percocet": "oxycodone",
    "oxycontin": "oxycodone",
    "advil": "ibuprofen",
    "motrin": "ibuprofen",
    "aleve": "naproxen",
    "celebrex": "celecoxib",
    # Antipsychotics
    "abilify": "aripiprazole",
    "zyprexa": "olanzapine",
    "seroquel": "quetiapine",
    "risperdal": "risperidone",
    "geodon": "ziprasidone",
    "clozaril": "clozapine",
    "haldol": "haloperidol",
    # Misc
    "tegretol": "carbamazepine",
    "dilantin": "phenytoin",
    "lamictal": "lamotrigine",
    "depakote": "valproate",
    "imuran": "azathioprine",
    "purinethol": "6-mercaptopurine",
    "tamoxifen": "tamoxifen",
    "ondansetron": "ondansetron",
    "zofran": "ondansetron",
    "metoprolol": "metoprolol",
    "lopressor": "metoprolol",
    "toprol": "metoprolol",
    "coreg": "carvedilol",
    "tacrolimus": "tacrolimus",
    "prograf": "tacrolimus",
    "abacavir": "abacavir",
    "ziagen": "abacavir",
}


# Generic → list of (gene, pathway_note). Multi-gene drugs list all relevant
# metabolizers/transporters. Note: this table reflects CPIC-prioritised
# gene-drug pairs, not every metabolic pathway.
DRUG_GENES: Dict[str, List[Dict[str, str]]] = {
    # SSRIs / SNRIs
    "sertraline":      [{"gene": "CYP2C19", "note": "Primary metabolizer; PMs/IMs may need lower dose."}],
    "escitalopram":    [{"gene": "CYP2C19", "note": "Primary metabolizer; PMs accumulate, UMs underexposed."}],
    "citalopram":      [{"gene": "CYP2C19", "note": "Primary metabolizer; QT-prolongation concern at high doses in PMs."}],
    "fluoxetine":      [{"gene": "CYP2D6",  "note": "Metaboliser + strong 2D6 inhibitor itself."}],
    "paroxetine":      [{"gene": "CYP2D6",  "note": "Strong 2D6 substrate AND inhibitor."}],
    "fluvoxamine":     [{"gene": "CYP2D6",  "note": "2D6 substrate."}],
    "venlafaxine":     [{"gene": "CYP2D6",  "note": "Primary metaboliser; PMs accumulate parent drug."}],
    "duloxetine":      [{"gene": "CYP2D6",  "note": "2D6 substrate; CYP1A2 also contributes."}],
    "bupropion":       [{"gene": "CYP2B6",  "note": "Primary metaboliser; PMs higher seizure risk."}],
    "amitriptyline":   [{"gene": "CYP2D6", "note": "Primary metaboliser."}, {"gene": "CYP2C19", "note": "Demethylation."}],
    "imipramine":      [{"gene": "CYP2D6", "note": "Primary."}, {"gene": "CYP2C19", "note": "Demethylation."}],
    "desipramine":     [{"gene": "CYP2D6", "note": "Primary."}],
    "mirtazapine":     [{"gene": "CYP2D6", "note": "Minor 2D6 contribution."}],

    # Stimulants
    "lisdexamfetamine":[{"gene": "CYP2D6", "note": "Modest role in d-amphetamine clearance."},
                        {"gene": "COMT",   "note": "Dopamine clearance — dose-response modifier."},
                        {"gene": "ANKK1/DRD2", "note": "Reward-pathway response."}],
    "amphetamine":     [{"gene": "CYP2D6", "note": "Modest role."},
                        {"gene": "COMT",   "note": "Dopamine clearance modifier."},
                        {"gene": "ANKK1/DRD2", "note": "Stimulant response."}],
    "dextroamphetamine":[{"gene": "CYP2D6", "note": "Modest role."},
                         {"gene": "COMT",   "note": "Dopamine clearance modifier."}],
    "methylphenidate": [{"gene": "COMT",   "note": "Slow COMT (Met/Met) may respond differently."},
                        {"gene": "ANKK1/DRD2", "note": "Stimulant response variation."}],
    "dexmethylphenidate":[{"gene": "COMT", "note": "Dopamine clearance modifier."}],
    "atomoxetine":     [{"gene": "CYP2D6", "note": "Primary clearance — PMs need lower dose."}],
    "guanfacine":      [{"gene": "CYP3A4",  "note": "Primary metaboliser."}],
    "clonidine":       [{"gene": "CYP2D6",  "note": "Minor role."}],

    # Statins
    "atorvastatin":    [{"gene": "SLCO1B1", "note": "Hepatic uptake; reduced function increases myopathy."},
                        {"gene": "CYP3A4",  "note": "Primary metaboliser."}],
    "simvastatin":     [{"gene": "SLCO1B1", "note": "STRONGLY affected by *5 (rs4149056) — myopathy risk; avoid 80 mg."}],
    "rosuvastatin":    [{"gene": "SLCO1B1", "note": "Affected but less so than simvastatin."}],
    "pravastatin":     [{"gene": "SLCO1B1", "note": "Affected but less so than simvastatin."}],
    "pitavastatin":    [{"gene": "SLCO1B1", "note": "Affected by SLCO1B1 variants."}],
    "lovastatin":      [{"gene": "SLCO1B1", "note": "Affected."},
                        {"gene": "CYP3A4",  "note": "Primary metaboliser."}],
    "fluvastatin":     [{"gene": "CYP2C9",  "note": "Primary metaboliser."}],

    # Anticoagulants / antiplatelets
    "warfarin":        [{"gene": "VKORC1",  "note": "Target enzyme — A allele = lower dose requirement."},
                        {"gene": "CYP2C9",  "note": "Primary metaboliser; *2/*3 reduce clearance."}],
    "clopidogrel":     [{"gene": "CYP2C19", "note": "REQUIRED for activation; PMs get minimal effect — use prasugrel/ticagrelor."}],
    "prasugrel":       [{"gene": "—",        "note": "Not heavily CYP2C19-dependent; safer alternative to clopidogrel for PMs."}],
    "ticagrelor":      [{"gene": "—",        "note": "Not CYP2C19-activated; safer alternative for CYP2C19 PMs."}],
    "apixaban":        [{"gene": "CYP3A4",  "note": "Modest role; major drug interaction concerns are CYP3A4 inhibitors."}],
    "rivaroxaban":     [{"gene": "CYP3A4",  "note": "Modest role."}],
    "dabigatran":      [{"gene": "—",        "note": "Not heavily CYP-dependent; mostly renal."}],
    "edoxaban":        [{"gene": "—",        "note": "Limited CYP involvement."}],

    # PPIs
    "omeprazole":      [{"gene": "CYP2C19", "note": "Primary metaboliser; *17 RM/UM may underrespond."}],
    "esomeprazole":    [{"gene": "CYP2C19", "note": "Primary metaboliser."}],
    "lansoprazole":    [{"gene": "CYP2C19", "note": "Primary metaboliser."}],
    "pantoprazole":    [{"gene": "CYP2C19", "note": "Primary metaboliser."}],
    "rabeprazole":     [{"gene": "CYP2C19", "note": "Modest role; mostly non-enzymatic."}],

    # Pain / opioids
    "codeine":         [{"gene": "CYP2D6", "note": "Pro-drug — PMs no analgesia; UMs RESPIRATORY DEPRESSION risk."}],
    "tramadol":        [{"gene": "CYP2D6", "note": "Pro-drug — PMs less effective; UMs higher AE risk."}],
    "hydrocodone":     [{"gene": "CYP2D6", "note": "Partial pro-drug metabolism."}],
    "oxycodone":       [{"gene": "CYP2D6", "note": "Conversion to oxymorphone (more potent)."}],
    "acetaminophen":   [{"gene": "—",        "note": "Not strongly PGx-dependent."},
                        {"gene": "CYP2E1",  "note": "Induced by chronic alcohol — hepatotoxicity risk."}],
    "ibuprofen":       [{"gene": "CYP2C9",  "note": "Metaboliser; PMs higher exposure."}],
    "naproxen":        [{"gene": "CYP2C9",  "note": "Modest role."}],
    "celecoxib":       [{"gene": "CYP2C9",  "note": "Primary metaboliser; PMs need dose reduction."}],

    # Antipsychotics
    "aripiprazole":    [{"gene": "CYP2D6", "note": "Major metaboliser — PMs reduce dose 50%."},
                        {"gene": "CYP3A4",  "note": "Secondary."}],
    "olanzapine":      [{"gene": "CYP1A2",  "note": "Primary metaboliser; smoking induces CYP1A2."},
                        {"gene": "CYP2D6",  "note": "Minor role."}],
    "quetiapine":      [{"gene": "CYP3A4",  "note": "Primary metaboliser."}],
    "risperidone":     [{"gene": "CYP2D6", "note": "Primary; PMs higher exposure to risperidone vs 9-OH."}],
    "ziprasidone":     [{"gene": "CYP3A4",  "note": "Primary."}],
    "clozapine":       [{"gene": "CYP1A2",  "note": "Primary; smoking strongly induces — dose adjustments at quit."},
                        {"gene": "HLA-DQB1", "note": "Agranulocytosis association (rare)."}],
    "haloperidol":     [{"gene": "CYP2D6", "note": "Major metaboliser."}],

    # Misc
    "carbamazepine":   [{"gene": "HLA-B*15:02", "note": "Severe SJS risk in carriers (Asian ancestry)."},
                        {"gene": "CYP3A4",  "note": "Auto-induces own metabolism."}],
    "phenytoin":       [{"gene": "CYP2C9",  "note": "Primary metaboliser; *3 PMs need dose reduction."},
                        {"gene": "HLA-B*15:02", "note": "SJS risk."}],
    "lamotrigine":     [{"gene": "HLA-B*15:02", "note": "SJS risk."}],
    "valproate":       [{"gene": "—",        "note": "Not strongly PGx-dependent."}],
    "azathioprine":    [{"gene": "TPMT",    "note": "PMs need DRAMATIC dose reduction — life-threatening toxicity."},
                        {"gene": "NUDT15",  "note": "Similar — particularly in East Asian patients."}],
    "6-mercaptopurine":[{"gene": "TPMT",    "note": "Same as azathioprine."},
                        {"gene": "NUDT15",  "note": "Same — critical in East Asian patients."}],
    "tamoxifen":       [{"gene": "CYP2D6", "note": "Activates to endoxifen; PMs may have reduced efficacy."}],
    "ondansetron":     [{"gene": "CYP2D6", "note": "UMs may have reduced antiemetic efficacy."}],
    "metoprolol":      [{"gene": "CYP2D6", "note": "Primary metaboliser; PMs higher exposure."}],
    "carvedilol":      [{"gene": "CYP2D6", "note": "Primary metaboliser."}],
    "tacrolimus":      [{"gene": "CYP3A5", "note": "Expressors metabolise rapidly — need ~2× the dose."}],
    "abacavir":        [{"gene": "HLA-B*57:01", "note": "CONTRAINDICATED in positive carriers — severe hypersensitivity."}],
}


def normalize_drug(name: str) -> str:
    """Strip whitespace, lowercase, brand → generic. Returns generic or input."""
    s = (name or "").strip().lower()
    return BRAND_TO_GENERIC.get(s, s)


def lookup_medication(drug_input: str, pgx_result: Dict) -> Dict:
    """For a single drug, find applicable genes + the user's phenotype +
    the CPIC recommendation if present."""
    generic = normalize_drug(drug_input)
    genes_info = DRUG_GENES.get(generic)

    if not genes_info:
        return {
            "input": drug_input,
            "generic": generic,
            "status": "unknown_drug",
            "message": (
                f"'{drug_input}' is not in the PGx lookup table. This does not "
                "mean it has no genetic implications — just that this report "
                "does not currently catalogue it."
            ),
            "findings": [],
        }

    findings: List[Dict] = []
    per_gene_results = (pgx_result or {}).get("per_gene", {})
    for entry in genes_info:
        gene = entry["gene"]
        if gene in ("—", "HLA-B*15:02", "CYP1A2", "CYP2B6", "ANKK1/DRD2", "COMT"):
            # Genes we don't currently phenotype via pgx.py — surface the note
            # but not a phenotype call.
            findings.append({
                "gene": gene,
                "pathway_note": entry["note"],
                "phenotype": None,
                "phenotype_class": "pheno-na",
                "activity_score": None,
                "recommendation": None,
                "covered_by_pgx": False,
            })
            continue
        gene_result = per_gene_results.get(gene)
        if not gene_result:
            findings.append({
                "gene": gene,
                "pathway_note": entry["note"],
                "phenotype": "Not phenotyped on this chip",
                "phenotype_class": "pheno-na",
                "activity_score": None,
                "recommendation": None,
                "covered_by_pgx": False,
            })
            continue
        # Find drug-specific recommendation if present
        drug_rec = None
        for d in gene_result.get("drug_recs", []):
            if any(generic.lower() in d.get("drug", "").lower() for _ in [0]) and \
               (generic.lower() in d.get("drug", "").lower() or
                drug_input.lower() in d.get("drug", "").lower()):
                phen_code = gene_result.get("phenotype_code")
                if phen_code and phen_code in d:
                    drug_rec = d[phen_code]
                    break
        findings.append({
            "gene": gene,
            "pathway_note": entry["note"],
            "phenotype": gene_result.get("phenotype"),
            "phenotype_class": gene_result.get("phenotype_class"),
            "activity_score": gene_result.get("activity_score"),
            "phenotype_code": gene_result.get("phenotype_code"),
            "guideline": gene_result.get("cpic_guideline"),
            "recommendation": drug_rec,
            "covered_by_pgx": True,
        })

    return {
        "input": drug_input,
        "generic": generic,
        "status": "ok",
        "findings": findings,
    }


def analyze_medications(drug_list: List[str], pgx_result: Dict) -> Dict:
    """Run lookup_medication over a user-provided list."""
    reviews = [lookup_medication(d, pgx_result) for d in drug_list]
    return {
        "n_input": len(drug_list),
        "n_known": sum(1 for r in reviews if r["status"] == "ok"),
        "n_unknown": sum(1 for r in reviews if r["status"] == "unknown_drug"),
        "reviews": reviews,
    }
