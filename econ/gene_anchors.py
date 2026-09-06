"""One gene → economics mapping, and one stated contract per field.

WHY THIS MODULE EXISTS
----------------------
Two tables used to answer "what is a pathogenic variant in this gene worth":
``ACMG_GENE_ECONOMICS`` in ``health_economics`` (17 genes) and ``_gene_to_econ``
in ``value_of_information`` (8 genes). They disagreed on every gene they shared,
by a consistent factor of roughly two.

The disagreement was not two opinions about one number. ``qaly_gain`` in the
ACMG table had **no documented semantics** — the field-documentation block in
``health_economics`` describes ``PGX_ECONOMICS``, not this one — so each consumer
invented its own contract and each was internally consistent:

* the PGx consumer multiplied it through ``p_rx x p_adr x rrr``;
* the ACMG consumer passed it straight to ``add()`` with nothing applied;
* the value-of-information path multiplied its own copy by ``p_event x rrr``.

Each table had been tuned to its own consumer, so the stated values differed by
about the ``rrr`` that one consumer applied and the other did not. Unifying the
tables without unifying the consumers would have produced the mirror-image
error.

A field name is not a contract. Every field below states what it means and who
applies it, and ``test_every_gene_anchor_field_has_a_stated_contract`` fails if a
new one appears without one.

FIELD CONTRACTS
---------------
``coi_key``
    The condition anchor. Selects BOTH the cost of illness and the QALY
    decrement from the parameter registry via ``engine.COI_KEY_TO_PARAM``, so
    the decrement is registry-backed and tiered rather than a module literal.
    ``None`` means the gene is deliberately not monetized — see NOT_VALUED_GENES.

``penetrance``
    Probability that a carrier of a pathogenic variant in this gene develops
    the condition, from clinically ascertained families. NOT a population
    frequency. The consumer is expected to correct it for ascertainment before
    use (``analyze_penetrance_posterior``).

``rrr``
    Relative risk reduction achieved by the surveillance or prophylaxis the
    finding indicates. Applied by the consumer, never pre-multiplied here.

``cost`` / ``outcome_value`` / ``clinical_benefit`` / ``src``
    Curated display economics for the report: what acting costs, the
    order-of-magnitude value of the averted outcome, the plain-language action,
    and the citation. ``outcome_value`` is a lifetime figure for one affected
    person; the consumer applies penetrance and effectiveness.

NOTE ON WHAT IS ABSENT
----------------------
There is deliberately no ``qaly_gain`` here. The QALY decrement comes from the
registry through ``coi_key``. A per-gene QALY literal is what allowed the two
tables to drift apart in the first place, and it never entered the provenance
count or the tornado.
"""
from __future__ import annotations

# Genes with a defensible registry-backed condition anchor.
#
# `penetrance` and `rrr` are carried over from the value-of-information table,
# which is the one whose consumer applied them correctly. `cost`,
# `outcome_value`, `clinical_benefit` and `src` are carried over from the ACMG
# table, which is the one that held curated display economics. Neither table
# "won": each supplied the quantity it actually had a contract for.
GENE_ANCHORS: dict[str, dict] = {
    "BRCA1": {
        "coi_key": "BreastOvarian", "penetrance": 0.60, "rrr": 0.45,
        "cost": 2_000, "outcome_value": 150_000,
        "clinical_benefit": "Enhanced screening + risk-reducing surgery option",
        "src": "Manchanda et al. (2015) J Clin Oncol",
        "finding": "BRCA1 pathogenic variant",
    },
    "BRCA2": {
        "coi_key": "BreastOvarian", "penetrance": 0.55, "rrr": 0.45,
        "cost": 2_000, "outcome_value": 140_000,
        "clinical_benefit": "Enhanced screening + risk-reducing surgery",
        "src": "Manchanda et al. (2015) J Clin Oncol",
        "finding": "BRCA2 pathogenic variant",
    },
    "PALB2": {
        "coi_key": "BreastOvarian", "penetrance": 0.40, "rrr": 0.40,
        "cost": 2_000, "outcome_value": 120_000,
        "clinical_benefit": "Enhanced screening; surgical discussion",
        "src": "Antoniou et al. (2014) NEJM — PALB2 breast-cancer risk",
        "finding": "PALB2 pathogenic variant",
    },
    "MLH1": {
        "coi_key": "Colorectal", "penetrance": 0.50, "rrr": 0.50,
        "cost": 1_500, "outcome_value": 100_000,
        "clinical_benefit": "Annual colonoscopy from age 25 + aspirin",
        "src": "Ladabaum et al. (2011) Ann Intern Med",
        "finding": "MLH1 Lynch syndrome variant",
    },
    "MSH2": {
        "coi_key": "Colorectal", "penetrance": 0.45, "rrr": 0.50,
        "cost": 1_500, "outcome_value": 100_000,
        "clinical_benefit": "Annual colonoscopy + gynecologic surveillance",
        "src": "Ladabaum et al. (2011) Ann Intern Med",
        "finding": "MSH2 Lynch syndrome variant",
    },
    "MSH6": {
        "coi_key": "Colorectal", "penetrance": 0.30, "rrr": 0.45,
        "cost": 1_500, "outcome_value": 80_000,
        "clinical_benefit": "Enhanced colonoscopy + endometrial screening",
        "src": "Ladabaum et al. (2011) Ann Intern Med",
        "finding": "MSH6 Lynch syndrome variant",
    },
    # LOOKUP, not a judgment: PMS2 is a Lynch mismatch-repair gene managed by
    # the same colonoscopy surveillance as MLH1/MSH2/MSH6, so it takes the
    # anchor those three already use. Its penetrance is set below the other MMR
    # genes on ten Broeke et al. (2015) J Clin Oncol, which found PMS2 carrier
    # risk substantially lower than MLH1/MSH2 — that is why its curated
    # outcome_value was already the smallest of the four. `src` stays the exact
    # Ladabaum string the provenance resolution map keys on; the penetrance
    # citation is recorded here rather than appended to it, because that map
    # matches source strings literally.
    "PMS2": {
        "coi_key": "Colorectal", "penetrance": 0.20, "rrr": 0.45,
        "cost": 1_500, "outcome_value": 60_000,
        "clinical_benefit": "Colonoscopy surveillance program",
        "src": "Ladabaum et al. (2011) Ann Intern Med",
        "finding": "PMS2 Lynch syndrome variant",
    },
    "LDLR": {
        "coi_key": "CAD", "penetrance": 0.50, "rrr": 0.50,
        "cost": 500, "outcome_value": 200_000,
        "clinical_benefit": "High-intensity statin + cascade screening",
        "src": "Nherera et al. (2011) Heart — FH CEA",
        "finding": "LDLR familial hypercholesterolemia",
    },
    "APOB": {
        "coi_key": "CAD", "penetrance": 0.40, "rrr": 0.45,
        "cost": 500, "outcome_value": 180_000,
        "clinical_benefit": "High-intensity statin therapy",
        "src": "Nherera et al. (2011) Heart — FH CEA",
        "finding": "APOB familial hypercholesterolemia",
    },
    # LOOKUP, not a judgment: PCSK9 gain-of-function causes the same
    # autosomal-dominant familial hypercholesterolaemia phenotype as LDLR and
    # APOB, treated on the same lipid-lowering pathway, so it takes the anchor
    # both already use.
    "PCSK9": {
        "coi_key": "CAD", "penetrance": 0.40, "rrr": 0.50,
        "cost": 3_000, "outcome_value": 200_000,
        "clinical_benefit": "PCSK9 inhibitor + cascade screening",
        "src": "Kazi et al. (2017) JAMA Cardiol",
        "finding": "PCSK9 familial hypercholesterolemia",
    },
}


# Genes detected and reported clinically, deliberately NOT priced.
#
# Each of these has a real curated action and a real clinical meaning. What none
# of them has is a registry-backed cost-of-illness anchor for its condition, and
# the registry holds no entry for endocrine neoplasia, cardiac arrhythmia,
# hypertrophic cardiomyopathy, pediatric retinoblastoma, or a multi-site
# cancer syndrome.
#
# Assigning them an approximate anchor would have been easy and would have
# raised the headline number: RET carries the largest QALY figure in the old
# table. That is exactly why it is not done here. Substituting `coi_mace` for
# hypertrophic cardiomyopathy because both are cardiac, or picking the larger
# of two single-site cancer anchors for a multi-site syndrome, is the
# assumption-laundering the parameter registry exists to prevent — and doing it
# on the genes with the biggest numbers would put the model's least defensible
# figures in its most prominent position.
#
# The finding is still reported. Only the dollar figure is withheld.
NOT_VALUED_GENES: dict[str, dict] = {
    "RET": {
        "finding": "RET MEN2 variant",
        "clinical_benefit": "Prophylactic thyroidectomy + calcitonin monitoring",
        "reason": "No registry anchor exists for endocrine neoplasia. Pricing "
                  "this would mean inventing a cost of illness and a baseline "
                  "event probability for the gene that previously carried the "
                  "largest QALY figure in the model.",
        "src": "Wells et al. (2015) Thyroid; Brandi (2001) JCEM",
    },
    "SCN5A": {
        "finding": "SCN5A channelopathy",
        "clinical_benefit": "Cardiac monitoring + beta-blocker/ICD assessment",
        "reason": "No registry anchor exists for cardiac arrhythmia. `coi_mace` "
                  "is atherosclerotic and is a different disease.",
        "src": "Kaufman et al. (2014) Circ Cardiovasc Genet",
    },
    "KCNQ1": {
        "finding": "KCNQ1 long-QT syndrome",
        "clinical_benefit": "Beta-blocker + activity restriction",
        "reason": "No registry anchor exists for cardiac arrhythmia.",
        "src": "Kaufman et al. (2014) Circ Cardiovasc Genet",
    },
    "KCNH2": {
        "finding": "KCNH2 long-QT syndrome",
        "clinical_benefit": "Beta-blocker + QT-prolonging drug avoidance",
        "reason": "No registry anchor exists for cardiac arrhythmia.",
        "src": "Kaufman et al. (2014) Circ Cardiovasc Genet",
    },
    "MYH7": {
        "finding": "MYH7 hypertrophic cardiomyopathy",
        "clinical_benefit": "Echo surveillance + exercise restriction",
        "reason": "Hypertrophic cardiomyopathy is not atherosclerotic MACE. "
                  "The same order of magnitude is not the same disease.",
        "src": "Maron et al. (2014) Circulation — HCM guidelines",
    },
    "MYBPC3": {
        "finding": "MYBPC3 hypertrophic cardiomyopathy",
        "clinical_benefit": "Echo surveillance + cascade screening",
        "reason": "As MYH7 — no cardiomyopathy anchor exists.",
        "src": "Maron et al. (2014) Circulation",
    },
    "RB1": {
        "finding": "RB1 retinoblastoma variant",
        "clinical_benefit": "Pediatric eye exams + family screening",
        "reason": "`qaly_loss_cancer` is an adult, stage-weighted decrement. "
                  "Applying it to a pediatric ocular cancer would borrow a "
                  "figure derived from a different population.",
        "src": "Soliman et al. (2016) J AAPOS",
    },
    "TP53": {
        "finding": "TP53 Li-Fraumeni variant",
        "clinical_benefit": "Annual whole-body MRI (Toronto protocol)",
        "reason": "Li-Fraumeni is a multi-site cancer syndrome, so neither the "
                  "colorectal nor the breast/ovarian cost anchor describes it. "
                  "Choosing the larger of the two to avoid leaving a gap is the "
                  "failure this withholding exists to prevent.",
        "src": "Villani et al. (2016) Lancet Oncol",
    },
}


def anchor_for(gene: str) -> dict | None:
    """The economics for `gene`, or None if it has no defensible anchor.

    None is a real answer, not a missing one: callers must report the finding
    and withhold the value rather than substitute a generic bucket. See
    NOT_VALUED_GENES for why each withheld gene is withheld.
    """
    return GENE_ANCHORS.get((gene or "").upper())


def not_valued_reason(gene: str) -> str:
    """Why `gene` carries no dollar figure, or "" if it does carry one."""
    return (NOT_VALUED_GENES.get((gene or "").upper()) or {}).get("reason", "")


def known_gene(gene: str) -> bool:
    g = (gene or "").upper()
    return g in GENE_ANCHORS or g in NOT_VALUED_GENES
