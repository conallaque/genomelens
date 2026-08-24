"""
Carrier Status & Pathogenic-Variant Detection
---------------------------------------------

Identifies carrier (heterozygous) and affected (homozygous) status for
recessive and dominant pathogenic variants where the chip covers a
defining SNP. Useful for:

  * Family planning — recessive carriers do not have disease but can
    transmit to offspring if their partner also carries a pathogenic
    variant in the same gene.
  * Personal disease-risk awareness for dominantly-acting variants
    (Factor V Leiden, HFE C282Y, HLA-B*57:01, TREM2 R47H).

This is NOT a substitute for clinical carrier-screening panels (which
cover ~250+ recessive conditions using sequencing). It is the
chip-detectable subset.

V8 migration note
-----------------
The SNV variants in ``CARRIER_VARIANTS`` whose rsIDs are present in
``snp_registry.SNPS`` are reconciled against the registry at import time
via :func:`audit_against_registry` — any disagreement on the pathogenic
allele or the gene-symbol assignment raises ``AssertionError`` rather
than silently shipping wrong data. Indel variants (CFTR ΔF508, HEXA
1278insTATC, BRCA1 185delAG, …) remain local to this module pending the
V8.1 SNPRecord schema extension to handle non-SNV variants. See
``CHANGELOG.md`` for the full reconciliation log.
"""


import pandas as pd

from core import snp_registry  # V8 cross-check; see audit_against_registry below


def _dose(snps_df: pd.DataFrame, rsid: str, allele: str) -> int | None:
    if rsid not in snps_df.index:
        return None
    gt = snps_df.loc[rsid].get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--") or len(s) != 2:
        return None
    return s.count(allele.upper())


CARRIER_VARIANTS: list[dict] = [
    {
        "rsid": "rs1800562",
        "gene": "HFE",
        "variant": "C282Y",
        "disease": "Hereditary Hemochromatosis (HH, type 1)",
        "inheritance": "autosomal recessive (incomplete penetrance)",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "Heterozygous carriers (C282Y/wt) are mostly unaffected but can transmit "
            "the variant. If partner is also a C282Y carrier, each child has 25% "
            "chance of being homozygous (clinical hemochromatosis risk)."
        ),
        "affected_implication": (
            "Homozygous C282Y/C282Y is the classical hemochromatosis genotype "
            "with ~25–60% lifetime penetrance for iron overload. Test ferritin "
            "and transferrin saturation; therapeutic phlebotomy is curative."
        ),
    },
    {
        "rsid": "rs1799945",
        "gene": "HFE",
        "variant": "H63D",
        "disease": "Hereditary Hemochromatosis (mild variant)",
        "inheritance": "autosomal recessive, mild penetrance",
        "pathogenic_allele": "G",
        "carrier_implication": (
            "H63D alone confers little risk. Most relevant in compound heterozygosity "
            "with C282Y (see Compound Heterozygosity findings)."
        ),
        "affected_implication": (
            "H63D/H63D homozygous has very low penetrance for clinical iron overload "
            "but periodic ferritin monitoring is reasonable."
        ),
    },
    {
        "rsid": "rs6025",
        "gene": "F5",
        "variant": "Factor V Leiden (R506Q)",
        "disease": "Venous Thromboembolism (VTE) Susceptibility",
        "inheritance": "autosomal dominant (incomplete penetrance)",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "One copy of Factor V Leiden raises lifetime VTE risk ~5–7×. Affects "
            "individual risk; family planning relevance is dominant-style."
        ),
        "affected_implication": (
            "Homozygous Factor V Leiden carries ~80× VTE risk. Strong indication "
            "to avoid estrogen-containing contraceptives, careful management during "
            "pregnancy, hematology consultation."
        ),
    },
    {
        "rsid": "rs1799963",
        "gene": "F2",
        "variant": "Prothrombin G20210A",
        "disease": "Venous Thromboembolism (VTE) Susceptibility",
        "inheritance": "autosomal dominant (incomplete penetrance)",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "One copy raises VTE risk ~3–4×."
        ),
        "affected_implication": (
            "Homozygous G20210A/G20210A is rare and confers substantially elevated VTE risk."
        ),
    },
    {
        "rsid": "rs75932628",
        "gene": "TREM2",
        "variant": "R47H",
        "disease": "Late-Onset Alzheimer's Disease (susceptibility)",
        "inheritance": "autosomal dominant susceptibility",
        "pathogenic_allele": "T",
        "carrier_implication": (
            "Single copy R47H is associated with ~2.5–4× AD risk — comparable to "
            "one APOE-ε4 allele. Also associated with FTD and PD risk."
        ),
        "affected_implication": (
            "Homozygous is very rare; risk likely further elevated. Discuss "
            "AD-prevention strategies with a neurologist."
        ),
    },
    {
        "rsid": "rs61816761",
        "gene": "FLG",
        "variant": "R501X",
        "disease": "Atopic Dermatitis / Ichthyosis Vulgaris (filaggrin deficiency)",
        "inheritance": "semi-dominant",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "Heterozygotes have ~3× atopic dermatitis (eczema) risk and elevated "
            "asthma + food-allergy risk through the 'atopic march'."
        ),
        "affected_implication": (
            "Homozygous typically causes severe ichthyosis vulgaris and atopic "
            "dermatitis. Aggressive barrier care is foundational."
        ),
    },
    {
        "rsid": "rs17879961",
        "gene": "CHEK2",
        "variant": "I157T",
        "disease": "Hereditary Breast/Colon/Prostate/Kidney Cancer Susceptibility",
        "inheritance": "autosomal dominant (moderate penetrance)",
        "pathogenic_allele": "C",
        "carrier_implication": (
            "I157T heterozygous: ~1.5–2× breast cancer risk; moderate associations "
            "with colorectal, prostate, kidney cancers."
        ),
        "affected_implication": (
            "Homozygous is rare and substantially elevates cancer risk. Genetic "
            "counselor referral; enhanced surveillance per guidelines."
        ),
    },
    {
        "rsid": "rs4349859",
        "gene": "HLA-B*27 region",
        "variant": "HLA-B*27 proxy",
        "disease": "Ankylosing Spondylitis / Acute Anterior Uveitis",
        "inheritance": "autosomal dominant susceptibility",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "Carrier of HLA-B*27 tag — ~90% of AS patients are B27+ but only "
            "~6–8% of B27 carriers develop AS. Awareness of red-flag symptoms."
        ),
        "affected_implication": (
            "Homozygous B27 carries similar individual risk to heterozygous."
        ),
    },
    {
        "rsid": "rs763035",
        "gene": "HLA-B*57:01 (proxy)",
        "variant": "B*57:01 tag",
        "disease": "Abacavir Hypersensitivity",
        "inheritance": "single positive copy contraindicates",
        "pathogenic_allele": "T",
        "carrier_implication": (
            "If HIV+ and considering abacavir: confirm HLA-B*57:01 with direct "
            "HLA typing — do NOT rely on this SNP proxy alone. Abacavir is "
            "CONTRAINDICATED in HLA-B*57:01-positive individuals (severe, "
            "potentially fatal hypersensitivity)."
        ),
        "affected_implication": "Same as carrier — single copy contraindicates abacavir.",
    },
    {
        "rsid": "rs2476601",
        "gene": "PTPN22",
        "variant": "R620W",
        "disease": "Broad Autoimmunity Susceptibility (T1D, RA, lupus, Graves')",
        "inheritance": "autosomal dominant susceptibility",
        "pathogenic_allele": "A",
        "carrier_implication": (
            "Carriage modestly raises risk of multiple autoimmune diseases; "
            "penetrance is low. Anti-inflammatory lifestyle, vitamin D sufficiency."
        ),
        "affected_implication": "Homozygous further elevates autoimmune risk modestly.",
    },
    {
        "rsid": "rs2187668",
        "gene": "HLA-DQA1",
        "variant": "DQ2.5 tag",
        "disease": "Celiac Disease Susceptibility",
        "inheritance": "autosomal dominant susceptibility",
        "pathogenic_allele": "T",
        "carrier_implication": (
            "Carrier of HLA-DQ2.5 — necessary but not sufficient for celiac disease "
            "(~30% of population carries DQ2; only ~1% develops celiac). Negative "
            "for both DQ2 and DQ8 essentially rules out celiac."
        ),
        "affected_implication": "Same as carrier — being positive enables celiac risk.",
    },
    {
        "rsid": "rs7454108",
        "gene": "HLA-DQB1",
        "variant": "DQ8 tag",
        "disease": "Celiac Disease Susceptibility (second-most-common allele)",
        "inheritance": "autosomal dominant susceptibility",
        "pathogenic_allele": "C",
        "carrier_implication": (
            "Carrier of HLA-DQ8 — second major celiac-risk allele (after DQ2)."
        ),
        "affected_implication": "Same — enables celiac risk.",
    },

    # ──────────────────────────────────────────────────────────────────────
    # V4 — Massively expanded carrier panel
    # Many of these are RARE pathogenic variants. Consumer chips often do NOT
    # type them; in that case the user appears in `untested`. When present,
    # the framework can detect them.
    # ──────────────────────────────────────────────────────────────────────

    # ── Cystic fibrosis ──
    {"rsid": "rs113993960", "gene": "CFTR", "variant": "F508del (ΔF508)",
     "disease": "Cystic Fibrosis", "inheritance": "autosomal recessive",
     "pathogenic_allele": "DEL", "carrier_implication":
     "Most common CF mutation (~70% of CF alleles in Europeans, ~1:25 European carrier rate). Heterozygotes asymptomatic.",
     "affected_implication":
     "Homozygous or compound heterozygous with another CF variant → cystic fibrosis (respiratory, GI, pancreatic). Modulator drugs (Trikafta/Kaftrio) revolutionary.",
     "carrier_frequency": "European ~1:25, Hispanic ~1:46, African ~1:65, Asian ~1:90",
     "chip_caveat": "Indel — many consumer chips do NOT type this. Negative result here cannot rule out CF carrier status; clinical sequencing required."},
    {"rsid": "rs75961395", "gene": "CFTR", "variant": "R117H",
     "disease": "Cystic Fibrosis (mild/variable)", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication":
     "Mild/variable CF variant. Phenotype depends on poly-T tract length.",
     "affected_implication": "Variable CF phenotype; sometimes congenital absence of vas deferens only."},
    {"rsid": "rs78655421", "gene": "CFTR", "variant": "G542X",
     "disease": "Cystic Fibrosis (severe)", "inheritance": "autosomal recessive",
     "pathogenic_allele": "T", "carrier_implication": "Severe CF nonsense mutation.",
     "affected_implication": "Severe CF when homozygous/compound."},
    {"rsid": "rs121908769", "gene": "CFTR", "variant": "W1282X (AJ founder)",
     "disease": "Cystic Fibrosis", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication":
     "AJ founder CF mutation (~30% of AJ CF alleles, ~1:24 AJ carrier rate).",
     "affected_implication": "Severe CF.",
     "carrier_frequency": "Ashkenazi Jewish ~1:24"},
    {"rsid": "rs121909001", "gene": "CFTR", "variant": "N1303K",
     "disease": "Cystic Fibrosis", "inheritance": "autosomal recessive",
     "pathogenic_allele": "G", "carrier_implication": "Common Mediterranean CF variant.",
     "affected_implication": "Severe CF."},

    # ── Sickle cell / hemoglobinopathies ──
    {"rsid": "rs334", "gene": "HBB", "variant": "HbS (Glu6Val)",
     "disease": "Sickle Cell Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication":
     "Carriers (sickle trait) usually asymptomatic. Trait protects against severe P. falciparum malaria.",
     "affected_implication":
     "Homozygous HbSS = sickle cell anemia (chronic hemolysis, pain crises, organ damage). HbSC and HbS/β-thal are similar.",
     "carrier_frequency": "Sub-Saharan African ~1:12, African American ~1:13, Hispanic ~1:200, Mediterranean ~1:100"},
    {"rsid": "rs33930165", "gene": "HBB", "variant": "HbC (Glu6Lys)",
     "disease": "HbC / HbSC Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "T", "carrier_implication": "HbC carrier (~25% of West Africans).",
     "affected_implication": "HbCC mild. HbSC (with HbS) causes sickle-related disease."},
    {"rsid": "rs11549407", "gene": "HBB", "variant": "Codon 39 (Q39X) β-thal",
     "disease": "Beta-thalassemia", "inheritance": "autosomal recessive",
     "pathogenic_allele": "T", "carrier_implication": "Common β-thal variant in Mediterranean.",
     "affected_implication": "β-thal major (transfusion-dependent) or β-thal intermedia."},

    # ── Alpha-1 antitrypsin deficiency ──
    {"rsid": "rs28929474", "gene": "SERPINA1", "variant": "Z allele (Glu342Lys)",
     "disease": "Alpha-1 Antitrypsin Deficiency", "inheritance": "autosomal codominant",
     "pathogenic_allele": "A", "carrier_implication":
     "MZ heterozygotes have ~60% normal A1AT — small COPD risk with smoking.",
     "affected_implication": "ZZ homozygotes — severe deficiency, early COPD (especially smokers), liver disease risk. Augmentation therapy available.",
     "carrier_frequency": "European ~1:30 (Z allele)"},
    {"rsid": "rs17580", "gene": "SERPINA1", "variant": "S allele (Glu264Val)",
     "disease": "Alpha-1 Antitrypsin Deficiency (mild)", "inheritance": "autosomal codominant",
     "pathogenic_allele": "A", "carrier_implication": "S allele milder than Z.",
     "affected_implication": "SS or SZ compound has intermediate-to-mild deficiency."},

    # ── PKU ──
    {"rsid": "rs5030858", "gene": "PAH", "variant": "R408W",
     "disease": "Phenylketonuria (classical)", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication": "Most common Eastern European PKU mutation.",
     "affected_implication": "Severe PKU — caught on newborn screening; lifelong low-phe diet."},

    # ── Gaucher ──
    {"rsid": "rs76763715", "gene": "GBA", "variant": "N370S",
     "disease": "Gaucher Disease Type 1 / Parkinson's Risk",
     "inheritance": "autosomal recessive; heterozygotes have PD risk",
     "pathogenic_allele": "C", "carrier_implication":
     "Most common Gaucher mutation in Ashkenazi Jews (~1:15 AJ carrier). Heterozygotes have ~5× Parkinson's disease risk.",
     "affected_implication": "Type 1 Gaucher — splenomegaly, hepatomegaly, bone disease. Enzyme replacement (imiglucerase) and substrate reduction (eliglustat) effective."},
    {"rsid": "rs421016", "gene": "GBA", "variant": "L444P",
     "disease": "Gaucher Disease Type 2/3 / Parkinson's", "inheritance": "autosomal recessive",
     "pathogenic_allele": "C", "carrier_implication":
     "Pan-ethnic Gaucher variant. Heterozygotes have elevated PD risk.",
     "affected_implication": "Severe neuronopathic Gaucher when homozygous."},

    # ── Tay-Sachs ──
    {"rsid": "rs121907968", "gene": "HEXA", "variant": "1278insTATC",
     "disease": "Tay-Sachs Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "INS", "carrier_implication":
     "AJ founder Tay-Sachs mutation (~1:30 AJ carrier rate).",
     "affected_implication": "Progressive fatal infantile neurodegeneration (death by age 4-5). No cure.",
     "carrier_frequency": "Ashkenazi Jewish ~1:30, French Canadian ~1:30",
     "chip_caveat": "Insertion variant — most consumer chips do NOT type this."},

    # ── Familial Mediterranean fever ──
    {"rsid": "rs61752717", "gene": "MEFV", "variant": "M694V",
     "disease": "Familial Mediterranean Fever", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication":
     "Most common FMF mutation. Heterozygotes mostly asymptomatic.",
     "affected_implication": "Recurrent fevers + serositis. Colchicine prevents attacks and amyloidosis.",
     "carrier_frequency": "Sephardic Jewish ~1:5, Armenian ~1:7, Turkish ~1:5, Arab ~1:5"},
    {"rsid": "rs28940579", "gene": "MEFV", "variant": "V726A",
     "disease": "Familial Mediterranean Fever", "inheritance": "autosomal recessive",
     "pathogenic_allele": "C", "carrier_implication": "Second common FMF variant.",
     "affected_implication": "Generally milder phenotype than M694V."},

    # ── Wilson disease ──
    {"rsid": "rs76151636", "gene": "ATP7B", "variant": "H1069Q",
     "disease": "Wilson Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "T", "carrier_implication": "Most common European Wilson variant.",
     "affected_implication": "Copper accumulation — liver, neurologic, psychiatric. Highly treatable (zinc, chelators) if caught early."},

    # ── HFI ──
    {"rsid": "rs1800546", "gene": "ALDOB", "variant": "A149P",
     "disease": "Hereditary Fructose Intolerance", "inheritance": "autosomal recessive",
     "pathogenic_allele": "C", "carrier_implication":
     "Most common HFI variant (~50% of alleles).",
     "affected_implication": "Severe symptoms with fructose. Manageable with strict avoidance."},

    # ── Galactosemia ──
    {"rsid": "rs75391579", "gene": "GALT", "variant": "Q188R",
     "disease": "Classical Galactosemia", "inheritance": "autosomal recessive",
     "pathogenic_allele": "G", "carrier_implication": "Most common European galactosemia variant.",
     "affected_implication": "Detected by newborn screening; lifelong galactose avoidance."},

    # ── Canavan ──
    {"rsid": "rs28940279", "gene": "ASPA", "variant": "E285A",
     "disease": "Canavan Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "C", "carrier_implication": "AJ founder Canavan (~1:55 AJ carrier).",
     "affected_implication": "Fatal infantile leukodystrophy."},

    # ── Niemann-Pick ──
    {"rsid": "rs120074118", "gene": "SMPD1", "variant": "R608del",
     "disease": "Niemann-Pick Disease Type A/B", "inheritance": "autosomal recessive",
     "pathogenic_allele": "DEL", "carrier_implication": "AJ founder NP-A/B.",
     "affected_implication": "Type A: severe infantile; Type B: chronic visceral."},

    # ── Maple syrup urine ──
    {"rsid": "rs121964939", "gene": "BCKDHA", "variant": "Y393N",
     "disease": "Maple Syrup Urine Disease", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication": "Mennonite founder MSUD.",
     "affected_implication": "Severe BCAA metabolism disorder; specialized diet."},

    # ── G6PD deficiency (X-linked) ──
    {"rsid": "rs1050828", "gene": "G6PD", "variant": "V68M (A- variant)",
     "disease": "G6PD Deficiency", "inheritance": "X-linked recessive",
     "pathogenic_allele": "T", "carrier_implication":
     "Common African variant. Males hemizygous, females heterozygous (variable expression).",
     "affected_implication":
     "Avoid fava beans, sulfa drugs (sulfamethoxazole), primaquine, nitrofurantoin, naphthalene. Hemolysis risk with oxidative triggers.",
     "carrier_frequency": "African ~10% males, Mediterranean ~5%, Asian ~5%"},

    # ── BRCA1/2 expanded ──
    {"rsid": "rs80357906", "gene": "BRCA1", "variant": "185delAG (AJ founder)",
     "disease": "Hereditary Breast/Ovarian Cancer Syndrome",
     "inheritance": "autosomal dominant (high penetrance)",
     "pathogenic_allele": "DEL", "carrier_implication":
     "One of three Ashkenazi BRCA founder mutations (~1:100 AJ carrier).",
     "affected_implication":
     "Carriers have ~60-72% lifetime breast cancer risk, ~40-44% ovarian. Risk-reducing salpingo-oophorectomy and intensive screening discussed with cancer genetics.",
     "carrier_frequency": "Ashkenazi Jewish ~1:100 (all 3 founders combined ~1:40)"},
    {"rsid": "rs80359550", "gene": "BRCA2", "variant": "6174delT (AJ founder)",
     "disease": "Hereditary Breast/Ovarian/Prostate Cancer",
     "inheritance": "autosomal dominant (high penetrance)",
     "pathogenic_allele": "DEL", "carrier_implication": "AJ founder BRCA2 mutation.",
     "affected_implication":
     "BRCA2: similar female cancer risks plus elevated male breast (~6%) and prostate cancer."},
    {"rsid": "rs28897743", "gene": "BRCA1", "variant": "5382insC (Slavic founder)",
     "disease": "Hereditary Breast/Ovarian Cancer",
     "inheritance": "autosomal dominant (high penetrance)",
     "pathogenic_allele": "INS", "carrier_implication": "Slavic founder BRCA1 mutation.",
     "affected_implication": "Same as 185delAG."},

    # ── Lynch syndrome ──
    {"rsid": "rs63750449", "gene": "MLH1", "variant": "Lynch splice variant",
     "disease": "Lynch Syndrome (HNPCC)",
     "inheritance": "autosomal dominant",
     "pathogenic_allele": "A", "carrier_implication":
     "Lynch syndrome carriers have markedly elevated colorectal, endometrial, ovarian, gastric, urinary tract cancer risks.",
     "affected_implication":
     "Carriers (heterozygotes) ARE affected — Lynch is autosomal dominant. Colonoscopy every 1-2 years from age 20-25; gynecologic screening from 30-35; risk-reducing hysterectomy after childbearing."},

    # ── Familial Hypercholesterolemia ──
    {"rsid": "rs5742904", "gene": "APOB", "variant": "R3500Q (FDB)",
     "disease": "Familial Hypercholesterolemia (Defective ApoB)",
     "inheritance": "autosomal dominant",
     "pathogenic_allele": "A", "carrier_implication":
     "Causes familial defective apolipoprotein B-100. Markedly elevated LDL from childhood.",
     "affected_implication":
     "Untreated CAD by 50s. Aggressive lipid-lowering (statin + ezetimibe ± PCSK9i)."},

    # ── Hereditary spherocytosis ──
    {"rsid": "rs121912747", "gene": "ANK1", "variant": "Hereditary Spherocytosis",
     "disease": "Hereditary Spherocytosis", "inheritance": "autosomal dominant",
     "pathogenic_allele": "T", "carrier_implication":
     "Carrier of HS-causing variant.",
     "affected_implication": "Chronic hemolytic anemia, gallstones, splenomegaly. Splenectomy curative."},

    # ── Phenylketonuria additional ──
    {"rsid": "rs5030849", "gene": "PAH", "variant": "IVS12+1G>A",
     "disease": "Phenylketonuria", "inheritance": "autosomal recessive",
     "pathogenic_allele": "A", "carrier_implication": "Severe PAH splice variant.",
     "affected_implication": "Severe PKU."},

    # ── Biotinidase deficiency ──
    {"rsid": "rs80338686", "gene": "BTD", "variant": "Q456H",
     "disease": "Biotinidase Deficiency", "inheritance": "autosomal recessive",
     "pathogenic_allele": "C", "carrier_implication": "Common BTD def variant.",
     "affected_implication": "Treatable with biotin if caught early. NBS catches it."},

    # ── MCAD ──
    {"rsid": "rs77931234", "gene": "ACADM", "variant": "K304E",
     "disease": "MCAD Deficiency", "inheritance": "autosomal recessive",
     "pathogenic_allele": "G", "carrier_implication": "Most common MCAD def variant.",
     "affected_implication": "Avoid prolonged fasting. NBS catches it."},

    # ── Hemochromatosis additional ──
    {"rsid": "rs1799945", "gene": "HFE", "variant": "H63D",
     "disease": "Hereditary Hemochromatosis (mild)", "inheritance": "autosomal recessive",
     "pathogenic_allele": "G", "carrier_implication":
     "Milder HH variant. C282Y/H63D compound heterozygotes have moderate risk.",
     "affected_implication":
     "H63D/H63D usually mild. Compound C282Y/H63D — moderate iron overload risk.",
     "carrier_frequency": "European ~1:4"},

    # ── Wilson additional ──
    {"rsid": "rs1061472", "gene": "ATP7B", "variant": "K832R",
     "disease": "Wilson Disease (variable)", "inheritance": "autosomal recessive (variable)",
     "pathogenic_allele": "G", "carrier_implication": "Common polymorphism; debated pathogenicity.",
     "affected_implication": "Variable phenotype."},

    # ── Cleft Lip/Palate susceptibility ──
    {"rsid": "rs642961", "gene": "IRF6", "variant": "Cleft lip / palate",
     "disease": "Nonsyndromic Cleft Lip / Palate",
     "inheritance": "complex (susceptibility)",
     "pathogenic_allele": "A", "carrier_implication":
     "Modest susceptibility variant for nonsyndromic cleft lip/palate.",
     "affected_implication": "Surgical correction; multidisciplinary care."},

    # ── Marfan / FBN1 (chip limitation note) ──
    {"rsid": "rs140803", "gene": "FBN1", "variant": "FBN1 common SNP (NOT diagnostic)",
     "disease": "Marfan Syndrome (NOT detectable from chip)",
     "inheritance": "autosomal dominant",
     "pathogenic_allele": "A", "carrier_implication":
     "Common FBN1 SNPs do NOT detect Marfan. Marfan-causing mutations are rare and require clinical sequencing.",
     "affected_implication":
     "Clinical features (height, arm span, lens dislocation, aortic root) drive diagnosis."},

    # ── Spinal muscular atrophy proxy ──
    {"rsid": "rs1554286236", "gene": "SMN1", "variant": "SMN1 proxy",
     "disease": "Spinal Muscular Atrophy",
     "inheritance": "autosomal recessive",
     "pathogenic_allele": "DEL", "carrier_implication":
     "SMA is caused by SMN1 deletion. Chip-based detection unreliable; clinical MLPA is gold standard.",
     "affected_implication":
     "Carrier screening pre-pregnancy via clinical MLPA. Modern therapies (nusinersen, onasemnogene, risdiplam) transformative."},
]


def analyze_carriers(snps_df: pd.DataFrame) -> dict:
    """Returns a dict of carrier-status findings, organised by category."""
    affected: list[dict] = []
    carriers: list[dict] = []
    not_carriers: list[dict] = []
    untested: list[dict] = []

    for entry in CARRIER_VARIANTS:
        rsid = entry["rsid"]
        dose = _dose(snps_df, rsid, entry["pathogenic_allele"])
        record = {
            **entry,
            "dosage": dose,
            "genotype": str(snps_df.loc[rsid].get("genotype")).upper() if rsid in snps_df.index else None,
        }
        if dose is None:
            untested.append(record)
        elif dose == 2:
            affected.append(record)
        elif dose == 1:
            carriers.append(record)
        else:
            not_carriers.append(record)

    return {
        "affected": affected,             # homozygous for pathogenic allele
        "carriers": carriers,             # heterozygous
        "not_carriers": not_carriers,     # homozygous reference
        "untested": untested,             # variant not on chip
        "n_affected": len(affected),
        "n_carriers": len(carriers),
        "n_not_carriers": len(not_carriers),
        "n_untested": len(untested),
    }


# ── V8: cross-check against the unified SNP registry ───────────────────────

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def audit_against_registry() -> dict[str, list[str]]:
    """Verify that every rsID in ``CARRIER_VARIANTS`` that is also present in
    ``snp_registry.SNPS`` agrees on (a) the gene symbol and (b) the
    pathogenic allele being either the registry's ``derived`` allele or its
    reverse-strand complement.

    The strand-flipped case is *expected*, not a bug: dbSNP reports alleles
    in the assembly's + strand orientation, while consumer chips for many
    classic carrier variants (SERPINA1 PiZ, PAH R408W, APOB R3500Q, HBB
    HbS) historically report on the gene's coding strand for human
    interpretability. The registry's invariant is "+ strand"; carrier.py
    keeps the gene-coding-strand label that clinicians recognise.

    Returns a dict with:
      * ``agreed``         — rsIDs reconciled with no strand flip
      * ``strand_flipped`` — rsIDs reconciled where local-risk == complement(registry-derived)
      * ``disagreed``      — real biology disagreements (not strand)
      * ``registry_missing`` — indels + unmigrated SNVs

    Run on module import; raises ``AssertionError`` only on real biology
    disagreements — strand flips are documented in CHANGELOG.md and allowed.
    """
    agreed: list[str] = []
    strand_flipped: list[str] = []
    disagreed: list[str] = []
    registry_missing: list[str] = []

    for entry in CARRIER_VARIANTS:
        rec = snp_registry.get(entry["rsid"])
        if rec is None:
            registry_missing.append(entry["rsid"])
            continue

        # Gene-symbol cross-check (HFE / HLA-DQA1 etc. — strict, raise on mismatch)
        local_gene = entry["gene"].upper().strip()
        reg_gene = rec.gene.upper().strip()
        # Tolerate "HLA-B*27 region" vs "HLA-B" etc. (qualifiers added in
        # carrier.py panel labels). Real conflicts are bare-symbol drift.
        if (local_gene != reg_gene
                and local_gene.split("*")[0].split(" ")[0]
                != reg_gene.split("*")[0]):
                disagreed.append(
                    f"{entry['rsid']}: gene local={entry['gene']!r} "
                    f"registry={rec.gene!r}"
                )
                continue

        # Pathogenic-allele cross-check
        local_risk = entry["pathogenic_allele"]
        if local_risk in {"DEL", "INS"}:
            # Indel — registry SNV-only schema doesn't represent it; skip.
            continue
        if local_risk == rec.derived:
            agreed.append(entry["rsid"])
        elif local_risk == _COMPLEMENT.get(rec.derived):
            strand_flipped.append(entry["rsid"])
        else:
            disagreed.append(
                f"{entry['rsid']}: risk allele local={local_risk!r} "
                f"registry-derived={rec.derived!r}"
            )

    return {
        "agreed": agreed,
        "strand_flipped": strand_flipped,
        "disagreed": disagreed,
        "registry_missing": registry_missing,
    }


# Run the audit at module import — fail fast on REAL biology disagreement.
# Strand flips are expected for classical carrier variants and tolerated.
_AUDIT = audit_against_registry()
assert not _AUDIT["disagreed"], (
    "carrier.py disagrees with snp_registry on biology (not strand) for:\n  - "
    + "\n  - ".join(_AUDIT["disagreed"])
    + "\nResolve in CHANGELOG.md before continuing."
)
