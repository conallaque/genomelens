"""
Immunogenetics — Viral / Bacterial / Parasitic Resistance & Historical Selection
=================================================================================

A comprehensive genotype-based screen for genetic resistance and susceptibility
to human pathogens, plus a Historical Selection Timeline that maps the user's
protective variants to the historical pandemics that selected them.

Coverage
--------

**Viral:**
  • HIV — CCR5-Δ32 (rs333), the most consequential viral-resistance variant
    known in humans. Δ32/Δ32 = near-total resistance to R5-tropic HIV; Δ32/+ =
    ~50% slower AIDS progression.
  • Noroviruses & rotaviruses — FUT2 secretor status (rs601338). Non-secretors
    are essentially resistant to GII.4 norovirus (the "48-hour stomach flu")
    and rotavirus P[8].
  • Hepatitis B — HLA-DPB1 rs9277535 (clearance vs chronicity).
  • Hepatitis C — IL28B/IFNL4 rs12979860 (spontaneous clearance).
  • SARS-CoV-2 — OAS1 rs2660 (antiviral response), 3p21.31 Neanderthal
    haplotype (severity, Zeberg & Pääbo 2020).
  • Influenza — IFITM3 rs12252 (severe H1N1/H7N9), MX1 rs469390.
  • Prion diseases — PRNP codon 129 (rs1799990). Heterozygotes (MV) resisted
    kuru during the epidemic and are relatively protected against CJD.

**Bacterial:**
  • Yersinia pestis (plague) — ERAP2 rs2549794 (Klunk 2022 Nature), the allele
    under strong positive selection during the Black Death.
  • Sepsis / RSV — TLR4 rs4986790, rs4986791.
  • Autoinflammatory receptors — TLR1 rs4833095.

**Parasitic:**
  • P. vivax malaria — Duffy (DARC/ACKR1) rs2814778 (Duffy-null = near-total
    resistance; near-fixed in West/Central Africans).
  • P. falciparum malaria — sickle-cell HbS (HBB rs334) and G6PD deficiency
    (rs1050828). Heterozygote advantage for both against severe malaria.

**Autoimmune trade-offs (HLA & non-HLA):**
  • PTPN22 R620W (rs2476601) — general autoimmunity risk locus.
  • STAT4 rs7574865 — lupus / RA.
  • IRF5 rs4833095 — general autoimmunity.
  • NOD2 tag (Crohn's).
  • HLA-DR/DQ trade-offs (checked via the existing HLA module).

**Historical Selection Timeline** — a narrative + rank list mapping the user's
detected protective variants to the pandemics/environmental pressures that
selected them (Black Death, malaria endemic zones, early dairying, endemic
noroviruses, etc.). Some of the strongest positive-selection signals in the
human genome sit on this list.

Citations
---------
Samson 1996 CCR5-Δ32; Dean 1996; Novembre 2005 Δ32 European gradient.
Lindesmith 2003 (FUT2 & norovirus); Kelly 1995 secretor status.
Ge 2009 & Prokunina-Olsson 2013 (IL28B / IFNL4 hep-C clearance).
Mead 2003 & 2009 (PRNP codon 129 & kuru resistance).
Klunk et al. Nature 2022 (ERAP2 & Black Death positive selection).
Zeberg & Pääbo Nature 2020 (3p21.31 Neanderthal COVID haplotype).
Everitt 2012 (IFITM3 rs12252 & severe H1N1).
Miller 1976 & Zimmerman 1976 (Duffy-null & P. vivax malaria).
Allison 1954 (HbS heterozygote advantage & P. falciparum malaria).
"""

from __future__ import annotations

from typing import Dict, List, Optional
import pandas as pd


CAT_VIRAL = "Viral Resistance & Susceptibility"
CAT_BACTERIAL = "Bacterial / Sepsis"
CAT_PARASITIC = "Parasitic (Malaria)"
CAT_AUTOIMMUNE = "Autoimmune Trade-offs"


def _gt(df: pd.DataFrame, rsid: str) -> Optional[str]:
    if rsid not in df.index:
        return None
    row = df.loc[rsid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    g = row.get("genotype")
    if g is None:
        return None
    s = str(g).upper().replace(" ", "").replace("-", "")
    return s or None


def _finding(category, name, gene, rsid, genotype, impact, verdict, mechanism,
             action, confidence, citation, historical=None):
    return {
        "category": category, "name": name, "gene": gene, "rsid": rsid,
        "genotype": genotype or "—", "impact": impact, "verdict": verdict,
        "mechanism": mechanism, "action": action, "confidence": confidence,
        "citation": citation, "historical": historical,
    }


# ─── VIRAL ────────────────────────────────────────────────────────────────────

def _ccr5_delta32(df):
    gt = _gt(df, "rs333")
    if gt is None:
        return None
    n_del = gt.count("D") if "D" in gt else (2 if gt == "--" else 0)
    if n_del >= 2:
        return _finding(
            CAT_VIRAL, "HIV resistance — CCR5-Δ32 homozygous", "CCR5", "rs333", gt,
            "protective",
            "NEAR-TOTAL RESISTANCE to R5-tropic HIV",
            "You have two copies of the 32-bp CCR5 deletion, which removes the "
            "co-receptor that most HIV strains use to enter T-cells. "
            "Homozygotes are functionally resistant to HIV infection — this is "
            "the genotype behind Timothy Ray Brown's HIV cure via bone-marrow "
            "transplant. ~1% of Northern Europeans.",
            "Standard STI prevention still applies (other STIs, X4-tropic HIV "
            "strains rarely). But your genotype is remarkable.",
            "high", "Samson 1996; Dean 1996; Novembre 2005",
            historical="Selected during pre-Holocene epidemics of R5-tropic pathogens; "
                       "steep frequency cline in Europe with peaks in the Baltic (Novembre 2005). "
                       "Popular hypothesis links it to Black Death survival — debated.",
        )
    if n_del == 1:
        return _finding(
            CAT_VIRAL, "HIV protection — CCR5-Δ32 heterozygous", "CCR5", "rs333", gt,
            "protective",
            "~50% slower HIV progression + partial infection resistance",
            "You carry ONE 32-bp deletion in CCR5. Heterozygotes who become "
            "infected progress to AIDS ~50% more slowly and show reduced viral "
            "replication. ~10% of Northern Europeans; near-absent outside Europe. "
            "One of the strongest positive-selection signals in the human genome.",
            "Materially reduces (not eliminates) HIV risk. Standard prevention "
            "still applies.",
            "high", "Samson 1996; Dean 1996; Novembre 2005",
            historical="Selected during pre-Holocene epidemics; strong European gradient. "
                       "Popular hypothesis links it to plague survival — debated.",
        )
    return _finding(
        CAT_VIRAL, "CCR5-Δ32", "CCR5", "rs333", gt, "neutral",
        "No CCR5 deletion — standard HIV susceptibility",
        "You carry two copies of the standard CCR5 receptor. Standard HIV risk.",
        "Standard STI prevention.", "high", "Samson 1996", None)


def _fut2_norovirus(df):
    gt = _gt(df, "rs601338")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A == 2:
        return _finding(
            CAT_VIRAL, "Norovirus & rotavirus resistance — FUT2 non-secretor",
            "FUT2", "rs601338", gt, "protective",
            "NEAR-COMPLETE RESISTANCE to GII.4 norovirus + strong rotavirus P[8] protection",
            "Homozygous FUT2 W143X → you don't express ABO/H antigens on gut "
            "epithelium. GII.4 norovirus (the 'stomach-flu cruise-ship pathogen') "
            "requires these antigens to bind and infect you. ~20% of Europeans "
            "share this genotype. Trade-off: modestly higher H. pylori and "
            "Crohn's susceptibility.",
            "You are practically the person on the cruise ship who doesn't get sick.",
            "high", "Lindesmith 2003; Kelly 1995",
            historical="Non-secretor allele has been maintained by balancing selection — "
                       "protects against gut viral pandemics but costs slight infection risk elsewhere.",
        )
    return _finding(
        CAT_VIRAL, "FUT2 secretor status", "FUT2", "rs601338", gt,
        "susceptible" if n_A == 0 else "intermediate",
        "Normal secretor — standard norovirus / rotavirus susceptibility"
        if n_A == 0 else "Heterozygous carrier — still a secretor phenotypically",
        "You express ABO antigens on gut epithelium — norovirus and rotavirus "
        "can bind and infect. Standard risk.",
        "Standard food-safety hygiene during outbreaks.",
        "high", "Lindesmith 2003", None)


def _il28b_hepc(df):
    gt = _gt(df, "rs12979860")
    if gt is None:
        return None
    n_C = gt.count("C")
    if n_C == 2:
        return _finding(
            CAT_VIRAL, "Hepatitis C spontaneous clearance — IL28B/IFNL4",
            "IL28B (IFNL4)", "rs12979860", gt, "protective",
            "~3× odds of spontaneous HCV clearance without treatment",
            "The IFNL4 rs12979860 CC genotype produces a more effective type-III "
            "interferon response to hepatitis C. CC individuals spontaneously "
            "clear an acute HCV infection ~55% of the time (vs ~15-20% for TT). "
            "This was also the strongest pre-treatment predictor of HCV therapy "
            "response before direct-acting antivirals.",
            "If you were ever HCV-exposed, your body is genetically primed to "
            "clear it. Modern DAAs cure >95% regardless of genotype now, but "
            "this is a genuine innate advantage.",
            "high", "Ge 2009 Nature; Prokunina-Olsson 2013 Nature Genet.",
            historical="One of the strongest interferon-locus polymorphisms in human "
                       "populations — reflects long co-evolution with hepatotropic viruses.",
        )
    if n_C == 1:
        verdict = "Intermediate HCV clearance capacity"
    else:
        verdict = "Reduced HCV clearance likelihood"
    return _finding(
        CAT_VIRAL, "IL28B/IFNL4 hepatitis C response", "IL28B (IFNL4)",
        "rs12979860", gt, "intermediate" if n_C == 1 else "susceptible",
        verdict,
        "IFNL4 rs12979860 T-carriers have a less effective type-III interferon "
        "response to HCV; spontaneous clearance rates are lower.",
        "Not action-guiding today (DAAs cure regardless), but relevant if you "
        "have known HCV exposure history.",
        "high", "Ge 2009; Prokunina-Olsson 2013", None)


def _hla_dpb1_hbv(df):
    gt = _gt(df, "rs9277535")
    if gt is None:
        return None
    n_G = gt.count("G")
    if n_G == 2:
        return _finding(
            CAT_VIRAL, "Chronic hepatitis B risk — HLA-DPB1", "HLA-DPB1",
            "rs9277535", gt, "susceptible",
            "Higher risk of chronic HBV persistence (if ever exposed)",
            "The HLA-DPB1 rs9277535 G/G genotype is associated with reduced HBV "
            "clearance and higher chance of chronic HBV. Modest per-allele effect "
            "but well-replicated (Kamatani 2009). This is only clinically relevant "
            "if you have been HBV-exposed and are unvaccinated.",
            "Ensure HBV vaccination is documented; check anti-HBs titer if in doubt.",
            "moderate", "Kamatani 2009 Nature Genet.", None)
    return _finding(
        CAT_VIRAL, "HLA-DPB1 hepatitis B response", "HLA-DPB1", "rs9277535",
        gt, "protective" if n_G == 0 else "neutral",
        "Favorable HBV clearance profile" if n_G == 0 else "Standard HBV clearance profile",
        "Your HLA-DPB1 genotype favors HBV clearance / lower chronicity risk.",
        "Standard hepatitis-B vaccination remains recommended.",
        "moderate", "Kamatani 2009", None)


def _oas1(df):
    gt = _gt(df, "rs2660")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A == 2:
        return _finding(
            CAT_VIRAL, "SARS-CoV-2 & RNA-virus response — OAS1", "OAS1",
            "rs2660", gt, "protective",
            "Better interferon-driven antiviral RNA cleavage",
            "OAS1 activates RNase L to degrade viral RNA. Your A/A genotype is "
            "linked to a more effective interferon response against SARS-CoV-2 "
            "(and other RNA viruses in general). Effect size is modest but "
            "replicated across large COVID GWAS.",
            "Standard viral precautions still apply — this shifts susceptibility "
            "distributions, not individual outcomes.",
            "moderate", "Zhou 2021 Nat. Med. (COVID host-genetics initiative)",
            historical="The favorable OAS1 haplotype was itself introgressed from archaic "
                       "humans and has been positively selected — a rare Neanderthal legacy that helps us.",
        )
    return None


def _covid_3p21(df):
    ancestral = 0
    typed = 0
    for rs, prot in [("rs35044562","A"), ("rs17713054","G"), ("rs13098911","C")]:
        gt = _gt(df, rs)
        if gt:
            typed += 1
            if all(a == prot for a in gt if a in "ACGT"):
                ancestral += 1
    if typed == 0:
        return None
    if ancestral == typed:
        return _finding(
            CAT_VIRAL, "Severe COVID protection — 3p21.31 Neanderthal haplotype absent",
            "LZTFL1 / 3p21.31", "rs35044562+rs17713054+rs13098911",
            f"all ancestral × {typed}", "protective",
            "Free of the Neanderthal severe-COVID risk haplotype (~2× lower severe risk)",
            "Zeberg & Pääbo (Nature 2020) identified a ~50 kb Neanderthal-derived "
            "haplotype at 3p21.31 that ~doubles risk of severe COVID / ICU / death. "
            "Your genome carries NONE of it (all three tag markers ancestral). "
            "This is a severity marker, not a susceptibility one.",
            "If infected, expect a milder course from this locus alone. Doesn't "
            "change infection prevention.",
            "high", "Zeberg & Pääbo 2020 Nature", None)
    return None


def _ifitm3_flu(df):
    gt = _gt(df, "rs12252")
    if gt is None:
        return None
    # rs12252 T/C SNP — many chips report on the antisense strand where C→G, T→A.
    # A "GG" or "GA" antisense genotype maps to CC/CT sense (risk).
    sense_gt = gt.translate(str.maketrans("ACGT", "TGCA")) if set(gt) <= {"A","G"} else gt
    n_C = sense_gt.count("C")
    if n_C >= 2:
        return _finding(
            CAT_VIRAL, "Severe influenza risk — IFITM3", "IFITM3", "rs12252",
            f"{gt} (sense-strand CC)", "susceptible",
            "Raised risk of severe H1N1 / H7N9 avian influenza",
            "IFITM3 restricts viral entry to endosomes. The rs12252-C variant "
            "produces a truncated protein with reduced antiviral activity. "
            "Homozygotes have materially higher rates of severe hospitalization "
            "during H1N1 pandemics (Everitt 2012).",
            "Take annual seasonal flu vaccination seriously — this is a locus "
            "where vaccination materially matters for you.",
            "moderate", "Everitt 2012 Nature; Xuan 2013 Cell",
            historical="rs12252-C is a Southeast-Asian-enriched allele that has been "
                       "under recent purifying selection in Europe.")
    if n_C == 1:
        return _finding(
            CAT_VIRAL, "IFITM3 flu response (heterozygous)", "IFITM3",
            "rs12252", gt, "intermediate",
            "Modest increase in severe-flu risk",
            "One copy of the truncating IFITM3 variant.",
            "Standard flu vaccination.", "moderate",
            "Everitt 2012; Xuan 2013", None)
    return None


def _prnp_prion(df):
    gt = _gt(df, "rs1799990")
    if gt is None:
        return None
    letters = set(gt)
    if letters == {"A","G"} or letters == {"C","T"}:
        # heterozygous MV
        return _finding(
            CAT_VIRAL, "Prion disease resistance — PRNP codon 129 MV heterozygous",
            "PRNP", "rs1799990", gt, "protective",
            "Robust protection against prion diseases (CJD, kuru)",
            "PRNP codon 129 encodes methionine (M) or valine (V). Heterozygotes "
            "(M/V) are strongly under-represented among people with sporadic and "
            "variant Creutzfeldt-Jakob disease AND were the primary survivors of "
            "the New Guinea kuru cannibalism-related prion epidemic. "
            "You are heterozygous — this is one of the clearest examples of "
            "heterozygote advantage documented in humans.",
            "Peace of mind. Sporadic CJD is rare regardless.",
            "high", "Mead 2003 Science (kuru); Palmer 1991; Mead 2009 NEJM",
            historical="The strongest genetic signature of a historical prion epidemic "
                       "found in modern populations. MV heterozygotes disproportionately "
                       "survived the kuru epidemic in the Fore people (Papua New Guinea).")
    # Homozygous MM or VV
    return _finding(
        CAT_VIRAL, "PRNP codon 129 homozygous", "PRNP", "rs1799990", gt,
        "susceptible",
        "Standard prion susceptibility — MM or VV homozygous",
        "Homozygotes at PRNP codon 129 have higher sporadic CJD susceptibility "
        "than heterozygotes. Sporadic CJD is still extremely rare (~1 in a million).",
        "No specific action — the baseline risk is very low.",
        "high", "Mead 2003; Palmer 1991", None)


def _mx1_flu(df):
    gt = _gt(df, "rs469390")
    if gt is None:
        return None
    return _finding(
        CAT_VIRAL, "MX1 antiviral response", "MX1", "rs469390", gt,
        "informational",
        f"MX1 rs469390 genotype {gt}",
        "MX1 is an interferon-induced antiviral GTPase that restricts influenza "
        "A and other RNA viruses. Variants at this locus have modest effects on "
        "flu susceptibility.",
        "General flu-vaccination recommendations apply.",
        "low", "Ciancanelli 2015; Haller 2015", None)


# ─── BACTERIAL ────────────────────────────────────────────────────────────────

def _erap2_plague(df):
    gt = _gt(df, "rs2549794")
    if gt is None:
        return None
    n_C = gt.count("C")
    if n_C >= 1:
        return _finding(
            CAT_BACTERIAL, "Black Death survivor allele — ERAP2", "ERAP2",
            "rs2549794", gt, "protective",
            f"Carries {n_C}× C — allele under strong positive selection during the Black Death",
            "Klunk et al. (Nature 2022) analysed ancient DNA from pre- and "
            "post-plague London and Denmark cemeteries and found rs2549794-C at "
            "much higher frequency in survivors — an ~40% survival advantage per "
            "copy, one of the strongest positive-selection events documented in "
            "recent humans. ERAP2 trims peptides for MHC-I antigen presentation. "
            "Trade-off: modestly raises Crohn's disease risk.",
            "Peace of mind. Also flags that your immune system is genetically "
            "biased toward strong intracellular-pathogen presentation — good in "
            "most contexts.",
            "high", "Klunk et al. 2022 Nature",
            historical="rs2549794-C rose from ~40% to ~50% in ~4 generations during the "
                       "Black Death (1348-50). One of the strongest documented episodes "
                       "of natural selection in modern human history.")
    return None


def _tlr4_sepsis(df):
    gt790 = _gt(df, "rs4986790")
    gt791 = _gt(df, "rs4986791")
    if not gt790 and not gt791:
        return None
    n_G = gt790.count("G") if gt790 else 0
    n_T = gt791.count("T") if gt791 else 0
    if n_G >= 1 or n_T >= 1:
        return _finding(
            CAT_BACTERIAL, "TLR4 endotoxin response — reduced", "TLR4",
            "rs4986790+rs4986791", f"{gt790 or '-'} / {gt791 or '-'}",
            "susceptible",
            "Blunted LPS response — raised risk of severe gram-negative sepsis / RSV",
            "TLR4 Asp299Gly / Thr399Ile carriers have reduced signalling in "
            "response to bacterial lipopolysaccharide. Modestly higher severe "
            "sepsis and RSV bronchiolitis rates.",
            "Awareness only. Standard antibiotic care applies.",
            "moderate", "Arbour 2000; Awomoyi 2007", None)
    return _finding(
        CAT_BACTERIAL, "TLR4 endotoxin response — normal", "TLR4",
        "rs4986790+rs4986791", f"{gt790 or '-'} / {gt791 or '-'}", "neutral",
        "Normal TLR4 signalling — standard sepsis / RSV response",
        "Both TLR4 variants are wild-type. Standard endotoxin sensing.",
        "None.", "moderate", "Arbour 2000", None)


# ─── PARASITIC ────────────────────────────────────────────────────────────────

def _duffy_malaria(df):
    gt = _gt(df, "rs2814778")
    if gt is None:
        return None
    n_C = gt.count("C")
    if n_C >= 2:
        return _finding(
            CAT_PARASITIC, "P. vivax malaria — near-total resistance (Duffy-null)",
            "DARC (ACKR1)", "rs2814778", gt, "protective",
            "Near-complete resistance to Plasmodium vivax malaria",
            "The Duffy-null genotype prevents P. vivax from binding to red-cell "
            "surfaces. Duffy-null is near-fixed in West/Central Africa (>95%) "
            "and near-absent in Europe.",
            "Peace of mind if visiting P. vivax endemic areas.",
            "high", "Miller 1976 NEJM; Zimmerman 1976",
            historical="One of the strongest positive-selection signals in the entire "
                       "human genome — driven by millennia of P. vivax pressure in Africa.")
    return _finding(
        CAT_PARASITIC, "Duffy antigen — malaria-susceptible", "DARC (ACKR1)",
        "rs2814778", gt, "susceptible",
        "Duffy-positive — standard P. vivax malaria susceptibility",
        "You have functional Duffy antigen on red cells — P. vivax can invade.",
        "Standard malaria chemoprophylaxis if visiting endemic areas.",
        "high", "Miller 1976", None)


def _hbs_sickle(df):
    gt = _gt(df, "rs334")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A == 1:
        return _finding(
            CAT_PARASITIC, "Sickle-cell trait (HbAS) — malaria protection",
            "HBB", "rs334", gt, "protective",
            "Heterozygote advantage — significant P. falciparum malaria resistance",
            "One copy of HbS. Confers ~90% reduction in severe P. falciparum "
            "malaria (classic heterozygote advantage — Allison 1954). Trade-off: "
            "risk of sickling under extreme hypoxic stress (very high altitude, "
            "sub-optimal anesthesia). Two copies = sickle-cell disease.",
            "Notify anesthesiologists before surgery. Careful with extreme "
            "altitude / dehydration.",
            "high", "Allison 1954; Aidoo 2002 Lancet", None)
    if n_A == 2:
        return _finding(
            CAT_PARASITIC, "Sickle-cell disease (HbSS) — homozygous",
            "HBB", "rs334", gt, "susceptible",
            "Homozygous HbS — sickle-cell disease",
            "Two copies of HbS. Chronic haemolysis, vaso-occlusive crises, "
            "organ damage. Would already be under haematology care.",
            "Haematology follow-up.", "high", "Allison 1954", None)
    return None   # HbAA is majority; not worth flagging.


def _g6pd_malaria(df):
    gt = _gt(df, "rs1050828")
    if gt is None:
        return None
    n_T = gt.count("T")
    # X-linked; hemizygous in males counts as 1 with C→T
    if n_T >= 1:
        return _finding(
            CAT_PARASITIC, "G6PD A- deficiency — malaria protection", "G6PD",
            "rs1050828", gt, "protective",
            "Modest P. falciparum protection; oxidant-drug caution",
            "Reduced G6PD activity confers ~40% protection against severe P. "
            "falciparum malaria (heterozygote advantage). Trade-off: hemolysis "
            "risk with oxidant drugs (primaquine, dapsone, sulfa) and fava beans.",
            "AVOID fava beans, primaquine, rasburicase, high-dose vitamin C. "
            "Flag G6PD status to prescribers.",
            "high", "Tishkoff 2001; Ruwende 1995 Nature",
            historical="Positive selection in African/Mediterranean populations by "
                       "P. falciparum malaria — one of the earliest and clearest examples "
                       "of gene×pathogen selection.")
    return None


# ─── AUTOIMMUNE TRADE-OFFS ────────────────────────────────────────────────────

def _ptpn22(df):
    gt = _gt(df, "rs2476601")
    if gt is None:
        return None
    n_A = gt.count("A")
    if n_A >= 1:
        return _finding(
            CAT_AUTOIMMUNE, "PTPN22 R620W — general autoimmunity risk",
            "PTPN22", "rs2476601", gt, "susceptible",
            f"Carrier of the general-autoimmunity risk allele ({n_A}× A)",
            "PTPN22 620W is one of the strongest common autoimmunity risk "
            "variants known — associated with type-1 diabetes, rheumatoid "
            "arthritis, lupus, Grave's, and juvenile idiopathic arthritis "
            "(OR ~1.5-2.0 per allele).",
            "Awareness of autoimmune-symptom red flags. Not preventive.",
            "high", "Bottini 2004 Nature Genet.", None)
    return None


def _stat4(df):
    gt = _gt(df, "rs7574865")
    if gt is None:
        return None
    n_T = gt.count("T")
    if n_T >= 1:
        return _finding(
            CAT_AUTOIMMUNE, "STAT4 — lupus / RA risk", "STAT4", "rs7574865", gt,
            "susceptible", f"{n_T}× T — modest lupus / RA risk allele",
            "STAT4 rs7574865-T is a replicated risk allele for SLE, rheumatoid "
            "arthritis and primary Sjögren's syndrome.",
            "None specific; part of an autoimmunity-risk profile.",
            "moderate", "Remmers 2007 NEJM", None)
    return None


def _irf5(df):
    gt = _gt(df, "rs4833095")
    if gt is None:
        return None
    return _finding(
        CAT_AUTOIMMUNE, "IRF5 / TLR1 region", "TLR1/IRF5", "rs4833095", gt,
        "informational",
        "Common autoimmunity / mycobacterial-susceptibility locus",
        "This region tags variants in TLR1 (leprosy susceptibility) and "
        "nearby IRF5 (autoimmunity).",
        "None specific.", "low", "multiple GWAS", None)


# ─── Historical selection narrative ──────────────────────────────────────────

def build_selection_timeline(findings: List[Dict]) -> List[Dict]:
    """Map protective findings to historical selection events, sorted by
    approximate epoch (oldest first)."""
    events: List[Dict] = []
    for f in findings:
        if not f.get("historical"):
            continue
        # Assign a rough epoch tag by keyword
        text = f["historical"].lower()
        if "black death" in text or "plague" in text:
            epoch = ("1348-1350 CE", "Bacterial", "Yersinia pestis (Black Death)")
        elif "malaria" in text or "falciparum" in text or "vivax" in text:
            epoch = ("~5,000-100,000 ya", "Parasitic", "Endemic malaria (P. falciparum / P. vivax)")
        elif "kuru" in text or "prion" in text:
            epoch = ("Recorded 20th c.; ancient origin", "Prion / dietary",
                     "Prion epidemics (kuru; ancient cannibalism)")
        elif "dairying" in text or "lactase" in text:
            epoch = ("~7,000 ya", "Dietary", "Neolithic dairy farming")
        elif "neanderthal" in text or "introgress" in text:
            epoch = ("~50,000 ya", "Archaic introgression", "Neanderthal → modern-human gene flow")
        elif "norovirus" in text or "gut" in text:
            epoch = ("Neolithic onwards", "Viral / gut pathogens",
                     "Gut-virus co-evolution (norovirus, rotavirus)")
        elif "hepatotrop" in text or "hepatitis" in text:
            epoch = ("Ancient", "Viral", "Hepatotropic viruses")
        else:
            epoch = ("Ancient", "General", "Balancing / positive selection")
        events.append({
            "epoch": epoch[0], "category": epoch[1], "driver": epoch[2],
            "finding": f["name"], "verdict": f["verdict"],
            "impact": f["impact"], "narrative": f["historical"],
        })
    return events


# ─── Master analyzer ─────────────────────────────────────────────────────────

def analyze_immunogenetics(df: pd.DataFrame) -> Dict:
    """Full immunogenetics work-up + Historical Selection Timeline."""
    analyzers = [
        # Viral
        _ccr5_delta32, _fut2_norovirus, _il28b_hepc, _hla_dpb1_hbv,
        _oas1, _covid_3p21, _ifitm3_flu, _prnp_prion, _mx1_flu,
        # Bacterial
        _erap2_plague, _tlr4_sepsis,
        # Parasitic
        _duffy_malaria, _hbs_sickle, _g6pd_malaria,
        # Autoimmune
        _ptpn22, _stat4, _irf5,
    ]
    findings: List[Dict] = []
    for a in analyzers:
        try:
            r = a(df)
        except Exception:
            continue
        if r is None:
            continue
        findings.append(r)

    by_category: Dict[str, List[Dict]] = {}
    for f in findings:
        by_category.setdefault(f["category"], []).append(f)

    n_protective = sum(1 for f in findings if f["impact"] == "protective")
    n_susceptible = sum(1 for f in findings if f["impact"] == "susceptible")
    n_intermediate = sum(1 for f in findings if f["impact"] == "intermediate")

    # Pull out headline resistances (only "protective" with "near-total" or
    # explicit strong wording) for the hero card.
    headlines = [f for f in findings if f["impact"] == "protective"
                 and any(kw in f["verdict"].lower() for kw in
                         ("near-total", "near-complete", "spontaneous",
                          "robust", "significant"))]

    return {
        "available": bool(findings),
        "n_findings": len(findings),
        "n_protective": n_protective,
        "n_susceptible": n_susceptible,
        "n_intermediate": n_intermediate,
        "findings": findings,
        "by_category": by_category,
        "categories": [CAT_VIRAL, CAT_BACTERIAL, CAT_PARASITIC, CAT_AUTOIMMUNE],
        "headlines": headlines,
        "historical_timeline": build_selection_timeline(findings),
    }
