"""
Unified SNP Registry — single source of truth for variant metadata.

Why this exists
---------------
Before V7 the codebase had **seven** modules each maintaining their own dict of
rsID → (position, ancestral, derived, gene) tuples. 46 % of rsIDs in the
project appeared in two or more such dicts — sometimes with *contradictory*
metadata (one module reading a SNP on the + strand, another on the − strand,
neither agreeing on the GRCh37 position). The strand-handling bug fixed in
`supplements._risk_dose` was one symptom; the disease is the duplication.

This module is the **only** place rsID metadata may be hard-coded going
forward. Every consumer module (`carrier`, `prs`, `pgx`, `traits`, …) should
migrate to ``snp_registry.get(rsid)`` instead of maintaining its own dict.

Design choices
--------------

* **`@dataclass(frozen=True)` not Pydantic.** We don't need network-deserialise
  semantics; we want speed and zero dependencies for what is essentially a
  read-mostly in-process lookup. ``frozen=True`` makes records hashable and
  prevents accidental mutation across modules.

* **Multi-build coordinates.** Every record carries both GRCh37 and GRCh38
  positions where known. Position-based fallback lookups (used when an rsID
  is reported under a non-standard label by a chip vendor) can match either.

* **Strand-aware dose helper lives here.** Same logic that landed in
  ``supplements._risk_dose`` — but at the registry level it gets a single
  authoritative implementation. ``supplements`` (and every other module) will
  re-export from here in the migration phase.

* **Provenance fields are first-class.** Every record carries the GWAS source
  citation and a `last_verified` date so a stale variant can be detected by
  the audit script (see `audit_registry()` at the bottom of this module).
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


# ── Core record ──────────────────────────────────────────────────────────────

# Single canonical ancestral/derived assignment per rsID, ALWAYS reported on
# the + strand. Chip-level strand differences are resolved by `risk_dose`.
@dataclass(frozen=True, slots=True)
class SNPRecord:
    """Canonical metadata for one rsID. Frozen → hashable + thread-safe."""

    rsid: str
    gene: str
    chrom: str
    # GRCh37/hg19 — the build most consumer chips still report against.
    pos_grch37: int | None
    # GRCh38/hg38 — populated where verifiable. Consumers that detect the
    # chip's build at parse time can use whichever applies.
    pos_grch38: int | None
    # All allele letters refer to the + strand of the reference assembly.
    ancestral: str
    derived: str
    # One-line description, used by the supplements / phewas / wellness
    # narratives to surface "what this SNP does" without each module
    # redefining it.
    description: str
    # GWAS / consortium / paper this record's effect direction is sourced
    # from. Free text, not parsed — purely for audit.
    source: str = ""
    # ISO-8601 date this record was last reconciled with the source. The
    # `audit_registry()` helper flags records older than `STALE_AFTER_DAYS`.
    last_verified: str = "2026-05-17"
    # Optional fields used by specific consumer modules — kept here so the
    # consumer doesn't need its own parallel dict.
    aliases: tuple[str, ...] = field(default_factory=tuple)
    clinical_significance: str = ""

    def __post_init__(self) -> None:  # validation only — frozen so no mutation
        if not self.rsid.startswith("rs"):
            raise ValueError(f"rsid must start with 'rs': {self.rsid!r}")
        if self.ancestral not in {"A", "T", "C", "G"}:
            raise ValueError(
                f"{self.rsid}: ancestral must be A/T/C/G on + strand, "
                f"got {self.ancestral!r}"
            )
        if self.derived not in {"A", "T", "C", "G"}:
            raise ValueError(
                f"{self.rsid}: derived must be A/T/C/G on + strand, "
                f"got {self.derived!r}"
            )
        if self.ancestral == self.derived:
            raise ValueError(f"{self.rsid}: ancestral == derived")


# ── Strand-aware dose helper (single canonical implementation) ──────────────

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement(allele: str) -> str:
    """Single-base complement on the reverse strand."""
    return _COMPLEMENT.get(allele.upper(), allele.upper())


def risk_dose(
    genotype: str | None,
    *,
    risk_allele: str,
    ref_allele: str,
) -> int | None:
    """
    Strand-aware risk-allele dosage (0/1/2) from a 2-char genotype string.

    Auto-detects whether the chip reports on the + or − strand by checking
    which of {risk, ref, complement(risk), complement(ref)} appear in the
    genotype, then counts the right side. Returns ``None`` for missing or
    no-call genotypes.

    The legacy ``count(risk) + count(complement(risk))`` antipattern would
    double-count a hypothetical mixed-strand "AT" genotype as homozygous;
    this implementation explicitly rejects mixed-strand calls.
    """
    if not genotype or len(genotype) != 2 or genotype in {"--", "00", "NN"}:
        return None

    r = risk_allele.upper()
    a = ref_allele.upper()
    r_c, a_c = complement(r), complement(a)

    plus_alleles = {r, a}
    minus_alleles = {r_c, a_c}
    gt_alleles = set(genotype.upper())

    if gt_alleles.issubset(plus_alleles):
        return genotype.upper().count(r)
    if gt_alleles.issubset(minus_alleles):
        return genotype.upper().count(r_c)
    # Mixed alleles (palindromic A/T or C/G SNPs, or chip noise) — fall back
    # to + strand interpretation rather than raise.
    return genotype.upper().count(r)


def risk_dose_from_df(
    snps_df: pd.DataFrame | None,
    rsid: str,
    *,
    risk_allele: str | None = None,
    ref_allele: str | None = None,
) -> int | None:
    """
    Same as :func:`risk_dose` but reads the genotype straight from a
    ``parse_dna_file``-style DataFrame indexed by rsID.

    If ``risk_allele`` / ``ref_allele`` are omitted, they are taken from the
    canonical :class:`SNPRecord` in :data:`SNPS`.
    """
    if snps_df is None or rsid not in snps_df.index:
        return None

    raw = snps_df.loc[rsid].get("genotype")
    if raw is None:
        return None
    gt = str(raw).strip().upper().replace(" ", "").replace("-", "")
    if not gt:
        return None
    # Collapse "AA" → handled fine by 2-char code; preserve hemizygous "A" → "AA"
    if len(gt) == 1:
        gt = gt + gt

    if risk_allele is None or ref_allele is None:
        rec = SNPS.get(rsid)
        if rec is None:
            return None
        risk_allele = risk_allele or rec.derived
        ref_allele = ref_allele or rec.ancestral

    return risk_dose(gt, risk_allele=risk_allele, ref_allele=ref_allele)


# ── The registry itself ──────────────────────────────────────────────────────
#
# Records below are the seed set for V7. Subsequent module migrations will
# add to this dict; once every consumer module imports from here we will
# delete the duplicate per-module dicts. We intentionally start small and
# migrate by domain — see `MIGRATION_PLAN.md` (or the comments in this file)
# for the order.

_RECORDS: list[SNPRecord] = [
    # ── Methylation cycle ────────────────────────────────────────────────
    SNPRecord(
        rsid="rs1801133", gene="MTHFR", chrom="1",
        pos_grch37=11_856_378, pos_grch38=11_796_321,
        ancestral="C", derived="T",
        description="MTHFR C677T — reduces methylenetetrahydrofolate reductase "
                    "activity ~30%/70% per allele copy.",
        source="Frosst 1995 (NEJM); CPIC 2019",
        aliases=("C677T",),
    ),
    SNPRecord(
        rsid="rs1801131", gene="MTHFR", chrom="1",
        pos_grch37=11_854_476, pos_grch38=11_794_419,
        ancestral="A", derived="C",
        description="MTHFR A1298C — reduces BH4 regeneration affecting "
                    "neurotransmitter & methylation pathways.",
        source="van der Put 1998",
        aliases=("A1298C",),
    ),
    SNPRecord(
        rsid="rs4680", gene="COMT", chrom="22",
        pos_grch37=19_951_271, pos_grch38=19_963_748,
        ancestral="G", derived="A",
        description="COMT Val158Met — A allele (Met) reduces enzyme activity "
                    "200-300%; slow catechol clearance.",
        source="Lachman 1996; Tunbridge 2006",
        aliases=("Val158Met",),
    ),

    # ── Vitamin D metabolism ─────────────────────────────────────────────
    SNPRecord(
        rsid="rs10741657", gene="CYP2R1", chrom="11",
        pos_grch37=14_914_878, pos_grch38=14_893_332,
        ancestral="A", derived="G",
        description="CYP2R1 — 25-hydroxylase. G allele = reduced 25(OH)D.",
        source="Manousaki 2017 (Plos Med)",
    ),
    SNPRecord(
        rsid="rs2282679", gene="GC", chrom="4",
        pos_grch37=72_618_323, pos_grch38=72_752_040,
        ancestral="T", derived="C",
        description="GC / vitamin-D binding protein. C allele = lower serum 25(OH)D.",
        source="Manousaki 2017",
        aliases=("VDBP",),
    ),
    SNPRecord(
        rsid="rs2228570", gene="VDR", chrom="12",
        pos_grch37=48_272_895, pos_grch38=47_879_112,
        ancestral="C", derived="T",
        description="VDR FokI — T allele produces a 3-aa-shorter, more active VDR.",
        source="Uitterlinden 2004",
        aliases=("FokI",),
    ),

    # ── B12 / iron ───────────────────────────────────────────────────────
    SNPRecord(
        rsid="rs602662", gene="FUT2", chrom="19",
        pos_grch37=49_206_985, pos_grch38=48_703_728,
        ancestral="G", derived="A",
        description="FUT2 — secretor status; A allele reduces B12 absorption.",
        source="Hazra 2009",
    ),
    SNPRecord(
        rsid="rs1801198", gene="TCN2", chrom="22",
        pos_grch37=31_019_043, pos_grch38=30_623_054,
        ancestral="C", derived="G",
        description="TCN2 Pro259Arg — transcobalamin variant affecting B12 transport.",
        source="Afman 2002",
    ),
    SNPRecord(
        rsid="rs855791", gene="TMPRSS6", chrom="22",
        pos_grch37=37_462_936, pos_grch38=37_066_896,
        ancestral="G", derived="A",
        description="TMPRSS6 — A allele = elevated hepcidin → lower serum iron/MCV.",
        source="Benyamin 2009; Astle 2016",
    ),
    SNPRecord(
        rsid="rs1800562", gene="HFE", chrom="6",
        pos_grch37=26_093_141, pos_grch38=26_092_913,
        ancestral="G", derived="A",
        description="HFE C282Y — homozygous → classical hereditary hemochromatosis.",
        source="Feder 1996",
        aliases=("C282Y",),
        clinical_significance="Pathogenic (when homozygous)",
    ),
    SNPRecord(
        rsid="rs1799945", gene="HFE", chrom="6",
        pos_grch37=26_091_179, pos_grch38=26_090_951,
        ancestral="C", derived="G",
        description="HFE H63D — mild iron-overload variant; compound het with "
                    "C282Y may be clinically relevant.",
        source="Feder 1996",
        aliases=("H63D",),
    ),

    # ── Fatty-acid metabolism / inflammation ─────────────────────────────
    SNPRecord(
        rsid="rs174547", gene="FADS1", chrom="11",
        pos_grch37=61_570_783, pos_grch38=61_803_311,
        ancestral="C", derived="T",
        description="FADS1 — T allele reduces ALA→EPA/DHA conversion rate.",
        source="Schaeffer 2006",
    ),
    SNPRecord(
        rsid="rs1800795", gene="IL6", chrom="7",
        pos_grch37=22_766_645, pos_grch38=22_727_026,
        ancestral="G", derived="C",
        description="IL6 -174 promoter — G allele = higher basal IL-6 transcription "
                    "(higher CRP).",
        source="Fishman 1998; Ligthart 2018",
    ),
    SNPRecord(
        rsid="rs1695", gene="GSTP1", chrom="11",
        pos_grch37=67_352_689, pos_grch38=67_585_218,
        ancestral="A", derived="G",
        description="GSTP1 Ile105Val — G allele = reduced GST activity (slower xenobiotic conjugation).",
        source="Watson 1998",
        aliases=("Ile105Val",),
    ),
    SNPRecord(
        rsid="rs6721961", gene="NFE2L2", chrom="2",
        pos_grch37=178_098_964, pos_grch38=177_234_237,
        ancestral="G", derived="T",
        description="NRF2 / NFE2L2 promoter — T allele reduces antioxidant-response binding.",
        source="Marzec 2007",
        aliases=("NRF2",),
    ),

    # ── Drug metabolism ──────────────────────────────────────────────────
    SNPRecord(
        rsid="rs762551", gene="CYP1A2", chrom="15",
        pos_grch37=75_041_917, pos_grch38=74_749_576,
        ancestral="A", derived="C",
        description="CYP1A2 *1F (-163 A>C) — C allele = lower inducibility "
                    "(slow caffeine metaboliser).",
        source="Sachse 1999",
        aliases=("*1F",),
    ),

    # ── Athletic / connective tissue ─────────────────────────────────────
    SNPRecord(
        rsid="rs1815739", gene="ACTN3", chrom="11",
        pos_grch37=66_328_095, pos_grch38=66_560_624,
        ancestral="C", derived="T",
        description="ACTN3 R577X — T (X) allele introduces a premature stop; "
                    "TT = α-actinin-3 deficient (endurance-biased).",
        source="Yang 2003",
        aliases=("R577X",),
    ),
    SNPRecord(
        rsid="rs6265", gene="BDNF", chrom="11",
        pos_grch37=27_679_916, pos_grch38=27_658_369,
        ancestral="G", derived="A",
        description="BDNF Val66Met — A (Met) reduces activity-dependent BDNF secretion.",
        source="Egan 2003",
        aliases=("Val66Met",),
    ),
    SNPRecord(
        rsid="rs1801260", gene="CLOCK", chrom="4",
        pos_grch37=56_412_708, pos_grch38=55_546_586,
        ancestral="A", derived="G",
        description="CLOCK 3111T>C — historical T/C nomenclature; T allele = evening preference.",
        source="Katzenberg 1998",
        aliases=("3111T/C",),
    ),

    # ── V8 migration: carrier panel SNVs ────────────────────────────────
    # Added during the carrier.py registry migration. Each record's
    # ancestral/derived assignment was reconciled against ClinVar / dbSNP +
    # the original carrier.py panel; documented agreements in CHANGELOG.md.
    # Indel variants from carrier.py (CFTR ΔF508, HEXA 1278insTATC, etc.)
    # are NOT added here — the SNPRecord schema is SNV-only. Schema
    # extension for indels is a V8.1 task; until then carrier.py keeps
    # those entries in its local CARRIER_VARIANTS dict.
    SNPRecord(
        rsid="rs6025", gene="F5", chrom="1",
        pos_grch37=169_519_049, pos_grch38=169_549_811,
        ancestral="G", derived="A",
        description="F5 Factor V Leiden (R506Q) — A allele = activated-protein-C "
                    "resistance; ~5-7× lifetime VTE risk per copy.",
        source="Bertina 1994 (Nature); ClinVar 642",
        aliases=("R506Q", "Factor V Leiden"),
        clinical_significance="Pathogenic (autosomal-dominant susceptibility)",
    ),
    SNPRecord(
        rsid="rs1799963", gene="F2", chrom="11",
        pos_grch37=46_761_055, pos_grch38=46_739_505,
        ancestral="G", derived="A",
        description="F2 / Prothrombin G20210A — A allele elevates prothrombin "
                    "levels and VTE risk.",
        source="Poort 1996; ClinVar 13310",
        aliases=("G20210A", "Prothrombin G20210A"),
        clinical_significance="Pathogenic (autosomal-dominant susceptibility)",
    ),
    SNPRecord(
        rsid="rs75932628", gene="TREM2", chrom="6",
        pos_grch37=41_129_252, pos_grch38=41_161_514,
        ancestral="C", derived="T",
        description="TREM2 R47H — T allele = late-onset Alzheimer's "
                    "susceptibility (~3× risk).",
        source="Guerreiro 2013 (NEJM); Jonsson 2013 (NEJM)",
        aliases=("R47H",),
        clinical_significance="Pathogenic susceptibility allele",
    ),
    SNPRecord(
        rsid="rs17879961", gene="CHEK2", chrom="22",
        pos_grch37=29_121_087, pos_grch38=28_725_099,
        ancestral="T", derived="C",
        description="CHEK2 I157T — C allele = elevated breast / colon / "
                    "prostate / kidney cancer risk.",
        source="Cybulski 2004; ClinVar 5605",
        aliases=("I157T",),
        clinical_significance="Pathogenic (low-penetrance)",
    ),
    SNPRecord(
        rsid="rs2476601", gene="PTPN22", chrom="1",
        pos_grch37=114_377_568, pos_grch38=113_834_946,
        ancestral="G", derived="A",
        description="PTPN22 R620W — A allele = broad autoimmunity "
                    "susceptibility (T1D, RA, lupus, Graves').",
        source="Bottini 2004 (Nat Genet); Begovich 2004",
        aliases=("R620W",),
    ),
    SNPRecord(
        rsid="rs2187668", gene="HLA-DQA1", chrom="6",
        pos_grch37=32_605_884, pos_grch38=32_638_107,
        ancestral="C", derived="T",
        description="HLA-DQ2.5 tag SNP — T allele in LD with the DQ2.5 "
                    "haplotype (coeliac-disease risk).",
        source="Monsuur 2008 (PLoS One); Karell 2003",
        aliases=("DQ2.5 tag",),
    ),
    SNPRecord(
        rsid="rs7454108", gene="HLA-DQB1", chrom="6",
        pos_grch37=32_814_869, pos_grch38=32_847_092,
        ancestral="T", derived="C",
        description="HLA-DQ8 tag SNP — C allele in LD with DQ8 haplotype "
                    "(coeliac, T1D risk).",
        source="Monsuur 2008",
        aliases=("DQ8 tag",),
    ),
    SNPRecord(
        rsid="rs334", gene="HBB", chrom="11",
        pos_grch37=5_248_232, pos_grch38=5_227_002,
        ancestral="T", derived="A",
        description="HBB sickle-cell mutation (Glu6Val on coding strand; "
                    "T>A on + strand) — A homozygotes = sickle-cell disease; "
                    "heterozygotes = trait (malaria-protective).",
        source="ClinVar 15333; Pauling 1949",
        aliases=("HbS", "Glu6Val", "rs10500170"),
        clinical_significance="Pathogenic (recessive)",
    ),
    SNPRecord(
        rsid="rs28929474", gene="SERPINA1", chrom="14",
        pos_grch37=94_844_947, pos_grch38=94_378_610,
        ancestral="C", derived="T",
        description="SERPINA1 PiZ allele (Glu342Lys) — T homozygotes = "
                    "α1-antitrypsin deficiency (lung/liver disease).",
        source="ClinVar 17968",
        aliases=("PiZ", "Glu342Lys", "Z allele"),
        clinical_significance="Pathogenic (recessive)",
    ),
    SNPRecord(
        rsid="rs5030858", gene="PAH", chrom="12",
        pos_grch37=103_245_445, pos_grch38=102_851_667,
        ancestral="C", derived="T",
        description="PAH R408W — T homozygotes / compound heterozygotes = "
                    "phenylketonuria (PKU).",
        source="ClinVar 591",
        aliases=("R408W",),
        clinical_significance="Pathogenic (recessive)",
    ),
    SNPRecord(
        rsid="rs76763715", gene="GBA", chrom="1",
        pos_grch37=155_205_634, pos_grch38=155_235_843,
        ancestral="T", derived="C",
        description="GBA N370S — C variant = Gaucher disease type I "
                    "(carrier frequency ~1 in 14 Ashkenazi Jewish).",
        source="ClinVar 4288",
        aliases=("N370S",),
        clinical_significance="Pathogenic (recessive)",
    ),
    SNPRecord(
        rsid="rs1050828", gene="G6PD", chrom="X",
        pos_grch37=153_762_633, pos_grch38=154_535_443,
        ancestral="C", derived="T",
        description="G6PD V68M (A− variant) — T allele = X-linked G6PD "
                    "deficiency (favism, drug-induced hemolysis).",
        source="ClinVar 10410",
        aliases=("V68M", "A-"),
        clinical_significance="Pathogenic (X-linked)",
    ),
    SNPRecord(
        rsid="rs5742904", gene="APOB", chrom="2",
        pos_grch37=21_229_160, pos_grch38=21_006_288,
        ancestral="C", derived="T",
        description="APOB R3500Q — T allele = familial defective ApoB-100 "
                    "(autosomal-dominant hypercholesterolaemia).",
        source="ClinVar 17890; Innerarity 1990",
        aliases=("R3500Q", "FDB"),
        clinical_significance="Pathogenic (autosomal-dominant)",
    ),

    # ── V8 migration: wellness panel SNVs ───────────────────────────────
    # Added during the wellness.py registry migration. Sources cited per
    # entry; the wellness module does not store a structured "risk allele"
    # field so the cross-check is presence-only (every wellness rsID must
    # be registered) — see wellness.audit_against_registry().
    SNPRecord(
        rsid="rs4880", gene="SOD2", chrom="6",
        pos_grch37=160_113_872, pos_grch38=159_692_840,
        ancestral="T", derived="C",
        description="SOD2 Ala16Val — C allele = mitochondrial-targeting "
                    "favoured (Ala); T = cytoplasmic (Val), lower mtSOD activity.",
        source="Sutton 2003 (Pharmacogenetics)",
        aliases=("Ala16Val",),
    ),
    SNPRecord(
        rsid="rs429358", gene="APOE", chrom="19",
        pos_grch37=45_411_941, pos_grch38=44_908_684,
        ancestral="T", derived="C",
        description="APOE ε4 SNP — C allele defines ε4 (with rs7412 C); "
                    "elevated late-onset Alzheimer's risk.",
        source="Corder 1993 (Science)",
        aliases=("APOE ε4 SNP",),
    ),
    SNPRecord(
        rsid="rs7412", gene="APOE", chrom="19",
        pos_grch37=45_412_079, pos_grch38=44_908_822,
        ancestral="C", derived="T",
        description="APOE ε2 SNP — T allele defines ε2; reduces LDL "
                    "binding and lowers AD risk.",
        source="Corder 1993",
        aliases=("APOE ε2 SNP",),
    ),
    SNPRecord(
        rsid="rs8192678", gene="PPARGC1A", chrom="4",
        pos_grch37=23_815_662, pos_grch38=23_813_945,
        ancestral="G", derived="A",
        description="PPARGC1A Gly482Ser — G = Gly482 (better aerobic "
                    "trainability); A = Ser482.",
        source="Lucia 2005",
        aliases=("Gly482Ser",),
    ),
    SNPRecord(
        rsid="rs1800012", gene="COL1A1", chrom="17",
        pos_grch37=48_275_363, pos_grch38=50_198_002,
        ancestral="G", derived="T",
        description="COL1A1 Sp1-binding site — T allele = altered "
                    "collagen type-I composition; soft-tissue injury risk.",
        source="Posthumus 2009",
        aliases=("Sp1",),
    ),
    SNPRecord(
        rsid="rs12722", gene="COL5A1", chrom="9",
        pos_grch37=137_721_567, pos_grch38=134_829_415,
        ancestral="T", derived="C",
        description="COL5A1 3'-UTR — CC = stiffer tendon, increased "
                    "tendinopathy risk; TT = protective running economy.",
        source="Mokone 2006",
    ),
    SNPRecord(
        rsid="rs41423247", gene="NR3C1", chrom="5",
        pos_grch37=142_778_575, pos_grch38=143_398_988,
        ancestral="C", derived="G",
        description="NR3C1 BclI — G allele = higher cortisol sensitivity "
                    "to glucocorticoids.",
        source="van Rossum 2003",
        aliases=("BclI",),
    ),
    SNPRecord(
        rsid="rs5751876", gene="ADORA2A", chrom="22",
        pos_grch37=24_827_015, pos_grch38=24_424_828,
        ancestral="C", derived="T",
        description="ADORA2A — T allele = caffeine-induced anxiety "
                    "tendency (adenosine A2A receptor).",
        source="Childs 2008",
    ),
    SNPRecord(
        rsid="rs225014", gene="DIO2", chrom="14",
        pos_grch37=80_682_044, pos_grch38=80_215_701,
        ancestral="T", derived="C",
        description="DIO2 Thr92Ala — C (Ala) = reduced peripheral T4→T3 "
                    "conversion; some patients prefer T3 supplementation.",
        source="Panicker 2009",
        aliases=("Thr92Ala",),
    ),
    SNPRecord(
        rsid="rs2802292", gene="FOXO3", chrom="6",
        pos_grch37=108_587_315, pos_grch38=108_265_853,
        ancestral="T", derived="G",
        description="FOXO3 — G allele = longevity-associated; over-represented "
                    "in centenarian cohorts across populations.",
        source="Willcox 2008 (PNAS)",
    ),
    SNPRecord(
        rsid="rs1050450", gene="GPX1", chrom="3",
        pos_grch37=49_394_834, pos_grch38=49_357_401,
        ancestral="C", derived="T",
        description="GPX1 Pro198Leu — T allele = reduced glutathione "
                    "peroxidase activity, increased oxidative susceptibility.",
        source="Forsberg 2000",
        aliases=("Pro198Leu",),
    ),
    SNPRecord(
        rsid="rs2736100", gene="TERT", chrom="5",
        pos_grch37=1_286_516, pos_grch38=1_286_401,
        ancestral="C", derived="A",
        description="TERT — A allele linked to longer telomere length "
                    "(GWAS leukocyte-telomere length).",
        source="Codd 2013",
    ),

    # ── Gut health panel (gut_health.py) ─────────────────────────────────
    SNPRecord(
        rsid="rs4988235", gene="MCM6", chrom="2",
        pos_grch37=136_608_646, pos_grch38=135_851_076,
        ancestral="G", derived="A",
        description="MCM6/LCT −13910 — derived A (+ strand) is the lactase-"
                    "persistence allele. Often reported as −13910C>T on the "
                    "− strand; the registry is canonical + strand G/A.",
        source="Enattah 2002 (Nat Genet)",
        aliases=("-13910C>T", "LCT -13910"),
    ),
    SNPRecord(
        rsid="rs10156191", gene="AOC1", chrom="7",
        pos_grch37=150_553_605, pos_grch38=150_856_517,
        ancestral="C", derived="T",
        description="AOC1 (DAO) Thr16Met — derived T lowers diamine-oxidase "
                    "activity, implicated in dietary histamine intolerance.",
        source="Maintz 2011",
        aliases=("Thr16Met",),
    ),
    SNPRecord(
        rsid="rs2066844", gene="NOD2", chrom="16",
        pos_grch37=50_745_926, pos_grch38=50_712_015,
        ancestral="C", derived="T",
        description="NOD2 R702W — derived T is a Crohn's-disease "
                    "susceptibility allele.",
        source="Hugot 2001 (Nature)",
        aliases=("R702W",),
    ),
    SNPRecord(
        rsid="rs11209026", gene="IL23R", chrom="1",
        pos_grch37=67_705_958, pos_grch38=67_240_275,
        ancestral="G", derived="A",
        description="IL23R R381Q — derived A is protective against "
                    "inflammatory bowel disease.",
        source="Duerr 2006 (Science)",
        aliases=("R381Q",),
    ),

    # ── Metals & oxidative stress panel (metal_oxidative.py) ─────────────
    SNPRecord(
        rsid="rs34637584", gene="LRRK2", chrom="12",
        pos_grch37=40_734_202, pos_grch38=40_340_400,
        ancestral="G", derived="A",
        description="LRRK2 G2019S — derived A is the most common Mendelian "
                    "risk variant for Parkinson's disease.",
        source="Healy 2008 (Lancet Neurol)",
        aliases=("G2019S",),
        clinical_significance="Pathogenic",
    ),
    SNPRecord(
        rsid="rs1001179", gene="CAT", chrom="11",
        pos_grch37=34_460_231, pos_grch38=34_438_684,
        ancestral="C", derived="T",
        description="CAT −262C>T promoter — derived T alters catalase "
                    "expression and oxidative-stress handling.",
        source="Forsberg 2001",
        aliases=("-262C>T",),
    ),
    SNPRecord(
        rsid="rs1061472", gene="ATP7B", chrom="13",
        pos_grch37=52_524_488, pos_grch38=51_950_352,
        ancestral="T", derived="C",
        description="ATP7B K832R — copper-transport variant. Canonical "
                    "+ strand T/C; commonly reported A/G on the − strand "
                    "(− strand G = + strand C).",
        source="dbSNP; ATP7B copper metabolism",
        aliases=("K832R",),
    ),
    SNPRecord(
        rsid="rs8052394", gene="MT1A", chrom="16",
        pos_grch37=56_673_828, pos_grch38=56_639_916,
        ancestral="A", derived="G",
        description="MT1A metallothionein variant — derived G associated "
                    "with altered zinc/cadmium handling.",
        source="dbSNP; metallothionein literature",
    ),
    SNPRecord(
        rsid="rs28366003", gene="MT2A", chrom="16",
        pos_grch37=56_642_491, pos_grch38=56_608_579,
        ancestral="A", derived="G",
        description="MT2A −5A>G promoter metallothionein variant — derived G "
                    "associated with altered metal-binding capacity.",
        source="dbSNP; metallothionein literature",
        aliases=("-5A>G",),
    ),
    SNPRecord(
        rsid="rs13107325", gene="SLC39A8", chrom="4",
        pos_grch37=103_188_709, pos_grch38=102_267_552,
        ancestral="C", derived="T",
        description="SLC39A8 (ZIP8) A391T — derived T (Thr391) lowers "
                    "manganese/zinc transport; highly pleiotropic (blood "
                    "pressure, lipids, neuropsychiatric, IBD). Canonical "
                    "+ strand C/T; arrays often report on the − strand as "
                    "G/A (− strand A = + strand T = effect allele).",
        source="dbSNP; GWAS of Mn metabolism & pleiotropy",
        aliases=("A391T",),
    ),
    SNPRecord(
        rsid="rs896378", gene="SLC39A14", chrom="8",
        pos_grch37=22_262_321, pos_grch38=22_404_808,
        ancestral="T", derived="C",
        description="SLC39A14 (ZIP14) common coding variant — ZIP14 is a "
                    "manganese/zinc/iron importer; benign missense near "
                    "Leu33. Research-grade. Reference + strand T; common "
                    "allele C. (Ancestral state ambiguous; reports + strand "
                    "reference vs common.)",
        source="dbSNP rs896378",
    ),
    SNPRecord(
        rsid="rs1893590", gene="ABCG1", chrom="21",
        pos_grch37=43_619_595, pos_grch38=42_199_485,
        ancestral="A", derived="C",
        description="ABCG1 −204A>C promoter variant — ABCG1 mediates "
                    "cholesterol/sterol efflux; derived C associated with "
                    "lower HDL-C in candidate-gene studies. Research-grade.",
        source="dbSNP; PMID 25398214 (ABCG1/HDL)",
        aliases=("-204A>C",),
    ),
    SNPRecord(
        rsid="rs4147565", gene="GSTM1", chrom="1",
        pos_grch37=110_231_777, pos_grch38=109_689_155,
        ancestral="G", derived="A",
        description="GSTM1 within-gene marker (glutathione-S-transferase Mu "
                    "1). The functional GSTM1-null state is a whole-gene "
                    "DELETION (CNV), not this SNP — null shows up as a "
                    "no-call/hemizygous genotype, so this is only a "
                    "research-grade proxy; PCR/CNV assay is the gold standard.",
        source="dbSNP; GST copy-number literature (PMC6118300)",
    ),
    SNPRecord(
        rsid="rs4630", gene="GSTT1", chrom="22",
        pos_grch37=24_376_322, pos_grch38=None,
        ancestral="G", derived="A",
        description="GSTT1 within-gene coding marker (glutathione-S-transferase "
                    "Theta 1). The functional GSTT1-null state is a whole-gene "
                    "DELETION (CNV), not this SNP — null shows up as a "
                    "no-call/hemizygous genotype, so this is only a "
                    "research-grade proxy; PCR/CNV assay is the gold standard. "
                    "GRCh38 omitted: this CNV-prone region has no clean "
                    "primary-chromosome hg38 placement.",
        source="dbSNP; GST copy-number literature (PMC6118300)",
    ),

    # ── Detoxification & environmental resilience (smoke / PAH / metals) ──
    # Added for detox.py: Phase I bioactivation (CYP1A1/1B1/AHR), Phase II
    # conjugation (EPHX1, NAT2, NQO1), and heavy-metal handling (ALAD, AS3MT,
    # PON1). GRCh37 positions taken from 23andMe/AncestryDNA + dbSNP; a few
    # GRCh38 positions omitted where not cleanly verified.
    SNPRecord(
        rsid="rs1048943", gene="CYP1A1", chrom="15",
        pos_grch37=75_012_985, pos_grch38=74_720_644,
        ancestral="A", derived="G",
        description="CYP1A1 Ile462Val (m2) — Val462 raises inducibility/activity "
                    "toward polycyclic aromatic hydrocarbons from smoke; more "
                    "Phase I bioactivation.",
        source="dbSNP; CYP1A1 smoking/PAH literature",
        aliases=("Ile462Val", "CYP1A1*2C", "m2"),
    ),
    SNPRecord(
        rsid="rs4646903", gene="CYP1A1", chrom="15",
        pos_grch37=75_011_641, pos_grch38=74_719_300,
        ancestral="T", derived="C",
        description="CYP1A1 m1 (3' UTR MspI) — variant allele associated with "
                    "higher CYP1A1 inducibility, especially in smokers.",
        source="dbSNP; CYP1A1 MspI literature",
        aliases=("CYP1A1*2A", "m1", "MspI"),
    ),
    SNPRecord(
        rsid="rs1056836", gene="CYP1B1", chrom="2",
        pos_grch37=38_298_203, pos_grch38=38_071_060,
        ancestral="C", derived="G",
        description="CYP1B1 Leu432Val — Val432 increases activity toward PAHs "
                    "and estrogen substrates; a Phase I activation route.",
        source="dbSNP; CYP1B1 literature",
        aliases=("Leu432Val", "L432V"),
    ),
    SNPRecord(
        rsid="rs2066853", gene="AHR", chrom="7",
        pos_grch37=17_379_110, pos_grch38=17_339_484,
        ancestral="G", derived="A",
        description="AHR Arg554Lys — modulates aryl-hydrocarbon-receptor "
                    "signalling that induces the CYP1 enzymes on smoke/dioxin "
                    "exposure.",
        source="dbSNP; AhR signalling literature",
        aliases=("Arg554Lys", "R554K"),
    ),
    SNPRecord(
        rsid="rs1051740", gene="EPHX1", chrom="1",
        pos_grch37=226_019_633, pos_grch38=225_831_932,
        ancestral="T", derived="C",
        description="EPHX1 Tyr113His (exon 3) — His113 (C) lowers microsomal "
                    "epoxide-hydrolase activity ('slow'), reducing clearance of "
                    "PAH epoxides from smoke.",
        source="dbSNP; EPHX1 activity literature",
        aliases=("Tyr113His", "Y113H"),
    ),
    SNPRecord(
        rsid="rs2234922", gene="EPHX1", chrom="1",
        pos_grch37=226_026_406, pos_grch38=225_838_705,
        ancestral="A", derived="G",
        description="EPHX1 His139Arg (exon 4) — Arg139 (G) raises epoxide-"
                    "hydrolase activity ('fast'); combines with Y113H to predict "
                    "overall EPHX1 activity.",
        source="dbSNP; EPHX1 activity literature",
        aliases=("His139Arg", "H139R"),
    ),
    SNPRecord(
        rsid="rs1801280", gene="NAT2", chrom="8",
        pos_grch37=18_257_854, pos_grch38=18_400_343,
        ancestral="T", derived="C",
        description="NAT2*5 (I114T, T341C) — slow-acetylator allele; slower "
                    "clearance of aromatic amines from tobacco/combustion smoke.",
        source="dbSNP; NAT2 acetylator literature",
        aliases=("NAT2*5", "I114T"),
    ),
    SNPRecord(
        rsid="rs1799930", gene="NAT2", chrom="8",
        pos_grch37=18_258_103, pos_grch38=18_400_592,
        ancestral="G", derived="A",
        description="NAT2*6 (R197Q, G590A) — slow-acetylator allele for aromatic "
                    "amine detoxification.",
        source="dbSNP; NAT2 acetylator literature",
        aliases=("NAT2*6", "R197Q"),
    ),
    SNPRecord(
        rsid="rs1799931", gene="NAT2", chrom="8",
        pos_grch37=18_258_370, pos_grch38=18_400_859,
        ancestral="G", derived="A",
        description="NAT2*7 (G286E, G857A) — slow-acetylator allele contributing "
                    "to the NAT2 slow phenotype.",
        source="dbSNP; NAT2 acetylator literature",
        aliases=("NAT2*7", "G286E"),
    ),
    SNPRecord(
        rsid="rs1800566", gene="NQO1", chrom="16",
        pos_grch37=69_745_145, pos_grch38=69_711_242,
        ancestral="C", derived="T",
        description="NQO1 Pro187Ser (C609T) — Ser187 (T) greatly reduces NQO1 "
                    "activity, impairing quinone detoxification (benzene/smoke) "
                    "and antioxidant recycling.",
        source="dbSNP; NQO1 P187S literature",
        aliases=("Pro187Ser", "P187S", "NQO1*2"),
    ),
    SNPRecord(
        rsid="rs2071746", gene="HMOX1", chrom="22",
        pos_grch37=35_776_290, pos_grch38=35_380_301,
        ancestral="A", derived="T",
        description="HMOX1 -413A>T promoter — modulates expression of heme "
                    "oxygenase-1, a cytoprotective enzyme induced by particulate "
                    "/ oxidative stress.",
        source="dbSNP; HMOX1 promoter literature",
        aliases=("-413A>T",),
    ),
    SNPRecord(
        rsid="rs662", gene="PON1", chrom="7",
        pos_grch37=94_937_446, pos_grch38=95_308_134,
        ancestral="A", derived="G",
        description="PON1 Q192R — Arg192 (G) shifts paraoxonase-1 substrate "
                    "specificity (organophosphate hydrolysis vs oxidised-lipid "
                    "clearance).",
        source="dbSNP; PON1 Q192R literature",
        aliases=("Q192R", "Gln192Arg"),
    ),
    SNPRecord(
        rsid="rs1800435", gene="ALAD", chrom="9",
        pos_grch37=116_153_891, pos_grch38=113_391_612,
        ancestral="C", derived="G",
        description="ALAD K59N (ALAD2 allele) — alters lead binding in blood; "
                    "associated in some cohorts with higher blood-lead retention "
                    "for a given exposure.",
        source="dbSNP; ALAD lead-kinetics literature",
        aliases=("K59N", "ALAD2"),
    ),
    SNPRecord(
        rsid="rs11191439", gene="AS3MT", chrom="10",
        pos_grch37=104_638_723, pos_grch38=102_878_966,
        ancestral="T", derived="C",
        description="AS3MT Met287Thr — influences arsenite-methyltransferase "
                    "activity and the arsenic methylation/excretion profile.",
        source="dbSNP; AS3MT arsenic-methylation literature",
        aliases=("Met287Thr", "M287T"),
    ),
]


# Public lookup table — built once at import time.
SNPS: dict[str, SNPRecord] = {r.rsid: r for r in _RECORDS}


# ── Lookup API ───────────────────────────────────────────────────────────────

def get(rsid: str) -> SNPRecord | None:
    """Return the canonical record for an rsID, or ``None`` if not registered."""
    return SNPS.get(rsid)


def require(rsid: str) -> SNPRecord:
    """Like :func:`get` but raises ``KeyError`` if the rsID is not registered.
    Use in consumers where missing metadata is a programming error, not a
    runtime condition."""
    rec = SNPS.get(rsid)
    if rec is None:
        raise KeyError(f"rsID {rsid!r} not in snp_registry. Add it to _RECORDS.")
    return rec


def by_gene(gene: str) -> list[SNPRecord]:
    """All registered SNPs for one gene."""
    return [r for r in SNPS.values() if r.gene == gene]


def lookup_by_pos(chrom: str, pos: int, build: str = "grch37") -> SNPRecord | None:
    """Fallback when a chip reports a variant under a non-standard rsID label
    — match by chromosomal position instead."""
    pos_attr = f"pos_{build.lower()}"
    for rec in SNPS.values():
        if rec.chrom == chrom and getattr(rec, pos_attr, None) == pos:
            return rec
    return None


# ── Migration audit ──────────────────────────────────────────────────────────

STALE_AFTER_DAYS = 365


def audit_registry(today: _dt.date | None = None) -> dict[str, Any]:
    """
    Pre-flight check for consumer modules. Returns a summary plus a list of
    records that are old enough to deserve human re-verification against the
    cited source.
    """
    today = today or _dt.date.today()
    stale: list[str] = []
    for rec in SNPS.values():
        verified = _dt.date.fromisoformat(rec.last_verified)
        if (today - verified).days > STALE_AFTER_DAYS:
            stale.append(rec.rsid)
    return {
        "n_records": len(SNPS),
        "n_genes": len({r.gene for r in SNPS.values()}),
        "n_stale": len(stale),
        "stale_rsids": stale,
        "checked": today.isoformat(),
    }


def dump_json(path: Path) -> None:
    """Serialise the registry to disk — useful for downstream non-Python tools
    (e.g. the Mermaid registry diagram in the docs, or external validation)."""
    Path(path).write_text(
        json.dumps([asdict(r) for r in SNPS.values()], indent=2)
    )
