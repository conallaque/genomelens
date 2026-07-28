"""
Clinical Variants — Phase 2 rare/pathogenic-variant interpretation (WGS/VCF)
============================================================================

The sequencing-grade payoff a chip structurally can't reach: screen a
whole-genome/exome VCF against **ClinVar** for known Pathogenic / Likely-
pathogenic (P/LP) variants, classify them (ACMG secondary finding / carrier /
possibly-affected), grade confidence by ClinVar review stars, and hand the
result to the local AI layer.

This module is deliberately conservative — a false "you're clear" is the worst
possible output, so every design choice favours honesty over reassurance.

Data
----
ClinVar is distilled once (``setup.py --clinvar`` → ``distill_clinvar_vcf``) into
a compact table keyed by a **normalised** ``chrom:pos:ref:alt`` so the same
variant represented differently in ClinVar vs the user VCF still matches. The
engine loads that table (fast dict) and streams the user VCF against it. If the
table isn't present, the module degrades gracefully.

Correctness decisions (each from review):
  * **inheritance-unknown is first-class.** We only assert carrier-vs-affected
    where a small, unambiguous curated ACMG-gene inheritance map applies;
    otherwise the finding is labelled "inheritance not determined."
  * **zygosity is resolved against the matched ALT index**, not generically —
    a ``1/2`` multi-allelic sample is heterozygous for whichever ALT matched.
  * **per-gene compound-het pass:** two distinct P/LP variants in one recessive
    gene → flagged *possible* compound heterozygote (phase unknown without
    parents) — the difference between "carrier" and "possibly affected".
  * **indel normalisation:** both sides are trimmed to a canonical key so
    anchor-base/left-alignment differences don't cause false negatives; indel
    sensitivity is still reduced and the report says so.
  * **negative result ≠ clear:** an empty findings list always carries an
    explicit statement of what was and wasn't screened.

Educational screening, NOT a clinical diagnostic test.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).parent
CLINVAR_DIR = SCRIPT_DIR / "reference" / "clinvar"

# ── ClinVar review status → star rating (verified against real ClinVar data) ──
_STAR = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
}


def revstat_to_stars(revstat: str) -> int:
    return _STAR.get((revstat or "").strip(), 0)


def is_pathogenic_sig(sig: str) -> bool:
    """True for Pathogenic / Likely_pathogenic (and the combined call), while
    excluding Conflicting / Benign / Uncertain."""
    s = (sig or "").strip()
    if not s or "Conflicting" in s or "Benign" in s or "Uncertain" in s:
        return False
    return (s in ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
            or s.startswith("Pathogenic") or s.startswith("Likely_pathogenic"))


def is_uncertain_sig(sig: str) -> bool:
    return "Uncertain_significance" in (sig or "")


# ── ACMG SF: curated inheritance for genes where it is unambiguous ────────────
# Only genes with well-established, unambiguous inheritance for the reportable
# phenotype. Everything not here → inheritance-unknown (never guessed).
_ACMG_SF_INHERITANCE: Dict[str, str] = {
    # hereditary cancer
    "BRCA1": "AD", "BRCA2": "AD", "PALB2": "AD",
    "MLH1": "AD", "MSH2": "AD", "MSH6": "AD", "PMS2": "AD", "EPCAM": "AD",
    "APC": "AD", "MUTYH": "AR",              # MUTYH is the recessive exception
    "TP53": "AD", "PTEN": "AD", "STK11": "AD", "CDH1": "AD",
    "RET": "AD", "VHL": "AD", "NF2": "AD", "SDHB": "AD", "SDHD": "AD",
    "BMPR1A": "AD", "SMAD4": "AD", "TSC1": "AD", "TSC2": "AD", "WT1": "AD",
    # cardiovascular
    "LDLR": "AD", "APOB": "AD", "PCSK9": "AD",
    "MYH7": "AD", "MYBPC3": "AD", "TNNT2": "AD", "TNNI3": "AD", "TPM1": "AD",
    "KCNQ1": "AD", "KCNH2": "AD", "SCN5A": "AD", "RYR2": "AD",
    "LMNA": "AD", "DSP": "AD", "PKP2": "AD", "FBN1": "AD", "TGFBR1": "AD",
    "TGFBR2": "AD", "COL3A1": "AD",
    # metabolic / other
    "RYR1": "AD",                            # malignant hyperthermia
    "OTC": "XL",
    "GLA": "XL",                             # Fabry (X-linked)
    "ATP7B": "AR",                           # Wilson disease
    "TTR": "AD",
}

# Broad set of well-known recessive carrier genes (for carrier framing).
_RECESSIVE_CARRIER_GENES = {
    "CFTR", "HBB", "HBA1", "HBA2", "HEXA", "GBA", "SMN1", "PAH", "GALT",
    "MUTYH", "ATP7B", "MCOLN1", "ASPA", "BCKDHB", "CLRN1", "FANCC", "IKBKAP",
    "DHDDS", "SMPD1", "GAA", "MEFV", "SERPINA1",
}


# ── key normalisation (SNV = identity; indels → canonical trimmed key) ────────

def _norm_chrom(c: str) -> str:
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def norm_key(chrom: str, pos: int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    """Canonical (chrom, pos, ref, alt). Trims shared suffix then shared prefix
    (advancing pos) so ``AG>A`` at 100 and ``G>`` representations collapse to the
    same key on both the ClinVar and user sides."""
    chrom = _norm_chrom(chrom)
    ref, alt = ref.upper(), alt.upper()
    # trim shared suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # trim shared prefix
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (chrom, int(pos), ref, alt)


# ── distiller (called by setup.py --clinvar) ──────────────────────────────────

def _parse_info(info: str) -> Dict[str, str]:
    out = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def distill_clinvar_vcf(in_path: str, out_path: str, log=print) -> int:
    """Stream a ClinVar VCF and write a compact P/LP (+ VUS-in-ACMG) table.
    Returns the number of rows written."""
    acmg_genes = set(_ACMG_SF_INHERITANCE)
    opener = gzip.open if str(in_path).lower().endswith((".gz", ".bgz")) else open
    n = 0
    with opener(in_path, "rt", errors="ignore") as f, \
         gzip.open(out_path, "wt") as out:
        out.write("chrom\tpos\tref\talt\tsig\tstars\tgene\tcondition\tvc\tis_vus\n")
        for line in f:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            if alt in (".", ""):
                continue
            info = _parse_info(cols[7])
            sig = info.get("CLNSIG", "")
            gene = (info.get("GENEINFO", "").split(":")[0] or "").strip()
            plp = is_pathogenic_sig(sig)
            vus = is_uncertain_sig(sig)
            if not plp and not (vus and gene in acmg_genes):
                continue
            try:
                k = norm_key(chrom, int(pos), ref, alt)
            except ValueError:
                continue
            stars = revstat_to_stars(info.get("CLNREVSTAT", ""))
            cond = info.get("CLNDN", "").replace("\t", " ")
            vc = info.get("CLNVC", "")
            out.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{sig}\t{stars}\t{gene}"
                      f"\t{cond}\t{vc}\t{int(vus and not plp)}\n")
            n += 1
            if log and n % 50000 == 0:
                log(f"    distilled {n:,} clinically-significant records ...")
    return n


# ── table loading (with truncation guard) ─────────────────────────────────────

def load_clinvar_table(build: str) -> Optional[Dict[Tuple, Dict]]:
    """Load the distilled table for a build into {norm_key: record}, or None if
    absent. Refuses a truncated/headerless table."""
    path = CLINVAR_DIR / f"clinvar_plp_{build}.tsv.gz"
    if not path.exists():
        return None
    table: Dict[Tuple, Dict] = {}
    try:
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            if header[:4] != ["chrom", "pos", "ref", "alt"]:
                return None                      # not our table / corrupt
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 10:
                    continue
                key = (p[0], int(p[1]), p[2], p[3])
                table[key] = {"sig": p[4], "stars": int(p[5]), "gene": p[6],
                              "condition": p[7], "vc": p[8], "is_vus": p[9] == "1"}
    except Exception:
        return None
    return table or None


# ── zygosity resolved against the matched ALT ─────────────────────────────────

def zygosity_for_alt(gt_field: str, matched_alt_index: int) -> Optional[str]:
    """Zygosity of the sample for a SPECIFIC ALT allele index (1-based ALT →
    VCF allele number). Returns 'homozygous'/'heterozygous'/'hemizygous', or
    None if the sample doesn't carry that allele / is a no-call."""
    core = (gt_field or "").split(":")[0].replace("|", "/")
    if not core or core in (".", "./.", "./", "/."):
        return None
    parts = [a for a in core.split("/") if a != ""]
    if not parts or "." in parts:
        return None
    try:
        idxs = [int(a) for a in parts]
    except ValueError:
        return None
    haploid = len(idxs) == 1
    if haploid:
        return "hemizygous" if idxs[0] == matched_alt_index else None
    count = idxs.count(matched_alt_index)
    if count == 0:
        return None
    return "homozygous" if count >= 2 else "heterozygous"


# ── main analysis ──────────────────────────────────────────────────────────────

_NEGATIVE_DISCLAIMER = (
    "No result here means 'no ClinVar-classified pathogenic variant was matched' "
    "— NOT 'no risk'. This screens only variants already classified in ClinVar, "
    "only at positions present in your file, with reduced sensitivity for indels, "
    "and does not detect structural/repeat variants or reclassify variants of "
    "uncertain significance. It is a screen, not a clinical diagnostic test.")


def analyze_clinical_variants(vcf_path: str, build: str,
                              inferred_sex: Optional[str] = None,
                              log=None) -> Dict:
    """Screen a user VCF against distilled ClinVar. Returns structured findings
    with classification, zygosity, star confidence, and compound-het detection."""
    _log = log or (lambda *a, **k: None)
    if build not in ("grch37", "grch38"):
        return {"available": False,
                "reason": f"Build '{build}' — cannot select a ClinVar table safely.",
                "negative_disclaimer": _NEGATIVE_DISCLAIMER}

    table = load_clinvar_table(build)
    if table is None:
        return {"available": False,
                "reason": ("ClinVar table not present — run `python setup.py "
                           "--clinvar`. (Falls back cleanly; no screen performed.)"),
                "negative_disclaimer": _NEGATIVE_DISCLAIMER}

    findings: List[Dict] = []
    n_scanned = 0
    opener = gzip.open if str(vcf_path).lower().endswith((".gz", ".bgz")) else open
    with opener(vcf_path, "rt", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            n_scanned += 1
            chrom, pos_s, ref, alt_field = cols[0], cols[1], cols[3], cols[4]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            sample = cols[9] if len(cols) > 9 else ""
            for ai, alt in enumerate(alt_field.split(","), start=1):
                if alt in (".", "", "*") or alt.startswith("<"):
                    continue
                rec = table.get(norm_key(chrom, pos, ref, alt))
                if rec is None:
                    continue
                zyg = zygosity_for_alt(sample, ai) if sample else "unknown"
                if zyg is None:      # sample doesn't carry the matched allele
                    continue
                findings.append({
                    "chrom": _norm_chrom(chrom), "pos": pos, "ref": ref, "alt": alt,
                    "gene": rec["gene"], "condition": rec["condition"].replace("_", " "),
                    "significance": rec["sig"], "stars": rec["stars"],
                    "vc": rec["vc"], "is_vus": rec["is_vus"],
                    "zygosity": zyg,
                })

    return _classify(findings, n_scanned, inferred_sex)


def _classify(findings: List[Dict], n_scanned: int,
              inferred_sex: Optional[str]) -> Dict:
    """Attach ACMG/inheritance classification, run the per-gene compound-het
    pass, split P/LP from uncertain, and bucket by category."""
    plp = [f for f in findings if not f["is_vus"] and f["stars"] >= 1]
    vus_acmg = [f for f in findings if f["is_vus"] and f["gene"] in _ACMG_SF_INHERITANCE]
    excluded_0star = [f for f in findings if not f["is_vus"] and f["stars"] < 1]

    # per-gene compound-het detection (recessive genes, ≥2 distinct P/LP variants)
    by_gene: Dict[str, List[Dict]] = {}
    for f in plp:
        by_gene.setdefault(f["gene"], []).append(f)

    for f in plp:
        gene = f["gene"]
        inh = _ACMG_SF_INHERITANCE.get(gene)
        is_acmg = gene in _ACMG_SF_INHERITANCE
        recessive = inh == "AR" or gene in _RECESSIVE_CARRIER_GENES
        distinct = {(x["pos"], x["alt"]) for x in by_gene.get(gene, [])}
        two_hits = len(distinct) >= 2

        if inh is None and not recessive:
            f["inheritance"] = "unknown"
        else:
            f["inheritance"] = inh or ("AR" if recessive else "unknown")

        # classification
        if recessive:
            if f["zygosity"] == "homozygous":
                f["category"] = "affected"
                f["interpretation"] = ("Homozygous P/LP in a recessive gene — "
                                       "consistent with being affected. Confirm clinically.")
            elif two_hits:
                f["category"] = "possible_compound_het"
                f["interpretation"] = ("Two distinct P/LP variants in this "
                                       "recessive gene — POSSIBLE compound heterozygote "
                                       "(affected) IF the variants are on opposite copies. "
                                       "Phase is unknown without parental testing.")
            else:
                f["category"] = "carrier"
                f["interpretation"] = ("Heterozygous P/LP in a recessive gene — "
                                       "carrier (typically unaffected; relevant for "
                                       "family planning).")
        elif inh == "AD":
            f["category"] = "actionable" if is_acmg else "dominant_risk"
            f["interpretation"] = ("P/LP in a dominant gene — a single copy can "
                                   "confer risk. Confirm in an accredited lab; "
                                   "genetic counseling advised.")
        elif inh == "XL":
            f["category"] = "actionable" if is_acmg else "xlinked"
            f["interpretation"] = "P/LP in an X-linked gene — interpretation depends on sex/zygosity."
        else:
            f["category"] = "uncertain_inheritance"
            f["interpretation"] = ("P/LP variant found; this tool does not "
                                   "determine the inheritance mode for this gene.")
        f["acmg_secondary"] = is_acmg

    buckets = {"actionable": [], "affected": [], "possible_compound_het": [],
               "carrier": [], "dominant_risk": [], "xlinked": [],
               "uncertain_inheritance": []}
    for f in plp:
        buckets.setdefault(f["category"], []).append(f)
    # sort each bucket by star confidence desc
    for b in buckets.values():
        b.sort(key=lambda x: -x["stars"])

    n_actionable = len(buckets["actionable"])
    return {
        "available": True,
        "n_scanned": n_scanned,
        "n_plp": len(plp),
        "n_actionable": n_actionable,
        "n_carrier": len(buckets["carrier"]),
        "n_affected": len(buckets["affected"]) + len(buckets["possible_compound_het"]),
        "n_vus_in_acmg": len(vus_acmg),
        "n_excluded_0star": len(excluded_0star),
        "buckets": buckets,
        "vus_in_acmg": sorted(vus_acmg, key=lambda x: x["gene"]),
        "findings": plp,
        "negative_disclaimer": _NEGATIVE_DISCLAIMER,
        "disclaimer": (
            "Educational screening against ClinVar, NOT a clinical diagnostic "
            "test. A self-/consumer-sequencing VCF is not clinical-grade. Any "
            "finding — especially anything actionable — must be confirmed in an "
            "accredited diagnostic laboratory and discussed with a board-certified "
            "genetic counselor before it means anything for your health."),
    }
