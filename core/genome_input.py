"""
Genome input layer — accept EITHER a consumer chip file OR a whole-genome VCF
=============================================================================

The pipeline was built chip-first: every module looks a variant up by rsID on a
DataFrame produced by the ``snps`` library. That library *can* read a VCF, but
with two consequences that matter for whole-genome sequencing (WGS):

  1. It **skips every variant that has no rsID** (the ID column is ``.``) — which
     is the entire rare/novel WGS payoff, and also any *curated* registry
     variant that happens to be position-only in the file.
  2. It skips indels and loads everything into memory.

This module adds a thin, dependency-free layer so the program cleanly runs on
either input:

  * ``looks_like_vcf(path)`` — detect input type.
  * ``vcf_gt_to_genotype(gt, ref, alt)`` — convert a VCF ``GT`` + ``REF``/``ALT``
    into the two-character, strand-as-reported, sorted genotype string the
    module layer expects. Conversion rules (tested one case each):

        VCF                              → emitted        reason
        ------------------------------------------------------------------
        0/1  REF=A ALT=G                 → "AG"           het SNV (sorted)
        1|1  REF=A ALT=G   (phased)      → "GG"           hom ALT
        0/0  REF=A ALT=G                 → "AA"           hom REF
        1    REF=A ALT=G   (haploid X/Y) → "GG"           haploid → homozygous
        1/2  REF=A ALT=G,T (multiallelic)→ "GT"           both ALTs (sorted)
        ./.  or partial "."              → None           no-call → skip
        REF=AT ALT=A       (indel)       → None           not a 2-char SNV → skip
        ALT=<DEL> / symbolic             → None           structural → skip
        allele index out of range        → None           malformed → skip

  * ``enrich_and_profile_vcf(...)`` — a SINGLE streaming pass (stdlib ``gzip``,
    memory-safe on a 100 GB file) that (a) **back-fills** curated registry
    positions the rsID reader dropped, gated on the *detected build* — and
    (b) profiles what the VCF contains that the pipeline can't yet interpret
    (total variants, how many lack rsIDs, how many fall in ACMG actionable
    genes), quantifying the value of the not-yet-built Phase-2 rare-variant
    engine instead of asserting it.

**Build safety:** back-fill selects ``pos_grch37`` vs ``pos_grch38`` from the
detected build and *refuses* (logs + skips) when the build is ``mixed`` or
``unknown`` — a silently wrong genotype at a curated clinical variant is worse
than an absent one.

**Strand:** VCF REF/ALT are plus-strand relative to the reference build; the
registry's ancestral/derived are documented plus-strand. We emit the observed
genotype verbatim (no flip); downstream compares it to plus-strand registry
alleles, so alignment holds.
"""

from __future__ import annotations

import gzip
import re

import pandas as pd

from . import snp_registry as _reg

# ── Input-type detection ──────────────────────────────────────────────────────

def looks_like_vcf(path: str) -> bool:
    """True if the file is a VCF (by extension, else by sniffing a ``##fileformat=VCF``
    header in the first lines — handles gzip)."""
    p = str(path).lower()
    if p.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        return True
    try:
        opener = gzip.open if p.endswith((".gz", ".bgz")) else open
        with opener(path, "rt", errors="ignore") as f:
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
                if line.startswith("##fileformat=VCF") or (
                        line.startswith("#CHROM") and "POS" in line):
                    return True
    except Exception:
        pass
    return False


# ── Build detection from the VCF header (no rsIDs required) ────────────────────
# Distinctive per-contig lengths: (GRCh37, GRCh38). A handful of unambiguous
# chromosomes is enough to call the build straight from the header — essential for
# whole-genome callsets (e.g. GIAB benchmarks) whose variants have no rsIDs, so the
# probe-position detector in provenance.py returns "unknown".
_CONTIG_BUILD_LEN: dict[str, tuple[int, int]] = {
    "1": (249250621, 248956422),
    "2": (243199373, 242193529),
    "3": (198022430, 198295559),
    "7": (159138663, 159345973),
    "X": (155270560, 156040895),
}


def detect_build_from_vcf_header(path: str, max_header_lines: int = 5000) -> str:
    """Detect ``grch37`` vs ``grch38`` from a VCF's ``##contig`` lengths (falling
    back to a ``##reference`` hint). Returns ``"grch37"``, ``"grch38"``, or
    ``"unknown"``. Does not read variant rows, so it is O(header)."""
    p = str(path)
    opener = gzip.open if p.lower().endswith((".gz", ".bgz")) else open
    votes = {"grch37": 0, "grch38": 0}
    ref_hint: str | None = None
    try:
        with opener(p, "rt", errors="ignore") as f:
            for i, line in enumerate(f):
                if not line.startswith("#"):
                    break                      # past the header
                if i > max_header_lines:
                    break
                low = line.lower()
                if low.startswith("##reference"):
                    if any(t in low for t in ("grch38", "hg38", "b38")):
                        ref_hint = "grch38"
                    elif any(t in low for t in ("grch37", "hg19", "b37")):
                        ref_hint = "grch37"
                elif line.startswith("##contig"):
                    m_id = re.search(r"ID=([^,>]+)", line)
                    m_len = re.search(r"length=(\d+)", line)
                    if m_id and m_len:
                        pair = _CONTIG_BUILD_LEN.get(_norm_chrom(m_id.group(1)))
                        if pair:
                            length = int(m_len.group(1))
                            if length == pair[0]:
                                votes["grch37"] += 1
                            elif length == pair[1]:
                                votes["grch38"] += 1
    except Exception:
        return "unknown"
    if votes["grch37"] or votes["grch38"]:
        return "grch37" if votes["grch37"] >= votes["grch38"] else "grch38"
    return ref_hint or "unknown"


# ── Genotype conversion ───────────────────────────────────────────────────────

_ACGT = set("ACGT")


def vcf_gt_to_genotype(gt: str, ref: str, alt: str) -> str | None:
    """Convert a VCF sample ``GT`` (+ REF/ALT) to a sorted 2-char SNV genotype,
    or None to skip (no-call, indel, symbolic, malformed). See module docstring
    for the full rules table."""
    if not gt:
        return None
    field = gt.split(":")[0].replace("|", "/")   # GT is the first FORMAT sub-field
    if field in (".", "./.", "./", "/."):
        return None
    parts = [a for a in field.split("/") if a != ""]
    if not parts:
        return None
    idxs: list[int] = []
    for a in parts:
        if a == ".":
            return None                    # partial no-call → skip
        try:
            idxs.append(int(a))
        except ValueError:
            return None
    if len(idxs) == 1:                     # haploid (e.g. male X/Y) → homozygous
        idxs = [idxs[0], idxs[0]]
    alts = alt.split(",") if alt not in (".", "") else []
    alleles: list[str] = []
    for i in idxs[:2]:
        if i == 0:
            seq = ref
        elif 1 <= i <= len(alts):
            seq = alts[i - 1]
        else:
            return None                    # index out of range → malformed
        seq = seq.upper()
        if len(seq) != 1 or seq not in _ACGT:
            return None                    # indel / MNV / symbolic → skip
        alleles.append(seq)
    return "".join(sorted(alleles))


# ── ACMG actionable-gene windows (telemetry only; approximate, curated subset) ──
# A representative subset of ACMG SF v3 secondary-findings genes with published
# gene coordinates. Used ONLY to COUNT how many VCF variants fall in actionable
# genes (a value-quantifier for the not-yet-built Phase-2 screen) — never for a
# clinical call. Coordinates are gene-body windows, deliberately approximate.
_ACMG_GENES: list[tuple[str, str, int, int, int, int]] = [
    # gene, chrom, g37_start, g37_end, g38_start, g38_end
    ("BRCA1",  "17", 41196312, 41277500, 43044295, 43125483),
    ("BRCA2",  "13", 32889611, 32973805, 32315474, 32400266),
    ("MLH1",   "3",  37034841, 37092337, 36993350, 37050846),
    ("MSH2",   "2",  47630206, 47710367, 47403067, 47483228),
    ("APC",    "5",  112043195, 112181936, 112707498, 112846239),
    ("MUTYH",  "1",  45794835, 45806142, 45329170, 45340477),
    ("TP53",   "17", 7565097, 7590856, 7661779, 7687538),
    ("PTEN",   "10", 89622870, 89731687, 87863113, 87971930),
    ("LDLR",   "19", 11200038, 11244506, 11089362, 11133830),
    ("MYH7",   "14", 23881948, 23904870, 23412739, 23435661),
    ("MYBPC3", "11", 47352957, 47374253, 47331406, 47352702),
    ("KCNQ1",  "11", 2466221, 2870339, 2444986, 2849110),
    ("SCN5A",  "3",  38589553, 38691164, 38548062, 38649687),
    ("RYR2",   "1",  237205701, 238004523, 237041851, 237833988),
    ("PCSK9",  "1",  55505221, 55530525, 55039548, 55064852),
]


def _build_acmg_index(build: str):
    """chrom → list of (start, end, gene) for the detected build."""
    si, ei = (2, 3) if build == "grch37" else (4, 5)
    idx: dict[str, list[tuple[int, int, str]]] = {}
    for row in _ACMG_GENES:
        gene, chrom = row[0], row[1]
        idx.setdefault(chrom, []).append((row[si], row[ei], gene))
    return idx


# ── Registry position map (for build-gated back-fill) ─────────────────────────

def _registry_position_map(build: str) -> dict[tuple[str, int], object]:
    """{(chrom, pos): SNPRecord} for the detected build's coordinates."""
    out: dict[tuple[str, int], object] = {}
    for r in _reg._RECORDS:
        pos = r.pos_grch37 if build == "grch37" else r.pos_grch38
        if pos:
            out[(str(r.chrom), int(pos))] = r
    return out


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def enrich_and_profile_vcf(snps_df: pd.DataFrame, path: str, build: str,
                           log=None) -> tuple[pd.DataFrame, dict]:
    """Single streaming pass over a VCF: build-gated back-fill of curated
    registry positions the rsID reader dropped, plus a profile of what the file
    contains that the pipeline can't yet interpret.

    Returns (possibly-enriched snps_df, profile dict)."""
    _log = log or (lambda *a, **k: None)
    profile = {"total_variants": 0, "without_rsid": 0, "acmg_gene_variants": 0,
               "backfilled": 0, "build": build, "backfill_gated": False}

    do_backfill = build in ("grch37", "grch38")
    if not do_backfill:
        _log(f"  VCF build is '{build}' — refusing position back-fill "
             "(cannot safely pick GRCh37 vs 38 coordinates).")
        profile["backfill_gated"] = True

    reg_pos = _registry_position_map(build) if do_backfill else {}
    have_rsids = set(snps_df.index) if snps_df is not None else set()
    acmg_idx = _build_acmg_index(build) if do_backfill else {}
    new_rows: dict[str, dict] = {}

    opener = gzip.open if str(path).lower().endswith((".gz", ".bgz")) else open
    with opener(path, "rt", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            profile["total_variants"] += 1
            chrom = _norm_chrom(cols[0])
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            vid, ref, alt = cols[2], cols[3], cols[4]
            if vid in (".", ""):
                profile["without_rsid"] += 1

            # ACMG-gene membership (telemetry only)
            for (s, e, _gene) in acmg_idx.get(chrom, ()):
                if s <= pos <= e:
                    profile["acmg_gene_variants"] += 1
                    break

            # Build-gated registry back-fill for positions the rsID reader dropped
            rec = reg_pos.get((chrom, pos)) if do_backfill else None
            if rec is not None and rec.rsid not in have_rsids and rec.rsid not in new_rows:
                sample = cols[9] if len(cols) > 9 else ""
                geno = vcf_gt_to_genotype(sample, ref, alt)
                if geno:
                    new_rows[rec.rsid] = {"chrom": chrom, "pos": pos,
                                          "genotype": geno, "source": "wgs"}

    if new_rows:
        add = pd.DataFrame.from_dict(new_rows, orient="index")
        add.index.name = snps_df.index.name or "rsid"
        # align columns to the existing frame
        for c in snps_df.columns:
            if c not in add.columns:
                add[c] = None
        snps_df = pd.concat([snps_df, add[snps_df.columns]])
        profile["backfilled"] = len(new_rows)
        _log(f"  VCF back-fill: added {len(new_rows)} curated registry "
             f"variant(s) that had no rsID in the file (build {build}).")

    _log(f"  VCF profile: {profile['total_variants']:,} variants · "
         f"{profile['without_rsid']:,} without rsID · "
         f"~{profile['acmg_gene_variants']:,} in ACMG actionable genes "
         "(not yet interpreted — Phase-2 rare-variant engine).")
    return snps_df, profile
