"""
Novel / Rare-Variant Pathogenicity Interpretation — Phase 3 (WGS/VCF)
=====================================================================

Phase 2 (``clinical_variants.py``) screens a VCF against *known* ClinVar variants.
It says nothing about the millions of variants ClinVar has not classified. Phase 3
fills that gap: for every protein-altering variant the sample carries that is **not**
already a ClinVar P/LP hit, look it up in offline, tabix-indexed **computational
pathogenicity predictors** and surface the *predicted-damaging, rare* ones — always
labelled as **computational predictions, never clinical calls**.

Predictor registry (pluggable, each independently graceful):
  * **AlphaMissense** (CC BY 4.0 — commercial-OK) — primary missense predictor.
  * **gnomAD** (open) — population allele frequency → rarity gate.
  * **REVEL / CADD / SpliceAI** — non-commercial; used when present, dropped by
    ``--commercial-safe``.

Everything is queried by ``chrom:pos:ref:alt`` via ``pysam``/tabix — no Ensembl VEP.
If ``pysam`` or a predictor table is missing, the module degrades gracefully (it
never blocks the chip path, which has no need of it).

Educational screening, NOT a clinical diagnostic test.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

# Reuse the strand/indel-normalisation + zygosity helpers from Phase 2 so the two
# modules agree on what "the same variant" and "carries this allele" mean.
from .clinical_variants import _norm_chrom, norm_key, zygosity_for_alt

# .parent.parent: moved into the risk package, but reference/ stays at the
# repository root. Left at .parent these resolve inside risk/ and every
# consumer degrades quietly on a missing directory — so nothing raises
# and the screen silently reports no findings.
SCRIPT_DIR = Path(__file__).parent.parent
REFERENCE_DIR = SCRIPT_DIR / "reference"

try:                                    # pysam ships arm64 wheels; optional dep.
    import pysam
    _HAVE_PYSAM = True
except Exception:                       # pragma: no cover - import guard
    pysam = None
    _HAVE_PYSAM = False

# ── thresholds (documented, adjustable) ───────────────────────────────────────
_REVEL_DAMAGING = 0.5           # REVEL ≥ 0.5 ≈ likely pathogenic
_CADD_DAMAGING = 20.0           # CADD PHRED ≥ 20 = top 1% most deleterious
_SPLICEAI_DELTA = 0.5           # SpliceAI Δ ≥ 0.5 = likely splice-altering
_RARE_AF = 0.01                 # ≤ 1% = rare
_VERY_RARE_AF = 0.001           # ≤ 0.1% = very rare
_MAX_CARRIED_SNVS = 3_000_000   # safety cap on per-variant tabix queries

_NEGATIVE_DISCLAIMER = (
    "No result here means 'no predictor flagged a carried variant as damaging' — "
    "NOT 'no risk'. This is a COMPUTATIONAL screen: predictions (AlphaMissense, "
    "REVEL, CADD, SpliceAI) estimate deleteriousness, they do not diagnose. It "
    "covers only single-nucleotide missense/splice variants present in your file, "
    "excludes anything already in the ClinVar screen, and does not confirm any "
    "finding. It is a hypothesis-generating screen, not a clinical test.")

_DISCLAIMER = (
    "Computational predictions only. A predicted-pathogenic call is a statistical "
    "estimate from a machine-learning model, NOT a clinical determination — many "
    "predicted-damaging variants are benign in reality. Nothing here should change "
    "medical decisions until confirmed in an accredited diagnostic laboratory and "
    "discussed with a board-certified genetic counselor.")



# ── gene assignment for predicted variants ────────────────────────────────────
#
# A predicted variant used to reach the economics with no gene symbol, so every
# one of them routed to the same generic bucket regardless of where it landed.
# The gene is recoverable two ways, and they are NOT equally good:
#
#   1. UniProt accession, carried in the AlphaMissense row itself. This is a
#      real assignment: the predictor scored a specific protein and says which.
#   2. Coordinate falling inside a gene's span. This is a GUESS. Gene spans are
#      not exon models, so an intronic position, or one inside an overlapping
#      gene, is assigned confidently and wrongly.
#
# Both are recorded, and which one answered is recorded with them, because a
# reader deciding whether to trust a gene label needs to know whether it came
# from the predictor or from arithmetic on an interval.
#
# The accession map is DERIVED, not transcribed. It comes from joining this
# repository's own ClinVar and AlphaMissense tables on coordinates and taking
# the most frequent accession per gene over ~60 supporting positions — first-hit
# picked a non-canonical MLH1 isoform, majority vote gives P40692. Every entry
# below reproduces from local data; none is from memory.
UNIPROT_TO_GENE: dict[str, str] = {
    'P04114': 'APOB',
    'P38398': 'BRCA1',
    'P51587': 'BRCA2',
    'Q12809': 'KCNH2',
    'P51787': 'KCNQ1',
    'P01130': 'LDLR',
    'P40692': 'MLH1',
    'P43246': 'MSH2',
    'P52701': 'MSH6',
    'Q14896': 'MYBPC3',
    'P12883': 'MYH7',
    'Q86YC2': 'PALB2',
    'Q8NBP7': 'PCSK9',
    'P54278': 'PMS2',
    'P06400': 'RB1',
    'P07949': 'RET',
    'Q14524': 'SCN5A',
    'P04637': 'TP53',
}


def gene_for_uniprot(uniprot: str) -> str:
    """Gene symbol for a UniProt accession, or "" if it is not one we anchor."""
    return UNIPROT_TO_GENE.get((uniprot or "").strip().upper(), "")


def _window_gene(build: str, chrom: str, pos: int) -> str:
    """Gene whose span contains this position, or "". A guess, not an assignment."""
    try:
        from core.genome_input import _build_acmg_index
    except Exception:
        return ""
    idx = _build_acmg_index(build)
    c = chrom[3:] if str(chrom).lower().startswith("chr") else str(chrom)
    for (start, end, gene) in idx.get(c, ()):
        if start <= int(pos) <= end:
            return gene
    return ""


# ── predictor registry ────────────────────────────────────────────────────────

# Filename tokens that identify which reference build a predictor table is keyed
# on. Every table this module reads is coordinate-keyed, so a table for the wrong
# build does not degrade gracefully — it answers a different question. A GRCh37
# position looked up in an hg38 table returns None wherever the coordinate is
# unused (silent under-detection) and, where the coordinate happens to be valid
# in hg38 for some other variant, returns THAT variant's score. Nothing about the
# result says it is misattributed.
#
# Worked example from this repository's own tables — APOE rs429358, whose GRCh37
# and GRCh38 positions are both real coordinates:
#     lookup 19:45411941 (its GRCh37 position) in AlphaMissense_hg38 -> None
#     lookup 19:44908684 (its GRCh38 position) in AlphaMissense_hg38 -> 0.0365
# Before this map existed, ``resolve()`` took the first glob hit and the ``build``
# argument was checked for membership in ("grch37","grch38") and then never used,
# so a GRCh37 whole genome was scored against whichever table happened to be on
# disk. The only table setup.py downloads by default is hg38.
_BUILD_TOKENS: dict[str, tuple[str, ...]] = {
    "grch37": ("hg19", "grch37", "b37"),
    "grch38": ("hg38", "grch38", "b38"),
}


def table_build(path: str) -> str | None:
    """Which build a predictor table filename declares, or None if it is silent.

    A table that names no build cannot be checked, which is itself worth
    distinguishing from one that names a conflicting build: the first is
    unverifiable, the second is wrong.
    """
    low = str(path).lower()
    for build, tokens in _BUILD_TOKENS.items():
        if any(t in low for t in tokens):
            return build
    return None


class Predictor:
    """One pluggable pathogenicity predictor backed by a tabix-indexed table."""

    def __init__(self, name: str, axis: str, license: str, commercial_ok: bool,
                 parse: Callable[[list[list[str]], str, str], float | None],
                 subdir: str = "", glob: str = "", remote_url: str = ""):
        self.name = name
        self.axis = axis                 # 'missense' | 'splice' | 'rarity'
        self.license = license
        self.commercial_ok = commercial_ok
        self._parse = parse
        self.subdir = subdir
        self.glob = glob
        self.remote_url = remote_url
        self._tabix = None               # (TabixFile, chrom_fmt) once opened
        self.source = ""
        self.build = ""                  # build of the resolved table, if declared
        self.reject_reason = ""          # why resolve() refused, for the caller

    # -- table resolution + open ------------------------------------------------
    def resolve(self, build: str = "") -> bool:
        """Locate + open the table for ``build``, refusing a mismatched one.

        ``build`` is the reference build of the *input genome*. A table keyed on
        a different build is not a degraded answer, it is an answer about other
        variants, so this refuses instead of scoring and records why in
        ``reject_reason``. Passing no build keeps the old first-hit behaviour,
        which only the licence manifest needs.
        """
        self.reject_reason = ""
        if not _HAVE_PYSAM:
            self.reject_reason = "pysam not installed"
            return False
        path = None
        d = REFERENCE_DIR / self.subdir
        if d.is_dir() and self.glob:
            hits = sorted(str(h) for h in d.glob(self.glob))
            if build and hits:
                # Prefer a table whose filename declares the input's build. A
                # table that declares nothing stays eligible — unverifiable is
                # not the same as wrong — but a declared conflict is fatal.
                matching = [h for h in hits if table_build(h) == build]
                silent = [h for h in hits if table_build(h) is None]
                if matching:
                    path = matching[0]
                elif silent:
                    path = silent[0]
                else:
                    wrong = sorted({table_build(h) or "?" for h in hits})
                    self.reject_reason = (
                        f"{self.name} table on disk is keyed on "
                        f"{'/'.join(wrong)} but the input genome is {build}; "
                        f"refusing to score coordinates against the wrong "
                        f"build (run `python setup.py --predictors "
                        f"--build {build}`)")
                    return False
            elif hits:
                path = hits[0]
        if path is None and self.remote_url:
            if build and table_build(self.remote_url) not in (None, build):
                self.reject_reason = (
                    f"{self.name} remote table is keyed on "
                    f"{table_build(self.remote_url)}, input genome is {build}")
                return False
            path = self.remote_url          # pysam/htslib can open http(s) + .tbi
        if path is None:
            self.reject_reason = f"no {self.name} table found"
            return False
        self._tabix = _open_tabix(path)
        self.source = path
        self.build = table_build(path) or ""
        if self._tabix is None:
            self.reject_reason = f"could not open {path}"
            return False
        return True

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> float | None:
        if self._tabix is None:
            return None
        rows = _fetch(self._tabix, chrom, pos)
        if not rows:
            return None
        try:
            return self._parse(rows, ref.upper(), alt.upper())
        except Exception:
            return None


def _open_tabix(path: str):
    """Open a tabix file (local path or http(s) URL) and return
    ``(TabixFile, chrom_fmt)`` where ``chrom_fmt`` maps a plain chrom to the
    table's own naming ('1' vs 'chr1'), or None on failure."""
    try:
        tbx = pysam.TabixFile(str(path))
        contigs = set(tbx.contigs)
    except Exception:
        return None

    def fmt(chrom: str) -> str:
        c = _norm_chrom(chrom)
        if c in contigs:
            return c
        if f"chr{c}" in contigs:
            return f"chr{c}"
        if c == "MT" and "chrM" in contigs:
            return "chrM"
        return c
    return tbx, fmt


def _fetch(tabix_pair, chrom: str, pos: int) -> list[list[str]]:
    tbx, fmt = tabix_pair
    try:
        return [row.split("\t") for row in tbx.fetch(fmt(chrom), pos - 1, pos)]
    except Exception:
        return []


def _info_field(info: str, key: str) -> str | None:
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
        if kv == key:
            return ""
    return None


# -- per-predictor row parsers -------------------------------------------------

def _parse_alphamissense(rows, ref, alt):
    """AlphaMissense TSV: CHROM POS REF ALT genome uniprot transcript prot_var
    am_pathogenicity am_class. Multiple transcript rows per locus → keep MAX
    pathogenicity (returned as a (score, class) tuple)."""
    best_s: float | None = None
    best_c = ""
    best_u = ""
    for c in rows:
        if len(c) < 10 or c[2].upper() != ref or c[3].upper() != alt:
            continue
        try:
            s = float(c[8])
        except ValueError:
            continue
        if best_s is None or s > best_s:
            # The UniProt accession travels with the score. It was being read
            # past and dropped, which is why a predicted variant reached the
            # economics anonymous: the predictor knew which protein it had
            # scored and nothing asked.
            best_s, best_c, best_u = s, c[9], c[5]
    if best_s is None:
        return None
    return (best_s, best_c, best_u)


def _parse_revel(rows, ref, alt):
    """REVEL table (bgzipped CSV→TSV): ... ref alt ... REVEL ... Match ref/alt,
    return MAX REVEL score. Column layout after setup normalisation:
    chr pos ref alt REVEL (5-col slim table produced by setup_revel)."""
    best = None
    for c in rows:
        if len(c) < 5 or c[2].upper() != ref or c[3].upper() != alt:
            continue
        try:
            v = float(c[4])
        except ValueError:
            continue
        best = v if best is None or v > best else best
    return best


def _parse_cadd(rows, ref, alt):
    """CADD whole_genome_SNVs.tsv: Chrom Pos Ref Alt RawScore PHRED. Return PHRED."""
    for c in rows:
        if len(c) < 6 or c[2].upper() != ref or c[3].upper() != alt:
            continue
        try:
            return float(c[5])
        except ValueError:
            return None
    return None


def _parse_spliceai(rows, ref, alt):
    """SpliceAI precomputed VCF: INFO SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|
    DS_DL|.... Return max delta (DS_*) for the matching ALT."""
    best = None
    for c in rows:
        if len(c) < 8 or c[3].upper() != ref:
            continue
        raw = _info_field(c[7], "SpliceAI")
        if not raw:
            continue
        for entry in raw.split(","):
            parts = entry.split("|")
            if len(parts) < 6 or parts[0].upper() != alt:
                continue
            try:
                delta = max(float(parts[2]), float(parts[3]),
                            float(parts[4]), float(parts[5]))
            except ValueError:
                continue
            best = delta if best is None or delta > best else best
    return best


def _parse_gnomad(rows, ref, alt):
    """gnomAD sites VCF: return AF for the matching ALT (handles multiallelic)."""
    for c in rows:
        if len(c) < 8 or c[3].upper() != ref:
            continue
        alts = [a.upper() for a in c[4].split(",")]
        if alt not in alts:
            continue
        af = _info_field(c[7], "AF")
        if af is None:
            return None
        parts = af.split(",")
        idx = alts.index(alt)
        try:
            return float(parts[idx]) if idx < len(parts) else float(parts[0])
        except ValueError:
            return None
    return None


def _registry() -> list[Predictor]:
    return [
        Predictor("AlphaMissense", "missense", "CC BY 4.0", True,
                  _parse_alphamissense, subdir="alphamissense",
                  glob="AlphaMissense_*.tsv.gz"),
        Predictor("gnomAD", "rarity", "open (no restrictions)", True,
                  _parse_gnomad, subdir="gnomad", glob="gnomad*.vcf.*gz"),
        Predictor("REVEL", "missense", "non-commercial", False,
                  _parse_revel, subdir="revel", glob="revel*.tsv.gz"),
        Predictor("SpliceAI", "splice", "CC BY-NC 4.0", False,
                  _parse_spliceai, subdir="spliceai", glob="spliceai*.vcf.*gz"),
        Predictor("CADD", "meta", "non-commercial", False,
                  _parse_cadd, subdir="cadd", glob="*whole_genome_SNVs*.tsv.gz",
                  remote_url=_CADD_REMOTE_GRCH38),
    ]


# Public CADD v1.7 whole-genome SNV tables (tabix-indexed) — used for a remote
# demo query on PUBLIC data only (never the user's private genome).
_CADD_REMOTE_GRCH38 = ("https://krishna.gs.washington.edu/download/CADD/v1.7/"
                       "GRCh38/whole_genome_SNVs.tsv.gz")


# ── main analysis ─────────────────────────────────────────────────────────────

def analyze_novel_variants(vcf_path: str, build: str,
                           clinvar_result: dict | None = None,
                           inferred_sex: str | None = None,
                           commercial_safe: bool = False,
                           allow_remote: bool = False,
                           log=None) -> dict:
    """Screen a user VCF for predicted-damaging novel/rare variants. Mirrors the
    ``clinical_variants`` result-dict contract for drop-in wiring."""
    _log = log or (lambda *a, **k: None)

    if not _HAVE_PYSAM:
        return _unavailable("pysam not installed — run `pip install pysam` "
                            "(optional 'wgs' extra) to enable Phase-3 predictors.")
    if build not in ("grch37", "grch38"):
        return _unavailable(f"Build '{build}' — cannot select predictor tables "
                            "safely (try --assume-build grch38).")

    # Resolve which predictor tables are actually available FOR THIS BUILD.
    # A table keyed on the other build is refused rather than used: see
    # _BUILD_TOKENS for why a wrong-build lookup is worse than no lookup.
    predictors: list[Predictor] = []
    dropped_nc: list[str] = []
    dropped_build: list[str] = []
    for p in _registry():
        if commercial_safe and not p.commercial_ok:
            dropped_nc.append(p.name)
            continue
        pdir = REFERENCE_DIR / p.subdir
        local_exists = pdir.is_dir() and any(pdir.glob(p.glob))
        # Only fall back to a remote table (CADD demo) when explicitly allowed.
        if not local_exists and not (p.remote_url and allow_remote):
            continue
        if p.resolve(build):
            predictors.append(p)
            _log(f"  Predictor ready: {p.name} ({p.license})"
                 + (f" [{p.build}]" if p.build else ""))
        elif p.reject_reason:
            dropped_build.append(p.reject_reason)
            _log(f"  Predictor REFUSED: {p.reject_reason}")

    if not predictors:
        if dropped_build:
            return _unavailable(
                "No predictor table matches this genome's reference build. "
                + " · ".join(dropped_build))
        return _unavailable(
            "No Phase-3 predictor tables found — run `python setup.py --predictors` "
            "(AlphaMissense + REVEL + gnomAD). Falls back cleanly; no screen performed.")

    have = {p.name for p in predictors}
    gnomad_present = "gnomAD" in have

    # Variants already flagged by the ClinVar screen — skip to add only new signal.
    known: set = set()
    for f in (clinvar_result or {}).get("findings") or []:
        known.add(norm_key(f["chrom"], f["pos"], f["ref"], f["alt"]))

    findings: list[dict] = []
    n_scanned = 0
    n_queried = 0
    import gzip as _gz
    opener = _gz.open if str(vcf_path).lower().endswith((".gz", ".bgz")) else open
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
            if len(ref) != 1 or ref.upper() not in "ACGT":
                continue                          # SNV predictors only
            sample = cols[9] if len(cols) > 9 else ""
            for ai, alt in enumerate(alt_field.split(","), start=1):
                alt = alt.upper()
                if len(alt) != 1 or alt not in "ACGT":
                    continue
                zyg = zygosity_for_alt(sample, ai) if sample else "unknown"
                if zyg is None:                   # sample doesn't carry this allele
                    continue
                if norm_key(chrom, pos, ref, alt) in known:
                    continue                      # already in ClinVar screen
                if n_queried >= _MAX_CARRIED_SNVS:
                    continue
                n_queried += 1
                rec = _score_variant(predictors, chrom, pos, ref, alt)
                if rec is None:
                    continue
                rec.update({"chrom": _norm_chrom(chrom), "pos": pos,
                            "ref": ref.upper(), "alt": alt, "zygosity": zyg})
                findings.append(rec)

    result = _classify(findings, n_scanned, n_queried, gnomad_present, build)
    result["predictors_used"] = [{"name": p.name, "license": p.license,
                                  "commercial_ok": p.commercial_ok,
                                  "table_build": p.build,
                                  "table": p.source} for p in predictors]
    result["commercial_safe"] = commercial_safe
    result["dropped_noncommercial"] = dropped_nc
    # Which predictors were refused for being keyed on the wrong build. Reported
    # rather than dropped silently: a screen that skipped half its predictors is
    # a different screen, and the reader has to be able to see that it did.
    result["input_build"] = build
    result["dropped_build_mismatch"] = dropped_build
    # Name only the pathogenicity predictors actually used in the disclaimer, so a
    # --commercial-safe report never even mentions the non-commercial ones.
    used = ", ".join(p.name for p in predictors if p.axis != "rarity") or "AlphaMissense"
    result["negative_disclaimer"] = _NEGATIVE_DISCLAIMER.replace(
        "(AlphaMissense, REVEL, CADD, SpliceAI)", f"({used})")
    return result


def _score_variant(predictors: list[Predictor], chrom, pos, ref, alt) -> dict | None:
    """Query every predictor; return a finding dict only if the ensemble flags the
    variant as damaging or splice-altering. Otherwise None."""
    am_score = am_class = None
    uniprot = ""
    revel = cadd = splice = af = None
    for p in predictors:
        v = p.lookup(chrom, pos, ref, alt)
        if v is None:
            continue
        if p.name == "AlphaMissense":
            am_score, am_class, uniprot = v
        elif p.name == "REVEL":
            revel = v
        elif p.name == "CADD":
            cadd = v
        elif p.name == "SpliceAI":
            splice = v
        elif p.name == "gnomAD":
            af = v

    missense_hits = []
    if am_class == "likely_pathogenic":
        missense_hits.append("AlphaMissense")
    if revel is not None and revel >= _REVEL_DAMAGING:
        missense_hits.append("REVEL")
    if cadd is not None and cadd >= _CADD_DAMAGING:
        missense_hits.append("CADD")
    splice_hit = splice is not None and splice >= _SPLICEAI_DELTA
    ambiguous = am_class == "ambiguous"

    if not missense_hits and not splice_hit and not ambiguous:
        return None

    return {
        "am_score": am_score, "am_class": am_class, "uniprot": uniprot,
        "revel": revel, "cadd_phred": cadd, "spliceai_delta": splice,
        "gnomad_af": af, "consensus": missense_hits,
        "n_missense_hits": len(missense_hits), "splice_hit": splice_hit,
    }


def _rarity(af: float | None) -> str:
    if af is None:
        return "unknown"
    if af <= _VERY_RARE_AF:
        return "very_rare"
    if af <= _RARE_AF:
        return "rare"
    if af <= 0.05:
        return "uncommon"
    return "common"


def _confidence(rec: dict) -> str:
    n = rec["n_missense_hits"]
    am = rec.get("am_score")
    if rec["splice_hit"] or n >= 2 or (am is not None and am >= 0.9):
        return "higher"
    if n == 1 or (am is not None and am >= 0.7):
        return "moderate"
    return "low"


def _classify(findings: list[dict], n_scanned: int, n_queried: int,
              gnomad_present: bool, build: str = "") -> dict:
    buckets: dict[str, list[dict]] = {
        "predicted_splice_disrupting": [], "predicted_pathogenic_rare": [],
        "predicted_pathogenic_uncommon": [], "predicted_pathogenic_common": [],
        "ambiguous": [],
    }
    for f in findings:
        rar = _rarity(f.get("gnomad_af"))
        f["rarity"] = rar
        f["confidence"] = _confidence(f)
        # GENE ASSIGNMENT, WITH ITS PROVENANCE. The UniProt accession comes
        # from the predictor and names the protein it actually scored; a
        # coordinate falling inside a gene's span is a guess, because spans are
        # not exon models. Both are usable, they are not equally trustworthy,
        # and which one answered is recorded so a reader can tell.
        if not f.get("gene"):
            g = gene_for_uniprot(f.get("uniprot", ""))
            if g:
                f["gene"], f["gene_basis"] = g, "uniprot_accession"
            else:
                g = _window_gene(build, f.get("chrom", ""), f.get("pos", 0))
                if g:
                    f["gene"], f["gene_basis"] = g, "coordinate_window"
        f["gene"] = f.get("gene", "")
        f.setdefault("gene_basis", "unassigned")
        parts = []
        if f.get("am_class"):
            parts.append(f"AlphaMissense {f['am_class'].replace('_', ' ')} "
                         f"({f['am_score']:.2f})")
        if f.get("revel") is not None:
            parts.append(f"REVEL {f['revel']:.2f}")
        if f.get("cadd_phred") is not None:
            parts.append(f"CADD {f['cadd_phred']:.0f}")
        if f.get("spliceai_delta") is not None:
            parts.append(f"SpliceAI Δ{f['spliceai_delta']:.2f}")
        af = f.get("gnomad_af")
        parts.append(f"gnomAD AF {af:.2e}" if af is not None else "gnomAD AF n/a")
        f["evidence"] = "; ".join(parts)
        f["interpretation"] = (
            "Computational prediction — a machine-learning model estimates this "
            "variant is damaging. NOT a clinical diagnosis; confirm in an "
            "accredited lab before it means anything.")

        if f["splice_hit"]:
            buckets["predicted_splice_disrupting"].append(f)
        elif f["am_class"] == "ambiguous" and f["n_missense_hits"] == 0:
            buckets["ambiguous"].append(f)
        elif rar in ("rare", "very_rare", "unknown"):
            buckets["predicted_pathogenic_rare"].append(f)
        elif rar == "uncommon":
            buckets["predicted_pathogenic_uncommon"].append(f)
        else:
            buckets["predicted_pathogenic_common"].append(f)

    conf_rank = {"higher": 0, "moderate": 1, "low": 2}
    for b in buckets.values():
        b.sort(key=lambda x: (conf_rank.get(x["confidence"], 3),
                              -(x.get("am_score") or 0)))

    n_rare_damaging = (len(buckets["predicted_pathogenic_rare"])
                       + len(buckets["predicted_splice_disrupting"]))
    return {
        "available": True,
        "n_scanned": n_scanned,
        "n_queried": n_queried,
        "n_predicted_pathogenic": sum(len(v) for k, v in buckets.items()
                                      if k != "ambiguous"),
        "n_rare_damaging": n_rare_damaging,
        "n_ambiguous": len(buckets["ambiguous"]),
        "gnomad_present": gnomad_present,
        "buckets": buckets,
        "findings": [f for b in buckets.values() for f in b],
        "negative_disclaimer": _NEGATIVE_DISCLAIMER,
        "disclaimer": _DISCLAIMER,
    }


def _unavailable(reason: str) -> dict:
    return {"available": False, "reason": reason,
            "negative_disclaimer": _NEGATIVE_DISCLAIMER, "disclaimer": _DISCLAIMER}
