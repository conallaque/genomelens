"""
TNRC18 rs117910193 marker lookup  (novelty / "for fun" feature)
===============================================================

A single-variant genotype lookup at:

    Gene            TNRC18 (Trinucleotide Repeat-Containing Gene 18 Protein)
    Chromosome      7
    Position        5401412
    dbSNP           rs117910193
    Variant         G>A  (reference allele G, minor allele A)

Classification requested:
    homozygous  G/G   → "Target Trait Marker (Wild Type)"
    heterozygous G/A  → "Standard Marker"

This is a plain genotype read-out. It does NOT assert any health, trait, disease, or
ancestry meaning — no such association is claimed or validated here. It simply reports
which allele pair you carry at this position and applies the labels above. Treat it as a
novelty lookup, not a result to act on.

Design notes:
  * Looks the variant up by rsID first, then falls back to a chromosome:position match,
    so it works on both chip files (rsID-keyed) and whole-genome VCFs (position-keyed).
  * Strand-aware: a genotype reported on the opposite strand (C/C or C/T) is recognized
    and normalized, so a strand flip cannot silently mis-call the marker.
"""

from __future__ import annotations

RSID = "rs117910193"
GENE = "TNRC18"
CHROM = "7"
POSITION = 5401412
REF_ALLELE = "G"       # wild-type / reference
ALT_ALLELE = "A"       # minor / variant

# Complement, for strand normalization (the variant may be reported on either strand).
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _norm_chrom(value: object) -> str:
    """Normalize a chromosome label to a bare number/letter ('chr7' -> '7')."""
    s = str(value).strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    return s.upper() if s in ("x", "y", "mt", "m") else s


def _clean_genotype(gt: object) -> str | None:
    """Extract the two called bases from a genotype string, or None if not a clean
    biallelic SNP call. Accepts 'GA', 'G/A', 'G|A', etc.; rejects indels/no-calls."""
    if gt is None:
        return None
    s = str(gt).upper().replace("/", "").replace("|", "").replace(" ", "").strip()
    if len(s) != 2 or any(b not in "ACGT" for b in s):
        return None
    return s


def _orient_to_ref_strand(alleles: str) -> str:
    """Return the genotype expressed on the reference (G/A) strand.

    The call may come on the opposite strand as C/C (=G/G) or C/T (=G/A). If the
    bases already belong to the {G, A} alphabet, they are returned as-is; if they
    belong to the complementary {C, T} alphabet, they are complemented. Mixed or
    ambiguous alphabets are returned unchanged (handled as 'other' downstream)."""
    ref_alphabet = {REF_ALLELE, ALT_ALLELE}                 # {'G','A'}
    comp_alphabet = {_COMPLEMENT[REF_ALLELE], _COMPLEMENT[ALT_ALLELE]}  # {'C','T'}
    bases = set(alleles)
    if bases <= ref_alphabet:
        return alleles
    if bases <= comp_alphabet:
        return "".join(_COMPLEMENT[b] for b in alleles)
    return alleles                                          # not on this locus's axis


def _locate(snps_df) -> dict | None:
    """Find the marker row by rsID, then by chrom:pos. Returns a small dict or None."""
    if snps_df is None or getattr(snps_df, "empty", True):
        return None
    # 1) rsID lookup (chip files and most exports)
    try:
        if RSID in snps_df.index:
            row = snps_df.loc[RSID]
            if hasattr(row, "iloc") and getattr(row, "ndim", 1) > 1:
                row = row.iloc[0]                            # duplicate index guard
            return {"how": "rsID", "genotype": row.get("genotype")}
    except Exception:
        pass
    # 2) positional fallback (whole-genome VCF without rsIDs)
    try:
        if {"chrom", "pos"} <= set(snps_df.columns):
            hit = snps_df[(snps_df["chrom"].map(_norm_chrom) == CHROM)
                          & (snps_df["pos"].astype("int64", errors="ignore") == POSITION)]
            if len(hit):
                return {"how": "chrom:pos", "genotype": hit.iloc[0].get("genotype")}
    except Exception:
        pass
    return None


def analyze_tnrc18_marker(snps_df) -> dict:
    """Look up rs117910193 and classify the genotype.

    Returns a dict with availability, the raw and strand-oriented genotype, a
    zygosity classification, and the requested marker label. Always returns cleanly
    (never raises) so it is safe to drop into a pipeline.
    """
    base = {
        "available": False,
        "rsid": RSID, "gene": GENE, "chrom": CHROM, "position": POSITION,
        "variant": f"{REF_ALLELE}>{ALT_ALLELE}",
        "disclaimer": ("Novelty genotype lookup only — no validated health, trait, or "
                       "ancestry association is claimed. Not medical advice."),
    }
    found = _locate(snps_df)
    if not found:
        base["reason"] = (f"{RSID} (chr{CHROM}:{POSITION}) was not typed in this file "
                          "— common on genotyping chips, which cover a sparse subset "
                          "of positions.")
        return base

    raw = _clean_genotype(found["genotype"])
    if raw is None:
        base["reason"] = (f"{RSID} present but not a clean biallelic SNP call "
                          f"(got {found['genotype']!r}).")
        base["matched_by"] = found["how"]
        return base

    oriented = _orient_to_ref_strand(raw)
    n_alt = sum(1 for b in oriented if b == ALT_ALLELE)
    n_ref = sum(1 for b in oriented if b == REF_ALLELE)

    # Classify per the requested logic.
    if n_ref == 2:
        zygosity = "homozygous reference"
        marker = "Target Trait Marker (Wild Type)"
        detail = (f"Homozygous {REF_ALLELE}{REF_ALLELE} — carries the target "
                  "wild-type marker at this position.")
    elif n_ref == 1 and n_alt == 1:
        zygosity = "heterozygous"
        marker = "Standard Marker"
        detail = (f"Heterozygous {REF_ALLELE}{ALT_ALLELE} — one reference and one "
                  "minor allele (the standard marker).")
    elif n_alt == 2:
        zygosity = "homozygous minor"
        marker = "Neither (homozygous minor allele)"
        detail = (f"Homozygous {ALT_ALLELE}{ALT_ALLELE} — two copies of the minor "
                  "allele; matches neither requested category.")
    else:
        zygosity = "off-axis genotype"
        marker = "Undetermined"
        detail = (f"Genotype {raw} does not sit on this locus's {REF_ALLELE}/"
                  f"{ALT_ALLELE} axis; cannot be classified.")

    base.update({
        "available": True,
        "matched_by": found["how"],
        "genotype_raw": raw,
        "genotype_oriented": oriented,
        "strand_flipped": (raw != oriented),
        "zygosity": zygosity,
        "alt_allele_count": n_alt,
        "marker": marker,
        "is_target_wild_type": (marker == "Target Trait Marker (Wild Type)"),
        "detail": detail,
    })
    return base


def build_tnrc18_html(result: dict | None) -> str:
    """Minimal, self-contained HTML card. Returns '' if there's nothing to show."""
    if not result:
        return ""
    if not result.get("available"):
        return (
            '<section class="tnrc18-card" id="tnrc18-marker" '
            'style="border:1px solid #e3e7ec;border-radius:10px;padding:14px 16px;'
            'margin:12px 0"><div style="font-weight:700">TNRC18 marker '
            f'({result["rsid"]})</div><div style="font-size:.85em;color:#8a94a3;'
            f'margin-top:4px">{result.get("reason","Not available.")}</div></section>')

    is_target = result.get("is_target_wild_type")
    accent = "#22683f" if is_target else "#5b6673"
    return (
        f'<section class="tnrc18-card" id="tnrc18-marker" '
        f'style="border:1px solid #e3e7ec;border-left:4px solid {accent};'
        f'border-radius:10px;padding:14px 16px;margin:12px 0">'
        f'<div style="font-weight:700;color:{accent}">TNRC18 marker · {result["rsid"]}</div>'
        f'<div style="font-size:1.4em;font-weight:700;margin:6px 0">{result["marker"]}</div>'
        f'<div style="font-size:.9em;color:#48545f">Genotype '
        f'<strong>{result["genotype_oriented"]}</strong> ('
        f'{result["zygosity"]}{" · strand-corrected" if result.get("strand_flipped") else ""}) '
        f'at chr{result["chrom"]}:{result["position"]} · matched by {result["matched_by"]}.</div>'
        f'<div style="font-size:.82em;color:#6a7683;margin-top:6px">{result["detail"]}</div>'
        f'<div style="font-size:.75em;color:#b0868a;margin-top:8px">{result["disclaimer"]}</div>'
        f'</section>')
