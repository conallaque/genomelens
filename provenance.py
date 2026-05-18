"""
Variant provenance + reference-build detection.

V8 introduces two parse-time annotations on ``snps_df`` that downstream
consumers can opt into:

  * ``snps_df["source"]`` — one of ``"chip"``, ``"imp_high_r2"``,
    ``"imp_low_r2"``. Every variant read from the raw chip file is
    tagged ``"chip"``; the imputation module is expected to overwrite
    imputed rows with one of the ``"imp_*"`` tags based on Beagle's
    DR2 quality score. Consumers (PRS, PheWAS, etc.) can read this
    column to weight or filter variants — the schema is added now;
    the actual reweighting in scoring modules is V8.1 work.

  * ``snps_df.attrs["build"]`` — one of ``"grch37"``, ``"grch38"``,
    or ``"unknown"``. Detected by sampling a small set of registered
    SNPs whose positions are known in both builds. Surfaced as a QC
    field so the report can warn users when a chip is on a less
    common assembly. Position-based fallback lookups in any module
    can use this to pick the right coordinate column.

Both annotations are *additive* — they extend the DataFrame without
breaking existing consumers that ignore them.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

import snp_registry

logger = logging.getLogger(__name__)


# ── Provenance tagging ──────────────────────────────────────────────────────

CHIP = "chip"
IMP_HIGH_R2 = "imp_high_r2"
IMP_LOW_R2 = "imp_low_r2"

# DR2 threshold separating "high-confidence" imputed variants from "low".
# Below this, scoring modules should down-weight or skip the variant.
# Default matches Beagle 5.4's commonly used --r2 filter; adjustable.
DEFAULT_HIGH_R2_THRESHOLD = 0.8


def tag_chip_source(snps_df: pd.DataFrame) -> pd.DataFrame:
    """Add a ``source`` column to a freshly-parsed chip DataFrame, marking
    every row as ``"chip"``. Idempotent — re-running on a tagged DataFrame
    is a no-op.

    This is called once at the end of ``parse_dna_file``. Imputation
    runs *afterwards* and overwrites the appropriate rows with
    ``"imp_high_r2"`` / ``"imp_low_r2"``.
    """
    if "source" in snps_df.columns:
        return snps_df
    snps_df = snps_df.copy()
    snps_df["source"] = CHIP
    return snps_df


def tag_imputed_rows(
    snps_df: pd.DataFrame,
    chip_rsids: set[str],
    dr2_by_rsid: Optional[dict[str, float]] = None,
    high_r2_threshold: float = DEFAULT_HIGH_R2_THRESHOLD,
) -> pd.DataFrame:
    """Mark every variant that is *not* in ``chip_rsids`` as imputed.

    Parameters
    ----------
    snps_df : DataFrame indexed by rsID with a ``source`` column.
    chip_rsids : the set of rsIDs that came from the raw chip file
        (i.e. were present *before* the imputation merge).
    dr2_by_rsid : optional rsID → Beagle DR2 mapping. When supplied,
        imputed variants with DR2 ≥ ``high_r2_threshold`` are tagged
        ``"imp_high_r2"``; the rest ``"imp_low_r2"``. Without DR2,
        every imputed row gets the conservative ``"imp_low_r2"`` tag.

    Returns the updated DataFrame.
    """
    if "source" not in snps_df.columns:
        snps_df = tag_chip_source(snps_df)
    else:
        snps_df = snps_df.copy()

    dr2_by_rsid = dr2_by_rsid or {}
    imputed_mask = ~snps_df.index.isin(chip_rsids)
    for rsid in snps_df.index[imputed_mask]:
        dr2 = dr2_by_rsid.get(rsid)
        if dr2 is not None and dr2 >= high_r2_threshold:
            snps_df.at[rsid, "source"] = IMP_HIGH_R2
        else:
            snps_df.at[rsid, "source"] = IMP_LOW_R2
    return snps_df


def provenance_summary(snps_df: pd.DataFrame) -> dict:
    """Return per-source counts for the QC card / report header.

    Always safe to call — returns zeros if the ``source`` column is absent.
    """
    if "source" not in snps_df.columns:
        return {CHIP: len(snps_df), IMP_HIGH_R2: 0, IMP_LOW_R2: 0,
                "total": len(snps_df), "tagged": False}
    counts = snps_df["source"].value_counts().to_dict()
    return {
        CHIP: int(counts.get(CHIP, 0)),
        IMP_HIGH_R2: int(counts.get(IMP_HIGH_R2, 0)),
        IMP_LOW_R2: int(counts.get(IMP_LOW_R2, 0)),
        "total": len(snps_df),
        "tagged": True,
    }


# ── Reference-build auto-detection ──────────────────────────────────────────

# Build the probe set lazily so we don't pay the cost at module load.
_PROBE_RSIDS: list[str] | None = None


def _build_probe_set() -> list[str]:
    """Pick a set of registry SNPs that have *both* GRCh37 and GRCh38
    positions and are reasonably likely to be on consumer chips. Sampling
    is biased toward common SNPs (MTHFR, COMT, APOE, FADS1, ACTN3, …)
    because those are present on essentially every consumer panel.
    """
    global _PROBE_RSIDS
    if _PROBE_RSIDS is None:
        _PROBE_RSIDS = [
            r.rsid for r in snp_registry.SNPS.values()
            if r.pos_grch37 is not None and r.pos_grch38 is not None
            and r.pos_grch37 != r.pos_grch38   # must distinguish builds
        ]
    return _PROBE_RSIDS


def detect_build(snps_df: pd.DataFrame, min_probes: int = 3) -> str:
    """Detect whether ``snps_df`` is on GRCh37 or GRCh38 by checking probe
    positions against the registry's dual coordinates.

    Returns one of ``"grch37"``, ``"grch38"``, ``"mixed"``, or
    ``"unknown"``:

      * ``grch37`` / ``grch38`` — clear majority match for one build
      * ``mixed`` — probes match different builds (the TellmeGen-style
        case where per-probe coordinates drift); position-based lookups
        for this chip should fall back to rsID matching only
      * ``unknown`` — fewer than ``min_probes`` probes resolved on either
        build (rare; a very sparse chip or a non-standard format)

    The chip's ``pos`` column is expected to be the second column of the
    DataFrame; we look up by rsID index and compare ``df.loc[rsid, "pos"]``
    against the two registry positions.
    """
    if "pos" not in snps_df.columns:
        return "unknown"

    probes = _build_probe_set()
    if not probes:
        return "unknown"

    n_grch37 = 0
    n_grch38 = 0
    for rsid in probes:
        if rsid not in snps_df.index:
            continue
        rec = snp_registry.SNPS[rsid]
        try:
            chip_pos = int(snps_df.loc[rsid, "pos"])
        except (ValueError, TypeError):
            continue
        if chip_pos == rec.pos_grch37:
            n_grch37 += 1
        elif chip_pos == rec.pos_grch38:
            n_grch38 += 1

    total_resolved = n_grch37 + n_grch38
    if total_resolved < min_probes:
        return "unknown"

    # Require a clear majority (>= 75% of resolved probes on one build).
    if n_grch37 >= 0.75 * total_resolved:
        return "grch37"
    if n_grch38 >= 0.75 * total_resolved:
        return "grch38"
    return "mixed"


def annotate_build(snps_df: pd.DataFrame) -> pd.DataFrame:
    """Detect the chip's reference build and stash it on ``df.attrs``.

    Stored under ``df.attrs["build"]`` so it survives copy operations
    (per pandas docs ``attrs`` is preserved across many DataFrame ops,
    though not all — consumers needing it should re-detect after heavy
    transformations).
    """
    build = detect_build(snps_df)
    snps_df.attrs["build"] = build
    if build == "mixed":
        logger.warning(
            "Detected mixed GRCh37/GRCh38 chip coordinates — position-based "
            "fallback lookups may be unreliable. rsID matching is unaffected."
        )
    elif build == "unknown":
        logger.warning(
            "Could not detect reference build from probe positions. "
            "Position-based fallback lookups will assume GRCh37."
        )
    return snps_df


# ── Combined parse-time annotation ──────────────────────────────────────────

def annotate_parsed(snps_df: pd.DataFrame) -> pd.DataFrame:
    """One-shot helper: ``parse_dna_file`` calls this at the end so every
    downstream consumer sees a DataFrame tagged with source + build.

    Equivalent to:
        snps_df = tag_chip_source(snps_df)
        snps_df = annotate_build(snps_df)
    """
    snps_df = tag_chip_source(snps_df)
    snps_df = annotate_build(snps_df)
    return snps_df
