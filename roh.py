"""
Runs of Homozygosity (ROH) Analysis
===================================

Detects contiguous stretches of homozygous genotypes across the autosomes.
Long ROH reveal population history (founder effects, drift, isolation)
and — at the longest scales — recent parental relatedness (consanguinity).

Algorithm (PLINK-like sliding window):
    1. For each autosome, walk the SNPs in physical order.
    2. Extend the run while a trailing window of `WINDOW_SNPS` SNPs holds no
       more than `max_hets` heterozygotes. The window is what slides — a het
       tolerated early does not consume the budget for the rest of the run.
    3. Break the run across any inter-SNP gap wider than `MAX_GAP_MB`, since an
       uncalled stretch (centromere, coverage hole) is absence of evidence, not
       evidence of homozygosity.
    4. Emit a run when it reaches `min_snps` SNPs AND `min_length_mb` Mb.

Output:
    * Per-run records (chr, start, end, n_snps, length_mb).
    * Classification: short (1-2 Mb), medium (2-8 Mb), long (>8 Mb).
    * F_ROH coefficient = total ROH length / autosomal genome length.
    * F_ROH_long — the same ratio over long (>8 Mb) runs only.
    * SVG ideogram visualisation of ROH locations.

Population context (illustrative, not used for hard cutoffs). Note these tiers
are keyed on F_ROH_LONG, not total F_ROH: only long ROH track recent pedigree
relatedness (Kirin et al. 2010, PLoS ONE 5:e13996 — r≈0.87 for sum-of-ROH
>10 Mb). Short and medium runs largely reflect ancient shared ancestry, so
including them mislabels founder-population background as consanguinity.
    * F_ROH_long ~ 0      — no evidence of recent parental relatedness
    * F_ROH_long to 0.02  — distant relatedness or strong founder background
    * F_ROH_long to 0.05  — consistent with second-cousin parents
    * F_ROH_long > 0.05   — consistent with first-cousin or closer

Framing: consanguinity is culturally normal in many populations. This is
informational, not judgmental. The clinically relevant consequence is
elevated rare-recessive-disease risk for offspring of related parents.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Approximate autosomal genome length in GRCh37/38 (Mb)
AUTOSOMAL_LENGTH_MB = 2881.0
# Sliding-window width (SNPs) over which `max_hets` is enforced. PLINK's
# --homozyg-window-snp default; the het budget applies per window, not per run.
WINDOW_SNPS = 50
# Maximum tolerated distance between consecutive genotyped SNPs inside one run.
# An uncalled stretch wider than this (centromere, assay coverage hole) carries
# no homozygosity evidence, so a run must not be extended across it. PLINK's
# --homozyg-gap default is 1 Mb; 1.0 Mb is kept here for the same reason.
MAX_GAP_MB = 1.0
# Length (Mb) at or above which a run counts as "long" — the only class that
# tracks recent pedigree relatedness (Kirin 2010).
LONG_ROH_MB = 8.0
# Approximate per-chromosome lengths (Mb) for ideogram rendering
CHROM_LENGTHS_MB = {
    "1": 249.3, "2": 243.2, "3": 198.0, "4": 191.2, "5": 180.9, "6": 171.1,
    "7": 159.1, "8": 146.4, "9": 141.2, "10": 135.5, "11": 135.0, "12": 133.9,
    "13": 115.2, "14": 107.3, "15": 102.5, "16": 90.4, "17": 81.2, "18": 78.1,
    "19": 59.1, "20": 63.0, "21": 48.1, "22": 51.3,
}

# Disease-relevant gene regions (chr, start_Mb, end_Mb, gene, condition)
# Light-touch annotation so we can flag ROH overlapping known disease genes
DISEASE_GENE_REGIONS: list[tuple] = [
    ("1", 11.0, 12.0, "MTHFR", "Folate metabolism"),
    ("1", 11.6, 11.9, "PCSK9", "Lipids"),
    ("4", 88.0, 89.0, "PKD2", "Polycystic kidney disease"),
    ("4", 100.4, 100.5, "ADH1B / ALDH2 region", "Alcohol metabolism"),
    ("5", 35.8, 36.0, "IL7R", "Multiple sclerosis"),
    ("6", 26.0, 33.0, "HLA / MHC", "Immune function — broad"),
    ("6", 161.0, 161.1, "LPA", "Lp(a) / coronary disease"),
    ("7", 117.0, 117.4, "CFTR", "Cystic fibrosis"),
    ("9", 21.9, 22.2, "9p21 / CDKN2A/B", "CAD + cancer + T2D"),
    ("9", 136.1, 136.2, "ABO", "Blood group / VTE"),
    ("11", 5.2, 5.3, "HBB", "Sickle cell / beta-thalassemia"),
    ("12", 6.0, 6.1, "VWF", "von Willebrand disease"),
    ("13", 32.8, 33.0, "BRCA2", "Breast/ovarian cancer"),
    ("14", 94.3, 94.4, "SERPINA1", "Alpha-1 antitrypsin"),
    ("16", 53.7, 54.2, "FTO", "Body weight"),
    ("17", 41.1, 41.3, "BRCA1", "Breast/ovarian cancer"),
    ("19", 45.4, 45.5, "APOE", "Alzheimer's / lipids"),
    ("20", 4.6, 4.8, "PRNP", "Prion disease"),
]


def _is_homozygous(gt: object) -> int:
    """Return 1 if homozygous (AA, GG, CC, TT, II, DD), 0 if heterozygous,
    -1 if uninterpretable (no-call, indel, missing)."""
    if gt is None:
        return -1
    s = str(gt).upper().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--", "00") or len(s) != 2:
        return -1
    a, b = s[0], s[1]
    if a not in "ACGTID" or b not in "ACGTID":
        return -1
    return 1 if a == b else 0


def _detect_roh_one_chrom(positions: np.ndarray, hom: np.ndarray,
                          min_snps: int, min_length_mb: float,
                          max_hets: int,
                          window_snps: int = WINDOW_SNPS,
                          max_gap_mb: float = MAX_GAP_MB) -> list[dict]:
    """Detect ROH on one chromosome.

    positions: SNP positions in bp (int)
    hom: 1=hom, 0=het, -1=uncalled (same length as positions)

    The het budget is enforced over a TRAILING WINDOW of ``window_snps`` called
    SNPs, so an isolated het early in a long segment does not terminate it — the
    behaviour the module docstring has always described. A whole-run counter
    instead fragments genuine IBD segments into several short runs, which both
    depresses F_ROH and moves length mass out of the long class that
    consanguinity inference depends on.

    A run also breaks across any inter-SNP gap wider than ``max_gap_mb``: an
    uncalled stretch is absence of evidence, and bridging it manufactures a
    single long run out of two unrelated homozygous stretches.
    """
    runs: list[dict] = []
    n = len(positions)
    max_gap_bp = max_gap_mb * 1e6

    def _emit(start_idx: int, end_idx: int) -> None:
        # Trim trailing non-homozygous calls so a run starts and ends on a hom.
        while end_idx > start_idx and hom[end_idx] != 1:
            end_idx -= 1
        if end_idx <= start_idx:
            return
        n_snps = end_idx - start_idx + 1
        length_mb = int(positions[end_idx] - positions[start_idx]) / 1e6
        if n_snps >= min_snps and length_mb >= min_length_mb:
            runs.append({
                "start_bp": int(positions[start_idx]),
                "end_bp": int(positions[end_idx]),
                "length_mb": round(length_mb, 3),
                "n_snps": int(n_snps),
            })

    i = 0
    while i < n:
        if hom[i] != 1:
            i += 1
            continue
        start_idx = i
        het_idx: list[int] = []      # positions of hets inside the current run
        j = i
        while j < n:
            # Gap guard — never extend across an uncalled stretch.
            if j > start_idx and (positions[j] - positions[j - 1]) > max_gap_bp:
                break
            if hom[j] == 0:
                het_idx.append(j)
                # Enforce max_hets over the trailing window only.
                window_start = max(start_idx, j - window_snps + 1)
                if sum(1 for h in het_idx if h >= window_start) > max_hets:
                    break
            j += 1
        _emit(start_idx, j - 1)
        # Resume after the run; always advance to guarantee termination.
        i = max(j, start_idx + 1)
    return runs


def detect_roh(snps_df: pd.DataFrame,
               min_snps: int = 50,
               min_length_mb: float = 1.0,
               max_hets: int = 1) -> dict:
    """Detect ROH genome-wide. Returns aggregated results with classification."""
    if snps_df.empty or "chrom" not in snps_df.columns:
        # Schema-consistent empty result: MUST carry the same keys as the full
        # return below, or downstream (pipeline log line, build_roh_html) hits a
        # KeyError and takes the whole report down with it.
        return {"runs": [], "f_roh": 0.0, "n_runs": 0, "total_roh_mb": 0.0,
                "long_roh_mb": 0.0, "f_roh_long": 0.0,
                "short": [], "medium": [], "long": [],
                "n_short": 0, "n_medium": 0, "n_long": 0,
                "population_context": ("No autosomal genotype data available to scan "
                                       "for runs of homozygosity."),
                "context_tier": "unavailable"}

    runs_all: list[dict] = []
    for c in CHROM_LENGTHS_MB:
        sub = snps_df[snps_df["chrom"].astype(str) == c]
        if len(sub) < min_snps:
            continue
        sub = sub.sort_values("pos")
        positions = sub["pos"].to_numpy(dtype=np.int64)
        hom = np.fromiter((_is_homozygous(g) for g in sub["genotype"]),
                          dtype=np.int8, count=len(sub))
        chrom_runs = _detect_roh_one_chrom(positions, hom, min_snps, min_length_mb, max_hets)
        for r in chrom_runs:
            r["chrom"] = c
            r = _annotate_roh_with_genes(r)
            runs_all.append(r)

    # Classify
    short = [r for r in runs_all if r["length_mb"] < 2.0]
    medium = [r for r in runs_all if 2.0 <= r["length_mb"] < LONG_ROH_MB]
    long_ = [r for r in runs_all if r["length_mb"] >= LONG_ROH_MB]

    total_mb = sum(r["length_mb"] for r in runs_all)
    f_roh = total_mb / AUTOSOMAL_LENGTH_MB
    # Only long ROH track recent pedigree relatedness (Kirin 2010): short and
    # medium runs are dominated by ancient shared ancestry, so blending them in
    # labels founder-population background as consanguinity.
    long_mb = sum(r["length_mb"] for r in long_)
    f_roh_long = long_mb / AUTOSOMAL_LENGTH_MB

    # Population context narrative — keyed on F_ROH_LONG, not total F_ROH.
    _founder_note = (" Total F_ROH here is "
                     f"{f_roh:.4f}; the short and medium runs making up the "
                     "difference reflect ancient shared ancestry (common in "
                     "founder populations such as Ashkenazi Jewish, Finnish, "
                     "Sardinian, French Canadian or Old Order Amish) rather "
                     "than recent parental relatedness.")
    if f_roh_long <= 0.0:
        context = ("No long (>8 Mb) runs of homozygosity were detected, so there is no "
                   "evidence of recent parental relatedness." + _founder_note)
        context_tier = "no_recent_relatedness"
    elif f_roh_long < 0.02:
        context = ("A small amount of long-ROH burden can reflect distant relatedness or a "
                   "strong founder-population background. It is not on its own evidence of "
                   "close parental relatedness." + _founder_note)
        context_tier = "distant_or_founder_background"
    elif f_roh_long < 0.05:
        context = ("This long-ROH burden is consistent with second-cousin parental "
                   "relatedness. Many cultures practice cousin marriage and this is "
                   "informational, not judgmental. Clinically, offspring of related parents "
                   "have elevated rare-recessive-disease risk; pre-pregnancy expanded carrier "
                   "screening can be valuable.")
        context_tier = "second_cousin"
    else:
        context = ("This long-ROH burden is consistent with first-cousin or closer parental "
                   "relatedness. Discuss with a genetic counsellor — expanded carrier "
                   "screening is particularly valuable for family-planning decisions.")
        context_tier = "first_cousin_or_closer"

    return {
        "runs": runs_all,
        "n_runs": len(runs_all),
        "total_roh_mb": round(total_mb, 2),
        "f_roh": round(f_roh, 5),
        "long_roh_mb": round(long_mb, 2),
        "f_roh_long": round(f_roh_long, 5),
        "short": short,
        "medium": medium,
        "long": long_,
        "n_short": len(short),
        "n_medium": len(medium),
        "n_long": len(long_),
        "population_context": context,
        "context_tier": context_tier,
        "parameters": {
            "min_snps": min_snps,
            "min_length_mb": min_length_mb,
            "max_hets": max_hets,
            "window_snps": WINDOW_SNPS,
            "max_gap_mb": MAX_GAP_MB,
        },
    }


def _annotate_roh_with_genes(roh: dict) -> dict:
    """Note any disease-relevant gene regions overlapping the ROH."""
    chrom = roh["chrom"]
    start_mb = roh["start_bp"] / 1e6
    end_mb = roh["end_bp"] / 1e6
    overlaps = []
    for c, gs, ge, gene, cond in DISEASE_GENE_REGIONS:
        if c != chrom:
            continue
        if gs < end_mb and ge > start_mb:
            overlaps.append({"gene": gene, "condition": cond,
                             "start_mb": gs, "end_mb": ge})
    if overlaps:
        roh["disease_genes"] = overlaps
    return roh


# ─── SVG ideogram rendering ──────────────────────────────────────────────────

def render_ideogram_svg(roh_runs: list[dict], width: int = 720) -> str:
    """Render a chromosome ideogram SVG with ROH bands highlighted."""
    chroms = list(CHROM_LENGTHS_MB.keys())
    n_chroms = len(chroms)
    chrom_height = 14
    chrom_gap = 8
    left_margin = 36
    right_margin = 20
    top_margin = 24
    bottom_margin = 20

    max_length = max(CHROM_LENGTHS_MB.values())
    scale = (width - left_margin - right_margin) / max_length

    height = top_margin + n_chroms * (chrom_height + chrom_gap) + bottom_margin

    rects = []
    labels = []
    for i, c in enumerate(chroms):
        y = top_margin + i * (chrom_height + chrom_gap)
        chrom_len = CHROM_LENGTHS_MB[c]
        bar_width = chrom_len * scale
        # Chromosome backbone
        rects.append(
            f'<rect x="{left_margin}" y="{y}" width="{bar_width:.1f}" '
            f'height="{chrom_height}" rx="6" fill="#e6e6e6" stroke="#bbb"/>'
        )
        labels.append(
            f'<text x="{left_margin-6}" y="{y + chrom_height - 3}" '
            f'text-anchor="end" font-size="11" font-family="sans-serif" fill="#888">{c}</text>'
        )
        # ROH overlays on this chromosome
        for run in roh_runs:
            if run["chrom"] != c:
                continue
            rx = left_margin + (run["start_bp"] / 1e6) * scale
            rw = max((run["length_mb"]) * scale, 1.0)
            # Colour by ROH length tier
            if run["length_mb"] >= 8:
                fill = "#f85149"  # long — recent relatedness
            elif run["length_mb"] >= 2:
                fill = "#d29922"  # medium
            else:
                fill = "#58a6ff"  # short
            rects.append(
                f'<rect x="{rx:.1f}" y="{y}" width="{rw:.1f}" '
                f'height="{chrom_height}" fill="{fill}" fill-opacity="0.85"/>'
            )

    legend_y = top_margin - 14
    legend = (
        f'<rect x="{left_margin}" y="{legend_y}" width="10" height="10" fill="#58a6ff"/>'
        f'<text x="{left_margin+14}" y="{legend_y+9}" font-size="10" font-family="sans-serif" '
        f'fill="#666">Short (1-2 Mb)</text>'
        f'<rect x="{left_margin+110}" y="{legend_y}" width="10" height="10" fill="#d29922"/>'
        f'<text x="{left_margin+124}" y="{legend_y+9}" font-size="10" font-family="sans-serif" '
        f'fill="#666">Medium (2-8 Mb)</text>'
        f'<rect x="{left_margin+230}" y="{legend_y}" width="10" height="10" fill="#f85149"/>'
        f'<text x="{left_margin+244}" y="{legend_y+9}" font-size="10" font-family="sans-serif" '
        f'fill="#666">Long (>8 Mb, recent relatedness)</text>'
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        f'{legend}'
        + "".join(rects) + "".join(labels) +
        '</svg>'
    )
    return svg
