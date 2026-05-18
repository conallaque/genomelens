"""
Runs of Homozygosity (ROH) Analysis
===================================

Detects contiguous stretches of homozygous genotypes across the autosomes.
Long ROH reveal population history (founder effects, drift, isolation)
and — at the longest scales — recent parental relatedness (consanguinity).

Algorithm (PLINK-like sliding window):
    1. For each autosome, walk the SNPs in physical order.
    2. Track current run length: extend when homozygous, allow up to
       `max_hets` heterozygotes per window.
    3. Emit a run when it reaches `min_snps` SNPs AND `min_length_mb` Mb.

Output:
    * Per-run records (chr, start, end, n_snps, length_mb).
    * Classification: short (1-2 Mb), medium (2-8 Mb), long (>8 Mb).
    * F_ROH coefficient = total ROH length / autosomal genome length.
    * SVG ideogram visualisation of ROH locations.

Population context (illustrative, not used for hard cutoffs):
    * F_ROH ~ 0.005 — outbred general population
    * F_ROH 0.01-0.03 — founder populations (Ashkenazi, Finnish, Sardinian)
    * F_ROH 0.03-0.08 — first/second-cousin parents likely
    * F_ROH > 0.08 — closer-than-first-cousin parents

Framing: consanguinity is culturally normal in many populations. This is
informational, not judgmental. The clinically relevant consequence is
elevated rare-recessive-disease risk for offspring of related parents.
"""

from __future__ import annotations

from typing import Dict, List, Tuple
import pandas as pd
import numpy as np


# Approximate autosomal genome length in GRCh37/38 (Mb)
AUTOSOMAL_LENGTH_MB = 2881.0
# Approximate per-chromosome lengths (Mb) for ideogram rendering
CHROM_LENGTHS_MB = {
    "1": 249.3, "2": 243.2, "3": 198.0, "4": 191.2, "5": 180.9, "6": 171.1,
    "7": 159.1, "8": 146.4, "9": 141.2, "10": 135.5, "11": 135.0, "12": 133.9,
    "13": 115.2, "14": 107.3, "15": 102.5, "16": 90.4, "17": 81.2, "18": 78.1,
    "19": 59.1, "20": 63.0, "21": 48.1, "22": 51.3,
}

# Disease-relevant gene regions (chr, start_Mb, end_Mb, gene, condition)
# Light-touch annotation so we can flag ROH overlapping known disease genes
DISEASE_GENE_REGIONS: List[Tuple] = [
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
                          max_hets: int) -> List[Dict]:
    """Detect ROH on one chromosome.

    positions: SNP positions in bp (int)
    hom: 1=hom, 0=het, -1=uncalled (same length as positions)
    """
    runs: List[Dict] = []
    n = len(positions)
    i = 0
    while i < n:
        if hom[i] != 1:
            i += 1
            continue
        # Start of a candidate run
        start_idx = i
        hets_in_window = 0
        j = i
        while j < n:
            if hom[j] == 1:
                j += 1
            elif hom[j] == 0:
                hets_in_window += 1
                if hets_in_window > max_hets:
                    break
                j += 1
            else:  # uncalled — skip without counting
                j += 1
        end_idx = j - 1
        # Trim trailing het/uncalled
        while end_idx > start_idx and hom[end_idx] != 1:
            end_idx -= 1

        n_snps = end_idx - start_idx + 1
        length_bp = int(positions[end_idx] - positions[start_idx])
        length_mb = length_bp / 1e6
        if n_snps >= min_snps and length_mb >= min_length_mb:
            runs.append({
                "start_bp": int(positions[start_idx]),
                "end_bp": int(positions[end_idx]),
                "length_mb": round(length_mb, 3),
                "n_snps": int(n_snps),
            })
        i = end_idx + 1
    return runs


def detect_roh(snps_df: pd.DataFrame,
               min_snps: int = 50,
               min_length_mb: float = 1.0,
               max_hets: int = 1) -> Dict:
    """Detect ROH genome-wide. Returns aggregated results with classification."""
    if snps_df.empty or "chrom" not in snps_df.columns:
        return {"runs": [], "f_roh": 0.0, "n_runs": 0, "total_roh_mb": 0.0,
                "short": [], "medium": [], "long": []}

    runs_all: List[Dict] = []
    for c in CHROM_LENGTHS_MB.keys():
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
    medium = [r for r in runs_all if 2.0 <= r["length_mb"] < 8.0]
    long_ = [r for r in runs_all if r["length_mb"] >= 8.0]

    total_mb = sum(r["length_mb"] for r in runs_all)
    f_roh = total_mb / AUTOSOMAL_LENGTH_MB

    # Population context narrative
    if f_roh < 0.005:
        context = ("F_ROH below ~0.005 is typical for an outbred individual with no recent "
                   "parental relatedness.")
        context_tier = "outbred"
    elif f_roh < 0.015:
        context = ("F_ROH in the 0.005-0.015 range is common in modern populations and may "
                   "reflect ancient population bottlenecks or founder-population background "
                   "(e.g. Ashkenazi Jewish, Finnish, Sardinian, French Canadian, Old Order "
                   "Amish, some Middle Eastern populations).")
        context_tier = "founder_background"
    elif f_roh < 0.04:
        context = ("F_ROH 0.015-0.04 suggests either strong founder-population background or "
                   "more distant consanguinity (second cousins or beyond). Common in some "
                   "endogamous communities.")
        context_tier = "endogamous_or_distant_consanguinity"
    elif f_roh < 0.08:
        context = ("F_ROH 0.04-0.08 is consistent with first- or second-cousin parental "
                   "relatedness. Many cultures practice cousin marriage and this is "
                   "informational, not judgmental. Clinically, offspring of related parents "
                   "have elevated rare-recessive-disease risk; pre-pregnancy expanded carrier "
                   "screening can be valuable.")
        context_tier = "first_second_cousin"
    else:
        context = ("F_ROH above 0.08 is consistent with closer-than-first-cousin parental "
                   "relatedness. Discuss with a genetic counsellor — expanded carrier "
                   "screening is particularly valuable for family-planning decisions.")
        context_tier = "very_close_relatedness"

    return {
        "runs": runs_all,
        "n_runs": len(runs_all),
        "total_roh_mb": round(total_mb, 2),
        "f_roh": round(f_roh, 5),
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
        },
    }


def _annotate_roh_with_genes(roh: Dict) -> Dict:
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

def render_ideogram_svg(roh_runs: List[Dict], width: int = 720) -> str:
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
