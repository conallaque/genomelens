"""
Local Ancestry Deconvolution (Chromosome Painting)
==================================================

Estimates the ancestral origin of segments of the genome ("chromosome
painting") instead of just global ancestry proportions. The state of the
art (RFMix, ChromoPainter) uses HMMs on phased haplotypes against panels
of thousands of phased reference samples. Without a full reference panel
this module implements a simplified per-window approach:

  1. Sweep each autosome in fixed-size windows (default 5 Mb).
  2. Within each window, collect ancestry-informative markers (AIMs)
     with known per-superpopulation allele frequencies.
  3. For each window, compute the log-likelihood of the user's genotypes
     under each of five 1000G superpopulations (EUR, AFR, EAS, SAS, AMR).
  4. Assign the window to its most likely population (with confidence
     from the log-likelihood gap).
  5. Flag windows whose call deviates strongly from the genome-wide
     ancestry estimate — these are the candidate ancestry-switch
     segments / admixture relics.
  6. Render an SVG "chromosome painting" using superpopulation colours.

Limitations (called out in the report):
  * The AIM panel is small (~14 markers); window calls in regions with
    few AIMs will be low-confidence or unassigned.
  * Without phasing we cannot resolve maternal vs paternal segments.
  * For real research-grade local ancestry, run RFMix against a 1000G
    reference panel.
"""

from __future__ import annotations

import math

import pandas as pd

# Reuse the AIMs panel + superpop allele frequencies from ancestry_pca.py
try:
    from ancestry_pca import AIMS_PRIORS, SUPERPOP_LONG, SUPERPOPS
except Exception:
    AIMS_PRIORS = {}
    SUPERPOPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
    SUPERPOP_LONG = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
                      "SAS": "South Asian", "AMR": "Admixed American"}


# Superpopulation colours for chromosome painting
SUPERPOP_COLORS = {
    "EUR": "#3b82f6",   # blue
    "AFR": "#f59e0b",   # orange
    "EAS": "#10b981",   # green
    "SAS": "#a855f7",   # purple
    "AMR": "#ef4444",   # red
    "UNK": "#cbd5e1",   # grey — unassigned (no AIM in window)
}

# Approximate chromosome lengths (Mb) for SVG sizing
CHROM_LENGTHS_MB = {
    "1": 249.3, "2": 243.2, "3": 198.0, "4": 191.2, "5": 180.9, "6": 171.1,
    "7": 159.1, "8": 146.4, "9": 141.2, "10": 135.5, "11": 135.0, "12": 133.9,
    "13": 115.2, "14": 107.3, "15": 102.5, "16": 90.4, "17": 81.2, "18": 78.1,
    "19": 59.1, "20": 63.0, "21": 48.1, "22": 51.3,
}


# Disease/trait genes by approximate chromosomal position — used to annotate
# interesting ancestry-discordant segments.
GENE_ANNOTATIONS: list[tuple] = [
    ("4", 100, 101, "ADH1B / ALDH2 region — alcohol metabolism"),
    ("6", 26, 33, "HLA / MHC — immune function"),
    ("7", 117, 118, "CFTR — cystic fibrosis"),
    ("9", 21, 23, "9p21 — CAD/T2D/cancer"),
    ("11", 5, 6, "HBB — sickle cell / β-thal"),
    ("15", 28, 29, "OCA2 / HERC2 — eye colour"),
    ("17", 41, 42, "BRCA1"),
    ("13", 32, 33, "BRCA2"),
    ("19", 45, 46, "APOE — Alzheimer's / lipids"),
]


def _dose(genotype: object, effect_allele: str) -> int | None:
    if genotype is None:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--") or len(gt) != 2:
        return None
    return gt.count(effect_allele)


def _binomial_loglik(dose: int, af: float) -> float:
    af = min(max(af, 1e-4), 1 - 1e-4)
    if dose == 0:
        return 2 * math.log(1 - af)
    if dose == 1:
        return math.log(2 * af * (1 - af))
    return 2 * math.log(af)


def _aim_by_chrom_pos(snps_df: pd.DataFrame) -> dict[str, list[dict]]:
    """Group the AIM panel by chromosome and attach the user's genotype."""
    by_chrom: dict[str, list[dict]] = {}
    for rsid, info in AIMS_PRIORS.items():
        if rsid not in snps_df.index:
            continue
        row = snps_df.loc[rsid]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        chrom = str(row.get("chrom"))
        pos = row.get("pos")
        if not chrom or pos is None:
            continue
        try:
            pos_mb = int(pos) / 1e6
        except Exception:
            continue
        dose = _dose(row.get("genotype"), info["effect_allele"])
        if dose is None:
            continue
        by_chrom.setdefault(chrom, []).append({
            "rsid": rsid, "pos_mb": pos_mb, "dose": dose, "info": info,
        })
    # Sort per-chrom
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: x["pos_mb"])
    return by_chrom


def _annotate_genes(chrom: str, start_mb: float, end_mb: float) -> list[str]:
    notes: list[str] = []
    for c, gs, ge, name in GENE_ANNOTATIONS:
        if c == chrom and gs < end_mb and ge > start_mb:
            notes.append(name)
    return notes


def analyze_local_ancestry(snps_df: pd.DataFrame,
                            window_mb: float = 5.0,
                            global_proportions: dict[str, float] | None = None
                            ) -> dict:
    """Compute per-window local ancestry calls.

    Returns dict with windows (list of window dicts), summary stats, and
    deviant segments (windows whose call differs from the global mode).
    """
    aim_by_chrom = _aim_by_chrom_pos(snps_df)
    if not aim_by_chrom:
        return {
            "available": False,
            "reason": "No AIM SNPs found on this chip — local ancestry cannot be estimated.",
            "windows": [],
        }

    all_windows: list[dict] = []
    for chrom, length in CHROM_LENGTHS_MB.items():
        aims = aim_by_chrom.get(chrom, [])
        # Window iteration
        n_windows = max(1, int(math.ceil(length / window_mb)))
        for w in range(n_windows):
            wstart = w * window_mb
            wend = min((w + 1) * window_mb, length)
            in_window = [a for a in aims if wstart <= a["pos_mb"] < wend]
            window = {
                "chrom": chrom,
                "start_mb": wstart,
                "end_mb": wend,
                "n_aims": len(in_window),
                "call": "UNK",
                "confidence": "n/a",
                "log_liks": {sp: 0.0 for sp in SUPERPOPS},
            }
            if in_window:
                loglik = {sp: 0.0 for sp in SUPERPOPS}
                for a in in_window:
                    for sp in SUPERPOPS:
                        af = a["info"][sp]
                        loglik[sp] += _binomial_loglik(a["dose"], af)
                # Find best superpop + gap to second-best
                ordered = sorted(loglik.items(), key=lambda x: -x[1])
                best_pop, best_ll = ordered[0]
                second_ll = ordered[1][1]
                gap = best_ll - second_ll
                # Confidence based on gap and number of AIMs
                if len(in_window) >= 2 and gap > 2.0:
                    conf = "high"
                elif gap > 1.0:
                    conf = "moderate"
                else:
                    conf = "low"
                window["call"] = best_pop
                window["confidence"] = conf
                window["log_liks"] = {sp: round(loglik[sp], 2) for sp in SUPERPOPS}
                window["gap"] = round(gap, 2)
            all_windows.append(window)

    # Summary stats
    by_call: dict[str, int] = {}
    for w in all_windows:
        if w["call"] != "UNK":
            by_call[w["call"]] = by_call.get(w["call"], 0) + 1
    total_called = sum(by_call.values())

    # Deviant segments: those whose call is not the genome-wide mode (and confident)
    deviant: list[dict] = []
    mode_pop = max(by_call.items(), key=lambda x: x[1])[0] if by_call else None
    for w in all_windows:
        if (w["call"] != mode_pop and w["call"] != "UNK"
                and w["confidence"] in ("high", "moderate")):
            w_ann = dict(w)
            w_ann["genes"] = _annotate_genes(w["chrom"], w["start_mb"], w["end_mb"])
            deviant.append(w_ann)

    return {
        "available": True,
        "window_mb": window_mb,
        "windows": all_windows,
        "n_windows_called": total_called,
        "n_windows_total": len(all_windows),
        "by_call": by_call,
        "mode_population": mode_pop,
        "deviant_segments": deviant,
        "n_deviant": len(deviant),
        "limitations": (
            "Based on a small AIM panel (~14 markers). Most windows will be unassigned "
            "(insufficient AIMs in that window). Without phased haplotypes we cannot "
            "separate maternal vs paternal chromosomes. For research-grade local ancestry, "
            "run RFMix against a phased 1000 Genomes reference panel."
        ),
    }


# ─── SVG chromosome painting ─────────────────────────────────────────────────

def render_chromosome_painting_svg(local_ancestry: dict, width: int = 760) -> str:
    """Render a chromosome painting SVG with windows coloured by local ancestry call."""
    chroms = list(CHROM_LENGTHS_MB.keys())
    n = len(chroms)
    bar_h = 14
    gap = 8
    left = 40
    right = 16
    top = 38
    bottom = 24

    max_len = max(CHROM_LENGTHS_MB.values())
    scale = (width - left - right) / max_len

    height = top + n * (bar_h + gap) + bottom

    parts = []
    # Legend
    legend_y = 12
    x_off = left
    for sp in SUPERPOPS + ["UNK"]:
        col = SUPERPOP_COLORS.get(sp, "#888")
        label = "Unassigned (no AIM)" if sp == "UNK" else f"{sp} ({SUPERPOP_LONG.get(sp, sp)})"
        parts.append(
            f'<rect x="{x_off}" y="{legend_y}" width="10" height="10" fill="{col}"/>'
            f'<text x="{x_off+13}" y="{legend_y+9}" font-size="10" '
            f'font-family="sans-serif" fill="#666">{label}</text>'
        )
        x_off += 120

    windows_by_chrom: dict[str, list[dict]] = {}
    for w in local_ancestry.get("windows", []):
        windows_by_chrom.setdefault(w["chrom"], []).append(w)

    for i, c in enumerate(chroms):
        y = top + i * (bar_h + gap)
        chrom_len = CHROM_LENGTHS_MB[c]
        bar_w = chrom_len * scale
        # Backbone
        parts.append(
            f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="{bar_h}" '
            f'rx="6" fill="#eee" stroke="#ccc"/>'
        )
        # Label
        parts.append(
            f'<text x="{left-6}" y="{y + bar_h - 3}" text-anchor="end" '
            f'font-size="11" font-family="sans-serif" fill="#888">{c}</text>'
        )
        # Window overlays
        for w in windows_by_chrom.get(c, []):
            wx = left + w["start_mb"] * scale
            ww = max((w["end_mb"] - w["start_mb"]) * scale, 0.5)
            col = SUPERPOP_COLORS.get(w["call"], "#888")
            opacity = "0.95" if w.get("confidence") == "high" else (
                "0.7" if w.get("confidence") == "moderate" else "0.4")
            parts.append(
                f'<rect x="{wx:.2f}" y="{y}" width="{ww:.2f}" height="{bar_h}" '
                f'fill="{col}" fill-opacity="{opacity}"/>'
            )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        + "".join(parts) +
        '</svg>'
    )
