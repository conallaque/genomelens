"""
PCA-Based Ancestry Estimation
=============================

Projects the user's genotype onto a principal-components space defined by
1000 Genomes superpopulation reference samples (EUR, AFR, EAS, SAS, AMR)
and estimates ancestry proportions.

Two operating modes:

  1. Full PCA (requires 1000G genotype data at the AIMs panel):
       a. Load 1000G AIM genotypes + superpop labels.
       b. Fit a PCA on the reference.
       c. Project the user's genotype onto the same PC axes.
       d. Run a soft k-NN classifier in PC space against the reference
          superpopulations to estimate proportions.
       e. Render a PNG scatter (user point overlaid on reference clouds).

  2. AIMs-only heuristic (fallback when only the curated AIM list is present):
       Uses the curated set of high-information markers and a simple Bayesian
       likelihood-by-superpopulation estimate. Approximate but informative
       even from chip data alone.

This is rough compared to commercial ancestry products (which use much
larger reference panels and proprietary algorithms). It is transparent and
local.
"""

from __future__ import annotations

import base64
import io
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


SCRIPT_DIR = Path(__file__).parent
REF_DIR = SCRIPT_DIR / "reference" / "ancestry"


def _log(msg: str) -> None:
    print(f"[ancestry] {msg}", flush=True)


# ── AIMs-only heuristic fallback ──────────────────────────────────────────────
# Each AIM has reference allele frequencies in five 1000G superpopulations.
# When 1000G genotype data isn't installed, we use these to compute a
# log-likelihood of the user's genotype under each superpop.
#
# Frequencies are derived from 1000 Genomes Phase 3 (EUR, AFR, EAS, SAS, AMR).
# Values are the frequency of the listed `effect_allele`.
#
# `ld_block`: markers sharing a block are in strong linkage disequilibrium (or
# are mirror tags of the same locus) and carry essentially the same ancestry
# signal. The heuristic counts only ONE marker per block so a single locus
# (e.g. eye-colour HERC2) is not double- or triple-counted. Treating correlated
# markers as independent in a naive-Bayes likelihood was the main driver of
# spurious AMR ("Admixed American") calls for Southern European samples.
#
# Caveat on AMR: the 1000G AMR group is recently admixed (predominantly
# European + Native American + African), so these are pooled frequencies that
# do not describe a single panmictic, Hardy-Weinberg population. Under a
# naive likelihood AMR tends to act as a "central" distribution that is never
# penalised hard on any one marker, which inflates its score for intermediate
# genotypes. We therefore (a) never let this heuristic exceed "moderate"
# confidence and (b) downgrade to "low" / flag the call as ambiguous whenever
# the top population's evidence margin over the runner-up is small.

AIMS_PRIORS: Dict[str, Dict] = {
    "rs3827760":  {"effect_allele": "G", "EUR": 0.005, "AFR": 0.003, "EAS": 0.93, "SAS": 0.005, "AMR": 0.35, "gene": "EDAR", "ld_block": "EDAR"},
    "rs1426654":  {"effect_allele": "A", "EUR": 1.00,  "AFR": 0.07,  "EAS": 0.005, "SAS": 0.55,  "AMR": 0.45, "gene": "SLC24A5"},
    "rs16891982": {"effect_allele": "C", "EUR": 0.97,  "AFR": 0.05,  "EAS": 0.005, "SAS": 0.32,  "AMR": 0.45, "gene": "SLC45A2"},
    "rs12913832": {"effect_allele": "G", "EUR": 0.63,  "AFR": 0.02,  "EAS": 0.005, "SAS": 0.05,  "AMR": 0.22, "gene": "HERC2", "ld_block": "HERC2_eye"},
    "rs1805007":  {"effect_allele": "T", "EUR": 0.09,  "AFR": 0.005, "EAS": 0.001, "SAS": 0.005, "AMR": 0.02, "gene": "MC1R"},
    "rs1042602":  {"effect_allele": "A", "EUR": 0.30,  "AFR": 0.001, "EAS": 0.001, "SAS": 0.04,  "AMR": 0.10, "gene": "TYR"},
    "rs17822931": {"effect_allele": "T", "EUR": 0.04,  "AFR": 0.001, "EAS": 0.85,  "SAS": 0.05,  "AMR": 0.30, "gene": "ABCC11"},
    "rs1129038":  {"effect_allele": "T", "EUR": 0.37,  "AFR": 0.98,  "EAS": 0.99,  "SAS": 0.95,  "AMR": 0.78, "gene": "HERC2", "ld_block": "HERC2_eye"},
    "rs2402130":  {"effect_allele": "G", "EUR": 0.05,  "AFR": 0.01,  "EAS": 0.62,  "SAS": 0.10,  "AMR": 0.30, "gene": "EDAR", "ld_block": "EDAR"},
    "rs671":      {"effect_allele": "A", "EUR": 0.002, "AFR": 0.002, "EAS": 0.21,  "SAS": 0.005, "AMR": 0.02, "gene": "ALDH2"},
    "rs1229984":  {"effect_allele": "A", "EUR": 0.03,  "AFR": 0.005, "EAS": 0.70,  "SAS": 0.07,  "AMR": 0.18, "gene": "ADH1B"},
    "rs2814778":  {"effect_allele": "C", "EUR": 0.005, "AFR": 0.99,  "EAS": 0.005, "SAS": 0.005, "AMR": 0.05, "gene": "DARC/ACKR1"},
    "rs4988235":  {"effect_allele": "T", "EUR": 0.51,  "AFR": 0.02,  "EAS": 0.001, "SAS": 0.18,  "AMR": 0.30, "gene": "LCT", "ld_block": "LCT_MCM6"},
    "rs182549":   {"effect_allele": "T", "EUR": 0.51,  "AFR": 0.02,  "EAS": 0.001, "SAS": 0.18,  "AMR": 0.30, "gene": "MCM6", "ld_block": "LCT_MCM6"},
}

# Number of independent ancestry signals in the panel (distinct LD blocks +
# standalone markers) — the realistic "expected" denominator for coverage.
N_INDEPENDENT_AIMS = len({p.get("ld_block", rsid) for rsid, p in AIMS_PRIORS.items()})

SUPERPOPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
SUPERPOP_LONG = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Admixed American",
}


def _dosage(genotype: object, effect_allele: str) -> Optional[int]:
    if genotype is None:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--") or len(gt) != 2:
        return None
    return gt.count(effect_allele)


def _binomial_loglik(dose: int, af: float) -> float:
    """Log-likelihood of observing `dose` effect alleles under Hardy-Weinberg
    with effect-allele frequency `af`. Two independent Bernoulli draws."""
    af = min(max(af, 1e-4), 1 - 1e-4)
    if dose == 0:
        return 2 * np.log(1 - af)
    if dose == 1:
        return np.log(2 * af * (1 - af))
    return 2 * np.log(af)


# Minimum independent AIMs below which the heuristic is not trustworthy at all.
MIN_AIMS_FOR_CALL = 4
# Log-likelihood margin (nats) of the top population over the runner-up below
# which the two cannot be reliably distinguished on this small panel.
# ~2.3 nats ≈ a 10:1 likelihood ratio.
AMBIGUOUS_MARGIN_NATS = 2.3


def estimate_ancestry_heuristic(snps_df: pd.DataFrame) -> Dict:
    """Rough single-population affinity from the curated AIM panel.

    This is NOT an admixture estimate. It computes, under a uniform prior, the
    posterior probability that the genotype best matches each *single* 1000G
    superpopulation, then reports the best match with a confidence that is
    driven by (a) how many independent AIMs were typed and (b) how decisively
    the top population beats the runner-up. LD-correlated markers are counted
    once per block (see AIMS_PRIORS).
    """
    loglik = {sp: 0.0 for sp in SUPERPOPS}
    used_aims: List[Dict] = []
    seen_blocks: set = set()
    n_redundant = 0
    for rsid, info in AIMS_PRIORS.items():
        if rsid not in snps_df.index:
            continue
        row = snps_df.loc[rsid]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        dose = _dosage(row.get("genotype"), info["effect_allele"])
        if dose is None:
            continue
        block = info.get("ld_block", rsid)
        if block in seen_blocks:
            # Correlated with a marker we already used — record it for the
            # report but do NOT add its (redundant) likelihood again.
            n_redundant += 1
            used_aims.append({
                "rsid": rsid, "gene": info["gene"],
                "genotype": str(row.get("genotype")).upper(),
                "effect_allele": info["effect_allele"], "dosage": dose,
                "ld_block": block, "counted": False,
            })
            continue
        seen_blocks.add(block)
        for sp in SUPERPOPS:
            loglik[sp] += _binomial_loglik(dose, info[sp])
        used_aims.append({
            "rsid": rsid, "gene": info["gene"],
            "genotype": str(row.get("genotype")).upper(),
            "effect_allele": info["effect_allele"], "dosage": dose,
            "ld_block": block, "counted": True,
        })

    n_independent = len(seen_blocks)
    if n_independent == 0:
        return {
            "available": False,
            "reason": "No AIM SNPs detected on chip — cannot estimate ancestry.",
            "n_aims_used": 0,
            "n_aims_independent": 0,
            "n_aims_expected": N_INDEPENDENT_AIMS,
            "confidence": "none",
        }

    # Convert log-likelihoods to probabilities (uniform prior). These are
    # relative affinities to a *single* population, not admixture fractions.
    mx = max(loglik.values())
    weights = {sp: np.exp(ll - mx) for sp, ll in loglik.items()}
    total = sum(weights.values())
    proportions = {sp: weights[sp] / total for sp in SUPERPOPS}

    sorted_props = sorted(proportions.items(), key=lambda x: -x[1])
    primary = sorted_props[0][0]

    # Evidence margin: difference in log-likelihood between best and 2nd-best
    # population. This — not the softmax probability — is the honest measure of
    # how confident the call is, because softmax exaggerates small gaps on a
    # tiny panel.
    sorted_ll = sorted(loglik.values(), reverse=True)
    margin_nats = float(sorted_ll[0] - sorted_ll[1]) if len(sorted_ll) > 1 else float("inf")
    runner_up = sorted_props[1][0] if len(sorted_props) > 1 else None

    ambiguous = margin_nats < AMBIGUOUS_MARGIN_NATS
    # Confidence is CAPPED at "moderate": a ~10-marker pigmentation/diet panel
    # can never give a high-confidence ancestry call.
    if n_independent < MIN_AIMS_FOR_CALL:
        confidence = "low"
        confidence_note = (
            f"Only {n_independent} independent AIM(s) typed (min "
            f"{MIN_AIMS_FOR_CALL} for any call) — estimate is unreliable."
        )
    elif ambiguous:
        confidence = "low"
        confidence_note = (
            f"{SUPERPOP_LONG[primary]} and {SUPERPOP_LONG[runner_up]} are nearly "
            f"equally likely on this panel (evidence margin "
            f"{margin_nats:.1f} nats) — they cannot be reliably distinguished."
        )
    else:
        confidence = "moderate"
        confidence_note = (
            f"Best single-population match on {n_independent} independent AIMs. "
            "This heuristic panel is too small for high confidence or for "
            "admixture proportions."
        )

    return {
        "available": True,
        "method": ("AIMs-only heuristic (best single-population match by "
                   "likelihood under 1000G superpop allele-frequency priors; "
                   "LD-correlated markers counted once)"),
        "is_admixture_estimate": False,
        "n_aims_used": len(used_aims),
        "n_aims_independent": n_independent,
        "n_aims_redundant": n_redundant,
        "n_aims_expected": N_INDEPENDENT_AIMS,
        "proportions": proportions,
        "sorted_proportions": sorted_props,
        "primary_population": primary,
        "runner_up_population": runner_up,
        "evidence_margin_nats": round(margin_nats, 2),
        "ambiguous": ambiguous,
        "confidence": confidence,
        "confidence_note": confidence_note,
        "used_aims": used_aims,
    }


# ── Full PCA mode (requires 1000G genotype data) ──────────────────────────────
def _try_full_pca(snps_df: pd.DataFrame) -> Optional[Dict]:
    """Attempts a full PCA-based ancestry projection. Looks for the 1000G
    reference genotype matrix at reference/ancestry/kgp_aims.npz.

    The expected file layout (created by user-supplied setup tooling):
        kgp_aims.npz: contains
            X        — (n_samples, n_aims) int8 dosage matrix
            samples  — (n_samples,) sample IDs
            superpop — (n_samples,) superpop labels ('EUR','AFR','EAS','SAS','AMR')
            rsids    — (n_aims,) rsID list (matched to X columns)
            effect_alleles — (n_aims,) effect allele characters

    Returns None if data not available — caller falls back to heuristic.
    """
    npz_path = REF_DIR / "kgp_aims.npz"
    if not npz_path.exists():
        return None
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        _log("scikit-learn not installed; full PCA disabled. `pip install scikit-learn`")
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        _log("matplotlib not installed; PCA plot will be skipped")
        plt = None

    data = np.load(npz_path, allow_pickle=True)
    X = data["X"].astype(float)            # ref dosages
    rsids = list(data["rsids"])
    effect_alleles = list(data["effect_alleles"])
    superpops = list(data["superpop"])

    # Build user dosage vector aligned to rsids
    user_dose = np.full(len(rsids), np.nan)
    for i, (rsid, ea) in enumerate(zip(rsids, effect_alleles)):
        if rsid in snps_df.index:
            row = snps_df.loc[rsid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            d = _dosage(row.get("genotype"), ea)
            if d is not None:
                user_dose[i] = d
    callable_mask = ~np.isnan(user_dose)
    if callable_mask.sum() < 100:
        _log(f"Only {callable_mask.sum()} AIMs callable; full PCA needs >=100. Falling back.")
        return None

    # Restrict to callable SNPs for both ref and user
    X_sub = X[:, callable_mask]
    user_sub = user_dose[callable_mask]

    # Mean-center on ref to match PCA conventions
    mean = X_sub.mean(axis=0)
    X_c = X_sub - mean
    user_c = user_sub - mean

    pca = PCA(n_components=4)
    ref_pcs = pca.fit_transform(X_c)
    user_pcs = pca.transform(user_c.reshape(1, -1))[0]

    # Soft k-NN classifier in PC space against superpop centroids
    centroids = {}
    for sp in SUPERPOPS:
        idx = [i for i, s in enumerate(superpops) if s == sp]
        if idx:
            centroids[sp] = ref_pcs[idx].mean(axis=0)
    distances = {sp: float(np.linalg.norm(user_pcs - c)) for sp, c in centroids.items()}
    inv = {sp: 1.0 / (d + 1e-6) for sp, d in distances.items()}
    s = sum(inv.values())
    proportions = {sp: inv[sp] / s for sp in SUPERPOPS}

    # Plot
    plot_b64 = ""
    if plt is not None:
        fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=130)
        colors = {"EUR": "#3b82f6", "AFR": "#10b981", "EAS": "#f59e0b",
                  "SAS": "#ef4444", "AMR": "#a855f7"}
        for sp in SUPERPOPS:
            idx = [i for i, s in enumerate(superpops) if s == sp]
            if idx:
                ax.scatter(ref_pcs[idx, 0], ref_pcs[idx, 1],
                           s=8, alpha=0.55, c=colors[sp], label=f"{sp} ({SUPERPOP_LONG[sp]})")
        ax.scatter([user_pcs[0]], [user_pcs[1]],
                   s=160, c="black", marker="*",
                   edgecolors="white", linewidths=1.5, label="You", zorder=10)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title("Ancestry PCA — 1000 Genomes superpopulations")
        ax.legend(loc="best", fontsize=8, framealpha=0.85)
        ax.grid(alpha=0.2)
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        plot_b64 = base64.b64encode(buf.getvalue()).decode("ascii")

    sorted_props = sorted(proportions.items(), key=lambda x: -x[1])
    n_callable = int(callable_mask.sum())
    return {
        "available": True,
        "method": "Full PCA against 1000 Genomes Phase 3 superpopulations",
        "is_admixture_estimate": False,
        "n_aims_used": n_callable,
        "n_aims_independent": n_callable,
        "n_aims_expected": len(rsids),
        "proportions": proportions,
        "sorted_proportions": sorted_props,
        "primary_population": sorted_props[0][0],
        "runner_up_population": sorted_props[1][0] if len(sorted_props) > 1 else None,
        "ambiguous": (sorted_props[0][1] - sorted_props[1][1] < 0.15) if len(sorted_props) > 1 else False,
        "confidence": "high",
        "confidence_note": f"PCA + k-NN on {n_callable} reference AIMs (full 1000G panel).",
        "pc_coords": {"PC1": float(user_pcs[0]), "PC2": float(user_pcs[1]),
                       "PC3": float(user_pcs[2]) if len(user_pcs) > 2 else None,
                       "PC4": float(user_pcs[3]) if len(user_pcs) > 3 else None},
        "variance_explained": {
            "PC1": float(pca.explained_variance_ratio_[0]),
            "PC2": float(pca.explained_variance_ratio_[1]),
        },
        "plot_png_b64": plot_b64,
    }


def analyze_ancestry(snps_df: pd.DataFrame) -> Dict:
    """Public entry point. Try full PCA; fall back to heuristic."""
    full = _try_full_pca(snps_df)
    if full is not None:
        return full
    return estimate_ancestry_heuristic(snps_df)
