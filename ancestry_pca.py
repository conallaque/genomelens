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

# Each AIM now declares BOTH + strand alleles (`effect_allele` + `other_allele`)
# so the dosage read can be strand-aware (see `_dosage`). Consumer chips report
# on either strand; the previous naive `count(effect_allele)` silently returned
# dosage 0 whenever the chip's strand differed from the table's, which flipped
# strong European markers (e.g. LCT rs4988235 read as 0 on an A/G-strand chip).
#
# `palindromic: True` marks A/T or C/G SNPs whose two alleles are reverse
# complements. Their strand cannot be recovered from genotype alone, so they are
# DISPLAYED for transparency but EXCLUDED from the likelihood. This is what
# previously mis-called European samples: SLC45A2 rs16891982 is a C/G SNP whose
# European (light-skin) allele is G at ~98%, but the table hard-coded "C", so a
# European GG genotype scored as "zero European alleles" and collapsed the EUR
# posterior. Excluding palindromes removes that entire failure mode.
#
# Frequencies are the effect-allele frequency in 1000 Genomes Phase 3 (EUR, AFR,
# EAS, SAS, AMR), rounded. `selection` flags loci under strong recent selection
# (diet/pigmentation) whose frequency tracks environment as well as ancestry —
# kept in the panel but the report labels them so they aren't over-read.

AIMS_PRIORS: Dict[str, Dict] = {
    # ── Skin pigmentation (strong continental discriminators) ──────────────
    "rs1426654":  {"effect_allele": "A", "other_allele": "G", "EUR": 0.999, "AFR": 0.065, "EAS": 0.006, "SAS": 0.55, "AMR": 0.48, "gene": "SLC24A5", "note": "Ala111Thr — derived A = lighter skin, near-fixed in Europe."},
    "rs16891982": {"effect_allele": "G", "other_allele": "C", "EUR": 0.98,  "AFR": 0.05,  "EAS": 0.01,  "SAS": 0.30, "AMR": 0.55, "gene": "SLC45A2", "palindromic": True, "note": "C/G palindrome — displayed only, cannot be strand-oriented from genotype (this was the marker that mis-called European samples)."},
    "rs1042602":  {"effect_allele": "A", "other_allele": "C", "EUR": 0.36,  "AFR": 0.003, "EAS": 0.001, "SAS": 0.05, "AMR": 0.12, "gene": "TYR", "note": "S192Y — A allele enriched in Europeans."},
    "rs1805007":  {"effect_allele": "T", "other_allele": "C", "EUR": 0.09,  "AFR": 0.005, "EAS": 0.001, "SAS": 0.005, "AMR": 0.02, "gene": "MC1R", "note": "R151C red-hair allele — European-specific, low frequency."},
    "rs12913832": {"effect_allele": "G", "other_allele": "A", "EUR": 0.635, "AFR": 0.03,  "EAS": 0.005, "SAS": 0.05, "AMR": 0.22, "gene": "HERC2", "ld_block": "HERC2_eye", "note": "Blue-eye allele — European-enriched."},
    # ── East Asian / Native American markers ───────────────────────────────
    "rs3827760":  {"effect_allele": "G", "other_allele": "A", "EUR": 0.005, "AFR": 0.003, "EAS": 0.90,  "SAS": 0.02, "AMR": 0.35, "gene": "EDAR", "ld_block": "EDAR", "note": "V370A — hair/tooth morphology, near-diagnostic for East Asian / Native American ancestry."},
    "rs17822931": {"effect_allele": "T", "other_allele": "C", "EUR": 0.05,  "AFR": 0.005, "EAS": 0.85,  "SAS": 0.05, "AMR": 0.30, "gene": "ABCC11", "note": "Dry-earwax allele — East Asian enriched."},
    "rs671":      {"effect_allele": "A", "other_allele": "G", "EUR": 0.002, "AFR": 0.002, "EAS": 0.21,  "SAS": 0.005, "AMR": 0.02, "gene": "ALDH2", "selection": True, "note": "ALDH2*2 alcohol-flush allele — essentially East-Asian-restricted."},
    "rs1229984":  {"effect_allele": "A", "other_allele": "G", "EUR": 0.03,  "AFR": 0.005, "EAS": 0.70,  "SAS": 0.07, "AMR": 0.18, "gene": "ADH1B", "selection": True, "note": "Arg48His — high in East Asia."},
    # ── African marker ─────────────────────────────────────────────────────
    "rs2814778":  {"effect_allele": "C", "other_allele": "T", "EUR": 0.005, "AFR": 0.97,  "EAS": 0.005, "SAS": 0.01, "AMR": 0.06, "gene": "DARC/ACKR1", "note": "Duffy-null — near-fixed in sub-Saharan Africa, near-absent elsewhere."},
    # ── Diet-selected (lactase persistence) ────────────────────────────────
    "rs4988235":  {"effect_allele": "A", "other_allele": "G", "EUR": 0.51,  "AFR": 0.10,  "EAS": 0.005, "SAS": 0.30, "AMR": 0.30, "gene": "LCT", "ld_block": "LCT_MCM6", "selection": True, "note": "-13910C>T lactase persistence — European/pastoralist enriched."},
    "rs182549":   {"effect_allele": "T", "other_allele": "C", "EUR": 0.51,  "AFR": 0.10,  "EAS": 0.005, "SAS": 0.30, "AMR": 0.30, "gene": "MCM6", "ld_block": "LCT_MCM6", "selection": True, "note": "MCM6 intronic — mirror tag of LCT persistence."},
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


_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _dosage(genotype: object, effect_allele: str,
            other_allele: Optional[str] = None) -> Optional[int]:
    """Strand-aware effect-allele dosage (0/1/2) from a 2-char genotype.

    When `other_allele` is supplied, detect whether the genotype is reported on
    the + or − strand (relative to the table's alleles) and count on the matching
    strand. This fixes the old `count(effect_allele)` behaviour, which returned 0
    whenever the chip's strand differed from the table's — silently zeroing out
    strong markers such as LCT rs4988235 on an A/G-strand consumer chip.

    Palindromic markers (effect/other are reverse complements) cannot be
    oriented and should be excluded upstream; if one slips through we still fall
    back to a plain + strand count rather than guessing.
    """
    if genotype is None:
        return None
    gt = str(genotype).upper().replace(" ", "").replace("-", "")
    if gt in ("", "NAN", "--") or len(gt) != 2:
        return None
    eff = effect_allele.upper()
    if other_allele is None:
        return gt.count(eff)

    oth = other_allele.upper()
    eff_c, oth_c = _COMPLEMENT.get(eff, eff), _COMPLEMENT.get(oth, oth)
    gt_alleles = set(gt)
    if gt_alleles <= {eff, oth}:            # reported on the table's strand
        return gt.count(eff)
    if gt_alleles <= {eff_c, oth_c}:        # reported on the opposite strand
        return gt.count(eff_c)
    return gt.count(eff)                    # mixed/noise — best-effort + strand


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
    n_palindromic = 0
    for rsid, info in AIMS_PRIORS.items():
        if rsid not in snps_df.index:
            continue
        row = snps_df.loc[rsid]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        dose = _dosage(row.get("genotype"), info["effect_allele"],
                       info.get("other_allele"))
        if dose is None:
            continue
        # Palindromic (A/T, C/G) markers can't be strand-oriented from genotype
        # alone — show them but never let them drive the likelihood.
        if info.get("palindromic"):
            n_palindromic += 1
            used_aims.append({
                "rsid": rsid, "gene": info["gene"],
                "genotype": str(row.get("genotype")).upper(),
                "effect_allele": info["effect_allele"], "dosage": dose,
                "ld_block": info.get("ld_block", rsid), "counted": False,
                "palindromic": True, "note": info.get("note", ""),
                "selection": info.get("selection", False),
            })
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
                "palindromic": False, "note": info.get("note", ""),
                "selection": info.get("selection", False),
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
            "palindromic": False, "note": info.get("note", ""),
            "selection": info.get("selection", False),
        })

    n_independent = len(seen_blocks)
    if n_independent == 0:
        return {
            "available": False,
            "reason": ("No scoreable AIM SNPs detected on chip — cannot estimate "
                       "ancestry (palindromic markers, if any, are excluded)."),
            "n_aims_used": len(used_aims),
            "n_aims_independent": 0,
            "n_aims_palindromic": n_palindromic,
            "n_aims_expected": N_INDEPENDENT_AIMS,
            "used_aims": used_aims,
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

    # DEMOTED: this heuristic must NOT emit a headline "you are population X".
    # Same PPV-collapse logic as the PheWAS demotion (Wilson & Jungner 1968): too few
    # informative markers cannot support a confident classification. docs/METHODS.md §22.
    # ancestry call. The panel is ~10 usable markers, several under strong natural
    # selection (pigmentation and diet genes: EDAR, ALDH2, ADH1B, LCT), whose
    # frequencies track local adaptation, not shared descent. A single selected
    # allele can flip the call — which is exactly how it mis-classified a European
    # sample as East Asian. A defensible superpopulation call requires the full
    # 1000G PCA path (reference/ancestry/ present); until then we present these as
    # what they actually are — individual pigmentation/trait markers — and set
    # primary_population to None so no downstream layer asserts an ancestry identity.
    return {
        "available": True,
        "method": ("Pigmentation & trait AIMs (individual markers shown for "
                   "interest only — NOT an ancestry classification; the panel is "
                   "too small and too selection-confounded to call ancestry)"),
        "is_admixture_estimate": False,
        "is_ancestry_call": False,
        "ancestry_call_suppressed": True,
        "suppression_reason": (
            "A superpopulation ancestry call needs the full 1000 Genomes PCA "
            "reference (reference/ancestry/). The fallback marker panel is too "
            "small and several markers are under natural selection, so it cannot "
            "reliably distinguish ancestry — it is shown as individual trait "
            "markers instead."),
        "headline": "Pigmentation & trait markers (not an ancestry call)",
        "n_aims_used": len(used_aims),
        "n_aims_independent": n_independent,
        "n_aims_redundant": n_redundant,
        "n_aims_palindromic": n_palindromic,
        "n_aims_expected": N_INDEPENDENT_AIMS,
        # Kept for the transparent per-marker table, clearly labelled as affinity
        # only — NOT a population identity.
        "marker_best_affinity": primary,
        "marker_affinity_proportions": proportions,
        "sorted_proportions": sorted_props,
        "primary_population": None,          # suppressed by design (see above)
        "runner_up_population": runner_up,
        "evidence_margin_nats": round(margin_nats, 2),
        "ambiguous": ambiguous,
        "confidence": "not applicable — no ancestry call made",
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


# ── Uniparental (Y-DNA / mtDNA) geographic cross-check ────────────────────────
#
# A whole-genome (autosomal) ancestry estimate and a haplogroup are different
# things: the autosomal panel samples DNA from every ancestor, while a Y-DNA or
# mtDNA haplogroup traces ONE line (strict paternal or strict maternal) back tens
# of thousands of years. They are not interchangeable — but they must be
# geographically *compatible*. When a small autosomal panel calls, say, "South
# Asian" for a man whose Y-DNA is a Near-Eastern/European/East-African lineage
# with essentially zero South-Asian majority, that contradiction is a red flag
# that the autosomal call is being driven by noise or a marker artefact, not by
# real ancestry. This module surfaces that check explicitly.
#
# Each haplogroup maps to (a) a human-readable geographic homeland and (b) a soft
# distribution over the five 1000G superpopulations describing where the lineage
# reaches appreciable frequency. Keys are matched by longest-prefix against the
# terminal haplogroup, so "T1a1a" resolves via the "T" entry, "R1b1a2" via "R1b".

_Y_GEOGRAPHY: Dict[str, Dict] = {
    # key: (region, {superpop: weight}, note)
    "R1b":  {"region": "Western & Atlantic Europe", "dist": {"EUR": 0.90, "AMR": 0.05, "SAS": 0.03, "AFR": 0.01, "EAS": 0.01}, "note": "The single most common paternal lineage in Western Europe."},
    "R1a":  {"region": "Eastern Europe, the Steppe, Central & South Asia", "dist": {"EUR": 0.55, "SAS": 0.35, "AMR": 0.05, "AFR": 0.025, "EAS": 0.025}, "note": "Spread with Bronze-Age steppe expansions; spans Slavic Europe and northern India."},
    "R1":   {"region": "Europe & South/Central Asia", "dist": {"EUR": 0.70, "SAS": 0.22, "AMR": 0.04, "AFR": 0.02, "EAS": 0.02}, "note": "Ancestor of R1a and R1b."},
    "I":    {"region": "Europe (indigenous Palaeolithic lineage)", "dist": {"EUR": 0.93, "AMR": 0.03, "SAS": 0.02, "AFR": 0.01, "EAS": 0.01}, "note": "The oldest Europe-specific Y lineage; peaks in Scandinavia (I1) and the Balkans/Sardinia (I2)."},
    "J":    {"region": "The Near East, Arabia, the Caucasus & the Mediterranean", "dist": {"EUR": 0.55, "SAS": 0.20, "AFR": 0.15, "AMR": 0.05, "EAS": 0.05}, "note": "Spread with Neolithic farming out of the Fertile Crescent."},
    "G":    {"region": "The Caucasus, Anatolia & Neolithic Europe", "dist": {"EUR": 0.72, "SAS": 0.16, "AFR": 0.06, "AMR": 0.03, "EAS": 0.03}, "note": "Carried by Europe's early farmers (e.g. Ötzi the Iceman)."},
    "T":    {"region": "The Near East, the Mediterranean & the Horn of Africa", "dist": {"EUR": 0.55, "AFR": 0.22, "SAS": 0.15, "AMR": 0.04, "EAS": 0.04}, "note": "An old, geographically scattered lineage of the Near East, Mediterranean Europe and East Africa. It is NOT an East- or South-East-Asian lineage — it reaches South Asia only at low frequency and is essentially absent as a majority there."},
    "L":    {"region": "South Asia & the Near East", "dist": {"SAS": 0.75, "EUR": 0.15, "AFR": 0.05, "AMR": 0.03, "EAS": 0.02}, "note": "Concentrated in the Indian subcontinent."},
    "E":    {"region": "Africa, the Levant & Mediterranean Europe", "dist": {"AFR": 0.62, "EUR": 0.28, "SAS": 0.05, "AMR": 0.03, "EAS": 0.02}, "note": "The most common African paternal lineage, with a Mediterranean reach."},
    "N":    {"region": "Northern Eurasia, Siberia & the Baltic", "dist": {"EUR": 0.50, "EAS": 0.42, "SAS": 0.04, "AMR": 0.02, "AFR": 0.02}, "note": "Common among Uralic-speaking and Siberian peoples."},
    "O":    {"region": "East & South-East Asia", "dist": {"EAS": 0.94, "SAS": 0.03, "AMR": 0.01, "EUR": 0.01, "AFR": 0.01}, "note": "The dominant paternal lineage of East Asia."},
    "Q":    {"region": "Siberia & the indigenous Americas", "dist": {"AMR": 0.70, "EAS": 0.24, "EUR": 0.03, "SAS": 0.02, "AFR": 0.01}, "note": "Crossed Beringia; the principal Native-American paternal lineage."},
    "C":    {"region": "Asia, Oceania & the Americas", "dist": {"EAS": 0.55, "AMR": 0.20, "SAS": 0.18, "EUR": 0.04, "AFR": 0.03}, "note": "A deep, widely scattered non-African lineage."},
}

_MT_GEOGRAPHY: Dict[str, Dict] = {
    "HV": {"region": "Europe & the Near East", "dist": {"EUR": 0.80, "SAS": 0.10, "AFR": 0.05, "AMR": 0.03, "EAS": 0.02}},
    "H":  {"region": "Europe (the most common European maternal lineage)", "dist": {"EUR": 0.85, "SAS": 0.07, "AFR": 0.04, "AMR": 0.02, "EAS": 0.02}},
    "V":  {"region": "Western Europe & the Mediterranean", "dist": {"EUR": 0.88, "AFR": 0.05, "SAS": 0.03, "AMR": 0.02, "EAS": 0.02}},
    "U":  {"region": "Europe, the Near East & South Asia", "dist": {"EUR": 0.65, "SAS": 0.22, "AFR": 0.06, "AMR": 0.04, "EAS": 0.03}},
    "K":  {"region": "Europe & the Near East", "dist": {"EUR": 0.80, "SAS": 0.10, "AFR": 0.05, "AMR": 0.03, "EAS": 0.02}},
    "J":  {"region": "Europe & the Near East", "dist": {"EUR": 0.72, "SAS": 0.15, "AFR": 0.07, "AMR": 0.03, "EAS": 0.03}},
    "T":  {"region": "Europe & the Near East", "dist": {"EUR": 0.72, "SAS": 0.15, "AFR": 0.07, "AMR": 0.03, "EAS": 0.03}},
    "W":  {"region": "Europe, the Near East & South Asia", "dist": {"EUR": 0.60, "SAS": 0.30, "AFR": 0.05, "AMR": 0.03, "EAS": 0.02}},
    "X":  {"region": "Europe, the Near East & (rarely) North America", "dist": {"EUR": 0.70, "SAS": 0.12, "AFR": 0.08, "AMR": 0.07, "EAS": 0.03}},
    "I":  {"region": "Europe & the Near East", "dist": {"EUR": 0.75, "SAS": 0.12, "AFR": 0.07, "AMR": 0.03, "EAS": 0.03}},
    "N":  {"region": "Eurasia (macro-lineage)", "dist": {"EUR": 0.55, "SAS": 0.20, "EAS": 0.15, "AFR": 0.05, "AMR": 0.05}},
    "L":  {"region": "Sub-Saharan Africa", "dist": {"AFR": 0.92, "EUR": 0.03, "SAS": 0.02, "AMR": 0.02, "EAS": 0.01}},
    "M":  {"region": "Asia (macro-lineage: East/South Asia & the Americas)", "dist": {"EAS": 0.45, "SAS": 0.35, "AMR": 0.12, "EUR": 0.05, "AFR": 0.03}},
    "D":  {"region": "East Asia & the Americas", "dist": {"EAS": 0.55, "AMR": 0.35, "SAS": 0.05, "EUR": 0.03, "AFR": 0.02}},
    "A":  {"region": "East Asia & the Americas", "dist": {"AMR": 0.50, "EAS": 0.42, "SAS": 0.04, "EUR": 0.02, "AFR": 0.02}},
    "B":  {"region": "East/South-East Asia & the Americas", "dist": {"EAS": 0.55, "AMR": 0.35, "SAS": 0.05, "EUR": 0.03, "AFR": 0.02}},
    "C":  {"region": "Siberia, East Asia & the Americas", "dist": {"AMR": 0.50, "EAS": 0.42, "SAS": 0.04, "EUR": 0.02, "AFR": 0.02}},
    "F":  {"region": "East & South-East Asia", "dist": {"EAS": 0.85, "SAS": 0.08, "AMR": 0.03, "EUR": 0.02, "AFR": 0.02}},
}


def _match_geography(terminal: str, path_labels: List[str], table: Dict[str, Dict]) -> Optional[Dict]:
    """Longest-prefix match of a haplogroup label against a geography table.

    Tries the terminal haplogroup first, then each node from deepest to
    shallowest, then a single-letter fallback. Longer keys win (so "R1b" beats
    "R1" beats a bare fallback)."""
    candidates = [terminal] + [p for p in reversed(path_labels or []) if p]
    keys_by_len = sorted(table.keys(), key=len, reverse=True)
    for cand in candidates:
        c = str(cand).strip().upper()
        for key in keys_by_len:
            if c == key.upper() or c.startswith(key.upper()):
                out = dict(table[key])
                out["matched_key"] = key
                return out
    # Fallback: bare first letter (e.g. an exotic subclade of a known macro-hg)
    if terminal:
        first = str(terminal).strip().upper()[:1]
        if first in table:
            out = dict(table[first])
            out["matched_key"] = first
            return out
    return None


def _argmax_dist(dist: Dict[str, float]) -> str:
    return max(dist.items(), key=lambda kv: kv[1])[0] if dist else ""


def haplogroup_geographic_prior(y_result: Optional[Dict],
                                mt_result: Optional[Dict]) -> Dict:
    """Turn Y-DNA and mtDNA calls into geographic expectations over the five
    1000G superpopulations, for use as an ancestry sanity-check."""
    out: Dict = {"paternal": None, "maternal": None}

    if y_result and y_result.get("terminal_haplogroup") not in (
            None, "Unknown", "Insufficient Y-chromosome SNPs"):
        term = y_result.get("terminal_haplogroup", "")
        labels = [n.get("haplogroup", "") for n in y_result.get("path", [])]
        geo = _match_geography(term, labels, _Y_GEOGRAPHY)
        if geo:
            out["paternal"] = {
                "haplogroup": term,
                "region": geo["region"],
                "dist": geo["dist"],
                "dominant": _argmax_dist(geo["dist"]),
                "note": geo.get("note", ""),
                "confidence": y_result.get("confidence", "unknown"),
            }

    if mt_result and mt_result.get("haplogroup") not in (None, "Unknown"):
        term = mt_result.get("haplogroup", "")
        geo = _match_geography(term, [], _MT_GEOGRAPHY)
        if geo:
            out["maternal"] = {
                "haplogroup": term,
                "region": geo["region"],
                "dist": geo["dist"],
                "dominant": _argmax_dist(geo["dist"]),
                "note": geo.get("note", ""),
                "confidence": mt_result.get("confidence", "unknown"),
            }
    return out


def cross_check_ancestry(autosomal: Dict,
                         y_result: Optional[Dict],
                         mt_result: Optional[Dict]) -> Optional[Dict]:
    """Compare the autosomal top call against the paternal/maternal lineages.

    Returns a verdict dict (or None if no uniparental data). The verdict is one
    of concordant / plausible / discordant / uninformative, plus a plain-English
    explanation. A discordant verdict is a strong signal that a small autosomal
    panel has misfired and should not be trusted at face value.
    """
    if not autosomal or not autosomal.get("available"):
        return None
    prior = haplogroup_geographic_prior(y_result, mt_result)
    pat, mat = prior.get("paternal"), prior.get("maternal")
    if not pat and not mat:
        return None

    # When the autosomal ancestry call is suppressed (fallback marker panel), the
    # cross-check falls back to the weak marker affinity, clearly flagged — it is a
    # sanity check against the haplogroups, not an ancestry determination.
    primary = autosomal.get("primary_population") or autosomal.get("marker_best_affinity")
    primary_long = SUPERPOP_LONG.get(primary, primary)

    lines: List[str] = []
    verdicts: List[str] = []

    def _assess(line: Dict, which: str) -> None:
        dist = line["dist"]
        weight_here = dist.get(primary, 0.0)
        dominant = line["dominant"]
        if weight_here >= 0.30 or primary == dominant:
            verdicts.append("concordant")
            lines.append(
                f"Your {which} lineage <strong>{line['haplogroup']}</strong> is "
                f"typical of {line['region']} — consistent with an autosomal "
                f"top match of {primary_long}."
            )
        elif weight_here >= 0.12:
            verdicts.append("plausible")
            lines.append(
                f"Your {which} lineage <strong>{line['haplogroup']}</strong> "
                f"({line['region']}) reaches only modest frequency in "
                f"{primary_long}; the two are compatible but not a clean match."
            )
        else:
            verdicts.append("discordant")
            lines.append(
                f"Your {which} lineage <strong>{line['haplogroup']}</strong> is a "
                f"lineage of {line['region']}, where {primary_long} ancestry is "
                f"minimal. {line.get('note','')} A {primary_long} autosomal call "
                f"is hard to reconcile with this lineage — on a panel this small, "
                f"that points to an unreliable autosomal estimate rather than real "
                f"{primary_long} ancestry."
            )

    if pat:
        _assess(pat, "paternal (Y-DNA)")
    if mat:
        _assess(mat, "maternal (mtDNA)")

    if "discordant" in verdicts:
        verdict = "discordant"
    elif "concordant" in verdicts and "plausible" not in verdicts:
        verdict = "concordant"
    elif "concordant" in verdicts or "plausible" in verdicts:
        verdict = "plausible"
    else:
        verdict = "uninformative"

    # A discordant lineage check, combined with a small/ambiguous panel, caps how
    # much the autosomal call should be trusted.
    suggested_conf = autosomal.get("confidence", "low")
    if verdict == "discordant":
        suggested_conf = "low"

    return {
        "verdict": verdict,
        "primary_population": primary,
        "primary_population_long": primary_long,
        "paternal": pat,
        "maternal": mat,
        "explanations": lines,
        "suggested_confidence": suggested_conf,
        "summary": {
            "concordant": "Your autosomal estimate agrees with your deep paternal/maternal lineages.",
            "plausible": "Your autosomal estimate is broadly compatible with your deep lineages.",
            "discordant": "Your autosomal estimate CONFLICTS with your deep paternal/maternal lineages — treat the small-panel autosomal call as unreliable.",
            "uninformative": "Your haplogroups don't add a strong geographic constraint here.",
        }[verdict],
    }


def analyze_ancestry(snps_df: pd.DataFrame,
                     y_result: Optional[Dict] = None,
                     mt_result: Optional[Dict] = None) -> Dict:
    """Public entry point. Try full PCA; fall back to heuristic. When Y-DNA and/or
    mtDNA results are supplied, attach a uniparental geographic cross-check that
    flags autosomal calls incompatible with the deep paternal/maternal lineages.
    """
    full = _try_full_pca(snps_df)
    result = full if full is not None else estimate_ancestry_heuristic(snps_df)

    if result.get("available") and (y_result or mt_result):
        cc = cross_check_ancestry(result, y_result, mt_result)
        if cc is not None:
            result["haplogroup_crosscheck"] = cc
            # A discordant deep-lineage check overrides an over-confident
            # small-panel call: never report better than "low" when the paternal
            # or maternal line flatly contradicts the autosomal top match.
            if cc["verdict"] == "discordant":
                result["confidence"] = "low"
                result["ambiguous"] = True
                extra = (
                    "  Uniparental cross-check: this autosomal call CONFLICTS "
                    "with your deep paternal/maternal lineage(s) and should be "
                    "treated as unreliable — see the Lineage Cross-Check below."
                )
                result["confidence_note"] = (result.get("confidence_note", "") or "") + extra
    return result
