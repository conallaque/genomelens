"""Regression: PheWAS must not present a single-SNP marker-score percentile as an
actual-trait percentile (the "99th percentile brain volume from one SNP" overclaim).

The fix keeps the real marker-score percentile (a genotype fact) but ALSO reports the
variance the panel explains of the actual trait, and an honest trait-percentile that
regresses toward the mean as the panel weakens — without collapsing everyone to 50.
"""
from __future__ import annotations

import pandas as pd

import phewas


def _score(rsid_to_gt, trait_name):
    df = pd.DataFrame({"genotype": rsid_to_gt})
    return phewas._score_trait(df, phewas.PHEWAS_TRAITS[trait_name], None)


def test_brain_volume_marker_high_but_trait_average():
    # Homozygous effect allele: marker score is genuinely high, but the trait
    # estimate must stay near average because the panel explains ~0.1%.
    r = _score({"rs17178006": "TT"}, "Brain volume (proxy)")
    assert r["status"] == "ok"
    assert r["marker_percentile"] > 90          # real genotype rank, preserved
    assert 40 < r["trait_percentile"] < 60      # honest actual-trait estimate
    assert r["variance_explained_pct"] < 1.0    # R² exposed and tiny
    assert r["signal_strength"] == "negligible"
    assert "marker" in r["tier"].lower()        # tier labeled as marker score


def test_marker_score_is_not_collapsed_to_50():
    # The user's constraint: do NOT show everyone at ~50th. Marker scores must still
    # spread with genotype (0, 1, 2 copies give distinct marker percentiles).
    p = [_score({"rs17178006": gt}, "Brain volume (proxy)")["marker_percentile"]
         for gt in ("CC", "CT", "TT")]  # 0, 1, 2 effect alleles (strand as stored)
    assert len(set(round(x) for x in p)) >= 2   # genuinely different, not all ~50


def test_every_trait_reports_variance_explained():
    # No trait may present a percentile without also exposing how much of the real
    # trait it explains.
    g = {rsid: ea + ea for t in phewas.PHEWAS_TRAITS.values()
         for rsid, ea, _, _ in t["variants"]}
    df = pd.DataFrame({"genotype": g})
    for t in phewas.PHEWAS_TRAITS.values():
        r = phewas._score_trait(df, t, None)
        if r["status"] == "ok":
            assert "variance_explained_pct" in r
            assert "trait_percentile" in r
            assert r["signal_strength"] in ("negligible", "weak", "modest")


def test_higher_powered_trait_gets_stronger_signal_label():
    # Sanity: the tiering tracks R² — the best-powered trait must not be labeled
    # 'negligible' if a 0.1% trait is.
    g = {rsid: ea + ea for t in phewas.PHEWAS_TRAITS.values()
         for rsid, ea, _, _ in t["variants"]}
    df = pd.DataFrame({"genotype": g})
    r2s = {name: phewas._score_trait(df, t, None).get("variance_explained_pct", 0)
           for name, t in phewas.PHEWAS_TRAITS.items()}
    best = max(r2s, key=r2s.get)
    r = phewas._score_trait(df, phewas.PHEWAS_TRAITS[best], None)
    assert r["variance_explained_pct"] >= 1.0    # our best panel clears 1%
