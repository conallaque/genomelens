

# ── low-coverage suppression ────────────────────────────────────────────────
def test_low_coverage_is_suppressed_not_rendered_with_a_warning():
    """Coverage below the declared threshold yields no percentile at all.

    The gate used to be `n_used == 0`, so one callable variant out of 1.7M
    produced a full percentile and tier carrying a banner that called the same
    percentile unreliable. A real genome scored 28% coverage for coronary
    artery disease and rendered "High, 100th percentile".
    """
    import pandas as pd

    from risk.pgs_catalog import PGS_MODERATE_COVERAGE_PCT, calculate_pgs

    variants = [{"rsid": f"rs{i}", "chrom": "1", "pos": str(i),
                 "effect_allele": "A", "other_allele": "G",
                 "weight": 0.1, "af": 0.3} for i in range(100)]

    def _typed(n):
        """n variants whose mean dose matches 2*af=0.6, so only coverage varies.

        An all-heterozygous set carries dose 1.0 against a predicted 0.6 and
        trips the informative-missingness gate instead — which is correct
        behavior, and exactly why this fixture mixes hom-reference calls in.
        """
        gts = {f"rs{i}": ("AG" if i % 5 < 3 else "GG") for i in range(n)}
        return pd.DataFrame({"genotype": list(gts.values())}, index=list(gts))

    r = calculate_pgs(_typed(20), variants)          # 20% -> under threshold
    assert r["status"] == "insufficient_data"
    assert "percentile" not in r
    assert str(PGS_MODERATE_COVERAGE_PCT) in r["reason"]

    r2 = calculate_pgs(_typed(60), variants)         # 60% -> above threshold
    assert r2["status"] == "computed"
    assert r2["percentile"] is not None


def test_a_callset_with_no_hom_reference_calls_is_not_given_a_percentile():
    """Absence means hom-ref on a gVCF, so the callable set is not random.

    A real 30x genome had zero 0/0 calls among 502,690 callable scoring
    variants; every retained variant carried an alt allele, and five unrelated
    traits pinned at exactly the 100th percentile. The bias is not correctable
    from the callable set, so the percentile is withheld.
    """
    import pandas as pd

    from risk.pgs_catalog import calculate_pgs

    variants = [{"rsid": f"rs{i}", "chrom": "1", "pos": str(i),
                 "effect_allele": "A", "other_allele": "G",
                 "weight": 0.1, "af": 0.3} for i in range(100)]
    # Every callable variant carries the effect allele: dose 1.0 vs 0.6 predicted.
    gts = {f"rs{i}": "AG" for i in range(80)}
    df = pd.DataFrame({"genotype": list(gts.values())}, index=list(gts))
    r = calculate_pgs(df, variants)
    assert r["status"] == "insufficient_data"
    assert r["informative_missingness"] is True
    assert r["dose_ratio"] > 1.5
    assert "percentile" not in r


def test_prs_and_pgs_disagreement_is_flagged_not_hidden():
    """Two modules scoring one condition must not disagree silently."""
    from risk.pgs_catalog import reconcile_prs_pgs

    prs = {"panels": {"Alzheimer's Disease (Late-onset)": {"result": {
        "status": "computed", "percentile": 96.0, "tier": "High",
        "confidence": "moderate"}}}}
    pgs = {"panels": {"alzheimers_disease": {"pgs_id": "PGS000334", "result": {
        "status": "computed", "percentile": 13.0, "tier": "Below Average",
        "confidence": "moderate"}}}}
    rows = reconcile_prs_pgs(prs, pgs)
    assert len(rows) == 1
    assert rows[0]["divergent"] is True
    assert rows[0]["percentile_gap"] == 83.0
    assert rows[0]["tiers_agree"] is False


def test_a_suppressed_score_is_not_reported_as_a_disagreement():
    """A withheld score is missing information, not conflicting information."""
    from risk.pgs_catalog import reconcile_prs_pgs

    prs = {"panels": {"Type 2 Diabetes": {"result": {
        "status": "computed", "percentile": 90.0, "tier": "High"}}}}
    pgs = {"panels": {"type_2_diabetes": {"result": {
        "status": "insufficient_data", "reason": "coverage"}}}}
    assert reconcile_prs_pgs(prs, pgs) == []
