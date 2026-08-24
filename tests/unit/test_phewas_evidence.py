"""PheWAS producer now returns the per-trait evidence it computes.

Locks the additive contract change: _score_trait includes used_variants
(rsid/effect_allele/dose/beta/af) and missing_variants so the renderer can
surface an auditable breakdown instead of only a coverage count.
"""

import pandas as pd

from risk import phewas


def _df_for_trait(trait_name):
    trait = phewas.PHEWAS_TRAITS[trait_name]
    idx = [v[0] for v in trait["variants"]]
    # homozygous effect allele for every typed variant
    gts = [v[1] + v[1] for v in trait["variants"]]
    return pd.DataFrame({"genotype": gts}, index=idx)


def test_score_trait_returns_used_and_missing_variants():
    name = next(iter(phewas.PHEWAS_TRAITS))
    res = phewas.analyze_phewas(_df_for_trait(name), sex="male")
    t = res["traits"][name]["result"]
    assert t["status"] == "ok"
    assert t.get("used_variants"), "used_variants missing"
    uv = t["used_variants"][0]
    for k in ("rsid", "effect_allele", "dose", "beta", "af"):
        assert k in uv, f"used_variant missing key {k}"
    assert "missing_variants" in t  # present (possibly empty) for typed traits
