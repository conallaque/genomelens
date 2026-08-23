"""Unit tests for multi_person_module.py — KING-robust relationship inference."""

from __future__ import annotations

import random

import pandas as pd

import multi_person_module as mp


def _df(g: dict) -> pd.DataFrame:
    return pd.DataFrame({"genotype": g})


def _random_genome(n=3000, seed=0):
    random.seed(seed)
    alleles = ["AA", "AG", "GG", "CC", "CT", "TT"]
    return {f"rs{i}": random.choice(alleles) for i in range(n)}


def test_identical_samples_are_duplicate_twin() -> None:
    g = _random_genome(seed=1)
    r = mp.analyze_multi_person(_df(g), _df(dict(g)))
    assert abs(r["king"]["kinship"] - 0.5) < 0.02
    assert r["relationship"]["degree"] == "Duplicate / MZ twin"
    assert r["king"]["concordance"] == 1.0


def test_unrelated_random_genomes_classified_unrelated() -> None:
    a = _random_genome(seed=2)
    b = _random_genome(seed=99)
    r = mp.analyze_multi_person(_df(a), _df(b))
    assert r["relationship"]["degree"] == "Unrelated"


def test_parent_child_simulation_refined_to_parent_child() -> None:
    parent = _random_genome(seed=3)
    random.seed(7)
    child = {}
    for rsid, g in parent.items():
        inherited = random.choice(g)            # one allele from parent
        other = random.choice("ACGT")           # other parent random
        child[rsid] = "".join(sorted(inherited + other))
    r = mp.analyze_multi_person(_df(parent), _df(child))
    assert r["relationship"]["degree"] == "1st-degree"
    # IBS0 ≈ 0 → refined toward parent-child
    assert r["king"]["ibs0_rate"] < 0.0025
    assert "parent" in r["relationship"]["label"].lower()


def test_too_few_markers_is_indeterminate() -> None:
    small = _df({f"rs{i}": "AG" for i in range(50)})
    r = mp.analyze_multi_person(small, small)
    assert r["relationship"]["degree"] == "Indeterminate"


def test_kinship_formula_denominator_zero_safe() -> None:
    # All homozygous identical → no hets → denom 0 → kinship None, no crash
    g = {f"rs{i}": "AA" for i in range(500)}
    r = mp.king_kinship(_df(g), _df(dict(g)))
    assert r["kinship"] is None


def test_sex_chromosomes_excluded_when_chrom_present() -> None:
    df_a = pd.DataFrame({
        "chrom": ["1", "X", "Y", "MT"],
        "genotype": ["AG", "AG", "A", "GG"],
    }, index=["rs1", "rs2", "rs3", "rs4"])
    df_b = df_a.copy()
    r = mp.king_kinship(df_a, df_b)
    # only the autosomal rs1 (AG) is counted
    assert r["n_shared_snps"] == 1


def test_duplicate_rsids_deduplicated_not_inflated() -> None:
    """parse_dna_file output can contain duplicate rsIDs; king_kinship must
    dedup (first wins) so the shared-SNP count isn't inflated and .loc doesn't
    fan out."""
    g = _random_genome(n=1000, seed=11)
    df = _df(g)
    df_dup = pd.concat([df, df.iloc[:200]])  # 200 duplicate rows
    assert not df_dup.index.is_unique
    r = mp.king_kinship(df_dup, _df(dict(g)))
    assert r["n_shared_snps"] == 1000  # not 1200


def test_chip_scale_runtime_is_vectorised() -> None:
    """A whole-chip-scale comparison must complete quickly (vectorised), not
    loop per-marker. 200k markers should be well under a couple seconds."""
    import time
    g_a = _random_genome(n=200000, seed=21)
    g_b = _random_genome(n=200000, seed=22)
    t = time.time()
    r = mp.king_kinship(_df(g_a), _df(g_b))
    elapsed = time.time() - t
    assert r["n_shared_snps"] == 200000
    assert elapsed < 5.0, f"king_kinship too slow ({elapsed:.1f}s) — not vectorised"


def test_render_standalone_html_has_no_crash_and_key_fields() -> None:
    g = _random_genome(seed=5)
    r = mp.analyze_multi_person(_df(g), _df(dict(g)))
    html = mp.render_multi_person_html(r)
    assert "Genome Comparison" in html
    assert "KING kinship" in html
