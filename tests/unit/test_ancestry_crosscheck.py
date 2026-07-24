"""Ancestry: strand-aware dosage, palindrome safety, and the Y-DNA/mtDNA
geographic cross-check.

These lock in the fix for the real-world bug where a European man with Y-DNA
T1a1a was mis-called South Asian: the SLC45A2 rs16891982 C/G palindrome had a
flipped effect allele, and `_dosage` was not strand-aware (so LCT rs4988235,
reported A/G, always read as dosage 0)."""

from __future__ import annotations

import pandas as pd

import ancestry_pca as ap


def test_dosage_is_strand_aware() -> None:
    # LCT persistence allele A, reported on the A/G strand — must read 2, not 0.
    assert ap._dosage("AA", "A", "G") == 2
    # Same locus reported on the opposite strand (T/C) — must still read 2.
    assert ap._dosage("TT", "A", "G") == 2
    assert ap._dosage("TC", "A", "G") == 1
    assert ap._dosage("GG", "A", "G") == 0


def test_palindromic_slc45a2_excluded_from_scoring() -> None:
    # A homozygous European SLC45A2 genotype must not be scored (palindrome).
    df = pd.DataFrame({"genotype": {"rs16891982": "GG"}})
    res = ap.estimate_ancestry_heuristic(df)
    assert res["n_aims_palindromic"] == 1
    assert all(a["counted"] is False for a in res["used_aims"]
               if a["rsid"] == "rs16891982")


def _european_profile() -> pd.DataFrame:
    """Genotypes matching the real T1a1a sample that used to mis-call as SAS."""
    return pd.DataFrame({"genotype": {
        "rs1426654": "AA",    # SLC24A5 light — near-fixed European
        "rs16891982": "GG",   # SLC45A2 European (palindrome — excluded)
        "rs12913832": "AG",   # HERC2
        "rs1042602": "AC",    # TYR
        "rs1805007": "CC",    # MC1R
        "rs3827760": "AA",    # EDAR — not East Asian
        "rs17822931": "CC",   # ABCC11 — not East Asian
        "rs671": "GG",        # ALDH2 — not East Asian
        "rs2814778": "TT",    # Duffy — not African
        "rs4988235": "AA",    # LCT persistence (strand-aware!)
    }})


def test_european_sample_calls_european_not_asian() -> None:
    res = ap.estimate_ancestry_heuristic(_european_profile())
    assert res["primary_population"] == "EUR"
    assert res["proportions"]["EUR"] > res["proportions"]["SAS"]
    assert res["proportions"]["EUR"] > 0.8


def test_crosscheck_concordant_for_t1a1a_and_european_autosomal() -> None:
    auto = ap.estimate_ancestry_heuristic(_european_profile())
    y = {"terminal_haplogroup": "T1a1a", "confidence": "high",
         "path": [{"haplogroup": h} for h in
                  ["CT", "F", "K", "LT", "T", "T1a", "T1a1", "T1a1a"]]}
    mt = {"haplogroup": "T", "confidence": "high"}
    result = ap.analyze_ancestry(_european_profile(), y_result=y, mt_result=mt)
    cc = result["haplogroup_crosscheck"]
    assert cc["verdict"] == "concordant"
    assert cc["paternal"]["haplogroup"] == "T1a1a"
    assert cc["paternal"]["dominant"] == "EUR"


def test_crosscheck_flags_discordance_and_caps_confidence() -> None:
    # A (hypothetical) East-Asian autosomal call against T1a1a must be flagged
    # discordant and forced to low confidence.
    auto = {"available": True, "primary_population": "EAS", "confidence": "moderate"}
    y = {"terminal_haplogroup": "T1a1a", "confidence": "high", "path": []}
    cc = ap.cross_check_ancestry(auto, y, None)
    assert cc["verdict"] == "discordant"
    assert cc["suggested_confidence"] == "low"
