"""Unit tests for the V35 analyze_family_planning() upgrade."""

from __future__ import annotations

import family_planning as fp


def _carrier(gene, variant, disease, inheritance, status="carrier"):
    return {"gene": gene, "variant": variant, "disease": disease,
            "inheritance": inheritance, "dosage": 1 if status == "carrier" else 2}


HFE = _carrier("HFE", "C282Y", "Hereditary Hemochromatosis (HH, type 1)",
               "autosomal recessive (incomplete penetrance)")
F5 = _carrier("F5", "Factor V Leiden (R506Q)",
              "Venous Thromboembolism (VTE) Susceptibility",
              "autosomal dominant (incomplete penetrance)")
FLG = _carrier("FLG", "R501X",
               "Atopic Dermatitis / Ichthyosis Vulgaris (filaggrin deficiency)",
               "semi-dominant")


def test_recessive_hardy_weinberg_math() -> None:
    r = fp.analyze_family_planning({"carriers": [HFE], "affected": []})
    it = next(x for x in r["recessive_items"] if x["gene"] == "HFE")
    # European HFE carrier freq 0.11 × 0.25 = 0.0275
    assert abs(it["child_two_copy_risk"] - 0.0275) < 1e-6
    # clinical = two-copy × penetrance (0.25, 0.60)
    lo, hi = it["child_clinical_risk"]
    assert abs(lo - 0.0275 * 0.25) < 1e-6
    assert abs(hi - 0.0275 * 0.60) < 1e-6


def test_dominant_transmission_and_penetrance_kept_separate() -> None:
    r = fp.analyze_family_planning({"carriers": [F5], "affected": []})
    it = next(x for x in r["dominant_items"] if x["gene"] == "F5")
    assert it["transmission"] == 0.50
    # penetrance is a separate descriptor, never merged into transmission
    assert "penetrance" in it["penetrance_text"].lower()


def test_mtdna_sex_gated_male_does_not_transmit() -> None:
    r = fp.analyze_family_planning({"carriers": [], "affected": []},
                                   inferred_sex="M")
    assert "not" in r["mtdna"]["transmission_note"].lower()
    assert r["mtdna"]["sex"] == "M"


def test_mtdna_sex_gated_female_transmits_all() -> None:
    r = fp.analyze_family_planning({"carriers": [], "affected": []},
                                   inferred_sex="F")
    assert "all" in r["mtdna"]["transmission_note"].lower()
    assert r["mtdna"]["sex"] == "F"


def test_semidominant_is_recessive_branch() -> None:
    r = fp.analyze_family_planning({"carriers": [FLG], "affected": []})
    it = next(x for x in r["recessive_items"] if x["gene"] == "FLG")
    assert it["semidominant"] is True
    assert it["child_two_copy_risk"] is not None


def test_hereditary_cancer_flags_mutyh_as_partner_relevant() -> None:
    r = fp.analyze_family_planning({"carriers": [], "affected": []})
    mutyh = next(c for c in r["hereditary_cancer"] if "MUTYH" in c["genes"])
    assert mutyh["partner_relevance"] == "HIGH"
    brca = next(c for c in r["hereditary_cancer"] if "BRCA1" in c["genes"])
    assert brca["partner_relevance"] == "low"


def test_affected_parent_uses_half_transmission() -> None:
    aff = _carrier("HFE", "C282Y", "Hereditary Hemochromatosis (HH, type 1)",
                   "autosomal recessive (incomplete penetrance)", status="affected")
    r = fp.analyze_family_planning({"carriers": [], "affected": [aff]})
    it = next(x for x in r["recessive_items"] if x["gene"] == "HFE")
    # affected parent: child two-copy = partner_freq × 0.5 = 0.11 × 0.5
    assert abs(it["child_two_copy_risk"] - 0.055) < 1e-6
