"""Regression: recessive carrier findings must produce a numeric offspring risk.

``family_planning._RECESSIVE_CARRIER_FREQ`` previously held two genes (HFE and
FLG). Every other recessive condition ``carrier.py`` can detect — cystic
fibrosis, sickle cell, SMA, Tay-Sachs and the rest, i.e. the conditions carrier
screening exists for — fell through with ``child_two_copy_risk = None``, so the
one quantitative reproductive figure in the report covered hemochromatosis and
eczema only.
"""
from __future__ import annotations

import carrier
import family_planning as fp


def _recessive_genes_detectable():
    """Genes carrier.py can call that are recessive (or semi-dominant), i.e.
    the ones for which partner carrier frequency is the relevant quantity."""
    out = set()
    for v in carrier.CARRIER_VARIANTS:
        inh = (v.get("inheritance") or "").lower()
        if ("recessive" in inh and "x-linked" not in inh) or "semi-dominant" in inh or "codominant" in inh:
            out.add(v["gene"])
    return out


def test_every_detectable_recessive_gene_has_a_frequency():
    missing = sorted(_recessive_genes_detectable() - set(fp._RECESSIVE_CARRIER_FREQ))
    assert not missing, (
        f"carrier.py can report these recessive genes but family_planning has no "
        f"carrier frequency for them, so offspring risk silently comes back None: "
        f"{missing}")


def test_frequencies_are_plausible_probabilities():
    for gene, spec in fp._RECESSIVE_CARRIER_FREQ.items():
        for anc, f in spec["freq"].items():
            assert 0.0 < f < 0.5, f"{gene}/{anc} carrier frequency {f} is implausible"


def test_every_entry_has_european_fallback():
    # The consumer falls back to the European key when ancestry is unknown, so a
    # missing one silently yields None.
    for gene, spec in fp._RECESSIVE_CARRIER_FREQ.items():
        assert "European" in spec["freq"], f"{gene} lacks the European fallback key"


def test_penetrance_is_an_ordered_pair_in_range():
    for gene, spec in fp._RECESSIVE_CARRIER_FREQ.items():
        lo, hi = spec["penetrance"]
        assert 0.0 < lo <= hi <= 1.0, f"{gene} penetrance {(lo, hi)} is malformed"


def test_every_entry_has_partner_test_and_note():
    for gene, spec in fp._RECESSIVE_CARRIER_FREQ.items():
        assert spec.get("partner_test"), f"{gene} missing partner_test"
        assert len(spec.get("note", "")) > 40, f"{gene} note is too thin"


def test_headline_screening_conditions_are_covered():
    # The conditions expanded carrier screening is actually built around.
    for gene in ("CFTR", "HBB", "SMN1", "HEXA"):
        assert gene in fp._RECESSIVE_CARRIER_FREQ


def test_dominant_genes_are_not_given_carrier_frequencies():
    # Partner carrier frequency is meaningless for a dominant condition; adding
    # one would imply the wrong reproductive model.
    for gene in ("BRCA1", "BRCA2", "MLH1", "FBN1"):
        assert gene not in fp._RECESSIVE_CARRIER_FREQ


def test_ancestry_specific_frequencies_differ_where_known():
    # Sanity that the table carries real ancestry structure rather than one
    # number copied across keys.
    cftr = fp._RECESSIVE_CARRIER_FREQ["CFTR"]["freq"]
    assert cftr["European"] > cftr["East Asian"]
    hbb = fp._RECESSIVE_CARRIER_FREQ["HBB"]["freq"]
    assert hbb["African American"] > hbb["European"]
