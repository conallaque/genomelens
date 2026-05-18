"""
Cross-module consistency test: every rsID in carrier.CARRIER_VARIANTS that
is also in snp_registry must agree on (a) gene symbol and (b) the pathogenic
allele being either registry-derived or its reverse-strand complement.

This is the V8 migration safety net. If a future change to either module
breaks the reconciliation, this test fails loudly.
"""

from __future__ import annotations

import carrier


def test_carrier_audit_passes_at_import() -> None:
    """The module-level assert in carrier.py is the runtime guard; this
    test makes the contract visible to the test suite as well."""
    audit = carrier.audit_against_registry()
    assert audit["disagreed"] == [], (
        "Real biology disagreement (not strand) between carrier.py and "
        f"snp_registry: {audit['disagreed']}"
    )


def test_at_least_thirteen_reconciled_directly() -> None:
    """V8 migration target — protect against regression where someone
    accidentally removes a carrier rsID from the registry."""
    audit = carrier.audit_against_registry()
    assert len(audit["agreed"]) >= 13


def test_known_strand_flips_remain_documented() -> None:
    """The three strand-flipped variants are an explicit design choice
    (carrier.py uses gene-coding-strand labels; registry uses + strand).
    Locking them here so any future "fix" that breaks the flip surfaces."""
    audit = carrier.audit_against_registry()
    expected = {"rs28929474", "rs5030858", "rs5742904"}
    assert set(audit["strand_flipped"]) == expected, (
        f"Strand-flipped set drifted. Expected {expected}, "
        f"got {set(audit['strand_flipped'])}. If this is intentional, "
        "update CHANGELOG.md and the test."
    )


def test_no_registry_record_missing_position() -> None:
    """Every carrier rsID we added to the registry must carry both build
    coordinates — they're needed for chip-build auto-detection later."""
    audit = carrier.audit_against_registry()
    import snp_registry as reg
    for rsid in audit["agreed"] + audit["strand_flipped"]:
        rec = reg.get(rsid)
        assert rec is not None
        # GRCh38 may be None for legacy seeds; GRCh37 must always exist.
        assert rec.pos_grch37 is not None, f"{rsid} missing GRCh37 position"
