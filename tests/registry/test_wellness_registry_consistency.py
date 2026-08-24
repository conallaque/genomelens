"""Cross-module consistency for wellness/core.py vs snp_registry."""

from __future__ import annotations

from wellness import core as wellness


def test_audit_returns_well_structured_dict() -> None:
    audit = wellness.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}
    assert isinstance(audit["registered"], list)
    assert isinstance(audit["missing"], list)


def test_at_least_25_wellness_rsids_registered() -> None:
    """Lock the V8 migration baseline so nobody accidentally pulls a
    registered rsID out (which would silently regress to wellness.py's old
    local-knowledge-only state)."""
    audit = wellness.audit_against_registry()
    assert len(audit["registered"]) >= 25


def test_known_deferred_rsids_match_changelog() -> None:
    """The wellness migration deferred rare antioxidant / longevity variants
    to V8.1 per CHANGELOG.md. Lock that set. (rs25531 was migrated into the
    registry in V6.12 for the neurochemistry module's 5-HTTLPR proxy.)"""
    audit = wellness.audit_against_registry()
    expected_deferred = {
        "rs12934922", "rs174546", "rs1799750", "rs1799752",
        "rs73598374", "rs7501331", "rs77086077", "rs9420907",
    }
    assert set(audit["missing"]) == expected_deferred, (
        f"Deferred-to-V8.1 set drifted from CHANGELOG. "
        f"Got {set(audit['missing'])}, expected {expected_deferred}. "
        "If a new rsID was added to wellness.py without registering it, "
        "add it to snp_registry._RECORDS. If a rsID was migrated, update "
        "this test + CHANGELOG together."
    )
