"""Cross-module consistency for traits.py vs snp_registry.

V8 status: 16 / 48 traits.py rsIDs registered (via earlier modules that
share their SNPs). The remaining 32 are documented in CHANGELOG.md under
V8.1 follow-ups — each is a well-characterised phenotype SNP (eye/hair
colour, taste, earwax) that needs literature-cited ancestral/derived +
GRCh37/38 coordinates before joining the registry.
"""

from __future__ import annotations

import traits


def test_audit_shape() -> None:
    audit = traits.audit_against_registry()
    assert set(audit.keys()) == {"registered", "missing"}


def test_at_least_16_traits_rsids_registered() -> None:
    """V8 baseline — protect against regression where a future change
    de-registers a rsID that traits.py relies on."""
    audit = traits.audit_against_registry()
    assert len(audit["registered"]) >= 16


def test_total_rsids_is_48() -> None:
    """If this count changes, traits.py gained / lost a rule; the V8.1
    deferred-list in CHANGELOG.md should be updated to match."""
    audit = traits.audit_against_registry()
    total = len(audit["registered"]) + len(audit["missing"])
    assert total == 48, (
        f"traits.py now references {total} rsIDs (was 48 at V8 cut). "
        "Update CHANGELOG.md V8.1 follow-ups + this test together."
    )
