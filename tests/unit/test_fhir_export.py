"""
Unit tests for fhir_export.py.

The V7 upgrade wrapped the existing Observations in a proper CG-IG
DiagnosticReport and added a Provenance resource. Lock those invariants here.
"""

from __future__ import annotations

import json
from pathlib import Path

import fhir_export as fx


# ── Stub upstream-module data ────────────────────────────────────────────────

_PGX_STUB = {
    "per_gene": {
        "CYP2D6": {
            "long_name": "Cytochrome P450 2D6",
            "phenotype": "Intermediate Metaboliser", "phenotype_code": "IM",
            "activity_score": 1.25,
            "callable_variants": 3, "total_variants": 5, "callability_pct": 60.0,
            "variant_calls": [
                {"rsid": "rs1065852", "called": True, "dosage": 1, "genotype": "GA"},
                {"rsid": "rs3892097", "called": True, "dosage": 0, "genotype": "CC"},
            ],
            "cpic_guideline": "CPIC v2024.1",
            "is_binary": False,
        },
    },
    "actionable_findings": [{
        "gene": "CYP2D6", "phenotype": "Intermediate Metaboliser",
        "phenotype_class": "pheno-im", "drug": "codeine",
        "recommendation": "Use alternative analgesic",
        "guideline": "CPIC v2024.1",
    }],
}

_CARRIER_STUB = {
    "affected": [],
    "carriers": [{
        "rsid": "rs6025", "gene": "F5",
        "variant": "Factor V Leiden (R506Q)", "disease": "VTE",
        "dosage": 1, "genotype": "GA",
        "carrier_implication": "1 copy raises VTE risk ~5-7×.",
    }],
}

_HLA_STUB = {"alleles": [
    {"allele": "HLA-B*27", "status": "carrier (heterozygous)",
     "confidence": "moderate", "dosage": 1},
]}


# ── Bundle structure ─────────────────────────────────────────────────────────

def test_bundle_has_diagnostic_report_envelope() -> None:
    """CG-IG mandates DiagnosticReport as the top-level genomic-report resource."""
    result = fx.build_fhir_bundle(
        pgx_result=_PGX_STUB, carrier_result=_CARRIER_STUB,
        hla_result=_HLA_STUB, apoe_genotype="ε3/ε4",
    )
    types = [e["resource"]["resourceType"] for e in result["bundle"]["entry"]]
    assert "DiagnosticReport" in types
    assert types.index("Patient") < types.index("DiagnosticReport"), \
        "DiagnosticReport should come immediately after Patient"


def test_bundle_has_provenance() -> None:
    result = fx.build_fhir_bundle(apoe_genotype="ε3/ε3", file_label="test.csv")
    types = [e["resource"]["resourceType"] for e in result["bundle"]["entry"]]
    assert types[-1] == "Provenance", "Provenance should be the last entry"


def test_diagnostic_report_references_all_observations() -> None:
    result = fx.build_fhir_bundle(
        pgx_result=_PGX_STUB, carrier_result=_CARRIER_STUB,
        hla_result=_HLA_STUB, apoe_genotype="ε3/ε4",
    )
    bundle = result["bundle"]
    obs_count = sum(1 for e in bundle["entry"]
                    if e["resource"]["resourceType"] == "Observation")
    diag = next(e["resource"] for e in bundle["entry"]
                if e["resource"]["resourceType"] == "DiagnosticReport")
    assert len(diag["result"]) == obs_count


def test_provenance_has_clinical_scope_policy() -> None:
    """The Provenance.policy[] should explicitly state what *isn't* in the bundle."""
    result = fx.build_fhir_bundle(apoe_genotype="ε3/ε3")
    prov = next(e["resource"] for e in result["bundle"]["entry"]
                if e["resource"]["resourceType"] == "Provenance")
    assert prov["policy"]
    assert any("clinically actionable" in p.lower() for p in prov["policy"])
    assert any("excluded" in p.lower() or "limit" in p.lower() for p in prov["policy"])


def test_summary_counts_match_emitted_resources() -> None:
    result = fx.build_fhir_bundle(
        pgx_result=_PGX_STUB, carrier_result=_CARRIER_STUB,
        hla_result=_HLA_STUB, apoe_genotype="ε3/ε4",
    )
    s = result["summary"]
    bundle = result["bundle"]
    assert s["n_pgx_observations"] == 1
    assert s["n_pgx_actionable"] == 1
    assert s["n_carrier_observations"] == 1
    assert s["n_hla_observations"] == 1
    assert s["apoe_included"] is True
    assert s["has_diagnostic_report"] is True
    assert s["has_provenance"] is True
    # Total entries: Patient + DiagnosticReport + 1 APOE + 1 PGx + 1 PGx-drug
    #                + 1 carrier + 1 HLA + Provenance = 8
    assert s["n_total_entries"] == len(bundle["entry"]) == 8


def test_carrier_with_zero_findings_omits_obs_but_keeps_report() -> None:
    result = fx.build_fhir_bundle(apoe_genotype="ε3/ε3")
    types = [e["resource"]["resourceType"] for e in result["bundle"]["entry"]]
    # Patient + DiagnosticReport + APOE Observation + Provenance
    assert types.count("Observation") == 1
    assert "DiagnosticReport" in types
    assert "Provenance" in types


# ── On-disk export ───────────────────────────────────────────────────────────

def test_export_fhir_writes_valid_json(tmp_path) -> None:
    out = tmp_path / "fhir.json"
    summary = fx.export_fhir(
        out_path=out,
        pgx_result=_PGX_STUB,
        apoe_genotype="ε3/ε4",
    )
    assert out.exists()
    data = json.loads(out.read_text())
    assert data["resourceType"] == "Bundle"
    assert data["type"] == "collection"
    assert summary["has_diagnostic_report"] is True
