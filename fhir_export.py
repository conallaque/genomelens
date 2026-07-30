"""
HL7 FHIR R4 Clinical Export
===========================

Exports clinically *actionable finding types* in FHIR R4 Bundle format —
format-compatible with EHR ingestion for review by a genetic counsellor or
physician, NOT a certified clinical record. Chip/genome-derived; requires
confirmatory clinical testing before any clinical action. Specifically:

  • PGx phenotypes        → Observation (LOINC + CPIC coding)
  • Carrier status        → Observation (genomic, ACMG-style coded variants)
  • HLA risk alleles      → Observation (HLA serotype)
  • APOE genotype         → Observation (LOINC 48006-2 / 48005-3)

We deliberately exclude PRS, trait predictions, ancestry, wellness, exercise,
nutrition, supplement, and bloodwork output — none of those are clinically
validated for individual decision-making per CPIC/ACMG/ESHG guidance.

Output is a single FHIR R4 Bundle (type=collection) JSON file suitable for
upload via the FHIR Bulk Import endpoint of any compliant EHR (Epic Cerner,
Athena, Allscripts, etc.).

Reference:
  HL7 FHIR R4 https://hl7.org/fhir/R4/
  CG-IG Clinical Genomics implementation guide https://hl7.org/fhir/uv/genomics-reporting/
"""

from __future__ import annotations

import datetime
import json
import uuid
from pathlib import Path
from typing import Dict, List, Optional


FHIR_VERSION = "4.0.1"
CG_IG_VERSION = "2.0.0"  # Clinical Genomics IG


# ── Helpers ──────────────────────────────────────────────────────────────────

def _uid() -> str:
    """FHIR-style URN UUID."""
    return f"urn:uuid:{uuid.uuid4()}"


def _now_iso() -> str:
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _coding(system: str, code: str, display: str) -> Dict:
    return {"system": system, "code": code, "display": display}


def _codeable(*codings: Dict, text: str = "") -> Dict:
    cc: Dict = {"coding": list(codings)}
    if text:
        cc["text"] = text
    return cc


def _patient(patient_id: str, sex: Optional[str]) -> Dict:
    res: Dict = {
        "resourceType": "Patient",
        "id": patient_id,
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/Patient-uv-genomics-reporting"
        ]},
        "active": True,
    }
    if sex:
        res["gender"] = "male" if sex.lower().startswith("m") else "female"
    return res


# ── PGx Observation builders ────────────────────────────────────────────────

def _pgx_observation(patient_ref: str, gene_name: str, gene_data: Dict) -> Dict:
    """Build a Genotype/Phenotype Observation for one PGx gene per CG-IG."""
    obs_id = _uid()

    # Map common pheno labels to CPIC standard terms where possible
    phenotype = gene_data.get("phenotype", "")
    activity = gene_data.get("activity_score")
    is_binary = gene_data.get("is_binary", False)

    components: List[Dict] = []

    # Activity score component (where applicable)
    if activity is not None and not is_binary:
        components.append({
            "code": _codeable(
                _coding("http://loinc.org", "82120-4", "Genetic activity score"),
                text="CPIC activity score",
            ),
            "valueQuantity": {
                "value": activity, "unit": "score",
                "system": "http://unitsofmeasure.org", "code": "1",
            },
        })

    # Per-variant genotype components
    for call in gene_data.get("variant_calls", []):
        if not call.get("called"):
            continue
        components.append({
            "code": _codeable(
                _coding("http://loinc.org", "53034-5", "Allelic state"),
                _coding("http://www.ncbi.nlm.nih.gov/snp", call["rsid"], call["rsid"]),
                text=f"{call['rsid']} genotype",
            ),
            "valueCodeableConcept": _codeable(
                _coding("http://loinc.org", "LA6705-3",
                        "Heterozygous" if call.get("dosage", 0) == 1 else
                        ("Homozygous" if call.get("dosage", 0) == 2 else "Wild-type")),
                text=call.get("genotype", ""),
            ),
        })

    res = {
        "resourceType": "Observation",
        "id": obs_id.replace("urn:uuid:", ""),
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genotype"
        ]},
        "status": "final",
        "category": [_codeable(
            _coding("http://terminology.hl7.org/CodeSystem/observation-category",
                    "laboratory", "Laboratory"),
            _coding("http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                    "GE", "Genetic"),
        )],
        "code": _codeable(
            _coding("http://loinc.org", "84413-4", "Genotype display name"),
            text=f"{gene_name} phenotype: {phenotype}",
        ),
        "subject": {"reference": patient_ref},
        "effectiveDateTime": _now_iso(),
        "valueCodeableConcept": _codeable(
            _coding("https://cpicpgx.org/lookup",
                    gene_data.get("phenotype_code", "UNK"), phenotype),
            text=phenotype,
        ),
        "component": components,
    }

    # Add cpic guideline as note
    if gene_data.get("cpic_guideline"):
        res["note"] = [{"text": f"CPIC: {gene_data['cpic_guideline']}"}]

    return res


def _pgx_drug_recommendation_observations(patient_ref: str,
                                          actionable: List[Dict]) -> List[Dict]:
    """One MedicationStatement-style Observation per actionable drug rec."""
    out: List[Dict] = []
    for f in actionable:
        out.append({
            "resourceType": "Observation",
            "id": uuid.uuid4().hex,
            "meta": {"profile": [
                "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/medication-assessed"
            ]},
            "status": "final",
            "category": [_codeable(_coding(
                "http://terminology.hl7.org/CodeSystem/observation-category",
                "laboratory", "Laboratory"))],
            "code": _codeable(
                _coding("http://loinc.org", "51963-7",
                        "Medication assessed [ID]"),
                text=f"{f['gene']} → {f['drug']} ({f['phenotype']})",
            ),
            "subject": {"reference": patient_ref},
            "effectiveDateTime": _now_iso(),
            "valueString": f["recommendation"],
            "note": [{"text": f"Guideline: {f.get('guideline','CPIC')}"}],
        })
    return out


# ── Carrier Observation builder ─────────────────────────────────────────────

def _carrier_observation(patient_ref: str, record: Dict, status: str) -> Dict:
    """
    status: 'affected' (homozygous), 'carrier' (heterozygous), or
            'not-carrier' (homozygous reference).
    """
    dosage = record.get("dosage")
    if dosage == 2:
        zygosity_code, zygosity_display = "LA6705-3", "Homozygous"
        clin_sig = "Pathogenic — biallelic"
    elif dosage == 1:
        zygosity_code, zygosity_display = "LA6706-1", "Heterozygous"
        clin_sig = "Pathogenic — carrier (one allele)"
    else:
        zygosity_code, zygosity_display = "LA6704-6", "Homozygous reference"
        clin_sig = "No pathogenic allele"

    return {
        "resourceType": "Observation",
        "id": uuid.uuid4().hex,
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"
        ]},
        "status": "final",
        "category": [_codeable(_coding(
            "http://terminology.hl7.org/CodeSystem/observation-category",
            "laboratory", "Laboratory"))],
        "code": _codeable(
            _coding("http://loinc.org", "69548-6", "Genetic variant assessment"),
            text=f"{record.get('gene','')} {record.get('variant','')} — {record.get('disease','')}",
        ),
        "subject": {"reference": patient_ref},
        "effectiveDateTime": _now_iso(),
        "valueCodeableConcept": _codeable(
            _coding("http://loinc.org", "LA9633-4", "Present"),
            text=clin_sig,
        ) if dosage and dosage > 0 else _codeable(
            _coding("http://loinc.org", "LA9634-2", "Absent"),
            text=clin_sig,
        ),
        "component": [
            {
                "code": _codeable(_coding("http://loinc.org", "48018-6",
                                          "Gene studied [ID]")),
                "valueCodeableConcept": _codeable(text=record.get("gene", "")),
            },
            {
                "code": _codeable(_coding("http://loinc.org", "48005-3",
                                          "Amino acid change [Type] in p.HGVS")),
                "valueCodeableConcept": _codeable(text=record.get("variant", "")),
            },
            {
                "code": _codeable(_coding("http://loinc.org", "53034-5",
                                          "Allelic state")),
                "valueCodeableConcept": _codeable(_coding(
                    "http://loinc.org", zygosity_code, zygosity_display)),
            },
            {
                "code": _codeable(_coding("http://loinc.org", "81252-9",
                                          "Discrete genetic variant")),
                "valueCodeableConcept": _codeable(_coding(
                    "http://www.ncbi.nlm.nih.gov/snp",
                    record.get("rsid", ""), record.get("rsid", ""))),
            },
        ],
        "note": [{"text": (
            record.get("affected_implication") if dosage == 2 else
            record.get("carrier_implication", "")
        )}],
    }


# ── HLA Observation builder ─────────────────────────────────────────────────

def _hla_observation(patient_ref: str, allele_result: Dict) -> Dict:
    return {
        "resourceType": "Observation",
        "id": uuid.uuid4().hex,
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"
        ]},
        "status": "final",
        "category": [_codeable(_coding(
            "http://terminology.hl7.org/CodeSystem/observation-category",
            "laboratory", "Laboratory"))],
        "code": _codeable(
            _coding("http://loinc.org", "81244-6", "HLA allele [Type]"),
            text=f"HLA {allele_result['allele']} (tag-SNP imputed)",
        ),
        "subject": {"reference": patient_ref},
        "effectiveDateTime": _now_iso(),
        "valueCodeableConcept": _codeable(
            _coding("http://loinc.org", "LA9633-4", "Present"),
            text=f"{allele_result['allele']} status: {allele_result.get('status')}",
        ),
        "note": [{"text": (
            f"Confidence: {allele_result.get('confidence', 'unknown')}. "
            f"Tag-SNP imputation, not direct typing. "
            f"For transplant matching or definitive abacavir screening, "
            f"order serotype-level HLA typing."
        )}],
    }


# ── APOE genotype Observation ───────────────────────────────────────────────

def _apoe_observation(patient_ref: str, apoe_genotype: str) -> Dict:
    return {
        "resourceType": "Observation",
        "id": uuid.uuid4().hex,
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genotype"
        ]},
        "status": "final",
        "category": [_codeable(_coding(
            "http://terminology.hl7.org/CodeSystem/observation-category",
            "laboratory", "Laboratory"))],
        "code": _codeable(
            _coding("http://loinc.org", "48006-2",
                    "APO E gene targeted mutation analysis"),
            text="APOE Alzheimer's risk genotype",
        ),
        "subject": {"reference": patient_ref},
        "effectiveDateTime": _now_iso(),
        "valueCodeableConcept": _codeable(text=apoe_genotype),
        "note": [{"text": (
            "APOE genotype is a clinically recognised risk modifier for "
            "late-onset Alzheimer's disease. ε4/ε4 carries highest lifetime "
            "risk; ε2/ε2 and ε2/ε3 may be protective."
        )}],
    }


# ── DiagnosticReport (top-level CG-IG resource) ─────────────────────────────

def _diagnostic_report(
    patient_ref: str,
    observation_refs: List[Dict],
    summary_text: str,
) -> Dict:
    """
    Per the HL7 Clinical Genomics IG, a genomic report is a DiagnosticReport
    that *contains* (via `result` references) the per-finding Observations.
    Previous bundle emitted the Observations directly into a collection; the
    DiagnosticReport is the envelope a clinician's EHR expects to ingest.
    """
    return {
        "resourceType": "DiagnosticReport",
        "id": uuid.uuid4().hex,
        "meta": {"profile": [
            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/"
            "diagnosticreport-uv-genomics-reporting"
        ]},
        "status": "final",
        "category": [_codeable(
            _coding("http://terminology.hl7.org/CodeSystem/v2-0074",
                    "GE", "Genetics"),
        )],
        "code": _codeable(
            _coding("http://loinc.org", "81247-9",
                    "Master HL7 genetic variant reporting panel"),
            text="Consumer-chip-derived clinical genomics report",
        ),
        "subject": {"reference": patient_ref},
        "effectiveDateTime": _now_iso(),
        "issued": _now_iso(),
        "performer": [{
            "display": "Local DNA Analysis Tool (consumer-chip pipeline)",
        }],
        "result": observation_refs,
        "conclusion": summary_text,
        "conclusionCode": [_codeable(
            _coding("http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                    "GS", "Genetic susceptibility"),
        )],
    }


def _provenance(
    target_refs: List[Dict],
    file_label: str,
) -> Dict:
    """
    A Provenance resource records the chain of custody / derivation: where the
    data came from (the chip raw file), what software produced the
    Observations, and when. EHRs use this to attribute findings.
    """
    return {
        "resourceType": "Provenance",
        "id": uuid.uuid4().hex,
        "target": target_refs,
        "recorded": _now_iso(),
        "occurredDateTime": _now_iso(),
        "agent": [{
            "type": _codeable(_coding(
                "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
                "assembler", "Assembler",
            )),
            "who": {"display": "Local DNA Analysis Tool"},
        }],
        "entity": [{
            "role": "source",
            "what": {
                "display": (
                    f"Consumer-chip raw genotype file "
                    f"({file_label or 'unspecified'}) — autosomal SNP array "
                    f"output, GRCh37/hg19 coordinates"
                ),
            },
        }],
        "policy": [
            # Statement of clinical scope
            "Findings limited to clinically actionable finding TYPES (not "
            "clinically validated results): PGx (CPIC), carrier status "
            "(autosomal recessive screening panels), HLA tag-SNP imputation, "
            "and APOE genotype. Chip/genome-derived; requires confirmatory "
            "clinical testing before any clinical action. PRS, ancestry, "
            "traits, wellness, and lifestyle modules are excluded from this "
            "clinical bundle.",
        ],
    }


# ── Bundle builder ──────────────────────────────────────────────────────────

def build_fhir_bundle(
    pgx_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    hla_result: Optional[Dict] = None,
    apoe_genotype: Optional[str] = None,
    patient_id: str = "patient-1",
    inferred_sex: Optional[str] = None,
    file_label: str = "",
) -> Dict:
    """Assemble a FHIR R4 Bundle of clinically-validated findings only."""
    patient_ref = f"Patient/{patient_id}"
    entries: List[Dict] = []
    obs_refs: List[Dict] = []  # collected for DiagnosticReport.result

    # Patient resource
    patient_res = _patient(patient_id, inferred_sex)
    entries.append({
        "fullUrl": f"urn:uuid:{patient_id}",
        "resource": patient_res,
    })

    def _add(res: Dict) -> None:
        full_url = f"urn:uuid:{res['id']}"
        entries.append({"fullUrl": full_url, "resource": res})
        obs_refs.append({"reference": f"{res['resourceType']}/{res['id']}"})

    # APOE
    if apoe_genotype:
        _add(_apoe_observation(patient_ref, apoe_genotype))

    # PGx genotype + phenotype observations
    n_pgx = 0
    n_pgx_actionable = 0
    if pgx_result and pgx_result.get("per_gene"):
        for gene_name, gene_data in pgx_result["per_gene"].items():
            if gene_data.get("callable_variants", 0) == 0:
                continue
            _add(_pgx_observation(patient_ref, gene_name, gene_data))
            n_pgx += 1

        # Actionable drug recommendations
        for drug_obs in _pgx_drug_recommendation_observations(
            patient_ref, pgx_result.get("actionable_findings", [])
        ):
            _add(drug_obs)
            n_pgx_actionable += 1

    # Carrier (affected + carrier only — not_carrier and untested skipped)
    n_carrier = 0
    if carrier_result:
        for rec in carrier_result.get("affected", []):
            _add(_carrier_observation(patient_ref, rec, "affected"))
            n_carrier += 1
        for rec in carrier_result.get("carriers", []):
            _add(_carrier_observation(patient_ref, rec, "carrier"))
            n_carrier += 1

    # HLA: only carrier/homozygous called alleles
    n_hla = 0
    if hla_result and hla_result.get("alleles"):
        for allele in hla_result["alleles"]:
            if allele.get("status") in ("carrier (heterozygous)", "homozygous"):
                _add(_hla_observation(patient_ref, allele))
                n_hla += 1

    # Build a human-readable conclusion paragraph for the DiagnosticReport
    conclusion_parts: List[str] = []
    if apoe_genotype:
        conclusion_parts.append(f"APOE genotype: {apoe_genotype}.")
    if n_pgx:
        conclusion_parts.append(
            f"{n_pgx} pharmacogenomic genes evaluated; "
            f"{n_pgx_actionable} actionable drug-level guidance entries."
        )
    if n_carrier:
        conclusion_parts.append(
            f"{n_carrier} pathogenic-allele carrier/affected finding(s)."
        )
    if n_hla:
        conclusion_parts.append(
            f"{n_hla} HLA risk allele(s) imputed via tag-SNP method "
            f"(not direct typing)."
        )
    if not conclusion_parts:
        conclusion_parts.append("No clinically actionable finding types reported.")
    conclusion = " ".join(conclusion_parts)

    # DiagnosticReport — proper top-level CG-IG resource
    diag_report = _diagnostic_report(patient_ref, list(obs_refs), conclusion)
    # Place DiagnosticReport immediately after Patient for ingestion order
    entries.insert(1, {
        "fullUrl": f"urn:uuid:{diag_report['id']}",
        "resource": diag_report,
    })

    # Provenance — chain of custody for ALL resources we just emitted
    all_refs = [
        {"reference": f"{e['resource']['resourceType']}/{e['resource']['id']}"}
        for e in entries
        if e['resource']['resourceType'] != 'Provenance'
        and e['resource'].get('id')
    ]
    prov = _provenance(all_refs, file_label)
    entries.append({"fullUrl": f"urn:uuid:{prov['id']}", "resource": prov})

    bundle = {
        "resourceType": "Bundle",
        "id": uuid.uuid4().hex,
        "meta": {
            "lastUpdated": _now_iso(),
            "profile": [
                "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/"
                "diagnosticreport-uv-genomics-reporting"
            ],
            "tag": [{
                "system": "https://hl7.org/fhir/R4",
                "code": "consumer-genomics",
                "display": "Generated from consumer chip raw data — clinically "
                           "actionable finding types only, not clinically "
                           "validated results",
            }],
        },
        "type": "collection",
        "timestamp": _now_iso(),
        "entry": entries,
    }

    summary = {
        "fhir_version": FHIR_VERSION,
        "cg_ig_version": CG_IG_VERSION,
        "n_pgx_observations": n_pgx,
        "n_pgx_actionable": n_pgx_actionable,
        "n_carrier_observations": n_carrier,
        "n_hla_observations": n_hla,
        "apoe_included": bool(apoe_genotype),
        "n_total_entries": len(entries),
        "has_diagnostic_report": True,
        "has_provenance": True,
        "file_label": file_label,
    }
    return {"bundle": bundle, "summary": summary}


def export_fhir(
    out_path: Path,
    pgx_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    hla_result: Optional[Dict] = None,
    apoe_genotype: Optional[str] = None,
    inferred_sex: Optional[str] = None,
    file_label: str = "",
) -> Dict:
    """Write the FHIR bundle to disk; return the summary dict."""
    result = build_fhir_bundle(
        pgx_result=pgx_result,
        carrier_result=carrier_result,
        hla_result=hla_result,
        apoe_genotype=apoe_genotype,
        inferred_sex=inferred_sex,
        file_label=file_label,
    )
    out_path = Path(out_path)
    out_path.write_text(json.dumps(result["bundle"], indent=2))
    return result["summary"]
