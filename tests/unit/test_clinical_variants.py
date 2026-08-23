"""Unit tests for clinical_variants.py (Phase-2 ClinVar screen)."""

from __future__ import annotations

import pytest

import clinical_variants as cv

# ── pure helpers ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("sig,expected", [
    ("Pathogenic", True), ("Likely_pathogenic", True),
    ("Pathogenic/Likely_pathogenic", True),
    ("Uncertain_significance", False), ("Benign", False),
    ("Likely_benign", False),
    ("Conflicting_classifications_of_pathogenicity", False),
    ("", False),
])
def test_is_pathogenic_sig(sig, expected):
    assert cv.is_pathogenic_sig(sig) is expected


@pytest.mark.parametrize("rev,stars", [
    ("practice_guideline", 4),
    ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
    ("criteria_provided,_single_submitter", 1),
    ("no_assertion_criteria_provided", 0),
    ("nonsense", 0),
])
def test_revstat_to_stars(rev, stars):
    assert cv.revstat_to_stars(rev) == stars


def test_norm_key_snv_identity():
    assert cv.norm_key("chr17", 100, "C", "T") == ("17", 100, "C", "T")


def test_norm_key_indel_trims_anchor_base():
    # AG>A (del) and a right-anchored representation collapse to the same key
    assert cv.norm_key("1", 100, "AG", "A") == cv.norm_key("1", 100, "AGG", "AG")


@pytest.mark.parametrize("gt,alt_idx,expected", [
    ("0/1", 1, "heterozygous"),
    ("1/1", 1, "homozygous"),
    ("1/2", 1, "heterozygous"),   # matched ALT1 in a multiallelic sample
    ("1/2", 2, "heterozygous"),
    ("1", 1, "hemizygous"),
    ("0/0", 1, None),             # doesn't carry the allele
    ("./.", 1, None),
    ("0/1:35", 1, "heterozygous"),
])
def test_zygosity_for_alt(gt, alt_idx, expected):
    assert cv.zygosity_for_alt(gt, alt_idx) == expected


# ── end-to-end distill + analyze against synthetic fixtures ───────────────────

_CLINVAR_ROWS = [
    ("17", "41250000", "C", "T", "Pathogenic",
     "criteria_provided,_multiple_submitters,_no_conflicts", "BRCA1", "HBOC", "single_nucleotide_variant"),
    ("7", "117199644", "A", "G", "Likely_pathogenic",
     "criteria_provided,_single_submitter", "CFTR", "Cystic_fibrosis", "single_nucleotide_variant"),
    ("7", "117199700", "G", "T", "Pathogenic",
     "reviewed_by_expert_panel", "CFTR", "Cystic_fibrosis", "single_nucleotide_variant"),
    ("15", "72346580", "C", "A", "Pathogenic",
     "criteria_provided,_multiple_submitters,_no_conflicts", "HEXA", "Tay-Sachs", "single_nucleotide_variant"),
    ("1", "45797000", "AG", "A", "Pathogenic",
     "criteria_provided,_single_submitter", "MUTYH", "MAP", "Deletion"),
    ("17", "7670000", "G", "A", "Uncertain_significance",
     "criteria_provided,_single_submitter", "TP53", "Li-Fraumeni", "single_nucleotide_variant"),
    ("2", "1000", "C", "G", "Pathogenic",
     "no_assertion_criteria_provided", "FOO", "x", "single_nucleotide_variant"),
    ("2", "2000", "C", "G", "Benign",
     "criteria_provided,_single_submitter", "BAR", "y", "single_nucleotide_variant"),
]


@pytest.fixture()
def clinvar_table(tmp_path, monkeypatch):
    raw = tmp_path / "clinvar_raw.vcf"
    with open(raw, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a, sig, rev, gene, dn, vc in _CLINVAR_ROWS:
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t"
                    f"CLNSIG={sig};CLNREVSTAT={rev};GENEINFO={gene}:1;CLNDN={dn};CLNVC={vc}\n")
    cvar_dir = tmp_path / "clinvar"
    cvar_dir.mkdir()
    dist = cvar_dir / "clinvar_plp_grch38.tsv.gz"
    n = cv.distill_clinvar_vcf(str(raw), str(dist), log=None)
    monkeypatch.setattr(cv, "CLINVAR_DIR", cvar_dir)
    return n


def _write_user_vcf(tmp_path, variants):
    p = tmp_path / "user.vcf"
    with open(p, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for c, pos, r, a, gt in variants:
            f.write(f"{c}\t{pos}\t.\t{r}\t{a}\t.\t.\t.\tGT\t{gt}\n")
    return str(p)


def test_distill_keeps_plp_and_acmg_vus_only(clinvar_table):
    # 5 P/LP + 1 VUS-in-ACMG(TP53) + 1 zero-star P/LP(FOO) = 7; Benign dropped
    assert clinvar_table == 7


def test_full_screen_classifies_every_branch(clinvar_table, tmp_path):
    user = _write_user_vcf(tmp_path, [
        ("17", "41250000", "C", "T", "0/1"),   # BRCA1 het -> actionable
        ("7", "117199644", "A", "G", "0/1"),    # CFTR het #1
        ("7", "117199700", "G", "T", "0/1"),    # CFTR het #2 -> compound-het
        ("15", "72346580", "C", "A", "1/1"),    # HEXA hom -> affected
        ("1", "45797000", "AG", "A", "0/1"),    # MUTYH indel het -> carrier
        ("17", "7670000", "G", "A", "0/1"),     # TP53 VUS
        ("2", "1000", "C", "G", "0/1"),         # 0-star -> excluded
        ("9", "5000", "A", "G", "0/1"),         # no match
    ])
    r = cv.analyze_clinical_variants(user, "grch38", inferred_sex="M")
    assert r["available"]
    b = r["buckets"]
    assert len(b["actionable"]) == 1 and b["actionable"][0]["gene"] == "BRCA1"
    assert len(b["affected"]) == 1 and b["affected"][0]["gene"] == "HEXA"
    assert len(b["possible_compound_het"]) == 2  # both CFTR variants flagged
    assert len(b["carrier"]) == 1 and b["carrier"][0]["gene"] == "MUTYH"
    assert r["n_vus_in_acmg"] == 1               # TP53 VUS counted, not a finding
    assert r["n_excluded_0star"] == 1            # FOO excluded
    assert r["negative_disclaimer"]


def test_build_gate_refuses_unknown(clinvar_table, tmp_path):
    user = _write_user_vcf(tmp_path, [("17", "41250000", "C", "T", "0/1")])
    r = cv.analyze_clinical_variants(user, "unknown")
    assert r["available"] is False


def test_missing_table_degrades_gracefully(tmp_path, monkeypatch):
    monkeypatch.setattr(cv, "CLINVAR_DIR", tmp_path / "nonexistent")
    user = _write_user_vcf(tmp_path, [("1", "1", "C", "T", "0/1")])
    r = cv.analyze_clinical_variants(user, "grch38")
    assert r["available"] is False
    assert "setup.py" in r["reason"]
    assert r["negative_disclaimer"]


def test_freshness_fields_surface_and_flag_stale(clinvar_table, tmp_path, monkeypatch):
    import datetime
    import json
    # write a 60-day-old meta sidecar next to the distilled table
    meta = cv.CLINVAR_DIR / "clinvar_plp_grch38.meta.json"
    old = (datetime.datetime.now() - datetime.timedelta(days=60)).isoformat(timespec="seconds")
    meta.write_text(json.dumps({"distilled": old, "source_last_modified": "x", "rows": 7}))
    user = _write_user_vcf(tmp_path, [("17", "41250000", "C", "T", "0/1")])
    r = cv.analyze_clinical_variants(user, "grch38")
    assert r["available"]
    assert r["clinvar_date"] == old
    assert r["clinvar_stale"] is True     # >45 days → stale


def test_negative_result_disclaimer_always_present(clinvar_table, tmp_path):
    user = _write_user_vcf(tmp_path, [("9", "999", "A", "G", "0/1")])  # no matches
    r = cv.analyze_clinical_variants(user, "grch38")
    assert r["available"] and r["n_plp"] == 0
    assert "not" in r["negative_disclaimer"].lower()
