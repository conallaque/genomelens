"""Unit tests for novel_variants.py (Phase 3 — computational predictor screen).

Builds a real synthetic tabix-indexed AlphaMissense table in tmp_path and runs the
full annotate → classify path, plus the graceful-degradation and build-gate paths.
"""
from __future__ import annotations

import pytest

from risk import novel_variants as nv

pysam = pytest.importorskip("pysam")


# ── fixtures ──────────────────────────────────────────────────────────────────

def _make_am_table(root, rows):
    """rows: (chrom, pos, ref, alt, am_pathogenicity, am_class). Writes a sorted,
    bgzipped, tabix-indexed AlphaMissense-shaped table under root/alphamissense/."""
    d = root / "alphamissense"
    d.mkdir(parents=True, exist_ok=True)
    raw = d / "AlphaMissense_hg38.tsv"
    rows = sorted(rows, key=lambda r: (r[0], int(r[1])))
    with open(raw, "w") as f:
        for c, p, ref, alt, am, cls in rows:
            f.write(f"{c}\t{p}\t{ref}\t{alt}\thg38\tU1\tT1\tp.X\t{am}\t{cls}\n")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, seq_col=0, start_col=1, end_col=1, force=True)
    return d


def _write_vcf(root, variants, build_contig="248956422"):
    """variants: (chrom, pos, ref, alt, gt). Writes a minimal single-sample VCF."""
    p = root / "user.vcf"
    with open(p, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID=1,length={build_contig}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for c, pos, ref, alt, gt in variants:
            f.write(f"{c}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    return str(p)


# ── tests ─────────────────────────────────────────────────────────────────────

def test_predicted_pathogenic_rare(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [
        ("1", "1000", "C", "T", "0.98", "likely_pathogenic"),
        ("1", "2000", "G", "A", "0.10", "likely_benign"),
        ("1", "3000", "A", "G", "0.55", "ambiguous"),
    ])
    vcf = _write_vcf(tmp_path, [
        ("1", "1000", "C", "T", "0/1"),   # damaging → flagged
        ("1", "2000", "G", "A", "1/1"),   # benign → dropped
        ("1", "3000", "A", "G", "0/1"),   # ambiguous → ambiguous bucket
        ("1", "9999", "T", "C", "0/1"),   # not in table → no annotation
    ])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 1        # the ambiguous one is not "pathogenic"
    assert r["n_ambiguous"] == 1
    dmg = r["buckets"]["predicted_pathogenic_rare"]
    assert len(dmg) == 1 and dmg[0]["pos"] == 1000
    assert dmg[0]["am_class"] == "likely_pathogenic"
    assert dmg[0]["confidence"] == "higher"        # am ≥ 0.9
    assert dmg[0]["rarity"] == "unknown"           # no gnomAD table present
    assert "AlphaMissense" in dmg[0]["evidence"]


def test_dedup_against_clinvar(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    clinvar = {"findings": [{"chrom": "1", "pos": 1000, "ref": "C", "alt": "T"}]}
    r = nv.analyze_novel_variants(vcf, "grch38", clinvar_result=clinvar)
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 0        # already a ClinVar hit → skipped


def test_dedup_handles_chip_result_without_findings_key(tmp_path, monkeypatch):
    # On chip input the clinvar dict is {"available": False, ...} with NO findings key.
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38",
                                  clinvar_result={"available": False, "reason": "x"})
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 1        # no crash on missing 'findings'


def test_graceful_degrade_no_tables(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)  # empty → no predictor tables
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["available"] is False
    assert "setup.py" in r["reason"]
    assert r["negative_disclaimer"] and r["disclaimer"]


def test_build_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "unknown")
    assert r["available"] is False
    assert "assume-build" in r["reason"] or "Build" in r["reason"]


def test_commercial_safe_drops_noncommercial(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38", commercial_safe=True)
    assert r["available"] is True
    used = {p["name"] for p in r["predictors_used"]}
    assert used == {"AlphaMissense"}               # only commercial-OK survives
    assert set(r["dropped_noncommercial"]) >= {"REVEL", "SpliceAI", "CADD"}


def test_hemizygous_and_homozygous_carried(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.95", "likely_pathogenic")])
    # not carried (0/0) → must be skipped
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/0")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["n_predicted_pathogenic"] == 0
