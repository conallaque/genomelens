"""Unit tests for genome_input.py — dual chip/VCF input, genotype conversion,
build-gated back-fill, and VCF profiling."""

from __future__ import annotations

import os
import tempfile

import pandas as pd
import pytest

import genome_input as gi
import snp_registry as reg

# ── genotype conversion rules (one case per rule) ─────────────────────────────

@pytest.mark.parametrize("gt,ref,alt,expected", [
    ("0/1", "A", "G", "AG"),        # het SNV, sorted
    ("1|1", "A", "G", "GG"),        # phased hom ALT
    ("0/0", "A", "G", "AA"),        # hom REF
    ("1", "A", "G", "GG"),          # haploid → homozygous
    ("1/2", "A", "G,T", "GT"),      # multiallelic, both ALTs
    ("0/1:35:99", "A", "G", "AG"),  # extra FORMAT sub-fields ignored
    ("./.", "A", "G", None),        # no-call
    ("0/.", "A", "G", None),        # partial no-call
    (".", "A", "G", None),          # single no-call
    ("0/1", "AT", "A", None),       # deletion (indel) → skip
    ("0/1", "A", "ATG", None),      # insertion → skip
    ("0/1", "A", "<DEL>", None),    # symbolic ALT → skip
    ("3/3", "A", "G", None),        # allele index out of range → skip
])
def test_gt_conversion(gt, ref, alt, expected):
    assert gi.vcf_gt_to_genotype(gt, ref, alt) == expected


# ── input-type detection ──────────────────────────────────────────────────────

def test_looks_like_vcf_by_extension():
    assert gi.looks_like_vcf("/x/genome.vcf")
    assert gi.looks_like_vcf("/x/genome.vcf.gz")
    assert not gi.looks_like_vcf("/x/genome.txt")


def test_looks_like_vcf_by_header():
    f = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False)
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n"); f.close()
    try:
        assert gi.looks_like_vcf(f.name)
    finally:
        os.unlink(f.name)


# ── build-gated back-fill + profiling ─────────────────────────────────────────

def _write_vcf(records_with_rsid, records_position_only, extra_lines=()):
    lines = ["##fileformat=VCFv4.2",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    for r, gt in records_with_rsid:
        lines.append(f"{r.chrom}\t{r.pos_grch37}\t{r.rsid}\t{r.ancestral}\t{r.derived}\t.\t.\t.\tGT\t{gt}")
    for r, gt in records_position_only:
        lines.append(f"{r.chrom}\t{r.pos_grch37}\t.\t{r.ancestral}\t{r.derived}\t.\t.\t.\tGT\t{gt}")
    lines.extend(extra_lines)
    f = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    f.write("\n".join(lines) + "\n"); f.close()
    return f.name


def test_backfill_recovers_position_only_registry_variant():
    r0, r1 = reg._RECORDS[0], reg._RECORDS[1]
    vcf = _write_vcf([(r0, "0/1")], [(r1, "1/1")])
    base = pd.DataFrame(
        {"chrom": [str(r0.chrom)], "pos": [r0.pos_grch37],
         "genotype": ["".join(sorted(r0.ancestral + r0.derived))], "source": ["chip"]},
        index=[r0.rsid])
    base.index.name = "rsid"
    try:
        enr, prof = gi.enrich_and_profile_vcf(base, vcf, "grch37")
        assert r1.rsid in enr.index                      # position-only variant recovered
        assert enr.loc[r1.rsid, "genotype"] == "".join(sorted(r1.derived + r1.derived))
        assert prof["backfilled"] == 1
        assert prof["without_rsid"] >= 1
    finally:
        os.unlink(vcf)


def test_build_gate_refuses_on_unknown_build():
    r1 = reg._RECORDS[1]
    vcf = _write_vcf([], [(r1, "1/1")])
    base = pd.DataFrame(columns=["chrom", "pos", "genotype", "source"])
    base.index.name = "rsid"
    try:
        enr, prof = gi.enrich_and_profile_vcf(base, vcf, "unknown")
        assert prof["backfill_gated"] is True
        assert prof["backfilled"] == 0          # never guesses coordinates
        assert r1.rsid not in enr.index
    finally:
        os.unlink(vcf)


def test_profile_counts_total_and_acmg():
    vcf = _write_vcf([], [], extra_lines=[
        "17\t41250000\t.\tC\tT\t.\t.\t.\tGT\t0/1",   # BRCA1 window (ACMG)
        "7\t123\t.\tA\tG\t.\t.\t.\tGT\t0/1",         # non-ACMG, no rsID
        "1\t456\trs999999\tA\tG\t.\t.\t.\tGT\t0/1",  # has rsID
    ])
    base = pd.DataFrame(columns=["chrom", "pos", "genotype", "source"]); base.index.name = "rsid"
    try:
        _, prof = gi.enrich_and_profile_vcf(base, vcf, "grch37")
        assert prof["total_variants"] == 3
        assert prof["without_rsid"] == 2
        assert prof["acmg_gene_variants"] == 1
    finally:
        os.unlink(vcf)
