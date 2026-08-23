"""Tests for provenance.py — V8 variant-source tagging + build detection."""

from __future__ import annotations

import pandas as pd

import provenance
import snp_registry as reg

# ── tag_chip_source ──────────────────────────────────────────────────────────

def test_tag_chip_source_adds_column() -> None:
    df = pd.DataFrame({"genotype": ["AA", "CT"]}, index=["rs1", "rs2"])
    tagged = provenance.tag_chip_source(df)
    assert "source" in tagged.columns
    assert (tagged["source"] == provenance.CHIP).all()


def test_tag_chip_source_is_idempotent() -> None:
    df = pd.DataFrame({"genotype": ["AA"], "source": [provenance.CHIP]},
                      index=["rs1"])
    out = provenance.tag_chip_source(df)
    assert (out["source"] == provenance.CHIP).all()


def test_tag_chip_source_does_not_mutate_caller() -> None:
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs1"])
    _ = provenance.tag_chip_source(df)
    assert "source" not in df.columns


# ── tag_imputed_rows ─────────────────────────────────────────────────────────

def test_tag_imputed_with_high_r2_marks_high() -> None:
    df = pd.DataFrame({"genotype": ["AA", "CT", "GG"]},
                      index=["rs_chip", "rs_imp_good", "rs_imp_bad"])
    df = provenance.tag_chip_source(df)
    dr2 = {"rs_imp_good": 0.95, "rs_imp_bad": 0.4}
    out = provenance.tag_imputed_rows(
        df, chip_rsids={"rs_chip"}, dr2_by_rsid=dr2,
    )
    assert out.at["rs_chip", "source"] == provenance.CHIP
    assert out.at["rs_imp_good", "source"] == provenance.IMP_HIGH_R2
    assert out.at["rs_imp_bad", "source"] == provenance.IMP_LOW_R2


def test_tag_imputed_without_dr2_defaults_to_low() -> None:
    """When DR2 isn't supplied (e.g. legacy imputation runs), every imputed
    row is conservatively tagged as low-confidence."""
    df = pd.DataFrame({"genotype": ["AA", "CT"]},
                      index=["rs_chip", "rs_imp"])
    df = provenance.tag_chip_source(df)
    out = provenance.tag_imputed_rows(df, chip_rsids={"rs_chip"})
    assert out.at["rs_imp", "source"] == provenance.IMP_LOW_R2


def test_provenance_summary_counts() -> None:
    df = pd.DataFrame(
        {"genotype": ["A"] * 4,
         "source": [provenance.CHIP, provenance.CHIP,
                    provenance.IMP_HIGH_R2, provenance.IMP_LOW_R2]},
        index=["a", "b", "c", "d"],
    )
    s = provenance.provenance_summary(df)
    assert s["chip"] == 2
    assert s["imp_high_r2"] == 1
    assert s["imp_low_r2"] == 1
    assert s["total"] == 4
    assert s["tagged"] is True


def test_provenance_summary_on_untagged_df() -> None:
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs1"])
    s = provenance.provenance_summary(df)
    assert s["tagged"] is False
    assert s["chip"] == 1   # graceful degradation — assume chip-only
    assert s["imp_high_r2"] == 0


# ── Build auto-detection ─────────────────────────────────────────────────────

def test_detect_build_grch37_via_synthetic_chip(synthetic_snps_df) -> None:
    """The conftest fixture genome uses GRCh37 positions for every SNP."""
    # Ensure the fixture has a `pos` column (it does — see conftest)
    assert "pos" in synthetic_snps_df.columns
    assert provenance.detect_build(synthetic_snps_df) == "grch37"


def test_detect_build_grch38_when_positions_match() -> None:
    """Build a probe DataFrame with GRCh38 positions for known registry SNPs."""
    probe_rsids = provenance._build_probe_set()
    assert len(probe_rsids) >= 3, "registry must have ≥3 dual-build probes"
    rows = []
    for rsid in probe_rsids[:6]:
        rec = reg.SNPS[rsid]
        rows.append({"rsid": rsid, "pos": rec.pos_grch38, "genotype": "AA"})
    df = pd.DataFrame(rows).set_index("rsid")
    assert provenance.detect_build(df) == "grch38"


def test_detect_build_mixed_when_some_probes_drift() -> None:
    """A chip that reports half its probes on GRCh37, half on GRCh38 is
    'mixed' — common with re-mapped Illumina arrays."""
    probe_rsids = provenance._build_probe_set()
    assert len(probe_rsids) >= 4
    rows = []
    half = len(probe_rsids) // 2
    for rsid in probe_rsids[:half]:
        rec = reg.SNPS[rsid]
        rows.append({"rsid": rsid, "pos": rec.pos_grch37, "genotype": "AA"})
    for rsid in probe_rsids[half:half * 2]:
        rec = reg.SNPS[rsid]
        rows.append({"rsid": rsid, "pos": rec.pos_grch38, "genotype": "AA"})
    df = pd.DataFrame(rows).set_index("rsid")
    assert provenance.detect_build(df) == "mixed"


def test_detect_build_unknown_without_pos_column() -> None:
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs1"])
    assert provenance.detect_build(df) == "unknown"


def test_detect_build_unknown_with_too_few_probes() -> None:
    df = pd.DataFrame({"pos": [1234], "genotype": ["AA"]}, index=["rs999999"])
    assert provenance.detect_build(df) == "unknown"


def test_annotate_build_sets_attrs(synthetic_snps_df) -> None:
    df = provenance.annotate_build(synthetic_snps_df.copy())
    assert df.attrs.get("build") == "grch37"


def test_annotate_parsed_does_both(synthetic_snps_df) -> None:
    df = provenance.annotate_parsed(synthetic_snps_df.copy())
    assert "source" in df.columns
    assert (df["source"] == provenance.CHIP).all()
    assert df.attrs.get("build") == "grch37"
