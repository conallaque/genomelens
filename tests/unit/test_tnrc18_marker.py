"""Tests for the TNRC18 rs117910193 novelty marker lookup."""
from __future__ import annotations

import pandas as pd

from risk import tnrc18_marker as tm


def _rsid(gt):
    return pd.DataFrame({"chrom": ["7"], "pos": [5401412], "genotype": [gt]},
                        index=["rs117910193"])


def _positional(gt):
    return pd.DataFrame({"chrom": ["chr7"], "pos": [5401412], "genotype": [gt]},
                        index=["7:5401412"])


def test_homozygous_gg_is_target_wild_type():
    r = tm.analyze_tnrc18_marker(_rsid("GG"))
    assert r["available"] and r["marker"] == "Target Trait Marker (Wild Type)"
    assert r["is_target_wild_type"] is True
    assert r["zygosity"] == "homozygous reference"


def test_heterozygous_ga_is_standard_marker():
    r = tm.analyze_tnrc18_marker(_rsid("GA"))
    assert r["marker"] == "Standard Marker"
    assert r["is_target_wild_type"] is False


def test_allele_order_does_not_matter():
    assert tm.analyze_tnrc18_marker(_rsid("AG"))["marker"] == "Standard Marker"


def test_homozygous_minor_is_neither():
    r = tm.analyze_tnrc18_marker(_rsid("AA"))
    assert r["zygosity"] == "homozygous minor"
    assert "Neither" in r["marker"]


def test_opposite_strand_is_normalised():
    # CC == GG (wild type), CT == GA (het), on the complementary strand.
    gg = tm.analyze_tnrc18_marker(_rsid("CC"))
    assert gg["marker"] == "Target Trait Marker (Wild Type)"
    assert gg["strand_flipped"] is True
    assert tm.analyze_tnrc18_marker(_rsid("CT"))["marker"] == "Standard Marker"


def test_slash_and_pipe_formats_parse():
    assert tm.analyze_tnrc18_marker(_rsid("G/A"))["marker"] == "Standard Marker"
    assert tm.analyze_tnrc18_marker(_rsid("G|G"))["marker"] == "Target Trait Marker (Wild Type)"


def test_positional_fallback_when_no_rsid():
    r = tm.analyze_tnrc18_marker(_positional("GA"))
    assert r["available"] and r["matched_by"] == "chrom:pos"


def test_not_typed_is_graceful():
    r = tm.analyze_tnrc18_marker(pd.DataFrame({"chrom": ["1"], "pos": [1],
                                               "genotype": ["AA"]}, index=["rs1"]))
    assert r["available"] is False and "not" in r["reason"].lower()


def test_empty_and_malformed_are_graceful():
    assert tm.analyze_tnrc18_marker(pd.DataFrame())["available"] is False
    assert tm.analyze_tnrc18_marker(_rsid("--"))["available"] is False
    assert tm.analyze_tnrc18_marker(_rsid("I"))["available"] is False   # indel/short


def test_never_raises_and_html_builds():
    for gt in ("GG", "GA", "AA", "CC", "--", "XY", "G"):
        r = tm.analyze_tnrc18_marker(_rsid(gt))
        assert isinstance(tm.build_tnrc18_html(r), str)      # renderer never crashes
    assert tm.build_tnrc18_html(None) == ""


def test_makes_no_health_claim():
    # The disclaimer must be present — this is a novelty lookup, not a trait/health call.
    r = tm.analyze_tnrc18_marker(_rsid("GG"))
    assert "no validated" in r["disclaimer"].lower()


def test_registry_entry_enables_wgs_backfill():
    # REGRESSION: on a whole-genome VCF without rsIDs, snps_df is empty at analysis
    # time, so the marker is only reachable if its coordinate is in the SNP registry
    # (which drives the VCF coordinate back-fill). Without this entry the WGS path
    # silently returned "not typed" even when the variant was in the genome.
    from core import snp_registry as reg
    rec = [r for r in reg._RECORDS if r.rsid == tm.RSID]
    assert rec, f"{tm.RSID} must be in the registry for WGS coordinate back-fill"
    assert rec[0].chrom == tm.CHROM
    assert rec[0].pos_grch38 == tm.POSITION      # the coordinate the back-fill matches
