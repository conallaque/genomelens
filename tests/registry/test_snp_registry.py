"""
Tests for the unified SNP registry.

The registry is the data foundation V7 builds on — these tests must stay
green or the migration of consumer modules will silently propagate bad data.
"""

from __future__ import annotations

import pandas as pd
import pytest

from core import snp_registry as reg

# ── Record-level invariants ──────────────────────────────────────────────────

def test_rsid_must_start_with_rs() -> None:
    with pytest.raises(ValueError, match="must start with"):
        reg.SNPRecord(
            rsid="rs1801133bad".replace("rs", "x"),
            gene="MTHFR", chrom="1",
            pos_grch37=1, pos_grch38=1,
            ancestral="C", derived="T", description="x",
        )


def test_ancestral_must_be_actg_on_plus_strand() -> None:
    with pytest.raises(ValueError, match="ancestral"):
        reg.SNPRecord(
            rsid="rs1", gene="X", chrom="1",
            pos_grch37=1, pos_grch38=1,
            ancestral="N", derived="T", description="x",
        )


def test_derived_must_differ_from_ancestral() -> None:
    with pytest.raises(ValueError, match="ancestral == derived"):
        reg.SNPRecord(
            rsid="rs1", gene="X", chrom="1",
            pos_grch37=1, pos_grch38=1,
            ancestral="C", derived="C", description="x",
        )


def test_record_is_frozen() -> None:
    """frozen=True → modules can't accidentally mutate canonical metadata."""
    rec = reg.get("rs1801133")
    assert rec is not None
    with pytest.raises(AttributeError):
        rec.gene = "ELSEWHERE"  # type: ignore[misc]


def test_record_is_hashable() -> None:
    """Important for using records in sets / as dict keys."""
    a = reg.get("rs1801133")
    b = reg.get("rs1801133")
    assert hash(a) == hash(b)


# ── Strand-aware dose helper ─────────────────────────────────────────────────

@pytest.mark.parametrize("genotype,risk,ref,expected", [
    # MTHFR C677T — risk T (+ strand) / A (− strand); ref C / G.
    ("TT", "T", "C", 2),
    ("AA", "T", "C", 2),     # − strand homozygous risk
    ("CT", "T", "C", 1),
    ("GA", "T", "C", 1),     # − strand het
    ("CC", "T", "C", 0),
    ("GG", "T", "C", 0),
    # Edge cases
    ("--", "T", "C", None),
    ("00", "T", "C", None),
    ("", "T", "C", None),
])
def test_risk_dose_canonical_cases(genotype, risk, ref, expected) -> None:
    assert reg.risk_dose(genotype, risk_allele=risk, ref_allele=ref) == expected


def test_risk_dose_from_df_uses_registry_alleles_when_omitted() -> None:
    """If the caller doesn't supply risk/ref alleles, take them from SNPS."""
    df = pd.DataFrame({"genotype": ["CT"]}, index=["rs1801133"])
    # MTHFR C677T canonical: ancestral=C, derived=T → 1 risk allele
    assert reg.risk_dose_from_df(df, "rs1801133") == 1


def test_risk_dose_from_df_unknown_rsid_returns_none() -> None:
    df = pd.DataFrame({"genotype": ["AA"]}, index=["rs999999"])
    assert reg.risk_dose_from_df(df, "rs999999") is None


def test_risk_dose_from_df_handles_hemizygous_y_chrom() -> None:
    """Y-chromosome SNPs sometimes report as a single allele 'T' instead of 'TT'."""
    df = pd.DataFrame({"genotype": ["T"]}, index=["rs1801133"])  # using known rsID
    # Single-char → duplicated to 'TT' internally → risk dose 2
    assert reg.risk_dose_from_df(df, "rs1801133") == 2


# ── Lookup API ───────────────────────────────────────────────────────────────

def test_get_returns_none_for_missing() -> None:
    assert reg.get("rs99999999") is None


def test_require_raises_for_missing() -> None:
    with pytest.raises(KeyError, match="not in snp_registry"):
        reg.require("rs99999999")


def test_by_gene_groups_correctly() -> None:
    mthfr = reg.by_gene("MTHFR")
    rsids = {r.rsid for r in mthfr}
    assert {"rs1801133", "rs1801131"}.issubset(rsids)


def test_lookup_by_pos_grch37() -> None:
    rec = reg.lookup_by_pos(chrom="1", pos=11_856_378, build="grch37")
    assert rec is not None
    assert rec.rsid == "rs1801133"


def test_lookup_by_pos_grch38() -> None:
    """The same SNP lookup-able via the other build's coordinates."""
    rec = reg.lookup_by_pos(chrom="1", pos=11_796_321, build="grch38")
    assert rec is not None
    assert rec.rsid == "rs1801133"


# ── Audit ────────────────────────────────────────────────────────────────────

def test_audit_summary_shape() -> None:
    audit = reg.audit_registry()
    assert audit["n_records"] >= 19   # current seed
    assert audit["n_genes"] >= 10
    assert isinstance(audit["stale_rsids"], list)


def test_audit_flags_stale_records() -> None:
    """Run audit far in the future — all records should be stale."""
    import datetime
    audit = reg.audit_registry(today=datetime.date(2099, 1, 1))
    assert audit["n_stale"] == audit["n_records"]
