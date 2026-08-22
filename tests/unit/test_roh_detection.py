"""Regression tests for three scientifically-wrong ROH detection behaviours.

1) HET BUDGET IS PER SLIDING WINDOW, NOT PER RUN — a whole-run counter ends a
   segment at the first het overrun, fragmenting genuine IBD tracts into several
   short runs. That both depresses F_ROH and moves length out of the long class
   that consanguinity inference depends on.

2) RUNS DO NOT BRIDGE LARGE GAPS — an uncalled stretch (centromere, coverage
   hole) carries no homozygosity evidence, so extending across it manufactures a
   single long run from two unrelated homozygous blocks.

3) CONSANGUINITY TIERS KEY ON LONG ROH ONLY — only long (>8 Mb) runs track
   recent pedigree relatedness (Kirin 2010). Blending in short/medium runs, which
   reflect ancient shared ancestry, labels founder-population background as
   cousin parentage.

These assert structural properties, not memorised outputs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import roh


def _hom_array(n):
    return np.ones(n, dtype=np.int8)


# ── 1) sliding window ────────────────────────────────────────────────────────

def test_scattered_hets_do_not_fragment_one_long_segment():
    # 800 SNPs over 20 Mb with two well-separated hets. With a whole-run budget
    # of max_hets=1 the second het would truncate the run; with a sliding window
    # the segment survives intact.
    n = 800
    pos = np.arange(n, dtype=np.int64) * 25_000
    hom = _hom_array(n)
    hom[71] = 0
    hom[619] = 0                      # >window_snps apart
    runs = roh._detect_roh_one_chrom(pos, hom, min_snps=50, min_length_mb=1.0,
                                     max_hets=1)
    assert len(runs) == 1
    assert runs[0]["length_mb"] > 15.0


def test_hets_inside_one_window_still_break_the_run():
    # The budget must still bite when hets are genuinely dense: two hets a few
    # SNPs apart exceed max_hets=1 within a single window.
    n = 800
    pos = np.arange(n, dtype=np.int64) * 25_000
    hom = _hom_array(n)
    hom[400] = 0
    hom[402] = 0
    runs = roh._detect_roh_one_chrom(pos, hom, min_snps=50, min_length_mb=1.0,
                                     max_hets=1)
    assert len(runs) >= 2, "dense hets must still terminate a run"


def test_window_tolerance_increases_with_max_hets():
    n = 400
    pos = np.arange(n, dtype=np.int64) * 25_000
    hom = _hom_array(n)
    for i in (100, 102, 104):
        hom[i] = 0
    strict = roh._detect_roh_one_chrom(pos, hom, 50, 1.0, max_hets=1)
    loose = roh._detect_roh_one_chrom(pos, hom, 50, 1.0, max_hets=5)
    assert max(r["length_mb"] for r in loose) >= max(r["length_mb"] for r in strict)


# ── 2) gap guard ─────────────────────────────────────────────────────────────

def test_large_uncalled_gap_is_not_bridged():
    # Two 5 Mb homozygous blocks separated by a ~40 Mb hole must stay separate.
    left = np.arange(60, dtype=np.int64) * 83_000
    right = 45_000_000 + np.arange(60, dtype=np.int64) * 83_000
    pos = np.concatenate([left, right])
    runs = roh._detect_roh_one_chrom(pos, _hom_array(120), 50, 1.0, 1)
    assert all(r["length_mb"] < 30 for r in runs), "gap was bridged into a false long run"


def test_small_gap_still_merges():
    # A sub-threshold gap is ordinary chip sparsity and must not split a run.
    pos = np.concatenate([np.arange(60, dtype=np.int64) * 83_000,
                          5_500_000 + np.arange(60, dtype=np.int64) * 83_000])
    runs = roh._detect_roh_one_chrom(pos, _hom_array(120), 50, 1.0, 1)
    assert len(runs) == 1


def test_gap_threshold_is_configurable_and_monotone():
    pos = np.concatenate([np.arange(60, dtype=np.int64) * 83_000,
                          8_000_000 + np.arange(60, dtype=np.int64) * 83_000])
    tight = roh._detect_roh_one_chrom(pos, _hom_array(120), 50, 1.0, 1,
                                      max_gap_mb=0.5)
    wide = roh._detect_roh_one_chrom(pos, _hom_array(120), 50, 1.0, 1,
                                     max_gap_mb=10.0)
    assert len(tight) > len(wide)


# ── 3) long-ROH-only consanguinity gating ────────────────────────────────────

def _synthetic_genome(segments, het_bg=0.35, seed=1):
    rng = np.random.default_rng(seed)
    frames = []
    for c, clen in roh.CHROM_LENGTHS_MB.items():
        n = int(clen * 40)
        pos = np.sort(rng.choice(int(clen * 1e6), size=n, replace=False))
        gts = ["AA" if (any(cc == c and s * 1e6 <= p <= (s + L) * 1e6
                            for cc, s, L in segments) or rng.random() > het_bg)
               else "AG" for p in pos]
        frames.append(pd.DataFrame({"chrom": c, "pos": pos, "genotype": gts}))
    return pd.concat(frames, ignore_index=True)


def test_founder_background_is_not_called_consanguineous():
    # Many short/medium runs, zero long runs: high TOTAL F_ROH but no evidence
    # of recent parental relatedness. This is the false-positive the old
    # blended-F_ROH tiers produced.
    segs = [(c, 10 + i * 12, 2.5)
            for c in ("1", "2", "3", "4", "5", "6", "7", "8") for i in range(6)]
    r = roh.detect_roh(_synthetic_genome(segs, seed=11))
    assert r["n_long"] == 0
    assert r["f_roh"] > 0.02, "test premise: total F_ROH is high enough to have misfired"
    assert r["f_roh_long"] == 0.0
    assert "cousin" not in r["context_tier"]


def test_genuine_long_roh_burden_is_still_detected():
    segs = [("1", 20, 40), ("2", 30, 35), ("3", 15, 30), ("4", 25, 28), ("5", 40, 25)]
    r = roh.detect_roh(_synthetic_genome(segs, seed=12))
    assert r["n_long"] >= 4
    assert "cousin" in r["context_tier"]


def test_f_roh_long_never_exceeds_total_f_roh():
    segs = [("1", 20, 30), ("2", 10, 3), ("3", 5, 1.5)]
    r = roh.detect_roh(_synthetic_genome(segs, seed=5))
    assert r["f_roh_long"] <= r["f_roh"]


def test_empty_genome_carries_the_long_roh_keys():
    # Schema parity: the degraded return must expose the same keys, or the
    # renderer and pipeline log hit a KeyError.
    r = roh.detect_roh(pd.DataFrame())
    for k in ("f_roh", "f_roh_long", "long_roh_mb", "context_tier"):
        assert k in r
