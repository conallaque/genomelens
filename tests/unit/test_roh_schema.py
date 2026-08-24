"""Regression: the ROH empty/degenerate result must be schema-consistent with the
full result, or the pipeline log line and build_roh_html crash the whole report.

This locks in the fix for a real crash: a near-empty genome made detect_roh take an
early return that omitted context_tier / n_short / n_medium / n_long, which then
raised KeyError in pipeline.py and renderers.build_roh_html.
"""
from __future__ import annotations

import pandas as pd
import pytest

from risk import roh

_FULL_KEYS = {"runs", "n_runs", "total_roh_mb", "f_roh", "short", "medium", "long",
              "n_short", "n_medium", "n_long", "population_context", "context_tier"}


def test_empty_input_result_has_full_schema():
    r = roh.detect_roh(pd.DataFrame())
    assert set(r.keys()) >= _FULL_KEYS, f"missing: {_FULL_KEYS - set(r.keys())}"
    # the keys the downstream log line and renderer index directly
    assert r["context_tier"] == "unavailable"
    assert r["n_short"] == 0 and r["n_medium"] == 0 and r["n_long"] == 0


def test_no_chrom_column_result_has_full_schema():
    r = roh.detect_roh(pd.DataFrame({"genotype": ["AA", "GG"]}))
    assert set(r.keys()) >= _FULL_KEYS


def test_build_roh_html_survives_empty_result():
    renderers = pytest.importorskip("report.renderers")
    # must not raise, even on the degenerate result
    html = renderers.build_roh_html(roh.detect_roh(pd.DataFrame()))
    assert isinstance(html, str)


def test_build_roh_html_skips_partial_result_without_crashing():
    renderers = pytest.importorskip("report.renderers")
    # A malformed partial dict (missing core stats) must degrade, not crash.
    html = renderers.build_roh_html({"runs": [], "f_roh": 0.0})
    assert html == ""
