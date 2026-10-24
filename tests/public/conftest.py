"""Shared loaders for the public verification suite.

Everything here reads only files committed to this repository. There is no
dependency on the GenomeLens engine, no network access and no private fixture.
A fresh clone can run the whole suite.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# The historical GIAB concordance record (an earlier engine build): explicit truth calls per sample and their aggregate. Every number here also
# appears in README.md and in docs/VALIDATION.md; the tests in test_historical_concordance.py exist to keep those places from drifting apart.
PUBLISHED = {"HG002": {"truth_calls": 169}, "HG003": {"truth_calls": 187}, "HG004": {"truth_calls": 170}}
AGGREGATE = 526
SAMPLES = tuple(PUBLISHED)


def readme() -> str:
    return (ROOT / "README.md").read_text(encoding="utf-8")


def validation() -> str:
    return (ROOT / "docs" / "VALIDATION.md").read_text(encoding="utf-8")
