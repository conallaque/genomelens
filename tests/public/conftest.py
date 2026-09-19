"""Shared loaders for the public verification suite.

Everything here reads only files committed to this repository. There is no
dependency on the GenomeLens engine, no network access and no private fixture.
A fresh clone can run the whole suite.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "artifacts" / "public-benchmarks"
SAMPLES = ("HG002", "HG003", "HG004")

# The published headline claims. Every number here also appears in README.md,
# in each summary.json and in each report-public.html; the cross-artifact tests
# exist to prove those three never drift apart.
PUBLISHED = {
    "HG002": {"truth_calls": 169, "nmb": 957, "canonical": 1239},
    "HG003": {"truth_calls": 187, "nmb": 2295, "canonical": 4127},
    "HG004": {"truth_calls": 170, "nmb": 937, "canonical": 1304},
}
AGGREGATE = 526


def summary(sample: str) -> dict:
    return json.loads((BENCH / sample / "summary.json").read_text())


def report_html(sample: str) -> str:
    return (BENCH / sample / "report-public.html").read_text()


def readme() -> str:
    return (ROOT / "README.md").read_text()


@pytest.fixture(params=SAMPLES)
def sample(request) -> str:
    return request.param


@pytest.fixture
def s(sample) -> dict:
    return summary(sample)
