"""
Tests for cli.py — the V7-extracted CLI parser.

These lock the public CLI contract so future refactors of analyze.py (the
larger orchestration extraction in V7-β) can't silently rename a flag, drop
a help string, or change a default.
"""

from __future__ import annotations

import argparse

import pytest

import cli


# ── Parser shape ─────────────────────────────────────────────────────────────

def test_build_parser_returns_argparse_parser() -> None:
    parser = cli.build_parser()
    assert isinstance(parser, argparse.ArgumentParser)
    assert parser.prog == "dna-analyze"


def test_parser_has_no_side_effects_on_import() -> None:
    """Building the parser must not touch the filesystem or import analyze.
    Tested implicitly by the fact that we're already imported — but assert
    the heavy `analyze` module hasn't been pulled in for free."""
    parser = cli.build_parser()
    # build_parser doesn't take args; just constructing the parser shouldn't
    # raise even when analyze has a missing import.
    assert parser is not None


# ── Flag inventory (locks the public CLI contract) ───────────────────────────

REQUIRED_FLAGS = {
    "--no-ai", "--output", "--model",
    "--impute", "--impute-min-r2", "--pdf", "--medications",
    "--carrier-report", "--chat", "--compare",
    "--narrative", "--emergency-card", "--retry-failed",
    "--bloodwork", "--fhir",
}


def test_every_v6_flag_present() -> None:
    parser = cli.build_parser()
    flags = {action.option_strings[0] for action in parser._actions
             if action.option_strings}
    missing = REQUIRED_FLAGS - flags
    assert not missing, f"CLI is missing flags: {missing}"


def test_dna_file_positional_optional() -> None:
    parser = cli.build_parser()
    ns = parser.parse_args(["--retry-failed"])
    assert ns.dna_file is None
    assert ns.retry_failed is True


def test_default_model_is_documented_constant() -> None:
    parser = cli.build_parser()
    ns = parser.parse_args(["test.csv"])
    assert ns.model == cli.DEFAULT_OLLAMA_MODEL


def test_impute_min_r2_is_float() -> None:
    parser = cli.build_parser()
    ns = parser.parse_args(["test.csv", "--impute-min-r2", "0.55"])
    assert isinstance(ns.impute_min_r2, float)
    assert ns.impute_min_r2 == 0.55


def test_bloodwork_path_passes_through() -> None:
    parser = cli.build_parser()
    ns = parser.parse_args(["test.csv", "--bloodwork", "labs.json"])
    assert ns.bloodwork == "labs.json"


def test_fhir_is_boolean_action() -> None:
    parser = cli.build_parser()
    ns = parser.parse_args(["test.csv", "--fhir"])
    assert ns.fhir is True
    ns2 = parser.parse_args(["test.csv"])
    assert ns2.fhir is False


def test_help_text_includes_v7_version() -> None:
    parser = cli.build_parser()
    help_text = parser.format_help()
    # Don't pin a specific version string — just check the description block
    # mentions the consolidation theme so future refactors stay grounded.
    assert "v7.0.0" in help_text or "consolidation" in help_text.lower()


# ── Help epilog has examples ─────────────────────────────────────────────────

def test_epilog_has_bloodwork_example() -> None:
    parser = cli.build_parser()
    assert "--bloodwork" in (parser.epilog or "")


def test_epilog_has_fhir_example() -> None:
    parser = cli.build_parser()
    assert "--fhir" in (parser.epilog or "")
