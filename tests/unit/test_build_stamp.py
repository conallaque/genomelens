"""Every generated artefact must say which code produced it.

THE INCIDENT THIS EXISTS FOR. A PDF was reviewed against figures that four
commits had changed, and the reviewer found the old values. The code was
correct, the tests passed, and the delivered file was in fact current — but
there was no way to establish that by looking at the file, because artefacts
carried no link to a commit. A stale copy and a current one were
indistinguishable, so the disagreement could not be settled by inspection.

These tests make the link exist and keep it there.
"""
from __future__ import annotations

import subprocess

import pytest

import build_stamp as bs


def _head() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True, cwd=bs._REPO)
    return out.stdout.strip()


def test_the_stamp_names_the_current_commit():
    s = bs.build_stamp()
    head = _head()
    if not head:                      # not a checkout; nothing to compare
        pytest.skip("not a git checkout")
    assert s["commit"] == head, (
        f"stamp says {s['commit']} but HEAD is {head}")


def test_the_marker_is_greppable_and_round_trips():
    # It has to survive being pulled out of HTML, out of extracted PDF text, or
    # read from PDF metadata — so it is one line with a fixed prefix.
    m = bs.build_stamp()["marker"]
    assert m.startswith(bs.MARKER_PREFIX)
    assert bs.commit_of(f"noise before\n{m}\nnoise after") == bs.build_stamp()["commit"]
    assert bs.marker_in("nothing here") == ""
    assert bs.commit_of("nothing here") == ""


def test_a_dirty_tree_is_declared_not_hidden():
    # An artefact built with uncommitted edits reflects code in no commit. That
    # is the failure mode that made the original mismatch unresolvable, so it
    # must be visible in the artefact rather than inferred.
    dirty = subprocess.run(["git", "status", "--porcelain"],
                           capture_output=True, text=True, cwd=bs._REPO).stdout.strip()
    assert bs.build_stamp()["dirty"] is bool(dirty)
    assert bs.is_dirty_build("GENOMELENS-BUILD: abc1234+dirty 2026-01-01T00:00:00")
    assert not bs.is_dirty_build("GENOMELENS-BUILD: abc1234 2026-01-01T00:00:00")


def test_the_stamp_degrades_instead_of_raising():
    # A report is worth more than a provenance line; losing the line must not
    # lose the report.
    import build_stamp
    orig = build_stamp._REPO
    try:
        build_stamp._REPO = "/nonexistent-path-for-this-test"
        s = build_stamp.build_stamp()
        assert s["commit"] == "unknown"
        assert s["marker"].startswith(build_stamp.MARKER_PREFIX)
    finally:
        build_stamp._REPO = orig


def test_a_rendered_report_carries_the_marker():
    # THE ACTUAL GUARD. Render the report header and confirm the marker is in
    # the output — not just that build_stamp() works in isolation.
    import renderers
    html = renderers.build_html_report(
        tier1_results=[], apoe_genotype=None, ai_results={},
        exec_summary=None, dna_filepath="test_genome.txt",
        no_ai=True, model="")
    assert bs.MARKER_PREFIX in html, "the report shipped with no build provenance"
    assert bs.commit_of(html) == bs.build_stamp()["commit"]
