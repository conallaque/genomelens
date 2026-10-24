"""How the repository is presented to a first-time visitor.

These tests check paths and wiring only. They deliberately do not fetch any external URL: a suite that fails because an external host is slow is
a suite people learn to ignore.
"""
from conftest import ROOT, readme


def test_pages_serves_files_verbatim():
    """Without .nojekyll, GitHub Pages runs the site through Jekyll, which can transform or omit files. Files must be served exactly as committed."""
    assert (ROOT / ".nojekyll").is_file(), ".nojekyll missing; Pages would process the files rather than serving them"


def test_the_readme_links_nothing_that_was_removed():
    """The earlier trio's per-genome reports left the tree; a link to one would send a visitor to a missing page."""
    text = readme()
    for gone in ("artifacts/public-benchmarks", "report-public.pdf", "report-public.html"):
        assert gone not in text, gone
