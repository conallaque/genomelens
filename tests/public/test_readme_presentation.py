"""How the three demonstrations are presented to a first-time visitor.

GitHub's file viewer shows a committed .html file as source code, so linking a
report's .html path from the README sends a reader to markup rather than to a
report. The primary call to action is therefore the PDF, which GitHub renders in
its own viewer, with a GitHub Pages URL offering the same report as a web page.

These tests check paths and wiring only. They deliberately do not fetch the
Pages URLs: a suite that fails because an external host is slow is a suite people
learn to ignore.
"""
import re
from conftest import BENCH, ROOT, SAMPLES

README = (ROOT / "README.md").read_text()
PAGES = "https://conallaque.github.io/genomelens/artifacts/public-benchmarks"
REL = "artifacts/public-benchmarks"


def test_preview_images_are_embedded_not_merely_linked(sample):
    """A linked image is a file listing. An embedded one is the proof."""
    assert f'<img src="{REL}/{sample}/preview.png"' in README, (
        f"{sample} preview is not embedded as an image in README")
    assert (BENCH / sample / "preview.png").is_file()


def test_the_primary_report_link_is_the_pdf(sample):
    assert f"[Report (PDF)]({REL}/{sample}/report-public.pdf)" in README
    assert (BENCH / sample / "report-public.pdf").is_file()


def test_no_readme_link_sends_a_reader_to_html_source(sample):
    """The failure this file exists to prevent: a relative .html link renders
    as source code in GitHub's file viewer."""
    bad = f"]({REL}/{sample}/report-public.html)"
    assert bad not in README, (
        f"{sample}: README links the .html path directly, which GitHub shows as "
        "source code. Link the PDF or the Pages URL instead.")


def test_a_rendered_web_report_is_offered(sample):
    assert f"{PAGES}/{sample}/report-public.html" in README


def test_sample_links_never_cross_reference_another_sample(sample):
    """A copy-paste error here would show a visitor the wrong genome."""
    others = [x for x in SAMPLES if x != sample]
    for block in re.findall(rf"\[[^\]]*\]\([^)]*{sample}[^)]*\)", README):
        for other in others:
            assert other not in block, f"{sample} link references {other}: {block}"


def test_pages_serves_files_verbatim():
    """Without .nojekyll, GitHub Pages runs the site through Jekyll, which can
    transform or omit files. The reports must be served exactly as committed."""
    assert (ROOT / ".nojekyll").is_file(), (
        ".nojekyll missing; Pages would process the reports rather than serving them")


def test_every_referenced_local_artifact_exists():
    missing = []
    for target in re.findall(r"\]\((artifacts/[^)]+)\)", README):
        if not (ROOT / target).exists():
            missing.append(target)
    assert not missing, f"README references missing artifacts: {missing}"
