"""Scan the committed release for accidental disclosure.

Two lessons from this project's own history are encoded here. First, a scanner
that matches nothing reports a clean result forever, so the controls below run
before any assertion and the suite fails if they do not behave. Second,
substring matching produces confident nonsense: an earlier scan flagged
"position" inside "disposition", and a case-insensitive search for CAC matched
inside "efficacy". Patterns are word-bounded and adjudicated in context.
"""
import re
from pathlib import Path

import pytest
from conftest import BENCH, ROOT, SAMPLES

TEXT_SUFFIXES = {".md", ".json", ".html", ".yml", ".yaml", ".py", ".txt", ".svg"}

PRIVATE_LOCATIONS = [r"/Users/", r"/home/", r"[A-Za-z]:\\Users\\"]
CREDENTIALS = [r"AKIA[0-9A-Z]{16}", r"gh[pousr]_[A-Za-z0-9]{20,}",
               r"-----BEGIN [A-Z ]*PRIVATE KEY", r"[Bb]earer\s+[A-Za-z0-9._\-]{24,}"]
COMMERCIAL = [r"gross\s+margin", r"per[- ]patient[ /]per[- ]month", r"\bPPPM\b",
              r"\bLTV\b", r"\b(TAM|SAM|SOM)\b", r"customer acquisition cost"]
WITHHELD_FIGURES = [r"\$36\.42", r"\b39\.8\s?%", r"\$56\.81", r"\b0\.774\s?%"]


# This module necessarily contains every pattern it searches for, so scanning
# itself produces guaranteed self-matches. It is excluded by path, and the
# exclusion is exactly one file so it cannot quietly widen into a way to hide a
# real finding. (The same shape appears in docs/PUBLIC-REPORT-SPEC.md, which
# names "payback" only to prohibit it -- adjudicated below rather than excluded.)
# Files that DEFINE disclosure patterns necessarily contain them, so scanning
# them yields guaranteed self-matches. Exactly two are exempt: this module and
# the release verifier. A private path written here is a search pattern; the
# same string in README.md or a published artifact would be a real leak, and
# those are still scanned.
#
# The exemption is pinned by a test below so it cannot quietly widen into a
# place to hide a finding.
PATTERN_DEFINING_FILES = {
    (ROOT / "tests" / "public" / "test_release_security.py").resolve(),
    (ROOT / "tools" / "verify_public_release.py").resolve(),
}


def public_text_files():
    for p in sorted(ROOT.rglob("*")):
        if not (p.is_file() and p.suffix in TEXT_SUFFIXES):
            continue
        if ".git" in p.parts or p.resolve() in PATTERN_DEFINING_FILES:
            continue
        yield p


def test_the_scan_exemption_covers_only_the_pattern_definitions():
    """Guards the exemption. If it widens, the suite says so."""
    scanned = {p.resolve() for p in public_text_files()}
    everything = {p.resolve() for p in ROOT.rglob("*")
                  if p.is_file() and p.suffix in TEXT_SUFFIXES and ".git" not in p.parts}
    skipped = everything - scanned
    assert skipped == PATTERN_DEFINING_FILES, (
        f"scan exemption has changed: {sorted(str(x) for x in skipped)}")
    assert len(PATTERN_DEFINING_FILES) == 2, "exemption set grew"


def test_the_scanner_can_actually_fire():
    """Positive and negative controls. Without these, every zero below is
    unfalsifiable."""
    readme = (ROOT / "README.md").read_text()
    assert re.search(r"GenomeLens", readme), "positive control failed"
    assert not re.search(r"zzqxwv-nonexistent-token", readme), "negative control failed"


def test_word_boundaries_are_what_prevent_the_false_positive():
    """The real defect was an unbounded substring search, not case-insensitivity.

    CAC is coronary artery calcium in clinical text and customer acquisition
    cost in commercial text. An early scan in this project searched for the bare
    substring case-insensitively and matched inside "effiCACy", reporting a
    commercial hit in prevention advice. Word boundaries fix it; dropping them
    reintroduces the bug even with case preserved.
    """
    probe = "efficacy of the CAC score"
    assert len(re.findall(r"\bCAC\b", probe)) == 1
    assert len(re.findall(r"\bCAC\b", probe, re.I)) == 1     # boundaries hold under -i
    assert len(re.findall(r"CAC", probe, re.I)) == 2           # unbounded: the old bug


@pytest.mark.parametrize("pattern", PRIVATE_LOCATIONS)
def test_no_local_filesystem_paths(pattern):
    hits = [str(p.relative_to(ROOT)) for p in public_text_files()
            if re.search(pattern, p.read_text(errors="replace"))]
    assert not hits, f"{pattern} appears in {hits}"


@pytest.mark.parametrize("pattern", CREDENTIALS)
def test_no_credentials(pattern):
    hits = [str(p.relative_to(ROOT)) for p in public_text_files()
            if re.search(pattern, p.read_text(errors="replace"))]
    assert not hits, f"credential-shaped string in {hits}"


@pytest.mark.parametrize("pattern", COMMERCIAL)
def test_no_commercial_terms(pattern):
    hits = [str(p.relative_to(ROOT)) for p in public_text_files()
            if re.search(pattern, p.read_text(errors="replace"), re.I)]
    assert not hits, f"commercial language {pattern} in {hits}"


@pytest.mark.parametrize("pattern", WITHHELD_FIGURES)
def test_withheld_figures_are_absent(pattern):
    """Numerical value-of-information examples are withheld from this release
    pending review of a treatment-effect parameter."""
    hits = [str(p.relative_to(ROOT)) for p in public_text_files()
            if re.search(pattern, p.read_text(errors="replace"))]
    assert not hits, f"withheld figure {pattern} published in {hits}"


# CAC is coronary artery calcium in clinical text and customer acquisition cost
# in commercial text. Banning the token would strip legitimate prevention advice,
# so each occurrence is adjudicated by context instead. These are the contexts
# in which the token is legitimate: clinical use, a prohibition list naming it,
# or documentation explaining the false-positive itself.
CAC_LEGITIMATE_CONTEXT = ("scan", "calcium", "prohibited", "efficacy",
                          "false positive", "unbounded", "adjudicat")


def test_cac_occurrences_are_clinical_not_commercial():
    """Adjudicated rather than banned: the term is legitimate clinical
    vocabulary and must not be stripped from prevention advice."""
    for p in public_text_files():
        text = p.read_text(errors="replace")
        for m in re.finditer(r"\bCAC\b", text):
            window = text[max(0, m.start() - 160): m.end() + 160].lower()
            assert any(w in window for w in CAC_LEGITIMATE_CONTEXT), (
                f"unadjudicated CAC usage in {p.relative_to(ROOT)}: "
                f"...{text[max(0, m.start()-60):m.end()+60]}...")


def test_payback_occurrences_are_health_economic_or_declarative():
    """'payback period' is cost-effectiveness vocabulary. The only other
    allowed use is naming it in a prohibition list."""
    for p in public_text_files():
        text = p.read_text(errors="replace")
        for m in re.finditer(r"\bpayback\b", text, re.I):
            window = text[max(0, m.start() - 160): m.end() + 160].lower()
            allowed = ("prohibited", "period", "margin", "adjudicated", "vocabulary")
            assert any(w in window for w in allowed), (
                f"unadjudicated payback usage in {p.relative_to(ROOT)}")


def test_pdf_text_is_scanned_not_just_its_bytes():
    """Raw-byte grep on a compressed PDF stream proves nothing. Text is
    extracted and scanned, or this test fails loudly rather than silently
    passing."""
    fitz = pytest.importorskip("fitz", reason="PyMuPDF not installed")
    for sample in SAMPLES:
        doc = fitz.open(BENCH / sample / "report-public.pdf")
        text = "\n".join(page.get_text() for page in doc)
        doc.close()
        assert len(text) > 500, f"{sample} PDF yielded no extractable text to scan"
        for pattern in PRIVATE_LOCATIONS + COMMERCIAL + WITHHELD_FIGURES:
            assert not re.search(pattern, text, re.I), f"{pattern} in {sample} PDF"
