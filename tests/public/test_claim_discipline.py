"""What the public artifacts are not allowed to claim.

The scoped concordance claim is defensible. The unscoped version of the same
sentence is not, and the distance between them is one word.
"""
import re
from conftest import ROOT, SAMPLES, report_html, summary

TEXT = {".md", ".html", ".json"}

# Each entry is (pattern, why it would be false).
OVERCLAIMS = [
    (r"100\s*%\s*(whole[- ]genome\s*)?accuracy", "concordance is scoped, not whole-genome accuracy"),
    (r"clinically\s+validated", "no clinical validation has been performed"),
    (r"GIAB[- ]certified", "GIAB certifies nothing"),
    (r"NIST[- ]certified", "NIST certifies nothing"),
    (r"FDA[- ](approved|cleared)", "no regulatory clearance is held"),
    (r"guaranteed\s+savings", "modeled economics are not realized savings"),
    (r"diagnostic\s+device", "explicitly disclaimed"),
]


def _public_text():
    for p in sorted(ROOT.rglob("*")):
        if p.is_file() and p.suffix in TEXT and ".git" not in p.parts:
            yield p, p.read_text(errors="replace")


def test_no_unscoped_accuracy_or_certification_claims():
    bad = []
    for p, text in _public_text():
        for pattern, why in OVERCLAIMS:
            for m in re.finditer(pattern, text, re.I):
                window = text[max(0, m.start() - 90): m.end() + 30]
                # A negated mention is the correct way to discuss these.
                if re.search(r"\bnot\b|\bno\b|never|prohibited|disclaim", window, re.I):
                    continue
                bad.append(f"{p.relative_to(ROOT)}: {m.group(0)!r} ({why})")
    assert not bad, "unscoped claims:\n  " + "\n  ".join(bad)


def test_the_scope_qualifier_travels_with_the_headline():
    """526/526 without its denominator is the number most likely to be quoted
    out of context, so the qualifier has to be adjacent, not merely elsewhere."""
    text = (ROOT / "README.md").read_text()
    i = text.find("526 / 526")
    assert i != -1
    window = text[max(0, i - 600): i + 900].lower()
    for token in ("v4.2.1", "grch38", "chr1", "not whole-genome accuracy"):
        assert token in window, f"scope qualifier {token!r} not adjacent to the aggregate"


def test_every_public_report_disclaims_clinical_status(sample):
    html = report_html(sample).lower()
    assert "not a diagnostic device" in html
    assert "not medical advice" in html
    assert "not clinical validation" in html


def test_reports_state_they_are_not_the_production_artifact(sample):
    html = report_html(sample).lower()
    assert "not the production" in html, (
        "the public/private boundary must be stated in the artifact itself, "
        "not only in the repository documentation")
