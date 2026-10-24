"""What the public artifacts are not allowed to claim.

The scoped concordance claim is defensible. The unscoped version of the same
sentence is not, and the distance between them is one word.
"""
import re

import fitz
from conftest import ROOT

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


def test_every_published_report_states_its_research_status_and_that_it_is_proprietary():
    pdfs = sorted((ROOT / "demonstrations" / "v1.11").rglob("*.pdf"))
    assert len(pdfs) >= 17, f"expected the demonstration reports, the dossier and the featured report; found {len(pdfs)}"
    for p in pdfs:
        text = " ".join(" ".join(pg.get_text() for pg in fitz.open(p)).split())
        assert re.search(r"not (medical|clinical) or financial advice", text), f"{p.name} does not state its research status"
        assert "not open-source software" in text, f"{p.name} does not state that it is proprietary"
