"""Controls for the generic disclosure scanner, using synthetic canaries only.

The public suite does not contain, and must not contain, any withheld value:
a deny-list of exact numbers publishes the numbers, and a failing assertion
that echoes its pattern or its match publishes them again. What this module
checks instead is value-free. It proves the scanner can find what a rule
describes (positive controls), does not find what it does not (negative
controls), reads every surface of a PDF, and cannot echo matched text.

Exact withheld-value screening happens in a private release gate that this
repository does not contain. Nothing here can confirm that any particular value
is absent, and a passing run must not be read as saying so.

Canaries are assembled at run time so that this file never holds a figure-shaped
token next to value-of-information vocabulary (the scan below would flag its own
source), and they are chosen so that they cannot coincide with any real figure.
"""
import re

import pytest
from conftest import ROOT
from scan_engine import (VOI_EXAMPLE_RULE, EmptyRuleSet, Finding, ProximityRule, Rule,
                         render, scan_pdf, scan_pdfs, scan_text, scan_tree)

LEAD = "The expected value of perfect information was "


def usd_canary() -> str:
    return "$" + "777777" + "." + "77"


def pct_canary() -> str:
    return "77" + "." + "777" + "%"


def canary_rule() -> Rule:
    return Rule("CTL-CANARY", r"(?<![\d.,])" + re.escape(usd_canary()) + r"(?![\d])")


def make_pdf(path, text="", link=None, title=None):
    fitz = pytest.importorskip("fitz", reason="PyMuPDF not installed")
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 72), text or "control page")
    if link:
        page.insert_link({"kind": fitz.LINK_URI, "from": fitz.Rect(72, 60, 300, 80), "uri": link})
    if title:
        doc.set_metadata({"title": title})
    doc.save(path)
    doc.close()


# ── the generic value-of-information rule: positive and negative controls ─────

def test_voi_rule_fires_on_a_synthetic_currency_figure():
    found = scan_text(LEAD + usd_canary() + " per person.", [VOI_EXAMPLE_RULE], "synthetic")
    assert [(f.rule_id, f.count) for f in found] == [("GEN-VOI-FIGURE", 1)]


def test_voi_rule_fires_on_a_synthetic_percent_figure():
    found = scan_text("EVSI was estimated at " + pct_canary() + " of the budget.", [VOI_EXAMPLE_RULE], "synthetic")
    assert len(found) == 1


@pytest.mark.parametrize("text", [
    "Value-of-information analysis asks whether more evidence could change the decision.",
    "The bound 0 <= EVSI <= EVPI holds for any decision problem.",
    "A threshold of " + "$" + "100,000 per QALY is a decision-maker's parameter.",
    "Coverage was " + "100" + "% of the rows checked.",
], ids=["vocabulary-only", "inequality", "figure-without-voi-vocabulary", "percent-without-voi-vocabulary"])
def test_voi_rule_ignores_benign_text(text):
    assert scan_text(text, [VOI_EXAMPLE_RULE], "synthetic") == []


def test_voi_rule_ignores_a_figure_beyond_its_window():
    far = LEAD + "stated in words only. " + ("Unrelated filler sentence. " * 30) + usd_canary()
    assert scan_text(far, [VOI_EXAMPLE_RULE], "synthetic") == []


@pytest.fixture(scope="module")
def scanner_is_live():
    """Controls run before the tree scan that depends on them. If the scanner
    cannot fire, its silence on the real tree proves nothing."""
    assert scan_text(LEAD + usd_canary(), [VOI_EXAMPLE_RULE], "positive-control"), "positive control failed"
    assert not scan_text("No figures here.", [VOI_EXAMPLE_RULE], "negative-control"), "negative control failed"


def test_public_text_files_have_no_numeric_value_of_information_example(scanner_is_live):
    """Value-free by design. Numerical value-of-information examples are withheld
    from this release (docs/RELEASE-STATUS.md). This looks for the category, by
    shape and context, in every text file in the tree. It cannot say whether a
    specific withheld value is absent."""
    findings = scan_tree(ROOT, [VOI_EXAMPLE_RULE], pdfs=False)
    assert not findings, render(findings)


def test_public_pdfs_have_no_numeric_value_of_information_example(scanner_is_live):
    """The same check on every surface of every PDF in the tree."""
    pytest.importorskip("fitz", reason="PyMuPDF not installed")
    findings = scan_pdfs(ROOT, [VOI_EXAMPLE_RULE])
    assert not findings, render(findings)


# ── the engine: surfaces, hygiene, refusals ───────────────────────────────────

def test_tree_scan_reads_text_files_and_pdfs(tmp_path):
    (tmp_path / "note.md").write_text("Synthetic: " + usd_canary() + "\n")
    make_pdf(tmp_path / "doc.pdf", text="Synthetic: " + usd_canary())
    found = scan_tree(tmp_path, [canary_rule()])
    assert {f.where for f in found} == {"note.md", "doc.pdf [text]"}


def test_pdf_scan_reads_text_link_targets_and_metadata_separately(tmp_path):
    canary = usd_canary()
    make_pdf(tmp_path / "text.pdf", text="x " + canary)
    make_pdf(tmp_path / "link.pdf", link="https://example.invalid/?q=" + canary)
    make_pdf(tmp_path / "meta.pdf", title="x " + canary)
    seen = {}
    for name in ("text", "link", "meta"):
        seen[name] = {f.where.split("[")[1].rstrip("]") for f in scan_pdf(tmp_path / f"{name}.pdf", [canary_rule()], name)}
    assert seen == {"text": {"text"}, "link": {"link-targets"}, "meta": {"metadata"}}


def test_a_clean_pdf_yields_no_findings(tmp_path):
    make_pdf(tmp_path / "clean.pdf", text="Nothing to see here.", link="https://example.invalid/", title="Clean")
    assert scan_pdf(tmp_path / "clean.pdf", [canary_rule()], "clean") == []


def test_near_misses_of_a_canary_are_not_findings():
    rule = canary_rule()
    lookalikes = ["$" + "777777.78", "$" + "1777777.77", "$" + "777777.775", "777777.77", "x" + usd_canary()[1:]]
    for text in lookalikes:
        assert scan_text(text, [rule], "synthetic") == [], "a boundary was not respected"
    assert scan_text("cost " + usd_canary() + ".", [rule], "synthetic")      # sentence-final period still matches


def test_findings_carry_no_matched_text():
    """The reason a failure cannot echo a protected value: there is no field to hold it."""
    assert set(Finding.__dataclass_fields__) == {"rule_id", "where", "count"}
    found = scan_text("leak " + usd_canary(), [canary_rule()], "somewhere")
    assert usd_canary() not in render(found) and usd_canary() not in repr(found)


def test_an_empty_rule_set_is_refused():
    with pytest.raises(EmptyRuleSet):
        scan_text("anything", [], "synthetic")
    with pytest.raises(EmptyRuleSet):
        scan_tree(ROOT, [])


def test_an_uncompilable_rule_does_not_echo_its_pattern():
    distinctive = "(unclosed-" + "pattern-marker"
    with pytest.raises(ValueError) as err:
        scan_text("anything", [Rule("CTL-BAD", distinctive)], "synthetic")
    assert "CTL-BAD" in str(err.value) and distinctive not in str(err.value)
    with pytest.raises(ValueError) as err2:
        scan_text("anything", [ProximityRule("CTL-BAD2", "ok", distinctive)], "synthetic")
    assert distinctive not in str(err2.value)


# ── regressions this change exists to prevent ─────────────────────────────────

# Shapes of a figure written out as a regular-expression literal. Detected by
# shape, so the check holds no value and cannot itself become a deny-list.
FIGURE_REGEX_LITERAL_SHAPES = (r"\\\$\d+\\\.\d+", r"\d+\\\.\d+(?:\\s\??)?\s?%")


def test_public_code_holds_no_figure_deny_list():
    hits = []
    for p in sorted((ROOT / "tests").rglob("*.py")) + sorted((ROOT / "tools").rglob("*.py")):
        text = p.read_text(errors="replace")
        if any(re.search(shape, text) for shape in FIGURE_REGEX_LITERAL_SHAPES):
            hits.append(str(p.relative_to(ROOT)))
    assert not hits, "a figure written as a regex literal appears in " + ", ".join(hits)


def test_no_collected_test_id_embeds_a_figure(request):
    """Parametrized test ids are built from their arguments and are printed in
    verbose runs, failure summaries and CI logs. An id must never carry a figure."""
    shape = re.compile(r"\$\s?\d+\.\d+|\d+\.\d+(?:s\??|\s)?%")
    # pytest escapes backslashes in parametrized ids, so compare with them removed
    bad = sum(1 for item in request.session.items if shape.search(item.nodeid.replace("\\", "")))
    assert bad == 0, f"{bad} collected test id(s) embed a figure-shaped token"
