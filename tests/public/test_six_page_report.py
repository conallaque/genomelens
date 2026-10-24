"""Checks on the featured public six-page report (demonstrations/v1.11/featured/): what it may and may not contain. They read only files committed to
this repository. The number check is the rule the report was built under: every significant number must already be public, or be calculated from
printed terms and shown with its calculation."""
import pathlib
import re

import fitz

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEMO = ROOT / "demonstrations" / "v1.11"
FEATURED = DEMO / "featured" / "GenomeLens_public_six_page_HEOR_report_hu4AE0F2.pdf"
README = (ROOT / "README.md").read_text(encoding="utf-8")
NUM = re.compile(r"\d[\d,]*(?:\.\d+)?")


def _text(path) -> str:
    return "\n".join(p.get_text() for p in fitz.open(path))


def _toks(t: str) -> set:
    return {x.rstrip(",").replace(",", "") for x in NUM.findall(t.replace("−", "-"))}


def _significant(tok: str) -> bool:
    if "." in tok:
        return True
    n = int(tok)
    return n >= 100 and not 1990 <= n <= 2035          # a year is a fact about a cited paper


def _known(tok: str, public: set) -> bool:
    if tok in public:
        return True
    return "." in tok and any("." in p and len(p) >= 4 and p[:4] == tok[:4] and (p.startswith(tok) or tok.startswith(p)) for p in public)


def _public_corpus() -> str:
    out = "\n".join(_text(p) for p in (DEMO / "GenomeLens_public_demonstration_hu4AE0F2.pdf", DEMO / "GenomeLens_public_demonstration_hu4AE0F2_age30.pdf"))
    for f in [ROOT / "README.md", ROOT / "reproduction" / "EXPECTED_OUTPUT.txt", ROOT / "reproduction" / "README.md"] + sorted((ROOT / "docs").glob("*.md")):
        out += "\n" + f.read_text(encoding="utf-8")
    return out


def _without_axis_labels(text: str) -> str:
    return "\n".join(l for l in text.splitlines() if not re.fullmatch(r"\s*(?:\d+(?:\.\d+)?k|0\.\d{1,2}|1(?:\.0)?)\s*", l))


def test_the_featured_report_exists_and_the_readme_links_it_with_its_preview():
    assert FEATURED.is_file()
    rel = "demonstrations/v1.11/featured/GenomeLens_public_six_page_HEOR_report_hu4AE0F2.pdf"
    assert README.count(f"]({rel})") >= 2, "the README should link the report from its first screen and from its featured section"
    assert "demonstrations/v1.11/previews/featured_six_page_page1.png" in README
    assert (DEMO / "previews" / "featured_six_page_page1.png").is_file()


def test_the_featured_report_has_exactly_six_pages_and_only_its_title_as_metadata():
    d = fitz.open(FEATURED)
    assert d.page_count == 6
    meta = {k: v for k, v in (d.metadata or {}).items() if v and k not in ("format", "title", "encryption")}
    assert not meta, meta
    assert not d.get_xml_metadata()
    # the only links are the citation links of its PMIDs (a PubMed page); no other target, no file link
    links = [l for p in d for l in p.get_links()]
    assert all(re.fullmatch(r"https://pubmed\.ncbi\.nlm\.nih\.gov/\d{6,9}/?", str(l.get("uri") or "")) and int(l.get("kind", -1)) == fitz.LINK_URI for l in links), links
    assert not [a for p in d for a in (p.annots() or ())]
    text = "\n".join(p.get_text() for p in d)
    for bad in ("file:" + "///", "/Us" + "ers/", "/private/" + "tmp", "Draft " + "notice", "{" + "{"):
        assert bad not in text, bad


def test_the_featured_report_states_that_it_is_hypothetical_proprietary_and_not_validated():
    text = " ".join(_text(FEATURED).split())
    assert "ASSUMPTION-DRIVEN SCENARIO" in text and "Model basis: hypothetical" in text
    assert "Proprietary." in text and "not open-source software" in text
    assert "none of the decision models has been validated against observed outcomes" in text
    assert "The model's inputs, structure and registry are proprietary and are not published" in text


def test_the_featured_report_discloses_that_its_comparator_is_an_assumption_that_overstates_the_increment():
    text = " ".join(_text(FEATURED).split())
    assert "an unidentified carrier would receive average-risk screening, so surveillance against none overstates the increment" in text
    assert "its comparator is assumed" in text


def test_every_significant_number_in_the_featured_report_is_public_or_calculated_from_printed_terms():
    public = _toks(_public_corpus())
    text = _without_axis_labels(_text(FEATURED))
    mine = {t for t in _toks(text) if _significant(t)}
    derived = set()
    m = re.search(r"([\d.]+) QALYs [×x] \$([\d,]+) [−-] \$([\d,]+)", text)
    if m:       # the sensitivity case at a second threshold is arithmetic on two printed terms
        dq, wtp, dc = m.group(1), m.group(2).replace(",", ""), m.group(3).replace(",", "")
        assert _known(dq, public) and _known(dc, public), "the calculation uses terms the public reports do not print"
        derived = {str(round(float(dq) * float(wtp) - float(dc))), wtp}
    unknown = sorted(t for t in mine if not _known(t, public) and t not in derived)
    assert not unknown, f"numbers in the featured report that no public report, document or recorded output prints: {unknown}"
    assert len(mine) > 20          # the check read the figures it is meant to guard
