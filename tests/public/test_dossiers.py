"""Checks on the public HEOR technical dossier (demonstrations/v1.11/dossiers/): what they may and may not contain. They read only files committed to
this repository. The number check is the rule the dossiers were built under: every significant number must already be public."""
import glob
import pathlib
import re

import fitz
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEMO = ROOT / "demonstrations" / "v1.11"
DOSSIERS = DEMO / "dossiers"
README = (ROOT / "README.md").read_text(encoding="utf-8")
PAIRS = {"GenomeLens_public_HEOR_dossier_hu4AE0F2_age30.pdf": "GenomeLens_public_demonstration_hu4AE0F2_age30.pdf"}
IDS = [k.split("dossier_")[-1][:-4] for k in PAIRS]
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


def _public_corpus(report: str) -> str:
    # the genome's public reports at both stated ages
    base = report.replace("_age30.pdf", ".pdf")
    out = "\n".join(_text(DEMO / r) for r in (base, base[:-4] + "_age30.pdf") if (DEMO / r).exists())
    for f in [ROOT / "README.md", ROOT / "reproduction" / "EXPECTED_OUTPUT.txt", ROOT / "reproduction" / "README.md"] + sorted((ROOT / "docs").glob("*.md")):
        out += "\n" + f.read_text(encoding="utf-8")
    return out


def _without_axis_labels(text: str) -> str:
    """The plane's axis ticks (a lone 0.5 or 10k on a line) are scale marks, not figures."""
    return "\n".join(l for l in text.splitlines() if not re.fullmatch(r"\s*(?:\d+(?:\.\d+)?k|0\.\d{1,2}|1(?:\.0)?)\s*", l))


def test_the_dossiers_exist_and_the_readme_links_each_one():
    assert sorted(p.name for p in DOSSIERS.glob("*.pdf")) == sorted(PAIRS)
    for name in PAIRS:
        assert f"demonstrations/v1.11/dossiers/{name}" in README, name


@pytest.mark.parametrize("name", list(PAIRS), ids=IDS)
def test_a_dossier_carries_only_its_title_and_no_link_annotation_or_path(name):
    d = fitz.open(DOSSIERS / name)
    meta = {k: v for k, v in (d.metadata or {}).items() if v and k not in ("format", "title", "encryption")}
    assert not meta, meta
    assert not d.get_xml_metadata()
    assert d.page_count == 8
    text = "\n".join(p.get_text() for p in d)
    assert not [l for p in d for l in p.get_links()]
    for bad in ("file:" + "///", "/Us" + "ers/", "/private/" + "tmp", "Draft " + "notice", "{" + "{"):
        assert bad not in text, bad


@pytest.mark.parametrize("name", list(PAIRS), ids=IDS)
def test_a_dossier_says_it_is_proprietary_and_that_the_inputs_are_not_published(name):
    text = " ".join(_text(DOSSIERS / name).split())
    assert "Proprietary." in text and "not a description sufficient to reproduce the analysis" in text
    assert "inputs and structure are not published" in text
    assert "No model here has been validated against observed outcomes" in text


@pytest.mark.parametrize("name", list(PAIRS), ids=IDS)
def test_every_significant_number_in_a_dossier_is_already_public(name):
    public = _toks(_public_corpus(PAIRS[name]))
    mine = {t for t in _toks(_without_axis_labels(_text(DOSSIERS / name))) if _significant(t)}
    unknown = sorted(t for t in mine if not _known(t, public))
    assert not unknown, f"numbers in the dossier that no public report, document or recorded output prints: {unknown}"
    assert len(mine) > 20          # the check read the figures it is meant to guard
