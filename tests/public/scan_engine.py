"""Generic disclosure-scan engine for the public verification suite.

This module knows how to search. It does not know what any withheld value is,
and it must never be given one: a list of withheld numbers in public code
publishes the numbers. Exact withheld-value screening is done before release by
a private gate that is not part of this repository (see docs/TESTING.md).

Two properties are deliberate, and test_scan_engine_controls.py pins both:

* A finding carries a rule id, a location and a count. It never carries the
  text that matched, so a failure message cannot echo what a rule exists to
  protect.
* A scan with no rules, or a rule that does not compile, is an error. A scanner
  that matches nothing reports a clean result forever.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Union

TEXT_SUFFIXES = frozenset({".md", ".json", ".html", ".yml", ".yaml", ".py", ".txt", ".svg"})


class EmptyRuleSet(ValueError):
    """A scan that cannot find anything is refused rather than reported clean."""


@dataclass(frozen=True)
class Rule:
    """Matches wherever `pattern` matches."""
    rule_id: str
    pattern: str
    flags: int = 0


@dataclass(frozen=True)
class ProximityRule:
    """Matches an `anchor` that has a `target` within `window` characters on either side.

    For categories recognisable by shape and context rather than by value."""
    rule_id: str
    anchor: str
    target: str
    window: int = 240
    flags: int = 0


@dataclass(frozen=True)
class Finding:
    """Where a rule matched and how often. Deliberately has no field for the matched text."""
    rule_id: str
    where: str
    count: int


AnyRule = Union[Rule, ProximityRule]


def _compile(pattern: str, flags: int, rule_id: str) -> re.Pattern:
    try:
        return re.compile(pattern, flags)
    except re.error:
        # re.error messages quote the pattern; report the rule id only.
        raise ValueError(f"rule {rule_id}: pattern does not compile") from None


def scan_text(text: str, rules: Iterable[AnyRule], where: str) -> list[Finding]:
    rules = list(rules)
    if not rules:
        raise EmptyRuleSet("no rules supplied")
    out: list[Finding] = []
    for rule in rules:
        if isinstance(rule, ProximityRule):
            anchor = _compile(rule.anchor, rule.flags, rule.rule_id)
            target = _compile(rule.target, rule.flags, rule.rule_id)
            n = sum(1 for m in anchor.finditer(text)
                    if target.search(text, max(0, m.start() - rule.window), m.end() + rule.window))
        else:
            n = sum(1 for _ in _compile(rule.pattern, rule.flags, rule.rule_id).finditer(text))
        if n:
            out.append(Finding(rule.rule_id, where, n))
    return out


def pdf_surfaces(path: Path) -> Iterator[tuple[str, str]]:
    """Every string a reader or a crawler can extract from a PDF, labelled by where it lives.

    Extracted text does not include a link's target, and neither includes the
    document properties or attachment names; each is its own surface."""
    import fitz  # PyMuPDF
    doc = fitz.open(path)
    try:
        yield "text", "\n".join(page.get_text() for page in doc)
        yield "link-targets", "\n".join(
            str(v) for page in doc for link in page.get_links()
            for k, v in link.items() if k in ("uri", "file", "name") and v)
        yield "metadata", "\n".join(str(v) for v in doc.metadata.values() if v)
        yield "attachments", "\n".join(doc.embfile_names())
    finally:
        doc.close()


def scan_pdf(path: Path, rules: Iterable[AnyRule], where: str) -> list[Finding]:
    rules = list(rules)
    out: list[Finding] = []
    for surface, text in pdf_surfaces(path):
        out += scan_text(text, rules, f"{where} [{surface}]")
    return out


def iter_text_files(root: Path, suffixes: frozenset[str] = TEXT_SUFFIXES,
                    skip: frozenset[Path] = frozenset()) -> Iterator[Path]:
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix in suffixes and ".git" not in p.parts and p.resolve() not in skip:
            yield p


def scan_tree(root: Path, rules: Iterable[AnyRule], suffixes: frozenset[str] = TEXT_SUFFIXES,
              skip: frozenset[Path] = frozenset(), pdfs: bool = True) -> list[Finding]:
    """Scan every text file and, unless `pdfs` is False, every PDF under `root` on all PDF surfaces."""
    rules = list(rules)
    out: list[Finding] = []
    for p in iter_text_files(root, suffixes, skip):
        out += scan_text(p.read_text(errors="replace"), rules, str(p.relative_to(root)))
    if pdfs:
        out += scan_pdfs(root, rules)
    return out


def scan_pdfs(root: Path, rules: Iterable[AnyRule]) -> list[Finding]:
    """Scan every PDF under `root` on all PDF surfaces. Needs PyMuPDF."""
    rules = list(rules)
    out: list[Finding] = []
    for p in sorted(root.rglob("*.pdf")):
        if ".git" not in p.parts:
            out += scan_pdf(p, rules, str(p.relative_to(root)))
    return out


def render(findings: Iterable[Finding]) -> str:
    """Failure text: rule ids, locations and counts only."""
    return "; ".join(f"{f.rule_id} x{f.count} in {f.where}" for f in findings) or "no findings"


# A category defined by shape and context, not by value: a currency or percent
# figure close to value-of-information vocabulary. It holds no withheld number
# and cannot confirm the absence of any particular one.
_CURRENCY = r"(?:US\$|USD|\$)\s?\d"
_FIGURE_WITH_UNIT = r"\d[\d,]*(?:\.\d+)?\s?(?:%|percent\b|per\s?cent\b|dollars?\b|USD\b)"
VOI_EXAMPLE_RULE = ProximityRule(
    "GEN-VOI-FIGURE",
    anchor=(r"\bexpected\s+value\s+of\s+(?:partial\s+)?(?:perfect|sample)\s+(?:parameter\s+)?information\b"
            r"|\bvalue[- ]of[- ]information\b|\bEVPI\b|\bEVSI\b|\bEVPPI\b|\bVOI\b"),
    target=_CURRENCY + "|" + _FIGURE_WITH_UNIT,
    window=240,
    flags=re.I,
)
