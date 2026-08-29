"""Stable semantic identifiers for economic pathways.

THE PROBLEM. One genomic finding is described three times by three code paths
that share no key:

* the **curated per-finding table** names the *action* — "Avoid clopidogrel
  non-response / stent thrombosis (MACE)";
* the **value-of-information parameterisation** names the *finding* — "CYP2C19
  Intermediate Metabolizer (IM) (clopidogrel)";
* the **condition pool** names the *disease* — "CAD".

Those strings describe one pathway and cannot be matched to each other, so the
report printed two different dollar figures for the same finding with no way for
a reader or a test to see that they referred to the same thing.

IDENTIFIERS ARE NOT DISPLAY TEXT. The first version of this module slugified the
``finding`` string, which works until someone improves the wording — at which
point every id silently changes and every reconciliation, test and cached
comparison breaks. Ids are built here from **semantic components** the producing
code already knows: gene, drug, phenotype, condition. Display strings can then
be rewritten freely.

Four levels, narrowest to broadest::

    finding_id            pgx:cyp2c19:intermediate_metabolizer
    action_id             pgx:clopidogrel:avoid_nonresponse
    condition_id          mace
    economic_pathway_id   pgx:cyp2c19:clopidogrel:mace

``economic_pathway_id`` is the join key: it identifies one genomic result acting
on one condition through one intervention, which is the unit both pricing paths
are trying to value.

MIGRATION FALLBACK. Extractors that do not yet supply semantic components get a
deterministic slug of the finding text, flagged ``legacy:``. Those ids are
explicitly unstable, and a test asserts that no canonical monetised record
depends on one.
"""
from __future__ import annotations

import re
import unicodedata

__all__ = [
    "LEGACY_PREFIX", "action_id", "condition_id", "economic_pathway_id",
    "finding_id", "gene_tokens", "is_legacy", "legacy_pathway_id", "slug",
]

LEGACY_PREFIX = "legacy:"

_GENE_RE = re.compile(
    r"\b(?:CYP[0-9][A-Z][0-9]+|SLCO1B1|TPMT|NUDT15|DPYD|UGT1A1|VKORC1|MTHFR|"
    r"APOE|PTPN22|BRCA[12]|MLH1|MSH[26]|PMS2|TP53|LDLR|PCSK9|APOB|KCNQ1|"
    r"HFE|F5|F2|SERPINA1|G6PD|HLA-[A-Z0-9*:]+|FTO|TCF7L2|LPA|CETP|COMT|"
    r"OPRM1|DRD2|ANKK1|CHRNA5|FAAH|PON1|GST[MPT][0-9]?|NAT2|ALDH2|ADH1B)\b")

# Phenotype text -> a short stable token. Free text varies ("Intermediate
# Metabolizer (IM)", "Intermediate metabolizer"); the token must not.
_PHENOTYPE_TOKENS = (
    ("poor metabol", "poor_metabolizer"),
    ("intermediate metabol", "intermediate_metabolizer"),
    ("rapid metabol", "rapid_metabolizer"),
    ("ultrarapid", "ultrarapid_metabolizer"),
    ("ultra-rapid", "ultrarapid_metabolizer"),
    ("normal metabol", "normal_metabolizer"),
    ("decreased function", "decreased_function"),
    ("intermediate function", "intermediate_function"),
    ("poor function", "poor_function"),
    ("increased function", "increased_function"),
    ("carrier", "carrier"),
    ("homozygous", "homozygous"),
    ("heterozygous", "heterozygous"),
)


def gene_tokens(text: str) -> list[str]:
    """Gene symbols appearing in a string, in order, de-duplicated."""
    seen: list[str] = []
    for m in _GENE_RE.finditer((text or "").upper()):
        if m.group(0) not in seen:
            seen.append(m.group(0))
    return seen


def slug(text: str, *, max_len: int = 48) -> str:
    """A stable, ASCII, lower-case token. Deterministic across platforms."""
    s = unicodedata.normalize("NFKD", str(text or ""))
    s = s.encode("ascii", "ignore").decode("ascii").lower()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return (s[:max_len].rstrip("_") or "unknown")


def _phenotype_token(phenotype: str) -> str:
    p = (phenotype or "").lower()
    for needle, token in _PHENOTYPE_TOKENS:
        if needle in p:
            return token
    return slug(phenotype, max_len=32) if p else "unspecified"


def finding_id(*, kind: str, gene: str = "", phenotype: str = "",
               variant: str = "") -> str:
    """What was found: ``pgx:cyp2c19:intermediate_metabolizer``."""
    parts = [slug(kind, max_len=16)]
    if gene:
        parts.append(slug(gene, max_len=24))
    if variant:
        parts.append(slug(variant, max_len=32))
    parts.append(_phenotype_token(phenotype))
    return ":".join(p for p in parts if p)


def action_id(*, kind: str, drug: str = "", action: str = "") -> str:
    """What could be done about it: ``pgx:clopidogrel:avoid_nonresponse``."""
    parts = [slug(kind, max_len=16)]
    if drug:
        parts.append(slug(drug, max_len=24))
    if action:
        parts.append(slug(action, max_len=32))
    return ":".join(p for p in parts if p)


def condition_id(condition: str) -> str:
    """The modelled outcome: ``mace``, ``myopathy``, ``alzheimer``."""
    return slug(condition, max_len=32)


def economic_pathway_id(*, kind: str, gene: str = "", drug: str = "",
                        condition: str = "", phenotype: str = "",
                        variant: str = "") -> str:
    """The join key: one finding, one intervention, one condition.

    ``pgx:cyp2c19:clopidogrel:mace`` — the unit both pricing paths value, and
    the level at which their disagreement is meaningful.
    """
    parts = [slug(kind, max_len=16)]
    for v, n in ((gene, 24), (drug, 24), (condition, 32)):
        if v:
            parts.append(slug(v, max_len=n))
    if not gene and not drug and not condition:
        parts.append(_phenotype_token(phenotype) if phenotype
                     else slug(variant, max_len=32))
    return ":".join(p for p in parts if p)


def legacy_pathway_id(finding: str, category: str = "",
                      pool_hint: str = "") -> str:
    """Migration fallback for extractors with no semantic components yet.

    Deterministic, but derived from display text and therefore **unstable**:
    rewording a finding changes the id. Flagged so callers and tests can tell
    the difference.
    """
    stem = slug(finding)
    if pool_hint:
        return f"{LEGACY_PREFIX}{slug(pool_hint, max_len=32)}:{stem}"
    genes = gene_tokens(finding)
    if genes:
        return f"{LEGACY_PREFIX}{slug(genes[0], max_len=24)}:{stem}"
    if category:
        return f"{LEGACY_PREFIX}{slug(category, max_len=24)}:{stem}"
    return f"{LEGACY_PREFIX}{stem}"


def is_legacy(pathway: str) -> bool:
    return str(pathway or "").startswith(LEGACY_PREFIX)
