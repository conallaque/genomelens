"""Screen the most-commonly-prescribed medications against this genome.

Takes the curated `data/top_prescribed_drugs.json` reference menu and, for each drug,
attaches its pharmacogenomic relevance computed from REAL bundled data:

  * `data/drug_database.json`  — 217 CPIC/DPWG drug↔gene↔marker↔dosing records.
  * `cpic_data/drugs.tsv` — ClinPGx/PharmGKB per-drug metadata (Top CPIC pair
    level, top clinical-annotation level).
  * the per-gene metabolizer phenotypes from `pgx.analyze_pgx` (the user's
    actual star-allele calls).

Nothing about a drug's pharmacogenomics is taken from the prescribed-drugs
file — that only supplies names/classes. A drug with no entry in the CPIC /
PharmGKB data is reported as "no known pharmacogenomic interaction," never as a
fabricated finding. Direction/dosing for actionable genes remains the authority
of the dedicated Pharmacogenomics (pgx.py) section.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pandas as pd

# .parent.parent: this module moved into the pgx package, but the data file
# it reads stays at the repository root. Left at .parent it would resolve
# inside pgx/ and the loader degrades quietly on a missing file.
_DIR = Path(__file__).resolve().parent.parent
_TOP_PATH = _DIR / "data/top_prescribed_drugs.json"
_CPIC_DB_PATH = _DIR / "data/drug_database.json"
_PGKB_TSV_PATH = _DIR / "cpic_data" / "drugs.tsv"

# Metabolizer phenotypes that can alter standard dosing (from pgx.py).
_ACTIONABLE_CODES = {"PM", "IM", "UM", "RM", "POS"}

_CACHE: dict[str, object] = {}


def _base_name(name: str) -> str:
    """Normalise a drug name for matching: lowercase, drop combination parts
    and parenthetical qualifiers (e.g. 'duloxetine DR' -> 'duloxetine')."""
    n = name.lower().split("/")[0].split("(")[0].strip()
    # strip common formulation suffixes
    for suf in (" dr", " er", " xr", " la", " odt", " ophthalmic", " topical",
                " nasal", " transdermal", " sodium", " fumarate", " chewable"):
        if n.endswith(suf):
            n = n[: -len(suf)].strip()
    return n


def _load_reference() -> list[dict]:
    if "top" not in _CACHE:
        _CACHE["top"] = (
            json.loads(_TOP_PATH.read_text()).get("drugs", [])
            if _TOP_PATH.exists() else []
        )
    return _CACHE["top"]  # type: ignore[return-value]


def _load_cpic_db() -> dict[str, dict]:
    """base drug name -> merged CPIC record (genes, markers, dosing)."""
    if "cpic" in _CACHE:
        return _CACHE["cpic"]  # type: ignore[return-value]
    out: dict[str, dict] = {}
    if _CPIC_DB_PATH.exists():
        for e in json.loads(_CPIC_DB_PATH.read_text()).get("drugs", []):
            key = _base_name(e.get("drug_name", ""))
            if not key:
                continue
            rec = out.setdefault(key, {"genes": set(), "markers": set(), "dosing": ""})
            rec["genes"].update(e.get("genes", []) or [])
            rec["markers"].update(e.get("snp_markers", []) or [])
            if not rec["dosing"] and e.get("dosing_recommendation"):
                rec["dosing"] = e["dosing_recommendation"]
    _CACHE["cpic"] = out
    return out


def _load_pgkb() -> dict[str, dict]:
    """base drug name -> {cpic_level, clin_level} from PharmGKB drugs.tsv."""
    if "pgkb" in _CACHE:
        return _CACHE["pgkb"]  # type: ignore[return-value]
    out: dict[str, dict] = {}
    if _PGKB_TSV_PATH.exists():
        try:
            csv.field_size_limit(sys.maxsize)
        except (OverflowError, ValueError):
            csv.field_size_limit(2**31 - 1)
        with _PGKB_TSV_PATH.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                names = set()
                for f in ("Name", "Generic Names"):
                    for n in (r.get(f) or "").split(","):
                        b = _base_name(n)
                        if b:
                            names.add(b)
                meta = {
                    "cpic_level": (r.get("Top CPIC Pairs Level") or "").strip(),
                    "clin_level": (r.get("Top Clinical Annotation Level") or "").strip(),
                }
                for n in names:
                    if n not in out or (meta["cpic_level"] and not out[n]["cpic_level"]):
                        out[n] = meta
    _CACHE["pgkb"] = out
    return out


def analyze_top_drugs(
    snps_df: pd.DataFrame | None,
    pgx_result: dict | None = None,
) -> dict:
    """Classify every reference drug by its pharmacogenomic relevance to this
    genome. Returns tiered lists + counts.
    """
    ref = _load_reference()
    if not ref:
        return {"available": False, "reason": "data/top_prescribed_drugs.json not present."}

    cpic = _load_cpic_db()
    pgkb = _load_pgkb()
    per_gene = (pgx_result or {}).get("per_gene", {}) if pgx_result else {}
    index = snps_df.index if snps_df is not None else pd.Index([])

    actionable: list[dict] = []       # your genotype may change dosing
    typed_normal: list[dict] = []     # PGx drug, your relevant gene typed & normal
    pgx_relevant: list[dict] = []     # PGx data exists, your gene unresolved
    no_pgx: list[dict] = []           # no PGx data in bundled databases

    for d in ref:
        key = _base_name(d["generic"])
        cp = cpic.get(key)
        pk = pgkb.get(key)
        entry = {
            "generic": d["generic"],
            "brand": d.get("brand", ""),
            "class": d.get("class", ""),
            "genes": sorted(cp["genes"]) if cp else [],
            "cpic_level": (pk or {}).get("cpic_level", ""),
            "clin_level": (pk or {}).get("clin_level", ""),
            "dosing": cp["dosing"] if cp else "",
        }

        if not cp and not pk:
            no_pgx.append(entry)
            continue

        # User phenotype across this drug's CPIC genes.
        gene_phenos = []
        any_actionable = any_typed = False
        for g in entry["genes"]:
            res = per_gene.get(g)
            if not res:
                continue
            code = res.get("phenotype_code")
            gene_phenos.append({"gene": g, "phenotype": res.get("phenotype", ""),
                                "code": code})
            if code and code not in ("IND", "NC"):
                any_typed = True
            if code in _ACTIONABLE_CODES:
                any_actionable = True
        # Also count the relevant markers typed on the chip.
        markers_typed = sum(1 for m in (cp["markers"] if cp else set()) if m in index)
        entry["gene_phenotypes"] = gene_phenos
        entry["markers_typed"] = markers_typed

        if any_actionable:
            actionable.append(entry)
        elif any_typed:
            typed_normal.append(entry)
        else:
            pgx_relevant.append(entry)

    # Sort: strongest CPIC evidence first within each tier.
    def _lvl_key(e):
        order = {"A": 0, "A/B": 1, "B": 2, "B/C": 3, "C": 4, "D": 5, "": 9}
        return order.get(e.get("cpic_level", ""), 9)

    for lst in (actionable, typed_normal, pgx_relevant):
        lst.sort(key=lambda e: (_lvl_key(e), e["generic"]))
    no_pgx.sort(key=lambda e: e["generic"])

    return {
        "available": True,
        "source": "Most-prescribed medications screened against CPIC + ClinPGx/PharmGKB",
        "n_screened": len(ref),
        "actionable": actionable,
        "typed_normal": typed_normal,
        "pgx_relevant": pgx_relevant,
        "no_pgx": no_pgx,
        "n_actionable": len(actionable),
        "n_typed_normal": len(typed_normal),
        "n_pgx_relevant": len(pgx_relevant),
        "n_no_pgx": len(no_pgx),
        "n_with_pgx": len(actionable) + len(typed_normal) + len(pgx_relevant),
    }
