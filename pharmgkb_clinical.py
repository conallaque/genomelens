"""PharmGKB / ClinPGx clinical-variant annotations.

Surfaces the `cpic_data/clinicalVariants.tsv` dataset (downloaded from
ClinPGx/PharmGKB) that no other module reads. For every rsID in that table
that is typed on the user's chip, we report the annotation the database
carries: gene, the drug(s) it is annotated for, the ClinPGx clinical-
annotation evidence level, and the associated phenotype.

IMPORTANT — accuracy framing (do not weaken this in the renderer):
  * The evidence levels (1A/1B/2A/2B/3/4) are the **ClinPGx / PharmGKB
    clinical-annotation** level system (see
    https://www.clinpgx.org/page/clinAnnLevels), NOT CPIC guideline strength.
    Level 1A/1B are strong/replicated; Level 3/4 are weak/unreplicated
    (single studies, case reports, in-vitro, or non-significant).
  * The source table has NO risk allele and NO direction of effect. So this
    module reports "you carry a genotype at an annotated position" — it does
    NOT make a genotype→phenotype call. Direction/dosing for the clinically
    actionable genes comes from the dedicated `pgx.py` star-allele module and
    the curated 217-drug database; this panel is complementary breadth, and
    the dedicated PGx section remains the authority for those genes.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

_TSV_PATH = Path(__file__).resolve().parent / "cpic_data" / "clinicalVariants.tsv"
_RSID_RE = re.compile(r"rs\d+")

# ClinPGx clinical-annotation level ordering (lower rank = stronger evidence).
_LEVEL_RANK = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "4": 5}
_HIGH_LEVELS = {"1A", "1B", "2A", "2B"}

# Genes with a dedicated star-allele analysis in pgx.py — surfaced here too for
# completeness, but flagged so the user knows the PGx section is authoritative.
_PGX_GENES = {
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "TPMT", "NUDT15",
    "SLCO1B1", "VKORC1", "UGT1A1", "DPYD",
}

_TSV_CACHE: Optional[List[Dict]] = None


def _load_table() -> List[Dict]:
    """Parse clinicalVariants.tsv into rows keyed by a single clean rsID.

    Rows whose `variant` is not a single rsID (star alleles, HGVS, multi-allele
    haplotypes) are skipped — those are handled by the star-allele PGx module.
    """
    global _TSV_CACHE
    if _TSV_CACHE is not None:
        return _TSV_CACHE
    rows: List[Dict] = []
    if not _TSV_PATH.exists():
        _TSV_CACHE = rows
        return rows
    with _TSV_PATH.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            variant = (r.get("variant") or "").strip()
            if not _RSID_RE.fullmatch(variant):
                continue
            gene = (r.get("gene") or "").strip()
            if not gene:
                continue
            level = (r.get("level of evidence") or "").strip()
            drugs = [c.strip() for c in (r.get("chemicals") or "").split(",") if c.strip()]
            phenos = [p.strip() for p in (r.get("phenotypes") or "").split(",") if p.strip()]
            rows.append({
                "rsid": variant,
                "gene": gene,
                "type": (r.get("type") or "").strip(),
                "level": level,
                "drugs": drugs,
                "phenotypes": phenos,
            })
    _TSV_CACHE = rows
    return rows


def analyze_pharmgkb_clinical(snps_df: Optional[pd.DataFrame]) -> Dict:
    """Report ClinPGx/PharmGKB clinical annotations for typed rsIDs.

    Returns a dict with `high` (Level 1A/1B/2A/2B) and `low` (Level 3/4)
    variant lists; each entry groups all annotation rows for one typed rsID.
    """
    table = _load_table()
    if snps_df is None or snps_df.empty or not table:
        return {"available": False, "reason": "No genotype data or dataset not present."}

    index = snps_df.index
    by_rsid: Dict[str, Dict] = {}
    for row in table:
        rsid = row["rsid"]
        if rsid not in index:
            continue
        entry = by_rsid.get(rsid)
        if entry is None:
            loc = snps_df.loc[rsid]
            # A duplicate rsID in the chip file makes .loc return a DataFrame;
            # take the first typed call in that case.
            if isinstance(loc, pd.DataFrame):
                loc = loc.iloc[0]
            gt = loc.get("genotype")
            gt = str(gt).upper().strip() if gt is not None else ""
            if gt in ("", "NAN", "--"):
                continue
            entry = {
                "rsid": rsid,
                "gene": row["gene"],
                "genotype": gt,
                "annotations": [],
                "pgx_gene": row["gene"] in _PGX_GENES,
            }
            by_rsid[rsid] = entry
        entry["annotations"].append({
            "drugs": row["drugs"],
            "level": row["level"],
            "phenotypes": row["phenotypes"],
            "type": row["type"],
        })

    # Best (strongest) evidence level per rsID drives tiering + sort order.
    def _best_rank(entry: Dict) -> int:
        return min(
            (_LEVEL_RANK.get(a["level"], 99) for a in entry["annotations"]),
            default=99,
        )

    high, low = [], []
    for entry in by_rsid.values():
        entry["annotations"].sort(key=lambda a: _LEVEL_RANK.get(a["level"], 99))
        entry["best_level"] = entry["annotations"][0]["level"] if entry["annotations"] else ""
        (high if _best_rank(entry) <= _LEVEL_RANK["2B"] else low).append(entry)

    high.sort(key=lambda e: (_best_rank(e), e["gene"]))
    low.sort(key=lambda e: (_best_rank(e), e["gene"]))

    n_rows = sum(len(e["annotations"]) for e in by_rsid.values())
    drugs = {d for e in by_rsid.values() for a in e["annotations"] for d in a["drugs"]}
    return {
        "available": True,
        "source": "ClinPGx / PharmGKB clinical variant annotations",
        "high": high,
        "low": low,
        "n_typed_variants": len(by_rsid),
        "n_high": len(high),
        "n_low": len(low),
        "n_annotation_rows": n_rows,
        "n_drugs": len(drugs),
        "n_dataset_rsid_rows": len(table),
    }
