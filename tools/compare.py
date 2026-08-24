"""
Report Comparison / Changelog
=============================

Diffs two tier1_results.json snapshots and emits a structured changelog
showing what's changed between two runs:

  * Tool / database version bumps
  * New variants matched (i.e. the SNP database grew)
  * Variants no longer matched (dropped from DB, or chip change)
  * Risk-level changes (significance up/down)
  * Risk-allele dosage changes (rare unless chip changed)
  * Polygenic risk score tier shifts
  * PGx phenotype changes
  * New / removed categories

Used via `python analyze.py file.txt --compare path/to/prev_tier1_results.json`.
"""

from __future__ import annotations

import json
from pathlib import Path


def _load(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def _variants_by_rsid(data: dict) -> dict[str, dict]:
    return {v["rsid"]: v for v in data.get("variants", [])}


def diff_runs(prev_path: Path, curr_path: Path) -> dict:
    """Compute the diff between two tier1_results.json files."""
    try:
        prev = _load(prev_path)
    except Exception as e:
        return {"error": f"Could not load previous report: {e}"}
    try:
        curr = _load(curr_path)
    except Exception as e:
        return {"error": f"Could not load current report: {e}"}

    prev_v = _variants_by_rsid(prev)
    curr_v = _variants_by_rsid(curr)

    new_rsids = sorted(set(curr_v) - set(prev_v))
    removed_rsids = sorted(set(prev_v) - set(curr_v))
    common_rsids = sorted(set(prev_v) & set(curr_v))

    new_variants: list[dict] = []
    for r in new_rsids:
        v = curr_v[r]
        new_variants.append({
            "rsid": r,
            "gene": v.get("gene"),
            "variant_name": v.get("variant_name"),
            "category": v.get("category"),
            "risk_copies": v.get("risk_copies"),
            "significance": v.get("significance"),
        })

    removed_variants: list[dict] = []
    for r in removed_rsids:
        v = prev_v[r]
        removed_variants.append({
            "rsid": r,
            "gene": v.get("gene"),
            "variant_name": v.get("variant_name"),
            "category": v.get("category"),
        })

    changed_significance: list[dict] = []
    changed_category: list[dict] = []
    changed_risk_copies: list[dict] = []
    for r in common_rsids:
        a, b = prev_v[r], curr_v[r]
        if a.get("significance") != b.get("significance"):
            changed_significance.append({
                "rsid": r, "gene": b.get("gene"),
                "from": a.get("significance"), "to": b.get("significance"),
            })
        if a.get("category") != b.get("category"):
            changed_category.append({
                "rsid": r, "gene": b.get("gene"),
                "from": a.get("category"), "to": b.get("category"),
            })
        if a.get("risk_copies") != b.get("risk_copies"):
            changed_risk_copies.append({
                "rsid": r, "gene": b.get("gene"),
                "from": a.get("risk_copies"), "to": b.get("risk_copies"),
            })

    # PRS tier shifts
    prs_prev = prev.get("prs_summary") or {}
    prs_curr = curr.get("prs_summary") or {}
    prs_changes: list[dict] = []
    for name in sorted(set(prs_prev) | set(prs_curr)):
        a = (prs_prev.get(name) or {}).get("tier")
        b = (prs_curr.get(name) or {}).get("tier")
        if a != b:
            prs_changes.append({"name": name, "from": a, "to": b,
                                "prev_pct": (prs_prev.get(name) or {}).get("percentile"),
                                "curr_pct": (prs_curr.get(name) or {}).get("percentile")})

    # PGx phenotype shifts
    pgx_prev = prev.get("pgx_summary") or {}
    pgx_curr = curr.get("pgx_summary") or {}
    pgx_changes: list[dict] = []
    for gene in sorted(set(pgx_prev) | set(pgx_curr)):
        a = (pgx_prev.get(gene) or {}).get("phenotype")
        b = (pgx_curr.get(gene) or {}).get("phenotype")
        if a != b:
            pgx_changes.append({"gene": gene, "from": a, "to": b})

    # Category set changes
    prev_cats = {v.get("category") for v in prev.get("variants", []) if v.get("category")}
    curr_cats = {v.get("category") for v in curr.get("variants", []) if v.get("category")}
    new_cats = sorted(curr_cats - prev_cats)
    dropped_cats = sorted(prev_cats - curr_cats)

    return {
        "prev_generated": prev.get("generated"),
        "curr_generated": curr.get("generated"),
        "prev_version": prev.get("report_version"),
        "curr_version": curr.get("report_version"),
        "prev_file_hash": prev.get("file_hash"),
        "curr_file_hash": curr.get("file_hash"),
        "prev_total": prev.get("total_matched"),
        "curr_total": curr.get("total_matched"),
        "new_variants": new_variants,
        "removed_variants": removed_variants,
        "changed_significance": changed_significance,
        "changed_category": changed_category,
        "changed_risk_copies": changed_risk_copies,
        "prs_changes": prs_changes,
        "pgx_changes": pgx_changes,
        "new_categories": new_cats,
        "dropped_categories": dropped_cats,
    }


# ── Pretty printer (terminal) ─────────────────────────────────────────────────
def render_diff_text(diff: dict) -> str:
    if "error" in diff:
        return diff["error"]

    out: list[str] = []
    def header(s: str) -> None:
        out.append(f"\n\033[1m{s}\033[0m")
        out.append("─" * len(s))

    out.append("\033[1mDNA Report Changelog\033[0m")
    out.append(f"\033[2mPrev: {diff['prev_generated']}  ({diff['prev_version']})  hash={diff['prev_file_hash']}\033[0m")
    out.append(f"\033[2mCurr: {diff['curr_generated']}  ({diff['curr_version']})  hash={diff['curr_file_hash']}\033[0m")
    out.append(f"\033[2mTotal variants:  {diff['prev_total']}  ->  {diff['curr_total']}\033[0m")

    if diff["prev_file_hash"] != diff["curr_file_hash"]:
        out.append("\n\033[33mNote: file hash differs — comparison is across two different chip files.\033[0m")

    if diff["new_variants"]:
        header(f"New variants matched ({len(diff['new_variants'])})")
        for v in diff["new_variants"][:40]:
            out.append(f"  + {v['rsid']:>12}  {v['gene']:14s}  [{v['category']}]  "
                       f"sig={v['significance']}  risk_copies={v['risk_copies']}")
        if len(diff["new_variants"]) > 40:
            out.append(f"  ... and {len(diff['new_variants'])-40} more")

    if diff["removed_variants"]:
        header(f"Variants no longer matched ({len(diff['removed_variants'])})")
        for v in diff["removed_variants"][:40]:
            out.append(f"  - {v['rsid']:>12}  {v['gene']:14s}  [{v['category']}]")

    if diff["changed_significance"]:
        header(f"Significance changes ({len(diff['changed_significance'])})")
        for c in diff["changed_significance"]:
            out.append(f"  ~ {c['rsid']:>12}  {c['gene']:14s}  {c['from']} -> {c['to']}")

    if diff["changed_risk_copies"]:
        header(f"Risk-copy changes ({len(diff['changed_risk_copies'])})")
        for c in diff["changed_risk_copies"]:
            out.append(f"  ~ {c['rsid']:>12}  {c['gene']:14s}  {c['from']} -> {c['to']}")

    if diff["changed_category"]:
        header(f"Category re-categorisations ({len(diff['changed_category'])})")
        for c in diff["changed_category"]:
            out.append(f"  ~ {c['rsid']:>12}  {c['gene']:14s}  {c['from']} -> {c['to']}")

    if diff["prs_changes"]:
        header(f"Polygenic risk score tier shifts ({len(diff['prs_changes'])})")
        for c in diff["prs_changes"]:
            a = c["from"] or "—"
            b = c["to"] or "—"
            pa = c.get("prev_pct")
            pb = c.get("curr_pct")
            extra = f" ({pa}th -> {pb}th)" if pa is not None and pb is not None else ""
            out.append(f"  ~ {c['name']:34s}  {a:14s} -> {b:14s}{extra}")

    if diff["pgx_changes"]:
        header(f"Pharmacogenomic phenotype changes ({len(diff['pgx_changes'])})")
        for c in diff["pgx_changes"]:
            out.append(f"  ~ {c['gene']:10s}  {c['from']} -> {c['to']}")

    if diff["new_categories"]:
        header("New report categories")
        for c in diff["new_categories"]:
            out.append(f"  + {c}")
    if diff["dropped_categories"]:
        header("Dropped categories")
        for c in diff["dropped_categories"]:
            out.append(f"  - {c}")

    if not any([
        diff["new_variants"], diff["removed_variants"], diff["changed_significance"],
        diff["changed_risk_copies"], diff["changed_category"], diff["prs_changes"],
        diff["pgx_changes"], diff["new_categories"], diff["dropped_categories"],
    ]):
        out.append("\n\033[32mNo material changes detected between runs.\033[0m")

    return "\n".join(out)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Diff two tier1_results.json files")
    ap.add_argument("prev", help="Previous tier1_results.json")
    ap.add_argument("curr", help="Current tier1_results.json")
    args = ap.parse_args()
    d = diff_runs(Path(args.prev), Path(args.curr))
    print(render_diff_text(d))
