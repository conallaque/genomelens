#!/usr/bin/env python3
"""Verify this public release. One command, public repository content only.

    python tools/verify_public_release.py

No network, no credentials, no GenomeLens engine, no private paths. A fresh
clone can run it.

SINGLE SOURCE OF TRUTH. The checks themselves live in `tests/public/`. This
script reads the published artifacts to print the headline figures, then runs
that suite and reports its result. It deliberately does not reimplement the
assertions: two copies of the same logic drift, and the copy nobody runs is the
one that rots.

What this does NOT do: reproduce the analysis. Regenerating these results
requires the production engine, which is not published. This verifies that what
is published here is internally consistent, correctly scoped, and free of the
disclosure classes the release excludes.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "artifacts" / "public-benchmarks"
SAMPLES = ("HG002", "HG003", "HG004")
WIDTH = 52


def dots(label: str, value: str) -> str:
    pad = max(1, WIDTH - len(label) - len(value))
    return f"  {label} {'.' * pad} {value}"


def main() -> int:
    print("GenomeLens public release verification\n")

    # ── artifacts present ────────────────────────────────────────────────
    required = ("summary.json", "report-public.html", "report-public.pdf", "preview.png")
    print("Artifacts")
    missing: list[str] = []
    for s in SAMPLES:
        absent = [n for n in required if not (BENCH / s / n).is_file()]
        missing += [f"{s}/{n}" for n in absent]
        print(dots(s, "PASS" if not absent else "FAIL"))
    if missing:
        print(f"\nRESULT: FAIL — missing artifacts: {', '.join(missing)}")
        return 1

    summaries = {s: json.loads((BENCH / s / "summary.json").read_text()) for s in SAMPLES}

    # ── benchmark ────────────────────────────────────────────────────────
    print("\nGIAB benchmark")
    total = mism = norm = 0
    for s, d in summaries.items():
        v = d["validation"]
        total += v["explicit_truth_calls"]
        mism += v["mismatches"]
        norm += v["normalization_failures"]
        print(dots(s, f"{v['concordant']}/{v['explicit_truth_calls']}"))
    print(dots("Aggregate", f"{total}/526"))
    print(dots("Mismatches", str(mism)))
    print(dots("Normalization failures", str(norm)))
    scope = {(d["dataset"]["benchmark_release"], d["dataset"]["reference_build"],
              d["dataset"]["benchmark_scope"]) for d in summaries.values()}
    print(dots("Scope", " · ".join(sorted(scope)[0])))

    # ── economics ────────────────────────────────────────────────────────
    print("\nEconomics")
    for s, d in summaries.items():
        e = d["economics"]
        icer = e["icer"] if isinstance(e["icer"], str) else f"{e['icer']:,}"
        print(dots(f"{s}  reference-case / canonical / ICER",
                   f"${e['reference_case_nmb_usd']:,} / "
                   f"${e['canonical_expected_nmb_usd']:,} / {icer}"))

    # ── the suite ────────────────────────────────────────────────────────
    print("\nPublic validation suite")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/public", "-q", "--no-header",
         "-p", "no:cacheprovider"],
        cwd=ROOT, capture_output=True, text=True)
    tail = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()]
    summary_line = tail[-1] if tail else "(no output)"

    if "no tests ran" in summary_line or not tail:
        print(dots("suite executed", "FAIL"))
        print("\nRESULT: FAIL — the suite did not run. An exit code from a suite\n"
              "that never executed is not evidence of anything.")
        print(proc.stdout[-1500:] or proc.stderr[-1500:])
        return 1

    print(dots("result", summary_line))
    if " skipped" in summary_line:
        # Report WHY, rather than implying a missing dependency. A skip whose
        # reason is unstated is indistinguishable from a check that was quietly
        # dropped.
        reasons = subprocess.run(
            [sys.executable, "-m", "pytest", "tests/public", "-q", "--no-header",
             "-rs", "-p", "no:cacheprovider"],
            cwd=ROOT, capture_output=True, text=True).stdout
        for ln in reasons.splitlines():
            if ln.startswith("SKIPPED"):
                print(f"  skipped: {ln.split('] ', 1)[-1].strip()}")

    if proc.returncode != 0:
        print("\nFailing checks:")
        for ln in proc.stdout.splitlines():
            if ln.startswith(("FAILED", "ERROR")):
                print(f"  {ln}")
        print("\nRESULT: FAIL")
        return proc.returncode

    print("\nRESULT: PASS")
    print("\nThis verifies the published artifacts. It does not re-run the analysis —\n"
          "that requires the production engine, which is not published. See\n"
          "docs/TESTING.md for what the public and private suites each cover.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
