#!/usr/bin/env python3
"""Verify this public release. One command, public repository content only.

    python tools/verify_public_release.py

No network, no credentials, no GenomeLens engine, no private paths. A fresh
clone can run it.

SINGLE SOURCE OF TRUTH. The checks themselves live in `tests/public/`. This
script reads the published files to print what they contain, then runs that
suite and reports its result. It deliberately does not reimplement the
assertions: two copies of the same logic drift, and the copy nobody runs is the
one that rots.

What this does NOT do: reproduce the analysis. Regenerating these results
requires the production engine, which is not published. This verifies that what
is published here is internally consistent, correctly scoped, and free of the
generic disclosure classes it screens for (local paths, credentials, commercial
language, numeric value-of-information examples). It does not screen for exact
withheld values: a list of those in public code would publish them, so that
check belongs to a private release gate this repository does not contain.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ("HG002", "HG003", "HG004")
WIDTH = 52


def dots(label: str, value: str) -> str:
    pad = max(1, WIDTH - len(label) - len(value))
    return f"  {label} {'.' * pad} {value}"


def main() -> int:
    print("GenomeLens public release verification\n")

    # -- published reports present ---------------------------------------
    print("Published reports")
    demo = ROOT / "demonstrations" / "v1.11"
    reports = sorted(demo.glob("GenomeLens_public_demonstration_*.pdf"))
    dossiers = sorted((demo / "dossiers").glob("*.pdf")) if (demo / "dossiers").is_dir() else []
    featured = sorted((demo / "featured").glob("*.pdf")) if (demo / "featured").is_dir() else []
    bench = sorted((ROOT / "benchmark" / "v1").glob("*.pdf"))
    print(dots("v1.11 demonstration reports", str(len(reports))))
    print(dots("featured six-page report", str(len(featured))))
    print(dots("HEOR technical dossier", str(len(dossiers))))
    print(dots("Benchmark v1 documents", str(len(bench))))
    if not (reports and bench):
        print("\nRESULT: FAIL - the published reports or the benchmark documents are missing")
        return 1

    # -- historical GIAB concordance, as the README states it ------------
    print("\nHistorical GIAB concordance (an earlier engine build; not v1.11)")
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    total = 0
    for s in SAMPLES:
        row = re.search(rf"\|\s*{s}\s*\|[^\n]*", text[text.index("## Earlier GIAB trio concordance"):])
        m = re.search(r"\b(\d+)/(\d+)\b", row.group(0)) if row else None
        if not m:
            print(dots(s, "FAIL"))
            print(f"\nRESULT: FAIL - the README no longer states {s}'s concordance")
            return 1
        total += int(m.group(2))
        print(dots(s, f"{m.group(1)}/{m.group(2)}"))
    print(dots("Aggregate", f"{total}/526"))

    # -- the suite -------------------------------------------------------
    print("\nPublic validation suite")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/public", "-q", "--no-header",
         "-p", "no:cacheprovider"],
        cwd=ROOT, capture_output=True, text=True)
    tail = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()]
    summary_line = tail[-1] if tail else "(no output)"

    if "no tests ran" in summary_line or not tail:
        print(dots("suite executed", "FAIL"))
        print("\nRESULT: FAIL - the suite did not run. An exit code from a suite\n"
              "that never executed is not evidence of anything.")
        print(proc.stdout[-1500:] or proc.stderr[-1500:])
        return 1

    print(dots("result", summary_line))
    print(dots("exact withheld-value screening", "NOT PERFORMED HERE"))
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
    print("\nThis verifies the published files. It does not re-run the analysis -\n"
          "that requires the production engine, which is not published - and it\n"
          "does not check for exact withheld values, which a private release gate\n"
          "screens before publication. See docs/TESTING.md for what the public and\n"
          "private suites each cover.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
