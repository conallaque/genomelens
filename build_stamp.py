"""Build provenance for generated artefacts.

WHY THIS EXISTS. A report or PDF used to carry nothing that said which code
produced it. That makes a stale artefact indistinguishable from a current one by
inspection: a reviewer searching a PDF for figures that "should" have changed
cannot tell whether the fix is missing or the file is old, and neither can the
person who generated it. A passing test suite proves the code is right; it says
nothing about the bytes in a particular file.

So every generated artefact now embeds a marker naming the commit it came from,
whether the working tree was clean at the time, and when it was written. The
marker is a single machine-findable line so it can be grepped out of HTML, out
of extracted PDF text, or read from PDF metadata:

    GENOMELENS-BUILD: <short-sha>[+dirty] <iso-8601>

``+dirty`` matters as much as the hash. An artefact built from a working tree
with uncommitted edits reflects code that exists in no commit, which is exactly
how a report and its git history drift apart without either being wrong.
"""

from __future__ import annotations

import datetime as _dt
import os
import subprocess

MARKER_PREFIX = "GENOMELENS-BUILD:"
_REPO = os.path.dirname(os.path.abspath(__file__))


def _git(*args: str) -> str:
    try:
        out = subprocess.run(("git", "-C", _REPO, *args), capture_output=True,
                             text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def build_stamp() -> dict:
    """Commit, cleanliness and timestamp for the code generating this artefact.

    Degrades to ``commit="unknown"`` outside a git checkout rather than raising:
    an artefact with no provenance is worse than one with partial provenance,
    but it is not worth failing a report over.
    """
    sha = _git("rev-parse", "--short", "HEAD") or "unknown"
    dirty = bool(_git("status", "--porcelain"))
    return {
        "commit": sha,
        "dirty": dirty,
        "generated_at": _dt.datetime.now().astimezone().isoformat(
            timespec="seconds"),
        "marker": f"{MARKER_PREFIX} {sha}{'+dirty' if dirty else ''} "
                  f"{_dt.datetime.now().astimezone().isoformat(timespec='seconds')}",
    }


def marker_in(text: str) -> str:
    """The build marker found in an artefact's text, or "" if it carries none."""
    for line in text.splitlines():
        i = line.find(MARKER_PREFIX)
        if i >= 0:
            return line[i:].strip()
    return ""


def commit_of(text: str) -> str:
    """The short commit hash an artefact records, or "" if unmarked."""
    m = marker_in(text)
    if not m:
        return ""
    parts = m.split()
    return parts[1].split("+")[0] if len(parts) > 1 else ""


def is_dirty_build(text: str) -> bool:
    return "+dirty" in marker_in(text)
