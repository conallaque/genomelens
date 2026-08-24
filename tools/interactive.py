#!/usr/bin/env python3
"""
GenomeLens Interactive Q&A
==========================

An offline, terminal chat interface that lets you ask plain-language
questions about your own DNA analysis. Your tier1_results.json is loaded,
distilled into a readable context block, and used to ground answers from a
local Ollama model (qwen3:14b by default) — nothing leaves your machine.

Usage:
    python interactive.py tier1_results.json
    python interactive.py ~/dna-project/tier1_results.json --model qwen3:14b

In the chat:
    - Ask anything: "What's my CAD risk?", "Should I take supplements?"
    - 'context'  shows the genetic context the model is working from
    - 'history'  shows the conversation so far
    - 'exit' / 'quit' / Ctrl-D / Ctrl-C  leaves

Conversation history is kept in memory for the session only; it is never
written to disk.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit(
        "The 'requests' package is required. Install it with:\n"
        "    pip install requests"
    )


OLLAMA_HOST = "http://localhost:11434"
GENERATE_URL = f"{OLLAMA_HOST}/api/generate"
TAGS_URL = f"{OLLAMA_HOST}/api/tags"
DEFAULT_MODEL = "qwen3:14b"
REQUEST_TIMEOUT = 600  # seconds — a 14B model on CPU can be slow to first token

SYSTEM_PROMPT = (
    "You are a genetic health advisor. Answer questions about the user's "
    "personal DNA results, which are provided to you in the GENETIC CONTEXT "
    "block below. Ground every answer in their specific findings — refer to "
    "their actual genes, variants, risk tiers, and percentiles by name. "
    "Be honest about the limitations of consumer genetic testing: most traits "
    "are polygenic and strongly shaped by lifestyle and environment, "
    "genotyping arrays miss many variants, and a risk percentile is not a "
    "diagnosis. Do NOT provide medical advice or tell the user to start, stop, "
    "or change any medication — for any clinical decision, advise them to "
    "consult a doctor or a board-certified genetic counselor. Never invent "
    "findings that are not present in the GENETIC CONTEXT; if something was not "
    "tested or is unknown, say so plainly. Keep answers focused and practical "
    "(a few short paragraphs unless real depth is needed)."
)

SUGGESTED_QUESTIONS = [
    "What are my biggest health risks based on my genetics?",
    "What does my CAD risk mean?",
    "Should I worry about my ancestry results?",
    "Which medications should I be cautious about?",
    "What should I eat based on my genome?",
    "Should I take any supplements?",
    "What should I tell my doctor at my next visit?",
]


# ── Context building ──────────────────────────────────────────────────────────
def load_results(path: Path) -> dict:
    """Load and validate the tier1 results JSON.

    Raises FileNotFoundError or ValueError with a human-readable message.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find results file: {path}\n"
            "Generate it first (e.g. run analyze.py) or pass the correct path:\n"
            "    python interactive.py /path/to/tier1_results.json"
        )
    try:
        text = path.read_text()
    except OSError as e:
        raise ValueError(f"Could not read {path}: {e}") from e
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(
            f"{path} is not valid JSON (line {e.lineno}, column {e.colno}): "
            f"{e.msg}.\nThe file may be incomplete or corrupted — try "
            "regenerating it."
        ) from e
    if not isinstance(data, dict):
        raise ValueError(
            f"{path} did not contain a results object (got "
            f"{type(data).__name__}). Expected a tier1_results.json file."
        )
    return data


def _ordinal(n: int) -> str:
    """Return n with its ordinal suffix, e.g. 89 -> '89th', 63 -> '63rd'."""
    if 10 <= n % 100 <= 20:
        return f"{n}th"
    return f"{n}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th') }"


def build_context(data: dict) -> str:
    """Distill the results JSON into a compact, readable context block.

    Only includes what is actually present in the data — never fabricates
    fields. Real numbers (percentiles, tiers, genotypes) are kept verbatim so
    the model can quote them.
    """
    lines: list[str] = []

    # --- Identity / quality ---
    lines.append(f"Report version: {data.get('report_version', '—')}")
    lines.append(f"Data quality (QC) grade: {data.get('qc_grade', '—')}")
    lines.append(f"Curated variants matched: {data.get('total_matched', 0)}")

    apoe = data.get("apoe_genotype")
    lines.append(f"APOE genotype: {apoe if apoe else 'not determined on this chip'}")
    ydna = data.get("y_haplogroup")
    mt = data.get("mt_haplogroup")
    lines.append(f"Y-DNA haplogroup (paternal line): {ydna or 'not determined'}")
    lines.append(f"mtDNA haplogroup (maternal line): {mt or 'not determined'}")

    # --- Polygenic risk scores (the headline numbers) ---
    prs = data.get("prs_summary") or {}
    scored = {
        name: p for name, p in prs.items()
        if p and p.get("tier") is not None and p.get("percentile") is not None
    }
    if scored:
        lines.append("\nPOLYGENIC RISK SCORES (population percentile — higher = more genetic predisposition):")
        for name, p in scored.items():
            lines.append(f"  - {name}: {p['tier']} ({_ordinal(round(p['percentile']))} percentile)")
    not_scored = [name for name, p in prs.items() if not (p and p.get("tier"))]
    if not_scored:
        lines.append(
            "  (Not scored / insufficient data: " + ", ".join(not_scored) + ")"
        )

    # --- Pharmacogenomics ---
    pgx = data.get("pgx_summary") or {}
    if pgx:
        lines.append("\nPHARMACOGENOMICS (how genes affect drug metabolism/response):")
        for gene, p in pgx.items():
            phen = (p or {}).get("phenotype")
            if phen:
                lines.append(f"  - {gene}: {phen}")

    # --- Variant-level findings, grouped by category ---
    variants = data.get("variants") or []
    risk_carrying = [v for v in variants if v.get("risk_copies", 0) > 0]

    high = [v for v in risk_carrying if v.get("significance") == "high"]
    if high:
        lines.append(f"\nHIGH-SIGNIFICANCE FINDINGS WHERE YOU CARRY RISK ALLELE(S) ({len(high)} total, showing top 25):")
        for v in high[:25]:
            lines.append(
                f"  - {v.get('gene')} {v.get('variant_name')} "
                f"[{v.get('category')}] — genotype {v.get('my_genotype')}, "
                f"{v.get('risk_copies')} risk allele(s)"
            )

    # Highlight the lifestyle-actionable categories the user asks about most:
    # nutrition, athletic performance / exercise, longevity, ancestry.
    spotlight = {
        "Diet & Nutrition": "NUTRITION-RELATED VARIANTS",
        "Athletic Performance": "EXERCISE / ATHLETIC PERFORMANCE VARIANTS",
        "Longevity": "LONGEVITY VARIANTS",
        "Ancestry Informative Markers": "ANCESTRY-INFORMATIVE MARKERS",
    }
    for category, heading in spotlight.items():
        in_cat = [v for v in risk_carrying if v.get("category") == category]
        if in_cat:
            lines.append(f"\n{heading} (you carry the noted allele(s)):")
            for v in in_cat[:10]:
                rec = (v.get("recommendation") or "").strip()
                rec = (rec[:140] + "…") if len(rec) > 140 else rec
                line = (
                    f"  - {v.get('gene')} {v.get('variant_name')} "
                    f"(genotype {v.get('my_genotype')})"
                )
                if rec:
                    line += f": {rec}"
                lines.append(line)

    # Category coverage summary so the model knows the breadth available.
    from collections import Counter
    cat_counts = Counter(v.get("category") for v in risk_carrying)
    if cat_counts:
        summary = ", ".join(
            f"{cat} ({n})" for cat, n in cat_counts.most_common()
        )
        lines.append(
            "\nAll categories with risk-carrying variants (count): " + summary
        )

    return "\n".join(lines)


# ── Terminal rendering ────────────────────────────────────────────────────────
THINK_RE = re.compile(r"<think>.*?</think>", re.DOTALL)


def strip_think(text: str) -> str:
    """Remove qwen3 <think>…</think> reasoning blocks."""
    return THINK_RE.sub("", text).strip()


def render(text: str) -> str:
    """Light ANSI formatting for an assistant response."""
    text = re.sub(r"\*\*(.+?)\*\*", "\033[1m\\1\033[0m", text)          # bold
    text = re.sub(r"\*([^*\n]+)\*", "\033[3m\\1\033[0m", text)           # italic
    out = []
    for line in text.split("\n"):
        stripped = line.lstrip()
        if stripped.startswith(("- ", "* ", "• ")):
            indent = len(line) - len(stripped)
            out.append(" " * indent + "\033[36m•\033[0m " + stripped[2:])
        else:
            out.append(line)
    return "\n".join(out)


# ── Ollama plumbing ───────────────────────────────────────────────────────────
def check_ollama(model: str) -> None:
    """Verify Ollama is reachable and the model is available. Exits on failure."""
    try:
        resp = requests.get(TAGS_URL, timeout=5)
        resp.raise_for_status()
    except requests.exceptions.ConnectionError:
        sys.exit(
            f"\033[31mCannot reach Ollama at {OLLAMA_HOST}.\033[0m\n"
            "Start it in another terminal with:\n"
            "    ollama serve\n"
            f"…then make sure the model is pulled:\n"
            f"    ollama pull {model}"
        )
    except requests.exceptions.RequestException as e:
        sys.exit(f"\033[31mOllama health check failed: {e}\033[0m")

    try:
        tags = resp.json().get("models", [])
        names = {m.get("name", "") for m in tags}
    except (ValueError, AttributeError):
        names = set()
    # Match with or without an explicit ':latest' tag.
    if names and model not in names and f"{model}:latest" not in names:
        avail = ", ".join(sorted(n for n in names if n)) or "none"
        print(
            f"\033[33mWarning: model '{model}' is not in Ollama's list "
            f"(available: {avail}).\n"
            f"Trying anyway — pull it with 'ollama pull {model}' if this "
            f"fails.\033[0m"
        )


def build_prompt(context: str, history: list[tuple[str, str]], question: str) -> str:
    """Assemble the /api/generate prompt from context + clean history + question."""
    parts = ["--- GENETIC CONTEXT (the user's DNA results) ---", context, ""]
    if history:
        parts.append("--- CONVERSATION SO FAR ---")
        for user_q, assistant_a in history:
            parts.append(f"User: {user_q}")
            parts.append(f"Advisor: {assistant_a}")
        parts.append("")
    parts.append(f"User: {question}")
    parts.append("Advisor:")
    return "\n".join(parts)


def ask_ollama(model: str, prompt: str) -> str:
    """Send one prompt to /api/generate and return the cleaned response text.

    Raises requests exceptions; the caller handles them.
    """
    payload = {
        "model": model,
        "prompt": prompt,
        "system": SYSTEM_PROMPT,
        "stream": False,
        "options": {"temperature": 0.4, "num_ctx": 8192},
    }
    resp = requests.post(GENERATE_URL, json=payload, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    data = resp.json()
    if data.get("error"):
        raise RuntimeError(data["error"])
    return strip_think(data.get("response", "")) or "(no response)"


# ── REPL ──────────────────────────────────────────────────────────────────────
def chat_loop(model: str, context: str) -> None:
    history: list[tuple[str, str]] = []  # (question, cleaned_answer) in memory

    bar = "\033[1m" + "─" * 64 + "\033[0m"
    print(bar)
    print("\033[1mGenomeLens — Interactive DNA Q&A\033[0m")
    print(f"Model: \033[36m{model}\033[0m  ·  fully offline via local Ollama")
    print("Ask about your results. Commands: 'help', 'context', 'history', 'exit'.")
    print(bar)

    while True:
        try:
            question = input("\n\033[1myou> \033[0m").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye — take care of yourself.")
            return

        if not question:
            continue
        cmd = question.lower()
        if cmd in ("exit", "quit", ":q"):
            print("Goodbye — take care of yourself.")
            return
        if cmd in ("help", "?"):
            print("\033[2mTry asking:\033[0m")
            for q in SUGGESTED_QUESTIONS:
                print(f"  · {q}")
            continue
        if cmd == "context":
            print("\033[2m" + context + "\033[0m")
            continue
        if cmd == "history":
            if not history:
                print("\033[2m(no conversation yet)\033[0m")
            for i, (q, a) in enumerate(history, 1):
                print(f"\033[2m{i}. you:\033[0m {q}")
                print(f"\033[2m   advisor:\033[0m {a[:200]}{'…' if len(a) > 200 else ''}")
            continue

        prompt = build_prompt(context, history, question)
        print("\033[2m  thinking…\033[0m", end="", flush=True)
        t0 = time.time()
        try:
            answer = ask_ollama(model, prompt)
        except requests.exceptions.ConnectionError:
            print("\r" + " " * 16 + "\r", end="")
            print(
                "\033[31m  Lost connection to Ollama. Is 'ollama serve' still "
                "running?\033[0m"
            )
            continue
        except requests.exceptions.Timeout:
            print("\r" + " " * 16 + "\r", end="")
            print(
                f"\033[31m  Ollama did not respond within {REQUEST_TIMEOUT}s. "
                "The model may be loading or the question too large — try "
                "again or ask something shorter.\033[0m"
            )
            continue
        except (requests.exceptions.RequestException, RuntimeError, ValueError) as e:
            print("\r" + " " * 16 + "\r", end="")
            print(f"\033[31m  Error talking to Ollama: {e}\033[0m")
            continue

        print("\r" + " " * 16 + "\r", end="")  # erase 'thinking…'
        elapsed = time.time() - t0
        print("\033[2m─── Advisor ───\033[0m")
        print(render(answer))
        print(f"\033[2m  ({elapsed:.1f}s)\033[0m")

        history.append((question, answer))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Interactive, offline Q&A about your DNA results."
    )
    parser.add_argument(
        "results",
        nargs="?",
        default="tier1_results.json",
        help="Path to tier1_results.json (default: ./tier1_results.json)",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Ollama model to use (default: {DEFAULT_MODEL})",
    )
    args = parser.parse_args(argv)

    path = Path(args.results).expanduser()
    try:
        data = load_results(path)
    except (FileNotFoundError, ValueError) as e:
        print(f"\033[31m{e}\033[0m", file=sys.stderr)
        return 1

    context = build_context(data)
    check_ollama(args.model)
    chat_loop(args.model, context)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
