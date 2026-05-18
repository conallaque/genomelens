"""
Interactive Q&A Chat Mode
=========================

After analyze.py runs the standard pipeline, --chat drops you into a REPL
where you can ask plain-language questions about your genetic findings.
Each question is sent to local Ollama with a compact context block
summarising your APOE, Y/mtDNA, PRS, PGx, traits, compound interactions,
and carrier findings.

Type 'exit', 'quit', or Ctrl-D to leave.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional


SCRIPT_DIR = Path(__file__).parent
DEFAULT_TIER1 = SCRIPT_DIR / "tier1_results.json"


def _summarise_context(
    tier1_path: Path,
    pgx_result: Optional[Dict] = None,
    prs_result: Optional[Dict] = None,
    traits_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
) -> str:
    """Build a compact context block describing the user's genetic state."""
    lines: List[str] = []

    if tier1_path.exists():
        try:
            data = json.loads(tier1_path.read_text())
            lines.append(f"REPORT VERSION: {data.get('report_version','—')}")
            apoe = data.get("apoe_genotype")
            ydna = data.get("y_haplogroup")
            mt = data.get("mt_haplogroup")
            grade = data.get("qc_grade")
            lines.append(f"APOE genotype: {apoe or 'not determined'}")
            lines.append(f"Y-DNA haplogroup: {ydna or 'not determined'}")
            lines.append(f"mtDNA haplogroup: {mt or 'not determined'}")
            lines.append(f"Data QC grade: {grade or '—'}")
            total = data.get("total_matched", 0)
            lines.append(f"Curated variants matched: {total}")

            # Risk-carrying variants by category
            variants = data.get("variants", [])
            risk_carrying = [v for v in variants if v.get("risk_copies", 0) > 0]
            high_sig = [v for v in risk_carrying if v.get("significance") == "high"]
            mod_sig = [v for v in risk_carrying if v.get("significance") == "moderate"]
            if high_sig:
                lines.append("\nHIGH-significance risk findings:")
                for v in high_sig[:20]:
                    lines.append(
                        f"  - {v['gene']} {v['variant_name']} "
                        f"({v['risk_copies']} risk allele(s), category: {v['category']})"
                    )
            if mod_sig:
                lines.append(f"\nModerate-significance risk findings ({len(mod_sig)}):")
                for v in mod_sig[:20]:
                    lines.append(
                        f"  - {v['gene']} {v['variant_name']} "
                        f"({v['risk_copies']} risk allele(s), category: {v['category']})"
                    )

            # PRS summary
            prs_summary = data.get("prs_summary") or {}
            if prs_summary:
                lines.append("\nPolygenic Risk Score tiers:")
                for name, p in prs_summary.items():
                    tier = (p or {}).get("tier")
                    pct = (p or {}).get("percentile")
                    if tier:
                        lines.append(f"  - {name}: {tier} ({pct}th percentile)")

            # PGx summary
            pgx_summary = data.get("pgx_summary") or {}
            if pgx_summary:
                lines.append("\nPharmacogenomic phenotypes:")
                for gene, p in pgx_summary.items():
                    phen = (p or {}).get("phenotype")
                    if phen:
                        lines.append(f"  - {gene}: {phen}")
        except Exception as e:
            lines.append(f"(Could not load tier1_results.json: {e})")
    else:
        lines.append("(No tier1_results.json found — context is limited.)")

    # Optional live results from this session
    if traits_result and traits_result.get("predictions"):
        lines.append("\nTrait predictions:")
        for t in traits_result["predictions"]:
            if t.get("confidence") in ("high", "moderate"):
                lines.append(f"  - {t.get('trait')}: {t.get('result')}")

    if carrier_result:
        n_aff = carrier_result.get("n_affected", 0)
        n_car = carrier_result.get("n_carriers", 0)
        if n_aff or n_car:
            lines.append(f"\nCarrier status: {n_aff} affected, {n_car} carrier(s)")

    if interactions_result and interactions_result.get("findings"):
        lines.append("\nCompound / variant-interaction findings:")
        for f in interactions_result["findings"][:8]:
            lines.append(f"  - [{f['severity'].upper()}] {f['title']}")

    if ancestry_result and ancestry_result.get("proportions"):
        props = ancestry_result["proportions"]
        sorted_props = sorted(props.items(), key=lambda x: -x[1])
        top = ", ".join(f"{sp} {p*100:.0f}%" for sp, p in sorted_props[:3])
        lines.append(f"\nEstimated ancestry (rough): {top}")

    return "\n".join(lines)


# ── Pretty terminal rendering ─────────────────────────────────────────────────
def _render_response(text: str) -> str:
    """Light terminal formatting for Ollama responses."""
    # Strip <think>...</think> blocks
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
    # Bold via ANSI
    text = re.sub(r"\*\*(.+?)\*\*", "\033[1m\\1\033[0m", text)
    # Italic via ANSI
    text = re.sub(r"\*([^*\n]+)\*", "\033[3m\\1\033[0m", text)
    # Bullet aesthetics
    lines = []
    for line in text.split("\n"):
        stripped = line.lstrip()
        if stripped.startswith(("- ", "* ", "• ")):
            indent = len(line) - len(stripped)
            lines.append(" " * indent + "\033[36m•\033[0m " + stripped[2:])
        else:
            lines.append(line)
    return "\n".join(lines)


def _print_response(text: str) -> None:
    print()
    print("\033[2m─── Assistant ───\033[0m")
    print(_render_response(text))
    print()


def _print_user_echo(text: str) -> None:
    print(f"\033[2m> {text}\033[0m")


# ── REPL ──────────────────────────────────────────────────────────────────────
SYSTEM_PROMPT = (
    "You are a knowledgeable genetic counselor with deep expertise in "
    "pharmacogenomics, polygenic risk, and preventive genomics. You are "
    "having a conversation with a person about their personal DNA analysis. "
    "Use the provided GENETIC CONTEXT to ground every answer in their "
    "specific findings. Be concrete, practical, and personalised — refer "
    "to their specific genes, variants, and risk tiers by name. Keep "
    "answers focused (3–6 short paragraphs unless the question requires "
    "depth). Always remind them this is educational and that medication "
    "and clinical decisions belong with a physician or board-certified "
    "genetic counselor. Do NOT invent findings not present in the context."
)


SUGGESTED_QUESTIONS = [
    "What are my biggest health risks based on my genetics?",
    "What should I eat based on my genome?",
    "Which medications should I be cautious about?",
    "How would I respond to Vyvanse, Adderall, or other ADHD medications?",
    "What should I tell my doctor at my next visit?",
    "How does my body respond to caffeine and alcohol?",
    "What lifestyle changes would have the biggest impact for me?",
    "Should I see a genetic counselor?",
]


def run_chat(
    model: str = "qwen3:14b",
    tier1_path: Optional[Path] = None,
    pgx_result: Optional[Dict] = None,
    prs_result: Optional[Dict] = None,
    traits_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
) -> None:
    """Start an interactive chat session. Returns when the user exits."""
    import requests

    tier1_path = tier1_path or DEFAULT_TIER1
    context = _summarise_context(
        tier1_path,
        pgx_result=pgx_result, prs_result=prs_result,
        traits_result=traits_result, carrier_result=carrier_result,
        interactions_result=interactions_result, ancestry_result=ancestry_result,
    )

    print("\033[1m" + "─" * 64 + "\033[0m")
    print("\033[1mDNA Chat Mode\033[0m — ask questions about your analysis.")
    print(f"Model: \033[36m{model}\033[0m (local Ollama)")
    print("Type 'exit', 'quit', or Ctrl-D to leave. Type 'help' for suggestions.")
    print("\033[1m" + "─" * 64 + "\033[0m")

    history: List[Dict[str, str]] = [
        {"role": "system", "content": SYSTEM_PROMPT + "\n\n--- GENETIC CONTEXT ---\n" + context},
    ]

    while True:
        try:
            user_in = input("\033[1myou> \033[0m").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye.")
            break
        if not user_in:
            continue
        if user_in.lower() in ("exit", "quit", ":q"):
            print("Goodbye.")
            break
        if user_in.lower() in ("help", "?", "suggestions"):
            print("\033[2mSuggested questions:\033[0m")
            for q in SUGGESTED_QUESTIONS:
                print(f"  · {q}")
            print()
            continue
        if user_in.lower() == "context":
            print(context)
            continue

        history.append({"role": "user", "content": user_in})

        # Stream response from Ollama
        try:
            payload = {
                "model": model,
                "messages": history,
                "stream": False,
                "options": {"temperature": 0.4, "num_ctx": 8192},
            }
            t0 = time.time()
            print("\033[2m  thinking...\033[0m", end="", flush=True)
            resp = requests.post(
                "http://localhost:11434/api/chat",
                json=payload, timeout=420,
            )
            resp.raise_for_status()
            content = resp.json()["message"]["content"]
            # Erase 'thinking...' line
            print("\r" + " " * 20 + "\r", end="", flush=True)
            elapsed = time.time() - t0
            _print_response(content)
            print(f"\033[2m  ({elapsed:.1f}s)\033[0m\n")
            history.append({"role": "assistant", "content": content})
        except requests.exceptions.ConnectionError:
            print("\n\033[31m  Cannot connect to Ollama. Start it: ollama serve\033[0m")
            break
        except Exception as e:
            print(f"\n\033[31m  Error: {e}\033[0m")
            history.pop()  # don't keep the failed exchange


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="DNA chat mode")
    ap.add_argument("--model", default="qwen3:14b")
    ap.add_argument("--tier1", default=str(DEFAULT_TIER1),
                    help="Path to tier1_results.json")
    args = ap.parse_args()
    run_chat(model=args.model, tier1_path=Path(args.tier1))
