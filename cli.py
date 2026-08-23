"""
Command-line interface for the DNA analysis pipeline.

V7 extracted the argparse parser into this module so it can be tested
independently and so the `dna-analyze` console-script entry point declared in
``pyproject.toml`` can point at a small, focused file rather than the 4 700-
line ``analyze.py``.

The full ``main()`` orchestration still lives in ``analyze.py`` for now; the
follow-up step in the V7 plan is to extract:

  * ``pipeline.py``   — the per-stage orchestration body of ``main()``
  * ``renderers/``    — the 24 ``build_*_html`` functions
  * ``tier2_ai.py``   — the Ollama integration and prompt assembly

Today this module's job is narrow: build the parser, then hand control to
``analyze.main()``. That keeps the import graph one-way (``cli`` → ``analyze``,
never the reverse) and means existing ``python analyze.py …`` invocations
keep working through the back-compat shim at the bottom of ``analyze.py``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

# Importing ``analyze`` at function scope (inside ``main``) avoids paying the
# multi-second import cost when the parser is constructed for ``--help`` or for
# unit tests that only want to inspect the flags.

DEFAULT_OLLAMA_MODEL = "qwen3:14b"


def build_parser() -> argparse.ArgumentParser:
    """Construct the full argparse parser for the DNA analysis CLI.

    Returns a configured ``ArgumentParser`` ready for ``parse_args()`` — but
    deliberately does **not** call it. This lets tests assert on flag
    definitions, build help text, and inspect defaults without invoking the
    pipeline.
    """
    parser = argparse.ArgumentParser(
        prog="dna-analyze",
        description=(
            "Local, privacy-first DNA analysis tool — v7.0.0 (consolidation)\n\n"
            "Reads a raw chip DNA file (23andMe / AncestryDNA / TellmeGen etc.),\n"
            "produces a comprehensive HTML report including curated variants,\n"
            "polygenic risk scores, CPIC pharmacogenomic phenotypes, compound\n"
            "heterozygosity, carrier status, trait predictions, ancestry, and\n"
            "optional local-AI interpretation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  dna-analyze ~/Downloads/genome.csv\n"
            "  dna-analyze genome.txt --impute               # imputes via Beagle\n"
            "  dna-analyze genome.txt --pdf                  # also emit PDF\n"
            "  dna-analyze genome.txt --bloodwork labs.json  # compare labs vs predictions\n"
            "  dna-analyze genome.txt --fhir                 # emit clinical EHR bundle\n"
            "  dna-analyze genome.txt --carrier-report\n"
            "  dna-analyze genome.txt --chat\n"
            "  dna-analyze genome.txt --compare prev_results.json\n\n"
            "First-time setup (downloads ~3 GB once):\n"
            "  python setup.py --all\n"
        ),
    )

    # ── Positional / core ────────────────────────────────────────────────
    parser.add_argument(
        "dna_file",
        nargs="?",
        default=None,
        help=("Path to your raw DNA file — a consumer chip export (CSV/TXT from "
              "23andMe, AncestryDNA, TellmeGen, …) OR a whole-genome/exome VCF "
              "(.vcf / .vcf.gz). Optional when --retry-failed is used."),
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output HTML path (default: report.html in project directory)",
    )
    parser.add_argument(
        "--no-ai",
        action="store_true",
        help="Skip Tier 2 AI analysis (faster, works without Ollama)",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_OLLAMA_MODEL,
        help=f"Ollama model name (default: {DEFAULT_OLLAMA_MODEL})",
    )

    # ── v3 / pipeline-stage flags ─────────────────────────────────────────
    parser.add_argument(
        "--impute",
        action="store_true",
        help=("Run Beagle 5.4 imputation against the 1000 Genomes Phase 3 panel "
              "before analysis. Requires `python setup.py --beagle` once. "
              "Adds ~30-90 min on first run; subsequent runs use cache."),
    )
    parser.add_argument(
        "--impute-min-r2",
        type=float, default=0.3,
        help="Minimum DR2 (R²) to keep an imputed variant (default 0.3).",
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Also emit a paginated PDF report (requires weasyprint).",
    )
    parser.add_argument(
        "--medications",
        default=None,
        help=("Comma-separated list of medication names (brand or generic) to "
              "cross-reference against your PGx phenotypes. Adds a focused "
              "Medication Review section."),
    )
    parser.add_argument(
        "--carrier-report",
        action="store_true",
        help=("Generate a standalone Carrier Status / Family Planning HTML "
              "document organised by severity, alongside the main report."),
    )
    parser.add_argument(
        "--chat",
        action="store_true",
        help=("After analysis completes, drop into an interactive Q&A REPL "
              "backed by local Ollama with your full genetic context loaded."),
    )
    parser.add_argument(
        "--compare",
        default=None,
        help=("Path to a previous tier1_results.json. Prints a structured "
              "changelog of what changed between runs."),
    )
    parser.add_argument(
        "--latitude",
        type=float,
        default=40.0,
        help=("Your approximate latitude in degrees (e.g. 43 for Michigan, 51 "
              "for London). Tunes the Environmental Optimization section's "
              "vitamin-D seasonality. Default 40.0 (temperate N hemisphere)."),
    )
    parser.add_argument(
        "--age",
        type=int,
        default=None,
        help=("Your age in years. Highlights the current decade in the "
              "Life-Stage Playbook. If omitted, age is taken from --bloodwork "
              "labs if present; otherwise all decades are shown un-highlighted."),
    )
    parser.add_argument(
        "--no-module-ai",
        action="store_true",
        help=("Skip the per-module AI interpretation pass ('AI on all tiers'). "
              "The executive summary and cross-category synthesis still run. "
              "Use for a faster run — each module interpretation is a separate "
              "local LLM call."),
    )
    parser.add_argument(
        "--compare-genome",
        default=None,
        metavar="PATH",
        help=("Path to a SECOND person's raw DNA file. Writes a standalone "
              "genome_comparison.html with a KING-robust relationship estimate "
              "(kinship + IBS0), genotype concordance, and shared recessive "
              "carrier risk for a couple. The second genome is never written "
              "to the repo and never committed."),
    )

    # ── v5 flags ──────────────────────────────────────────────────────────
    parser.add_argument(
        "--narrative",
        action="store_true",
        help=("After analysis, send a comprehensive summary to local Ollama to "
              "generate a warmly-written 'genetic counsellor explaining results' "
              "narrative_report.html. Requires Ollama running."),
    )
    # Opt-OUT, not opt-in. Every other companion page (supplements, exercise,
    # nutrition, economic analysis, longevity, master plan) is written on every
    # run; the emergency card was the only one gated behind a flag, so a normal
    # run produced no card and it looked like the feature had been dropped.
    parser.add_argument(
        "--emergency-card",
        action="store_true",
        default=True,
        help=("Write the one-page emergency_card.html — clinically actionable "
              "findings only (drug hypersensitivities, clotting disorders, PGx "
              "extremes) for emergency clinicians. On by default."),
    )
    parser.add_argument(
        "--no-emergency-card",
        dest="emergency_card",
        action="store_false",
        help="Skip the emergency card.",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help=("Re-run AI only for categories listed in failed_categories.json "
              "(written automatically when an earlier run hit Ollama timeouts) "
              "and patch the fresh interpretations into the existing report.html "
              "without regenerating the rest. dna_file is not required when "
              "this flag is used."),
    )

    # ── v6 flags ──────────────────────────────────────────────────────────
    parser.add_argument(
        "--bloodwork",
        default=None,
        help=("Path to a JSON file of measured lab values (LDL, testosterone, "
              "vitamin D, CRP, glucose, ferritin, etc.). Generates a "
              "bloodwork.html comparing genetic predictions to actual labs."),
    )
    parser.add_argument(
        "--fhir",
        action="store_true",
        help=("Export clinically-validated findings (PGx, carrier, HLA, APOE) "
              "as an HL7 FHIR R4 Bundle JSON file (fhir_bundle.json) suitable "
              "for ingestion into a clinical EHR."),
    )

    # ── Phase-3 / WGS flags ───────────────────────────────────────────────
    parser.add_argument(
        "--assume-build",
        choices=["grch37", "grch38"],
        default=None,
        help=("Force the reference build instead of auto-detecting it. Useful "
              "for whole-genome VCFs whose variants lack rsIDs (e.g. GIAB "
              "benchmark callsets), where probe-based detection cannot resolve "
              "the build. Takes precedence over auto-detection."),
    )
    parser.add_argument(
        "--commercial-safe",
        action="store_true",
        help=("Restrict the Phase-3 novel-variant predictors to commercially "
              "licensed sources only (AlphaMissense + gnomAD). Disables REVEL, "
              "SpliceAI, and CADD, which are licensed for non-commercial use. "
              "Also tags the report footer accordingly."),
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Console-script entry point declared in pyproject.toml.

    Delegates to ``analyze.main()`` after the parser is built so the existing
    orchestration body keeps working unchanged. When the V7 decomposition
    extracts ``pipeline.py`` and ``renderers/``, this function will own the
    full orchestration directly and ``analyze.py`` will be deleted (or
    reduced to a back-compat shim).
    """
    import analyze  # local import — see module docstring
    analyze.main(argv=argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
