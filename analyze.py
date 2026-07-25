#!/usr/bin/env python3
"""
DNA Analysis Tool — Local, Privacy-First
Two-tier analysis: deterministic SNP lookup + local Ollama AI interpretation.
Usage: python analyze.py /path/to/raw_dna_file.csv [--no-ai] [--model MODEL]
"""

import sys
import json
import math
import re
import time
import argparse
import datetime
from pathlib import Path
from collections import defaultdict
from typing import Optional, Dict, List, Tuple

from y_haplogroup import analyze_y_haplogroup

try:
    from mt_haplogroup import analyze_mt_haplogroup
except ImportError:
    analyze_mt_haplogroup = None

# Professional-grade analysis modules
try:
    from prs import analyze_polygenic_scores
    from pgx import analyze_pgx
    from interactions import detect_interactions
    from carrier import analyze_carriers
    from traits import predict_traits
    from qc import run_qc
    from references import collect_references_used, get_reference, level_class
    from counseling import evaluate_counseling_triggers
    PROFESSIONAL_MODULES_LOADED = True
except ImportError as _e:
    PROFESSIONAL_MODULES_LOADED = False
    _MODULE_LOAD_ERROR = str(_e)

# v3 modules (graceful degradation when missing)
try:
    from imputation import impute_genotypes, imputation_available
except ImportError:
    impute_genotypes = None
    imputation_available = None
try:
    from pgs_catalog import analyze_expanded_pgs
except ImportError:
    analyze_expanded_pgs = None
try:
    from ancestry_pca import analyze_ancestry
except ImportError:
    analyze_ancestry = None
try:
    from pdf_export import html_to_pdf, weasyprint_available
except ImportError:
    html_to_pdf = None
    weasyprint_available = lambda: False
try:
    from medications import analyze_medications
except ImportError:
    analyze_medications = None
try:
    from family_planning import build_carrier_report, render_carrier_html
except ImportError:
    build_carrier_report = None
    render_carrier_html = None
try:
    from chat import run_chat
except ImportError:
    run_chat = None
try:
    from compare import diff_runs, render_diff_text
except ImportError:
    diff_runs = None
    render_diff_text = None
try:
    from wellness import analyze_wellness
except ImportError:
    analyze_wellness = None

# V5 modules
try:
    from hla import analyze_hla
except ImportError:
    analyze_hla = None
try:
    from roh import detect_roh, render_ideogram_svg
except ImportError:
    detect_roh = None
    render_ideogram_svg = None
try:
    from local_ancestry import analyze_local_ancestry, render_chromosome_painting_svg
except ImportError:
    analyze_local_ancestry = None
    render_chromosome_painting_svg = None
try:
    from phewas import analyze_phewas
except ImportError:
    analyze_phewas = None
try:
    from mendelian_randomization import analyze_mr
except ImportError:
    analyze_mr = None
try:
    from genetic_age import analyze_genetic_age
except ImportError:
    analyze_genetic_age = None
try:
    from pgx_simulation import simulate_pgx
except ImportError:
    simulate_pgx = None
try:
    from reproductive import analyze_reproductive
except ImportError:
    analyze_reproductive = None
try:
    from emergency_card import build_emergency_card
except ImportError:
    build_emergency_card = None
try:
    from narrative import generate_narrative_report
except ImportError:
    generate_narrative_report = None
try:
    from bloodwork import compare_bloodwork, render_bloodwork_html
except ImportError:
    compare_bloodwork = None
    render_bloodwork_html = None
try:
    from supplements import build_supplement_stack, render_supplements_html
except ImportError:
    build_supplement_stack = None
    render_supplements_html = None
try:
    from exercise import analyze_exercise, render_exercise_html
except ImportError:
    analyze_exercise = None
    render_exercise_html = None
try:
    from nutrition import analyze_nutrition, render_nutrition_html
except ImportError:
    analyze_nutrition = None
    render_nutrition_html = None
try:
    from fhir_export import export_fhir
except ImportError:
    export_fhir = None
try:
    from personalized_plan import build_personalized_plan, render_plan_html
except ImportError:
    build_personalized_plan = None
    render_plan_html = None
try:
    from health_economics import analyze_health_economics
except ImportError:
    analyze_health_economics = None
try:
    from health_economics import (analyze_personal_economics,
                                  render_economic_analysis_html)
except ImportError:
    analyze_personal_economics = None
    render_economic_analysis_html = None
try:
    from metal_oxidative import analyze_metal_oxidative
except ImportError:
    analyze_metal_oxidative = None
try:
    from detox import analyze_detox
except ImportError:
    analyze_detox = None

SCRIPT_DIR = Path(__file__).parent
DB_PATH = SCRIPT_DIR / "snp_database.json"
OLLAMA_URL = "http://localhost:11434/api/chat"
OLLAMA_MODEL = "qwen3:14b"
REPORT_VERSION = "6.6.1-premium"

CATEGORY_ORDER = [
    "Hereditary Conditions",
    "Cardiovascular",
    "Cancer Risk",
    "Alzheimer's/Neurological",
    "Neurological",
    "Blood Clotting",
    "Diabetes",
    "Iron Metabolism",
    "Drug Metabolism/Pharmacogenomics",
    "Stimulant & Medication Response",
    "Methylation & Folate",
    "Diet & Nutrition",
    "Athletic Performance",
    "Dopamine & Mental Health",
    "Testosterone & Hormones",
    "Male Fertility & Sperm",
    "Reproductive Health",
    "Detoxification",
    "Heavy Metal & Environmental Toxin Response",
    "Inflammation",
    "Autoimmune & HLA",
    "Sleep",
    "Skin & Aging",
    "Eye Health",
    "Immune System",
    "Respiratory",
    "Gastrointestinal",
    "Musculoskeletal",
    "Kidney",
    "Dental & Oral",
    "Longevity",
    "Wellness",
    "Ancestry Informative Markers",
]

# Category-specific extra questions appended to the base AI prompt.
# Each value is a string that will be injected before the standard 6 questions.
CATEGORY_PROMPTS: Dict[str, str] = {
    "Testosterone & Hormones": (
        "Focus your interpretation on:\n"
        "  • What this genotype profile implies for *natural* testosterone levels, "
        "androgen receptor sensitivity, and bioavailable (free) testosterone.\n"
        "  • DHT conversion (5-alpha reductase activity) — implications for hair "
        "loss, prostate, muscle/strength response.\n"
        "  • Aromatase activity and T:E2 ratio — implications for estrogen balance, "
        "gynecomastia susceptibility, mood/libido/bone effects.\n"
        "  • SHBG variants and how to estimate free vs total T from a lab panel.\n"
        "  • Lifestyle/training/dietary levers that meaningfully shift hormone "
        "profile: sleep, resistance training, body composition, zinc, vitamin D, "
        "magnesium, boron, alcohol, healthy fats.\n"
        "  • When the genetic context warrants discussing free T / SHBG / E2 "
        "testing with a physician.\n"
    ),
    "Male Fertility & Sperm": (
        "Focus your interpretation on:\n"
        "  • Sperm quality optimization based on this genetic profile — DNA "
        "fragmentation risk, motility, and concentration vulnerabilities.\n"
        "  • Oxidative stress impact on fertility — which antioxidants and "
        "cofactors matter most given these variants (CoQ10, vitamin E, vitamin C, "
        "NAC, selenium, zinc, glutathione precursors).\n"
        "  • Supplements with the strongest fertility evidence: zinc, CoQ10, "
        "folate (methylfolate if MTHFR), selenium, L-carnitine, omega-3.\n"
        "  • Lifestyle factors: heat exposure (hot tubs, laptops on lap, tight "
        "underwear, sauna), alcohol, smoking, sleep, body composition.\n"
        "  • When a semen analysis + reproductive endocrinology consult is "
        "appropriate given findings.\n"
    ),
    "Stimulant & Medication Response": (
        "Focus your interpretation on:\n"
        "  • Predicted response to lisdexamfetamine (Vyvanse), other amphetamines "
        "(Adderall, Dexedrine), methylphenidate (Ritalin/Concerta), and atomoxetine "
        "(Strattera) given CYP2D6, CYP2B6, CYP3A4, COMT, DRD2/4, DAT1, ADRA2A, "
        "SLC6A2, and ADHD-susceptibility variants present.\n"
        "  • Optimal dosing considerations — fast vs slow metabolism, ultra-rapid "
        "vs poor metabolizer phenotypes where determinable.\n"
        "  • Side-effect susceptibility: anxiety (COMT Met/Met, ADORA2A T/T), "
        "insomnia, appetite suppression, rebound, cardiovascular sensitivity.\n"
        "  • Which medication CLASS the genetics suggest exploring first "
        "(stimulant vs non-stimulant; amphetamine vs methylphenidate vs guanfacine).\n\n"
        "**CRITICAL DISCLAIMER**: This is genetic context for informed conversation "
        "with a psychiatrist. It is NOT medical advice and NOT prescriptive. "
        "ADHD diagnosis is clinical (DSM-5). Medication decisions — including which "
        "drug to try, at what dose, and how to titrate — must involve a licensed "
        "psychiatrist or qualified prescriber. Genetic factors are ONE input among "
        "many; clinical response and structured outcome monitoring (rating scales, "
        "follow-ups) remain primary.\n"
    ),
    "Heavy Metal & Environmental Toxin Response": (
        "Focus your interpretation on:\n"
        "  • Individual vulnerability to specific exposures: lead, mercury "
        "(methylmercury from fish, dental amalgam), arsenic (inorganic in "
        "rice/well water), cadmium (smoke, some grain), aluminum, and fluoride.\n"
        "  • Which water filtration system would matter most for this genetic "
        "profile (carbon block, reverse osmosis, distillation) and what it "
        "should target.\n"
        "  • Supplements that support detox pathways: N-acetylcysteine, "
        "selenium, zinc, glutathione precursors, sulforaphane, chlorella, "
        "vitamin C, magnesium — match the supplement to the pathway "
        "(GST, MT, GPX, methylation).\n"
        "  • Occupational/environmental considerations: smoking (cadmium, "
        "PAHs), traffic-pollution proximity, char-grilled meat, drinking-water "
        "source, fish consumption (mercury vs omega-3 trade-off).\n"
        "  • When provoked vs unprovoked heavy-metal urine testing is "
        "informative (and the risks of chelation without proper medical "
        "supervision).\n"
    ),
    "Longevity": (
        "Focus your interpretation on:\n"
        "  • What this profile suggests about pathways of healthy aging — "
        "stress resistance, inflammaging, telomere biology, lipid metabolism, "
        "FOXO/IGF-1 signalling.\n"
        "  • Evidence-based longevity levers most relevant to these variants: "
        "exercise (zone 2 + VO2max + resistance), Mediterranean / "
        "minimally-processed diet, time-restricted eating, sleep, social "
        "connection, polyphenols, stress management.\n"
        "  • Cardiovascular and dementia prevention focus areas given the "
        "genetic profile.\n"
        "  • Pitfalls — supplements/interventions with weak evidence vs "
        "strong evidence; avoid over-supplementation.\n"
    ),
    "Ancestry Informative Markers": (
        "Focus your interpretation on:\n"
        "  • What these markers suggest about continental and regional "
        "ancestry — paternal/maternal lineages (Y-DNA, mtDNA) and "
        "autosomal pigmentation/trait markers (SLC24A5, SLC45A2, EDAR, "
        "HERC2, MC1R, etc.).\n"
        "  • Migration history and population-genetic context — admixture "
        "events, founder populations, geographic distribution.\n"
        "  • Practical implications: vitamin D needs (skin pigmentation × "
        "latitude), Duffy-null neutropenia awareness, alcohol metabolism "
        "(ALDH2/ADH1B), lactase persistence.\n"
        "  • Make clear that these are population-average associations — "
        "individual ancestry is complex and not reducible to a few SNPs.\n"
    ),
}

APOE_INFO = {
    "E2/E2": {
        "risk_label": "Reduced Risk",
        "color_class": "apoe-green",
        "description": (
            "APOE E2/E2 (~1% of people) is associated with significantly reduced Alzheimer's "
            "risk and generally lower LDL cholesterol. However, this genotype is linked to a "
            "rare condition called Type III hyperlipoproteinemia in some individuals, where "
            "triglycerides can be markedly elevated. Discuss cardiovascular lipid monitoring "
            "with your doctor."
        ),
    },
    "E2/E3": {
        "risk_label": "Below Average Risk",
        "color_class": "apoe-green",
        "description": (
            "APOE E2/E3 confers below-average Alzheimer's risk and is typically associated "
            "with favorable cholesterol profiles. This is a protective genotype for "
            "neurodegeneration compared to the most common E3/E3."
        ),
    },
    "E3/E3": {
        "risk_label": "Average Risk",
        "color_class": "apoe-blue",
        "description": (
            "APOE E3/E3 is the most common genotype (~60–65% of the population) and is "
            "associated with average lifetime risk for Alzheimer's disease and cardiovascular "
            "disease. Standard prevention measures apply — the lifestyle factors below are "
            "meaningful even at average genetic risk."
        ),
    },
    "E2/E4": {
        "risk_label": "Mixed (Near Average)",
        "color_class": "apoe-yellow",
        "description": (
            "APOE E2/E4 carries opposing influences: E2 is protective while E4 increases risk. "
            "The net effect on Alzheimer's risk is roughly average, though the E4 allele may "
            "modestly increase cardiovascular risk. A heart-healthy lifestyle is particularly "
            "important. This is a rare genotype (~1% of people)."
        ),
    },
    "E3/E4": {
        "risk_label": "Above Average Risk",
        "color_class": "apoe-orange",
        "description": (
            "APOE E3/E4 (carried by ~20–25% of people) increases Alzheimer's risk approximately "
            "3-fold compared to E3/E3. It is also associated with higher LDL cholesterol and "
            "increased cardiovascular risk. Proactive lifestyle measures have the strongest "
            "evidence for risk reduction: daily aerobic exercise, Mediterranean/MIND diet, "
            "quality sleep (7–9 hours), cognitive engagement, and blood pressure control."
        ),
    },
    "E4/E4": {
        "risk_label": "High Risk",
        "color_class": "apoe-red",
        "description": (
            "APOE E4/E4 (~2% of the population) is the highest-risk APOE genotype, increasing "
            "Alzheimer's risk approximately 8–12 fold and significantly elevating cardiovascular "
            "risk. This does NOT mean Alzheimer's is certain — many E4/E4 carriers never develop "
            "it. However, preventive action is especially important: daily aerobic exercise "
            "(proven to reduce amyloid burden), Mediterranean or MIND diet, quality sleep "
            "(sleep clears amyloid beta), stress management, avoiding head trauma, and "
            "controlling blood pressure and cholesterol. Consider discussing risk-reduction "
            "strategies with a neurologist. Clinical trials for prevention in high-risk "
            "individuals are ongoing."
        ),
    },
}


def log(msg: str) -> None:
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def load_snp_database() -> Dict:
    if not DB_PATH.exists():
        log(f"ERROR: SNP database not found at {DB_PATH}")
        sys.exit(1)
    with open(DB_PATH) as f:
        return json.load(f)


def parse_dna_file(filepath: str):
    """Parse a raw DNA file using the snps library (auto-detects format)."""
    try:
        from snps import SNPs
    except ImportError:
        log("ERROR: snps library not installed. Run: pip install snps")
        sys.exit(1)

    path = Path(filepath)
    if not path.exists():
        log(f"ERROR: File not found: {filepath}")
        sys.exit(1)

    log(f"Parsing DNA file: {path.name} ...")
    try:
        # assign_par_snps=False avoids some internet requests
        s = SNPs(str(path), assign_par_snps=False)
    except Exception as e:
        log(f"ERROR: Failed to parse DNA file: {e}")
        sys.exit(1)

    if s.snps is None or s.snps.empty:
        log("ERROR: No SNPs found in file. Check the file format.")
        sys.exit(1)

    log(f"  Loaded {len(s.snps):,} SNPs | Format detected: {s.source}")

    # V8: tag every parsed variant with provenance + auto-detect chip build.
    # `source` column = "chip" for everything at this stage; imputation
    # overwrites imputed rows later. `attrs["build"]` is one of
    # grch37 / grch38 / mixed / unknown — surfaced by the QC card.
    try:
        from provenance import annotate_parsed
        snps_df = annotate_parsed(s.snps)
        log(f"  Build auto-detected: {snps_df.attrs.get('build', 'unknown')}")
    except ImportError:
        snps_df = s.snps
    return snps_df


def determine_apoe(rs429358_gt: Optional[str], rs7412_gt: Optional[str]) -> Optional[str]:
    """
    Determine APOE genotype from the two defining SNPs.
    rs429358 C allele = E4 haplotype
    rs7412   T allele = E2 haplotype
    """
    if not rs429358_gt or not rs7412_gt:
        return None

    gt_429358 = str(rs429358_gt).upper().replace(" ", "").replace("-", "")
    gt_7412 = str(rs7412_gt).upper().replace(" ", "").replace("-", "")

    if len(gt_429358) < 1 or len(gt_7412) < 1:
        return None

    e4_count = gt_429358.count("C")
    e2_count = gt_7412.count("T")
    e3_count = 2 - e4_count - e2_count

    if e3_count < 0:
        return None

    alleles = sorted(["E4"] * e4_count + ["E2"] * e2_count + ["E3"] * e3_count)
    if len(alleles) == 2:
        return f"{alleles[0]}/{alleles[1]}"
    return None


def tier1_lookup(
    snps_df, database: Dict
) -> Tuple[List[Dict], Optional[str]]:
    """Match the person's SNPs against the curated database."""
    log("Matching SNPs against database ...")

    results = []
    apoe_gts: Dict[str, str] = {}

    for rsid, entry in database.items():
        if rsid not in snps_df.index:
            continue

        row = snps_df.loc[rsid]
        genotype = row.get("genotype", None)
        if genotype is None or str(genotype).strip() in ("nan", "--", ""):
            continue

        genotype = str(genotype).upper().strip()
        risk_allele = entry["risk_allele"].upper()

        # Handle insertion/deletion markers
        if risk_allele in ("INS", "DEL", "I", "D"):
            risk_copies = 0  # Can't reliably count from raw genotype string
        else:
            risk_copies = genotype.count(risk_allele)

        # Provenance: distinguish a directly-typed call from a statistically
        # imputed one. When --impute is used, snps_df carries `source`
        # ('chip'/'imputed') and `r2` columns; surface them so imputed matches
        # (which can be many) are never mistaken for measured genotypes.
        source = row.get("source", "chip") if hasattr(row, "get") else "chip"
        r2 = row.get("r2") if hasattr(row, "get") else None
        try:
            r2 = round(float(r2), 2) if r2 is not None and str(r2) != "nan" else None
        except (TypeError, ValueError):
            r2 = None

        results.append(
            {
                "rsid": rsid,
                "gene": entry["gene"],
                "variant_name": entry["name"],
                "category": entry["category"],
                "my_genotype": genotype,
                "risk_allele": risk_allele,
                "risk_copies": risk_copies,
                "significance": entry["significance"],
                "summary": entry["summary"],
                "recommendation": entry["recommendation"],
                "cross_references": entry.get("cross_references", []),
                "chip_coverage_note": entry.get("chip_coverage_note"),
                "source": source if source in ("chip", "imputed") else "chip",
                "r2": r2,
            }
        )

        if rsid in ("rs429358", "rs7412"):
            apoe_gts[rsid] = genotype

    apoe_genotype = determine_apoe(
        apoe_gts.get("rs429358"), apoe_gts.get("rs7412")
    )
    log(
        f"  Matched {len(results)} variants | "
        f"APOE: {apoe_genotype or 'not determined'}"
    )
    return results, apoe_genotype


def call_ollama(prompt: str, model: str, timeout: int = 600) -> str:
    """Send a prompt to local Ollama and return the response text."""
    import requests

    payload = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "stream": False,
        "options": {"temperature": 0.3, "num_ctx": 8192},
    }

    try:
        resp = requests.post(OLLAMA_URL, json=payload, timeout=timeout)
        resp.raise_for_status()
        content = resp.json()["message"]["content"]
        # Strip <think>...</think> blocks from thinking models (e.g. qwen3)
        content = re.sub(r"<think>.*?</think>", "", content, flags=re.DOTALL).strip()
        return content
    except requests.exceptions.ConnectionError:
        raise ConnectionError(
            "Cannot connect to Ollama at localhost:11434. "
            "Start it with: ollama serve"
        )
    except requests.exceptions.Timeout:
        raise TimeoutError("Ollama request timed out.")
    except Exception as e:
        raise RuntimeError(f"Ollama error: {e}")


def build_category_map(tier1_results: List[Dict]) -> Dict[str, List[Dict]]:
    """Group results by category, expanding each result into every category
    listed in its cross_references. Cross-referenced entries are shallow
    copies with summary replaced by the cross-ref implication and an
    is_cross_ref flag set so the renderer can label them.
    """
    by_cat: Dict[str, List[Dict]] = defaultdict(list)
    for r in tier1_results:
        by_cat[r["category"]].append(r)
        for xref in r.get("cross_references", []) or []:
            xcat = xref.get("category")
            if not xcat:
                continue
            shadow = dict(r)
            shadow["summary"] = xref.get("implication", r["summary"])
            shadow["is_cross_ref"] = True
            shadow["primary_category"] = r["category"]
            by_cat[xcat].append(shadow)
    return by_cat


# Categories with more than this many variants are split into sub-batches when
# sent to the local AI, to keep each request small enough to finish within the
# Ollama timeout. Empirically tuned for qwen3:14b on consumer hardware.
AI_BATCH_MAX = 15


def _split_into_batches(items: List[Dict], target_max: int = AI_BATCH_MAX) -> List[List[Dict]]:
    """Split a list of variants into roughly equal batches, each no larger
    than target_max. For total counts > target_max the batches will each be
    in the 12-15 range whenever the math allows; smaller totals stay as one
    batch and pathologically small overflows (e.g. 16 -> 8+8) accept slightly
    smaller batches rather than exceed target_max.
    """
    n = len(items)
    if n <= target_max:
        return [items]
    n_batches = math.ceil(n / target_max)
    base, extras = divmod(n, n_batches)
    out: List[List[Dict]] = []
    start = 0
    for i in range(n_batches):
        size = base + (1 if i < extras else 0)
        out.append(items[start:start + size])
        start += size
    return out


def _build_category_prompt(
    category: str,
    snps: List[Dict],
    batch_idx: Optional[int] = None,
    total_batches: Optional[int] = None,
) -> str:
    """Build the AI prompt for one category (or one sub-batch of a category)."""
    lines = []
    for s in snps:
        if s["risk_copies"] == 0:
            status = "no risk allele (normal genotype)"
        elif s["risk_copies"] == 1:
            status = "1 copy of risk allele (heterozygous)"
        else:
            status = "2 copies of risk allele (homozygous)"
        xref_tag = ""
        if s.get("is_cross_ref"):
            xref_tag = f" [cross-referenced from {s.get('primary_category', '?')}]"
        lines.append(
            f"• {s['rsid']} | {s['gene']} {s['variant_name']}{xref_tag} | "
            f"Genotype: {s['my_genotype']} | Status: {status} | "
            f"Significance: {s['significance']}\n"
            f"    Notes: {s['summary']}"
        )

    category_focus = CATEGORY_PROMPTS.get(category, "")

    batch_note = ""
    if total_batches and total_batches > 1:
        batch_note = (
            f"NOTE: This is batch {batch_idx} of {total_batches} for this "
            f"category — the variants were split into smaller groups so each "
            f"request fits comfortably within the local model's context. "
            f"Interpret ONLY the variants shown below; other batches cover "
            f"the rest of the category.\n\n"
        )

    return (
        "You are a knowledgeable genetic counselor providing educational "
        "(non-diagnostic) information about DNA variants.\n\n"
        + batch_note
        + f"=== CATEGORY: {category} ===\n"
        "Variants found in this person's DNA:\n\n"
        + "\n".join(lines)
        + "\n\n"
        + (category_focus + "\n" if category_focus else "")
        + "Please provide:\n\n"
        "**1. Combined Interpretation**\n"
        "What do these variants mean together? Note any important "
        "gene-gene or variant-variant interactions within this category.\n\n"
        "**2. Prioritized Actionable Findings**\n"
        "Which findings are most immediately actionable and relevant to "
        "daily life? Rank them by priority.\n\n"
        "**3. Risk Interactions**\n"
        "Do any combinations amplify or reduce each other's risk? "
        "(e.g. two pro-inflammatory variants together vs. one alone)\n\n"
        "**4. Specific Recommendations**\n"
        "Concrete lifestyle, dietary, supplement, or monitoring "
        "recommendations. Be specific (e.g. '1000mg/day of X' rather "
        "than 'supplement with X').\n\n"
        "**5. Evidence Quality**\n"
        "For key findings, note: well-established (replicated in large "
        "trials) vs. emerging evidence vs. preliminary (single studies).\n\n"
        "**6. Discuss With Your Doctor**\n"
        "Specific tests, referrals, or conversations this person should "
        "have with their healthcare provider.\n\n"
        "Be focused and practical. Remind the reader this is educational, "
        "not diagnostic."
    )


def _run_category_ai_with_batching(
    category: str, snps: List[Dict], model: str
) -> Tuple[str, bool]:
    """Run the per-category AI interpretation, splitting into sub-batches if
    the category has more than AI_BATCH_MAX variants. Returns
    (combined_ai_text, had_failure). If had_failure is True the caller should
    record this category for retry; the returned text still contains the
    best partial output we managed to obtain.
    """
    batches = _split_into_batches(snps, target_max=AI_BATCH_MAX)

    if len(batches) == 1:
        try:
            return (
                call_ollama(_build_category_prompt(category, batches[0]), model=model),
                False,
            )
        except Exception as e:
            log(f"    WARNING: AI analysis failed for {category}: {e}")
            return f"*AI analysis unavailable: {e}*", True

    sizes = ", ".join(str(len(b)) for b in batches)
    log(f"    Splitting {len(snps)} variants into {len(batches)} sub-batches ({sizes})")
    parts: List[str] = []
    had_failure = False
    offset = 0
    for i, batch in enumerate(batches, start=1):
        first = offset + 1
        last = offset + len(batch)
        log(f"    Batch {i}/{len(batches)} (variants {first}-{last}) ...")
        try:
            text = call_ollama(
                _build_category_prompt(category, batch, i, len(batches)),
                model=model,
            )
            parts.append(
                f"**Batch {i} of {len(batches)} (variants {first}–{last})**\n\n{text}"
            )
        except Exception as e:
            log(f"    WARNING: Batch {i}/{len(batches)} failed for {category}: {e}")
            had_failure = True
            parts.append(
                f"**Batch {i} of {len(batches)} (variants {first}–{last})** — "
                f"*AI analysis unavailable: {e}*"
            )
        offset += len(batch)
        time.sleep(0.5)

    return "\n\n---\n\n".join(parts), had_failure


def tier2_analysis(
    tier1_results: List[Dict], apoe_genotype: Optional[str], model: str
) -> Tuple[Dict[str, str], str, List[Dict]]:
    """Run per-category AI interpretation then an executive summary.

    Returns (ai_results, exec_summary, failed_categories). failed_categories
    is a list of {"category", "snps"} dicts for any categories whose AI
    analysis didn't fully succeed; callers persist this list so the user can
    re-run just those categories via --retry-failed.
    """
    categories = build_category_map(tier1_results)

    ai_results: Dict[str, str] = {}
    failed_categories: List[Dict] = []

    for category in sorted(categories.keys()):
        snps = categories[category]
        n = len(snps)
        suffix = ", will batch" if n > AI_BATCH_MAX else ""
        log(f"  AI analyzing: {category} ({n} variants{suffix}) ...")

        ai_text, had_failure = _run_category_ai_with_batching(category, snps, model)
        ai_results[category] = ai_text
        if had_failure:
            failed_categories.append({
                "category": category,
                "snps": list(snps),
            })

        time.sleep(0.5)

    # Executive summary
    log("  AI generating executive summary ...")

    high_risk = [
        r for r in tier1_results
        if r["significance"] == "high" and r["risk_copies"] > 0
    ]
    moderate_risk = [
        r for r in tier1_results
        if r["significance"] == "moderate" and r["risk_copies"] > 0
    ]

    summary_lines = []
    if apoe_genotype:
        apoe_info = APOE_INFO.get(apoe_genotype, {})
        summary_lines.append(
            f"APOE Genotype: {apoe_genotype} "
            f"({apoe_info.get('risk_label', 'see report')})"
        )
    for r in high_risk:
        summary_lines.append(
            f"HIGH SIGNIFICANCE: {r['gene']} {r['variant_name']} — "
            f"{r['risk_copies']} risk allele(s)"
        )
    for r in moderate_risk[:20]:
        summary_lines.append(
            f"MODERATE: {r['gene']} {r['variant_name']} — "
            f"{r['risk_copies']} risk allele(s)"
        )

    exec_prompt = (
        "You are a genetic counselor providing an executive summary of a "
        "comprehensive DNA analysis. Based on the findings below, write "
        "a 4–6 paragraph executive summary.\n\n"
        "Key findings:\n"
        + "\n".join(summary_lines)
        + f"\n\nTotal variants analyzed: {len(tier1_results)}\n"
        f"Categories covered: "
        f"{', '.join(sorted(set(r['category'] for r in tier1_results)))}\n\n"
        "Your summary should:\n"
        "1. Highlight the most important genetic findings and their overall "
        "meaning\n"
        "2. Identify the top 3–5 highest-priority actionable areas\n"
        "3. Note cross-category patterns (e.g. compounding inflammation "
        "variants, methylation + detox interactions)\n"
        "4. Provide a balanced, grounded perspective — not alarmist, not "
        "dismissive\n"
        "5. Give clear next steps: what to prioritize, monitor, or discuss "
        "with a doctor\n\n"
        "Use plain language for a scientifically literate non-specialist. "
        "This is educational information only, not a medical diagnosis."
    )

    try:
        exec_summary = call_ollama(exec_prompt, model=model)
    except Exception as e:
        log(f"  WARNING: Executive summary failed: {e}")
        exec_summary = f"*Executive summary unavailable: {e}*"

    return ai_results, exec_summary, failed_categories


def cross_category_synthesis(
    tier1_results: List[Dict], model: str
) -> str:
    """Generate a synthesis of how findings across categories interact.

    This is a separate AI call that focuses on cross-category compounding —
    e.g. methylation + heavy-metal detox, testosterone + fertility, dopamine +
    stimulant response — rather than within-category interpretation.
    """
    log("  AI generating Cross-Category Interactions synthesis ...")

    # Group by category, prioritise risk-carrying variants and any with cross-refs
    by_cat: Dict[str, List[Dict]] = defaultdict(list)
    cross_ref_pairs: List[Tuple[str, str, Dict]] = []
    for r in tier1_results:
        if r["risk_copies"] > 0 or r.get("cross_references"):
            by_cat[r["category"]].append(r)
        for xref in r.get("cross_references", []) or []:
            cross_ref_pairs.append((r["category"], xref.get("category", "?"), r))

    cat_blocks = []
    for cat in sorted(by_cat.keys()):
        items = by_cat[cat]
        if not items:
            continue
        lines = []
        for s in items[:12]:  # cap per category to keep prompt reasonable
            lines.append(
                f"  • {s['gene']} {s['variant_name']} ({s['rsid']}) — "
                f"{s['risk_copies']} risk allele(s), "
                f"{s['significance']} significance"
            )
        cat_blocks.append(f"**{cat}**\n" + "\n".join(lines))

    xref_lines = []
    for primary, secondary, snp in cross_ref_pairs:
        xref_lines.append(
            f"  • {snp['gene']} {snp['variant_name']} bridges "
            f"{primary} ↔ {secondary}"
        )

    if not cat_blocks:
        return "*No risk-carrying or cross-referenced variants were found, so no cross-category synthesis was generated.*"

    prompt = (
        "You are a genetic counselor synthesising findings ACROSS multiple "
        "categories of a comprehensive DNA report. Per-category interpretations "
        "have already been generated separately. Your job here is different: "
        "identify how findings from DIFFERENT categories COMPOUND or INTERACT.\n\n"
        "Findings by category (risk-carrying or cross-referenced variants only):\n\n"
        + "\n\n".join(cat_blocks)
        + "\n\nVariants that already have explicit cross-category bridges in the "
        "database:\n"
        + ("\n".join(xref_lines) if xref_lines else "  (none)")
        + "\n\n"
        "Please write a Cross-Category Interactions analysis covering 4–7 of the "
        "most clinically meaningful compounding patterns for this person. For "
        "each pattern:\n\n"
        "**Pattern name** (e.g. 'Methylation × Heavy-Metal Detox')\n"
        "- Which variants from which categories combine\n"
        "- What the combined biology means (mechanism in 1–2 sentences)\n"
        "- The PRACTICAL implication — what this person should prioritise that "
        "they would NOT prioritise from a single-category view\n\n"
        "Patterns to consider when present:\n"
        "  • Methylation + Heavy-Metal Detox (folate/B12 cycle drives metal excretion)\n"
        "  • Testosterone + Fertility (hormone milieu shapes spermatogenesis)\n"
        "  • Dopamine + Stimulant Response (same receptors matter for ADHD med choice)\n"
        "  • Inflammation + Cardiovascular (chronic inflammation drives atherosclerosis)\n"
        "  • Inflammation + Fertility (oxidative stress harms sperm)\n"
        "  • Detoxification + Cancer Risk (slow detox amplifies carcinogen damage)\n"
        "  • Estrogen metabolism (CYP1B1, COMT) + Cancer Risk\n"
        "  • APOE/lipids + Cognition + Cardiovascular\n"
        "  • Iron (HFE) + Heavy Metal handling (transport-pathway competition)\n"
        "  • Insulin resistance × Lipids × Inflammation (metabolic syndrome)\n\n"
        "Only include patterns that are actually relevant to this person's variants. "
        "Be specific — name the variants, name the action.\n\n"
        "Reminder: educational, not diagnostic."
    )

    try:
        return call_ollama(prompt, model=model)
    except Exception as e:
        log(f"  WARNING: Cross-category synthesis failed: {e}")
        return f"*Cross-category synthesis unavailable: {e}*"


FAILED_CATEGORIES_PATH = SCRIPT_DIR / "failed_categories.json"


def write_failed_categories(
    failed: List[Dict],
    model: str,
    report_path: Path,
    dna_file: Optional[str] = None,
) -> None:
    """Persist the list of categories whose AI analysis didn't fully succeed
    so the user can re-run just those via --retry-failed. Always overwrites
    the file, so a clean run after a partial one produces an empty failed
    list rather than stale entries.
    """
    payload = {
        "timestamp": datetime.datetime.now().isoformat(),
        "model": model,
        "report_path": str(report_path),
        "dna_file": dna_file,
        "failed": failed,
    }
    FAILED_CATEGORIES_PATH.write_text(
        json.dumps(payload, indent=2, default=str),
        encoding="utf-8",
    )
    if failed:
        names = ", ".join(e["category"] for e in failed)
        log(f"  Logged {len(failed)} failed categor"
            f"{'y' if len(failed) == 1 else 'ies'} -> "
            f"{FAILED_CATEGORIES_PATH.name} ({names})")
        log(f"  Re-run with: python analyze.py --retry-failed")


def retry_failed_categories(model: str, override_report: Optional[Path] = None) -> int:
    """Read failed_categories.json, re-run AI for each entry, and patch the
    new interpretations into the existing report.html. Returns a process
    exit code (0 success, non-zero if something prevented the retry).
    """
    if not FAILED_CATEGORIES_PATH.exists():
        log(f"  --retry-failed: no {FAILED_CATEGORIES_PATH.name} found in "
            f"{SCRIPT_DIR}. Nothing to retry.")
        return 1

    try:
        payload = json.loads(FAILED_CATEGORIES_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        log(f"  --retry-failed: failed_categories.json is malformed: {e}")
        return 1

    failed = payload.get("failed", []) or []
    if not failed:
        log("  --retry-failed: no failed categories recorded. Nothing to do.")
        return 0

    report_path = override_report or Path(payload.get("report_path") or (SCRIPT_DIR / "report.html"))
    if not report_path.exists():
        log(f"  --retry-failed: report not found at {report_path}. Aborting.")
        return 1

    html = report_path.read_text(encoding="utf-8")
    retry_model = model or payload.get("model") or OLLAMA_MODEL

    log("=" * 60)
    log(f"Retrying AI for {len(failed)} failed categor"
        f"{'y' if len(failed) == 1 else 'ies'} (model: {retry_model})")
    log(f"  Patching: {report_path}")
    log("=" * 60)

    still_failed: List[Dict] = []
    patched = 0
    for entry in failed:
        category = entry.get("category")
        snps = entry.get("snps") or []
        if not category or not snps:
            log(f"  WARNING: skipping malformed entry: {entry!r}")
            still_failed.append(entry)
            continue

        n = len(snps)
        suffix = ", will batch" if n > AI_BATCH_MAX else ""
        log(f"  AI re-analyzing: {category} ({n} variants{suffix}) ...")

        ai_text, had_failure = _run_category_ai_with_batching(category, snps, retry_model)
        if had_failure:
            still_failed.append({"category": category, "snps": snps})

        cat_id = _cat_id(category)
        html, ok = _patch_ai_section_in_html(html, cat_id, md_to_html(ai_text))
        if ok:
            patched += 1
            log(f"    Patched '{category}' into {report_path.name}")
        else:
            log(f"    WARNING: could not locate <details id=\"{cat_id}\"> "
                f"in {report_path.name}; AI text not patched")

        time.sleep(0.5)

    report_path.write_text(html, encoding="utf-8")

    payload["failed"] = still_failed
    payload["timestamp"] = datetime.datetime.now().isoformat()
    payload["last_action"] = "retry"
    payload["model"] = retry_model
    FAILED_CATEGORIES_PATH.write_text(
        json.dumps(payload, indent=2, default=str),
        encoding="utf-8",
    )

    log("=" * 60)
    log(f"Retry complete. Patched {patched}/{len(failed)} categories in {report_path}.")
    if still_failed:
        names = ", ".join(e["category"] for e in still_failed)
        log(f"Still failing: {len(still_failed)} ({names}). "
            f"Re-run --retry-failed to try again.")
    else:
        log("All previously failed categories now succeeded.")
    log("=" * 60)
    return 0


# ── HTML helpers ──────────────────────────────────────────────────────────────

# ── Renderers ────────────────────────────────────────────────────────────────
# V8 extraction — the 3 038 lines that previously lived here (md_to_html,
# significance_badge, risk_indicator, Y_ANCIENT_DNA_REFS, 23 build_*_html
# functions, and build_html_report) now live in renderers.py.
#
# An eager re-export from this module was the source of a hard circular-
# import deadlock when analyze.py runs as a script — renderers needs
# analyze's CATEGORY_ORDER / APOE_INFO / build_category_map at module-load
# time, and `from renderers import …` here would re-enter analyze before
# those names existed in the second-load. Internal callers
# (pipeline.run_pipeline) now import directly from renderers; the
# back-compat shim below is a *lazy* attribute hook so legacy external
# code doing `from analyze import build_html_report` still works without
# the eager-load cycle.
_RENDERER_NAMES = frozenset({
    "Y_ANCIENT_DNA_REFS", "_cat_id", "_esc", "_patch_ai_section_in_html",
    "md_to_html", "risk_indicator", "significance_badge",
    "build_ancestry_html", "build_carrier_html", "build_counseling_html",
    "build_expanded_pgs_html", "build_genetic_age_html", "build_hla_html",
    "build_html_report", "build_imputation_html", "build_interactions_html",
    "build_economics_html",
    "build_local_ancestry_html", "build_medications_html", "build_mr_html",
    "build_mtdna_html", "build_pgx_html", "build_pgx_sim_html",
    "build_phewas_html", "build_prs_html", "build_qc_html",
    "build_references_html", "build_reproductive_html", "build_roh_html",
    "build_traits_html", "build_wellness_html",
    "build_ydna_ancient_block", "build_ydna_html",
})


def __getattr__(name):
    """Lazy back-compat: resolve renderer names on first access only."""
    if name in _RENDERER_NAMES:
        import renderers
        return getattr(renderers, name)
    raise AttributeError(f"module 'analyze' has no attribute {name!r}")



# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None) -> None:
    # V7: parser lives in cli.py — this lets the parser be tested and the
    # `dna-analyze` console-script entry point use a small, focused module.
    from cli import build_parser
    parser = build_parser()
    args = parser.parse_args(argv)

    # V8: orchestration body lives in pipeline.run_pipeline.
    # analyze.main() is now a thin shim — parser + delegation.
    from pipeline import run_pipeline
    rc = run_pipeline(args)
    if rc:
        sys.exit(rc)

if __name__ == "__main__":
    main()
