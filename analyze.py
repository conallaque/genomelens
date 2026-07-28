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
    from family_planning import (build_carrier_report, render_carrier_html,
                                  analyze_family_planning)
except ImportError:
    build_carrier_report = None
    render_carrier_html = None
    analyze_family_planning = None
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
try:
    from urologic import analyze_urologic
except ImportError:
    analyze_urologic = None
try:
    from deep_ancestry import analyze_deep_ancestry
except ImportError:
    analyze_deep_ancestry = None
try:
    from blood_type import analyze_blood_type
except ImportError:
    analyze_blood_type = None
try:
    from immunogenetics import analyze_immunogenetics
except ImportError:
    analyze_immunogenetics = None
try:
    from ancestral_story import analyze_ancestral_story
except ImportError:
    analyze_ancestral_story = None
try:
    from neurochemistry import analyze_neurochemistry
except ImportError:
    analyze_neurochemistry = None
try:
    from holistic_synthesis import analyze_holistic_synthesis
except ImportError:
    analyze_holistic_synthesis = None
try:
    from addiction_genetics import analyze_addiction_genetics
except ImportError:
    analyze_addiction_genetics = None
try:
    from polygenic_traits import analyze_polygenic_traits
except ImportError:
    analyze_polygenic_traits = None
try:
    from environmental_optimization import analyze_environmental_optimization
except ImportError:
    analyze_environmental_optimization = None
try:
    from life_stage_playbook import analyze_life_stage_playbook
except ImportError:
    analyze_life_stage_playbook = None
try:
    from clinical_variants import analyze_clinical_variants
except ImportError:
    analyze_clinical_variants = None
try:
    from multi_person_module import (analyze_multi_person,
                                      render_multi_person_html)
except ImportError:
    analyze_multi_person = None
    render_multi_person_html = None

SCRIPT_DIR = Path(__file__).parent
DB_PATH = SCRIPT_DIR / "snp_database.json"
OLLAMA_URL = "http://localhost:11434/api/chat"
OLLAMA_MODEL = "qwen3:14b"
REPORT_VERSION = "6.22.0-premium"

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


# ONE context size for EVERY Ollama call. Changing num_ctx between calls forces
# Ollama to reload the model, and reloading a large model (e.g. the 30B on a
# 25GB-unified Mac) intermittently OOMs → HTTP 500. Keeping num_ctx constant
# means the model loads once and never reloads. 16384 fits the largest prompt
# (cross-category synthesis, capped at 48k chars ≈ 12k tokens + output).
AI_NUM_CTX = 16384


def call_ollama(prompt: str, model: str, timeout: int = 1800,
                stream: bool = True, num_ctx: int = AI_NUM_CTX,
                num_predict: int = 1024, think: Optional[bool] = None,
                retries: int = 2) -> str:
    """Send a prompt to local Ollama and return the response text.

    Streaming (default) resets a 90s read-idle guard on every token so a long
    reasoning generation never trips the socket timeout. ``keep_alive`` holds
    the model resident between calls, and transient 5xx errors (model reload /
    momentary OOM) are retried with backoff — the whole tier-2 pass makes many
    sequential calls, so a single transient blip must not blank a section."""
    import requests

    payload = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "stream": stream,
        # Constant num_ctx (see AI_NUM_CTX) — avoids model-reload thrash.
        "options": {"temperature": 0.3, "num_ctx": num_ctx,
                    "num_predict": num_predict},
        # Keep the model loaded between the many sequential tier-2 calls.
        "keep_alive": "15m",
    }
    # Thinking-model control: for direct synthesis/interpretation tasks we
    # disable reasoning (think=False) so the whole num_predict budget produces
    # the visible answer rather than being consumed inside a <think> block —
    # which otherwise strips to an empty result. Left at the model default when
    # think is None.
    if think is not None:
        payload["think"] = think

    def _attempt() -> str:
        if not stream:
            resp = requests.post(OLLAMA_URL, json=payload, timeout=timeout)
            resp.raise_for_status()
            content = resp.json()["message"]["content"]
        else:
            content_parts: List[str] = []
            last_token_t = time.time()
            with requests.post(OLLAMA_URL, json=payload, stream=True,
                               timeout=(30, timeout)) as resp:
                resp.raise_for_status()
                for line in resp.iter_lines(decode_unicode=True):
                    if not line:
                        if time.time() - last_token_t > 90:
                            raise TimeoutError(
                                "Ollama stream idle for >90s (model stalled)")
                        continue
                    try:
                        obj = json.loads(line)
                    except Exception:
                        continue
                    piece = (obj.get("message") or {}).get("content", "")
                    if piece:
                        content_parts.append(piece)
                        last_token_t = time.time()
                    if obj.get("done"):
                        break
            content = "".join(content_parts)
        raw = content
        stripped = re.sub(r"<think>.*?</think>", "", content, flags=re.DOTALL).strip()
        if not stripped and raw.strip():
            after = re.sub(r"^.*</think>", "", raw, flags=re.DOTALL).strip()
            stripped = after or raw.replace("<think>", "").replace("</think>", "").strip()
        return stripped

    last_err: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            return _attempt()
        except requests.exceptions.ConnectionError:
            raise ConnectionError(
                "Cannot connect to Ollama at localhost:11434. "
                "Start it with: ollama serve")
        except requests.exceptions.HTTPError as e:
            status = getattr(getattr(e, "response", None), "status_code", None)
            last_err = e
            if status and 500 <= status < 600 and attempt < retries:
                time.sleep(3 * (attempt + 1))   # brief backoff, then retry
                continue
            raise RuntimeError(f"Ollama error: {e}")
        except (requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
            raise TimeoutError("Ollama request timed out.")
        except Exception as e:
            raise RuntimeError(f"Ollama error: {e}")
    raise RuntimeError(f"Ollama error after retries: {last_err}")


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
# Ollama read-idle window. Empirically tuned for qwen3:14b in reasoning mode
# on consumer hardware — an 8-variant prompt reliably completes in ~2–4 min,
# whereas 15 variants was hitting the wall clock on large panels.
AI_BATCH_MAX = 8
AI_BATCH_MIN = 2   # never split below this


def _split_evenly(items: List[Dict], max_size: int) -> List[List[Dict]]:
    """Split a list into batches of at most `max_size`, sizes as equal as
    possible (so the last batch isn't a lonely small one)."""
    n = len(items)
    if n <= max_size:
        return [items]
    import math as _m
    n_batches = _m.ceil(n / max_size)
    base, extras = divmod(n, n_batches)
    out: List[List[Dict]] = []
    start = 0
    for i in range(n_batches):
        size = base + (1 if i < extras else 0)
        out.append(items[start:start + size])
        start += size
    return out


def _split_into_batches(items: List[Dict], target_max: int = AI_BATCH_MAX) -> List[List[Dict]]:
    """Back-compat alias for :func:`_split_evenly` (identical behaviour)."""
    return _split_evenly(items, target_max)


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
    the category has more than AI_BATCH_MAX variants. On a timeout the batch
    is automatically halved and retried (down to AI_BATCH_MIN), so a single
    slow response no longer kills the entire category's output."""

    def _try_one(batch: List[Dict], idx: Optional[int] = None,
                 total: Optional[int] = None) -> Tuple[str, bool]:
        """Attempt one AI call. On TimeoutError, halve the batch and recurse
        (up to AI_BATCH_MIN). Returns (text, had_failure)."""
        try:
            t0 = time.time()
            text = call_ollama(
                _build_category_prompt(category, batch, idx, total),
                model=model,
            )
            log(f"      ✓ {len(batch)} variants in {time.time()-t0:.0f}s")
            return text, False
        except TimeoutError as te:
            if len(batch) <= AI_BATCH_MIN:
                log(f"    WARNING: {len(batch)}-variant batch timed out even at "
                    f"minimum size — recording for --retry-failed")
                return f"*AI analysis timed out for this batch: {te}*", True
            half = max(AI_BATCH_MIN, len(batch) // 2)
            log(f"    Batch of {len(batch)} timed out; auto-splitting to {half} and retrying")
            sub = _split_evenly(batch, half)
            sub_parts: List[str] = []
            any_fail = False
            for j, sb in enumerate(sub, start=1):
                t, f = _try_one(sb, j, len(sub))
                sub_parts.append(t); any_fail |= f
            return "\n\n".join(sub_parts), any_fail
        except Exception as e:
            log(f"    WARNING: AI call failed for {category}: {e}")
            return f"*AI analysis unavailable: {e}*", True

    batches = _split_evenly(snps, AI_BATCH_MAX)
    if len(batches) == 1:
        return _try_one(batches[0])

    sizes = ", ".join(str(len(b)) for b in batches)
    log(f"    Splitting {len(snps)} variants into {len(batches)} sub-batches ({sizes})")
    parts: List[str] = []
    had_failure = False
    offset = 0
    for i, batch in enumerate(batches, start=1):
        first = offset + 1
        last = offset + len(batch)
        log(f"    Batch {i}/{len(batches)} (variants {first}-{last}) ...")
        text, fail = _try_one(batch, i, len(batches))
        parts.append(
            f"**Batch {i} of {len(batches)} (variants {first}–{last})**\n\n{text}"
        )
        had_failure |= fail
        offset += len(batch)
        time.sleep(0.3)

    return "\n\n---\n\n".join(parts), had_failure


def _build_module_context(
    bloodwork_result: Optional[Dict] = None,
    immunogenetics_result: Optional[Dict] = None,
    neurochemistry_result: Optional[Dict] = None,
    addiction_genetics_result: Optional[Dict] = None,
    deep_ancestry_result: Optional[Dict] = None,
    blood_type_result: Optional[Dict] = None,
    family_planning_result: Optional[Dict] = None,
    polygenic_traits_result: Optional[Dict] = None,
    environmental_optimization_result: Optional[Dict] = None,
    holistic_synthesis_result: Optional[Dict] = None,
    detox_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
    clinical_variants_result: Optional[Dict] = None,
    y_result: Optional[Dict] = None,
    mt_result: Optional[Dict] = None,
    max_chars: int = 7000,
    per_section_chars: int = 900,
) -> str:
    """Compact, per-section-budgeted digest of the higher-order analysis
    modules — the 'whole person' beyond the tier-1 SNP category findings.

    Pure and None-safe: every section is optional and independently truncated
    to ``per_section_chars`` so no single module can crowd the others out, then
    the whole digest is capped at ``max_chars``. Shared by the executive
    summary, the cross-category synthesis, and the interactive chat assistant
    so all three reason over the same rich context.
    """
    sections: List[Tuple[str, str]] = []

    def _add(title: str, lines: List[str]) -> None:
        lines = [str(x) for x in lines if x]
        if not lines:
            return
        body = "\n".join(lines)
        if len(body) > per_section_chars:
            body = body[:per_section_chars].rstrip() + " …"
        sections.append((title, body))

    # ── Clinical variants (ClinVar P/LP — highest clinical weight, goes first) ──
    try:
        clv = clinical_variants_result or {}
        if clv.get("available") and clv.get("n_plp", 0) > 0:
            ls = [f"{clv['n_plp']} ClinVar pathogenic/likely-pathogenic variant(s) "
                  f"matched: {clv.get('n_actionable',0)} actionable, "
                  f"{clv.get('n_carrier',0)} carrier, {clv.get('n_affected',0)} "
                  f"affected-consistent. (Screening, not diagnosis.)"]
            for f in (clv.get("findings") or [])[:6]:
                ls.append(f"- {f.get('gene')} [{f.get('category')}] "
                          f"{f.get('condition','')[:60]} ({f.get('stars')}★, {f.get('zygosity')})")
            _add("CLINICAL VARIANTS (ClinVar)", ls)
    except Exception:
        pass

    # ── Holistic synthesis (leverage score — highest-signal framing) ──
    try:
        hs = holistic_synthesis_result or {}
        gl = hs.get("genome_leverage") or {}
        if gl.get("tier"):
            ls = [f"Genome Leverage Score: {gl.get('score')}/100 ({gl.get('tier')}). "
                  f"{gl.get('narrative','')}"]
            for ins in (hs.get("insights") or [])[:5]:
                ls.append(f"- {ins.get('title','')}: {ins.get('explanation','')[:160]}")
            _add("HOLISTIC SYNTHESIS", ls)
    except Exception:
        pass

    # ── Bloodwork (biological age + flagged markers) ──
    try:
        bw = bloodwork_result or {}
        clin = bw.get("clinical") or {}
        adv = clin.get("advanced") or {}
        bio = adv.get("biological_age") or {}
        ls = []
        if bio.get("phenoage") is not None:
            accel = bio.get("accel")
            ls.append(f"PhenoAge biological age {bio.get('phenoage'):.1f} "
                      f"(accel {accel:+.1f} yr)" if accel is not None
                      else f"PhenoAge {bio.get('phenoage'):.1f}")
        prevent = adv.get("prevent") or {}
        if prevent.get("risk10_total") is not None:
            ls.append(f"PREVENT 10-yr CVD risk ~{prevent['risk10_total']:.1f}%")
        flags = clin.get("flags") or []
        for f in flags[:8]:
            nm = f.get("name") or f.get("marker") or "?"
            val = f.get("value")
            note = (f.get("note") or f.get("interpretation") or "")[:80]
            ls.append(f"- FLAG {nm}={val}: {note}")
        _add("BLOODWORK (measured labs)", ls)
    except Exception:
        pass

    # ── Immunogenetics ──
    try:
        ig = immunogenetics_result or {}
        if ig.get("available"):
            ls = [f"{ig.get('n_protective',0)} protective / "
                  f"{ig.get('n_susceptible',0)} susceptible pathogen findings."]
            for h in (ig.get("headlines") or [])[:4]:
                ls.append(f"- PROTECTIVE {h.get('name','')} ({h.get('gene','')} "
                          f"{h.get('genotype','')}): {h.get('verdict','')}")
            for f in [x for x in (ig.get('findings') or []) if x.get('impact')=='susceptible'][:3]:
                ls.append(f"- SUSCEPTIBLE {f.get('name','')}: {f.get('verdict','')}")
            _add("IMMUNOGENETICS", ls)
    except Exception:
        pass

    # ── Neurochemistry ──
    try:
        nc = (neurochemistry_result or {}).get("composite") or {}
        if nc:
            ls = [f"COMT {nc.get('comt_class')} · MAOA {nc.get('maoa_class')} · "
                  f"BDNF {nc.get('bdnf_class')}.",
                  nc.get("stress_response_profile", ""),
                  f"Stimulants: {nc.get('stimulant_response','')[:120]}",
                  f"Caffeine: {nc.get('caffeine_protocol','')[:120]}"]
            _add("NEUROCHEMISTRY", ls)
    except Exception:
        pass

    # ── Addiction genetics ──
    try:
        ag = (addiction_genetics_result or {}).get("composite") or {}
        if ag:
            ls = [f"Alcohol tier: {ag.get('alcohol_tier')}. Overall: {ag.get('overall_tier')}."]
            for f in (ag.get("clinical_flags") or [])[:4]:
                ls.append(f"- {f.get('title','')}: {f.get('text','')[:100]}")
            _add("ADDICTION GENETICS", ls)
    except Exception:
        pass

    # ── Ancestry + haplogroups ──
    try:
        ls = []
        anc = ancestry_result or {}
        props = anc.get("proportions") or {}
        if props:
            top = sorted(props.items(), key=lambda kv: -kv[1])[:3]
            ls.append("Autosomal: " + ", ".join(f"{k} {v*100:.0f}%" for k, v in top))
        cc = anc.get("haplogroup_crosscheck") or {}
        if cc.get("verdict"):
            ls.append(f"Lineage cross-check: {cc['verdict']}")
        if (y_result or {}).get("terminal_haplogroup"):
            ls.append(f"Y-DNA {y_result['terminal_haplogroup']}")
        if (mt_result or {}).get("haplogroup"):
            ls.append(f"mtDNA {mt_result['haplogroup']}")
        _add("ANCESTRY & LINEAGE", ls)
    except Exception:
        pass

    # ── Deep ancestry ──
    try:
        dp = deep_ancestry_result or {}
        ls = []
        n = dp.get("neanderthal") or {}
        if n.get("available"):
            ls.append(f"Neanderthal ~{n.get('approx_pct')}% ({n.get('tier')})")
        ap = dp.get("ancient_populations") or {}
        if ap.get("available"):
            ls.append("Ancient-pop: " + ", ".join(
                f"{p['short']} {p['affinity']*100:.0f}%" for p in ap.get("populations", [])[:3]))
        ea = dp.get("european_axis") or {}
        if ea.get("available"):
            ls.append(f"N-S axis: {ea.get('lean')} (index {ea.get('index')})")
        _add("DEEP ANCESTRY", ls)
    except Exception:
        pass

    # ── Blood type ──
    try:
        bt = blood_type_result or {}
        if bt.get("available"):
            _add("BLOOD TYPE", [
                f"{bt.get('combined') or '—'} — ABO {(bt.get('abo') or {}).get('phenotype')} "
                f"(genotype {(bt.get('abo') or {}).get('genotype')}), "
                f"Rh {(bt.get('rhd') or {}).get('status')}. "
                f"FUT2: {(bt.get('secretor') or {}).get('secretor_status')}."])
    except Exception:
        pass

    # ── Family planning ──
    try:
        fp = family_planning_result or {}
        if fp.get("available"):
            _add("FAMILY PLANNING", [fp.get("summary", "")[:per_section_chars]])
    except Exception:
        pass

    # ── Trait genetics ──
    try:
        pt = polygenic_traits_result or {}
        if pt.get("available"):
            ls = [f"{f.get('trait')}: {f.get('call')}" for f in (pt.get("findings") or [])[:8]]
            _add("TRAIT GENETICS", ls)
    except Exception:
        pass

    # ── Environmental optimization ──
    try:
        eo = environmental_optimization_result or {}
        if eo.get("available"):
            ls = []
            if eo.get("circadian"):
                ls.append(f"Chronotype: {eo['circadian'].get('lean')}")
            if eo.get("exercise"):
                ls.append(f"Exercise fit: {eo['exercise'].get('lean')}")
            if eo.get("vitamin_d"):
                ls.append(f"Vitamin-D: {eo['vitamin_d'].get('tendency')}")
            _add("ENVIRONMENTAL OPTIMIZATION", ls)
    except Exception:
        pass

    # ── Detox ──
    try:
        dt = detox_result or {}
        tier = dt.get("smoke_resilience_tier") or dt.get("tier")
        if tier:
            _add("DETOXIFICATION", [f"Smoke/xenobiotic resilience tier: {tier}"])
    except Exception:
        pass

    if not sections:
        return ""

    out = "\n\n".join(f"[{title}]\n{body}" for title, body in sections)
    if len(out) > max_chars:
        out = out[:max_chars].rstrip() + "\n… [context truncated to fit]"
    return out


# Per-module AI ("AI on all tiers"): each new module section gets its own local
# AI interpretation. Keys are the report section anchor ids so the renderer can
# attach each interpretation to the right section.
_MODULE_AI_SPECS: Dict[str, Dict[str, str]] = {
    "clinical-variants": {
        "kwarg": "clinical_variants_result", "title": "Clinical Variants (ClinVar)",
        "focus": "Explain each pathogenic/likely-pathogenic finding plainly — what "
                 "the gene/condition is, what carrier vs affected vs actionable "
                 "means here, and the ClinVar star confidence. Be careful and "
                 "non-alarmist: stress this is a screen needing accredited-lab "
                 "confirmation + genetic counseling, and that absence of findings "
                 "is not absence of risk."},
    "holistic-synthesis": {
        "kwarg": "holistic_synthesis_result", "title": "Holistic Synthesis",
        "focus": "Explain what the Genome Leverage Score means for how much this "
                 "person's choices (vs their genes) drive their long-term outcome, "
                 "and translate the top cross-panel insights into a short priority order."},
    "immunogenetics": {
        "kwarg": "immunogenetics_result", "title": "Immunogenetics",
        "focus": "Explain the headline protective findings and what they mean day "
                 "to day, then any susceptibilities and their vaccination/behaviour "
                 "implications. Note how notable the overall combination is."},
    "neurochemistry": {
        "kwarg": "neurochemistry_result", "title": "Neurochemistry",
        "focus": "Translate the COMT/MAOA/BDNF profile into practical implications "
                 "for stress, focus, stimulant/caffeine response, and learning — as "
                 "tendencies, never as a fixed verdict on who they are."},
    "addiction-genetics": {
        "kwarg": "addiction_genetics_result", "title": "Addiction Genetics",
        "focus": "Explain the alcohol and overall susceptibility tiers plainly, and "
                 "surface the clinically-useful flags (never-smoke, naltrexone "
                 "response, opioid dosing). Stress that behaviour dominates."},
    "deep-ancestry": {
        "kwarg": "deep_ancestry_result", "title": "Deep Ancestry",
        "focus": "Tell the story of the Neanderthal, ancient-population, and "
                 "North-South axis findings — where these ancestral components came "
                 "from and what they say about deep origins."},
    "blood-type": {
        "kwarg": "blood_type_result", "title": "Blood Type",
        "focus": "Explain the predicted ABO/Rh type and secretor status, what the "
                 "hidden allele means for children, and the confidence caveats."},
    "family-planning": {
        "kwarg": "family_planning_result", "title": "Family Planning",
        "focus": "Explain the reproductive-relevant findings for future children — "
                 "keeping transmission separate from disease probability — and what, "
                 "if anything, is worth discussing with a partner or counselor."},
    "polygenic-traits": {
        "kwarg": "polygenic_traits_result", "title": "Trait Genetics",
        "focus": "Interpret the single-variant trait calls (taste, chronotype, "
                 "appearance) in a light, accurate way; reinforce why no polygenic "
                 "score is given for height/cognition/personality."},
    "environmental-optimization": {
        "kwarg": "environmental_optimization_result", "title": "Environmental Optimization",
        "focus": "Turn the chronotype, exercise-modality, and vitamin-D findings into "
                 "a concise, prioritized set of concrete behavioural recommendations."},
}


def ai_interpret_modules(model: str, log=None, **results) -> Dict[str, str]:
    """AI-interpret each higher-order module ('AI on all tiers'). Returns a dict
    keyed by report section-id → AI interpretation text. Each call is small and
    fully local; per-module failures are logged and skipped, never fatal."""
    _log = log or (lambda *a, **k: None)
    out: Dict[str, str] = {}
    for section_id, spec in _MODULE_AI_SPECS.items():
        res = results.get(spec["kwarg"])
        if not res:
            continue
        digest = _build_module_context(**{spec["kwarg"]: res},
                                       per_section_chars=1400, max_chars=2200)
        if not digest:
            continue
        prompt = (
            f"You are a genetic counselor. Interpret the following "
            f"{spec['title']} analysis for this individual. {spec['focus']}\n\n"
            f"{digest}\n\n"
            "Write 2–4 short, plain-language paragraphs. Be specific and grounded "
            "in the data shown; distinguish strong findings from weak tendencies. "
            "Do not just restate the numbers — explain what they mean and what to "
            "do. Educational only, not a medical diagnosis."
        )
        try:
            _log(f"  AI interpreting module: {spec['title']} ...")
            txt = call_ollama(prompt, model=model, num_ctx=AI_NUM_CTX,
                              num_predict=1100, think=False)
            if txt and txt.strip():
                out[section_id] = txt.strip()
        except Exception as e:
            _log(f"  WARNING: module AI for {section_id} failed: {e}")
    return out


def tier2_analysis(
    tier1_results: List[Dict], apoe_genotype: Optional[str], model: str,
    module_context: str = "",
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

    module_block = (
        "\n\nBeyond the SNP-level findings above, the report computed these "
        "higher-order analyses. Weave the most important of them into your "
        "summary — do not ignore them:\n\n" + module_context
        if module_context else "")

    exec_prompt = (
        "You are a genetic counselor providing an executive summary of a "
        "comprehensive DNA analysis. Based on ALL the findings below — the "
        "SNP-level findings AND the higher-order module analyses — write a "
        "5–7 paragraph executive summary that reads like one coherent picture "
        "of this whole person, not a list.\n\n"
        "Key SNP-level findings:\n"
        + "\n".join(summary_lines)
        + f"\n\nTotal variants analyzed: {len(tier1_results)}\n"
        f"Categories covered: "
        f"{', '.join(sorted(set(r['category'] for r in tier1_results)))}\n"
        + module_block
        + "\n\nYour summary should:\n"
        "1. Open with the single most important through-line about this person "
        "(use the Genome Leverage framing if present).\n"
        "2. Integrate genetics WITH measured labs where both are present "
        "(e.g. an APOE genotype next to the actual LDL; a PhenoAge result).\n"
        "3. Identify the top 3–5 highest-priority, highest-leverage actions.\n"
        "4. Note cross-domain patterns that only appear when modules combine "
        "(e.g. immunogenetics × inflammation labs; neurochemistry × addiction; "
        "ancestry-adapted diet × metabolic markers).\n"
        "5. Be balanced and grounded — neither alarmist nor dismissive; "
        "distinguish deterministic findings from probabilistic tendencies.\n"
        "6. End with clear next steps: what to prioritize, monitor, or discuss "
        "with a physician.\n\n"
        "Use plain language for a scientifically literate non-specialist. "
        "Educational information only, not a medical diagnosis."
    )

    try:
        # Larger context now that the exec summary reasons over every module.
        exec_summary = call_ollama(exec_prompt, model=model,
                                   num_ctx=16384, num_predict=2048, think=False)
    except Exception as e:
        log(f"  WARNING: Executive summary failed: {e}")
        exec_summary = f"*Executive summary unavailable: {e}*"

    return ai_results, exec_summary, failed_categories


def cross_category_synthesis(
    tier1_results: List[Dict], model: str, module_context: str = "",
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

    if not cat_blocks and not module_context:
        return "*No risk-carrying or cross-referenced variants were found, so no cross-category synthesis was generated.*"

    module_block = (
        "\n\nHigher-order module analyses (reason across these AND the "
        "categories above — the richest interactions bridge the two):\n\n"
        + module_context if module_context else "")

    prompt = (
        "You are a genetic counselor synthesising findings ACROSS multiple "
        "domains of a comprehensive DNA report. Per-category interpretations "
        "have already been generated separately. Your job here is different: "
        "identify how findings from DIFFERENT domains COMPOUND or INTERACT.\n\n"
        "Findings by category (risk-carrying or cross-referenced variants only):\n\n"
        + "\n\n".join(cat_blocks)
        + "\n\nVariants that already have explicit cross-category bridges in the "
        "database:\n"
        + ("\n".join(xref_lines) if xref_lines else "  (none)")
        + module_block
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
        "  • Insulin resistance × Lipids × Inflammation (metabolic syndrome)\n"
        "  • Immunogenetics × measured inflammation labs (e.g. FUT2 non-secretor "
        "microbiome baseline vs actual hs-CRP)\n"
        "  • Neurochemistry × Addiction (dopamine/COMT/reward × alcohol/nicotine tiers)\n"
        "  • Ancestry-adapted diet (LCT / steppe / farmer) × metabolic labs\n"
        "  • Genome Leverage tier × current lab trajectory (how much environment "
        "dominates this person's outcome)\n\n"
        "Only include patterns that are actually relevant to this person's data. "
        "Be specific — name the variants/modules, name the action.\n\n"
        "Reminder: educational, not diagnostic."
    )

    # This is the largest single prompt in the pipeline — a full risk-variant
    # snapshot across every category. Cap prompt size and give it a bigger
    # context window than per-category calls (which use the default 4096).
    if len(prompt) > 48_000:
        prompt = prompt[:48_000] + "\n\n[...truncated to fit context]"
    try:
        return call_ollama(prompt, model=model,
                           num_ctx=16384, num_predict=2560, think=False)
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

    # These renderer helpers are re-exported via this module's __getattr__ shim
    # for external callers, but bare-name lookups inside this function do not
    # trigger module __getattr__ — import them explicitly to avoid NameError.
    from renderers import _cat_id, _patch_ai_section_in_html, md_to_html

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
