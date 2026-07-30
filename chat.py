"""
Interactive Q&A Chat Mode (V25 — deep mode)
===========================================

After ``analyze.py`` runs the standard pipeline, ``--chat`` drops you into a
REPL where you can ask plain-language questions about your genetic findings.

V25 upgrades:
  • Deep-mode system prompt with a structured answering rubric.
  • Streaming responses (tokens print as the model generates them).
  • Vastly richer context — bloodwork (PhenoAge, PREVENT, flagged markers),
    detox/oxidative panel, longevity variants, personal economics, PRS, PGx,
    carrier, ancestry, interactions.
  • Follow-up suggestions after every answer, tailored to the topic.
  • REPL commands: /help, /context, /save <path>, /deep, /brief,
    /model <name>, /topic <name>, /reset, /new.

Type 'exit', 'quit', ':q', or Ctrl-D to leave.
"""

from __future__ import annotations

import datetime as _dt
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional


SCRIPT_DIR = Path(__file__).parent
DEFAULT_TIER1 = SCRIPT_DIR / "tier1_results.json"


# ── Context assembly ─────────────────────────────────────────────────────────

def _summarise_context(
    tier1_path: Path,
    pgx_result: Optional[Dict] = None,
    prs_result: Optional[Dict] = None,
    traits_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
    bloodwork_result: Optional[Dict] = None,
    detox_result: Optional[Dict] = None,
    economics_result: Optional[Dict] = None,
    personal_econ_result: Optional[Dict] = None,
    urologic_result: Optional[Dict] = None,
    deep_ancestry_result: Optional[Dict] = None,
    blood_type_result: Optional[Dict] = None,
    immunogenetics_result: Optional[Dict] = None,
    neurochemistry_result: Optional[Dict] = None,
    addiction_genetics_result: Optional[Dict] = None,
    holistic_synthesis_result: Optional[Dict] = None,
    clinical_variants_result: Optional[Dict] = None,
    novel_variants_result: Optional[Dict] = None,
    polygenic_traits_result: Optional[Dict] = None,
    environmental_optimization_result: Optional[Dict] = None,
    family_planning_result: Optional[Dict] = None,
    y_result: Optional[Dict] = None,
    mt_result: Optional[Dict] = None,
) -> str:
    """Build a rich context block describing the user's genetic + phenotypic
    state. Ordered so the highest-signal / most-actionable items come first."""
    lines: List[str] = []

    # Header ─────────────────────────────────────────────────────────────────
    if tier1_path.exists():
        try:
            data = json.loads(tier1_path.read_text())
            lines.append(f"REPORT VERSION: {data.get('report_version','—')}")
            lines.append(f"APOE genotype: {data.get('apoe_genotype') or 'not determined'}")
            ydna = data.get("y_haplogroup") or (y_result or {}).get("terminal_haplogroup")
            mtdna = data.get("mt_haplogroup") or (mt_result or {}).get("haplogroup")
            lines.append(f"Y-DNA haplogroup: {ydna or 'not determined'}")
            lines.append(f"mtDNA haplogroup: {mtdna or 'not determined'}")
            lines.append(f"Data QC grade: {data.get('qc_grade') or '—'}")
            lines.append(f"Curated variants matched: {data.get('total_matched', 0)}")

            variants = data.get("variants", [])
            risk_carrying = [v for v in variants if v.get("risk_copies", 0) > 0]
            high_sig = [v for v in risk_carrying if v.get("significance") == "high"]
            mod_sig = [v for v in risk_carrying if v.get("significance") == "moderate"]
            if high_sig:
                lines.append("\nHIGH-significance risk findings:")
                for v in high_sig[:25]:
                    gt = v.get("my_genotype", "?")
                    lines.append(f"  - {v['gene']} {v['variant_name']} "
                                 f"[{gt}, {v['risk_copies']} risk allele(s)] · {v['category']}")
                    if v.get("summary"):
                        lines.append(f"      {v['summary'][:180]}")
            if mod_sig:
                lines.append(f"\nModerate-significance risk findings ({len(mod_sig)}):")
                for v in mod_sig[:25]:
                    gt = v.get("my_genotype", "?")
                    lines.append(f"  - {v['gene']} {v['variant_name']} "
                                 f"[{gt}, {v['risk_copies']} risk allele(s)] · {v['category']}")

            prs_summary = data.get("prs_summary") or {}
            if prs_summary:
                lines.append("\nPolygenic Risk Score tiers:")
                for name, p in prs_summary.items():
                    tier = (p or {}).get("tier"); pct = (p or {}).get("percentile")
                    if tier:
                        lines.append(f"  - {name}: {tier} ({pct}th percentile)")

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

    # Traits / carriers / interactions / ancestry (unchanged)
    if traits_result and traits_result.get("predictions"):
        lines.append("\nTrait predictions (high/moderate confidence):")
        for t in traits_result["predictions"]:
            if t.get("confidence") in ("high", "moderate"):
                lines.append(f"  - {t.get('trait')}: {t.get('result')}")

    if carrier_result:
        n_aff = carrier_result.get("n_affected", 0)
        n_car = carrier_result.get("n_carriers", 0)
        if n_aff or n_car:
            lines.append(f"\nCarrier status: {n_aff} affected, {n_car} carrier(s)")
            for c in (carrier_result.get("carriers") or [])[:10]:
                lines.append(f"  - {c.get('gene','?')} {c.get('condition','')} "
                             f"({c.get('inheritance','')})")

    if interactions_result and interactions_result.get("findings"):
        lines.append("\nCompound / variant-interaction findings:")
        for f in interactions_result["findings"][:10]:
            lines.append(f"  - [{(f.get('severity') or '').upper()}] {f.get('title','')}")

    if ancestry_result and ancestry_result.get("proportions"):
        props = ancestry_result["proportions"]
        top = ", ".join(f"{sp} {p*100:.0f}%"
                        for sp, p in sorted(props.items(), key=lambda x: -x[1])[:3])
        conf = ancestry_result.get("confidence", "—")
        lines.append(f"\nEstimated ancestry: {top} (confidence: {conf})")
        cc = ancestry_result.get("haplogroup_crosscheck") or {}
        if cc.get("verdict"):
            lines.append(f"Y/mtDNA lineage cross-check: {cc['verdict']} — {cc.get('summary','')}")

    # ── Bloodwork (V6.x) — new richer context ────────────────────────────────
    if bloodwork_result:
        clinical = bloodwork_result.get("clinical") or {}
        adv = clinical.get("advanced") or {}
        bio = adv.get("biological_age") or {}
        if bio:
            lines.append("\nBLOOD WORK — biological age & risk:")
            lines.append(f"  - PhenoAge: {bio.get('phenoage')} yr "
                         f"({bio.get('accel'):+g} vs chronological "
                         f"{bio.get('chronological')})")
            if bio.get("mortality_10yr_pct") is not None:
                lines.append(f"  - Modeled 10-yr mortality risk: {bio['mortality_10yr_pct']}%")
            if bio.get("recoverable_years"):
                lines.append(f"  - Recoverable biological years if all "
                             f"markers optimized: {bio['recoverable_years']}")
            for lev in (bio.get("levers") or [])[:5]:
                lines.append(f"    · lever: {lev['marker']} "
                             f"{lev['current']} → {lev['ideal']} = −{lev['years_cost']} yr")

        for idx in adv.get("indices", []):
            if idx.get("id") == "prevent_ascvd":
                lines.append(f"  - PREVENT 10-yr ASCVD risk: {idx['value']}% "
                             f"({idx.get('status_label','')})")
                break

        flags = clinical.get("flags") or []
        if flags:
            lines.append("  Flagged biomarkers (out-of-range):")
            for f in flags[:12]:
                geno = f" · {f['genotype_note']}" if f.get("genotype_note") else ""
                lines.append(f"    - {f['name']} = {f['value']} {f['unit']} "
                             f"({f['status_label']}){geno}")

        gl = adv.get("genetic_longevity") or {}
        if gl.get("variants"):
            fav = gl.get("n_favorable", 0); adv_ct = gl.get("n_adverse", 0)
            lines.append(f"  Longevity variants: {fav} favorable / {adv_ct} adverse")
            for v in gl["variants"][:8]:
                tag = "FAV" if v["favorable"] else "adv"
                lines.append(f"    - {v['gene']} {v['rsid']} [{v['genotype']}] {tag}: {v['label']}")

    # ── Blood type (ABO + Rh + secretor) ─────────────────────────────────────
    if blood_type_result and blood_type_result.get("available"):
        bt = blood_type_result
        lines.append(f"\nBLOOD TYPE (predicted): {bt.get('combined') or '—'} "
                     f"(ABO {bt['abo']['phenotype']}, genotype {bt['abo']['genotype']}; "
                     f"Rh {bt['rhd']['status']}). "
                     f"FUT2: {bt['secretor'].get('secretor_status','?')}.")
        if bt["abo"].get("carries_hidden_O"):
            lines.append("  · Carries a hidden O allele (recessive)")

    # ── Deep ancestry (Neanderthal + ancient populations + N/S axis) ─────────
    if deep_ancestry_result and deep_ancestry_result.get("available"):
        n = deep_ancestry_result.get("neanderthal") or {}
        if n.get("available"):
            lines.append(f"\nDEEP ANCESTRY — Neanderthal affinity ~{n.get('approx_pct','?')}% "
                         f"({n.get('tier','?')}); {n.get('n_carrying',0)} of "
                         f"{n.get('n_typed',0)} Neanderthal-tagged loci carrying.")
            for v in n.get("variants", [])[:5]:
                if v["n_alleles"] > 0:
                    lines.append(f"  · Neanderthal allele: {v['gene']} {v['rsid']} "
                                 f"[{v['genotype']}, dose {v['n_alleles']}] — {v['trait'][:100]}")
        ap = deep_ancestry_result.get("ancient_populations") or {}
        if ap.get("available"):
            top3 = ap["populations"][:3]
            lines.append("Ancient-population affinity (Yamnaya/EEF/WHG): " +
                         ", ".join(f"{p['short']} {p['affinity']*100:.0f}%" for p in top3))
        ea = deep_ancestry_result.get("european_axis") or {}
        if ea.get("available"):
            lines.append(f"N-S Europe axis: {ea['lean']} (index {ea['index']})")
        tl = deep_ancestry_result.get("haplogroup_timeline") or {}
        if tl.get("y"):
            lines.append(f"Y-DNA timeline: {tl['y']['haplogroup']} — TMRCA ~"
                         f"{tl['y']['tmrca_kya']} kya · {tl['y']['origin']}")
        if tl.get("mt"):
            lines.append(f"mtDNA timeline: {tl['mt']['haplogroup']} — TMRCA ~"
                         f"{tl['mt']['tmrca_kya']} kya · {tl['mt']['origin']}")

    # ── Urologic panel ───────────────────────────────────────────────────────
    if urologic_result and urologic_result.get("available"):
        lines.append(f"\nUROLOGIC PANEL: {urologic_result['n_findings']} findings "
                     f"({urologic_result.get('n_flagged',0)} flagged) across "
                     f"{len(urologic_result.get('categories',[]))} sub-panels "
                     f"(OAB, BPH/prostate, kidney stones, testicular, androgen).")
        for f in urologic_result.get("findings", []):
            if f.get("impact") in ("higher-load", "reduced", "reduced-clearance"):
                lines.append(f"  - {f['gene']} {f['trait'][:60]} "
                             f"[{f['genotype']}, {f['impact']}]: {f['result'][:140]}...")

    # ── Detox / smoke resilience ─────────────────────────────────────────────
    if detox_result and detox_result.get("available"):
        sr = detox_result.get("smoke_resilience") or {}
        if sr:
            lines.append(f"\nDetox / smoke resilience: {sr.get('tier','?')} "
                         f"(activation {sr.get('activation_hits',0)}, "
                         f"clearance-deficit {sr.get('clearance_deficit_hits',0)}, "
                         f"antioxidant-deficit {sr.get('antioxidant_deficit_hits',0)})")
            if sr.get("activate_but_dont_clear"):
                lines.append("  ⚠ Activate-but-don't-clear genotype pattern detected.")

    # ── Personal economics ──────────────────────────────────────────────────
    if personal_econ_result and personal_econ_result.get("available"):
        lines.append(f"\nModeled 10-yr net economic value of acting: "
                     f"${personal_econ_result.get('total_net'):,} "
                     f"(ROI {personal_econ_result.get('roi')}× vs analysis cost)")
        for i in (personal_econ_result.get("items") or [])[:5]:
            lines.append(f"  - {i['category']}: {i['finding']} · net ${i['net']:,} "
                         f"({i['confidence']} confidence)")

    # ── Shared higher-order module digest (immunogenetics, neurochemistry,
    #    addiction, holistic synthesis, trait genetics, env-opt, family planning)
    #    so chat reasons over the same rich context as the report AI. ──
    try:
        from analyze import _build_module_context
        digest = _build_module_context(
            bloodwork_result=bloodwork_result,
            immunogenetics_result=immunogenetics_result,
            neurochemistry_result=neurochemistry_result,
            addiction_genetics_result=addiction_genetics_result,
            deep_ancestry_result=deep_ancestry_result,
            blood_type_result=blood_type_result,
            family_planning_result=family_planning_result,
            polygenic_traits_result=polygenic_traits_result,
            environmental_optimization_result=environmental_optimization_result,
            holistic_synthesis_result=holistic_synthesis_result,
            clinical_variants_result=clinical_variants_result,
            novel_variants_result=novel_variants_result,
            detox_result=detox_result,
            ancestry_result=ancestry_result,
            y_result=y_result, mt_result=mt_result,
            max_chars=9000,
        )
        if digest:
            lines.append("\n=== HIGHER-ORDER MODULE ANALYSES ===\n" + digest)
    except Exception:
        pass

    return "\n".join(lines)


# ── Terminal rendering ───────────────────────────────────────────────────────

def _render_response(text: str) -> str:
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
    text = re.sub(r"\*\*(.+?)\*\*", "\033[1m\\1\033[0m", text)
    text = re.sub(r"\*([^*\n]+)\*", "\033[3m\\1\033[0m", text)
    lines = []
    for line in text.split("\n"):
        stripped = line.lstrip()
        if stripped.startswith(("- ", "* ", "• ")):
            indent = len(line) - len(stripped)
            lines.append(" " * indent + "\033[36m•\033[0m " + stripped[2:])
        elif re.match(r"^#{1,6}\s+", stripped):
            n = len(stripped) - len(stripped.lstrip("#"))
            lines.append("\033[1;36m" + stripped.lstrip("# ").strip() + "\033[0m"
                         if n <= 2 else "\033[1m" + stripped.lstrip("# ") + "\033[0m")
        else:
            lines.append(line)
    return "\n".join(lines)


def _stream_and_render(chunk_iter) -> str:
    """Stream text to the terminal (light bold/italic on-the-fly) and return
    the full concatenated content."""
    buf = []
    in_think = False
    pending_bold = False
    print("\033[2m─── Assistant ───\033[0m")
    for chunk in chunk_iter:
        buf.append(chunk)
        # Filter <think>...</think> blocks live
        s = chunk
        if "<think>" in s or in_think:
            in_think = True
            if "</think>" in s:
                s = s.split("</think>", 1)[1]
                in_think = False
            else:
                continue
        # Very light live styling: bold on ** pairs
        i = 0
        while i < len(s):
            if s[i:i+2] == "**":
                pending_bold = not pending_bold
                print("\033[1m" if pending_bold else "\033[0m", end="", flush=True)
                i += 2
                continue
            print(s[i], end="", flush=True)
            i += 1
    print("\033[0m")
    return "".join(buf)


# ── System prompts ───────────────────────────────────────────────────────────

_DEEP_RUBRIC = """
When answering ANY substantive question, structure your reply with clearly-
labeled sections (skip a section if it truly doesn't apply):

  1. **Direct answer** — the bottom line in 1–3 sentences.
  2. **What your data actually shows** — cite the specific gene(s), variant(s),
     genotype(s), risk-copy counts, PRS tier(s), or biomarker values from the
     GENETIC CONTEXT that drive the answer. Never invent findings.
  3. **Mechanism / biology** — explain WHY those findings lead to your conclusion,
     in plain language.
  4. **Personal fit** — reconcile the finding with the person's other results
     (e.g. an APOE ε4 carrier with high LDL vs an ε2 carrier; a slow metabolizer
     stacked with high-risk PRS; a MTHFR carrier with elevated homocysteine).
  5. **Concrete action plan** — specific, prioritised steps (dose, food, test,
     habit) they can act on this month. Prefer specificity over "eat healthy".
  6. **Uncertainty & limits** — what the evidence base can and can't say,
     effect sizes (small/moderate/large), and where results are chip-inferred
     vs directly measured.
  7. **When to see a clinician** — the specific situation, specialty, or test
     that warrants a professional (physician, board-certified genetic
     counselor, cardiologist, endocrinologist, etc.).
  8. **Follow-ups** — 2–3 sharp questions they could ask next to go deeper.

Length: DEFAULT to depth (aim for 500–1200 words), unless the user explicitly
asks for a short answer, in which case give 3–6 sentences with no rubric.
Never end without at least one concrete action step and a reminder that this
is educational, not medical advice.
"""

_BRIEF_RUBRIC = """
Answer concisely: 3–6 sentences, direct answer first, one concrete action,
one caveat. No headings. Still ground every claim in the user's specific
findings from GENETIC CONTEXT.
"""

_BASE_PROMPT = (
    "You are a senior clinical geneticist / genetic counselor and functional-"
    "medicine-literate researcher with deep expertise across pharmacogenomics "
    "(CPIC), polygenic risk, biological aging (PhenoAge/GrimAge), "
    "cardiovascular risk (AHA PREVENT), metabolic health, methylation, "
    "hormones, detox, and preventive medicine. You are having a private "
    "conversation with an individual about THEIR OWN DNA + blood-work "
    "analysis. Every claim must be grounded in the GENETIC CONTEXT block "
    "below — refer to specific genes, variants, genotypes, biomarker values, "
    "and PRS/PGx tiers by name. If context is silent on something, say so "
    "rather than inventing findings.\n\n"
    "You may be direct and technically precise (this user is scientifically "
    "literate). Use bullet lists and headings freely. Cite evidence tier "
    "(well-established / emerging / preliminary) for major claims. Bring in "
    "up-to-date knowledge — CPIC guidelines, PREVENT, ADA/AHA thresholds, "
    "cutting-edge longevity science — but never contradict the specific "
    "genotype/biomarker data in front of you.\n\n"
    "Medical decisions, prescriptions, and dose changes always belong with a "
    "licensed physician or board-certified genetic counselor — remind the "
    "user at the end of every substantive answer."
)


def system_prompt(mode: str = "deep") -> str:
    return _BASE_PROMPT + "\n\n" + (_BRIEF_RUBRIC if mode == "brief" else _DEEP_RUBRIC)


# ── Suggested / follow-up questions ──────────────────────────────────────────

SUGGESTED_QUESTIONS = [
    "What are my biggest health risks based on my genetics and blood work?",
    "Walk me through my biological age — what's dragging it up, what's my top lever?",
    "Which medications should I be cautious about, based on my PGx phenotype?",
    "How should I eat given my genome, blood work, and detox profile?",
    "What should I tell my primary care doctor at my next visit?",
    "How do my longevity variants (FOXO3, APOE, etc.) affect my long-term outlook?",
    "What supplement stack would give me the biggest personalised return?",
    "Explain my PREVENT cardiovascular risk and what to do about it.",
    "Should I see a genetic counselor? What would I bring?",
]


_TOPIC_HINTS = [
    (r"medication|drug|dose|adderall|vyvanse|ssri|statin|warfarin|clopidog|codeine|opioid",
     ["Which specific drugs should I discuss with my doctor?",
      "How would my CYP2D6/CYP2C19 phenotype change SSRI dosing?",
      "Any anesthesia-relevant PGx findings I should carry on a card?"]),
    (r"cardio|heart|cvd|ascvd|prevent|lipids|ldl|apob|blood pressure|bp\b",
     ["What targets should I aim for — LDL, ApoB, SBP?",
      "How would statin therapy change my 10-year ASCVD risk?",
      "Is Lp(a) worth measuring given my family history?"]),
    (r"biolog|phenoage|aging|longev|lifespan|healthspan|foxo",
     ["Which two markers would move my biological age most?",
      "How does my FOXO3/APOE combo actually change what I should do?",
      "What retest cadence makes sense to track biological age?"]),
    (r"diet|nutrition|eat|food|carb|fat|protein|caffeine|alcohol",
     ["A specific weekly meal pattern for my methylation + APOE profile?",
      "Caffeine dose ceiling given my CYP1A2 status?",
      "Alcohol threshold given my ADH1B/ALDH2 and liver panel?"]),
    (r"supplement|vitamin|mineral|nac|omega|magnes|methylfolate",
     ["Rank my top 5 supplements by expected impact.",
      "What doses and forms for methylfolate/B12/D given my MTHFR + labs?",
      "Any interactions between my current recs and prescription meds?"]),
    (r"detox|smoke|wildfire|toxin|metal|glutathione|nrf2|gst",
     ["A concrete protocol for the next wildfire-smoke event.",
      "Do I need heavy-metal testing given my ALAD/AS3MT genotype?",
      "How does my NRF2/GST profile change everyday exposure calls?"]),
    (r"ancestry|haplogroup|foxo|y-dna|mtdna|europe|african|asian",
     ["What does T1a1a T1a1a subclade actually tell me practically?",
      "Are any of my risk estimates less reliable due to ancestry?",
      "What genealogy tests would go deeper than a consumer chip?"]),
    (r"diabetes|glucose|hba1c|insulin|homa|tyg|prediab",
     ["A 90-day glucose-reduction plan I could realistically follow?",
      "Which of my genes make me more sensitive to refined carbs?",
      "Would CGM help me given my prediabetes flag?"]),
    (r"carrier|reproductive|family|fertility|pregnan",
     ["Which findings matter most before family planning?",
      "How would a partner's carrier profile change our risk?",
      "When would you recommend a fertility work-up on this profile?"]),
]


def suggest_followups(user_q: str, answer: str) -> List[str]:
    text = (user_q + " " + answer).lower()
    for pat, qs in _TOPIC_HINTS:
        if re.search(pat, text):
            return qs
    # Fallback: general deep-dive prompts
    return [
        "What would meaningfully change your recommendation here?",
        "Which single action would I get the most return on this month?",
        "How would you know in 3 months whether it's working?",
    ]


# ── REPL ─────────────────────────────────────────────────────────────────────

_HELP_TEXT = """\
\033[1mCommands:\033[0m
  \033[36m/help\033[0m              show this help
  \033[36m/deep\033[0m              switch to deep, structured answers (default)
  \033[36m/brief\033[0m             switch to short, focused answers
  \033[36m/model <name>\033[0m      change Ollama model (e.g. huihui_ai/qwen3-abliterated:30b-a3b)
  \033[36m/context\033[0m           print the genetic-context block sent to the model
  \033[36m/topic <name>\033[0m      biases the next answer toward a topic (e.g. pgx, cardio, detox, longevity)
  \033[36m/save <path>\033[0m       save the current conversation to a markdown file
  \033[36m/reset\033[0m /\033[36m/new\033[0m       clear conversation history (keep context)
  \033[36m/suggest\033[0m           print example questions
  \033[36mexit\033[0m / \033[36mquit\033[0m / \033[36m:q\033[0m  leave the chat
"""


def _save_transcript(history: List[Dict[str, str]], path: Path, model: str) -> None:
    stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# DNA Chat Transcript — {stamp}\n\n_Model:_ `{model}`\n\n---\n\n")
        for m in history:
            if m["role"] == "system":
                continue
            head = "🧑 You" if m["role"] == "user" else "🧬 Assistant"
            f.write(f"### {head}\n\n{m['content']}\n\n---\n\n")


def run_chat(
    model: str = "qwen3:14b",
    tier1_path: Optional[Path] = None,
    pgx_result: Optional[Dict] = None,
    prs_result: Optional[Dict] = None,
    traits_result: Optional[Dict] = None,
    carrier_result: Optional[Dict] = None,
    interactions_result: Optional[Dict] = None,
    ancestry_result: Optional[Dict] = None,
    bloodwork_result: Optional[Dict] = None,
    detox_result: Optional[Dict] = None,
    economics_result: Optional[Dict] = None,
    personal_econ_result: Optional[Dict] = None,
    urologic_result: Optional[Dict] = None,
    deep_ancestry_result: Optional[Dict] = None,
    blood_type_result: Optional[Dict] = None,
    immunogenetics_result: Optional[Dict] = None,
    neurochemistry_result: Optional[Dict] = None,
    addiction_genetics_result: Optional[Dict] = None,
    holistic_synthesis_result: Optional[Dict] = None,
    clinical_variants_result: Optional[Dict] = None,
    novel_variants_result: Optional[Dict] = None,
    polygenic_traits_result: Optional[Dict] = None,
    environmental_optimization_result: Optional[Dict] = None,
    family_planning_result: Optional[Dict] = None,
    y_result: Optional[Dict] = None,
    mt_result: Optional[Dict] = None,
) -> None:
    """Start an interactive chat session. Returns when the user exits."""
    import requests

    tier1_path = tier1_path or DEFAULT_TIER1
    context = _summarise_context(
        tier1_path, pgx_result=pgx_result, prs_result=prs_result,
        traits_result=traits_result, carrier_result=carrier_result,
        interactions_result=interactions_result, ancestry_result=ancestry_result,
        bloodwork_result=bloodwork_result, detox_result=detox_result,
        economics_result=economics_result, personal_econ_result=personal_econ_result,
        urologic_result=urologic_result, deep_ancestry_result=deep_ancestry_result,
        blood_type_result=blood_type_result,
        immunogenetics_result=immunogenetics_result,
        neurochemistry_result=neurochemistry_result,
        addiction_genetics_result=addiction_genetics_result,
        holistic_synthesis_result=holistic_synthesis_result,
        clinical_variants_result=clinical_variants_result,
        novel_variants_result=novel_variants_result,
        polygenic_traits_result=polygenic_traits_result,
        environmental_optimization_result=environmental_optimization_result,
        family_planning_result=family_planning_result,
        y_result=y_result, mt_result=mt_result,
    )
    mode = "deep"
    topic_bias = None

    print("\033[1m" + "─" * 72 + "\033[0m")
    print("\033[1mDNA Chat — Deep Mode\033[0m")
    print(f"Model: \033[36m{model}\033[0m (local Ollama, streaming) · "
          f"context: {len(context.splitlines())} lines")
    print("Type your question, or \033[36m/help\033[0m for commands. \033[36mexit\033[0m to leave.")
    print("\033[1m" + "─" * 72 + "\033[0m")

    def _system_msg() -> Dict[str, str]:
        prompt = system_prompt(mode)
        if topic_bias:
            prompt += (f"\n\nThe user is currently focused on the topic: "
                       f"**{topic_bias}**. Weight your answer toward that lens "
                       f"while remaining grounded in their data.")
        return {"role": "system",
                "content": prompt + "\n\n--- GENETIC CONTEXT ---\n" + context}

    history: List[Dict[str, str]] = [_system_msg()]

    while True:
        try:
            user_in = input("\033[1myou> \033[0m").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye.")
            break
        if not user_in:
            continue
        low = user_in.lower()
        if low in ("exit", "quit", ":q"):
            print("Goodbye.")
            break
        if low in ("help", "/help", "?"):
            print(_HELP_TEXT); continue
        if low in ("suggest", "/suggest", "suggestions"):
            print("\033[2mSuggested questions:\033[0m")
            for q in SUGGESTED_QUESTIONS: print(f"  · {q}")
            print(); continue
        if low in ("/context", "context"):
            print(context); continue
        if low in ("/deep",):
            mode = "deep"; history[0] = _system_msg()
            print("\033[2m(deep mode)\033[0m"); continue
        if low in ("/brief",):
            mode = "brief"; history[0] = _system_msg()
            print("\033[2m(brief mode)\033[0m"); continue
        if low.startswith("/model"):
            parts = user_in.split(maxsplit=1)
            if len(parts) == 2:
                model = parts[1].strip()
                print(f"\033[2m(model → {model})\033[0m")
            else:
                print("Usage: /model <ollama-model-name>")
            continue
        if low.startswith("/topic"):
            parts = user_in.split(maxsplit=1)
            topic_bias = parts[1].strip() if len(parts) == 2 else None
            history[0] = _system_msg()
            print(f"\033[2m(topic → {topic_bias or 'cleared'})\033[0m"); continue
        if low.startswith("/save"):
            parts = user_in.split(maxsplit=1)
            path = Path(parts[1] if len(parts) == 2 else "dna_chat.md").expanduser()
            _save_transcript(history, path, model)
            print(f"\033[2m(saved → {path})\033[0m"); continue
        if low in ("/reset", "/new"):
            history = [_system_msg()]
            print("\033[2m(conversation reset)\033[0m"); continue

        history.append({"role": "user", "content": user_in})

        # Stream response from Ollama
        try:
            payload = {
                "model": model, "messages": history, "stream": True,
                "options": {"temperature": 0.35, "num_ctx": 16384,
                            "num_predict": 2048 if mode == "deep" else 512},
            }
            t0 = time.time()
            with requests.post("http://localhost:11434/api/chat",
                               json=payload, stream=True, timeout=900) as resp:
                resp.raise_for_status()
                def _chunks():
                    for line in resp.iter_lines():
                        if not line: continue
                        try:
                            obj = json.loads(line)
                        except Exception:
                            continue
                        msg = obj.get("message") or {}
                        piece = msg.get("content", "")
                        if piece:
                            yield piece
                        if obj.get("done"):
                            break
                content = _stream_and_render(_chunks())
            elapsed = time.time() - t0
            print(f"\n\033[2m  ({elapsed:.1f}s)\033[0m")
            history.append({"role": "assistant", "content": content})
            # Follow-up prompts
            followups = suggest_followups(user_in, content)
            if followups:
                print("\033[2m  next-question ideas:\033[0m")
                for q in followups:
                    print(f"    \033[36m·\033[0m {q}")
            print()
        except requests.exceptions.ConnectionError:
            print("\n\033[31m  Cannot connect to Ollama. Start it: ollama serve\033[0m")
            break
        except Exception as e:
            print(f"\n\033[31m  Error: {e}\033[0m")
            if history and history[-1]["role"] == "user":
                history.pop()


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="DNA chat mode")
    ap.add_argument("--model", default="qwen3:14b")
    ap.add_argument("--tier1", default=str(DEFAULT_TIER1),
                    help="Path to tier1_results.json")
    args = ap.parse_args()
    run_chat(model=args.model, tier1_path=Path(args.tier1))
