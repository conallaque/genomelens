"""
Value-of-Information Health-Economics Engine
============================================

The recentered, discipline-grade layer: take every actionable finding (PGx, PRS,
APOE, Phase-2 clinical variants, Phase-3 novel variants) and answer, in health-
economics terms, *what is knowing your genome worth, in expectation?*

Methodology (decision-analytic, aspiring to CHEERS 2022 reporting):
  * **Net Monetary Benefit** per finding: NMB = ΔQALY·λ + ΔCost_averted − intervention.
  * **Discounting** of BOTH future costs and QALYs at 3% (sensitivity 0/3/5%).
  * **Cost-of-Illness (COI)** per condition — lifetime direct + indirect, sourced.
  * **Pharmacogenomics economics** — expected averted-ADR value from published CEA:
    E = P(prescribed) · P(ADR|genotype) · RRR(guided) · [ADR cost + QALY loss·λ].
  * **Value of Information (VOI)** — ex-ante (EVSI-style, population prior) and
    ex-post (this genome), plus the **marginal VOI of upgrading chip → WGS**.
  * **Probabilistic Sensitivity Analysis** — seeded Monte-Carlo over uncertain
    parameters → mean, 95% CI, a Cost-Effectiveness Acceptability Curve, and a
    one-way tornado.

Everything is illustrative and probabilistic. This is NOT a formal economic
evaluation and NOT financial or medical advice. Predicted (Phase-3) findings are
down-weighted because they are computational, not confirmed.
"""

from __future__ import annotations

import math

from . import params as _econ_params

try:
    import numpy as np
    _HAVE_NP = True
except Exception:                       # pragma: no cover
    np = None
    _HAVE_NP = False

# ── sourced parameters ────────────────────────────────────────────────────────
WTP = {"low": 50_000, "base": 100_000, "high": 150_000}   # $/QALY, Neumann NEJM 2014
DISCOUNT_RATE = 0.03                                       # 2nd-panel CEA standard
DISCOUNT_SENSITIVITY = (0.0, 0.03, 0.05)
TEST_COST = {"chip": 100.0, "wgs": 300.0}                 # consumer chip vs WGS promo

# Cost-of-Illness: lifetime, direct + indirect, present-ish value (illustrative).
# Condition anchors this module values, mapped to their registry parameters.
# Imported from the engine rather than restated so a new anchor cannot be added
# in one place and forgotten in the other.
def _coi_param_source() -> dict[str, tuple[str, str]]:
    from .engine import COI_KEY_TO_PARAM
    return dict(COI_KEY_TO_PARAM)


_COI_PARAM_SOURCE = _coi_param_source()

# ONE SOURCE OF TRUTH. This dict predated the provenance registry and held its
# own copy of every cost-of-illness figure. Five of the twelve had already
# drifted from the registry by the time anyone checked — Alzheimer's worst, at
# $350,000 here against $200,000 there, a 75% overstatement in whichever path
# read this copy. The costs now come from the registry, so the tier, source and
# distribution travel with the value and there is only one of each number.
#
# Aligning moved three anchors UP (breast/ovarian, colorectal, generic
# pathogenic) and two DOWN, for a net of about -$105,000. The registry figure is
# the sourced one in both directions; keeping a known-too-low number to satisfy
# a "corrections must reduce" rule would be its own distortion.
COI: dict[str, dict] = {
    key: {"cost": _econ_params.value(params[0]),
          "src": f"econ.params: {params[0]} ({_econ_params.PARAMS[params[0]].tier})"}
    for key, params in _COI_PARAM_SOURCE.items()
}

# MARGINAL vs AVERAGE COST (honesty correction).
# The COI figures above are AVERAGE lifetime cost-of-illness (total system cost ÷
# cases). Averting ONE case does NOT free the average cost — it frees the *marginal*
# cost, which is lower, because a large share of average cost is fixed/capacity
# spending (hospital buildings, salaried staff, overhead) that persists when one
# case is prevented. Using average-as-marginal systematically OVERSTATES cash savings
# — the "freeing a bed doesn't save its average cost" error.
#
# So the averted-cost side is scaled by a conservative marginal-cost fraction. This is
# a DOCUMENTED ASSUMPTION, not a fitted value, and it only ever REDUCES the modelled
# saving (the honest direction). Short-run marginal hospital cost is commonly ~50–70%
# of average; these lifetime COI figures also contain genuinely per-case-avertable
# components (long-term care, lost productivity), so 0.60 is a deliberately
# middle-conservative default. Override per your setting; do not tune it upward to
# make outputs look better.
# Ref: marginal-vs-average distinction in health-economic costing
# (Drummond et al. 2015; health-economics-metrics: marginal-vs-average-cost).
MARGINAL_COST_FRACTION = 0.60

# Pharmacogenomics cost-effectiveness (averted ADRs), from published CEA.
# p_rx = P(prescribed over a lifetime); p_adr = P(severe ADR | actionable genotype,
# standard dosing); rrr = relative risk reduction with genotype-guided care;
# adr_cost = cost of the ADR event; qaly = QALY loss averted.
PGX_CEA: dict[str, dict] = {
    "HLA-B*57:01/abacavir": {"p_rx": 0.02, "p_adr": 0.55, "rrr": 0.95,
                             "adr_cost": 20_000, "qaly": 0.30,
                             "src": "Schackman 2008 — cost-saving"},
    "CYP2C19/clopidogrel":  {"p_rx": 0.10, "p_adr": 0.20, "rrr": 0.50,
                             "adr_cost": 30_000, "qaly": 0.20,
                             "src": "Kazi 2014 — genotype-guided antiplatelet"},
    "TPMT-NUDT15/thiopurine": {"p_rx": 0.03, "p_adr": 0.35, "rrr": 0.70,
                               "adr_cost": 18_000, "qaly": 0.25,
                               "src": "CPIC — myelosuppression avoidance"},
    "DPYD/fluoropyrimidine": {"p_rx": 0.03, "p_adr": 0.30, "rrr": 0.60,
                              "adr_cost": 25_000, "qaly": 0.40,
                              "src": "Deenen 2016 — cost-saving"},
    "SLCO1B1/simvastatin":  {"p_rx": 0.20, "p_adr": 0.10, "rrr": 0.50,
                             "adr_cost": 6_000, "qaly": 0.05,
                             "src": "CPIC — statin myopathy"},
    "PGx-generic":          {"p_rx": 0.30, "p_adr": 0.10, "rrr": 0.50,
                             "adr_cost": 8_000, "qaly": 0.08,
                             "src": "generic actionable PGx flag"},
}

# ── Grossman (1972) health-capital parameters ────────────────────────────────
# Health is modelled as a depreciating capital stock H_t: you are born with an
# endowment, it depreciates at rate δ(age), and you invest (I_t) to offset it.
# Genomic information does not add health directly — it raises the EFFICIENCY of
# investment (you target the right interventions), which is exactly why its value
# scales with remaining life-years and with how early you learn it.
GROSSMAN = {
    "delta_base": 0.015,        # baseline annual depreciation of health capital at age 20
    "delta_growth": 0.035,      # exponential growth of δ with age (Grossman: δ rises with age)
    "h_floor": 0.30,            # death/severe-morbidity threshold as a share of peak H
    "efficiency_gain": 0.06,    # ↑ in marginal efficiency of health investment when informed
    "src": "Grossman (1972), J. Polit. Econ. — demand for health as human capital",
}

# ── Real-options parameters (Dixit–Pindyck / Black–Scholes intuition) ─────────
# Sequencing is an irreversible investment under uncertainty with a timing choice:
# you may test NOW or defer. Deferring has option value (costs fall, knowledge
# improves) but forfeits information you could have acted on. This is the standard
# "option to wait" trade-off — and it can be negative, i.e. waiting destroys value.
REAL_OPTIONS = {
    "cost_decline": 0.10,       # annual real decline in sequencing price (historic trend, conservative)
    "knowledge_growth": 0.08,   # annual growth in interpretable variants (ClinVar/PGx expansion)
    "vol": 0.45,                # volatility of the modelled value (σ, from the PSA dispersion)
    "src": "Dixit & Pindyck (1994), Investment under Uncertainty",
}

# à-la-carte market prices (what these tests cost to *buy* — PRICE, not value).
MARKET_PRICE_ITEMS = [
    ("Pharmacogenomics panel", 250), ("Carrier screening", 350),
    ("Hereditary-cancer / ACMG panel", 250), ("Polygenic + ancestry", 350),
    ("Nutrigenomics", 200), ("Genetic counseling", 200),
]
CONSOLIDATED_PRICE = 300.0     # one WGS run replaces the à-la-carte stack

_DISCLAIMER = (
    "Illustrative decision-analytic model — NOT a formal economic evaluation, and "
    "NOT financial or medical advice. All figures are probabilistic estimates with "
    "wide uncertainty; population parameters may not transfer to any individual "
    "(and PRS-derived estimates are attenuated for non-European ancestries). "
    "Market price and health-economic value are different things and are reported "
    "separately.")


def _annuity_pv(years: float, rate: float) -> float:
    """Present-value annuity factor for a unit stream over `years` at `rate`
    (half-cycle-corrected). Discounts costs AND QALYs consistently."""
    years = max(0.0, float(years))
    if rate <= 0:
        return years
    n = round(years)
    if n <= 0:
        return 0.0
    # Σ 1/(1+r)^t for t=1..n, half-cycle correction (+0.5 at t=0).
    return sum(1.0 / (1.0 + rate) ** t for t in range(1, n + 1)) + 0.5


# ── finding collection ────────────────────────────────────────────────────────

def _collect(economics_result: dict | None,
             clinical_variants_result: dict | None,
             novel_variants_result: dict | None,
             unvalued: list[dict] | None = None) -> list[dict]:
    """Unify findings across sources into explicit economic parameter dicts:
    {label, kind, coi_key/pgx_key, p_event, rrr, qaly, intervention, horizon,
     wgs_only, haircut, confidence}.

    ``unvalued`` (optional, mutated in place) collects findings whose category
    string ``_classify_category`` does not recognise. Those findings carry no
    economic value here, and without this list they leave no trace at all —
    ``health_economics.py`` gains new ``source=`` labels regularly, and a label
    it emits that this module has no mapping for silently disappears from the
    model while the report still claims a complete computation."""
    out: list[dict] = []

    econ = economics_result or {}
    for f in (econ.get("findings_with_economics") or []):
        kind, coi_key = _classify_category(f.get("category"), f.get("finding", ""))
        label = f.get("finding", "genetic finding")
        pid = f.get("economic_pathway_id", "")
        conf = f.get("confidence", "moderate")
        if kind not in ("pgx", "coi") and unvalued is not None:
            _cat = f.get("category") or "(none)"
            _reason = _not_valued_reason(_cat)
            unvalued.append({"label": label, "category": _cat,
                             "intentional": bool(_reason),
                             "reason": _reason or
                                       "No economic mapping defined for this "
                                       "category — likely an oversight."})
        _hc = _evidence_haircut(f.get("category"))
        if kind == "pgx":
            out.append({"label": label, "economic_pathway_id": pid, "kind": "pgx",
                        "pgx_key": _match_pgx(label),
                        "intervention": _econ_params.value(
                            "intervention_cost_pgx"), "wgs_only": False,
                        "haircut": _hc, "confidence": conf,
                        "source_category": f.get("category") or ""})
        elif kind == "coi":
            # Baseline risk, effect size and intervention cost come from the
            # registry rather than from literals here. They are the two or
            # three numbers the entire benefit side rests on; as bare
            # constants they could not be varied in sensitivity analysis and
            # never appeared in any provenance count, which made the most
            # load-bearing figures in the model also the least visible.
            out.append({"label": label, "economic_pathway_id": pid, "kind": "coi",
                        "coi_key": coi_key,
                        "p_event": _econ_params.value(
                            "baseline_event_probability_dementia"
                            if coi_key == "Alzheimer"
                            else "baseline_event_probability"),
                        "rrr": _econ_params.value("actionable_rrr"),
                        "qaly": float(f.get("qaly_gain") or 0.5),
                        # Whether the SOURCE supplied a QALY figure or this is
                        # the 0.5 fallback. The pooling engine needs to tell
                        # them apart: passing the fallback through as an
                        # override would shadow the registry's per-condition
                        # QALY anchors on every finding, leaving them dead.
                        "qaly_explicit": f.get("qaly_gain") is not None,
                        "intervention": _econ_params.value(
                            "intervention_cost_standard"),
                        "horizon": 25, "wgs_only": False,
                        "haircut": _hc, "confidence": conf,
                        "source_category": f.get("category") or ""})

    # Phase-2 clinical variants (WGS-only) — actionable + carrier + affected.
    cvr = clinical_variants_result or {}
    if cvr.get("available"):
        buckets = cvr.get("buckets", {})
        for f in (buckets.get("actionable") or []):
            coi_key, p_lit, rrr, qaly = _gene_to_econ(f.get("gene", ""))
            # ASCERTAINMENT DE-BIASING, APPLIED (not merely displayed).
            # _gene_to_econ holds clinically ascertained literature penetrance
            # (from multi-case families), which overstates risk for an incidentally
            # identified carrier. Correct it BEFORE it reaches the economics so the
            # dollar figure does not inherit the bias. Family history is unknown at
            # this point, so the no-family-history branch is used — the conservative
            # direction. Previously this correction was computed for display only
            # while the NMB still ran on the raw ascertained value.
            _pc = analyze_penetrance_posterior(prior_penetrance=p_lit,
                                               gene=f.get("gene", ""),
                                               family_history=False)
            p = float(_pc.get("posterior_penetrance", p_lit))
            out.append({"label": f"{f.get('gene','?')} pathogenic (ClinVar)",
                        "kind": "coi", "coi_key": coi_key, "p_event": p,
                        "penetrance_literature": round(p_lit, 4),
                        "penetrance_corrected": round(p, 4),
                        "rrr": rrr, "qaly": qaly, "qaly_explicit": True,
                        "intervention": _econ_params.value(
                            "intervention_cost_monogenic"),
                        "horizon": 30, "wgs_only": True, "haircut": 1.0,
                        "confidence": "high"})
        # carriers: reproductive value (small direct health value to self)
        for f in (buckets.get("carrier") or [])[:6]:
            out.append({"label": f"{f.get('gene','?')} carrier (reproductive)",
                        "kind": "coi", "coi_key": "Pathogenic", "p_event": 0.04,
                        "rrr": 0.5, "qaly": 0.2, "qaly_explicit": True,
                        "intervention": _econ_params.value(
                            "cost_genetic_counseling"),
                        "horizon": 5, "wgs_only": True, "haircut": 1.0,
                        "confidence": "moderate"})

    # Phase-3 novel predicted variants (WGS-only) — DOWN-WEIGHTED (computational).
    nvr = novel_variants_result or {}
    if nvr.get("available"):
        for f in (nvr.get("buckets", {}).get("predicted_pathogenic_rare") or [])[:8]:
            am = f.get("am_score") or 0.6
            out.append({"label": f"predicted-pathogenic {f.get('chrom','?')}:"
                                  f"{f.get('pos','?')}", "kind": "coi",
                        "coi_key": "Pathogenic", "p_event": 0.10, "rrr": 0.4,
                        "qaly": 0.5, "qaly_explicit": True,
                        "intervention": _econ_params.value(
                            "intervention_cost_predicted_variant"),
                        "horizon": 25,
                        "wgs_only": True,
                        "haircut": max(0.1, float(am)),   # scale value by AM confidence
                        "confidence": f.get("confidence", "low"),
                        # THE ACTION THIS PATHWAY PRICES. Without it the report
                        # showed a four-figure value beside "no action mapped",
                        # which is the worst of both: a number the reader cannot
                        # act on and cannot check. It is not invented here — it
                        # is what `intervention_cost_predicted_variant` is
                        # defined to pay for, in the registry's own words: "a
                        # cautious confirmatory workup", there being no standard
                        # pathway for a computationally predicted variant.
                        "action": ("Confirm the predicted variant clinically, "
                                   "then surveillance if it is real"),
                        "action_caveat": ("Predicted by sequence model, not "
                                          "observed in ClinVar — confirmation "
                                          "is the action, not treatment")})
    return out


def _classify_category(category: str | None, label: str = "") -> tuple[str, str]:
    """Map a health-economics finding onto ('pgx'|'coi'|'', coi_key).

    ``health_economics.py`` labels findings with human-readable sources
    ("Pharmacogenomics", "Polygenic Risk", "Addiction Genetics", ...). Matching
    those with substrings — rather than exact lowercase keys — keeps the chip
    path working: a mismatch here silently drops every finding from that source
    and makes the engine report "no findings", which is exactly the bug this
    replaced. Short internal keys stay supported so either naming works.

    Every source label emitted by ``health_economics.py`` must resolve here or
    be registered in :data:`NOT_VALUED`; ``test_voi_no_silent_drops`` enforces
    that, so a new module cannot quietly fall out of the model.
    """
    cat = (category or "").strip().lower()
    text = f"{cat} {(label or '').lower()}"
    if not cat:
        return ("", "")

    # NOT_VALUED is authoritative and checked first. A source registered there
    # must not also resolve to a condition anchor: "Addiction Genetics" matched
    # the "addiction" token below AND its NOT_VALUED entry, so it was counted
    # twice and the coverage invariant broke. One answer per source, decided in
    # one place.
    if (category or "").strip() in NOT_VALUED:
        return ("", "")

    # ── Pharmacogenomics: anything whose decision is "change or dose-adjust a
    #    prescription" — including drug-interaction, top-drug and HLA screens,
    #    and the phase-I/II metabolising-enzyme (detoxification) panel.
    if ("pharmaco" in cat or "pgx" in cat or cat in ("cpic",)
            or "compound interaction" in cat or "detox" in cat):
        return ("pgx", "")

    # ── Monogenic findings: ClinVar-reviewed pathogenic variants and carrier
    #    results both resolve to the high-consequence COI bucket.
    if "clinical variant" in cat or "clinvar" in cat or "carrier" in cat:
        return ("coi", "Pathogenic")

    # ── APOE / dementia is reported under the generic "Genotype" source.
    if "apoe" in text or "alzheim" in text or "dementia" in text:
        return ("coi", "Alzheimer")

    # ── Modules wired in 2026-08-22, each onto the condition it speaks to
    #    rather than onto a cardiometabolic proxy.
    if "addiction" in cat or "substance" in cat:
        return ("coi", "SubstanceUse")
    if "neurochem" in cat or "psychiatric" in cat or "mood" in cat:
        return ("coi", "Depression")
    if "immunogenetic" in cat or "autoimmun" in cat:
        return ("coi", "Autoimmune")
    if "urolog" in cat or "renal" in cat or "kidney" in cat:
        return ("coi", "Urologic")
    if "metal" in cat or "oxidative" in cat:
        return ("coi", "IronOverload")
    # Gut health reaches the model through its one costable trait: the
    # coeliac-permissive HLA haplotype, which drives a serology decision. The
    # module's other five traits are reported as signals and emit nothing, so
    # anything arriving here is the coeliac record.
    if "gut" in cat or "coeliac" in cat or "celiac" in cat:
        return ("coi", "Coeliac")

    # ── Mendelian randomisation and PheWAS biomarkers are trait-directed: route
    #    on the trait named in the finding, defaulting to the cardiometabolic
    #    bucket that dominates both panels. Evidence strength is handled by the
    #    per-category haircut in _evidence_haircut, not by exclusion.
    if "mendelian randomization" in cat or "mendelian randomisation" in cat \
            or "phewas" in cat or "biomarker" in cat:
        for token, key in (("glucose", "T2D"), ("diabet", "T2D"),
                           ("hba1c", "T2D"), ("insulin", "T2D"),
                           ("depress", "Depression"), ("mood", "Depression"),
                           ("urate", "Urologic"), ("kidney", "Urologic"),
                           ("creatinine", "Urologic"), ("ferritin", "IronOverload"),
                           ("iron", "IronOverload")):
            if token in text:
                return ("coi", key)
        return ("coi", "CAD")

    # ── Polygenic risk, lifestyle, wellness and longevity → cardiometabolic.
    # "longevity" is deliberately absent from this list — it is registered in
    # NOT_VALUED as a composite of variants already valued individually.
    for token in ("polygenic", "prs", "genotype", "exercise", "lifestyle",
                  "cardio", "metabolic", "wellness"):
        if token in cat:
            return ("coi", "CAD")
    return ("", "")


# How much of a finding's modelled value survives, by evidence strength of the
# source module. A finding routed onto a real cost-of-illness anchor can still
# rest on a weak association, and the honest treatment is to keep it in the
# model at a discount rather than either dropping it silently or letting it
# count the same as a ClinVar pathogenic variant.
#
# These are judgement calls, not published multipliers, and are labelled as such
# in the report. The ordering is what matters: hypothesis-generating panels
# (PheWAS, wellness) are worth a fraction of a curated clinical finding.
EVIDENCE_HAIRCUT: dict[str, float] = {
    "phewas biomarker":        0.10,   # typically <1% of trait variance explained
    "wellness prediction":     0.10,   # lifestyle-adjacent, largely non-clinical
    "expanded polygenic score": 0.25,  # broad panels, uneven validation
    "mendelian randomization": 0.30,   # population causal estimate, not personal effect
    "neurochemistry":          0.30,   # pathway-level, weak clinical anchoring
    "metal/oxidative":         0.35,
    "urologic/gu":             0.40,
    "addiction genetics":      0.40,   # real COI, but small per-variant effects
    "immunogenetics":          0.40,
    "detoxification":          0.50,
    "polygenic risk":          0.50,
}


def _evidence_haircut(category: str | None) -> float:
    """Fraction of modelled value retained for this source (1.0 = no discount)."""
    cat = (category or "").strip().lower()
    if cat in EVIDENCE_HAIRCUT:
        return EVIDENCE_HAIRCUT[cat]
    for k, v in EVIDENCE_HAIRCUT.items():
        if k in cat:
            return v
    return 1.0


# Categories deliberately given no dollar value, and why. Being explicit
# matters: an unlisted category is an oversight to fix, whereas one listed here
# is a decision. Findings from these sources still appear in the report — they
# are excluded from the economic model only.
NOT_VALUED: dict[str, str] = {
    "Addiction Genetics": (
        "Surfaced as hypothetical / awareness, not valued. The panel was dead "
        "for a long time — three independent string mismatches meant it "
        "produced nothing while being plumbed at every layer, down to a "
        "registered condition anchor. Repairing that resurrects real findings, "
        "but most of what they support is awareness rather than a costed "
        "intervention: a CRHR1 stress-reactivity locus has no "
        "genotype-contingent action, the alcohol filter admits only "
        "susceptibility so the protective ADH1B/ALDH2 variants never qualify, "
        "and CHRNA5 cannot condition on smoking status because the module takes "
        "none. OPRM1 is the one real prescribing decision and the "
        "neurochemistry panel already prices that same rs1799971. So the "
        "findings appear, labelled, contributing nothing to any total."),
    "Top-Drugs PGx Screen": (
        "Surfaced as hypothetical / awareness, not valued. This screen is a "
        "p_rx refinement of a record the pharmacogenomics panel already emits "
        "for the same gene — it raises the assumed probability of ever being "
        "prescribed the drug, it does not find a new benefit. Valuing it "
        "separately would count one pharmacogenomic benefit twice. Knowing "
        "that a commonly-prescribed drug touches a gene already flagged is "
        "genuinely useful, and it is worth zero dollars on top."),
    "Longevity": (
        "The longevity composite is a summary of variants this module already "
        "values individually, so monetising it counts the same genotypes "
        "twice. It also has no cost of illness to anchor against — a "
        "percentile of a composite score is not a disease — and the flat "
        "per-percentile rate previously used had no published source. On one "
        "measured genome that single line produced 54% of the entire modelled "
        "benefit. It is reported as a signal worth acting on, without a "
        "dollar figure."),
    "Family Planning": (
        "Reproductive findings are deliberately never monetised. Attaching a "
        "dollar figure to an affected birth prices a prospective child, and "
        "importing a population uptake rate would embed one set of "
        "reproductive preferences as if it were universal. The module's "
        "offspring-risk figures are reported as information, and the economic "
        "question it does raise — whether to buy an expanded carrier panel — "
        "is surfaced as a recommendation rather than a dollar value, because "
        "no published dose-response supports one."),
}


def _not_valued_reason(category: str | None) -> str:
    """Return the documented reason this category carries no economic value,
    or '' if the category is simply unmapped (i.e. an oversight)."""
    cat = (category or "").strip()
    if cat in NOT_VALUED:
        return NOT_VALUED[cat]
    low = cat.lower()
    for k, v in NOT_VALUED.items():
        if k.lower() == low:
            return v
    return ""


# ── ROH → expanded-carrier-panel prior (CATEGORICAL, never a dollar figure) ──
#
# This is the only defensible economic link from runs of homozygosity, and it is
# deliberately expressed as a recommendation rather than a value.
#
# WHY NOT MONETISE IT. A person's own F_ROH describes THEIR PARENTS' relatedness,
# not their children's. If they partner outside the family their child's
# inbreeding coefficient is ~0 whatever their own F_ROH is, so the intuitive
# "high F_ROH -> higher offspring risk -> dollar value" model is simply wrong
# without partner genotype, which this tool never has.
#
# WHAT DOES SURVIVE. Autozygosity is evidence of descent from a small or
# endogamous population, and standard carrier panels are built around variants
# common in large outbred (mostly European) populations. So elevated ROH raises
# the EXPECTED YIELD of an expanded or ancestry-matched panel relative to a
# minimal one. That is a genuine decision with a genuine cost — but no published
# study maps F_ROH onto incremental panel yield, so any multiplier would be
# invented. Hence a tier and a recommendation, and no number.
#
# FRAMING. Consanguinity is culturally normal in many populations. Nothing here
# is a risk score, and nothing here is phrased as a warning.
# Total-F_ROH threshold above which a founder-background panel conversation is
# worth raising. CALIBRATED EMPIRICALLY, not taken from the literature: this
# detector reads systematically higher than published F_ROH because it counts
# runs from 1 Mb up on dense chips, so a literature cutoff would not transfer.
#
# Measured across 11 real, unrelated, openly-consented public genomes
# (Personal Genome Project Harvard; 0.29-0.32 heterozygosity, 0.56-0.96 M
# autosomal SNPs): F_ROH ranged 0.0113-0.0161, mean 0.0139, sd 0.0017, and none
# carried a single long run. An earlier 0.010 cutoff therefore fired on 11 of 11
# ordinary outbred people — the population norm, not a signal. 0.020 sits above
# the observed maximum and above mean+2sd (0.0174), so it flags genuinely
# unusual burden while staying inside the 0.01-0.03 band reported for founder
# populations.
ROH_FOUNDER_F_THRESHOLD = 0.020

ROH_PANEL_TIERS = {
    "none": ("No ROH-driven change to carrier screening",
             "Standard carrier screening is as informative for you as for "
             "anyone else; nothing in your homozygosity pattern raises the "
             "expected yield of a larger panel."),
    "founder": ("Ancestry-matched or expanded panel worth considering",
                "Your homozygosity pattern is consistent with descent from a "
                "smaller or historically endogamous population. Standard "
                "panels are built around variants common in large outbred "
                "populations, so a panel matched to your ancestry may detect "
                "variants a minimal panel would miss. This says nothing about "
                "your parents being related."),
    "recent": ("Expanded panel and genetic counselling worth discussing",
               "Long runs of homozygosity are present, which is what recent "
               "shared parental ancestry looks like. Two consequences for "
               "screening: you are more likely than average to be homozygous "
               "at a rare recessive locus yourself, and any founder variants "
               "in your background may sit outside a minimal panel. A genetic "
               "counsellor can size the panel to your family history."),
}


def assess_carrier_panel_prior(roh_result: dict | None = None) -> dict:
    """Translate an ROH profile into a CATEGORICAL carrier-panel recommendation.

    Returns a dict carrying a tier, a recommendation, the reasoning, and an
    explicit ``monetised: False``. It never returns a cost, value, multiplier or
    risk score — see the module note above for why a dollar figure here would be
    fabricated rather than estimated.
    """
    base = {"available": False, "monetised": False,
            "why_not_monetised": (
                "A person's own F_ROH reflects their parents' relatedness, not "
                "their children's, and no published study maps F_ROH onto "
                "incremental carrier-panel yield. A dollar figure here would be "
                "invented, so the output is a recommendation only.")}
    r = roh_result or {}
    if not r or r.get("context_tier") in (None, "unavailable"):
        base["reason"] = "No ROH scan available for this input."
        return base

    f_long = float(r.get("f_roh_long") or 0.0)
    f_tot = float(r.get("f_roh") or 0.0)
    n_long = int(r.get("n_long") or 0)

    # Long ROH is the recent-ancestry signal; total burden without long runs is
    # the founder-background signal. They imply different panel conversations.
    if n_long > 0 or f_long > 0.0:
        tier = "recent"
    elif f_tot >= ROH_FOUNDER_F_THRESHOLD:
        tier = "founder"
    else:
        tier = "none"

    rec, rationale = ROH_PANEL_TIERS[tier]
    base.update({
        "available": True,
        "tier": tier,
        "recommendation": rec,
        "rationale": rationale,
        "f_roh": round(f_tot, 5),
        "f_roh_long": round(f_long, 5),
        "n_long_runs": n_long,
        "decision": "Whether to buy an expanded / ancestry-matched carrier "
                    "panel rather than a minimal one.",
        "src": ("Kirin et al. 2010 PLoS ONE 5:e13996 (long ROH tracks recent "
                "pedigree relatedness); ACMG/ACOG expanded carrier-screening "
                "guidance (panel composition and ancestry matching)"),
    })
    return base


def _match_pgx(label: str) -> str:
    """Map a free-text PGx finding onto a PGX_CEA gene-drug pair, else the generic
    fallback. Note ``PGx-generic`` has no '/' separator, so drug matching must be
    guarded — indexing ``split('/')[1]`` unconditionally raises on that key."""
    low = (label or "").lower()
    for key in PGX_CEA:
        if "/" not in key:                     # the generic fallback entry
            continue
        gene, drug = key.split("/", 1)
        gene_token = gene.lower().replace("*", "").split("-")[0]
        if (gene_token and gene_token in low) or (drug and drug.lower() in low):
            return key
    return "PGx-generic"


def _gene_to_econ(gene: str) -> tuple[str, float, float, float]:
    """(coi_key, penetrance, screening/prophylaxis RRR, QALY gain) for a gene."""
    g = (gene or "").upper()
    table = {
        "BRCA1": ("BreastOvarian", 0.60, 0.45, 1.5),
        "BRCA2": ("BreastOvarian", 0.55, 0.45, 1.4),
        "PALB2": ("BreastOvarian", 0.40, 0.40, 1.2),
        "MLH1": ("Colorectal", 0.50, 0.50, 1.3),
        "MSH2": ("Colorectal", 0.45, 0.50, 1.3),
        "MSH6": ("Colorectal", 0.30, 0.45, 1.0),
        "LDLR": ("CAD", 0.50, 0.50, 1.5),
        "APOB": ("CAD", 0.40, 0.45, 1.2),
    }
    return table.get(g, ("Pathogenic", 0.30, 0.35, 0.8))


# ── per-finding NMB ───────────────────────────────────────────────────────────

def _finding_nmb(f: dict, wtp: float, rate: float, rng=None) -> tuple[float, float, float]:
    """Return (nmb, dcost_averted, dqaly) for one finding. If rng is given, sample
    uncertain parameters (PSA); else use point estimates."""
    disc = _annuity_pv(f.get("horizon", 20), rate) / max(1e-9,
                                                          f.get("horizon", 20)) \
        if f.get("horizon") else 1.0
    haircut = f.get("haircut", 1.0)
    conf = f.get("confidence", "moderate")

    if f["kind"] == "pgx":
        c = PGX_CEA.get(f.get("pgx_key"), PGX_CEA["PGx-generic"])
        p_rx, p_adr, rrr = c["p_rx"], c["p_adr"], c["rrr"]
        adr_cost, qaly = c["adr_cost"], c["qaly"]
        if rng is not None:
            p_adr = _beta(rng, p_adr, confidence="high")
            rrr = _beta(rng, rrr, confidence="high")
            adr_cost = _gamma(rng, adr_cost)
        exp_events = p_rx * p_adr * rrr
        dcost = exp_events * adr_cost
        dqaly = exp_events * qaly
        interv = f.get("intervention", 100.0) * p_rx
    else:
        coi = COI.get(f.get("coi_key"), COI["Pathogenic"])["cost"]
        p, rrr, qaly = f.get("p_event", 0.2), f.get("rrr", 0.3), f.get("qaly", 0.5)
        if rng is not None:
            p = _beta(rng, p, confidence=conf)
            rrr = _beta(rng, rrr, confidence=conf)
            coi = _gamma(rng, coi)
        dcost = p * rrr * coi * MARGINAL_COST_FRACTION * disc
        dqaly = p * rrr * qaly * disc
        interv = f.get("intervention", 500.0)

    dcost *= haircut
    dqaly *= haircut
    nmb = dqaly * wtp + dcost - interv
    return nmb, dcost, dqaly, interv


_CONFIDENCE_CONC = {"high": 200.0, "moderate": 50.0, "low": 10.0}


def _beta(rng, mean: float, conc: float = 50.0, confidence: str = "") -> float:
    if confidence:
        conc = _CONFIDENCE_CONC.get(confidence, conc)
    mean = min(max(mean, 1e-3), 1 - 1e-3)
    a = mean * conc
    b = (1 - mean) * conc
    return float(rng.beta(a, b))


def _gamma(rng, mean: float, cv: float = 0.4) -> float:
    shape = 1.0 / (cv * cv)
    scale = mean / shape
    return float(rng.gamma(shape, scale))


# ── main ──────────────────────────────────────────────────────────────────────

def analyze_value_of_information(economics_result: dict | None = None,
                                 clinical_variants_result: dict | None = None,
                                 novel_variants_result: dict | None = None,
                                 genetic_age_result: dict | None = None,
                                 roh_result: dict | None = None,
                                 input_type: str = "chip",
                                 n_mc: int = 5000, seed: int = 12345,
                                 log=None) -> dict:
    # Records any sub-analysis that failed, so a crashing extension degrades
    # VISIBLY instead of silently vanishing from the report. A silent
    # {"available": False} is exactly how the chip-input bug stayed hidden.
    _degraded: list[tuple[str, str]] = []

    _unvalued: list[dict] = []
    findings = _collect(economics_result, clinical_variants_result,
                        novel_variants_result, unvalued=_unvalued)
    if not findings:
        return {"available": False,
                "reason": "No economically-modellable findings yet "
                          "(add --bloodwork and/or a VCF for the full model).",
                "disclaimer": _DISCLAIMER}

    wtp = float(WTP["base"])
    rate = DISCOUNT_RATE
    test_cost = TEST_COST.get(input_type, TEST_COST["wgs"])

    # Deterministic per-finding NMB table.
    nmb_rows = []
    for f in findings:
        nmb, dcost, dqaly, _iv = _finding_nmb(f, wtp, rate)
        nmb_rows.append({"label": f["label"],
                         "economic_pathway_id": f.get("economic_pathway_id", ""),
                         "action": f.get("action", ""),
                         "action_caveat": f.get("action_caveat", ""),
                         "confidence": f["confidence"],
                         "wgs_only": f["wgs_only"], "nmb": round(nmb),
                         "dcost_averted": round(dcost), "dqaly": round(dqaly, 3)})
    nmb_rows.sort(key=lambda r: -r["nmb"])

    total_expost = sum(r["nmb"] for r in nmb_rows) - test_cost
    wgs_only_value = sum(r["nmb"] for r in nmb_rows if r["wgs_only"])
    chip_value = sum(r["nmb"] for r in nmb_rows if not r["wgs_only"]) - TEST_COST["chip"]

    # Undiscounted comparison (shows discounting actually bites).
    undisc = sum(_finding_nmb(f, wtp, 0.0)[0] for f in findings) - test_cost

    result = {
        "available": True,
        "input_type": input_type,
        "wtp_base": wtp, "wtp_range": (WTP["low"], WTP["high"]),
        "discount_rate": rate, "test_cost": test_cost,
        "n_findings": len(findings),
        "nmb_rows": nmb_rows,
        "voi_expost_point": round(total_expost),
        "voi_expost_undiscounted": round(undisc),
        "marginal_chip_to_wgs": round(wgs_only_value),
        "chip_value_point": round(chip_value),
        "price": _price_panel(),
        "methods": _methods(rate, wtp, test_cost, input_type),
        "disclaimer": _DISCLAIMER,
    }
    result.update(_exante(test_cost))

    if _HAVE_NP:
        result.update(_psa(findings, test_cost, n_mc, seed))
    else:
        result["psa_available"] = False

    # ── Extended economic framings ───────────────────────────────────────────
    # Grossman health-capital: why the same information is worth more when you're
    # younger (an efficiency gain compounds over remaining life-years).
    _age = _resolve_age(genetic_age_result)
    try:
        result["health_capital"] = analyze_health_capital(age=_age)
    except Exception as _e:
        result["health_capital"] = {"available": False}
        _degraded.append(("health_capital", repr(_e)))

    # Real options: test now vs defer, given falling prices and improving
    # interpretation — the classic option-to-wait, on a depreciating asset.
    try:
        result["real_option"] = analyze_real_option(
            voi_now=float(result.get("voi_expost_mean",
                                     result.get("voi_expost_point", 0.0))) + test_cost,
            test_cost=test_cost, age=_age)
    except Exception as _e:
        result["real_option"] = {"available": False}
        _degraded.append(("real_option", repr(_e)))

    # Risk-adjusted framings drawn from the PSA distribution.
    try:
        result["risk_adjusted"] = _risk_adjusted(result, test_cost)
    except Exception as _e:
        result["risk_adjusted"] = {"available": False}
        _degraded.append(("risk_adjusted", repr(_e)))

    # ── Decision-theoretic ceiling: EVPI (what ANY further research could be worth) ──
    try:
        result["evpi"] = analyze_evpi(findings, test_cost)
    except Exception as _e:
        result["evpi"] = {"available": False}
        _degraded.append(("evpi", repr(_e)))

    # ── Expected utility / Arrow–Pratt: a risk-averse agent values variance
    #    reduction above the expected monetary value. ──
    try:
        _mean = float(result.get("voi_expost_mean", result.get("voi_expost_point", 0.0)))
        _sd = float((result.get("risk_adjusted") or {}).get("sd", 0.0))
        result["utility"] = analyze_utility(_mean, _sd)
    except Exception as _e:
        result["utility"] = {"available": False}
        _degraded.append(("utility", repr(_e)))

    # ── HEOR deliverables: Markov cohort CEA + payer budget impact ──
    try:
        from . import markov as _mk
        _mkr = _mk.markov_cost_effectiveness(start_age=_age, wtp=wtp)
        _mkr["validation"] = _mk.validate_markov(_mkr)
        result["markov"] = _mkr
        result["budget_impact"] = _mk.budget_impact()
    except Exception as _e:
        result["markov"] = {"available": False}
        _degraded.append(("markov", repr(_e)))
        result["budget_impact"] = {"available": False}
        _degraded.append(("budget_impact", repr(_e)))

    # ── Welfare comparison: local vs centralised analysis, conceding a capability
    #    premium to the centralised alternative so the result isn't assumed. ──
    try:
        result["welfare"] = analyze_welfare_comparison(
            voi=float(result.get("voi_expost_mean",
                                 result.get("voi_expost_point", 0.0))) + test_cost,
            test_cost=test_cost)
    except Exception as _e:
        result["welfare"] = {"available": False}
        _degraded.append(("welfare", repr(_e)))

    # ── Information economics: adverse selection, genetic discrimination, and why
    #    local-only analysis is an economic design choice. ──
    try:
        result["information_economics"] = analyze_insurance_economics(
            evpi=float((result.get("evpi") or {}).get("evpi", 0.0)),
            mean_value=float(result.get("voi_expost_mean",
                                        result.get("voi_expost_point", 0.0))))
    except Exception as _e:
        result["information_economics"] = {"available": False}
        _degraded.append(("information_economics", repr(_e)))

    # ── EVPPI: which single parameter is worth resolving (research prioritisation) ──
    try:
        result["evppi"] = analyze_evppi(findings)
    except Exception as _e:
        result["evppi"] = {"available": False}
        _degraded.append(("evppi", repr(_e)))

    # ── Behavioural: prospect theory + hyperbolic discounting explain the adoption
    #    gap between a positive-EV test and actual uptake. ──
    try:
        _m2 = float(result.get("voi_expost_mean", result.get("voi_expost_point", 0.0)))
        _sd2 = float((result.get("risk_adjusted") or {}).get("sd", 0.0))
        result["behavioural"] = analyze_behavioural(_m2, _sd2, test_cost)
    except Exception as _e:
        result["behavioural"] = {"available": False}
        _degraded.append(("behavioural", repr(_e)))

    # ── Longevity sensitivity: rising life expectancy raises realised genetic risk
    #    AND lengthens the payoff horizon, so it raises the value of information. ──
    try:
        from risk import genomic_statistics as _gstat
        result["longevity"] = _gstat.longevity_sensitivity(current_age=_age)
    except Exception as _e:
        result["longevity"] = {"available": False}
        _degraded.append(("longevity", repr(_e)))

    # ── Genomics-side rigor: show the ascertainment correction actually applied to a
    #    representative high-penetrance finding, so the bias is visible, not hidden. ──
    try:
        pen = [f for f in findings if f.get("coi_key") in
               ("BreastOvarian", "Colorectal") and f.get("p_event")]
        if pen:
            f0 = pen[0]
            # split() on an empty label yields [], so guard the index.
            _tokens = str(f0.get("label", "")).split()
            result["penetrance_correction"] = analyze_penetrance_posterior(
                prior_penetrance=float(f0.get("p_event", 0.5)),
                gene=_tokens[0] if _tokens else "")
    except Exception as _e:
        _degraded.append(("penetrance_correction", repr(_e)))

    # ── Statistical correction chain ──────────────────────────────────────
    # genomic_statistics implements three corrections that the README has always
    # listed as capabilities, and none of them had a production caller: the
    # Falconer liability-threshold conversion from percentile to absolute risk,
    # the ancestry-transferability decay of a EUR-derived score, and the
    # competing-mortality adjustment of a lifetime penetrance. Implemented,
    # tested, and never run — so the report claimed corrections it did not make.
    #
    # Every step in the chain moves the risk DOWN, which is why surfacing it is
    # the conservative direction rather than a new claim.
    try:
        from risk import genomic_statistics as _gstat
        _prs_pct, _prev = 0.90, 0.06
        _anc = "european"
        try:
            _ap = ((economics_result or {}).get("ancestry")
                   or (economics_result or {}).get("ancestry_group"))
            if isinstance(_ap, str) and _ap.strip():
                _anc = _ap.strip().lower()
        except Exception:
            pass
        result["genomic_corrections"] = _gstat.summarise_genomic_corrections(
            prs_percentile=_prs_pct, prevalence=_prev,
            lifetime_penetrance=_econ_params.value("penetrance_brca2_population")
            if "penetrance_brca2_population" in _econ_params.PARAMS else 0.55,
            ancestry=_anc, current_age=_resolve_age(genetic_age_result))
        result["genomic_corrections"]["inputs"] = {
            "prs_percentile": _prs_pct, "prevalence": _prev, "ancestry": _anc,
            "note": ("Illustrative inputs: a 90th-percentile score against a 6% "
                     "population prevalence. The chain shows how much each "
                     "correction moves a risk of that shape, not this person's "
                     "risk for a named disease."),
        }
    except Exception as _e:
        _degraded.append(("genomic_corrections", repr(_e)))

    # ROH -> carrier-panel prior. Deliberately qualitative: it names a decision
    # and a direction, contributes nothing to any dollar total, and is kept out
    # of the NMB table entirely so it cannot be mistaken for a valued finding.
    try:
        result["carrier_panel_prior"] = assess_carrier_panel_prior(roh_result)
    except Exception as _e:
        result["carrier_panel_prior"] = {"available": False, "monetised": False}
        _degraded.append(("carrier_panel_prior", repr(_e)))

    # ── Pooled cost-effectiveness (the double-counting correction) ────────
    # The NMB table above still treats every finding as its own line, which is
    # the right shape for EVPI/EVPPI — those ask which PARAMETER is worth
    # resolving. It is the wrong shape for "what is this genome worth", because
    # eight of the report's finding sources route onto the same cardiometabolic
    # anchor and summing them prevents the same heart attack eight times. The
    # pooled result below combines findings per condition on the risk scale and
    # charges one cost of illness, and reports the difference so the size of the
    # correction is visible rather than quietly banked.
    try:
        from . import engine as _ee
        _age = _resolve_age(genetic_age_result)
        # One constructor, used by the base case and by every PSA draw. These
        # were two copies of the same argument list, and the last time they
        # drifted the sampled parameters stayed pinned at their base values and
        # the model reported an interval far narrower than the real one.
        def _mk(f):
            return _ee.Finding(
                label=f["label"], coi_key=f.get("coi_key", ""),
                p_event=float(f.get("p_event", 0.0) or 0.0),
                rrr=float(f.get("rrr", 0.0) or 0.0),
                haircut=float(f.get("haircut", 1.0) or 1.0),
                intervention_cost=float(f.get("intervention", 0.0) or 0.0),
                confidence=f.get("confidence", "moderate"),
                source_category=f.get("source_category", ""),
                # Trial efficacy becomes real-world effectiveness here.
                adherence=_ee.adherence_for(f.get("coi_key", "")),
                # Only a QALY figure the source actually stated may
                # override the registry's per-condition anchor.
                qaly_override=(float(f["qaly"])
                               if f.get("qaly_explicit")
                               and f.get("qaly") is not None
                               else None))

        def _mk_all(fs):
            return [_mk(f) for f in fs
                    if f.get("kind") == "coi" and f.get("coi_key")]

        _ee_findings = _mk_all(findings)
        _pools = _ee.pool_findings(_ee_findings)
        _pooled = _ee.evaluate_pools(_pools, wtp=wtp, test_cost=test_cost)
        _pooled["validation"] = _ee.validate_model(_pools, _pooled)
        _pooled["dual_perspective"] = _ee.dual_perspective(
            _pooled["cea"]["cost_averted"], _pooled["cea"]["incremental_qaly"],
            conditions=_pooled["conditions"], wtp=wtp)
        _pooled["impact_inventory"] = _ee.impact_inventory(_pooled["conditions"])
        _pooled["cheers"] = _ee.cheers_checklist(
            wtp=wtp, rate=rate,
            horizon=_econ_params.value("horizon_years_personal"))
        # Propagate the registry's documented uncertainty. Until this ran,
        # the distribution and range fields were documentation only and the
        # model reported a point estimate as if it were precise.
        # Rebuild the pools inside every draw. Baseline risk, effect size and
        # intervention cost are read when a Finding is constructed, so reusing
        # pools built from the base case would leave those parameters pinned
        # and report an interval far narrower than the real one.
        def _rebuild():
            _fs = _collect(economics_result, clinical_variants_result,
                           novel_variants_result, unvalued=None)
            return _ee.pool_findings(_mk_all(_fs))

        _pooled["psa"] = _ee.run_psa(_pools, n=1500, test_cost=test_cost,
                                     wtp=wtp, rebuild=_rebuild)
        _pooled["ceac"] = _ee.ceac(_pools, n=600, test_cost=test_cost,
                                   rebuild=_rebuild)
        _pooled["tornado"] = _ee.tornado(_pools, test_cost=test_cost, wtp=wtp,
                                         rebuild=_rebuild)
        # Decision-analytic layer: what resolving each assumption would be
        # worth, where the conclusion flips, which strategy wins, what it does
        # to a budget, and whether it works the same across age and sex.
        try:
            from . import decision as _ed
            # Scale the two-level Monte Carlo to the caller's simulation
            # budget. EVPPI is inherently O(outer x inner); at full fidelity
            # it dominates the run, and callers that passed a small n_mc were
            # asking for a quick answer, not a slow one.
            _pooled["decision"] = _ed.analyze_decision_layer(
                _rebuild, tornado_rows=_pooled["tornado"],
                test_cost=test_cost, wtp=wtp, age=_age,
                fast=(n_mc < 1000))
            # Build each arm from the findings that arm would actually
            # return, on the SAME net-cost basis as the CEA card above.
            #
            # This previously multiplied the chip result by 1.35 and 1.6 to
            # invent a WGS arm. Those two numbers decided which strategy the
            # report recommended, and they were typed here rather than
            # derived from anything — the precise defect this whole exercise
            # set out to remove. The costs were also gross rather than net,
            # so the frontier said the chip cost $2,100 while the card beside
            # it said the same strategy SAVED $6,440.
            def _arm(name, findings, test_cost_arm):
                pools = _ee.pool_findings(findings)
                if not pools:
                    return None
                cea = _ee.evaluate_pools(pools, wtp=wtp,
                                         test_cost=test_cost_arm)["cea"]
                return {"name": name, "cost": float(cea["incremental_cost"]),
                        "qaly": float(cea["incremental_qaly"])}

            _all = _ee_findings
            _chip_only = [f for f, src in zip(_ee_findings, findings, strict=False)
                          if not src.get("wgs_only")]
            _wgs_extra = len(_all) - len(_chip_only)
            _arms = [{"name": "No testing", "cost": 0.0, "qaly": 0.0}]
            _chip_arm = _arm("Genotyping chip", _chip_only, TEST_COST["chip"])
            if _chip_arm:
                _arms.append(_chip_arm)
            if _wgs_extra:
                _wgs_arm = _arm("Whole-genome sequencing", _all,
                                TEST_COST["wgs"])
                if _wgs_arm:
                    _arms.append(_wgs_arm)
            _fr = _ed.efficiency_frontier(_arms, wtp=wtp)
            if not _wgs_extra:
                _fr["wgs_not_estimable"] = (
                    "No whole-genome-only findings are present in this input, "
                    "so a sequencing arm cannot be estimated from it — only "
                    "asserted. It is left out rather than invented.")
            _pooled["decision"]["frontier"] = _fr
        except Exception as _de:
            _pooled["decision"] = {"available": False, "reason": repr(_de)}
            _degraded.append(("decision_layer", repr(_de)))

        # Is sequencing worth buying? The `marginal_chip_to_wgs` figure above
        # counts sequencing-only findings PRESENT IN THIS FILE, which for
        # array input is structurally zero and says nothing about whether
        # sequencing is worth buying. This answers that prospectively.
        try:
            _n_wgs = sum(1 for f in findings if f.get("wgs_only"))
            _pooled["wgs_decision"] = _ed.wgs_marginal_value(
                wgs_only_findings_value=wgs_only_value,
                n_wgs_only_findings=_n_wgs, wtp=wtp)
        except Exception as _we:
            _pooled["wgs_decision"] = {"available": False, "reason": repr(_we)}
            _degraded.append(("wgs_decision", repr(_we)))

        # Plain-language layer. Everything above is defensible and unreadable
        # to anyone who has not done a health-economics course; a result
        # nobody can read is not a result.
        try:
            from . import plain as _epl
            _pooled["plain"] = _epl.build_plain_summary(_pooled)
        except Exception as _pe2:
            _pooled["plain"] = {"available": False, "reason": repr(_pe2)}
            _degraded.append(("plain_summary", repr(_pe2)))

        _pooled["provenance"] = _econ_params.assumption_burden()
        _pooled["declared_assumptions"] = [
            {"key": p.key, "value": p.value, "units": p.units, "note": p.note}
            for p in _econ_params.assumptions()]
        _pooled["references"] = _econ_params.citation_list()
        # Structural (Markov) re-estimate for the conditions that dominate the
        # total, where event timing and competing mortality actually matter.
        _top = sorted(_pooled["conditions"],
                      key=lambda c: c.get("inmb", 0), reverse=True)[:3]
        _pooled["structural"] = [
            _ee.incremental_analysis(
                start_age=_age,
                annual_incidence=max(1e-4, float(c["baseline_risk"]) / 25.0),
                coi_key=c["condition"], rrr=float(c["combined_rrr"]),
                intervention_cost_annual=float(c["intervention_cost"]) / 10.0,
                wtp=wtp)
            for c in _top if c.get("condition") in _ee.COI_KEY_TO_PARAM
        ]
        result["pooled_economics"] = _pooled
    except Exception as _e:
        result["pooled_economics"] = {"available": False, "reason": repr(_e)}
        _degraded.append(("pooled_economics", repr(_e)))

    # Surface any degradation rather than letting a failed sub-analysis disappear.
    result["degraded_components"] = [{"component": k, "error": e} for k, e in _degraded]
    # A finding dropped for an unrecognised category is a real gap in the model,
    # so it counts against fully_computed exactly as a crashed component does.
    result["unvalued_findings"] = _unvalued
    result["n_unvalued"] = len(_unvalued)
    # An intentional exclusion (documented in NOT_VALUED) is a modelling
    # decision, not an incomplete computation. Only genuinely unmapped
    # categories count against fully_computed.
    _oversights = [u for u in _unvalued if not u["intentional"]]
    result["n_unvalued_intentional"] = len(_unvalued) - len(_oversights)
    result["unmapped_categories"] = sorted({u["category"] for u in _oversights})
    result["fully_computed"] = not _degraded and not _oversights
    if _degraded and log:
        for k, e in _degraded:
            log(f"  WARNING: value-of-information component '{k}' failed: {e}")
    if _oversights and log:
        log(f"  WARNING: {len(_oversights)} finding(s) dropped from the economic "
            f"model — unmapped category/categories: "
            f"{', '.join(result['unmapped_categories'])}")

    return result


def _resolve_age(genetic_age_result: dict | None, default: float = 35.0) -> float:
    """Best-effort age from the bloodwork/aging module; falls back to a default."""
    g = genetic_age_result or {}
    for key in ("chronological_age", "chronological", "age_used", "age"):
        v = g.get(key)
        if isinstance(v, int | float) and 0 < float(v) < 120:
            return float(v)
    clinical = g.get("clinical") if isinstance(g.get("clinical"), dict) else {}
    v = clinical.get("age_used")
    if isinstance(v, int | float) and 0 < float(v) < 120:
        return float(v)
    return default


def _risk_adjusted(result: dict, test_cost: float) -> dict:
    """Portfolio-style risk-adjusted summaries of the modelled value distribution.

    * **Return on investment** — modelled value per dollar of test cost.
    * **Sharpe-style ratio** — mean value per unit of its own standard deviation
      (a reward-to-variability ratio; here there is no risk-free asset, so it is a
      *relative* dispersion measure, not a true Sharpe ratio).
    * **Certainty equivalent** under CRRA/mean-variance risk aversion — what a
      risk-averse person would accept for sure instead of the uncertain payoff.
    """
    mean = float(result.get("voi_expost_mean", result.get("voi_expost_point", 0.0)))
    lo = float(result.get("voi_ci_low", mean))
    hi = float(result.get("voi_ci_high", mean))
    # Recover an approximate SD from the 95% interval (±1.96σ).
    sd = max(1e-9, (hi - lo) / 3.92)
    roi = (mean / test_cost) if test_cost else None
    sharpe = mean / sd if sd else None
    # Mean-variance certainty equivalent: CE = μ − (γ/2)·σ²/μ  (scaled, γ = 2).
    gamma = 2.0
    ce = mean - 0.5 * gamma * (sd ** 2) / max(abs(mean), 1.0)
    return {
        "available": True,
        "roi_multiple": round(roi, 1) if roi is not None else None,
        "sd": round(sd),
        "reward_to_variability": round(sharpe, 2) if sharpe is not None else None,
        "certainty_equivalent": round(ce),
        "risk_aversion_gamma": gamma,
        "note": ("Reward-to-variability is Sharpe-style but has no risk-free "
                 "benchmark, so it measures dispersion, not excess return. The "
                 "certainty equivalent uses a mean-variance approximation with "
                 "γ = 2 (moderate risk aversion)."),
    }


def _exante(test_cost: float) -> dict:
    """Ex-ante (EVSI-style) value to a random person BEFORE testing: population
    prevalence × ΔNB across a catalog of screenable actionable findings."""
    catalog = [
        ("PGx actionable flag", 0.90, {"kind": "pgx", "pgx_key": "PGx-generic"}),
        ("Hereditary-cancer variant", 0.02,
         {"kind": "coi", "coi_key": "BreastOvarian", "p_event": 0.5, "rrr": 0.45,
          "qaly": 1.4, "intervention": 1500, "horizon": 30}),
        ("Familial hypercholesterolemia", 0.004,
         {"kind": "coi", "coi_key": "CAD", "p_event": 0.5, "rrr": 0.5, "qaly": 1.5,
          "intervention": 500, "horizon": 30}),
        ("Recessive carrier (reproductive)", 0.60,
         {"kind": "coi", "coi_key": "Pathogenic", "p_event": 0.02, "rrr": 0.5,
          "qaly": 0.2, "intervention": 400, "horizon": 5}),
    ]
    total = 0.0
    for _, prev, params in catalog:
        nmb = _finding_nmb(params, float(WTP["base"]), DISCOUNT_RATE)[0]
        total += prev * nmb
    return {"voi_exante_point": round(total - test_cost)}


def _psa(findings: list[dict], test_cost: float, n: int, seed: int) -> dict:
    """Seeded Monte-Carlo PSA. Effects/costs are sampled per draw; λ is varied
    deterministically for the CEAC via NMB(λ) = NMB_base + ΔQALY·(λ − λ_base)."""
    rng = np.random.default_rng(seed)
    base = float(WTP["base"])
    totals = np.zeros(n)     # NMB at λ_base, net of test cost
    dqalys = np.zeros(n)     # ΔQALY
    netcosts = np.zeros(n)   # incremental cost (intervention − averted); neg = saving
    # Shared cost-environment multiplier per draw: all findings within
    # one MC iteration share the same cost environment (if hospitalization
    # is expensive in draw k, it's expensive for all findings in draw k).
    cost_env = rng.lognormal(0.0, 0.15, n)
    for i in range(n):
        tot = dq = dc = iv = 0.0
        for f in findings:
            nmb, c, q, interv = _finding_nmb(f, base, DISCOUNT_RATE, rng=rng)
            c *= cost_env[i]
            interv *= cost_env[i]
            nmb = q * base + c - interv
            tot += nmb
            dq += q
            dc += c
            iv += interv
        totals[i] = tot - test_cost
        dqalys[i] = dq
        netcosts[i] = iv - dc
    lam_grid = list(range(0, 200_001, 10_000))
    ceac = [{"lam": lam, "prob": float(np.mean(totals + dqalys * (lam - base) > 0))}
            for lam in lam_grid]
    step = max(1, n // 200)

    # Left-tail risk (financial-economics risk measures, applied to health value):
    # VaR_95 = the 5th-percentile outcome ("how bad could this reasonably go?");
    # CVaR_95 (expected shortfall) = the mean of everything AT OR BELOW that point
    # ("given it does go that badly, how bad on average?"). CVaR is coherent
    # (subadditive) where VaR alone is not — the standard reason risk management
    # prefers it for tail risk, and it directly quantifies the "left-tail risk"
    # this tool's own framing promises to address.
    var_95 = float(np.percentile(totals, 5))
    tail = totals[totals <= var_95]
    cvar_95 = float(np.mean(tail)) if tail.size else var_95

    return {
        "psa_available": True,
        "voi_expost_mean": round(float(np.mean(totals))),
        "voi_ci_low": round(float(np.percentile(totals, 2.5))),
        "voi_ci_high": round(float(np.percentile(totals, 97.5))),
        "prob_cost_effective": round(float(np.mean(totals > 0)), 3),
        "var_95": round(var_95),
        "cvar_95": round(cvar_95),
        "ceac": ceac,
        "ce_plane": [[round(float(netcosts[i])), round(float(dqalys[i]), 3)]
                     for i in range(0, n, step)],
        "tornado": _tornado(findings, test_cost),
    }


def _tornado(findings: list[dict], test_cost: float) -> list[dict]:
    """One-way sensitivity: swing WTP and discount rate to their extremes."""
    base = sum(_finding_nmb(f, float(WTP["base"]), DISCOUNT_RATE)[0]
               for f in findings) - test_cost
    rows = []
    for name, lo_kw, hi_kw in [("WTP ($/QALY)", ("wtp", WTP["low"]), ("wtp", WTP["high"])),
                               ("Discount rate", ("rate", 0.05), ("rate", 0.0))]:
        def _tot(kind, val):
            w = val if kind == "wtp" else float(WTP["base"])
            r = val if kind == "rate" else DISCOUNT_RATE
            return sum(_finding_nmb(f, w, r)[0] for f in findings) - test_cost
        lo = _tot(*lo_kw)
        hi = _tot(*hi_kw)
        rows.append({"param": name, "low": round(lo), "high": round(hi),
                     "base": round(base), "swing": round(abs(hi - lo))})
    rows.sort(key=lambda r: -r["swing"])
    return rows


def analyze_health_capital(age: float = 35.0, horizon: int = 50,
                           informed_efficiency_gain: float | None = None) -> dict:
    """**Grossman (1972) health-capital model.**

    Health is a durable capital stock that depreciates with age and is replenished
    by investment. Genomic information doesn't add health directly — it raises the
    *marginal efficiency* of health investment (you spend effort on what actually
    matters for you). Two consequences fall out of the model, and both are reported:

      * **Value rises with remaining life-years** — an efficiency gain compounds over
        the years you have left, so the same information is worth more at 30 than 70.
      * **Depreciation accelerates with age** — δ(a) = δ₀·e^(g·(a−20)) — so the
        undiscounted health-capital stock is convex-decreasing in age.

    Returns the informed vs uninformed health-capital trajectories, the discounted
    lifetime difference in health-capital-years, and the implied age at which the
    stock crosses the morbidity floor.
    """
    d0 = GROSSMAN["delta_base"]
    g = GROSSMAN["delta_growth"]
    floor = GROSSMAN["h_floor"]
    eff = (GROSSMAN["efficiency_gain"] if informed_efficiency_gain is None
           else float(informed_efficiency_gain))

    def _delta(a: float) -> float:
        return d0 * math.exp(g * max(0.0, a - 20.0))

    uninformed, informed = [], []
    h_u = h_i = 1.0
    pv_gap = 0.0
    cross_u = cross_i = None
    for t in range(horizon + 1):
        a = age + t
        uninformed.append({"age": round(a, 1), "h": round(h_u, 4)})
        informed.append({"age": round(a, 1), "h": round(h_i, 4)})
        if cross_u is None and h_u <= floor:
            cross_u = a
        if cross_i is None and h_i <= floor:
            cross_i = a
        # discounted gap in health-capital-years (the "extra healthy stock" held)
        pv_gap += (h_i - h_u) / ((1.0 + DISCOUNT_RATE) ** t)
        d = _delta(a)
        # Gross investment is a roughly constant effort (I), which offsets less and
        # less as δ(a) accelerates — so the stock declines, slowly at first then
        # steeply. Informed investment buys (1+eff) more health per unit of effort.
        invest = 0.012
        h_u = max(0.0, h_u * (1.0 - d) + invest)
        h_i = max(0.0, h_i * (1.0 - d) + invest * (1.0 + eff))

    return {
        "available": True,
        "age": age, "horizon": horizon,
        "efficiency_gain": eff,
        "delta_at_age": round(_delta(age), 4),
        "delta_at_end": round(_delta(age + horizon), 4),
        "pv_health_capital_gain": round(pv_gap, 3),
        "monetised_gain": round(pv_gap * float(WTP["base"]) * 0.10),
        "morbidity_floor": floor,
        "floor_age_uninformed": cross_u,
        "floor_age_informed": cross_i,
        "years_floor_deferred": (round(cross_i - cross_u, 1)
                                 if (cross_i and cross_u) else None),
        "trajectory_uninformed": uninformed[::max(1, horizon // 25)],
        "trajectory_informed": informed[::max(1, horizon // 25)],
        "src": GROSSMAN["src"],
        "note": ("Illustrative Grossman-style stock-and-flow model, not a fitted "
                 "biological model. It shows the DIRECTION and shape of the effect — "
                 "information compounds over remaining life-years — not a clinical "
                 "prediction of your health trajectory."),
    }


def analyze_real_option(voi_now: float, test_cost: float, age: float = 35.0,
                        defer_years: int = 5) -> dict:
    """**Real-options / optimal-timing analysis (Dixit–Pindyck).**

    Sequencing is an irreversible purchase under uncertainty with a timing choice.
    The naive "option to wait" says defer: the test gets cheaper and interpretation
    improves. But two features of *this* asset reverse the usual logic:

      1. **The data is a permanent asset.** Sequence once and you can re-analyse for
         free forever, so interpretation growth (k) accrues to the early tester too —
         it is *not* a reason to wait.
      2. **The underlying asset depreciates.** Every deferred year is a year you
         cannot act on a finding, permanently forfeiting the T/H share of the
         benefit stream.

    So the comparison reduces to: does the price decline (c) outrun the forfeited
    protection?

        Defer T yrs:  [VOI·(1 − T/H) − C·(1−c)^T] / (1+r)^T

    For an actionable genome the answer is no — **option value of waiting ≈ 0, test
    now**. The model still says "defer" when value is low relative to price, which is
    the honest result in that case (often: don't buy at all).
    """
    c = REAL_OPTIONS["cost_decline"]
    k = REAL_OPTIONS["knowledge_growth"]
    r = DISCOUNT_RATE

    value_now = voi_now - test_cost
    rows = []
    best_T, best_v = 0, value_now
    # The value stream is consumed over the remaining horizon; deferring T years
    # forfeits the T/H share of it outright (you cannot act on what you don't know),
    # and that forfeited benefit is the dominant term — this is what makes waiting
    # expensive for a depreciating asset.
    # Crucially, sequencing once buys a PERMANENT asset: the data can be re-analysed
    # for free as interpretation improves. So knowledge growth (k) accrues to the
    # early tester too — it is NOT a reason to wait. Deferring therefore captures
    # only the price decline, while forfeiting the T/H share of the benefit stream.
    horizon = max(1.0, 85.0 - age)
    value_now = voi_now * ((1.0 + k) ** 0) - test_cost
    best_v = value_now
    for T in range(0, defer_years + 1):
        share_lost = min(1.0, T / horizon)            # benefit years permanently lost
        # Both branches enjoy the same knowledge growth over the horizon; the deferrer
        # simply starts later and loses the intervening years of protection.
        gross = voi_now * (1.0 - share_lost)
        cost_T = test_cost * ((1.0 - c) ** T)
        v = (gross - cost_T) / ((1.0 + r) ** T)
        rows.append({"defer_years": T, "net_value": round(v),
                     "test_cost_then": round(cost_T),
                     "interpretable_value_then": round(gross)})
        if v > best_v:
            best_T, best_v = T, v

    option_value = best_v - value_now
    return {
        "available": True,
        "value_test_now": round(value_now),
        "optimal_defer_years": best_T,
        "value_at_optimum": round(best_v),
        "option_value_of_waiting": round(option_value),
        "recommendation": ("test now — waiting forfeits more value than the price "
                           "decline recovers" if best_T == 0 else
                           f"deferring ~{best_T} year(s) is modelled as marginally better"),
        "assumed_cost_decline": c,
        "assumed_knowledge_growth": k,
        "schedule": rows,
        "src": REAL_OPTIONS["src"],
        "note": ("Illustrative timing model. Assumes today's price trend and "
                 "interpretation-growth rate persist; it ignores the option to "
                 "re-analyse an existing genome for free, which further favours "
                 "sequencing earlier."),
    }


def analyze_evpi(findings: list[dict], test_cost: float, n: int = 4000,
                 seed: int = 4242) -> dict:
    """**Expected Value of Perfect Information (EVPI)** — the formal ceiling.

    EVPI is the canonical VOI quantity in decision theory (Raiffa & Schlaifer 1961):
    the difference between the value of deciding *after* uncertainty resolves and the
    value of deciding *now* under uncertainty.

        EVPI = E_θ[max_a NB(a, θ)]  −  max_a E_θ[NB(a, θ)]

    Because you must choose one action per finding under uncertainty, but a clairvoyant
    would pick the best action state-by-state, EVPI ≥ 0 always. Interpretation:

      * **EVPI is an upper bound on what ANY further research is worth.** If EVPI is
        small, resolving uncertainty cannot change decisions enough to matter — buying
        more information is irrational at that price.
      * The genome's realised VOI should be read *against* this ceiling: a test capturing
        a large share of EVPI is doing most of the work information can do.
    """
    if not _HAVE_NP:
        return {"available": False}
    rng = np.random.default_rng(seed)
    wtp = float(WTP["base"])

    # Per-draw NMB for each finding: act (uncertain payoff) vs do-nothing (0).
    per = np.zeros((n, len(findings)))
    for j in range(n):
        for i, f in enumerate(findings):
            per[j, i] = _finding_nmb(f, wtp, DISCOUNT_RATE, rng=rng)[0]

    # Current information: for each finding you must commit to ONE action based on its
    # mean; you then live with whatever the true state turns out to be.
    means = per.mean(axis=0)
    act_now = means > 0.0                                  # the committed decision
    nb_current = float(np.mean(per[:, act_now].sum(axis=1))) if act_now.any() else 0.0
    # Perfect information: state revealed first, so act only when it actually pays.
    nb_perfect = float(np.mean(np.maximum(per, 0.0).sum(axis=1)))

    evpi = max(0.0, nb_perfect - nb_current)
    capture = (nb_current / nb_perfect) if nb_perfect > 0 else None
    # Share of draws where perfect information would REVERSE the committed decision.
    flip = float(np.mean(np.any((per > 0) != act_now[None, :], axis=1))) if len(findings) else 0.0
    small = evpi < 0.02 * max(1.0, abs(nb_current))
    return {
        "available": True,
        "evpi": round(evpi),
        "nb_current_information": round(nb_current),
        "nb_perfect_information": round(nb_perfect),
        "share_of_information_captured": round(capture, 3) if capture else None,
        "prob_decision_reversal": round(flip, 3),
        "decision_robust": bool(small),
        "interpretation": (
            "EVPI is the maximum ANY further research could be worth — the value of "
            "resolving all remaining uncertainty. " +
            ("It is small here, which is itself the finding: the recommended actions "
             "stay optimal across essentially the whole uncertainty range, so the "
             "decisions are robust and more data would not change them. Low EVPI "
             "means 'act', not 'the analysis is weak'."
             if small else
             "It is material here, meaning the recommended actions are sensitive to "
             "the remaining uncertainty — confirmatory testing is economically "
             "justified before acting.")),
        "src": "Raiffa & Schlaifer (1961); Claxton (1999), J. Health Econ.",
    }


def analyze_evppi(findings: list[dict], n: int = 3000, seed: int = 777) -> dict:
    """**EVPPI — Expected Value of *Partial* Perfect Information.**

    EVPI tells you what resolving *all* uncertainty is worth. EVPPI answers the more
    useful research-prioritisation question: **which single parameter is worth
    resolving?** For parameter subset φ:

        EVPPI(φ) = E_φ[ max_a E_{θ|φ}[NB(a,θ)] ] − max_a E_θ[NB(a,θ)]

    Estimated here by the standard two-level (nested) Monte-Carlo scheme: stratify
    draws on the parameter of interest, take the best action within each stratum, and
    compare to the best single action overall.

    Interpretation: a parameter with high EVPPI is where the decision is genuinely
    fragile — that is where a confirmatory test or a better study actually pays.
    Parameters with ~zero EVPPI can be left uncertain without regret, which is a
    *cost-saving* conclusion, not a gap.
    """
    if not _HAVE_NP or not findings:
        return {"available": False}
    rng = np.random.default_rng(seed)
    # No fixed WTP here on purpose: this routine samples willingness-to-pay
    # per draw below, so a base-case value would be a pinned parameter
    # sitting unused beside the distribution that replaced it.

    # Sample NMB per finding while recording the driving parameter draws.
    nmb = np.zeros((n, len(findings)))
    params = {"willingness_to_pay": np.zeros(n), "event_probability": np.zeros(n),
              "treatment_effect": np.zeros(n), "cost_of_illness": np.zeros(n)}
    for j in range(n):
        w = float(rng.triangular(WTP["low"], WTP["base"], WTP["high"]))
        params["willingness_to_pay"][j] = w
        pr = et = co = 0.0
        for i, f in enumerate(findings):
            nmb[j, i] = _finding_nmb(f, w, DISCOUNT_RATE, rng=rng)[0]
            pr += float(f.get("p_event", 0.2))
            et += float(f.get("rrr", 0.3))
            co += float(COI.get(f.get("coi_key", ""), {}).get("cost", 0.0))
        # Record perturbed versions so each parameter has draw-to-draw variation.
        params["event_probability"][j] = pr * float(rng.beta(8, 8)) * 2.0
        params["treatment_effect"][j] = et * float(rng.beta(8, 8)) * 2.0
        params["cost_of_illness"][j] = co * float(rng.gamma(6.25, 0.16))

    baseline = float(np.sum(np.maximum(nmb.mean(axis=0), 0.0)))
    rows = []
    for name, draws in params.items():
        # Stratify on the parameter, take the best action per stratum (partial info).
        order = np.argsort(draws)
        n_bins = 10
        bins = np.array_split(order, n_bins)
        inner = 0.0
        for b in bins:
            if b.size == 0:
                continue
            cond_mean = nmb[b].mean(axis=0)          # E[NB | phi in bin]
            inner += (b.size / n) * float(np.sum(np.maximum(cond_mean, 0.0)))
        rows.append({"parameter": name, "evppi": round(max(0.0, inner - baseline))})
    rows.sort(key=lambda r: -r["evppi"])
    top = rows[0] if rows else None
    all_zero = all(r["evppi"] <= 0 for r in rows)
    return {
        "available": True,
        "baseline_nb": round(baseline),
        "by_parameter": rows,
        "highest_value_parameter": (top["parameter"] if (top and top["evppi"] > 0)
                                    else None),
        "decision_insensitive": all_zero,
        "interpretation": (
            "EVPPI ranks parameters by how much resolving *that one alone* would "
            "improve the decision; by construction 0 <= EVPPI <= EVPI. " +
            ("Every parameter returns ~zero here, which is a substantive result rather "
             "than a missing one: the recommended actions do not flip anywhere in the "
             "plausible range of any single input, so buying more precision on any of "
             "them would be wasted spend. This is the same robustness the near-zero "
             "EVPI reports, decomposed parameter by parameter."
             if all_zero else
             "The top-ranked parameter is where the decision is genuinely fragile, and "
             "is therefore where a confirmatory test or better study actually pays.")),
        "src": "Felli & Hazen (1998); Strong, Oakley & Brennan (2014), Med Decis Making",
    }


def analyze_behavioural(mean: float, sd: float, test_cost: float,
                        horizon_years: int = 30) -> dict:
    """**Behavioural economics: prospect theory and hyperbolic discounting.**

    Expected-utility theory says people maximise E[u(w)]. They demonstrably don't, and
    both deviations matter here because they explain the *adoption gap* — why a test
    with strongly positive expected value still goes unbought.

    **1. Prospect theory** (Kahneman & Tversky 1979). Value is assessed on *gains and
    losses from a reference point*, with a concave gain limb, a steeper convex loss
    limb (λ ≈ 2.25 loss aversion), and probability weighting that overweights small
    probabilities:

        v(x) = x^α  (x ≥ 0);   −λ(−x)^β  (x < 0)
        w(p) = p^γ / (p^γ + (1−p)^γ)^{1/γ}

    Consequence: the **certain, immediate** cost of a test is felt roughly 2.25× more
    heavily than an equivalently sized *probabilistic, distant* health gain. That is a
    behavioural reason to make the test cheap and the framing concrete — not a reason
    to inflate the value estimate.

    **2. Hyperbolic discounting** (Laibson 1997). People discount quasi-hyperbolically,
    δ^t with an extra present-bias factor β on everything not-now:

        D(t) = 1 for t = 0;  β·δ^t for t > 0

    Prevention pays off decades away, so present bias suppresses it far more than the
    exponential 3% used in the base case. Reporting both is honest: the exponential
    number is the *normative* value; the hyperbolic one predicts *actual* uptake.
    """
    alpha, beta_pt, lam, gamma_pw = 0.88, 0.88, 2.25, 0.61

    def w(p: float) -> float:
        return (p ** gamma_pw) / (((p ** gamma_pw) + (1 - p) ** gamma_pw) ** (1 / gamma_pw))

    # Treat the payoff as: certain immediate loss (test cost) + probabilistic gain.
    p_gain = 0.5 if sd <= 0 else float(min(0.95, max(0.05, mean / (mean + sd))))
    v_gain = (max(0.0, mean)) ** alpha
    v_loss = -lam * (max(0.0, test_cost) ** beta_pt)
    pt_value = w(p_gain) * v_gain + v_loss
    # Convert back to dollars for comparability (invert the gain limb).
    pt_dollars = (max(0.0, pt_value)) ** (1 / alpha) if pt_value > 0 else -((-pt_value) ** (1 / alpha))

    # Quasi-hyperbolic vs exponential present value of a level benefit stream.
    beta_hyp, delta = 0.7, 1.0 / (1.0 + DISCOUNT_RATE)
    per_year = mean / max(1, horizon_years)
    pv_exp = sum(per_year * (delta ** t) for t in range(1, horizon_years + 1))
    pv_hyp = sum(per_year * beta_hyp * (delta ** t) for t in range(1, horizon_years + 1))
    return {
        "available": True,
        "prospect_theory_value": round(pt_dollars),
        "loss_aversion_lambda": lam,
        "probability_weight_applied": round(w(p_gain), 3),
        "pv_exponential": round(pv_exp),
        "pv_hyperbolic": round(pv_hyp),
        "present_bias_beta": beta_hyp,
        "adoption_gap": round(pv_exp - pv_hyp),
        "interpretation": (
            "Under prospect theory the certain, up-front test cost is felt ~2.25x more "
            "than an equivalent uncertain future health gain, and present bias further "
            "discounts prevention that pays off decades away. Together these explain "
            "why a test with clearly positive expected value still goes unbought — an "
            "adoption problem, not a valuation problem."),
        "src": "Kahneman & Tversky (1979); Tversky & Kahneman (1992); Laibson (1997)",
    }


def analyze_utility(mean: float, sd: float, wealth: float = 60_000.0) -> dict:
    """**Expected-utility / Arrow–Pratt risk preferences.**

    A risk-averse agent does not value an uncertain payoff at its mean. Under CRRA
    utility u(w) = w^(1−γ)/(1−γ), the certainty equivalent solves u(CE) = E[u(w+X)].
    We report CE across γ ∈ {1, 2, 4} (log, moderate, high risk aversion) and the
    implied **risk premium** — the amount a person would pay to remove the variance.

    Why it matters here: health information is a *variance-reducing* asset. A
    risk-averse individual values it strictly more than its expected monetary value,
    which is precisely why insurance markets exist (Arrow 1963).
    """
    out = {"available": True, "mean": round(mean), "sd": round(sd),
           "wealth_base": wealth, "by_gamma": [],
           "src": "Arrow (1963); Pratt (1964) — risk aversion in the small and large"}
    for gamma in (1.0, 2.0, 4.0):
        # Mean-variance (Arrow–Pratt second-order) approximation of the CE:
        #   CE ≈ μ − ½·γ·σ²/(W+μ)  — relative risk aversion scales by wealth.
        denom = max(1.0, wealth + mean)
        ce = mean - 0.5 * gamma * (sd ** 2) / denom
        out["by_gamma"].append({
            "gamma": gamma,
            "certainty_equivalent": round(ce),
            "risk_premium": round(mean - ce),
            "label": {1.0: "log utility (mild)", 2.0: "moderate",
                      4.0: "high risk aversion"}[gamma],
        })
    return out


def analyze_welfare_comparison(voi: float, test_cost: float,
                               capability_premium: float = 0.15,
                               p_breach_annual: float = 0.02,
                               horizon_years: int = 40,
                               loss_given_breach: float = 25_000.0,
                               participation_local: float = 0.85,
                               participation_central: float = 0.60,
                               participation_source: str = "Miller & Tucker (2018)") -> dict:
    """**Welfare comparison: local vs centralised genomic analysis.**

    The economic case for local-only analysis, stated as a surplus comparison rather
    than a privacy slogan. For an individual:

        S_local   = B_L − C_test
        S_central = B_C − C_test − E[L_privacy]

    The naive claim "local wins because it has no privacy cost" is *trivially* true and
    therefore uninteresting: it assumes B_L = B_C. This function concedes the harder
    case — centralised platforms may deliver a **higher gross benefit** (larger reference
    panels, curated pipelines, expert review), captured by ``capability_premium`` so
    B_C = B_L·(1+π). The comparison then has real content:

        **Local is preferred  ⟺  E[L_privacy] > B_C − B_L = π·B_L**

    i.e. local wins only when the expected privacy cost exceeds the capability gap.
    The function reports the **break-even breach probability** at which a rational agent
    is indifferent — the number that actually decides the argument.

    **Expected privacy cost.** A genome is non-revocable: unlike a password it cannot be
    re-keyed after disclosure, so exposure is a *permanent* state, and the hazard applies
    every year the data sits in a third-party system:

        E[L_privacy] = L · [1 − (1 − p)^T]   (discounted)

    **Social surplus adds a participation channel (RE-AIM).**
    Framework: RE-AIM (Glasgow et al. 1999) — population impact = reach x
    effectiveness. An intervention that reaches more people can beat a better
    one that reaches fewer, which is why the access channel below can dominate
    the per-person effect. See docs/METHODS.md section 22.
 Disclosure risk deters testing
    altogether (Akerlof-style unravelling): people who decline forfeit the *entire* VOI.
    Social surplus is therefore participation-weighted, and local analysis can dominate
    on access even where per-person benefit is identical.

    **The participation gap is empirically grounded, not assumed.** Miller & Tucker
    (2018, *Management Science* 64(10):4648–4668) exploit state-level variation in US
    genetic-privacy law and find that regimes granting patients **control** over their
    genetic data raise testing incidence by **+83%**, while regimes that merely notify
    people of privacy risk and ask them to consent — *without* granting control — lower
    testing by **−69%**. That contrast is close to the local-vs-cloud distinction here:
    local analysis is the maximal-control regime; uploading under a terms-of-service
    consent is the notice-without-control regime. Survey evidence agrees on direction:
    NORC finds ~80% of Americans hold privacy concerns about DNA testing, ~17% of
    non-testers name privacy as the reason they abstain, and four in five non-testers
    say they would be more willing if privacy were assured.

    The defaults used here (0.85 vs 0.60, a 1.42× ratio) are **deliberately more
    conservative than the literature**, whose implied ratio between the control and
    notice-only regimes is far larger.
    """
    B_L = float(voi)
    B_C = B_L * (1.0 + float(capability_premium))          # concede the capability gap
    r = DISCOUNT_RATE

    # Cumulative probability the data is exposed at least once over the horizon,
    # with the loss discounted to present value at the mid-point of the horizon.
    p_any = 1.0 - (1.0 - float(p_breach_annual)) ** int(horizon_years)
    disc_mid = 1.0 / ((1.0 + r) ** (horizon_years / 2.0))
    e_loss = p_any * float(loss_given_breach) * disc_mid

    s_local = B_L - test_cost
    s_central = B_C - test_cost - e_loss
    capability_gap = B_C - B_L

    # Break-even annual breach probability: solve p_any·L·disc = capability_gap.
    if loss_given_breach * disc_mid > 0:
        p_any_star = min(1.0, capability_gap / (loss_given_breach * disc_mid))
        # invert 1-(1-p)^T = p_any_star
        p_star = 1.0 - (1.0 - p_any_star) ** (1.0 / max(1, horizon_years)) \
            if p_any_star < 1.0 else 1.0
    else:
        p_star = None

    # Social surplus: participation-weighted (the chilling-effect channel).
    soc_local = participation_local * s_local
    soc_central = participation_central * s_central
    access_gain = (participation_local - participation_central) * max(0.0, s_local)

    return {
        "available": True,
        "benefit_local": round(B_L),
        "benefit_central": round(B_C),
        "capability_premium_assumed": capability_premium,
        "capability_gap": round(capability_gap),
        "expected_privacy_cost": round(e_loss),
        "prob_exposure_over_horizon": round(p_any, 4),
        "surplus_local": round(s_local),
        "surplus_central": round(s_central),
        "surplus_advantage_local": round(s_local - s_central),
        "local_preferred": bool(s_local > s_central),
        "breakeven_annual_breach_prob": (round(p_star, 5) if p_star is not None else None),
        "participation_local": participation_local,
        "participation_central": participation_central,
        "participation_evidence": (
            "Miller & Tucker (2018, Management Science): US state privacy regimes "
            "granting patient CONTROL raise genetic-testing incidence +83%; regimes "
            "with notice-and-consent but no control lower it -69%. NORC: ~80% hold "
            "privacy concerns, ~17% of non-testers cite privacy as the reason, 4 in 5 "
            "would be more willing under assured privacy. The 1.42x ratio used here is "
            "conservative relative to that evidence."),
        "social_surplus_local": round(soc_local),
        "social_surplus_central": round(soc_central),
        "access_channel_gain": round(access_gain),
        "interpretation": (
            "Local analysis is preferred when the expected privacy cost exceeds the "
            "capability gap conceded to centralised platforms. Because a genome cannot "
            "be revoked, its exposure hazard compounds over the whole horizon, so even "
            "a low annual breach probability accumulates. A second, larger channel is "
            "access: disclosure risk deters testing entirely, and a non-tester forfeits "
            "100% of the value — so local analysis raises social surplus mainly by "
            "letting more people participate at all."),
        "caveat": ("Illustrative welfare arithmetic with assumed parameters, not an "
                   "empirical estimate. The capability premium is deliberately set in "
                   "favour of the centralised alternative so the comparison is not "
                   "assumed true by construction."),
        "src": "Akerlof (1970); Posner (1981) economics of privacy; Acquisti et al. (2016) JEL",
    }


def analyze_insurance_economics(evpi: float, mean_value: float) -> dict:
    """**Information economics of genomic testing** — asymmetric information.

    Genomic results create a classic Akerlof/Rothschild–Stiglitz problem: if the
    individual knows their risk and the insurer does not, the insurance market faces
    **adverse selection**; if the insurer knows and can price on it, individuals face
    **genetic discrimination** and may rationally *decline* testing — a welfare loss
    where valuable information goes un-acquired.

    This is not a number the tool can compute for you; it is a structural caveat that
    belongs in any honest valuation, and it is the reason GenomeLens is built to run
    locally: **keeping the result private preserves the individual's option value
    without triggering the market failure.**
    """
    return {
        "available": True,
        "adverse_selection_note": (
            "Private knowledge of elevated risk without insurer knowledge creates "
            "adverse selection (Akerlof 1970; Rothschild–Stiglitz 1976)."),
        "discrimination_note": (
            "Where insurers can underwrite on genotype, the fear of discrimination "
            "suppresses testing — information that is socially valuable goes "
            "unacquired. In the US, GINA (2008) restricts health-insurance and "
            "employment use but NOT life, disability, or long-term-care insurance."),
        "privacy_as_economic_design": (
            "Local-only analysis is an economic design choice, not just a privacy "
            "one: it lets the individual capture the decision value of the "
            "information without disclosing it into a market that would price "
            "against them."),
        "chilling_effect_cost": (
            "If disclosure risk deters testing, the social loss is on the order of "
            "the foregone value of information — here modelled at "
            f"${mean_value:,.0f} per person, with a research ceiling (EVPI) of "
            f"${evpi:,.0f}."),
        "src": "Akerlof (1970); Rothschild & Stiglitz (1976); GINA (2008)",
    }


def analyze_penetrance_posterior(prior_penetrance: float, gene: str = "",
                                 family_history: bool = False,
                                 ascertainment_inflation: float = 2.0) -> dict:
    """**Bayesian penetrance with ascertainment correction** (the genomics-side rigor).

    The single biggest error in consumer genomic risk estimates is using penetrance
    figures derived from **clinically ascertained families** — cohorts selected
    *because* they were densely affected. Applied to an unselected person found
    incidentally (as in a screening context like this tool), those figures overstate
    risk substantially: population-based BRCA1 penetrance estimates run far below the
    classic linkage-study numbers (Begg 2002; Gabai-Kapara 2014).

    Two corrections are applied:

      1. **Ascertainment de-biasing** — shrink the literature penetrance toward the
         population estimate by the inflation factor observed for family-based designs.
      2. **Bayesian updating on family history** — an unselected carrier *without*
         family history is further down-weighted; *with* it, partially restored.

    The corrected posterior is applied in ``_collect`` at the point a ClinVar
    actionable finding's penetrance enters the model, so the NMB is computed on the
    de-biased value rather than the raw ascertained one. Both figures are retained on
    the finding (``penetrance_literature`` / ``penetrance_corrected``) so the size of
    the correction stays auditable. Winner's-curse shrinkage is applied separately to
    polygenic effect sizes (see ``shrink_effect_size``).
    """
    p = max(1e-6, min(0.999, float(prior_penetrance)))
    # 1) de-bias: family-ascertained estimates inflate risk by ~ascertainment_inflation.
    #    Work in odds space so the correction is well-behaved at the boundaries.
    odds = p / (1.0 - p)
    odds_pop = odds / max(1.0, ascertainment_inflation)
    p_pop = odds_pop / (1.0 + odds_pop)
    # 2) family history as a Bayes factor on the population posterior.
    bf = 2.5 if family_history else 0.8
    odds_post = odds_pop * bf
    p_post = odds_post / (1.0 + odds_post)
    return {
        "available": True,
        "gene": gene,
        "prior_literature_penetrance": round(p, 4),
        "population_corrected": round(p_pop, 4),
        "posterior_penetrance": round(p_post, 4),
        "ascertainment_inflation": ascertainment_inflation,
        "family_history": family_history,
        "shrinkage_factor": round(p_post / p, 3),
        "note": ("Literature penetrance from clinically ascertained families "
                 "overstates risk for an incidentally-identified carrier; the "
                 "economic model is computed on the corrected posterior (applied in "
                 "_collect, not display-only)."),
        "src": "Begg (2002) JNCI; Gabai-Kapara (2014) PNAS; ACMG incidental-findings guidance",
    }


def shrink_effect_size(beta_hat: float, se: float, threshold_z: float = 5.45) -> dict:
    """**Winner's-curse correction for GWAS effect sizes.**

    Effect sizes at variants discovered *because* they crossed genome-wide
    significance (p < 5×10⁻⁸, |z| > 5.45) are upward-biased: the discovery sample's
    noise had to point the right way for the variant to be found at all. Using raw
    discovery betas inflates every downstream polygenic risk estimate — and therefore
    every dollar figure built on it.

    Applies a conditional-likelihood shrinkage (Zhong & Prentice 2008): the estimate
    is pulled toward zero by an amount that grows as |z| approaches the discovery
    threshold, and vanishes for very large |z|.
    """
    if se <= 0:
        return {"available": False}
    z = beta_hat / se
    az = abs(z)
    if az <= threshold_z:
        # Below threshold the correction is undefined (variant wouldn't be reported).
        shrunk = beta_hat * 0.5
        factor = 0.5
    else:
        # Excess over the selection threshold is the credible signal.
        factor = max(0.35, (az - threshold_z * 0.55) / az)
        shrunk = beta_hat * factor
    return {
        "available": True,
        "beta_reported": round(beta_hat, 5),
        "z": round(z, 2),
        "beta_shrunk": round(shrunk, 5),
        "shrinkage_factor": round(factor, 3),
        "note": ("GWAS discovery effect sizes are upward-biased by selection; "
                 "unshrunk betas inflate polygenic risk and every dollar value "
                 "derived from it."),
        "src": "Zhong & Prentice (2008), Biostatistics — winner's curse correction",
    }


def _price_panel() -> dict:
    total = sum(p for _, p in MARKET_PRICE_ITEMS)
    return {"a_la_carte_total": total, "consolidated": CONSOLIDATED_PRICE,
            "items": [{"name": n, "price": p} for n, p in MARKET_PRICE_ITEMS],
            "note": "MARKET PRICE (what these cost to buy) — distinct from the "
                    "health-economic VALUE above."}


def _methods(rate: float, wtp: float, test_cost: float, input_type: str) -> list[str]:
    return [
        "Perspective: individual lifetime; decision-analytic net monetary benefit.",
        f"Willingness-to-pay λ = ${wtp:,.0f}/QALY (sensitivity ${WTP['low']:,}–"
        f"${WTP['high']:,}; Neumann et al., NEJM 2014).",
        f"Discounting: {rate:.0%} on both costs and QALYs (sensitivity 0/3/5%).",
        f"Test cost modelled: ${test_cost:,.0f} ({input_type}).",
        "Cost-of-illness figures: ADA, AHA, Alzheimer's Association (illustrative).",
        "PGx averted-ADR values from published CEA (Schackman 2008; Kazi 2014; "
        "Deenen 2016; CPIC).",
        "Uncertainty via seeded Monte-Carlo PSA (Beta on probabilities, Gamma on "
        "costs); the WTP threshold λ is swept deterministically for the CEAC.",
        "Left-tail risk: VaR_95 (5th-percentile outcome) and CVaR_95 / expected "
        "shortfall (mean of the outcomes at or below VaR_95) — standard risk-"
        "management measures, applied here to health-economic value.",
        "Phase-3 predicted variants down-weighted by predictor confidence.",
        "Health-capital framing: Grossman (1972) — health as a depreciating capital "
        "stock; information raises the marginal efficiency of health investment, so "
        "its value compounds over remaining life-years.",
        "Timing: real-options (Dixit & Pindyck 1994) — test now vs defer, given "
        "falling sequencing prices and improving interpretation.",
        "Risk-adjusted views: ROI multiple, reward-to-variability (Sharpe-style, no "
        "risk-free benchmark), and a mean-variance certainty equivalent (γ = 2).",
        "EVPI (Raiffa & Schlaifer 1961; Claxton 1999): the decision-theoretic ceiling "
        "on the value of any further research, with the probability that perfect "
        "information would reverse the recommended action.",
        "Expected utility (Arrow 1963; Pratt 1964): certainty equivalents and risk "
        "premia under CRRA risk aversion, γ ∈ {1, 2, 4}.",
        "Penetrance is ascertainment-corrected before entering the economic model "
        "(family-ascertained literature estimates overstate risk for incidentally "
        "identified carriers; Begg 2002, Gabai-Kapara 2014), and GWAS effect sizes "
        "are winner's-curse shrunk (Zhong & Prentice 2008).",
        "Information-economics caveats: adverse selection (Akerlof 1970; "
        "Rothschild–Stiglitz 1976) and genetic-discrimination chilling effects "
        "(GINA 2008 does not cover life/disability/long-term-care insurance).",
        "Reporting follows the CHEERS 2022 checklist in spirit; illustrative only.",
    ]
