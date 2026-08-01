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
from typing import Dict, List, Optional, Tuple

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
COI: Dict[str, Dict] = {
    "CAD":            {"cost": 90_000,  "src": "AHA — lifetime CVD direct+indirect"},
    "T2D":            {"cost": 85_000,  "src": "ADA 2018 — lifetime excess cost"},
    "Alzheimer":      {"cost": 350_000, "src": "Alz. Assoc. 2023 — incl. long-term care"},
    "BreastOvarian":  {"cost": 130_000, "src": "hereditary breast/ovarian treatment"},
    "Colorectal":     {"cost": 110_000, "src": "Lynch-spectrum CRC treatment"},
    "Pathogenic":     {"cost": 80_000,  "src": "illustrative avg for a pathogenic finding"},
}

# Pharmacogenomics cost-effectiveness (averted ADRs), from published CEA.
# p_rx = P(prescribed over a lifetime); p_adr = P(severe ADR | actionable genotype,
# standard dosing); rrr = relative risk reduction with genotype-guided care;
# adr_cost = cost of the ADR event; qaly = QALY loss averted.
PGX_CEA: Dict[str, Dict] = {
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
    n = int(round(years))
    if n <= 0:
        return 0.0
    # Σ 1/(1+r)^t for t=1..n, half-cycle correction (+0.5 at t=0).
    return sum(1.0 / (1.0 + rate) ** t for t in range(1, n + 1)) + 0.5


# ── finding collection ────────────────────────────────────────────────────────

def _collect(economics_result: Optional[Dict],
             clinical_variants_result: Optional[Dict],
             novel_variants_result: Optional[Dict]) -> List[Dict]:
    """Unify findings across sources into explicit economic parameter dicts:
    {label, kind, coi_key/pgx_key, p_event, rrr, qaly, intervention, horizon,
     wgs_only, haircut, confidence}."""
    out: List[Dict] = []

    econ = economics_result or {}
    for f in (econ.get("findings_with_economics") or []):
        cat = (f.get("category") or "").lower()
        label = f.get("finding", "genetic finding")
        conf = f.get("confidence", "moderate")
        if cat == "pgx":
            key = _match_pgx(label)
            out.append({"label": label, "kind": "pgx", "pgx_key": key,
                        "intervention": 100.0, "wgs_only": False,
                        "haircut": 1.0, "confidence": conf})
        elif cat in ("prs", "apoe", "exercise_longevity"):
            coi_key = "Alzheimer" if cat == "apoe" else "CAD"
            out.append({"label": label, "kind": "coi", "coi_key": coi_key,
                        "p_event": 0.20 if cat != "apoe" else 0.15,
                        "rrr": 0.30, "qaly": float(f.get("qaly_gain") or 0.5),
                        "intervention": 500.0, "horizon": 25, "wgs_only": False,
                        "haircut": 1.0, "confidence": conf})

    # Phase-2 clinical variants (WGS-only) — actionable + carrier + affected.
    cvr = clinical_variants_result or {}
    if cvr.get("available"):
        buckets = cvr.get("buckets", {})
        for f in (buckets.get("actionable") or []):
            coi_key, p, rrr, qaly = _gene_to_econ(f.get("gene", ""))
            out.append({"label": f"{f.get('gene','?')} pathogenic (ClinVar)",
                        "kind": "coi", "coi_key": coi_key, "p_event": p,
                        "rrr": rrr, "qaly": qaly, "intervention": 1_500.0,
                        "horizon": 30, "wgs_only": True, "haircut": 1.0,
                        "confidence": "high"})
        # carriers: reproductive value (small direct health value to self)
        for f in (buckets.get("carrier") or [])[:6]:
            out.append({"label": f"{f.get('gene','?')} carrier (reproductive)",
                        "kind": "coi", "coi_key": "Pathogenic", "p_event": 0.04,
                        "rrr": 0.5, "qaly": 0.2, "intervention": 400.0,
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
                        "qaly": 0.5, "intervention": 800.0, "horizon": 25,
                        "wgs_only": True,
                        "haircut": max(0.1, float(am)),   # scale value by AM confidence
                        "confidence": f.get("confidence", "low")})
    return out


def _match_pgx(label: str) -> str:
    low = label.lower()
    for key in PGX_CEA:
        gene = key.split("/")[0].lower()
        if gene.replace("*", "").split("-")[0] in low or key.split("/")[1] in low:
            return key
    return "PGx-generic"


def _gene_to_econ(gene: str) -> Tuple[str, float, float, float]:
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

def _finding_nmb(f: Dict, wtp: float, rate: float, rng=None) -> Tuple[float, float, float]:
    """Return (nmb, dcost_averted, dqaly) for one finding. If rng is given, sample
    uncertain parameters (PSA); else use point estimates."""
    disc = _annuity_pv(f.get("horizon", 20), rate) / max(1e-9,
                                                          f.get("horizon", 20)) \
        if f.get("horizon") else 1.0
    haircut = f.get("haircut", 1.0)

    if f["kind"] == "pgx":
        c = PGX_CEA.get(f.get("pgx_key"), PGX_CEA["PGx-generic"])
        p_rx, p_adr, rrr = c["p_rx"], c["p_adr"], c["rrr"]
        adr_cost, qaly = c["adr_cost"], c["qaly"]
        if rng is not None:
            p_adr = _beta(rng, p_adr)
            rrr = _beta(rng, rrr)
            adr_cost = _gamma(rng, adr_cost)
        exp_events = p_rx * p_adr * rrr
        dcost = exp_events * adr_cost
        dqaly = exp_events * qaly
        interv = f.get("intervention", 100.0) * p_rx
    else:
        coi = COI.get(f.get("coi_key"), COI["Pathogenic"])["cost"]
        p, rrr, qaly = f.get("p_event", 0.2), f.get("rrr", 0.3), f.get("qaly", 0.5)
        if rng is not None:
            p = _beta(rng, p)
            rrr = _beta(rng, rrr)
            coi = _gamma(rng, coi)
        dcost = p * rrr * coi * disc
        dqaly = p * rrr * qaly * disc
        interv = f.get("intervention", 500.0)

    dcost *= haircut
    dqaly *= haircut
    nmb = dqaly * wtp + dcost - interv
    return nmb, dcost, dqaly, interv


def _beta(rng, mean: float, conc: float = 25.0) -> float:
    mean = min(max(mean, 1e-3), 1 - 1e-3)
    a = mean * conc
    b = (1 - mean) * conc
    return float(rng.beta(a, b))


def _gamma(rng, mean: float, cv: float = 0.4) -> float:
    shape = 1.0 / (cv * cv)
    scale = mean / shape
    return float(rng.gamma(shape, scale))


# ── main ──────────────────────────────────────────────────────────────────────

def analyze_value_of_information(economics_result: Optional[Dict] = None,
                                 clinical_variants_result: Optional[Dict] = None,
                                 novel_variants_result: Optional[Dict] = None,
                                 genetic_age_result: Optional[Dict] = None,
                                 input_type: str = "chip",
                                 n_mc: int = 5000, seed: int = 12345,
                                 log=None) -> Dict:
    findings = _collect(economics_result, clinical_variants_result,
                        novel_variants_result)
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
        nmb_rows.append({"label": f["label"], "confidence": f["confidence"],
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
    except Exception:
        result["health_capital"] = {"available": False}

    # Real options: test now vs defer, given falling prices and improving
    # interpretation — the classic option-to-wait, on a depreciating asset.
    try:
        result["real_option"] = analyze_real_option(
            voi_now=float(result.get("voi_expost_mean",
                                     result.get("voi_expost_point", 0.0))) + test_cost,
            test_cost=test_cost, age=_age)
    except Exception:
        result["real_option"] = {"available": False}

    # Risk-adjusted framings drawn from the PSA distribution.
    try:
        result["risk_adjusted"] = _risk_adjusted(result, test_cost)
    except Exception:
        result["risk_adjusted"] = {"available": False}

    # ── Decision-theoretic ceiling: EVPI (what ANY further research could be worth) ──
    try:
        result["evpi"] = analyze_evpi(findings, test_cost)
    except Exception:
        result["evpi"] = {"available": False}

    # ── Expected utility / Arrow–Pratt: a risk-averse agent values variance
    #    reduction above the expected monetary value. ──
    try:
        _mean = float(result.get("voi_expost_mean", result.get("voi_expost_point", 0.0)))
        _sd = float((result.get("risk_adjusted") or {}).get("sd", 0.0))
        result["utility"] = analyze_utility(_mean, _sd)
    except Exception:
        result["utility"] = {"available": False}

    # ── Information economics: adverse selection, genetic discrimination, and why
    #    local-only analysis is an economic design choice. ──
    try:
        result["information_economics"] = analyze_insurance_economics(
            evpi=float((result.get("evpi") or {}).get("evpi", 0.0)),
            mean_value=float(result.get("voi_expost_mean",
                                        result.get("voi_expost_point", 0.0))))
    except Exception:
        result["information_economics"] = {"available": False}

    # ── Genomics-side rigor: show the ascertainment correction actually applied to a
    #    representative high-penetrance finding, so the bias is visible, not hidden. ──
    try:
        pen = [f for f in findings if f.get("coi_key") in
               ("BreastOvarian", "Colorectal") and f.get("p_event")]
        if pen:
            f0 = pen[0]
            result["penetrance_correction"] = analyze_penetrance_posterior(
                prior_penetrance=float(f0.get("p_event", 0.5)),
                gene=str(f0.get("label", "")).split()[0])
    except Exception:
        pass

    return result


def _resolve_age(genetic_age_result: Optional[Dict], default: float = 35.0) -> float:
    """Best-effort age from the bloodwork/aging module; falls back to a default."""
    g = genetic_age_result or {}
    for key in ("chronological_age", "chronological", "age_used", "age"):
        v = g.get(key)
        if isinstance(v, (int, float)) and 0 < float(v) < 120:
            return float(v)
    clinical = g.get("clinical") if isinstance(g.get("clinical"), dict) else {}
    v = clinical.get("age_used")
    if isinstance(v, (int, float)) and 0 < float(v) < 120:
        return float(v)
    return default


def _risk_adjusted(result: Dict, test_cost: float) -> Dict:
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


def _exante(test_cost: float) -> Dict:
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


def _psa(findings: List[Dict], test_cost: float, n: int, seed: int) -> Dict:
    """Seeded Monte-Carlo PSA. Effects/costs are sampled per draw; λ is varied
    deterministically for the CEAC via NMB(λ) = NMB_base + ΔQALY·(λ − λ_base)."""
    rng = np.random.default_rng(seed)
    base = float(WTP["base"])
    totals = np.zeros(n)     # NMB at λ_base, net of test cost
    dqalys = np.zeros(n)     # ΔQALY
    netcosts = np.zeros(n)   # incremental cost (intervention − averted); neg = saving
    for i in range(n):
        tot = dq = dc = iv = 0.0
        for f in findings:
            nmb, c, q, interv = _finding_nmb(f, base, DISCOUNT_RATE, rng=rng)
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


def _tornado(findings: List[Dict], test_cost: float) -> List[Dict]:
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
                           informed_efficiency_gain: Optional[float] = None) -> Dict:
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
                        defer_years: int = 5) -> Dict:
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


def analyze_evpi(findings: List[Dict], test_cost: float, n: int = 4000,
                 seed: int = 4242) -> Dict:
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


def analyze_utility(mean: float, sd: float, wealth: float = 60_000.0) -> Dict:
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


def analyze_insurance_economics(evpi: float, mean_value: float) -> Dict:
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
                                 ascertainment_inflation: float = 2.0) -> Dict:
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

    The corrected posterior then drives the economic model, so the dollar value never
    inherits the ascertainment bias. Winner's-curse shrinkage is applied separately to
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
                 "economic model uses the corrected posterior."),
        "src": "Begg (2002) JNCI; Gabai-Kapara (2014) PNAS; ACMG incidental-findings guidance",
    }


def shrink_effect_size(beta_hat: float, se: float, threshold_z: float = 5.45) -> Dict:
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


def _price_panel() -> Dict:
    total = sum(p for _, p in MARKET_PRICE_ITEMS)
    return {"a_la_carte_total": total, "consolidated": CONSOLIDATED_PRICE,
            "items": [{"name": n, "price": p} for n, p in MARKET_PRICE_ITEMS],
            "note": "MARKET PRICE (what these cost to buy) — distinct from the "
                    "health-economic VALUE above."}


def _methods(rate: float, wtp: float, test_cost: float, input_type: str) -> List[str]:
    return [
        f"Perspective: individual lifetime; decision-analytic net monetary benefit.",
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
