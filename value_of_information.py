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
    return result


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
    return {
        "psa_available": True,
        "voi_expost_mean": round(float(np.mean(totals))),
        "voi_ci_low": round(float(np.percentile(totals, 2.5))),
        "voi_ci_high": round(float(np.percentile(totals, 97.5))),
        "prob_cost_effective": round(float(np.mean(totals > 0)), 3),
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
        "Phase-3 predicted variants down-weighted by predictor confidence.",
        "Reporting follows the CHEERS 2022 checklist in spirit; illustrative only.",
    ]
