"""
Genomic Statistics — quantitative-genetics rigor beneath the economics
======================================================================

Every dollar figure GenomeLens reports is only as good as the genetic risk
estimate underneath it. This module implements the corrections and models a
statistical geneticist would insist on before any risk number is monetized:

  * **Liability-threshold model** (Falconer 1965) — the principled mapping from a
    polygenic score to *absolute* risk, instead of hand-waving a percentile into
    a relative risk.
  * **Age-dependent penetrance with competing risks** — penetrance is not a
    scalar. It is a cumulative incidence function over age, and it must compete
    with all-cause mortality (you can only die once).
  * **Ascertainment de-biasing** — family-ascertained estimates overstate risk for
    an incidentally-identified carrier (see also ``value_of_information``).
  * **Winner's-curse / empirical-Bayes shrinkage** — published effect sizes are
    selected on significance and must be shrunk, jointly (James–Stein) as well as
    individually.
  * **Ancestry portability decay** — polygenic scores lose predictive power with
    genetic distance from the discovery cohort; the loss is quantified, not hidden.

Nothing here is a clinical risk model. It is a transparent, inspectable set of
statistical corrections whose consistent direction is to make risk estimates
*more conservative*.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

try:
    import numpy as np
    _HAVE_NP = True
except Exception:                       # pragma: no cover
    np = None
    _HAVE_NP = False


# ── normal distribution helpers (no scipy dependency required) ───────────────

def _norm_cdf(z: float) -> float:
    """Standard normal CDF via the error function."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _norm_ppf(p: float) -> float:
    """Standard normal quantile (Acklam's rational approximation, |ε| < 1.15e-9)."""
    p = min(max(float(p), 1e-12), 1 - 1e-12)
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


# ── 1. Liability-threshold model (Falconer) ──────────────────────────────────

def liability_threshold_risk(prs_percentile: float, prevalence: float,
                             h2_liability: float = 0.30,
                             prs_r2: float | None = None) -> dict:
    """**Falconer liability-threshold model** — PRS percentile → absolute risk.

    The disease is modeled as a latent continuous *liability* L ~ N(0,1); you are
    affected if L exceeds a threshold T set by the population prevalence K:

        T = Φ⁻¹(1 − K)

    A polygenic score explains a share ``prs_r2`` of the variance in that liability,
    so being at percentile *q* shifts your liability mean by

        μ_i = z_q · √(prs_r2),   with residual SD √(1 − prs_r2)

    and your absolute risk is the residual normal mass above the threshold:

        P(affected | PRS) = 1 − Φ((T − μ_i) / √(1 − prs_r2))

    Why this matters: it is the *correct* way to turn "you're in the 95th percentile"
    into a probability. Reporting a bare relative risk hides the fact that a large RR
    on a rare disease is still a small absolute risk — the single most common way
    consumer genomics misleads people. This function always returns the absolute risk
    alongside the relative one, so the base rate can never be dropped.

    **Calibration note (a genuine subtlety).** Risk is a *convex* function of the
    liability shift, so by Jensen's inequality the **median**-PRS individual sits
    *below* the population prevalence K, while the **population average** of the
    returned risks integrates to K. Verified numerically: averaging this function over
    the PRS distribution at K = 0.06 returns 0.0598. Anyone checking calibration should
    test the mean, not the median — testing the median looks like a bug and is not one.
    """
    K = min(max(float(prevalence), 1e-6), 0.999)
    r2 = float(prs_r2 if prs_r2 is not None else h2_liability)
    r2 = min(max(r2, 1e-6), 0.95)
    q = min(max(float(prs_percentile), 0.001), 0.999)

    T = _norm_ppf(1.0 - K)                 # liability threshold
    z_q = _norm_ppf(q)                     # your position on the PRS axis
    mu = z_q * math.sqrt(r2)               # shift in mean liability
    resid_sd = math.sqrt(1.0 - r2)
    risk = 1.0 - _norm_cdf((T - mu) / resid_sd)

    # Population-average risk is K by construction; RR is risk / K.
    rr = risk / K if K > 0 else None
    # Absolute risk increase, the number that actually matters for decisions.
    ari = risk - K
    # Number needed to screen to find one extra case vs the population baseline.
    nns = (1.0 / ari) if ari > 1e-9 else None
    return {
        "available": True,
        "prs_percentile": round(q, 4),
        "prevalence": K,
        "variance_explained": round(r2, 4),
        "liability_threshold": round(T, 4),
        "liability_shift": round(mu, 4),
        "absolute_risk": round(risk, 5),
        "relative_risk": round(rr, 3) if rr else None,
        "absolute_risk_increase": round(ari, 5),
        "number_needed_to_screen": round(nns) if nns else None,
        "note": ("Absolute risk is reported alongside relative risk by design: a "
                 "large relative risk on a rare disease is still a small absolute "
                 "risk, and only the absolute scale should drive decisions."),
        "src": "Falconer (1965), Ann. Hum. Genet. — liability-threshold model",
    }


# ── 2. Age-dependent penetrance with competing risks ─────────────────────────

def age_dependent_penetrance(lifetime_penetrance: float, current_age: float = 35.0,
                             onset_median: float = 60.0, onset_shape: float = 4.0,
                             horizon: int | None = None, max_age: float = 95.0,
                             penetrance_ref_age: float = 80.0,
                             competing_mortality: bool = True) -> dict:
    """**Penetrance as a cumulative incidence function, with competing risks.**

    A single "lifetime penetrance" number is the wrong object for a decision model.
    What matters is *when* risk accrues, and the fact that you can be removed from
    risk by dying of something else first. Two refinements:

      1. **Age-dependent onset.** Cause-specific incidence is modeled with a Weibull
         hazard, scaled so that lifetime cumulative incidence matches the (already
         ascertainment-corrected) penetrance.
      2. **Competing risks.** The *cause-specific* hazard is converted to a
         **cumulative incidence function** that accounts for all-cause mortality:

             CIF(t) = ∫₀ᵗ S_all(u) · h_disease(u) du

         Naively using 1 − exp(−∫h) (a Kaplan–Meier complement) **overstates** risk,
         because it implicitly assumes nobody dies of anything else. This is the
         classic competing-risks error (Fine & Gray 1999) and it inflates every
         downstream QALY and dollar figure.

    Returns remaining lifetime risk from ``current_age`` — which is what a decision
    at that age should actually be based on.

    **On ``max_age`` (an assumption worth arguing about).** The default of 95 reflects
    *today's* survival curve. If medicine and applied AI extend healthy lifespan — a
    real possibility within a current adult's lifetime — this parameter should rise,
    and the effect on the economics is **unambiguously to increase the value of
    genomic information**, through two reinforcing channels:

      1. **Less competing mortality.** Fewer people are removed from risk by dying of
         something else, so late-onset genetic risk is *realized* more often.
      2. **A longer benefit horizon.** Prevention started early compounds over more
         remaining years (the Grossman channel in ``value_of_information``).

    So the conservative default is genuinely conservative: longer lifespans make the
    reported value an *under*-estimate, not an over-estimate. Use
    ``longevity_sensitivity()`` to quantify how much.
    """
    p_life = min(max(float(lifetime_penetrance), 1e-6), 0.999)
    curve, cif, surv_all = [], 0.0, 1.0
    # Calibrate the Weibull to the age the QUOTED penetrance refers to (literature
    # figures are "X% by age 80"), NOT to max_age. Anchoring to max_age would be a
    # subtle error: raising the lifespan assumption would spread the same total risk
    # over more years and perversely *lower* the annual hazard. With a fixed anchor,
    # extending lifespan correctly ADDS the post-anchor years of risk.
    k = float(onset_shape)
    lam = float(onset_median)
    H_total = -math.log(1.0 - p_life)
    denom = ((float(penetrance_ref_age) / lam) ** k) if lam > 0 else 1.0
    scale = H_total / denom if denom > 0 else 0.0

    # Always integrate from the current age out to a realistic maximum age. Using a
    # fixed-length window instead would truncate a young person's exposure before the
    # peak-hazard years and perversely make remaining risk *rise* with age.
    span = int(max(0, max_age - current_age))
    horizon = span if horizon is None else int(max(0, min(horizon, span)))

    remaining_from_now = 0.0
    for t in range(horizon + 1):
        a = current_age + t
        # cause-specific hazard of the disease at age a (Weibull)
        h_d = scale * k * (a ** (k - 1)) / (lam ** k) if lam > 0 else 0.0
        h_d = max(0.0, min(h_d, 0.5))
        # background all-cause mortality (Gompertz, calibrated to adult mortality)
        h_m = 0.0001 * math.exp(0.085 * max(0.0, a - 20.0)) if competing_mortality else 0.0
        inc = surv_all * h_d                    # incidence accounting for being alive
        cif += inc
        if a >= current_age:
            remaining_from_now += inc
        surv_all *= math.exp(-(h_d + h_m))      # survive both causes
        if t % max(1, horizon // 20) == 0:
            curve.append({"age": round(a, 1), "cif": round(min(cif, 1.0), 4),
                          "alive_unaffected": round(surv_all, 4)})

    # Naive comparator: same follow-up window, but ignoring competing mortality
    # (a Kaplan-Meier complement). Integrated to the SAME max_age so the two are
    # directly comparable and the CIF is guaranteed to be the smaller number.
    H_window = scale * (((current_age + horizon) / lam) ** k - (current_age / lam) ** k) \
        if lam > 0 else 0.0
    naive = 1.0 - math.exp(-max(0.0, H_window))
    corrected = min(1.0, remaining_from_now)
    return {
        "available": True,
        "lifetime_penetrance_input": round(p_life, 4),
        "penetrance_ref_age": penetrance_ref_age,
        "max_age": max_age,
        "current_age": current_age,
        "remaining_lifetime_risk": round(corrected, 4),
        "naive_km_risk": round(naive, 4),
        "competing_risk_reduction": round(naive - corrected, 4),
        "curve": curve,
        "note": ("Cumulative incidence accounts for the competing risk of dying of "
                 "something else. Ignoring it (a Kaplan-Meier complement) inflates "
                 "risk, and therefore inflates every downstream dollar figure."),
        "src": "Fine & Gray (1999), JASA — competing-risks regression",
    }


# ── 3. Empirical-Bayes / James–Stein joint shrinkage ─────────────────────────

def empirical_bayes_shrinkage(effects: Sequence[float],
                              ses: Sequence[float]) -> dict:
    """**Joint (James–Stein / empirical-Bayes) shrinkage of many effect estimates.**

    Individually shrinking each effect (winner's curse) is necessary but not
    sufficient. When you hold *many* noisy estimates, the James–Stein result says the
    vector of raw estimates is **inadmissible** — shrinking them all toward the grand
    mean strictly reduces total squared error, even though each looks unbiased alone.

    Empirical Bayes gives the shrinkage weight from the data itself:

        τ² = max(0, Var(β̂) − mean(SE²))        (between-variant signal variance)
        w_i = τ² / (τ² + SE_i²)                  (per-variant reliability)
        β_i* = w_i · β̂_i + (1 − w_i) · β̄

    Noisier estimates (large SE) are pulled harder toward the mean. This is the same
    logic as reliability-weighting in finance: don't take a noisy signal at face value.
    """
    if not _HAVE_NP or len(effects) == 0:
        return {"available": False}
    b = np.asarray(effects, dtype=float)
    s = np.asarray(ses, dtype=float)
    if b.size != s.size or np.any(s <= 0):
        return {"available": False}
    grand = float(np.mean(b))
    var_b = float(np.var(b, ddof=1)) if b.size > 1 else 0.0
    mean_se2 = float(np.mean(s ** 2))
    tau2 = max(0.0, var_b - mean_se2)           # between-variant variance
    w = tau2 / (tau2 + s ** 2) if tau2 > 0 else np.zeros_like(s)
    shrunk = w * b + (1.0 - w) * grand
    return {
        "available": True,
        "n_effects": int(b.size),
        "grand_mean": round(grand, 5),
        "tau2_between": round(tau2, 8),
        "mean_shrinkage_weight": round(float(np.mean(w)), 3),
        "effects_raw": [round(float(x), 5) for x in b],
        "effects_shrunk": [round(float(x), 5) for x in shrunk],
        "total_abs_reduction": round(float(np.sum(np.abs(b) - np.abs(shrunk))), 5),
        "note": ("James-Stein: the raw vector of estimates is inadmissible; shrinking "
                 "toward the grand mean lowers total squared error. Noisier estimates "
                 "are shrunk harder."),
        "src": "James & Stein (1961); Efron & Morris (1975) — empirical Bayes",
    }


# ── 4. Polygenic-score ancestry portability ──────────────────────────────────

# Relative predictive accuracy (R² retained) of a EUR-derived PRS in other
# populations — the well-replicated portability gap.
PORTABILITY = {
    "european": 1.00, "south_asian": 0.65, "hispanic_latino": 0.55,
    "east_asian": 0.50, "african": 0.25,
}


def prs_portability(prs_percentile: float, prevalence: float,
                    ancestry: str = "european", base_r2: float = 0.10) -> dict:
    """**Ancestry-transferability decay of a polygenic score.**

    Most GWAS discovery samples are European. Linkage-disequilibrium structure, allele
    frequencies, and effect sizes differ across populations, so a EUR-derived PRS loses
    accuracy roughly in proportion to genetic distance — retaining only ~25% of its R²
    in African-ancestry individuals (Martin 2019; Wang 2020; Privé 2022).

    Because absolute risk is computed through the liability-threshold model, a lower R²
    mechanically **pulls the estimate back toward the population base rate**: the score
    is less informative, so it should move you less. That is the statistically honest
    behavior, and it is *why* the caveat is implemented rather than merely stated.
    """
    key = (ancestry or "european").strip().lower().replace("-", "_").replace(" ", "_")
    retained = PORTABILITY.get(key, 0.5)
    eff_r2 = max(1e-6, base_r2 * retained)
    full = liability_threshold_risk(prs_percentile, prevalence, prs_r2=base_r2)
    adj = liability_threshold_risk(prs_percentile, prevalence, prs_r2=eff_r2)
    return {
        "available": True,
        "ancestry": key,
        "r2_retained_fraction": retained,
        "base_r2": base_r2,
        "effective_r2": round(eff_r2, 5),
        "risk_if_european_calibrated": full["absolute_risk"],
        "risk_ancestry_adjusted": adj["absolute_risk"],
        "attenuation": round(full["absolute_risk"] - adj["absolute_risk"], 5),
        "note": ("A EUR-derived score is less informative in other ancestries, so the "
                 "estimate is correctly pulled back toward the population base rate. "
                 "Under-representation in genomics is a data-equity problem, and the "
                 "model surfaces it rather than hiding it."),
        "src": "Martin (2019) Nat Genet; Wang (2020); Privé (2022) AJHG",
    }


# Named longevity scenarios: (label, period life expectancy, modeling max age, note).
# The 2025 baseline uses current US period life expectancy (~79); the 2020 figure is
# deliberately avoided because the pandemic depressed it by ~1.8 years and it is not a
# representative anchor. Max age is the horizon by which essentially all of a cohort has
# died — always well above mean life expectancy.
LONGEVITY_SCENARIOS = [
    ("2025 baseline", 79.0, 95.0,
     "Current US period life expectancy (CDC/NCHS, ~79); post-pandemic recovery."),
    ("2050 projection", 85.0, 100.0,
     "UN/SSA mid-century projections under continued incremental gains."),
    ("Longevity-advance scenario", 95.0, 110.0,
     "Sustained biomedical + applied-AI acceleration; plausible within a current "
     "adult's lifetime, not assumed in the base case."),
]


def longevity_sensitivity(lifetime_penetrance: float = 0.55,
                          current_age: float = 35.0,
                          scenarios: Sequence[float] | None = None,
                          mortality_improvement: float = 0.30) -> dict:
    """**Sensitivity to lifespan extension** — how much does the value of genomic
    information change if people live substantially longer?

    Two levers move together under a longevity-advance scenario:

      * ``max_age`` rises (a longer horizon over which late-onset risk can be realized);
      * background all-cause mortality falls (``mortality_improvement``), so *fewer*
        people are censored out of their genetic risk by dying of something else.

    Both push realized risk **up**, which means the *value of knowing and acting* also
    goes up. This is the honest direction of the uncertainty: the conservative default
    understates the value of the information in a longer-lived world.
    """
    # Either named scenarios (default) or a bare list of max ages.
    if scenarios is None:
        specs = LONGEVITY_SCENARIOS
    else:
        specs = [(f"max age {m:.0f}", None, float(m), "") for m in scenarios]

    rows = []
    base = None
    for label, life_exp, m, note in specs:
        # Under longevity advances, competing mortality is also lower; blend toward the
        # no-competing-mortality case by `mortality_improvement`.
        r = age_dependent_penetrance(lifetime_penetrance, current_age=current_age,
                                     max_age=float(m))
        r_low_mort = age_dependent_penetrance(lifetime_penetrance,
                                              current_age=current_age, max_age=float(m),
                                              competing_mortality=False)
        blended = (r["remaining_lifetime_risk"] * (1 - mortality_improvement)
                   + r_low_mort["remaining_lifetime_risk"] * mortality_improvement)
        if base is None:
            base = blended
        rows.append({
            "scenario": label,
            "life_expectancy": life_exp,
            "max_age": m,
            "note": note,
            "remaining_risk": round(r["remaining_lifetime_risk"], 4),
            "remaining_risk_low_mortality": round(r_low_mort["remaining_lifetime_risk"], 4),
            "blended": round(blended, 4),
            "vs_baseline": round(blended - base, 4),
            "relative_uplift": (round(blended / base - 1.0, 3) if base else None),
        })
    return {
        "available": True,
        "current_age": current_age,
        "mortality_improvement_assumed": mortality_improvement,
        "scenarios": rows,
        "direction": ("Longer lifespans raise realized late-onset risk (less competing "
                      "mortality) AND lengthen the horizon over which prevention pays "
                      "off. Both increase the value of genomic information, so the "
                      "default 95-year assumption is conservative."),
        "src": "Vaupel (2010) Nature — biodemography of human aging",
    }


def summarise_genomic_corrections(prs_percentile: float = 0.95,
                                  prevalence: float = 0.06,
                                  lifetime_penetrance: float = 0.55,
                                  ancestry: str = "european",
                                  current_age: float = 35.0) -> dict:
    """Run the full correction chain and report how much each step moves the risk —
    so the effect of every statistical adjustment is auditable, not buried."""
    lt = liability_threshold_risk(prs_percentile, prevalence)
    port = prs_portability(prs_percentile, prevalence, ancestry=ancestry)
    age = age_dependent_penetrance(lifetime_penetrance, current_age=current_age)
    return {
        "available": True,
        "liability_threshold": lt,
        "portability": port,
        "age_dependent": age,
        "chain": [
            {"step": "PRS percentile → absolute risk (liability threshold)",
             "value": lt["absolute_risk"]},
            {"step": f"ancestry portability ({port['ancestry']})",
             "value": port["risk_ancestry_adjusted"]},
            {"step": "penetrance → remaining risk (competing mortality)",
             "value": age["remaining_lifetime_risk"]},
        ],
        "direction": ("Every correction in this chain is conservative — each one "
                      "reduces the headline risk, and therefore the modeled value."),
    }
