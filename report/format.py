"""One place where numbers become strings.

WHY THIS EXISTS. Before this module the repository defined five separate money
formatters — ``money()`` three times and ``_money()`` once in
``report/renderers.py``, plus ``_money()`` in ``econ/health_economics.py`` —
and they did not agree:

===========================  ==================  ================  ===============
where                        negative            missing value     rounding
===========================  ==================  ================  ===============
renderers.py:2132            ``-$216``           ``&mdash;``       ``round()``
renderers.py:2321            ``-$216``           ``—``             ``round()``
renderers.py:2570            ``-$216``           ``—``             ``round()``
renderers.py:5397            ``$-216``           ``—``             ``,.0f``
health_economics.py:2346     ``-$216``           *raises*          ``round()``
===========================  ==================  ================  ===============

Three real defects came out of that table, all of them visible in generated
output rather than theoretical:

1. **``$-216`` instead of ``-$216``.** The value-of-information per-finding
   table used the fourth formatter, so the one place in the report where a
   finding can carry a negative net monetary benefit is also the one place that
   renders the minus sign in the wrong position.
2. **Two rounding modes.** ``round()`` is banker's rounding (half to even);
   ``,.0f`` rounds half away from zero. The same half-dollar amount could
   therefore appear as ``$0`` in one table and ``$1`` in another. This module
   uses ``ROUND_HALF_UP`` throughout, which is the convention money is expected
   to follow and which no longer varies by call site.
3. **One formatter raises on ``None``** where the other four degrade to a dash.

FORMATTING IS PRESENTATION ONLY. Nothing here changes a stored value. Callers
pass full precision and receive a string; the payload keeps the number.
"""
from __future__ import annotations

from decimal import ROUND_HALF_UP, Decimal, InvalidOperation

# The dash shown when a quantity is genuinely absent. A literal character rather
# than an HTML entity: these strings are escaped by the renderer on the way out,
# and an entity that survives escaping renders as "&amp;mdash;".
MISSING = "—"

__all__ = [
    "MISSING",
    "count",
    "healthy_days",
    "money",
    "months",
    "one_in",
    "percentage",
    "probability",
    "qaly",
    "ratio",
    "signed_money",
    "simulation_count",
    "years",
]


def _dec(v) -> Decimal | None:
    """Coerce to Decimal, or None when the value is not a finite number."""
    if v is None or isinstance(v, bool):
        return None
    try:
        d = Decimal(str(v))
    except (InvalidOperation, ValueError, TypeError):
        return None
    return d if d.is_finite() else None


def _round(d: Decimal, places: int = 0) -> Decimal:
    return d.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP)


# ── money ─────────────────────────────────────────────────────────────────────

def money(v, *, missing: str = MISSING) -> str:
    """``$1,234`` / ``-$1,234``. The sign leads; the currency symbol follows it.

    Whole dollars. Cents are noise at every magnitude this report deals in, and
    showing them invites the reader to believe a modeled estimate is precise to
    the penny.
    """
    d = _dec(v)
    if d is None:
        return missing
    n = _round(d)
    return f"-${abs(n):,}" if n < 0 else f"${n:,}"


def signed_money(v, *, missing: str = MISSING) -> str:
    """``+$1,234`` / ``-$1,234``, for quantities that are *changes*.

    An incremental cost of ``$4,223`` and one of ``-$4,223`` mean opposite
    things, and the reader should not have to infer the sign from a caption. Use
    this for anything incremental; use :func:`money` for levels.
    """
    d = _dec(v)
    if d is None:
        return missing
    n = _round(d)
    if n == 0:
        return "$0"
    return f"-${abs(n):,}" if n < 0 else f"+${n:,}"


# ── health ────────────────────────────────────────────────────────────────────

def qaly(v, *, missing: str = MISSING) -> str:
    """Quality-adjusted life-years, at precision that cannot hide a real value.

    THE BUG THIS FIXES. Per-finding QALY gains were rendered at two decimals, so
    a finding contributing 0.00065 QALYs displayed as ``0.00`` — indistinguishable
    from a finding contributing nothing at all, in a table whose whole purpose is
    to separate those two cases. A reader could not tell "too small to matter"
    from "not modeled".

    So precision adapts: three decimals normally, more when three would round a
    genuinely non-zero quantity to zero, and an explicit ``<0.00001`` floor below
    the point where more digits stop being meaningful.
    """
    d = _dec(v)
    if d is None:
        return missing
    if d == 0:
        return "0.000"
    for places in (3, 4, 5):
        r = _round(d, places)
        if r != 0:
            return f"{r:.{places}f}"
    return "<0.00001" if d > 0 else ">-0.00001"


def healthy_days(v, *, missing: str = MISSING) -> str:
    """QALYs restated as days, because a QALY is not an everyday unit.

    ``0.032 QALYs`` means little on its own; "about 12 additional days of healthy
    life" is the same number in a unit a reader already owns. Rounded to whole
    days and always hedged by the caller's surrounding prose — this is a
    population-average expectation, not a promise to an individual.
    """
    d = _dec(v)
    if d is None:
        return missing
    days = _round(d * Decimal("365.25"))
    return f"{days:,} day" + ("" if abs(days) == 1 else "s")


# ── proportions ───────────────────────────────────────────────────────────────

def percentage(v, *, places: int = 1, missing: str = MISSING) -> str:
    """A proportion in 0..1 rendered as a percentage: ``0.327`` -> ``32.7%``."""
    d = _dec(v)
    if d is None:
        return missing
    return f"{_round(d * 100, places):.{places}f}%"


def probability(v, *, missing: str = MISSING) -> str:
    """A probability as whole percent, with the tails protected.

    ``100%`` and ``0%`` are claims of certainty. A probabilistic analysis that
    ran 1,500 iterations and saw no failure has evidence for ">99%", not for
    "100%", so values that round to the extremes are reported as bounds instead.
    """
    d = _dec(v)
    if d is None:
        return missing
    pct = _round(d * 100)
    if pct >= 100 and d < 1:
        return ">99%"
    if pct <= 0 and d > 0:
        return "<1%"
    return f"{pct:.0f}%"


def simulation_count(p_value, n_iterations, *, missing: str = MISSING) -> str:
    """``all 1,500 of 1,500 simulations`` — the count, not a rounded percentage.

    A probability of 0.9987 over 1,500 runs is not "100%", and reporting it that
    way claims a certainty the analysis does not have. Saying how many runs did
    what is both more precise and easier to check: a reader can see the
    denominator and judge it. :func:`probability` supplies the secondary label.
    """
    if not n_iterations:
        return missing
    d = _dec(p_value)
    if d is None:
        return missing
    n = int(n_iterations)
    hits = int(_round(d * Decimal(n)))
    if hits >= n:
        return f"all {n:,} of {n:,} simulations"
    return f"{hits:,} of {n:,} simulations"


def one_in(v, *, missing: str = MISSING) -> str:
    """``1 in 56`` from a probability, for frequencies a reader can picture.

    The caller is responsible for saying *one in fifty-six what*. This report
    carries at least two unrelated "1 in N" quantities — the chance sequencing
    finds something an array missed, and the number needed to screen to avert
    one case of a given condition — and they have been mistaken for each other.
    """
    d = _dec(v)
    if d is None or d <= 0:
        return missing
    return f"1 in {_round(Decimal(1) / d):,}"


# ── everything else ───────────────────────────────────────────────────────────

def ratio(numerator, denominator, *, missing: str = MISSING) -> str:
    """``2.4x``. Never rendered as ``N:1`` and never called a return on investment.

    The report's health-economic value already includes monetized health, so a
    ratio built from it is not money coming back to anyone. Callers must label
    what is being divided.
    """
    n, dd = _dec(numerator), _dec(denominator)
    if n is None or dd is None or dd == 0:
        return missing
    return f"{_round(n / dd, 1):.1f}x"


def months(v, *, missing: str = MISSING) -> str:
    d = _dec(v)
    if d is None:
        return missing
    r = _round(d, 1)
    return f"{r:.1f} month" + ("" if r == 1 else "s")


def years(v, *, places: int = 0, missing: str = MISSING) -> str:
    d = _dec(v)
    if d is None:
        return missing
    r = _round(d, places)
    return f"{r:.{places}f} year" + ("" if r == 1 else "s")


def count(v, *, missing: str = MISSING) -> str:
    d = _dec(v)
    if d is None:
        return missing
    return f"{_round(d):,}"
