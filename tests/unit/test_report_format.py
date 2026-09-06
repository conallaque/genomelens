"""The single formatting layer.

Each test here corresponds to a divergence found between the five money
formatters this module replaced, or to a rendering defect visible in generated
output. They are regressions, not hypotheticals.
"""
from __future__ import annotations

from report import format as fmt

# ── the divergences between the old formatters ────────────────────────────────

def test_negative_money_puts_the_sign_before_the_dollar():
    """REGRESSION. The value-of-information per-finding table used the one
    formatter that produced `$-216`, and it is the only table in the report
    where a finding can carry a negative net monetary benefit."""
    assert fmt.money(-216) == "-$216"
    assert fmt.money(-216) != "$-216"


def test_rounding_is_half_up_everywhere():
    """The old formatters mixed round() (banker's, half to even) with ,.0f
    (half away from zero), so the same half-dollar amount could render as $0 in
    one table and $1 in another."""
    assert fmt.money(0.5) == "$1"
    assert fmt.money(1.5) == "$2"
    assert fmt.money(2.5) == "$3"


def test_missing_values_degrade_instead_of_raising():
    """One of the five raised TypeError on None."""
    for bad in (None, "", "n/a", float("nan"), float("inf")):
        assert fmt.money(bad) == fmt.MISSING


def test_missing_marker_is_a_character_not_an_html_entity():
    """An entity that survives the renderer's escaping renders as
    `&amp;mdash;`."""
    assert "&" not in fmt.MISSING


def test_money_is_thousands_separated_and_whole_dollars():
    assert fmt.money(1234567.89) == "$1,234,568"
    assert fmt.money(0) == "$0"


# ── signed money ──────────────────────────────────────────────────────────────

def test_signed_money_marks_direction_explicitly():
    assert fmt.signed_money(4223) == "+$4,223"
    assert fmt.signed_money(-4223) == "-$4,223"
    assert fmt.signed_money(0) == "$0"


# ── QALY precision ────────────────────────────────────────────────────────────

def test_small_but_real_qaly_never_renders_as_zero():
    """REGRESSION. Per-finding QALYs rendered at two decimals, so a finding
    contributing 0.00065 displayed as `0.00` — indistinguishable from a finding
    modeled to contribute nothing, in the table whose job is to tell those
    apart."""
    assert fmt.qaly(0.00065) not in ("0.00", "0.000")
    assert float(fmt.qaly(0.00065)) > 0


def test_qaly_adds_digits_only_when_needed():
    assert fmt.qaly(0.032) == "0.032"
    assert fmt.qaly(0.0319) == "0.032"
    assert fmt.qaly(0.0004) == "0.0004"
    assert fmt.qaly(0.00004) == "0.00004"


def test_exact_zero_qaly_is_distinguishable_from_a_tiny_one():
    assert fmt.qaly(0) == "0.000"
    assert fmt.qaly(0.0000001) == "<0.00001"


def test_healthy_days_translates_qalys_into_an_everyday_unit():
    assert fmt.healthy_days(0.032) == "12 days"
    assert fmt.healthy_days(1 / 365.25) == "1 day"


# ── probability ───────────────────────────────────────────────────────────────

def test_probability_does_not_claim_certainty_it_does_not_have():
    """1500 simulations with no failure is evidence for '>99%', not for '100%'.
    The live fixture reports p_cost_effective = 0.9987, which the old renderer
    printed as 100%."""
    assert fmt.probability(0.9987) == ">99%"
    assert fmt.probability(1.0) == "100%"
    assert fmt.probability(0.0001) == "<1%"
    assert fmt.probability(0.0) == "0%"


def test_probability_rounds_to_whole_percent_in_between():
    assert fmt.probability(0.956) == "96%"
    assert fmt.probability(0.5) == "50%"


def test_percentage_keeps_one_decimal_by_default():
    assert fmt.percentage(0.628) == "62.8%"
    assert fmt.percentage(0.628, places=0) == "63%"
    assert fmt.percentage(None) == fmt.MISSING


# ── frequencies ───────────────────────────────────────────────────────────────

def test_one_in_n_from_a_probability():
    assert fmt.one_in(0.018) == "1 in 56"
    assert fmt.one_in(0) == fmt.MISSING
    assert fmt.one_in(None) == fmt.MISSING


# ── ratio ─────────────────────────────────────────────────────────────────────

def test_ratio_does_not_render_as_a_return_on_investment():
    """`26.7:1` reads as money returned. The quantity includes monetized health,
    so it is not."""
    assert fmt.ratio(8000, 300) == "26.7x"
    assert ":" not in fmt.ratio(8000, 300)
    assert fmt.ratio(1, 0) == fmt.MISSING


# ── everything is presentation only ───────────────────────────────────────────

def test_formatting_never_mutates_its_input():
    from decimal import Decimal
    v = Decimal("0.0319")
    before = str(v)
    fmt.qaly(v)
    fmt.money(v)
    assert str(v) == before


def test_every_formatter_survives_none():
    for f in (fmt.money, fmt.signed_money, fmt.qaly, fmt.healthy_days,
              fmt.percentage, fmt.probability, fmt.one_in, fmt.months,
              fmt.years, fmt.count):
        assert f(None) == fmt.MISSING
