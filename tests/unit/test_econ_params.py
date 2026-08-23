"""Provenance guarantees for the economic parameter registry.

The registry exists because an internal audit found roughly two thirds of the
economic parameters had no traceable source. These tests are what stop that
state from returning: a new constant cannot enter the model without either a
citation or an explicit statement that it is a judgement call.
"""

import re

import pytest

import econ_params as ep


# ── Registry integrity ────────────────────────────────────────────────────

def test_registry_validates_clean():
    assert ep.validate_registry() == []


def test_every_parameter_has_a_recognised_tier():
    for p in ep.PARAMS.values():
        assert p.tier in ep.TIERS, f"{p.key} has tier {p.tier!r}"


def test_sourced_parameters_carry_a_resolvable_citation():
    # A citation has to be something a reader can actually look up. PMID, DOI,
    # ISBN and https URLs all resolve outside this repository; "internal
    # estimate" does not, and is exactly what the registry exists to prevent.
    pattern = re.compile(r"(PMID:\d+|doi:10\.\S+|ISBN:[\dX-]+|https?://\S+)", re.I)
    for p in ep.PARAMS.values():
        if not p.sourced:
            continue
        assert pattern.search(p.citation), (
            f"{p.key}: citation {p.citation!r} has no PMID / DOI / ISBN / URL")


def test_assumptions_are_declared_and_uncited():
    # An assumption that cites a source is mis-tiered — either the source
    # supports it (published/derived) or it does not (assumption, no citation).
    for p in ep.assumptions():
        assert p.note.strip(), f"{p.key}: assumption without a note"
        assert not p.citation.strip(), f"{p.key}: assumption must not cite"
        assert len(p.note) > 60, (
            f"{p.key}: an assumption's note must explain the judgement, "
            f"not just label it")


def test_derived_parameters_state_their_derivation():
    for p in ep.by_tier("derived"):
        assert len(p.note) > 40, (
            f"{p.key}: tier=derived must record the arithmetic step taken "
            f"from the source, so the derivation can be disputed separately "
            f"from the source")


def test_ranges_are_ordered_and_bracket_the_value():
    for p in ep.PARAMS.values():
        if p.low is None or p.high is None:
            continue
        assert p.low <= p.high, f"{p.key}: low > high"
        assert p.low <= p.value <= p.high, (
            f"{p.key}: value {p.value} outside its own stated range "
            f"[{p.low}, {p.high}]")


def test_probabilities_and_fractions_stay_in_unit_interval():
    for p in ep.PARAMS.values():
        u = p.units.lower()
        if "probability" in u or "fraction" in u or "relative risk" in u \
                or "multiplier" in u or "utility" in u:
            assert 0.0 <= p.value <= 1.0, (
                f"{p.key} is a {p.units} but has value {p.value}")


def test_no_duplicate_keys():
    assert len(ep.PARAMS) == len(set(ep.PARAMS))


# ── The headline provenance claim ─────────────────────────────────────────

def test_most_of_the_model_is_sourced():
    # The audit that prompted this module put ~65% of parameters in the
    # "invented" bucket. This test pins the improvement so it cannot silently
    # erode: a new unsourced constant drops the percentage and fails here.
    burden = ep.assumption_burden()
    assert burden["pct_sourced"] >= 85.0, (
        f"only {burden['pct_sourced']}% of parameters carry a citation; "
        f"the registry exists to keep this high")
    assert burden["n_assumption"] <= 6, (
        "declared assumptions are allowed but should stay few and "
        "individually justified")


def test_assumption_burden_arithmetic_is_consistent():
    b = ep.assumption_burden()
    assert b["n_published"] + b["n_derived"] + b["n_assumption"] == b["n_parameters"]
    assert abs(b["pct_sourced"] + b["pct_assumption"] - 100.0) < 0.15


def test_citation_list_is_deduplicated_and_attributed():
    refs = ep.citation_list()
    assert refs, "a model with sourced parameters must produce a reference list"
    assert len({r["citation"] for r in refs}) == len(refs), "duplicate citations"
    for r in refs:
        assert r["source"] and r["params"]


# ── Access discipline ─────────────────────────────────────────────────────

def test_unknown_key_fails_loudly():
    # Silently returning a default for an unregistered key would reintroduce
    # exactly the untraceable-constant problem the registry prevents.
    with pytest.raises(KeyError) as e:
        ep.get("cost_of_something_nobody_registered")
    assert "registered in econ_params" in str(e.value)


def test_value_reads_the_registry():
    assert ep.value("discount_rate") == ep.get("discount_rate").value


def test_cite_renders_assumptions_distinctly():
    assumption = ep.assumptions()[0]
    assert "Declared assumption" in assumption.cite()
    sourced = ep.by_tier("published")[0]
    assert "Declared assumption" not in sourced.cite()


# ── Specific values other modules depend on ───────────────────────────────

def test_key_method_parameters_match_their_cited_conventions():
    assert ep.value("wtp_per_qaly") == 100_000
    assert ep.value("discount_rate") == 0.03
    # Second Panel sensitivity range must be carried, not just the point value.
    assert ep.get("discount_rate").low == 0.0
    assert ep.get("discount_rate").high == 0.05


def test_effect_sizes_match_their_trials():
    # CTT primary prevention and DPP are the two effect sizes the model leans
    # on hardest; a typo here changes every cardiometabolic dollar figure.
    assert ep.value("statin_rrr_primary") == pytest.approx(0.27, abs=0.005)
    assert ep.value("dpp_rrr") == pytest.approx(0.58, abs=0.005)


def test_marginal_cost_fraction_only_ever_reduces_claimed_savings():
    assert 0.0 < ep.value("marginal_cost_fraction") < 1.0


def test_retired_longevity_parameter_has_not_returned():
    # $10,000 per longevity percentile was the single largest line in the
    # report and had no source. If something like it reappears, it must go
    # through the registry and justify itself.
    for key in ep.PARAMS:
        assert "per_percentile" not in key, (
            f"{key}: pricing a percentile of a composite score needs an "
            f"explicit anchor; the previous one had none")


# ── The coverage claim must not overstate itself ──────────────────────────

def test_coverage_is_reported_against_the_whole_model_not_just_the_registry():
    # A true statement about a subset, phrased as a statement about the whole,
    # is the same species of error the registry was built to fix. The burden
    # report must carry the honest denominator.
    b = ep.assumption_burden()
    assert "n_unregistered" in b and "pct_of_model_registered" in b
    assert b["n_unregistered"] > 0, (
        "the curated per-finding tables in health_economics.py hold hundreds "
        "of load-bearing numbers; reporting zero here means the counter broke")
    assert b["n_total_known"] == b["n_parameters"] + b["n_unregistered"]
    assert b["pct_of_model_registered"] < b["pct_sourced"], (
        "registry coverage of the whole model is necessarily lower than the "
        "cited share within the registry; if these are equal the denominator "
        "is wrong")


def test_unregistered_counter_sees_the_curated_tables():
    n = ep.count_unregistered_parameters()
    assert n > 100, f"expected the curated econ tables to dominate, got {n}"


def test_burden_scope_names_what_is_not_covered():
    assert "not yet registered" in ep.assumption_burden()["scope"]
