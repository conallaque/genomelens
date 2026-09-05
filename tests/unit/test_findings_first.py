"""The findings-first pages.

These assert the *contract* the pages have with the payload — that the dollar
column is canonical, that legacy values never surface, that a negative value
renders and renders correctly — rather than asserting the HTML looks a
particular way. Layout can change; those guarantees cannot.
"""
from __future__ import annotations

import re

import pytest

from report.findings_first import (
    MAIN_REPORT_PAGES,
    render_findings_first,
    render_page_two,
)
from report.payload import (
    EconomicsReportPayload,
    FindingEconomics,
    PricingPath,
    ReferenceCase,
    ReportMetadata,
    Uncertainty,
)


def _fixture_like() -> EconomicsReportPayload:
    """The synthetic whole-genome fixture's canonical values."""
    def f(name, pid, nmb, conf, legacy=None, cat="", action=""):
        return FindingEconomics(
            display_name=name, economic_pathway_id=pid,
            pricing_path=PricingPath.VOI_PARAMETRIC,
            canonical_expected_nmb=nmb,
            economic_value_basis=("parametric_expected_nmb" if nmb is not None
                                  else "not_yet_standardised"),
            legacy_curated_value=legacy,
            evidence_confidence=conf, category=cat, action_summary=action)

    return EconomicsReportPayload(
        metadata=ReportMetadata(is_synthetic=True, willingness_to_pay=100_000.0,
                                input_label="synthetic_wgs.vcf"),
        reference_case=ReferenceCase(
            incremental_cost=-4223.0, incremental_qalys=0.0319, nmb=7410.0,
            wtp=100_000.0, dominance_status="dominant (more health, lower cost)",
            icer_note="not defined — strategy dominates"),
        uncertainty=Uncertainty(
            psa_available=True, psa_iterations=1500, psa_mean_nmb=7759.0,
            nmb_ci_low=1435.0, nmb_ci_high=20594.0,
            probability_cost_effective=0.9987, probability_cost_saving=0.956),
        findings=[
            f("CYP2C19 Intermediate Metabolizer (IM) (clopidogrel)",
              "pgx:cyp2c19:clopidogrel:mace", 490.0, "high", 4452.0,
              "Pharmacogenomic / genomic", "Avoid clopidogrel non-response"),
            f("CYP2D6 Intermediate Metabolizer (IM) (codeine)",
              "pgx:cyp2d6:codeine:opioid_toxicity", 210.0, "high", 2500.0,
              "Pharmacogenomic / genomic", "Avoid opioid toxicity"),
            f("SLCO1B1 Intermediate Function (IM-like) (simvastatin)",
              "pgx:slco1b1:simvastatin:myopathy", 90.0, "high", 1521.0,
              "Pharmacogenomic / genomic", "Prevent statin myopathy"),
            f("APOE ε4 carrier", "legacy:apoe:apoe_4", 3514.0, "moderate", 550.0,
              "Risk", "Intensive cardiovascular risk reduction"),
            f("PTPN22 R620W carrier", "legacy:ptpn22", 2509.0, "moderate", 5800.0,
              "Carrier Screening", "Awareness for autoimmune symptoms"),
            f("Folate Metabolism (MTHFR)", "legacy:mthfr", -216.0, "low", 2035.0,
              "Wellness", "Targeted supplementation"),
            f("1 actionable wellness variant(s)", "unkeyed:wellness", None, "low",
              1767.0, "Wellness", "Nutrient metabolism"),
        ])


@pytest.fixture
def html_out() -> str:
    return render_findings_first(_fixture_like())


def _text(html_str: str) -> str:
    """Visible text, with tags collapsed to nothing rather than to a space.

    A label split by <br> — "Not yet<br>standardised" — must still read as one
    phrase, or every assertion about wording depends on where the line happens
    to break. The break becomes a space, not nothing: deleting it welds the two
    words together and the assertion fails for a reason that has nothing to do
    with the report.
    """
    import html as _h
    stripped = re.sub(r"<style.*?</style>", " ", html_str, flags=re.S)
    stripped = re.sub(r"<br\s*/?>", " ", stripped)
    return re.sub(r"\s+", " ", _h.unescape(re.sub(r"<[^>]+>", " ", stripped)))


# ── the dollar column is canonical ────────────────────────────────────────────

@pytest.mark.parametrize("value", ["$490", "$210", "$90", "$3,514", "$2,509"])
def test_canonical_values_render(html_out, value):
    assert value in _text(html_out)


@pytest.mark.parametrize("legacy", ["4,452", "2,500", "1,521", "2,035", "18,625"])
def test_legacy_curated_values_never_appear(html_out, legacy):
    """They live in the payload audit fields and the reconciliation output, not
    beside the canonical figure."""
    assert legacy not in _text(html_out)


def test_canonical_per_finding_sum_is_not_shown_as_a_total(html_out):
    """$10,072 is a diagnostic, not a user-facing total: standalone expected
    values are not additive."""
    total = sum(f.canonical_expected_nmb or 0
                for f in _fixture_like().findings)
    assert f"{total:,.0f}" not in _text(html_out)


# ── the negative value, which is a feature ────────────────────────────────────

def test_mthfr_renders_negative_with_the_sign_in_front(html_out):
    txt = _text(html_out)
    assert "-$216" in txt
    assert "$-216" not in txt


def test_a_negative_finding_is_not_dropped(html_out):
    assert "Folate Metabolism (MTHFR)" in _text(html_out)


# ── uncertainty ───────────────────────────────────────────────────────────────

def test_uncertainty_reports_a_count_not_a_bare_hundred_percent(html_out):
    """The PROBABILITY must be a count, not a rounded 100%.

    Originally asserted that the string "100%" appeared nowhere in the report,
    which was a proxy for the real rule and too blunt to survive: the page-5
    pooling note now says "any value above 100% is the artifact this correction
    removes", where 100% is the name of a boundary, not a reported probability.
    Suppressing that sentence to satisfy a string match would have removed an
    explanation to protect a proxy for the thing it explains.

    Pinned on the uncertainty block specifically, which is where a rounded
    probability would actually mislead.
    """
    txt = _text(html_out)
    assert "1,498 of 1,500" in txt
    assert ">99%" in txt
    # No bare "100%" presented AS the probability of cost-effectiveness.
    for phrase in ("cost-effective (100%)", "cost-effective in 100%",
                   "100% of simulations", "cost-saving (100%)"):
        assert phrase not in txt, f"probability rendered as a bare {phrase!r}"


# ── reference case ────────────────────────────────────────────────────────────

def test_reference_case_is_labelled_and_unambiguous(html_out):
    txt = _text(html_out)
    assert "7,410" in txt
    assert "Reference-case NMB" in txt


def test_not_additive_statement_is_present(html_out):
    assert "not additive" in _text(html_out)


# ── ordering ──────────────────────────────────────────────────────────────────

def test_groups_appear_in_clinical_order(html_out):
    txt = _text(html_out)
    order = [txt.index(g) for g in ("MEDICATION & PRESCRIBING",
                                    "RISK & PREVENTION",
                                    "LOWER-CONFIDENCE & EXPLORATORY")
             if g in txt]
    assert order == sorted(order)


def test_high_confidence_pgx_precedes_a_larger_moderate_finding(html_out):
    """$90 SLCO1B1 (high) must appear above $3,514 APOE (moderate): clinical
    grouping outranks dollar value."""
    txt = _text(html_out)
    assert txt.index("SLCO1B1") < txt.index("APOE")


def test_within_a_group_ordering_is_by_canonical_value_descending(html_out):
    txt = _text(html_out)
    assert txt.index("CYP2C19") < txt.index("CYP2D6") < txt.index("SLCO1B1")


def test_unstandardised_findings_render_without_an_invented_value(html_out):
    txt = _text(html_out)
    assert "1 actionable wellness variant(s)" in txt
    assert "Not yet standardised" in txt
    assert "no registry-backed" in txt.lower()


# ── labelling rules ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("banned", ["ROI", "Return on investment", "your genome is worth"])
def test_forbidden_value_language_is_absent(html_out, banned):
    assert banned.lower() not in _text(html_out).lower()


def test_synthetic_input_is_labelled_from_metadata(html_out):
    assert "synthetic input" in _text(html_out).lower()

    p = _fixture_like()
    p.metadata.is_synthetic = False
    assert "synthetic input" not in _text(render_findings_first(p)).lower()


def test_evidence_badge_carries_no_economic_colour():
    """Confidence and economic direction are separate dimensions. A HIGH badge
    must not be styled with the positive-value colour."""
    out = render_page_two(_fixture_like())
    for m in re.finditer(r'<span class="badge[^"]*">', out):
        assert "good" not in m.group(0)


def test_empty_payload_renders_without_crashing():
    out = render_findings_first(EconomicsReportPayload())
    assert "<html" in out


# ── the eight-page constraint ─────────────────────────────────────────────────

def test_report_declares_eight_pages():
    from report.findings_first import MAIN_REPORT_PAGES, page_count
    assert MAIN_REPORT_PAGES == 8
    # MAIN_REPORT_PAGES is the body count; the report also renders a
    # cover sheet carrying the synthetic-input and not-medical-advice
    # statements, so a minimal payload renders one more than the body.
    assert page_count(_fixture_like()) == MAIN_REPORT_PAGES + 1


def test_every_page_renders_and_is_numbered_over_eight():
    """A section that silently disappears would still produce valid HTML, so
    the numbering is checked rather than the section count.

    The total counts the cover as well. It is a numbered sheet like any other,
    and a footer reading "11 / 10" is exactly what happens when the two
    disagree.
    """
    out = render_findings_first(_fixture_like())
    total = MAIN_REPORT_PAGES + 1          # body sheets plus the cover
    for n in range(1, total + 1):
        assert f"{n} / {total}" in out, f"page {n} footer missing"
    assert out.count('<section class="sheet') == total
    # The cover is a sheet of the document, not a bolted-on A4 page.
    assert '<section class="sheet cover">' in out


def test_page_order_matches_the_agreed_structure():
    out = _text(render_findings_first(_fixture_like()))
    titles = [
        "Findings-first economics",
        "Your findings at a glance",
        "Medication-genotype findings",
        "Risk, prevention & exploratory findings",
        "How the findings combine",
        "Uncertainty, evidence & model discipline",
        "The testing decision",
        "Methods & provenance",
    ]
    positions = [out.index(x) for x in titles]
    assert positions == sorted(positions), "pages are out of order"


def test_page_two_omits_the_detailed_economic_basis():
    """Page 2 scans, pages 3-4 explain. Putting the basis line on every row is
    what pushed this page onto a second sheet."""
    from report.findings_first import render_page_two
    txt = _text(render_page_two(_fixture_like()))
    assert "medical cost averted" not in txt.lower()
    assert "cost to act" not in txt.lower()


def test_pgx_page_carries_the_detail_page_two_dropped():
    from report.findings_first import render_page_three
    txt = _text(render_page_three(_fixture_like()))
    assert "Medical cost averted" in txt
    assert "Cost to act" in txt
    assert "an instruction to start, stop or change a medication" in txt


def test_methods_page_reports_both_provenance_denominators():
    from report.findings_first import render_page_eight
    p = _fixture_like()
    p.provenance.registry_n_parameters = 65
    p.provenance.registry_pct_sourced = 75.4
    p.provenance.model_pct_resolvable = 46.9
    p.provenance.model_pct_unsourced = 4.2
    txt = _text(render_page_eight(p))
    assert "75.4%" in txt and "46.9%" in txt and "4.2%" in txt


def test_testing_page_never_calls_the_observed_figure_expected():
    from report.findings_first import render_page_seven
    p = _fixture_like()
    p.testing_decision.observed_wgs_only_findings = 2
    p.testing_decision.observed_wgs_only_value = 3475.0
    p.testing_decision.prospective_gross_expected_value = 491.0
    txt = _text(render_page_seven(p))
    assert "Observed standalone contribution" in txt
    i = txt.index("Observed — in this genome")
    assert "expected" not in txt[i:i + 260].lower().replace(
        "not an expectation", "")


def test_build_stamp_does_not_add_a_ninth_page():
    """REGRESSION. The pipeline's stamp chokepoint appended a visible <div>
    after the last page box, which paginated as a blank ninth page on an
    eight-page report. Paged documents now take the marker as a comment."""
    from core import build_stamp as bs
    marker = bs.build_stamp()["marker"]
    stamped = render_findings_first(_fixture_like()).replace(
        "</body>", f"<!--\n{marker}\n-->" + "</body>", 1)
    assert bs.marker_in(stamped) == marker
    assert stamped.count('class="sheet"') == 8


def test_the_glance_pages_never_drop_or_duplicate_a_finding():
    """The invariant render_page_two's docstring claims, at every size.

    The sheet split was rewritten to spread rows evenly rather than filling
    each sheet to its cap — a fixed budget left the last sheet holding five
    rows and 55% white space once the phantom action-rows were removed. An
    even split is easy to get subtly wrong at a boundary, and getting it wrong
    means silently losing a finding, so every size up to a few sheets is
    checked rather than one representative case.
    """
    from report.findings_first import render_page_two
    from report.payload import EconomicsReportPayload, FindingEconomics

    for n in list(range(0, 40)) + [47, 48, 49, 60, 97]:
        p = EconomicsReportPayload()
        p.findings = [FindingEconomics(display_name=f"F{i}",
                                       canonical_expected_nmb=1.0,
                                       evidence_confidence="high")
                      for i in range(n)]
        html = render_page_two(p)
        shown = sum(html.count(f">F{i}<") for i in range(n))
        assert shown == n, f"{n} findings in, {shown} rendered"


def test_a_continuation_header_reports_the_remainder_not_the_total():
    """A group spanning sheets repeated its full count on every sheet.

    "Risk & prevention, 12 findings" appeared on the sheet showing twelve and
    again on the sheet showing the last four, so the same twelve read as
    twenty-four to anyone adding up the headers.
    """
    import re

    from report.findings_first import render_page_two
    from report.payload import EconomicsReportPayload, FindingEconomics

    p = EconomicsReportPayload()
    p.findings = [FindingEconomics(display_name=f"F{i}",
                                   canonical_expected_nmb=1.0,
                                   evidence_confidence="high")
                  for i in range(60)]
    heads = re.findall(r'<div class="grouplabel">(.*?)</div>',
                       render_page_two(p))
    assert len(heads) > 2, "fixture must span at least three sheets"
    conts = [h for h in heads if "continued" in h]
    assert conts, "a 60-finding group must continue"
    for h in conts:
        assert "of 60 findings" in h, f"continuation repeats the total: {h}"
    # Each continuation must say how many were already shown, and those
    # counts must strictly increase down the document.
    earlier = [int(m) for h in conts
               for m in re.findall(r"(\d+) shown earlier", h)]
    assert earlier == sorted(earlier) and len(set(earlier)) == len(earlier), (
        f"'shown earlier' counts must increase: {earlier}")


def test_detail_card_fields_are_bounded():
    """The medication sheet must not grow with its inputs.

    That sheet bounds how many cards it shows but bounded nothing about their
    text, and it renders 24px clear of the footer — so ~60 extra characters on
    each of three PGx action strings pushed it past, onto the next page's
    masthead. The card title was clipped at 62 characters; the field values
    were not clipped at all.
    """
    from report.findings_first import _DETAIL_FIELD_CHARS, _detail_card
    from report.payload import FindingEconomics

    long = "extremely verbose clinical qualifier " * 10
    html = _detail_card(FindingEconomics(display_name="G", evidence_confidence="high"),
                        fields=[("Potential decision", long)], source="s")
    assert long not in html, "an unbounded field value reached the card"
    # The clip keeps it within the budget, plus the ellipsis the clipper adds.
    body = html.split('class="fv">')[1].split("</div>")[0]
    assert len(body) <= _DETAIL_FIELD_CHARS + 2, len(body)
