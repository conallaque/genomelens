# Health-economics report architecture

Working document for the findings-first report redesign. Phase 1 (map the
existing implementation) is complete; this records what was found before any
code changed, so the "before" state is auditable.

Baseline artefact: `scripts/make_econ_sample.py` on the synthetic whole-genome
fixture. No human genome is involved at any step.

---

## 1. Where economics is calculated

| Module | Lines | Owns |
|---|---:|---|
| `econ/health_economics.py` | 3,208 | per-finding records, category→economics tables, cohort scaling, personal sheet |
| `econ/value_of_information.py` | 1,939 | VOI, PSA, tornado, EVPI/EVPPI, health capital, real options, risk preferences, privacy welfare |
| `econ/params.py` | 1,246 | parameter provenance registry (65 params, 16 assumptions, 75.4% sourced) |
| `econ/engine.py` | 1,210 | `Finding`/`ConditionPool`, pooling, adherence, life table, 3-state cohort, `run_psa` |
| `econ/decision.py` | 676 | decision layer, competing-mortality framing |
| `econ/frontier.py` | 582 | no-test / chip / WGS strategy frontier |
| `econ/plain.py` | 398 | plain-language translation of results |
| `econ/markov.py` | 323 | illustrative Markov cohort CEA + budget impact |

## 2. Where economics is formatted

| Module | Lines | Notes |
|---|---:|---|
| `report/renderers.py` | 7,403 | the great majority of econ HTML |
| `econ/health_economics.py` | (4 fns, ~47 markup lines) | **the only calculation module that emits HTML** |
| `pipeline.py` `_build_consolidated_econ_page` | ~150 | **regex-scrapes `<section>` blocks out of already-generated HTML** |

## 3. Are calculation and presentation coupled?

Mostly **no**, and this is the single most important Phase 1 finding — it
changes the plan.

Measured, not assumed:

- Recomputation of core totals inside `report/renderers.py`: **5 occurrences**,
  all of them counts or trivial sums (`sum(1 for ...)`, tornado swing total).
  No renderer independently derives NMB, ICER, incremental cost, or QALYs.
- HTML emission inside `econ/`: confined to `health_economics.py`.

So the spec's mandate — *one calculation layer → one normalised model → many
views* — is already ~80% true. A full `EconomicsReportResult` dataclass
retro-fitted over 9,606 lines of engine code would be a large typing exercise
that buys little against the stated goal.

**Two real violations remain, and they are where the risk actually lives:**

1. **`_build_consolidated_econ_page` parses HTML to build HTML.** It regex-lifts
   `<section id="value-of-information">` and `<section id="health-economics">`
   out of `report.html` and splices them beside the personal sheet's `<body>`.
   The lifting was a deliberate choice (one computation, two presentations) but
   the mechanism is string-scraping: any change to section ids, nesting, or
   `<style>` placement silently produces a partial page.
2. **Five separate money formatters.** `report/renderers.py` defines `money()`
   three times and `_money()` once; `econ/health_economics.py` defines `_money()`.
   Same quantity, five rounding paths.

## 4. Canonical model — the chosen shape

Rather than re-plumbing the engine, introduce a **normalising view at the report
boundary**:

```
build_report_payload(economics_result, voi_result, personal_econ, ...) -> dict
```

One function, assembled from the payloads the engine already returns, serialised
to JSON alongside the PDF. This yields the canonical-model benefit (every
displayed total has one source, and that source is inspectable) and the JSON
artefact, without changing a single economic calculation.

Section renderers then consume the normalised payload instead of the raw dicts,
and `_build_consolidated_econ_page` calls them directly instead of scraping.

## 5. Presentation bugs vs calculation bugs

Found on the synthetic fixture during Phase 1. **Nothing here has been changed
yet** — documented first, per the brief.

### Presentation (numbers are right; the report lets them collide)

- **P1 — two "Net monetary benefit" values, unqualified.** The deterministic
  reference case reports **$7,410**; the PSA reports **$7,759** (95% interval
  $1,435–$20,594). Both are labelled "Net monetary benefit". These are
  legitimately different quantities (point estimate vs mean of the simulated
  distribution) and the difference is small, but nothing on the page says so.
- **P2 — the Markov model sits beside the reference case at similar visual
  weight.** It reports **ΔQALY +0.436** and **NMB $46,301** against the
  reference case's **+0.032** and **$7,410** — a 13.6× and 6.2× difference. The
  Markov model is an *illustrative structural re-estimate* with its own
  parameters, not this genome's reference case, and the section text says so in
  prose. A reader scanning figures will not get that. This is the most
  misreadable thing in the current output.
- **P3 — per-finding values and pooled values for the same condition differ
  without adjacency.** APOE ε4 shows **$3,869** cost averted / **$3,514** NMB on
  the individual sheet and **$1,630** / **$1,516** in the pooled condition table.
  The gap is the pooling correction working correctly; the two figures are ~40
  pages apart.
- **P4 — ΔQALY rendered as `0.00`** in the per-finding table for values that are
  non-zero at three decimals (APOE, MTHFR). Precision chosen per-table.

### Calculation — traced

- **C1 — WGS retrospective vs prospective. RESOLVED: no contamination in the
  calculation; the defect is in the label.**

  `marginal_chip_to_wgs` = **$3,475** is computed at
  `value_of_information.py:711` as

  ```python
  wgs_only_value = sum(r["nmb"] for r in nmb_rows if r["wgs_only"])
  ```

  — a sum over findings **this genome actually contains**, i.e.
  `analysis_basis = observed_findings`. It is passed into `decision.py` as
  `wgs_only_findings_value` and re-emitted there as `retrospective_value`, which
  is why the two figures agree.

  The prospective quantities come from a different function and different
  inputs: `number_needed_to_sequence` = 56, `value_per_finding` = $25,189,
  `gross_expected_value` = $491, `net_expected_value` = $291, all derived from
  `wgs_yield_acmg_secondary x (1 - chip_detection_share_monogenic)`. These are
  `analysis_basis = prospective_expected_yield`.

  **The two bases are computed independently and never combined.** The renderer,
  however, labels the observed figure *"added expected value the WGS-only
  findings unlock"* — prospective language on a retrospective quantity. Fixed by
  tagging, not by changing arithmetic: `TestingDecision` now carries
  `prospective_*` and `observed_*` fields with explicit `AnalysisBasis`, and
  `_check_wgs_separation` reconstructs the observed sum from the findings so a
  future change that sources it from the prior fails a test.

- **C2 — MTHFR negative NMB. RESOLVED: correct, no change.**

  The identity at `value_of_information.py:651` is

  ```python
  nmb = dqaly * wtp + dcost - interv
  ```

  which is `lambda x dQALY - dCost` with `dCost = interv - dcost`. Sign
  convention is consistent with every other row. For MTHFR: `dcost` = $219,
  `dqaly` ~ 0.00065, `interv` = $500, giving `65 + 219 - 500 = -216`. The
  economics are unchanged; a test pins the identity.

  What C2 *did* surface is a display defect: `dqaly` of 0.00065 rendered as
  `0.00` at two decimals, indistinguishable from a finding modelled to
  contribute nothing. `report/format.py:qaly` now adds digits until the value is
  visible.

- **C3 — the same finding is priced by two independent code paths, and they
  disagree by up to 17x. NOT RESOLVED. Documented, not consolidated.**

  Both tables render in the same consolidated page:

  | finding | curated per-finding table | value-of-information | ratio |
  |---|---:|---:|---|
  | SLCO1B1 | $1,521 | $90 | 16.9x |
  | CYP2D6 | $2,500 | $210 | 11.9x |
  | CYP2C19 | $4,452 | $490 | 9.1x |
  | APOE e4 | $550 | $3,514 | 6.4x *(other direction)* |
  | PTPN22 | $5,626 | $2,509 | 2.2x |
  | MTHFR | $2,035 | -$216 | **sign flip** |

  `analyze_personal_economics` prices from curated per-category cost/value/QALY
  tables; `analyze_value_of_information._collect` builds its own finding list
  from `p_event x rrr x cost_of_illness x intervention`. The parameterisations
  were never reconciled.

  Per the brief — *"if the same quantity is calculated differently in multiple
  places, trace which is authoritative before consolidating"* — **nothing has
  been consolidated.** Which path is authoritative is an open question that
  needs a decision, not a merge.

  The two paths also **share no join key**: the curated table names an *action*
  ("Avoid clopidogrel non-response"), the VOI rows name a *finding* ("CYP2C19
  Intermediate Metabolizer"). Only PTPN22 matches by name, so cross-path
  agreement is not machine-checkable. The validator reports both the one
  detectable divergence and the structural fact that the rest cannot be checked.

  **Recommended fix (deferred):** emit a stable `finding_id` from both paths.
  Small, but it is a calculation-layer change and out of scope here.

### Not a bug — do not "fix"

- **Registry 75.4% sourced vs report 46.9% verified PMID/DOI.** Different
  denominators: registry-only (65 params) vs whole-model coverage (382). The
  test comment at `tests/unit/test_econ_params.py:88-92` states this explicitly.
  Both must stay, each labelled with its denominator. Collapsing them into one
  number would create a new inconsistency, not remove one.

## 6. Constraint carried into implementation

`tests/unit/test_econ_params.py:96` asserts `pct_sourced >= 75.0`. Live value is
**75.4%** (49 sourced / 65 total). **One** added unsourced parameter drops it to
74.2% and fails. This redesign must add no registry parameters. The spec's
"evidence foundation" weighted-provenance score would add some; the spec's own
escape clause (retain raw provenance categories, leave a design note) is taken.

## 7. Implementation order

1. Replace the regex scraper with direct payload-driven rendering.
2. Collapse five money formatters into one, diffing the bodies first.
3. Add the report consistency validator + tests (NMB identity, net-cash
   identity, dominance sign agreement, pooling/no-double-charge, WGS
   prospective-vs-observed separation, budget-impact peak PMPM).
4. Findings-first page order, pages 1–8.
5. Move advanced economics into the appendix.
6. Regenerate sample PDF + JSON payload; inspect page by page.

---

## 8. Methodology change record — canonicalisation (reporting change)

> The findings-first report now uses registry-backed **parametric expected NMB**
> as the canonical per-finding economic value. Legacy curated outcome values are
> retained for audit and calibration but are no longer presented as comparable
> NMB estimates.

This is a **reporting and canonicalisation change, not a change to the
economics.** No calculation was altered, no parameter was retuned, and the
pooled reference case is unchanged at **$7,410**.

### What the curated value actually is — traced, not assumed

For `CYP2C19`, `outcome_value` in the curated table is:

```
p_rx x p_adr x rrr x adr_cost  =  0.10 x 0.20 x 0.50 x 30,000  =  $300
```

That is already an **expected** cost averted, conditioned on exposure and
efficacy. `analyze_personal_economics` then multiplies it — and the QALY gain —
by `prevalence`:

```
avoided = 300 x 0.30 x 0.6 x 0.8626 x 0.5  = $23
qaly    = 0.35 x 0.30 x 0.5                = 0.0525   ->  $4,529
```

`prevalence` is the **population frequency of the genotype** — the probability
that a random person carries the variant. This person carries it; that is why
the finding is in the report. So two distinct defects sit in the curated path:

1. **Conditioning on a settled fact.** An already-expected value is multiplied
   by the probability of something already observed. Structurally the same error
   the sequencing section was fixed for: a forecast about a fact already
   settled.
2. **QALYs never conditioned on the event.** The parametric path conditions the
   QALY on `p_rx x p_adr x rrr` = 0.01. The curated path uses prevalence (0.30)
   instead — **30x** — so the QALY term is roughly 26x larger after the other
   factors, and it supplies about **98%** of the $4,452 curated figure.

The curated value is therefore neither a gross event value nor an expected NMB.
It is a **prevalence-weighted mixture containing a conditioning error**, and it
is named `legacy_curated_value` with basis
`curated_prevalence_weighted_mixture` for exactly that reason. It is retained
for audit and calibration; it is never rendered beside the canonical figure.

**Not claimed:** that every curated number is wrong in every context. The
per-category tables encode real literature. What is established is that the
quantity they produce after `analyze_personal_economics` is not an expected NMB
and must not be labelled as one.

### Why parametric is canonical

Decomposable into explicit probabilities, effect sizes, costs and utilities ·
participates in the PSA · parameters carry registry provenance · responds
correctly to willingness-to-pay, adherence and horizon changes · does not depend
on a second opaque valuation pathway.

### Reconciliation, synthetic whole-genome fixture

| Finding | Legacy curated | Canonical NMB | Diff | Ratio | Sign |
|---|---:|---:|---:|---|---|
| MTHFR folate metabolism | $2,035 | **-$216** | $2,251 | 9.4x | **opposite** |
| SLCO1B1 / simvastatin | $1,521 | **$90** | $1,431 | 16.9x | agree |
| CYP2D6 / codeine | $2,500 | **$210** | $2,290 | 11.9x | agree |
| CYP2C19 / clopidogrel | $4,452 | **$490** | $3,962 | 9.1x | agree |
| APOE e4 | $550 | **$3,514** | $2,964 | 6.4x | agree |
| PTPN22 R620W | $5,800 | **$2,509** | $3,291 | 2.3x | agree |

Reasons, in order of contribution: QALY conditioned on genotype prevalence
rather than event probability (dominant) · `outcome_value` multiplied by
prevalence a second time · evidence haircut applied by the parametric path only
· real-world adherence applied by the curated path only · marginal-cost fraction
applied to every curated row but only to the condition branch of the parametric
path.

APOE runs the other way because the curated path prices it as a statin
intervention while the parametric path prices the Alzheimer cost-of-illness
anchor.

### Totals

| Quantity | Value | What it is |
|---|---:|---|
| Legacy curated total | $18,625 | pre-canonicalisation personal sheet |
| Canonical expected NMB, summed | $10,072 | standalone per-finding values, **not additive** |
| Pooled reference case | **$7,410** | authoritative; overlapping signals pooled first |

No pro-rata allocation is invented to make the middle row reconcile to the
bottom one.

## 9. Further defects found during canonicalisation

- **D1 — unrecognised confidence vocabulary.** The novel-variants module emits
  `"higher"` for its most confident predictions (AlphaMissense 0.9963).
  `_CONFIDENCE_CONC` in `value_of_information.py` has keys `high`/`moderate`/
  `low` only, so `_beta` silently falls back to the **moderate** concentration
  for the run's single most confident computational finding. Display is
  normalised in the payload; the PSA spread is **not** changed, and the
  validator reports it.
- **D2 — curated sheet double-counts PTPN22.** Two curated records describe one
  pathway ($5,626 carrier row + $174 symptom-awareness row) and both render.
- **D3 — marginal-cost fraction applied inconsistently inside the parametric
  path.** The condition branch of `_finding_nmb` scales avoided cost by
  `MARGINAL_COST_FRACTION`; the pharmacogenomics branch does not.
- **D4 — `ConditionPool` models cost and adherence at condition level.**
  Benefits pool correctly on the risk scale, but `intervention_cost()` charges
  `max(cost)` and `adherence()` asserts one value per pool. Correct where one
  action covers every finding (the fixture's only multi-finding pool); wrong as
  soon as one condition maps to genuinely separate interventions. Latent on this
  fixture, so no number changes yet.
- **D5 — `baseline_risk() = max(p_event)`.** The docstring says this "keeps the
  strongest evidence of elevated risk"; the highest estimate is not the
  strongest evidence. Needs an explicit source hierarchy, and at minimum each
  pool should record which risk source won and why.
