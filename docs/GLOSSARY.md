# Glossary

The normative vocabulary. Where another document uses one of these terms, this
is what it means.

## Two different axes

Findings and economic arms are classified on **separate axes** and the two lists
are not interchangeable. A single finding can produce no arm, one arm, or
several arms with different beneficiaries.

### Finding disposition — what the evidence supports about the finding

| Term | Meaning |
|---|---|
| `supported` | The evidence supports reporting this finding. |
| `conditional` | Reportable, but its consequence depends on a condition not known to hold. |
| `unresolved` | Could not be evaluated. Not a negative result. |
| `refused` | The model declines to make a claim. Not a negative biological result. |

In a report's page footer, **UNRESOLVED** is the label on a computed figure whose applicability to this donor is not established (the donor's age is stated, not observed). It is a statement about the donor, not about the model: the figure was evaluated.

### Economic disposition — what the model did with a pathway

| Term | Meaning |
|---|---|
| `modeled` | A monetary result was produced, with its evidence basis stated. |
| `conditional` | Value is not booked because the triggering condition is not known to hold. |
| `unresolved` | Could not be evaluated. Never rendered as `0`. |
| `refused` | The model declines to value it. |
| `not economically applicable` | No economic question arises. |

### Evidence basis — what a monetary contribution rests on

`evidence-derived` · `model/scenario assumption` · `conditional` · `unresolved` ·
`refused` · `not economically applicable`

Every monetary contribution carries exactly one of these, reported with the number.

## Quantities

| Term | Meaning |
|---|---|
| **QALY** | Quality-adjusted life year — survival weighted by health-related quality of life. |
| **Incremental cost / QALYs** | The difference against a stated comparator, not a total. |
| **ICER** | `ΔCost / ΔQALY`. See [`ECONOMICS.md`](ECONOMICS.md#core-quantities) for when it is undefined or ambiguous. |
| **Dominance** | Less costly *and* more effective. Reported as a status, never as a ratio. |
| **NMB** | `(ΔQALY × λ) − ΔCost`, meaningless without λ stated. |
| **λ (WTP)** | Willingness-to-pay threshold per QALY. |
| **PSA** | Probabilistic sensitivity analysis — parameter uncertainty propagated, not a point estimate. |
| **EVPI** | Ceiling on what resolving all represented uncertainty could be worth. |
| **EVSI** | What one specific additional measurement is expected to be worth. |
| **Estimand** | The precise quantity a number is an estimate of — including whose it is and against what comparator. Two numbers with different estimands are not addable. |
| **Callability** | Whether a position could be interrogated at all — distinct from what was found there. |
| **Model basis** | What kind of model stands behind a value: a published model's printed result (restated by stated arithmetic where the source prints another unit), or a hypothetical, assumption-based construction built here from published components or assumptions where the best-ranked model of the decision for the person is not a published one. |
| **Primary decision-relevant comparison** | The decision's own value: acting on a returned finding against not having it, as the best-ranked model of that decision gives it. |
| **Secondary published-strategy comparison** | A published result of a different question (for example, which surveillance schedule is best among carriers who are already known), shown beside a decision with each comparator named; never the decision's value and never added to a total. |
| **Evidence-qualified** | A value for which every model input is sourced and the value does not rest on an age that is stated rather than observed. |
| **Assumption-driven scenario** | A value that rests on at least one assumption, commonly the stated age; shown as a scenario. |
| **Supported calculation at the assumed age** | The decisions whose every input except the stated age is sourced; a stricter test than the modeled headline, never a total. A value above zero stays a scenario at the stated age; a supported zero made only of measured nulls is the same true zero as those findings, because a measured null does not depend on the age. |
| **Hypothetical** | A value whose best-ranked model is one built here, from components and assumptions, shown as a scenario. It is labeled, never counted as a supported estimate. |
| **Stated age** | The standardized age used for the public demonstrations, because public genomes carry no age: 45 for every genome, except one demonstration of a published-model decision, which uses the age its source reports (30). It is an assumption, labeled on every page. |
| **Callset** | The file of calls a genome's analysis starts from. A variants-only callset records only where a person differs from the reference, so silence in it is not evidence of a reference sequence. |
| **Validation (three kinds)** | *Reproduction*: a model rebuilt from a source reproduces that source's printed outputs. *Cross-check*: results set beside a model developed by others. *External validation*: results compared with observed outcomes in a population. No GenomeLens model has had the third. |
