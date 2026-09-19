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
