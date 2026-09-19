# Economic methods

Standard health-economic definitions, stated so published results are interpretable. How each genomic finding is routed through the production engine is not published.

## Core quantities

**Incremental cost** — difference in expected cost between a strategy and its comparator.

**Incremental QALYs** — difference in quality-adjusted life years. A QALY weights survival by health-related quality of life.

**ICER** — incremental cost-effectiveness ratio, `ΔCost / ΔQALY`. Genuinely undefined only when `ΔQALY = 0`. Under dominance the ratio is finite but uninterpretable, so it is not reported. A positive ICER is also quadrant-ambiguous: more costly and more effective, and less costly and less effective, produce numerically indistinguishable ratios, so the signs of `ΔCost` and `ΔQALY` are reported alongside it — otherwise a favorable-looking ratio can mean health was traded away for savings.

**Extended dominance** — an option ruled out because a mixture of two others achieves more at the same cost. It is assessed before any option is called cost-effective; comparing options pairwise without it can favor an option that is never efficient.

**Dominance** — a strategy that is both less costly and more effective. Reported as a **status**, never as a ratio and never as `$0/QALY`.

**Net monetary benefit** — `NMB = (ΔQALY × λ) − ΔCost`, where λ is the willingness-to-pay threshold. NMB is reported with λ stated, because it is meaningless without it.

## Uncertainty

**Probabilistic sensitivity analysis** propagates parameter uncertainty through the model rather than reporting a point estimate alone. Parameters are drawn from registered distributions; results are summarized as distributions. Where joint structure between parameters is not represented, draws are independent, which understates decision uncertainty rather than overstating it.

**Decision uncertainty** asks a different question from parameter uncertainty: not "how precise is this number" but "how often does the preferred decision change". A wide interval that never changes the decision and a narrow one that straddles the threshold have very different consequences.

**EVPI** — expected value of perfect information: the ceiling on what resolving *all* represented uncertainty could be worth. Reported per decision for one individual, never scaled to a population; population and individual EVPI differ by orders of magnitude and answer different questions. It bounds any further research over the uncertainty actually modeled.

**EVSI** — expected value of sample information: what one *specific* additional measurement is expected to be worth.

**Net information value** — EVSI minus what the measurement costs to acquire. Treatment cost belongs inside the action's net benefit and is never subtracted a second time here.

The relations `0 ≤ EVSI ≤ EVPI` and `net = EVSI − cost` are enforced, not assumed.

Numerical value-of-information examples are [not published in this release](RELEASE-STATUS.md#methods-hold), pending review of one treatment-effect parameter.

## Estimands kept apart

**Reference-case** results use a stated reference configuration. **Canonical / standardized pathway** results route through one aggregation authority so a single health consequence is credited once, however many genomic routes reach it. These answer different questions and are never summed.

**Reproductive** value accrues to a prospective child and is never folded into the tested adult's total. **Cascade / family** effects reach relatives and are reported separately. **Testing-replacement** value concerns assays the analysis substitutes for. **Payer / budget impact** is a different perspective with a different time horizon.

Action value and information value can represent closely related economic consequences viewed from different decision perspectives, so they are not automatically added.

## Evidence basis

Every monetary contribution is classified — evidence-derived, model/scenario assumption, conditional, unresolved, refused, or not economically applicable — and the classification is reported with the number. Parameter provenance distinguishes published values, values derived by a stated step, and declared assumptions, and the model reports how much of its output rests on each.

## Bounds instead of invented prices

Where an intervention's benefit is estimable but its real-world cost is not adequately sourced, GenomeLens can report the **maximum model-consistent cost** at which the modeled decision would still hold. This is a bound, not an estimate of actual cost, not a recommended or fair price, and it is never added to any total. The alternative — treating an unknown cost as zero, or inventing a market price — would be worse than the bound.
