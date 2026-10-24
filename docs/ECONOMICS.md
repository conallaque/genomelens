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

**Probabilistic sensitivity analysis** propagates parameter uncertainty through the model rather than reporting a point estimate alone. Parameters are drawn from registered distributions; results are summarized as distributions. Where joint structure between parameters is not represented, draws are independent, and the direction of the error depends on the omitted correlations: the variance of net monetary benefit is λ²·Var(ΔQALY) + Var(ΔCost) − 2λ·Cov(ΔQALY, ΔCost), so independence understates it when the two are negatively correlated and overstates it when they are positively correlated.

**Decision uncertainty** asks a different question from parameter uncertainty: not "how precise is this number" but "how often does the preferred decision change". A wide interval that never changes the decision and a narrow one that straddles the threshold have very different consequences.

**EVPI** — expected value of perfect information: the ceiling on what resolving *all* represented uncertainty could be worth. Reported per decision for one individual, never scaled to a population; population and individual EVPI differ by orders of magnitude and answer different questions. It bounds any further research over the uncertainty actually modeled.

**EVSI** — expected value of sample information: what one *specific* additional measurement is expected to be worth.

**Net information value** — EVSI minus what the measurement costs to acquire. Treatment cost belongs inside the action's net benefit and is never subtracted a second time here.

The relations `0 ≤ EVSI ≤ EVPI` and `net = EVSI − cost` are enforced, not assumed.

Numerical value-of-information examples are [not published in this release](RELEASE-STATUS.md#methods-hold), pending review of one treatment-effect parameter.

## Choosing the model for a decision

GenomeLens values a decision, not a gene. Several models can describe one decision (a published analysis, a model rebuilt here from published components, an analysis of a related question); the report counts one of them and shows the others beside it, each with the reason it was not counted. The counted model is chosen by the validity of its evidence and its fit to the person, not by whether it is published: a published analysis of a *different* question is never the value (a comparison of screening schedules among carriers who are already known, for example, has no "not identified" arm, so it answers which schedule is best, not what the finding is worth, and is carried beside the decision as a **secondary published-strategy comparison** that is never the value and never in a total); a published analysis whose inputs are known to bias the decision ranks below a model that corrects them (a published value at the person's age whose penetrance is not corrected for ascertainment is shown but not counted); the counted model can carry biases of its own (at 45 the model built here compares with no surveillance, which overstates the increment and is printed as an assumption, and it credits colorectal cancer deaths averted only, with no morbidity or stage-specific treatment cost, which understates it; neither is sized, so the counted value is labeled hypothetical, not a supported estimate); and a source's value is never carried to another age as the value. Where the best-ranked model is one built here, its value is labeled **hypothetical and assumption-based**: a scenario, never a supported estimate. A published result is that source's result under its own population, comparator and horizon; where a source prints its result in another unit (per person screened, for example) the report restates it by stated arithmetic and says so, and reporting it is not a validation of the source.

Each headline value in the public demonstrations, and each modeled value in the private economic report and technical dossier, carries a **cost per QALY gained** calculated from its own change in cost and QALYs and checked against them, or states that the value dominates its comparator (more health at lower cost), in which case no ratio is printed. Intervals are **parameter uncertainty only**, and each printed interval is labeled so: the draws use the ranges a source prints or the registry holds (rescaled to a per-carrier value where the source prints another unit), and comparator uncertainty is not in them. Where no range is registered for an input, it is held at its point; an alternative model shown beside a counted value is not sampled and prints no interval; a distribution is never invented.

### How the age-30 value is derived

The published-model demonstration restates one published result per carrier. The inputs are the source's printed figures and GenomeLens's arithmetic on them; nothing else enters.

| Step | Quantity | Whose |
|---|---|---|
| 1 | Cost per person screened at age 30 with cascade testing, about $246 (2021 dollars); the scenario without cascade testing removes about $12 per person, leaving about $234 | The source's printed figures |
| 2 | ICER of the scenario without cascade testing, about $153,200 per QALY; carrier prevalence 0.3%; test $200 | The source's printed figures |
| 3 | Quality-adjusted life-year gain per person screened: $234 divided by $153,200, 0.00153 | GenomeLens (the source does not print it; it is backed out of the printed cost and ICER) |
| 4 | Remove the $200 test as sunk: about $34 per person screened | GenomeLens |
| 5 | Divide by the carrier prevalence: +0.5091 QALYs and about +$11,333 per carrier | GenomeLens |
| 6 | Restate 2021 dollars as 2019 dollars (medical-care price index 0.9489): +$10,754 | GenomeLens |
| 7 | Net monetary benefit: $100,000 x 0.5091 - $10,754 = $40,160; ICER $10,754 / 0.5091 = $21,121 per QALY | GenomeLens |

The arithmetic adds three assumptions. (a) The program's net effect other than the test accrues to carriers: non-carriers incur no downstream cost or benefit. (b) The four-gene pooled average applies to an MSH2 carrier: MSH2 risk is above the pooled average, so the average probably understates the carrier's benefit, by an unknown amount. (c) The printed figures are approximate ("about"), and the result inherits their rounding. The source did not report an MSH2-specific individual result and did not validate this transformation, so the value is a published result restated per carrier, not a result the source reported for an MSH2 carrier.

## Estimands kept apart

**Reference-case** results use a stated reference configuration. **Canonical / standardized pathway** results are built so that a single health consequence is credited once, however many genomic routes reach it. These answer different questions and are never summed.

**Reproductive** value accrues to a prospective child and is never folded into the tested adult's total. **Cascade / family** effects reach relatives and are reported separately. **Testing-replacement** value concerns assays the analysis substitutes for. **Payer / budget impact** is a different perspective with a different time horizon.

Action value and information value can represent closely related economic consequences viewed from different decision perspectives, so they are not automatically added.

## Evidence basis

Every monetary contribution is classified — evidence-derived, model/scenario assumption, conditional, unresolved, refused, or not economically applicable — and the classification is reported with the number. Parameter provenance distinguishes published values, values derived by a stated step, and declared assumptions, and the model reports how much of its output rests on each.

## Illustration

A short, fully synthetic worked example of these quantities, combining published clinical evidence with independently constructed economic inputs, is given in the [illustrative case study](../research/GenomeLens_TSC_Illustrative_HEOR_Case_Study.pdf). It is an educational calculation and not an output of the GenomeLens engine.

## Bounds instead of invented prices

Where an intervention's benefit is estimable but its real-world cost is not adequately sourced, GenomeLens can report the **maximum model-consistent cost** at which the modeled decision would still hold. This is a bound, not an estimate of actual cost, not a recommended or fair price, and it is never added to any total. The alternative — treating an unknown cost as zero, or inventing a market price — would be worse than the bound.
