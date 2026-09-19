# GenomeLens — technical overview for partners

## What GenomeLens adds

Most genomic pipelines answer **"what variants are present?"**

GenomeLens answers a different question: *what could these findings change, what evidence supports that change, what are the modeled health and economic consequences, how uncertain are they, and where should the model stop?*

It is a **decision and economic layer that sits downstream of sequencing** — not a variant caller, and not a competitor to one.

## For a sequencing company

```
your WGS output
     ↓
interpretation · actionability · evidence provenance
health economics · uncertainty · auditable report
     ↓
a decision/economic layer on top of data you already produce
```

The proof points are deliberately external where possible: a public GIAB trio, an external truth set, a published benchmark denominator, real generated reports, and explicit refusal behavior.

## For a payer or HTA reader

The economic layer is built around the distinctions payers care about and vendors usually blur:

- evidence-derived results separated from scenario assumptions
- conditional pathways that do not book value until the condition is met
- refusals that are reported rather than zero-filled
- perspective and estimand stated rather than implied
- budget-impact kept apart from cost-effectiveness
- value-of-information treated as its own question

GenomeLens does not claim payer validation and does not present modeled results as realized savings.

## For an HEOR reader

GenomeLens reports incremental cost, incremental QALYs, ICER/dominance and NMB with the relevant perspective, horizon, discounting and willingness-to-pay assumptions made explicit.

Uncertainty is carried alongside the point estimate. Parameter provenance distinguishes published evidence, derived inputs and declared assumptions, allowing the report to show not just the result but what the result depends on.

## What is public and what is not

The public repository exposes selected outputs, methods and validation behavior. Production GenomeLens contains additional evidence qualification, modality handling and economic decision logic that is intentionally not part of the public distribution.

The public material is designed to make the approach inspectable without reproducing the production engine.

## Status

Research software under active development. Not a diagnostic device, not clinically validated, not payer-approved. Benchmarked against public reference genomes within a stated scope. Engagement enquiries welcome.
