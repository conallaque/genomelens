# GenomeLens — technical overview

## What GenomeLens adds

Most genomic pipelines answer **"what variants are present?"**

GenomeLens addresses a different question: *what could these findings change, what evidence supports that change, what are the modeled health and economic consequences, how uncertain are they, and where should the model stop?*

It operates downstream of variant calling. It is not a variant caller, and it does not perform sequencing.

## Where it sits in a genomic workflow

```
variant callset (WGS)
     ↓
interpretation · actionability · evidence provenance
health economics · uncertainty · auditable report
```

The proof points are deliberately external where possible: a public GIAB trio (a historical result from an earlier engine build, not representative of v1.11 validation), an external truth set, a published benchmark denominator, and explicit refusal behavior. Sanitized reports for the current engine build are [not published](../docs/RELEASE-STATUS.md#trio-reports).

## For a health-economic reader

The economic layer is built around the distinctions that health technology assessment depends on:

- evidence-derived results separated from scenario assumptions
- conditional pathways that do not book value until the condition is met
- refusals that are reported rather than zero-filled
- perspective and estimand stated rather than implied
- budget impact kept apart from cost-effectiveness
- value of information treated as its own question

GenomeLens reports incremental cost, incremental QALYs, ICER or dominance, and net monetary benefit, with the relevant perspective, horizon, discounting and willingness-to-pay assumptions made explicit. Uncertainty is carried alongside the point estimate. Parameter provenance distinguishes published evidence, derived inputs and declared assumptions, so a report can show not only the result but what the result depends on.

GenomeLens does not claim payer validation and does not present modeled results as realized savings.

## What is public and what is not

The public repository exposes selected outputs, methods and validation behavior. The production implementation contains additional evidence qualification, modality handling and economic decision logic that is intentionally not part of the public distribution. The public material is designed to make the approach inspectable without reproducing the production engine.

## Status

Research software under active development. The current development version, v1.11, is private; the benchmark evidence published here describes the frozen v1 build. Not a diagnostic device, not clinically validated, not payer-approved. Benchmarked against public reference genomes within a stated scope.
