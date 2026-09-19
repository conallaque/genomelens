# Validation

## Real-data benchmark

GenomeLens is run end-to-end against the Genome in a Bottle Ashkenazim trio — HG002, HG003, HG004 — using the publicly released GIAB GRCh38 benchmark callsets and their defined confident regions.

```
526 / 526  explicit benchmarked calls concordant
  0        mismatches
  0        normalization failures
```

Per sample: HG002 169, HG003 187, HG004 170.

## The denominator, stated

The figure counts positions where the truth set **states a genotype**. Positions the truth set declines to state are not counted as agreement — they are reported separately as a callability result. Counting them toward concordance would inflate it.

This matters more than the percentage. A concordance rate whose denominator is unpublished is not a measurement.

**What the benchmark is not:** not whole-genome accuracy; not clinical validation; not an endorsement by NIST or GIAB. Technical concordance says a genotype matches an external truth call within a defined scope. It says nothing about clinical utility.

## Out of benchmark scope

Regions outside the confident-region definition; variant classes the truth set does not adjudicate; and any capability requiring read-level evidence rather than a variant callset. These are stated rather than silently folded into a pass.

## Provenance

Published artifacts carry the engine build that produced them, the input callset and truth-set version, and the region definitions applied. A result that cannot be attributed to a specific build is not publishable as evidence.

## Synthetic stress testing — PLANNED, NOT YET RUN

A synthetic whole-genome **variant-call** cohort is designed to test behavior at scale: failure modes, output distributions, economic coverage and systems behavior across many inputs.

It is explicitly **planned and not yet executed**, and when it runs it will answer a different question from the GIAB benchmark. Synthetic inputs test stability and coherence at scale; they are not evidence of real-world accuracy, disease prevalence, or clinical validity, and they do not substitute for real-genome validation.

`OWNER REVIEW` — do not describe this cohort as existing until it has run.
