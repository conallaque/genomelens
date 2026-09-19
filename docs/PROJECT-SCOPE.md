# Project scope

What this project set out to do, how it was designed and tested, and where its
boundaries are. Stage semantics, methods and validation design are public.
Production routing rules, qualification logic and parameter registries are not.

## Problem definition

Genomic interpretation and economic decision analysis are usually separate
disciplines with separate tooling. Interpretation pipelines answer *what
variants are present*. Health-economic models answer *is this intervention
worth doing*, typically for a population and a decision defined in advance.

Very little sits between them. A supported genomic finding does not carry, with
itself, the evidence grade behind it, whether the relevant decision is actually
in play for this person, what the modeled consequence would be, how uncertain
that consequence is, or whether the honest answer is that no claim can be made.

GenomeLens explores how to connect those stages explicitly — finding, evidence,
actionability, applicability, economics, uncertainty, refusal — **without
reducing a genome to a single price.** The refusal to produce that single number
is a design position, not a missing feature: constructing it requires adding
quantities that answer different questions, for different beneficiaries, on
different grades of evidence.

## Research design

The public conceptual sequence is:

```
finding → evidence → actionability → applicability
        → economics → uncertainty → refusal / reporting
```

Each stage is a place where the analysis can stop. Stopping is a first-class
outcome: a pathway that reaches "actionable but not applicable" is reported as
such rather than being valued at zero or quietly dropped. The distinction
between *modeled*, *conditional*, *unresolved*, *refused* and *not
economically applicable* is carried end to end.

Stage semantics are public. The rules that decide which findings take which
route through the production engine are not.

## Quantitative and health-economic methods

The economic layer is built on standard cost-effectiveness machinery, stated
publicly so that published results are interpretable: incremental cost,
incremental QALYs, ICER and dominance-as-a-status, and net monetary benefit
against an explicit willingness-to-pay threshold. Perspective, time horizon and
discount rate are stated rather than implied.

Uncertainty is treated as a separate question from the point estimate.
Probabilistic sensitivity analysis propagates registered parameter
distributions; decision uncertainty asks how often the preferred decision
changes rather than how precise the number is. Value-of-information methods —
EVPI as the ceiling, EVSI for a specific measurement, and net information value
after acquisition cost — are carried with their identities (`0 ≤ EVSI ≤ EVPI`,
`net = EVSI − cost`) enforced rather than assumed.

Seven estimands are kept apart and are not summed unless an explicit additive
attribution contract permits it: reference-case, canonical standardized
pathway, reproductive, cascade/family, testing-replacement, payer/budget-impact
and value-of-information. Where a benefit is estimable but its real-world cost
is not adequately sourced, the model can report a maximum model-consistent cost
bound instead of inventing a price or treating the unknown as zero.

See [`ECONOMICS.md`](ECONOMICS.md).

## Genomic and bioinformatic work

Inputs are whole-genome variant callsets on a stated reference build. The
analysis is explicit about what that modality can and cannot answer: copy
number, structural variants, repeat expansions, some star-allele architectures
and mitochondrial heteroplasmy are not claimed as production-supported from a
SNV/indel callset alone, and anything requiring read-level evidence is labeled
rather than presented at uniform confidence.

Callability is treated as its own result. *Explicitly called*, *callable and
reference*, *explicit no-call* and *never interrogated* are four different
states, and collapsing them into "negative" is how false reassurance gets
manufactured. Variant normalization and truth-scope handling are part of the
benchmark rather than a preprocessing detail.

See [`ARCHITECTURE.md`](ARCHITECTURE.md) and [`LIMITATIONS.md`](LIMITATIONS.md).

## Validation and QA

External benchmarking uses the Genome in a Bottle Ashkenazim trio — HG002,
HG003, HG004 — against the publicly released GIAB GRCh38 callsets and their
defined confident regions: **526 / 526** explicit benchmarked calls concordant,
zero mismatches, zero normalization failures, across 169 / 187 / 170 calls per
sample.

The denominator is published deliberately. The figure counts positions where
the truth set *states a genotype*; positions it declines to state are reported
separately as callability rather than counted as agreement, which would inflate
the rate. A concordance figure whose denominator is unpublished is not a
measurement — and the scoped claim is more defensible than the unscoped one it
replaces.

Validation design also covers a failure taxonomy, independent checking of
economic identities, fail-closed behavior on unknown inputs, and explicit
refusal states. A synthetic stress cohort is designed but **not yet run**, and
is described as planned rather than existing.

See [`VALIDATION.md`](VALIDATION.md).

## Product and system design

GenomeLens is positioned downstream of sequencing rather than in competition
with it: a decision and economic layer over variant output that a sequencing
provider already produces. Public and production outputs are deliberately
different artifacts — the public report is built from an explicit allowlist, so
that a field added to production tomorrow cannot become public by default. A
blacklist fails open; an allowlist fails closed.

Evidence qualification is treated as a product principle rather than a
reporting afterthought: the report shows not only the result but what the
result depends on, and which parts of it rest on published evidence, a derived
step, or a declared assumption.

See [`PUBLIC-REPORT-SPEC.md`](PUBLIC-REPORT-SPEC.md) and [`OVERVIEW.md`](../partner/OVERVIEW.md).

## Research communication

Every reported result is designed to carry, alongside the number: the evidence
basis, the assumptions in force, the uncertainty around it, the limitations
that bound it, and — where applicable — the refusal state that explains why no
number is given. A refusal is not a negative biological result, and a missing
value is never silently rendered as `$0`.

## Creator role

The project was independently conceived and directed. That direction covers the
project concept and scope, the research questions, the analytical and evidence
requirements, the system architecture, the health-economic framework and
estimand separation, the validation strategy and benchmark design, the
failure-and-refusal semantics, the public/private export boundary, testing and
audit requirements, iterative methodological correction, and product direction.

Implementation was carried out with substantial AI-assisted software
engineering. The design decisions, evidence standards, validation requirements
and methodological corrections above are the project's own work, and are what
the public material is intended to evidence.

[`OWNER REVIEW`](RELEASE-STATUS.md#ai-assistance) — the wording of the AI-assisted development disclosure in the
preceding paragraph is pending owner approval before publication.

## What is public, and what is not

Public: methodology, validation design and results, architecture at the level
of stage semantics, economic definitions and estimand separation, limitations,
the public report contract, and a synthetic schema example.

Not public: production source, internal economic parameter registries, routing
and qualification logic, evidence mappings, and internal state machinery.

The public material is intended to make the approach inspectable and
criticizable without reproducing the engine that implements it.
