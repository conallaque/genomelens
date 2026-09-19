# GenomeLens

**From whole-genome evidence to decision-relevant health economics.**

GenomeLens is an experimental genomics and decision-analysis platform. It asks not only what a genome contains, but what supported findings could change, what evidence stands behind that change, what the modeled health and economic consequences are, how uncertain those consequences remain — and when the system should refuse to make a claim at all.

```
WHOLE-GENOME INTERPRETATION → CLINICAL ACTIONABILITY → EVIDENCE PROVENANCE
→ HEALTH ECONOMICS → UNCERTAINTY / VALUE OF INFORMATION → AUDITABLE REPORTING
```

## Real whole-genome demonstrations

GenomeLens was run end to end on three publicly released Genome in a Bottle whole-genome callsets. The reports below are sanitized public outputs from those runs — real results, not mock-ups.

| Sample | GIAB concordance | Reference-case NMB | Canonical expected NMB | Report |
|---|---|---|---|---|
| **HG002** | **169 / 169** · 0 mismatches · 0 normalization failures | **$957** | **$1,239** | [View report](artifacts/public-benchmarks/HG002/report-public.html) · [summary](artifacts/public-benchmarks/HG002/summary.json) |
| **HG003** | **187 / 187** · 0 mismatches · 0 normalization failures | **$2,295** | **$4,127** | [View report](artifacts/public-benchmarks/HG003/report-public.html) · [summary](artifacts/public-benchmarks/HG003/summary.json) |
| **HG004** | **170 / 170** · 0 mismatches · 0 normalization failures | **$937** | **$1,304** | [View report](artifacts/public-benchmarks/HG004/report-public.html) · [summary](artifacts/public-benchmarks/HG004/summary.json) |

```
526 / 526  scoped explicit GIAB truth calls concordant
  0        mismatches
  0        normalization failures
```

Scope: GIAB **v4.2.1** · **GRCh38** · **chr1–22** · a fixed curated validation panel. This is agreement on the explicitly benchmarked calls inside that scope — not whole-genome accuracy, not clinical validation, and not an endorsement by NIST or GIAB. The denominator is published in [`docs/VALIDATION.md`](docs/VALIDATION.md#per-sample).

<p align="center">
  <a href="artifacts/public-benchmarks/HG002/report-public.html"><img src="artifacts/public-benchmarks/HG002/preview.png" alt="GenomeLens public report, HG002" width="32%"></a>
  <a href="artifacts/public-benchmarks/HG003/report-public.html"><img src="artifacts/public-benchmarks/HG003/preview.png" alt="GenomeLens public report, HG003" width="32%"></a>
  <a href="artifacts/public-benchmarks/HG004/report-public.html"><img src="artifacts/public-benchmarks/HG004/preview.png" alt="GenomeLens public report, HG004" width="32%"></a>
</p>

These public reports are built from an explicit allowlist, not by redacting a production report. The production GenomeLens report is a separate internal artifact and is not published.

## It runs, on real public genomes

GenomeLens is exercised end-to-end against the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) Ashkenazim trio — HG002, HG003, HG004 — using the publicly released GIAB GRCh38 benchmark callsets and defined confident regions. These are real public reference genomes with an external truth set, not mockups or simulations.

```
526 / 526  explicit benchmarked calls concordant
  0        mismatches
  0        normalization failures
```

Per sample: [HG002](https://www.nist.gov/programs-projects/genome-bottle) [169](docs/VALIDATION.md#hg002), [HG003](https://www.nist.gov/programs-projects/genome-bottle) [187](docs/VALIDATION.md#hg003), [HG004](https://www.nist.gov/programs-projects/genome-bottle) [170](docs/VALIDATION.md#hg004). The three are a father–mother–son trio, so the calls are not statistically independent.

**What that is not.** Not 100% whole-genome accuracy. Not clinical validation. Not an endorsement by NIST or GIAB. It is agreement on the explicitly benchmarked calls inside the defined truth and callability scope — a denominator we publish rather than hide.

## What a report contains

Each run produces a genomic and health-economic report containing: supported findings, carrier status, pharmacogenomics, callability and source state; evidence provenance per finding; reference-case and standardized pathway economics; incremental cost, incremental QALYs, ICER or dominance; probabilistic sensitivity analysis; explicit refusals and limitations.

Reports are published only through the allowlist export profile in [`PUBLIC-REPORT-SPEC.md`](docs/PUBLIC-REPORT-SPEC.md). A sanitized trio report is [not published in this release](docs/RELEASE-STATUS.md#trio-reports).

## It does not price a genome

GenomeLens deliberately produces no single "your genome is worth $X" number. That figure cannot be constructed honestly: it requires adding quantities that answer different questions, for different beneficiaries, on different evidence.

Seven estimands are kept distinct — reference-case · canonical standardized pathway · reproductive · cascade/family · testing-replacement · payer/budget-impact · value-of-information — and are not summed unless an explicit additive-attribution contract permits it.

Action value and information value can represent closely related economic consequences viewed from different decision perspectives, so GenomeLens does not automatically add them.

Headline trio economics are regenerated against each engine build and published only with the build identifier that produced them. None is [published in this release](docs/RELEASE-STATUS.md#trio-economics).

## Evidence, assumption, and refusal are different states

```
EVIDENCE-DERIVED · MODEL / SCENARIO ASSUMPTION · CONDITIONAL
UNRESOLVED · REFUSED · NOT ECONOMICALLY APPLICABLE
```

Findings and economic arms are classified on two separate axes; both vocabularies are defined in [`GLOSSARY.md`](docs/GLOSSARY.md).

A missing value is never silently converted to $0. A pathway that could not be evaluated and one evaluated at zero are different facts.

A refusal is not a negative biological result. Declining to monetize a finding says nothing about whether the finding is real.

An explicit economic-validation matrix requires every tracked pathway to terminate in a defined disposition rather than falling silently through the model. The claim is "every tracked row reaches an explicit terminal disposition" — coverage and accountability, not a claim that every row is empirically validated.

The matrix counts themselves are regenerated per build and are [not published in this release](docs/RELEASE-STATUS.md#matrix-counts).

## Uncertainty and value of information

GenomeLens carries probabilistic sensitivity analysis, decision-uncertainty summaries, and value-of-information methods for asking whether resolving additional clinical uncertainty could change a modeled decision — and whether that information is worth its acquisition cost.

In the current lipid value-of-information pathway, holding clinical state fixed and changing the genome does not change the result. This is a property of the current model, not a claim that genomes are irrelevant to lipid decisions — and because it derives from the same model as the withheld figures below, it is reported as a current model property rather than as an established finding.

Numerical value-of-information examples are [not published in this release](docs/RELEASE-STATUS.md#methods-hold), pending review of one treatment-effect parameter.

## Where GenomeLens stops

Will not monetize a pharmacogenomic finding merely because a genotype is guideline-actionable when the medication or indication is unknown.

Distinguishes an unsupported absence claim from a true negative.

Does not claim production support for every copy-number, structural, repeat-expansion, CYP2D6 star-allele or mitochondrial-heteroplasmy problem.

Does not read a GIAB benchmark as clinical validation. Concordance establishes agreement with the external truth calls within the benchmarked scope; it says nothing about clinical utility.

Does not treat economic arms as additive by default.

Where a benefit is estimable but its real-world cost is not adequately sourced, GenomeLens can report a maximum model-consistent cost bound rather than silently treating the unknown cost as zero or inventing an actual market price.

## Public repository boundary

This repository contains selected methodology, validation evidence, curated reports and integration examples. The production engine, internal economic parameter registries, routing logic and qualification machinery are maintained separately.

Public artifacts are generated through a deliberately constrained export layer so that validation evidence and interpretable results can be inspected without exposing the production engine.

**Research software. Not a diagnostic device. Not medical advice.**
