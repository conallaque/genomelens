# GenomeLens

## Health-economic decision analysis for genomic medicine

GenomeLens connects genomic evidence to the decisions it may change, and evaluates the potential health and economic consequences of those decisions while preserving uncertainty, missingness, and explicit refusal states.

```
GENOMIC EVIDENCE  →  DECISION / COMPARATOR  →  HEALTH CONSEQUENCE
                  →  ΔCOST + ΔQALY  →  ICER / NMB (where supported)  →  UNCERTAINTY / APPLICABILITY
```

**[View Benchmark v1](benchmark/v1/) · [Read the Public Evidence Brief](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Evidence_Brief.pdf) · [View the frozen results](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Results.pdf)**

<a href="benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Evidence_Brief.pdf">
  <img src="benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Cover.png" alt="GenomeLens Deterministic Benchmark v1 — Economic decision analysis for genomic medicine" width="520">
</a>

---

## Deterministic Benchmark v1

A frozen, public-data technical proof of concept: 30 primary public genomes and 6 reference/parent genomes, 36 in total, run end to end as a deterministic control (no AI interpretation) on pinned engine and tooling builds.

| | | | |
|---|---|---|---|
| **36 / 36**<br>public genomes completed | **7,723**<br>tests passed (26 skipped, 0 failed) | **9 / 9**<br>release gates | **0**<br>parity failures |
| **0**<br>arithmetic discrepancies | **301 / 301**<br>real-row assertions | **0 / 1,399**<br>zero-audit defects | **0 / 144**<br>forbidden-string hits across 144 rendered PDFs checked |

These figures matter for specific reasons. *Parity* means every dollar and percentage printed in a report traces to a structured result, so a report cannot say something the analysis did not. The *zero audit* examines every displayed zero and requires it to be a supported null or no-change state, never a missing value rendered as zero. *Arithmetic* recomputes each net-monetary-benefit value from its recorded components wherever those components are separable. *Real-row assertions* are checked against rows the pipeline actually produced, not fixtures.

Full record: [Public Evidence Brief](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Evidence_Brief.pdf) (narrative, 6 pages) and [Public Results](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Results.pdf) (compact technical record, 4 pages).

---

## Two analysis modes

| Individual genomic economics | Reproductive & carrier economics |
|---|---|
| **Unit of analysis:** one genome | **Unit of analysis:** a biologically paired parent pair |
| Genomic finding → evidence qualification → decision / comparator → health consequence → incremental cost and QALYs → ICER / NMB where supportable → uncertainty, conditionality or refusal | Parent A + Parent B → relevant findings → inheritance compatibility → offspring genotype / phenotype probability → reproductive decision → expected health and economic consequence |
| Decision contexts: screening, surveillance, treatment selection, pharmacogenomics, prevention, context-dependent future decisions | Decision pathways: resolving technical uncertainty, targeted testing, prenatal testing, PGT-M, or no change where evidence does not support intervention |

Pair-level reproductive value is a **different estimand**. It is not two individual reports added together and is never silently added to either parent's personal-health value.

---

## Decision economics, not genome monetization

- **Genotype alone is not economic value.** Value attaches to a decision, so a finding with no decision it could change has nothing to monetize.
- **Actionability needs a decision context**, and an economic estimate needs a defensible evidence chain from finding to outcome.
- **Unsupported pathways can be refused.** A refusal is a reported result with its reason, not a silent omission.
- **Missing is not zero.** Absent donor context, unresolved technical calls, missing intervention evidence, or economic terms that cannot be separated do not become "$0".
- **Distinct estimands stay distinct.** Present expected value, value if and when a decision arises, reproductive value, family or cascade value, conditional scenarios, alternative models and sensitivity analyses are reported separately and never silently combined.
- **Uncertainty stays visible** rather than being collapsed into a point estimate.

---

## What v1 establishes

- Reproducible execution across the intended public cohort, from pinned engine and tooling builds.
- Evidence and economic routing consistency, with no finding dropped from the economics without a recorded reason.
- Arithmetic integrity wherever economic terms are recomputable.
- Explicit missingness and zero semantics.
- Parity between structured results and rendered reports.
- Separate parent-pair reproductive economics.
- Resistance to unsupported external analyst suggestions entering the deterministic economic layer: 278 of 278 challenge items answered, none accepted into an authoritative result.

## What v1 does not establish

- Diagnostic sensitivity or specificity of whole-genome sequencing.
- Population prevalence or representativeness.
- Clinical efficacy for an individual patient.
- A universal cash value of a genome.
- Individualized medical or reproductive advice.
- Regulatory clearance or medical-device validation.

### Disclosed limitations of the release

- **Visual review: ACCEPTED P2 VISUAL EXCEPTION.** 49 of 9,482 rendered pages were flagged, all in the full-genome `report.pdf`; the economics PDFs were clean. The flags were a bounding-box false positive and minor ancestry-label overlap and clipping. This is not an unconditional pass.
- Mutation harnesses were not rerun because their source anchors were stale.
- No separate full manual adversarial-review pass was completed across every final artifact; one Lynch report was read end to end.
- 20 economic rows are explicitly classified as not independently recomputable from separable ΔQALY / ΔCost source terms, rather than being filled with invented components.

---

## Reproducibility record

| | |
|---|---|
| Benchmark | GenomeLens Deterministic Benchmark v1 (frozen) |
| Engine | `4eaf729` |
| Tooling | `a9b3d3a` |
| Cohort | 30 primary public genomes + 6 reference/parent genomes |
| Data | Public, open-consent research genomes only. No personal genome, private genomic data or bloodwork. |
| Archive | Frozen benchmark artifacts were archived and integrity-checked locally. Internal benchmark packages and raw run outputs are not publicly distributed. |

The production implementation, reference data assets and economic parameterization are maintained privately. This repository publishes what GenomeLens evaluates, what it outputs, the evidence that it behaves as described, the principles that govern its output, and where its evidence stops.

---

---

## Earlier reference demonstrations (superseded by Benchmark v1)

> **These are dated snapshots from an earlier engine build**, kept for the record. They are not the current benchmark: the figures below do not describe Deterministic Benchmark v1. Use [Benchmark v1](benchmark/v1/) for current results.

GenomeLens was run end to end on three publicly released Genome in a Bottle whole-genome callsets. The reports below are sanitized public outputs from those runs — real results, not mock-ups.

| Sample | Explicit truth calls | Concordance | Reference-case NMB | Canonical expected NMB | Incremental cost | Incremental QALY | ICER / status\* | Evidence basis | Report |
|---|---|---|---|---|---|---|---|---|---|
| **HG002** | 169 | **100.0%** (169/169, 0 mismatches) | **$158.72** | **$700.00** | $871.10 | 0.0103 | $84,588 | curated clinical, published association, published guideline + modeled economics | [Report (PDF)](artifacts/public-benchmarks/HG002/report-public.pdf) · [Web report](https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG002/report-public.html) · [JSON](artifacts/public-benchmarks/HG002/summary.json) |
| **HG003** | 187 | **100.0%** (187/187, 0 mismatches) | **$1,497.24** | **$3,588.00** | $388.53 | 0.0189 | $20,603 | curated clinical, published association, published guideline + modeled economics | [Report (PDF)](artifacts/public-benchmarks/HG003/report-public.pdf) · [Web report](https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG003/report-public.html) · [JSON](artifacts/public-benchmarks/HG003/summary.json) |
| **HG004** | 170 | **100.0%** (170/170, 0 mismatches) | **$138.56** | **$765.00** | $588.49 | 0.0073 | $80,943 | published association, published guideline + modeled economics | [Report (PDF)](artifacts/public-benchmarks/HG004/report-public.pdf) · [Web report](https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG004/report-public.html) · [JSON](artifacts/public-benchmarks/HG004/summary.json) |

> **\*** ICER shows the additional cost required to gain one additional
> quality-adjusted life year (QALY). Lower positive ICERs generally indicate better
> value at a given willingness-to-pay threshold. **Dominant** means better modeled
> outcomes at lower cost.

```
526 / 526  scoped explicit GIAB truth calls concordant
  0        mismatches
  0        normalization failures
```

Scope: GIAB **v4.2.1** · **GRCh38** · **chr1–22** · a fixed curated validation panel. This is agreement on the explicitly benchmarked calls inside that scope — not whole-genome accuracy, not clinical validation, and not an endorsement by NIST or GIAB. The denominator is published in [`docs/VALIDATION.md`](docs/VALIDATION.md#per-sample).

**What the columns mean**

- **Sample** — The public Genome in a Bottle (GIAB) genome analyzed by GenomeLens.
- **GIAB concordance** — Of the benchmark positions where the GIAB truth set states a
  definitive genotype (the *explicit truth calls*, the denominator shown), how many
  GenomeLens matched. A scoped validation measure over a fixed panel, not a claim of
  whole-genome accuracy.
- **Reference-case NMB** — *Net monetary benefit*: the modeled monetary value of
  expected health gains, minus the additional costs, under the reference-case
  assumptions. Higher positive values mean greater modeled net benefit at the stated
  willingness-to-pay threshold.
- **Canonical expected NMB** — A separate, standardized estimate of expected net benefit
  across the modeled pathways. It answers a different question from reference-case NMB,
  so it is not a component of that figure and the two are never added together.
- **Report** — The sanitized public output for that genome, as a rendered PDF, a web
  page, and the underlying JSON.

Reference-case and canonical expected NMB answer different modeling questions. They are
reported side by side rather than combined, because no single number is "the value of
this genome."

> **For nontechnical readers:** These results are a research demonstration of
> what GenomeLens can produce from real whole-genome data. The benchmark
> measures agreement only at the specific GIAB locations tested and is not a
> claim of whole-genome or clinical accuracy. The health-economic figures are
> modeled estimates based on published evidence and assumptions; they are not a
> diagnosis, medical advice, or a statement that a genome is literally “worth” a
> particular dollar amount.

<p align="center">
  <a href="https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG002/report-public.html"><img src="artifacts/public-benchmarks/HG002/preview.png" alt="GenomeLens public report, HG002" width="32%"></a>
  <a href="https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG003/report-public.html"><img src="artifacts/public-benchmarks/HG003/preview.png" alt="GenomeLens public report, HG003" width="32%"></a>
  <a href="https://conallaque.github.io/genomelens/artifacts/public-benchmarks/HG004/report-public.html"><img src="artifacts/public-benchmarks/HG004/preview.png" alt="GenomeLens public report, HG004" width="32%"></a>
</p>

Click any preview for the rendered report. **Report (PDF)** opens in GitHub's
own viewer; **Web report** is the same report rendered via GitHub Pages; **JSON**
is the machine-readable public summary. These reports are built from an explicit
allowlist, not by redacting a production report. The production GenomeLens report is a separate internal artifact and is not published.

### Genome in a Bottle trio concordance

GenomeLens is exercised end-to-end against the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) Ashkenazim trio — HG002, HG003, HG004 — using the publicly released GIAB GRCh38 benchmark callsets and defined confident regions. These are real public reference genomes with an external truth set, not mockups or simulations.

```
526 / 526  explicit benchmarked calls concordant
  0        mismatches
  0        normalization failures
```

Per sample: [HG002](https://www.nist.gov/programs-projects/genome-bottle) [169](docs/VALIDATION.md#hg002), [HG003](https://www.nist.gov/programs-projects/genome-bottle) [187](docs/VALIDATION.md#hg003), [HG004](https://www.nist.gov/programs-projects/genome-bottle) [170](docs/VALIDATION.md#hg004). The three are a father–mother–son trio, so the calls are not statistically independent.

**What that is not.** Not 100% whole-genome accuracy. Not clinical validation. Not an endorsement by NIST or GIAB. It is agreement on the explicitly benchmarked calls inside the defined truth and callability scope — a denominator we publish rather than hide.

---

## Validation and QA

[![Public validation](https://github.com/conallaque/genomelens/actions/workflows/public-validation.yml/badge.svg)](https://github.com/conallaque/genomelens/actions/workflows/public-validation.yml)

GenomeLens is developed against a private automated regression suite covering
genomic interpretation, evidence handling, health economics, uncertainty,
reproducibility and reporting. Current private-engine release run:
**7,723 passing · 26 skipped · 0 failed** for the frozen Deterministic Benchmark v1 build (engine `4eaf729`). The earlier release-build status that accompanied the trio demonstrations is recorded in [`docs/TESTING.md`](docs/TESTING.md).

That suite stays private, because a test asserting what the engine does under a
given input is a specification of the engine. So a **representative public
verification suite** ships here instead, written from the published artifacts
outward rather than by exporting private tests. It is purpose-built
around the public claims — not a sample drawn from the private suite — and lets
what this repository publishes be inspected independently. **Public validation:
161 checks — 161 passed.**

```bash
pytest tests/public                     # 161 checks: 161 passed
python tools/verify_public_release.py   # one-command release verification
```

Both read only files committed to this repository — no network, no credentials,
no engine dependency, so a fresh clone can run them.

[Public validation suite](tests/public) · [Testing methodology](docs/TESTING.md)

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

---

## Repository guide

| | |
|---|---|
| [`benchmark/v1/`](benchmark/v1/) | Deterministic Benchmark v1: evidence brief, results record, cover |
| [`docs/`](docs/) | Architecture, economics, glossary, limitations, scope, report spec, release status, testing, validation |
| [`tests/public/`](tests/public) · [`tools/verify_public_release.py`](tools/verify_public_release.py) | Public verification suite and one-command release check |
| [`artifacts/public-benchmarks/`](artifacts/public-benchmarks/) | Earlier trio reference reports (superseded by v1) |
| [`partner/OVERVIEW.md`](partner/OVERVIEW.md) · [`examples/`](examples/) | Partner overview and integration example |
