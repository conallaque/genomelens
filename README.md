# GenomeLens

### Applied Health Economics and Genomic Decision Analysis

GenomeLens is an independently developed research project on how genomic evidence can inform decision-relevant health-economic analysis. It follows the chain that conventional pipelines leave implicit: a genomic finding, the clinical decision it could change, the comparator that decision is judged against, the resulting health outcomes and costs, and the resources they consume. Where the chain is not supported by evidence, the analysis is designed to stop and say so.

The project does not assign an intrinsic monetary price to a genome. Economic value attaches to a decision, a comparator and a perspective. A genotype with no decision it could change has nothing to value.

**Documents:** [Research and Methodology Overview (PDF)](research/GenomeLens_v1.11_Research_and_Methodology_Overview.pdf) · [Illustrative HEOR case study (PDF)](research/GenomeLens_TSC_Illustrative_HEOR_Case_Study.pdf) · [Deterministic Benchmark v1](benchmark/v1/) · [Documentation](docs/)

---

## Research objective

Genomic interpretation and economic evaluation usually run as separate workflows. Interpretation answers *which variants are present*; economic evaluation answers *whether an intervention is worth its cost*, for a population and a decision defined in advance. Little connects the two.

GenomeLens investigates how to connect them while accounting explicitly for:

- clinical applicability to the person or population in question;
- the economic decision context, including comparator, perspective and time horizon;
- the quality of the evidence behind each step;
- incremental costs and incremental health outcomes;
- uncertainty, both parametric and structural;
- methodological limitations; and
- the circumstances in which an economic estimate is not supported and should not be produced.

This is an applied research investigation. It does not claim to have solved these problems in general.

## Current development status

| | |
|---|---|
| Current development version | **v1.11**: private, under active development, not a public software release |
| Public benchmark evidence | GenomeLens Deterministic Benchmark v1, frozen: engine `4eaf729`, tooling `a9b3d3a` |
| This repository | Selected research methodology, historical benchmark evidence and sanitized public examples |
| Not distributed | The production implementation and its associated research assets |

The benchmark results published here describe the frozen v1 build. They are results for that build only. They are not results for v1.11 and do not validate it. Work in progress in v1.11 is not described here as finalized.

## Research methods

GenomeLens uses standard health-economic methods, stated here at the conceptual level so that published results are interpretable.

| Method | What it provides |
|---|---|
| Cost-effectiveness and cost-utility analysis | Incremental cost (ΔC) and incremental QALYs (ΔQ) of a strategy against a stated comparator |
| Incremental cost-effectiveness ratio | `ICER = ΔC / ΔQ`, interpreted with the signs of ΔC and ΔQ; dominance is reported as a status, not as a ratio |
| Net monetary benefit | `NMB = λ·ΔQ − ΔC`, where λ is the willingness-to-pay threshold; meaningless unless λ is stated |
| Perspective | Whose costs and outcomes count: health-care sector, payer, societal or patient |
| Time horizon and discounting | The period over which consequences accrue, and the rate at which later costs and QALYs are discounted |
| One-way sensitivity analysis | The effect on the result of varying one input across a plausible range |
| Probabilistic sensitivity analysis | Parameter uncertainty propagated through the model, summarized as distributions rather than a point estimate |
| Value-of-information analysis | The ceiling on what resolving uncertainty could be worth (EVPI), and what one specific measurement could be worth (EVSI) |

Three questions are kept apart: whether a strategy lowers expected spending (**cost-saving**), whether the health it buys is worth its incremental cost at a stated threshold (**cost-effective**), and whether a budget-holder can pay for it at the scale proposed (**affordable**, a budget-impact question). A positive net monetary benefit does not imply cost savings, and a cost-effective strategy can still be unaffordable.

See [`docs/ECONOMICS.md`](docs/ECONOMICS.md) for definitions and the [Research and Methodology Overview](research/GenomeLens_v1.11_Research_and_Methodology_Overview.pdf) for a short exposition.

## Research applications

The project explores how this framework applies in several settings. These are research directions. They are not described as completed or validated production capabilities.

- Genomic medicine and precision health
- Rare-disease economics
- Pediatric genomic screening research
- Pharmacogenomics
- Reproductive and family decision economics
- Healthcare resource allocation
- Economic uncertainty and evidence appraisal

An [illustrative case study](research/GenomeLens_TSC_Illustrative_HEOR_Case_Study.pdf) shows the method on a rare disease using published clinical evidence and synthetic economic inputs. It is an educational calculation and not an output of the GenomeLens engine.

## Research integrity

These principles govern how results are interpreted and reported.

- **A computable result is not necessarily a defensible result.** A number that can be calculated still needs an evidence chain from finding to outcome.
- **Clinical actionability and individual applicability are different questions.** A finding can be actionable in principle and not applicable to this person now.
- **Missing evidence is not zero benefit.** An absent input is reported as absent. It is never rendered as `$0`.
- **Perspectives and beneficiaries must be read correctly.** Quantities that answer different questions, for different beneficiaries, are not added together.
- **Unsupported aggregation does not produce totals.** A sum is reported only where its components are comparable.
- **Uncertainty and assumptions stay visible.** The result carries the evidence basis and the assumptions it rests on.
- **Declining to estimate can be a legitimate result.** A refusal is reported with its reason. It is not a negative biological finding.

The vocabulary used to report these states (supported, conditional, unresolved, refused, not economically applicable) is defined in [`docs/GLOSSARY.md`](docs/GLOSSARY.md).

---

## Public benchmark and technical evidence

### Deterministic Benchmark v1

A frozen, public-data technical proof of concept: 30 primary public genomes and 6 reference or parent genomes, 36 in total, run end to end as a deterministic control (no AI interpretation) on pinned engine and tooling builds. The results below describe that specific historical build.

| Measure | Result |
|---|---|
| Public genomes completed | 36 / 36 |
| Tests passed (private suite, frozen build) | 7,723 passed · 26 skipped · 0 failed |
| Release gates | 9 / 9 |
| Parity failures | 0 |
| Arithmetic discrepancies | 0 |
| Real-row assertions | 301 / 301 |
| Zero-audit defects | 0 / 1,399 |
| Forbidden-string hits across rendered PDFs | 0 / 144 |

*Parity* means every dollar and percentage printed in a report traces to a structured result, so a report cannot say something the analysis did not. The *zero audit* examines every displayed zero and requires it to be a supported null or no-change state, never a missing value rendered as zero. *Arithmetic* recomputes each net-monetary-benefit value from its recorded components wherever those components are separable. *Real-row assertions* are checked against rows the pipeline actually produced, not fixtures.

Full record: [Public Evidence Brief](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Evidence_Brief.pdf) (narrative, 6 pages) and [Public Results](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Results.pdf) (compact technical record, 4 pages). The [benchmark cover](benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Cover.png) is kept with them.

### What v1 establishes

- Reproducible execution across the intended public cohort, from pinned engine and tooling builds.
- Evidence and economic routing consistency, with no finding dropped from the economics without a recorded reason.
- Arithmetic integrity wherever economic terms are recomputable.
- Explicit missingness and zero semantics.
- Parity between structured results and rendered reports.
- Separate parent-pair reproductive economics.
- Resistance to unsupported external analyst suggestions entering the deterministic economic layer: 278 of 278 challenge items answered, none accepted into an authoritative result.

### What v1 does not establish

- Diagnostic sensitivity or specificity of whole-genome sequencing.
- Population prevalence or representativeness.
- Clinical efficacy for an individual patient.
- A universal cash value of a genome.
- Individualized medical or reproductive advice.
- Regulatory clearance or medical-device validation.
- Any property of versions after v1, including the current private development version.

### Disclosed limitations of the release

- **Visual review: ACCEPTED P2 VISUAL EXCEPTION.** 49 of 9,482 rendered pages were flagged, all in the full-genome `report.pdf`; the economics PDFs were clean. The flags were a bounding-box false positive and minor ancestry-label overlap and clipping. This is not an unconditional pass.
- Mutation harnesses were not rerun because their source anchors were stale.
- No separate full manual adversarial-review pass was completed across every final artifact; one Lynch report was read end to end.
- 20 economic rows are explicitly classified as not independently recomputable from separable ΔQALY / ΔCost source terms, rather than being filled with invented components.

### Reproducibility record

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

## Earlier reference demonstrations (historical; superseded by Benchmark v1)

> **These are historical results from an earlier engine build**, kept for the record. They are not the current benchmark and are **not representative of v1.11 validation**: the figures below do not describe Deterministic Benchmark v1 or any later version, including the current private development version, v1.11. Use [Benchmark v1](benchmark/v1/) for current public results.

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

### Genome in a Bottle trio concordance (historical snapshot)

The concordance figures in this section come from the same earlier engine build as the snapshots above. They are historical and are not representative of v1.11 validation.

GenomeLens was exercised end-to-end against the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) Ashkenazim trio — HG002, HG003, HG004 — using the publicly released GIAB GRCh38 benchmark callsets and defined confident regions. These are real public reference genomes with an external truth set, not mockups or simulations.

```
526 / 526  explicit benchmarked calls concordant
  0        mismatches
  0        normalization failures
```

Per sample: [HG002](https://www.nist.gov/programs-projects/genome-bottle) [169](docs/VALIDATION.md#hg002), [HG003](https://www.nist.gov/programs-projects/genome-bottle) [187](docs/VALIDATION.md#hg003), [HG004](https://www.nist.gov/programs-projects/genome-bottle) [170](docs/VALIDATION.md#hg004). The three are a father–mother–son trio, so the calls are not statistically independent.

**What that is not.** Not 100% whole-genome accuracy. Not clinical validation. Not an endorsement by NIST or GIAB. It is agreement on the explicitly benchmarked calls inside the defined truth and callability scope — a denominator we publish rather than hide.

Sanitized reports for the current engine build are not published. The three earlier snapshots above are the only per-genome reports in this repository, and they are historical (see [`docs/RELEASE-STATUS.md`](docs/RELEASE-STATUS.md#trio-reports)).

---

## Public verification

[![Public validation](https://github.com/conallaque/genomelens/actions/workflows/public-validation.yml/badge.svg)](https://github.com/conallaque/genomelens/actions/workflows/public-validation.yml)

GenomeLens is developed against a private automated regression suite covering genomic interpretation, evidence handling, health economics, uncertainty, reproducibility and reporting. The figure published for the frozen Deterministic Benchmark v1 build (engine `4eaf729`) is **7,723 passing · 26 skipped · 0 failed**. The earlier release-build status that accompanied the trio demonstrations is recorded in [`docs/TESTING.md`](docs/TESTING.md). No test figure is published for any later version.

That suite stays private, because a test asserting what the engine does under a given input is a specification of the engine. A **representative public verification suite** ships here instead, written from the published artifacts outward rather than by exporting private tests. It is purpose-built around the public claims and lets what this repository publishes be inspected independently. **Public validation: 175 checks — 175 passed.**

```bash
pytest tests/public                     # 175 checks: 175 passed
python tools/verify_public_release.py   # one-command release verification
```

Both read only files committed to this repository — no network, no credentials, no engine dependency, so a fresh clone can run them.

[Public validation suite](tests/public) · [Testing methodology](docs/TESTING.md)

## What a public report contains

Each run produces a genomic and health-economic report containing supported findings, carrier status, pharmacogenomics, callability and source state; evidence provenance per finding; reference-case and standardized pathway economics; incremental cost, incremental QALYs, ICER or dominance; probabilistic sensitivity analysis; and explicit refusals and limitations. Reports are published only through the allowlist export profile in [`PUBLIC-REPORT-SPEC.md`](docs/PUBLIC-REPORT-SPEC.md).

## Where GenomeLens stops

- It does not monetize a pharmacogenomic finding merely because a genotype is guideline-actionable when the medication or indication is unknown.
- It distinguishes an unsupported absence claim from a true negative.
- It does not claim production support for every copy-number, structural, repeat-expansion, CYP2D6 star-allele or mitochondrial-heteroplasmy problem.
- It does not read a GIAB benchmark as clinical validation. Concordance establishes agreement with the external truth calls within the benchmarked scope; it says nothing about clinical utility.
- It does not treat economic arms as additive by default.
- Where a benefit is estimable but its real-world cost is not adequately sourced, it can report a maximum model-consistent cost bound rather than silently treating the unknown cost as zero or inventing an actual market price.

---

## Research limitations

GenomeLens is research software. Its economic outputs are model-dependent: they rest on assumed treatment effects, event costs, utility decrements, discount rate, time horizon and willingness-to-pay threshold, and where an assumption's plausible range changes the sign of a result, that is disclosed rather than resolved by choosing a convenient value.

Technical reproducibility does not establish clinical effectiveness, payer approval, or external health-economic validation. Public benchmark agreement is not clinical validity. Modeled results are not realized savings. Nothing in this repository has been independently validated beyond what is stated in [`docs/VALIDATION.md`](docs/VALIDATION.md).

**Research software. Not a diagnostic device. Not medical advice.** Full statement: [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).

## Creator role and AI-assisted development

The project was independently conceived and directed by its author: the research direction, the health-economic methodology, the model design, and the validation strategy and its oversight.

The software engineering and the preparation of documents were carried out with substantial AI assistance. That includes the code and the two research PDFs in [`research/`](research/), which carry signed Content Credentials recording that they were produced with an AI tool. The author set the requirements, reviewed the output and made the methodological corrections; neither the code nor the documents are presented as written by hand. The author's oversight is not independent scientific validation; what has been checked externally, and within what scope, is stated in [`docs/VALIDATION.md`](docs/VALIDATION.md). Details are in [`docs/PROJECT-SCOPE.md`](docs/PROJECT-SCOPE.md#creator-role).

## Public repository boundary

This repository contains selected methodology, validation evidence, curated reports and integration examples. The production engine, internal economic parameter registries, routing logic and qualification machinery are maintained separately. Public artifacts are generated through a deliberately constrained export layer so that validation evidence and interpretable results can be inspected without exposing the production engine.

---

## Repository guide

| | |
|---|---|
| [`research/`](research/) | Research and Methodology Overview and an illustrative HEOR case study (PDF) |
| [`benchmark/v1/`](benchmark/v1/) | Deterministic Benchmark v1: evidence brief, results record, cover |
| [`docs/`](docs/) | Architecture, economics, glossary, limitations, scope, report spec, release status, testing, validation |
| [`tests/public/`](tests/public) · [`tools/verify_public_release.py`](tools/verify_public_release.py) | Public verification suite and one-command release check |
| [`artifacts/public-benchmarks/`](artifacts/public-benchmarks/) | Earlier trio reference reports (historical; superseded by v1) |
| [`partner/OVERVIEW.md`](partner/OVERVIEW.md) · [`examples/`](examples/) | Technical overview and public schema example |
