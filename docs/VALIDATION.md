# Validation

> The concordance results on this page are historical: they come from an earlier engine build and are published as a dated snapshot. They are not representative of v1.11 validation, and they do not describe Deterministic Benchmark v1 (see [`benchmark/v1/`](../benchmark/v1/)) or any later version, including the current private development version.

## Real-data benchmark

GenomeLens was run end-to-end, in an earlier engine build, against the Genome in a Bottle Ashkenazim trio — HG002, HG003, HG004 — using the publicly released GIAB GRCh38 benchmark callsets and their defined confident regions.

```
526 / 526  explicit benchmarked calls concordant
  0        mismatches
  0        normalization failures
```

Per sample: HG002 169, HG003 187, HG004 170.

## Per sample

<a id="hg002"></a>

### HG002

**Sample ID** — HG002 / NA24385 (son), [GIAB Ashkenazim trio](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/).

**Benchmarked calls** — 169, all concordant. Zero mismatches, zero normalization failures.

**Truth and callability scope** — the GIAB GRCh38 benchmark callset for this sample, restricted to its defined confident regions.

The 169 counts **explicit calls compared within the GIAB benchmark scope** — positions where the truth set states a genotype. Positions the truth set declines to state are reported as callability, not as agreement.

<a id="hg003"></a>

### HG003

**Sample ID** — HG003 / NA24149 (father), [GIAB Ashkenazim trio](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/).

**Benchmarked calls** — 187, all concordant. Zero mismatches, zero normalization failures.

**Truth and callability scope** — the GIAB GRCh38 benchmark callset for this sample, restricted to its defined confident regions.

The 187 counts **explicit calls compared within the GIAB benchmark scope** — positions where the truth set states a genotype. Positions the truth set declines to state are reported as callability, not as agreement.

<a id="hg004"></a>

### HG004

**Sample ID** — HG004 / NA24143 (mother), [GIAB Ashkenazim trio](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/).

**Benchmarked calls** — 170, all concordant. Zero mismatches, zero normalization failures.

**Truth and callability scope** — the GIAB GRCh38 benchmark callset for this sample, restricted to its defined confident regions.

The 170 counts **explicit calls compared within the GIAB benchmark scope** — positions where the truth set states a genotype. Positions the truth set declines to state are reported as callability, not as agreement.

Per-sample counts sum to the 526 reported above (169 + 187 + 170).

---

## The denominator, stated

The figure counts positions where the truth set **states a genotype**. Positions the truth set declines to state are not counted as agreement — they are reported separately as a callability result. Counting them toward concordance would inflate it.

This matters more than the percentage. A concordance rate whose denominator is unpublished is not a measurement.

**Precision.** 526 comparisons with zero errors is not evidence of a 100% rate. By the rule of three, the one-sided 95% lower bound is approximately 99.4% pooled, and approximately 98.3% at the per-sample counts of roughly 170. The point estimate is reported with that bound in mind.

**Direction.** The figure measures agreement where both the truth set states a genotype and the pipeline reports one. It is a genotype-concordance measure, not a recall measure: a position where the truth set states a genotype and the pipeline emits nothing is a separate quantity, and `0 mismatches` cannot detect it.

**Independence.** HG002 is the son of HG003 and HG004. Roughly half of HG002's alleles are present in a parent by descent, so the 526 comparisons are not 526 independent observations, and the trio is one nuclear family of a single ancestry rather than a population sample.

**What the benchmark is not:** not whole-genome accuracy; not clinical validation; not an endorsement by NIST or GIAB. Technical concordance says a genotype matches an external truth call within a defined scope. It says nothing about clinical utility.

## Out of benchmark scope

**Chromosome scope.** The GIAB small-variant benchmark used here covers the autosomes, chr1–22. chrX, chrY and chrM therefore have **zero benchmark coverage** — a different statement from the modality limits in [`LIMITATIONS.md`](LIMITATIONS.md#modality), which concern what a variant callset can resolve at all. Nothing here is evidence about sex chromosomes or mitochondrial variants.

Regions outside the confident-region definition; variant classes the truth set does not adjudicate; and any capability requiring read-level evidence rather than a variant callset. These are stated rather than silently folded into a pass.

## Provenance

Published artifacts carry the engine build that produced them, the input callset and truth-set version, and the region definitions applied. A result that cannot be attributed to a specific build is not publishable as evidence.

**The trio figure does not yet meet that standard.** The truth-set release (GIAB v4.2.1), reference build and chromosome scope are published with it. The confident-region file identity, the engine build identifier, the retrieval date and the rule by which the benchmarked positions were selected are not. Until they are, the figure should be read as a reported result rather than as an independently re-derivable one. See [`RELEASE-STATUS.md`](RELEASE-STATUS.md#benchmark-manifest).

## Synthetic stress testing — PLANNED, NOT YET RUN

A synthetic whole-genome **variant-call** cohort is designed to test behavior at scale: failure modes, output distributions, economic coverage and systems behavior across many inputs.

It is explicitly **planned and not yet executed**, and when it runs it will answer a different question from the GIAB benchmark. Synthetic inputs test stability and coherence at scale; they are not evidence of real-world accuracy, disease prevalence, or clinical validity, and they do not substitute for real-genome validation.

It is [not run in this release](RELEASE-STATUS.md#synthetic-cohort), and is not described as existing.

## v1.11 release candidate

<a id="v111-release-candidate"></a>

The figures below describe the v1.11 release candidate built from engine commit `e89f545c`. They come from the final validation run of that build, after its last change. They are software and engineering results; none is evidence of clinical effectiveness, and none validates the health-economic models against observed outcomes.

| Check | Result |
|---|---|
| Automated software tests (full private suite) | 11,165 passed, 68 skipped, 1 failed (an environmental disk-space failure in a test that starts a pipeline from sequencing reads), full private suite on a clean snapshot of the final tree |
| Every economic document equals what its renderer makes from the saved payload alone | 300 of 300 documents byte-identical (the six-page report, the three editions, the findings-first sheets and the public report of each genome) |
| Every printed figure traces to a stored result | every printed figure traces to a stored result, in every format checked |
| Differences from the validated release candidate | 23 economics payloads compared field by field with the validated release candidate's; 0 unexplained differences. Every difference is traced to an approved change. |
| PDF rendering (collisions, clipping, unreadable text on every page) | 424 PDFs, 16888 pages rendered; 0 pages with defects |
| Public export audit | 16 public PDFs audited (15 at the standard stated age and 1 at the age its source reports); 0 findings |
| Build | `e89f545c` |

**The 6 published analyses behind the registry's published rows** (rows of the registry, not independent studies): BRCA2 carriers, prostate cancer screening (a cost-effectiveness analysis for the UK National Screening Committee), UK NSC / SENSS (ScHARR) 2026 (16 rows); genomic findings in healthy adults sequenced for another indication: Lynch syndrome, BRCA1/BRCA2, malignant hyperthermia, hypertrophic and arrhythmogenic cardiomyopathy, long QT (Bennette 2015), PMID 25394171 (6 rows); MSH6 and PMS2 carriers, gynecologic surgery, PMID 34499134 (2 rows); BRCA1/BRCA2, risk-reducing surgery and surveillance, PMID 33119106 (1 row); Lynch syndrome, colonoscopy surveillance after a returned finding (Guzauskas 2022), PMID 35227606 (1 row); familial hypercholesterolemia, statin treatment, PMID 35961838 (1 row).

**Kinds of checking, kept apart.**

- *Automated software tests* assert what the software does under given inputs.
- *Regression tests* pin results that were reviewed and accepted, so a later change that moves one is noticed.
- *Economic calculation reproduction*: one published model (a hemochromatosis cost-effectiveness model) reproduces its source's printed untreated-cohort outputs when rebuilt from the source's parameters (the treated arm is not tested by those outputs); that check is made in the private suite and has no public artifact here. Separately, the arithmetic of one published table (hereditary colorectal cancer surveillance, 20 printed strategies) is recomputed from its printed totals in a second language: which strategies are dominated, the efficiency frontier, each ICER, the optimal strategy at 100,000 USD per QALY and its net monetary benefit against two comparators. All 20 classifications and each gene's optimal strategy agree with the source; three of the four printed ICERs reproduce and one does not; and the printed QALYs are larger than the source's stated 3% discounting can give over its horizon, so its discounting convention is unresolved. The Python and the R scripts are in [`reproduction/`](../reproduction/) and print the same recorded output; the public suite checks it. This verifies arithmetic against a source and does not rebuild the source's model or show that it predicts outcomes.
- *Cross-check*: the platform's own construction for one decision (hereditary colorectal cancer surveillance against no surveillance) was set beside the published schedule comparison strategy by strategy. The two agree on the ordering of strategies. The platform's model credits a smaller surveillance effect (a fraction of the published quality-adjusted life-year gain, and higher cancer incidence in 16 of 20 comparisons); their horizons, risk sources and cost inputs differ, the questions differ, and the source's discounting convention is unresolved (see above). The gap is stated, not tuned away: the platform's inputs were not adjusted to match the published results.
- *Genome-processing validation*: the frozen Deterministic Benchmark v1 above describes the earlier build; the release candidate was run end to end on 50 public genomes.
- *PDF rendering validation*: every page of every generated report is checked for overflow, clipping and collisions, and each document for its required sections and its figures' state lines.
- *External validation against observed outcomes*: none. No model here has been validated against observed outcomes in a population.

## What this repository lets you verify

<a id="what-this-repository-lets-you-verify"></a>

**You can:** run both reproduction scripts on the published table (which you supply, as `reproduction/README.md` describes) and match their recorded output; recompute each net monetary benefit and cost-effectiveness ratio in the demonstration reports from the quality-adjusted life-years and the change in cost printed beside it; run the public verification suite; verify the signed Content Credentials of the two research PDFs; compare the printed results with the papers behind the cited PMIDs; and see the class of every result.

**You cannot:** reproduce any value from a genome (the engine, the registry, its inputs and the reference data are private); reproduce the hypothetical values (their sources and horizon are named; their structure, effect fractions, costs and probabilistic distributions are not); inspect the registry itself or the private test suite; reproduce the Deterministic Benchmark v1 run or the 50-genome run; or obtain any validation beyond arithmetic. A question of that kind is for the repository owner.
