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
