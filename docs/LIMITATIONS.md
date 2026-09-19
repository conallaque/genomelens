# Limitations

Published in detail because a system that states its boundaries is easier to evaluate than one that does not.

## Status

Research software under active development. **Not a diagnostic device. Not medical advice. Not clinically validated. Not payer-approved.** No regulatory clearance is held or claimed.

## Benchmark scope

The reported concordance covers positions where the truth set states a genotype, within defined confident regions. It is not whole-genome accuracy, and it is not an endorsement by NIST or GIAB. Technical concordance is not clinical validity.

## Modality

A SNV/indel variant callset cannot resolve every genomic event. Copy number, structural variants, repeat expansions, some star-allele architectures, HLA typing and mitochondrial heteroplasmy are **not claimed as production-supported**. Capabilities requiring read-level evidence are not exercised by a variant callset alone, and results are labeled accordingly rather than presented at uniform confidence.

## Difficult regions

Segmental duplications, homopolymers, high-homology gene families and repeat-rich regions are harder for any short-read pipeline. Callability is reported rather than assumed.

## Absence of evidence

A locus never interrogated is not a locus that came back normal. GenomeLens distinguishes explicit calls, callable reference, explicit no-calls and untested positions — but a negative result is only as strong as the assay's coverage of that position.

## Pharmacogenomic conditionality

Guideline actionability is not the same as a decision in play. Where the relevant medication or indication is unknown, the finding is reported and the economics are withheld. Published economics therefore **may** understate what would be realized if prescribing context were known. The direction is not signed by the model: a withheld pathway with negative net monetary benefit at the stated threshold would move the total the other way.

## Model dependence

Results depend on assumed treatment effects, event costs, utility decrements, discount rate, time horizon and willingness-to-pay threshold. Where an assumption's plausible range changes the *sign* of a result, that is disclosed rather than resolved by choosing a convenient value.

## Ancestry and population

The benchmark itself is one nuclear family of a single ancestry. It establishes nothing about callability in populations with a different variant spectrum or divergent haplotype structure.

Effect estimates and allele frequencies are not uniformly transportable across populations. Where an estimate is transferred from a source population to a different target, that transfer is an assumption.

## Synthetic versus real data

Synthetic inputs test behavior at scale. They are not evidence of accuracy, prevalence, or clinical validity, and they never substitute for real-genome validation.

## Economic uncertainty

Modeled economic results are not realized savings. They are consequences of a model, conditional on its inputs, and should be read as such.
