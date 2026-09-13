# Validation

Factual status. No clinical, payer or regulatory validation is claimed.

## Automated tests

| | |
|---|---|
| Passing | **3,620** |
| Skipped | 16 (optional public benchmark fixtures not installed) |
| Known-failing | **10** |

The ten are held red deliberately and the count has not moved. Nine pin
curated records whose stored risk allele conflicts with the GRCh38 reference —
including risk alleles equal to the reference base, which would fire on
essentially every genome, and single-base alleles stored against insertions.
The tenth pins a stop condition covering curated sources that would otherwise
bypass validation. Silencing any of them would hide the defect; they stay red
until the underlying data is corrected.

This fixed set is also the change-control baseline. Every change is measured
against it, and a run is only accepted when the failing set matches it exactly
— not merely when the count matches.

## What the tests check

**Representation equivalence.** One biology encoded as an array export, a
block-compressed callset and an all-sites callset must produce identical
results — genotypes, dosages, callability and downstream economics.

**Fail-closed behaviour.** Unsupported genome builds refuse. Un-orientable
alleles refuse. Multi-base reference alleles inside a reference block refuse,
because a block does not exclude an indel the caller was never asked about.

**Semantic separation.** Conditional payoffs cannot be constructed wearing an
expected-value label. Values from different beneficiaries cannot be summed.
Unresolved quantities cannot be coerced to zero.

**Renderer purity.** Report modules are parsed and checked structurally for
arithmetic. A figure that appears on a page but not in the computed payload
fails the build.

**Reconciliation.** Cohort inputs must equal successes plus partials plus
failures. Every pathway must appear in exactly one economic state.

## Recent verified improvements

Behaviour verified by the test suite and by mutation testing. Each is an
implementation-correctness result; see the limitations note at the end of this
page for what that does and does not mean.

**Sequenced homozygous-reference is distinguished from absent.** In a
block-compressed whole-genome callset a homozygous-reference locus carries no
genotype row — it is silence inside a reference block. Read literally, that is
indistinguishable from a locus never examined. Consumers across the trait,
nutrition, pharmacogenomic-interaction, metal-handling and oxidative panels now
resolve the two separately, so a locus that was sequenced and found
reference-matching is no longer reported as untested.

**Canonical evidence access across assay generations.** Access to genotype and
dosage evidence is routed through a single canonical layer rather than reading
an array-era genotype column directly. Representation conventions inherited from
array exports no longer cause valid whole-genome evidence to be dropped, and the
canonical layer refuses to express a non-substitution as an allele string, so a
sentinel cannot be counted as if it were a base.

**Scoring-path disclosure for polygenic panels.** A polygenic panel now records
whether it was scored by genomic coordinate or by rsID fallback, and when the
fallback was used it records the reason. The scoring path is part of the result
rather than an unlogged internal detail.

**Ascertainment bias in instrument scoring.** Instruments used for
Mendelian-randomization-style analysis now separate callable
homozygous-reference loci from loci with no call. Scoring only the loci that
carry explicit records measures variant ascertainment rather than burden, and
biases the result in a direction that depends on the input representation.

**Threshold-qualified runs-of-homozygosity language.** Output no longer asserts
short or medium runs when none were detected, and the wording names the
thresholds a run was measured against. Where the input cannot support the
analysis, the result is withheld rather than reported as a zero the algorithm
never computed.

**Model defaults are distinguished from observations.** Person-level attributes
carry explicit provenance, and an attribute supplied by a model default cannot
be rendered as though it were observed. Defaults are labeled as such and
reported separately from measured inputs.

**Whole-genome hardening is in progress.** Work toward feature parity between
whole-genome input and the earlier array-based report is active and incomplete.
The items listed under *Active validation and integration work* below are part
of that effort and are not claimed as finished.

## Active validation and integration work

Listed for completeness. None of the following is complete, and none should be
read as a validated capability:

- **ABO blood typing.** The mapping from the canonical whole-genome observation
  to functional versus O-associated alleles is supported by coordinate,
  population-frequency and assembly evidence, but has not been demonstrated
  directly. It is not production-valid and is treated as open.
- **Direct HLA typing.** Current handling is tag-based. Direct typing is not
  implemented.
- **Specialized copy-number and structural-variant calling.** Not implemented.
- **Repeat-expansion calling.** Not implemented.
- **Economic remediation.** Work on unresolved and unpriced result states is
  ongoing and unfinished.
- **Full whole-genome feature parity.** Not reached.
- **Externally licensed polygenic score integrations.** Not integrated;
  licensing is unresolved.
- **External data redistribution terms.** Unresolved for some reference assets,
  which is why those assets are not part of this repository.
- **Production readiness.** Not claimed.

## Mutation testing

A passing test is not evidence until it has been shown to fail. Guards are
verified by planting the defect they exist to catch and confirming the guard
goes red — then removing the mutation.

This has repeatedly caught tests that passed without checking anything: a
comparison over an empty collection, a source-text assertion that dropped the
region it was meant to scan, a reconciliation whose numerator and denominator
came from the same filtered list.

## Shareability scanning

Artifacts intended to leave the building are scanned for personal data, local
paths, internal module references, raw genotype content and unsupported
commercial claims. Each scanner carries positive controls proving it can fire,
and negative controls proving it does not flag legitimate content.

## Limitations of this evidence

These are software correctness tests. They establish that the model does what
it says, consistently and reproducibly. They do not establish that the model's
clinical or economic assumptions are correct, and they are not evidence of
patient benefit.

**Validated implementation behaviour is not clinical validation.** Everything on
this page describes software behaviour verified against tests and a public
benchmark. It is not a clinical validation study, not an analytical validation
under a laboratory standard, and not a regulatory clearance of any kind.
GenomeLens is not a medical device and no diagnostic claim is made for it.
