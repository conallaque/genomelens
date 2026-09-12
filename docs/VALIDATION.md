# Validation

Factual status. No clinical, payer or regulatory validation is claimed.

## Automated tests

| | |
|---|---|
| Passing | **3,297** |
| Skipped | 16 (optional public benchmark fixtures not installed) |
| Known-failing | **10** |

The ten are held red deliberately. Each pins an unresolved upstream data defect
— curated records whose stored risk allele equals the GRCh38 reference base,
which would fire on essentially every genome. Silencing them would hide the
defect; they stay red until the data is corrected.

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
