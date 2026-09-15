# Capability matrix

What the engine does today, at capability level. Marker panels, routing
rules, parameter registries and qualification gates are not published.

Status vocabulary:

| Status | Meaning |
|---|---|
| **VALIDATED** | implemented, asserted by test, and audited adversarially |
| **OPERATIONAL** | implemented and asserted by test |
| **PARTIAL** | implemented for a subset, with the boundary stated |
| **REFUSES** | deliberately declines to produce output |
| **OPEN** | known gap, tracked against a frozen acceptance spec |

---

## Genomic layer

| Capability | Status | Note |
|---|---|---|
| Variant callset parsing | OPERATIONAL | coordinate handling, allele orientation, ploidy, indel representation |
| Representation invariance | VALIDATED | array vs block-compressed vs all-sites callsets asserted to agree |
| Callability / reference-block semantics | VALIDATED | an absent explicit row is a typed state, not a missing genotype |
| Typed source states | VALIDATED | explicit call · callable reference · explicit no-call · no usable evidence · never interrogated |
| Structural / copy-number phenotypes | REFUSES | where a phenotype depends on copy number, no diplotype is claimed |
| Cross-modality observation identity | OPEN | one biological fact observed through two assay types must resolve to one canonical observation |

## Clinical layer

| Capability | Status | Note |
|---|---|---|
| Pharmacogenomic interpretation | OPERATIONAL | category level only |
| Hereditary / carrier interpretation | OPERATIONAL | category level only |
| Polygenic and absolute-risk context | PARTIAL | population-average parameters, not individual prediction |
| Clinical validation | REFUSES | none claimed; concordance with a truth set says genotypes are right, not that they are useful |

## Decision layer

| Capability | Status | Note |
|---|---|---|
| Decision-pathway mapping | OPERATIONAL | a finding with no intervention does not become a decision |
| Usual-care comparator requirement | VALIDATED | a decision usual care reaches anyway is not genomic value |
| Attribution against comparator | VALIDATED | comparator existence does not imply attribution |
| Transportability checks | OPERATIONAL | publication does not imply applicability to this population |

## Economic layer

| Capability | Status | Note |
|---|---|---|
| Terminal adjudication of every pathway | VALIDATED | 36 / 36 rows, 0 unresolved |
| One credit per health consequence | VALIDATED | single credit-attribution authority; 0 duplicate credit |
| Canonical economic identity at the producer | VALIDATED | identity authored where the decision model is created, never derived from a label downstream |
| Construct separation | VALIDATED | quantities of different economic kinds are never summed; no single scalar "genome value" exists |
| Absent is not zero | VALIDATED | absent totals serialize as absent |
| Conditional vs expected reporting | VALIDATED | both published, labelled, never interchanged |
| Beneficiary separation | VALIDATED | reproductive and family-cascade value never pooled into a personal total |
| Alternate-model handling | VALIDATED | a second model of one decision is reportable and non-additive |
| External source values | VALIDATED | a number a cited source reported is never presented as an engine output |
| Uncertainty reporting | OPERATIONAL | point estimates carry intervals where the source provides them |
| Headline monetary total | REFUSES | 4 distinct economic pathways currently have defensible numeric estimates within conditional/scenario frames; 0 are currently eligible for a portfolio headline total |
| Countervailing-harm representation | OPEN | one intervention with both a benefit and an opposing harm is not yet fully representable; mandatory before any pathway becomes creditable |

## Survival / lifetime layer

| Capability | Status | Note |
|---|---|---|
| Life-table survival engine | VALIDATED | vendored national period life table, checksummed, reconstructs published life expectancy |
| Measure semantics | VALIDATED | central death rate and annual probability are typed distinctly, so the conversion for one cannot be applied to the other |
| Integration convention | VALIDATED | the quadrature rule the published table was built with, not an arbitrary approximation |
| Competing risks | OPERATIONAL | composed on the hazard scale, never by adding annual probabilities |
| Disease-specific lifetime models | PARTIAL | archetype selection complete; first archetype parameterised, not yet released |

## Verification

| Capability | Status | Note |
|---|---|---|
| Automated test suite | VALIDATED | 4,864 passing at the September 15, 2026 validation snapshot · 10 deliberately red · 0 unexpected |
| Mutation testing | VALIDATED | every core invariant must fail at least one test when deliberately broken |
| Adversarial review | VALIDATED | independent reviewers instructed to assume the work is wrong |
| Fail-closed inputs | VALIDATED | unsupported inputs refuse rather than coerce |
| Report reconciliation | VALIDATED | rendered figures must trace to the computed payload; renderers may not compute |
| Production reachability auditing | VALIDATED | safety guards are proven reachable by import-graph trace and behavioural probe, not by reading a docstring |
| Shareability scanning | OPERATIONAL | artifacts scanned for personal data, private paths and unsupported claims, with positive controls |

## Privacy

| Capability | Status | Note |
|---|---|---|
| Local-first inference | OPERATIONAL | a model whose name marks it cloud-hosted is refused before genomic content is assembled |
| Execution-mode separation | OPERATIONAL | personal, public-benchmark and synthetic contexts are distinct typed modes |
| Privacy-policy enforcement across generated artifact types — hardening in progress | OPEN | enforcement is being extended so that policy is applied uniformly at artifact-generation time |
| Certification | REFUSES | no HIPAA, SOC 2, HITRUST, FDA or security certification claimed or implied |

---

**Three entries above say OPEN and are load-bearing.** Cross-modality
observation identity, countervailing-harm representation, and uniform
privacy-policy enforcement across artifact types are remaining
release-hardening items, named here rather than omitted. A capability
matrix with no OPEN rows would be a marketing document.
