# GenomeLens — Genomic Decision & Health-Economic Engine

**Turning sequence data into measurable decision and economic value —
and refusing to produce a number when the evidence behind it does not
hold.**

---

## Scorecard — 15 September 2026

```
36 / 36   economic pathway rows terminally adjudicated
0         unresolved rows
0         duplicate economic credit
0         fake-zero economic outputs
0         fabricated health outcomes in a published quantity
15        distinct economic decision pathways behind those 36 rows

Canonical economic-credit aggregation          VALIDATED
Remaining release-hardening items               IN PROGRESS

4,864     passing at the September 15, 2026 validation snapshot
10        held deliberately red (each encodes a known curated-data defect)
0         unexpected failures

Synthetic and public-reference validation only
Local-first / offline-capable architecture
```

Every figure above is produced by the engine's own audit tooling, not
transcribed by hand. Two of them deserve reading twice: **0 duplicate
economic credit** and **0 fake-zero outputs** are guarantees about
arithmetic, and the section below explains why they are the hard part.

**Release hardening is in progress and is reported as such.** Canonical
economic-credit aggregation is validated. Remaining release-hardening
items are open work, tracked internally against a frozen acceptance
specification. Readiness is not claimed beyond what is stated here.

---

## The chain

```
Genome
  ↓
Canonical genomic observations
  ↓
Clinical interpretation
  ↓
Actionable decision pathways
  ↓
Health outcomes
  ↓
Costs / QALYs / life-years
  ↓
Economic result
  ↓
Evidence + uncertainty + provenance
```

Each arrow is a gate that can refuse, and the monetary result is the
**end** of that chain rather than the beginning. A finding with no
intervention behind it does not become a decision. A decision usual care
would have reached anyway does not become *genomic* value. A pathway with
no sourced probability for its triggering event does not become an
expected value.

---

## What GenomeLens refuses to do

These are enforced by tests, not by policy documents. Each refusal exists
because the opposite behaviour is the easy, plausible-looking default —
and each one makes the headline number smaller.

- **It does not monetize unsupported pathways.** A pathway whose causal
  chain has not been established carries no monetary value, however
  complete its parameters are.
- **It does not treat missing evidence as zero.** `0` means *measured,
  and found to be nothing*. Absent means absent. Collapsing the second
  into the first invents precision that was never measured.
- **It does not double-count two genes informing one decision.** Two
  pharmacogenes governing a single drug dose are two genomic findings and
  **one** economic consequence. Three genes converging on one
  lipid-lowering decision are one decision, not three.
- **It does not merge conditional results into unconditional headline
  values.** "Worth $X *if* this drug decision arises" and "worth $Y
  today" are the same pathway answering different questions. They are
  reported separately and never silently interchanged.
- **It does not convert reproductive outcomes into personal net monetary
  benefit.** Value accruing to a prospective child is not the sequenced
  person's gain, and the two are never pooled into one figure attributed
  to one party.
- **It does not emit unsupported economic values.** Where a required
  input is absent the result is *unavailable*, with the specific missing
  evidence named — not rounded to zero, and not filled with a default.

The engine currently publishes **no** headline monetary total, because no
pathway has cleared the bar. That is the governance layer holding, not
the pipeline failing.

---

## Results on a real benchmark genome — snapshot, 12 September 2026

GenomeLens turns a variant callset into decision and economic intelligence, and
refuses to produce a number when the evidence behind it does not hold. It was
run end to end against a public, open-consent reference genome with a published
truth set, and cross-checked against a second, independent assay of the same
individual. The economic layer then published **nothing** — every result failed
at least one evidence gate. Both of those are the intended outcome.

| | |
|---|---|
| **Benchmark** | GIAB / NIST v4.2.1, GRCh38 |
| **Subject** | HG002 / NA24385 / huAA53E0 (public, open consent) |
| **Scope** | chr1–22, 481,622 high-confidence regions |
| **Runtime** | 190 s over the benchmark truth set |

### Genotype concordance

Two comparisons, reported separately because they are not equally strong
evidence and merging them would overstate the result.

**1 — Against the benchmark's explicit records.** Every locus where the truth
set carries a variant call, compared against what the engine independently
parsed from the same file:

```
169 explicit-record loci compared
  → 169 agree          100.00%
  → 0 mismatches · 0 normalization failures
```

This exercises coordinate handling, allele orientation, ploidy and indel
representation — where a resolver actually goes wrong.

**2 — Against a second, independent assay.** A consumer genotyping array of the
same individual, joined by rsID (the array is GRCh37 and the benchmark GRCh38,
so position is not a valid join key). This is the genuinely external check,
because the two measurements come from different technologies:

```
316 mutually scoreable loci
  → 315 agree          99.68%
  → 1 genuine discordance
  → 3 representation artefacts, classified as such rather than as mismatches
      (98.75% if all three were counted as mismatches instead)
```

The single discordance was **resolved against the array**. A second GIAB
benchmark built from a diploid assembly — not from short-read mapping, so not
subject to the same failure mode — makes the same call as the truth set. The
array probe sits in a region of high pseudogene similarity where
cross-hybridization is documented to produce spurious heterozygous calls. It is
disclosed rather than dropped.

### What the callability layer recovered

A variant-only callset carries no reference blocks, so an absent row is not
evidence of anything. Supplying the producer's confident-region definition let
**239** further loci resolve as confidently homozygous reference instead of
unresolved — the same claim the *WGS-native by design* section below makes
structurally, now measured on real data.

These 239 are reported as a **callability result, not as concordance.** Both the
engine and the truth definition decide "homozygous reference" from the same two
facts — no record, inside the confident region — so their agreement confirms the
region logic fires correctly; it does not independently verify a genotype.
Counting them toward a concordance rate would inflate it.

The 490 requested loci are fully accounted for: 169 explicit records, 239
confident-reference, 1 resolved as dosage-only (an indel carried as a copy
count rather than an allele string), 9 outside
truth scope, and 72 not assayed because they sit on chrX, chrY or chrM, which
this benchmark does not cover. 169 + 239 + 1 + 9 + 72 = 490.

### Pipeline output for the run

| | |
|---|---|
| Variants matched | 123 |
| Risk alleles | 77 |
| Clinically classified pathogenic / likely pathogenic | 1 |
| Actionable | 0 |
| Carrier | 0 |
| Findings reaching the economic engine | 21 |
| Report consistency | 0 errors · 4 warnings |

These are engine counts for this run, not a clinical reading of a person.

### Economic qualification for this run

Of the 36 results the qualification layer evaluated in this benchmark run,
none cleared the bar to be valued:

| State | Count |
|---|---|
| UNRESOLVED | 33 (91.7%) |
| REFUSED | 3 (8.3%) |
| EVIDENCE_QUALIFIED | 0 |
| PARTIALLY_SUPPORTED | 0 |
| ASSUMPTION_DRIVEN | 0 |
| **Bookable** | **0** |
| **Headline-eligible** | **0** |

A platform that wanted a number here could have produced one. The gates are
what stopped it:

```
finding present
  → no usual-care comparator registered
      → attributable value UNRESOLVED   (not $0)

finding present
  → no registered survival evidence
      → life-years gained UNAVAILABLE   (not 0 years)
```

Unresolved is not zero. Zero would assert that acting here is worth nothing;
unresolved says the evidence does not yet support a number, and names what
would unlock it. A run that publishes no valuations is the governance layer
holding, not the pipeline failing.

The dollar figures shown under *A worked example* below come from a synthetic
profile and are illustrative. No economic result from this benchmark run is
published, because none qualified.

### Current economic model status — 15 September 2026

The benchmark figures above are a dated snapshot. The economic layer has
since moved on, and its current status is reported in a different unit of
account: **result rows**, of which there are 36, behind **15 distinct
economic pathways**.

| | |
|---|---|
| Rows with a terminal adjudication | **36 / 36** |
| Unresolved rows | **0** |
| Rows dropped without accounting | **0** |
| Placeholder `$0` used to represent missing economics | **0** |
| Duplicate economic credit | **0** |
| Rows numeric within their own validated frame | **8** |
| Distinct pathways carrying a defensible numeric output | **4** |
| **Headline-eligible pathways** | **0** |

Stated precisely, because the distinction is the whole point:

> **4 distinct economic pathways currently have defensible numeric
> estimates within conditional/scenario frames; 0 are currently eligible
> for a portfolio headline total.**

A conditional estimate is real inside the frame that conditions it. A
headline is exactly where that condition would be dropped, which is why
the second number is zero and not a rounding of the first.

**A published count was wrong and is corrected here.** An earlier
revision of this page reported *17* distinct economic pathways and
"4 of 8 with published economics". Both figures came from counting
**rows** where the unit is **pathways** — several rows can share one
economic pathway, and one pathway appeared five times. The corrected
figures are **15** distinct pathways, of which **4** carry a defensible
numeric output and **8 rows** are numeric within a validated frame. The
audit that certified pathway coverage was itself miscounting, in the
reassuring direction, which is the direction worth being suspicious of.

**These two 36s are not the same 36.** The benchmark snapshot above counts
results from a public reference genome; this counts rows from a
development reference configuration. The equal totals are coincidental
and no ratio between them is meaningful.

**Adjudicated is not validated, and neither is monetized.** All 36 rows
carry a terminal decision. Zero are headline-eligible. Coverage is not
offered as a substitute for either.

Full detail: [benchmark results](docs/RESULTS.md) ·
[capability matrix](docs/CAPABILITY-MATRIX.md).

### One decision, one credit

The hardest problem in this layer is not computing a value. It is making
sure one health consequence is credited exactly **once**, however many
genomic routes reach it.

A coronary polygenic score, three familial-hypercholesterolaemia genes
and an ischaemic-stroke score are five different **findings**. They are
substantially **one** intervention — lower this person's LDL cholesterol
— preventing **one** class of event in a person who has one
cardiovascular system. Booking them independently sums to a large,
individually defensible, entirely fictitious number.

The engine now resolves this with a single credit-attribution authority:

- economic identity is authored where the **decision model** is created,
  never derived downstream from a gene symbol, a display label or a
  condition name
- one aggregation path produces every total; no renderer or serializer
  recomputes one
- a result that cannot be shown *not* to double count is **excluded and
  reported**, never admitted with a synthesised identity
- quantities of different economic constructs are never added, so there
  is deliberately **no single scalar "genome value"**
- absent totals serialize as absent, never as `$0`

Deduplication semantics, identity registries and aggregation internals
are not published.

**Technical validation is not clinical validation.** Concordance with a truth
set says the genotypes are right; it says nothing about clinical utility.

Full detail: [benchmark results](docs/RESULTS.md).

## What it does

A variant callset is evidence. It is not yet a decision, and a decision is not
yet a value. GenomeLens connects the chain that turns one into the other:

```
genomic evidence
  → clinical context
    → decision pathway
      → usual-care comparator
        → health / economic consequence
          → genomic-attributable value
            → uncertainty and provenance
```

Each arrow is a gate that can refuse. A finding with no intervention behind it
does not become a decision. A decision that usual care would have reached
anyway does not become *genomic* value. A pathway with no sourced probability
for its triggering event does not become an expected value. What survives every
gate is reported; what does not is reported as unresolved, with the specific
evidence that would unlock it.

## Why this is different

Most genomic reporting stops at the finding. The distinctions below are the
product:

| | |
|---|---|
| **A finding is not a value** | A variant is only worth something if some decision changes because of it. |
| **Management value is not genomic value** | Correct treatment is worth something however it was discovered. Only the part sequencing *added*, over the route that would have happened anyway, is attributable to the genome. |
| **A conditional payoff is not an expected present value** | "Worth $11,600 if this drug decision arises" and "worth $38.91 today" are the same pathway at different questions. They differ by the probability the decision arises at all. |
| **A QALY is not added lifespan** | Quality-adjusted life-years are a weighted construct, not months on a calendar. |
| **Unresolved is not zero** | Zero asserts that acting here is worth nothing. Unresolved says the evidence does not yet support a number. Collapsing the second into the first invents precision. |
| **Beneficiaries stay separate** | Value accruing to a payer, a relative, or a prospective child is not the sequenced person's gain, and is never pooled into one figure attributed to one party. |

## What it analyzes

Pharmacogenomics · prevention and risk · hereditary findings · carrier and
reproductive context · testing replacement · family cascade · callability and
coverage · economic qualification.

Category level only. Marker panels, routing rules and qualification gates are
not published.

## WGS-native by design

On a whole-genome callset, most loci have no variant row — the caller emitted a
reference block instead. Treating that silence as "untyped" reports a fully
sequenced genome as unassayed.

GenomeLens distinguishes, as typed states:

- an explicit call
- callable reference evidence covering the locus
- an explicit no-call
- genuinely no usable evidence
- a locus an array never interrogated

**An absent explicit row is not a missing genotype.** The same biology encoded
as an array export, a block-compressed callset, or an all-sites callset is
required to produce the same answer, and that equivalence is asserted by test
rather than assumed.

The resolver itself is not published.

## Economic modeling

Conditional payoff, expected future value, genomic attribution against a
usual-care comparator, beneficiary separation, and explicitly unresolved
pathways. Parameterization, eligibility gates and attribution formulas are not
published.

## A worked example (synthetic)

**DPYD / fluoropyrimidines**, from a synthetic genome — illustrating the
distinctions above. This is not a benchmark result:

| | |
|---|---|
| Identity | verified |
| Payoff if the drug decision arises | **$11,600** |
| Probability that decision arises | **1.05%** |
| Expected present value | **$38.91** |

Both numbers are true and they are not interchangeable. A platform that reports
the first as value delivered overstates its product by roughly two orders of
magnitude. GenomeLens reports both, labeled, and lets the reader choose the one
their question needs.

*(A larger pharmacogenomic figure exists in the same synthetic profile. It is
not used as the flagship, because that allele is reached by tag proxy rather
than direct typing — and an example you have to caveat is not an example.)*

## Validation

| | |
|---|---|
| Automated tests | **4,864 passing at the September 15, 2026 validation snapshot** |
| Known-failing | 10, held red deliberately — each encodes a known curated-data defect. Turning them green without fixing the cause would delete the only record that the defect exists. |
| Representation equivalence | array vs block-compressed vs all-sites callsets asserted to agree |
| Fail-closed | unsupported inputs refuse rather than coerce |
| Mutation testing | guards are re-verified by planting the defect they exist to catch and confirming they fail. Guards that survived their own mutation — and therefore proved nothing — have been found this way and given real tests |
| Adversarial review | finished work is re-audited by independent reviewers told to assume it is wrong and to find the defect. Most of the serious findings in the latest round were in code written during that same round |
| Aggregation invariants | permutation invariance, idempotence and byte-identical determinism are asserted on economic totals |
| Report reconciliation | rendered figures must trace to the computed payload; renderers may not compute |
| Shareability scanning | artifacts scanned for personal data, private paths and unsupported claims, with positive controls |

No clinical validation is claimed. See limitations.

## Privacy

GenomeLens supports privacy-sensitive local workflows in which personal genomic
data need not be sent to external LLM services. Analysis runs against a local
inference endpoint, and a model whose name marks it as cloud-hosted is refused
before genomic content is assembled into a prompt.

No HIPAA, SOC 2, HITRUST, FDA or security-certification status is claimed or
implied.

## Partner / pilot

For a sequencing provider the shape is:

```
sequencing output → GenomeLens → per-genome economic intelligence
                              → cohort-level decision and evidence summaries
```

Per genome: which decisions the sequence can change, who benefits, how well
evidenced each one is. Across a cohort: pathway prevalence, economic-state
distribution, coverage, and an evidence-gap map naming what would unlock each
unvalued pathway.

Synthetic demonstrations: [partner overview](assets/partner-overview.pdf) ·
[pilot concept](assets/pilot-onepager.pdf) ·
[example payload](examples/genomelens_partner_v1.example.json).

No current partner deployment is claimed.

## Limitations

- Research and pilot stage. **Not a diagnostic.**
- Not a payer submission and not a regulatory evaluation.
- Not proof of population return on investment. Individual conditional payoffs
  are not summed into a cohort total; they refer to different unrealized
  triggers and are not commensurable.
- Population-average parameters, not individual predictions.
- Many pathways remain deliberately unresolved pending sourced evidence.
- Structural-variant limitations remain in specific areas. Where a gene's
  phenotype depends on copy number — CYP2D6 is the clearest case — this
  analysis does not measure it and does not claim a resolved diplotype however
  completely the gene's SNP loci are called.
- Two sources are mixed here and are labeled as such. The worked example,
  partner overview and example payload are generated from **synthetic data**.
  The benchmark section at the top is measured on a **public, open-consent
  reference genome** against its published truth set.

## What is public and what is not

This repository is **documentation, benchmark results and synthetic
artifacts**. It does not contain the production engine.

Selected production methods, reference assets, routing logic and economic
parameterization are intentionally withheld from the public repository. The
intent is credibility without recipe: enough to evaluate whether the approach
is sound, not enough to reconstruct it.

## Documentation

[Benchmark results](docs/RESULTS.md) ·
[Capability matrix](docs/CAPABILITY-MATRIX.md) ·
[Architecture](docs/ARCHITECTURE.md) · [Methods](docs/METHODS.md) ·
[Validation](docs/VALIDATION.md) · [Partner pilot](docs/PARTNER-PILOT.md) ·
[Changelog](CHANGELOG.md)

## License

See [LICENSE](LICENSE). Documentation and synthetic artifacts only.
