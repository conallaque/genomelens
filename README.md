# GenomeLens

### Economic Intelligence for Genome Sequencing

**Turning sequence data into measurable decision and economic value.**

*Not just what your genome says — what it could change.*

---

## Results on a real benchmark genome

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
confident-reference, 1 no-call (a record present but unreconcilable), 9 outside
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

### Economic qualification: 0 of 36

Of the 36 results the qualification layer evaluated, none cleared the bar to
be valued:

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
| Automated tests | **3,430 passing** |
| Known-failing | 10, held red deliberately — each encodes a known curated-data defect. Turning them green without fixing the cause would delete the only record that the defect exists. |
| Representation equivalence | array vs block-compressed vs all-sites callsets asserted to agree |
| Fail-closed | unsupported inputs refuse rather than coerce |
| Mutation testing | guards are re-verified by planting the defect they exist to catch and confirming they fail |
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
[Architecture](docs/ARCHITECTURE.md) · [Methods](docs/METHODS.md) ·
[Validation](docs/VALIDATION.md) · [Partner pilot](docs/PARTNER-PILOT.md)

## License

See [LICENSE](LICENSE). Documentation and synthetic artifacts only.
