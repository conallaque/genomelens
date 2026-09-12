# GenomeLens

### Economic Intelligence for Genome Sequencing

**Turning sequence data into measurable decision and economic value.**

*Not just what your genome says — what it could change.*

---

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

## A worked example

**DPYD / fluoropyrimidines**, from a synthetic genome:

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
| Automated tests | **3,297 passing**, 16 skipped |
| Known-failing | 10, held red deliberately — each pins an unresolved data defect rather than being silenced |
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
- All figures and artifacts here are generated from synthetic data.

## What is public and what is not

This repository is **documentation and synthetic artifacts**. It does not
contain the production engine.

Selected production methods, reference assets, routing logic and economic
parameterization are intentionally withheld from the public repository. The
intent is credibility without recipe: enough to evaluate whether the approach
is sound, not enough to reconstruct it.

## Documentation

[Architecture](docs/ARCHITECTURE.md) · [Methods](docs/METHODS.md) ·
[Validation](docs/VALIDATION.md) · [Partner pilot](docs/PARTNER-PILOT.md)

## License

See [LICENSE](LICENSE). Documentation and synthetic artifacts only.
