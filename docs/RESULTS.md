# Results

Measured output on a public, open-consent reference genome. No raw genotype
content is reproduced here. No clinical, payer or regulatory validation is
claimed.

## Real-genome benchmark — snapshot, 12 September 2026

The figures in this section and in *Result qualification* below are a dated
benchmark snapshot, not a description of the engine's current state. They are
retained unchanged as a historical result. For current development status see
*Health-economic model status* below.

The benchmark input is the Genome in a Bottle / NIST **v4.2.1** high-confidence
truth set for **HG002 / NA24385 / huAA53E0**, on GRCh38 (no-alt analysis set).
The subject is a public reference genome released under PGP-Harvard open
consent.

| | |
|---|---|
| Truth set | GIAB / NIST v4.2.1 |
| Reference build | GRCh38, no-alt analysis set |
| Subject | HG002 / NA24385 / huAA53E0 |
| Consent | PGP-Harvard open consent, public release |
| Chromosomal scope | **chr1–22 only** — no chrX, chrY, chrM |
| High-confidence regions | **481,622** |

**Scope is corroborated, not assumed.** The chr1–22 restriction is not taken
from the release notes. It was confirmed independently by reading the tabix
index of the callset, which reports exactly **22** sequences. Every scope
statement in this document rests on that check.

**Why this genome.** It was selected because it carries benchmarking
infrastructure: a published truth set, an explicit high-confidence region
definition, and a second independent assay of the same individual. Those three
properties are what make a concordance rate meaningful rather than decorative.
The subject's ancestry is incidental provenance of a public dataset. It was
**not** a selection criterion, it is not a variable in any result below, and it
is not offered as evidence of generalisation across populations.

**What this input is, and what it is not.** A high-confidence benchmark truth
set is a curated set of variant calls within regions the producer is willing to
stand behind. It is **not** a complete whole-genome variant representation.
Regions the producer excluded are not represented at all — not as reference and
not as absent. For that reason this input is never called "full WGS" here or in
the outputs it summarises, and any rate computed from it is a rate over
benchmark-defined regions.

## Technical validation

**490** loci were requested against the benchmark. Results are reported as two
separate measurements, because they are not equally strong evidence and a single
blended rate would overstate the weaker one.

**Genotype concordance against explicit records.** Every locus where the truth
set carries a variant call, compared against what the engine independently
parsed from the same file.

| | |
|---|---|
| Explicit-record loci compared | **169** |
| Agreements | **169** |
| Concordance | **100.00%** |
| Mismatches | **0** |
| Normalization failures | **0** |

This is the measurement that exercises coordinate handling, allele orientation,
ploidy and indel representation.

**Confident-reference resolution — a callability result, not a concordance.**
A further **239** loci resolved as confidently homozygous reference rather than
unresolved once the producer's confident-region definition was supplied.

| Evidence type | Count |
|---|---|
| — single-nucleotide | 229 |
| — deletion | 4 |
| — insertion | 4 |
| — multiallelic | 2 |
| **Confident-reference total** | **239** |

These 239 are deliberately **not** folded into a concordance rate. The engine
emits *confidently homozygous reference* when there is no record at a position
inside the confident region; the truth definition asserts *homozygous
reference* from the same two facts. Their agreement therefore confirms that the
region logic fires where it should — it does not independently verify a
genotype, and counting it as concordance would inflate the figure.

**Independent cross-assay concordance.** The external check is a second,
independent measurement of the same individual by a different technology: a
consumer genotyping array, joined by rsID.

| | |
|---|---|
| Mutually scoreable loci | **316** |
| Agreements | **315** |
| Concordance | **99.68%** |
| Genuine discordances | **1** |
| Representation artefacts (not counted as mismatches) | **3** |
| Concordance if artefacts counted as mismatches | **98.75%** |

The scoreable loci here are exclusively single-nucleotide: every insertion,
deletion and multiallelic row resolved to no-probe, off-chromosome, or an
indel placeholder encoding the array could not express as sequence. **99.68% is
an SNV figure, not a whole-callset one.**

**Callability semantics.** A locus absent from a variant-only callset has not
thereby been shown to be reference. A variant-only callset carries **zero
reference blocks**: absence is the file's default state for everything the
producer did not call, including everything the producer never examined.
Absence alone is therefore not evidence.

The evidence that converts absence into a finding comes from elsewhere — the
producer's own high-confidence region definition. Inside those regions the
producer asserts it would have called a variant had one been present, so absence
becomes a positive statement of homozygous reference. Outside them, absence
remains uninformative.

GenomeLens keeps five distinct states and never collapses them into a single
"missing".

| Observed condition | Interpretation | Count |
|---|---|---|
| Present in callset | Explicit genotype evidence | **169** |
| Absent, inside confident region | Confidently homozygous reference | **239** |
| Absent, outside confident region | **Outside truth scope** — not reference, not absent, not a negative finding | **9** |
| Outside benchmark chromosomes | **Not assayed** | **72** |
| Indel resolved as dosage-only | **No allele string** | **1** |
| **Requested total** | | **490** |

The no-call is a single locus where a record exists in the source but could not
be reconciled to the curated representation. It is reported as a no-call rather
than resolved by preference.

**Effect of supplying the region definition.** Without the confident-region
input, only explicit records can be evaluated at all, and **239** curated loci
resolve as unresolved. Supplying the definition resolves them as confidently
homozygous reference. Nothing was relaxed to achieve this: the evidence that
licenses the call arrived in a separate file, and the engine refuses the call
without it.

**The rule inverts by file type, and the two look identical.** In an all-sites
callset, reference blocks carry the homozygous-reference evidence directly, so
an absent row is genuinely silent — the block already spoke for that position.
In a variant-only callset, an absent row is silent for the opposite reason: no
block spoke for it at all. At row level the two cases are indistinguishable, and
they mean opposite things. Applying the all-sites rule to a variant-only file
manufactures reference calls; applying the variant-only rule to an all-sites
file discards real ones. The resolution layer decides by file type, not by row.

**Second assay.** The same individual has an independent consumer genotyping
array. It was joined to the benchmark by **rsID**, not by coordinate: the array
is on GRCh37 and the benchmark on GRCh38, so position is not a valid join key
between them, and using it would produce alignment artefacts indistinguishable
from genotype disagreement.

| Category | Outcome |
|---|---|
| Genuine discordance | **1**, reported as discordant |
| Representation artefacts, correctly classified | **2** |

**One genuine discordance remains.** It falls in a gene with well-documented
pseudogene homology, a context in which short-read genotyping is known to be
difficult. It is reported as a discordance rather than normalised away, because
a rule written to absorb it would also absorb real disagreements elsewhere.

**Two further disagreements were representation artefacts, not genotype
mismatches** — one strand-orientation difference and one indel placeholder
encoding. Both were classified as such by the resolution layer rather than
counted as errors. Classifying an artefact correctly and suppressing a
discordance are different operations, and the count above separates them.

## End-to-end GenomeLens output

A full pipeline run was executed on the benchmark truth set.

| | |
|---|---|
| Runtime | **190 s** |
| Variants matched | **123** |
| Risk alleles | **77** |
| Clinically classified pathogenic / likely-pathogenic | **1** |
| Actionable | **0** |
| Carrier | **0** |
| Findings reaching the economic engine | **21** |
| Report consistency | **0 errors / 4 warnings** |

These are pipeline counts describing what the run produced. They are not a
clinical characterisation of the subject, and nothing here should be read as a
clinical finding about an individual.

**Y-DNA and mtDNA were not assayed.** The benchmark covers chr1–22 only, so
those analyses had no input. They are reported as not assayed, which is distinct
from reporting them as negative.

**What the run does and does not support.** No whole-genome-specific analysis
path was exercised: **0 of 21** findings required evidence available only from
whole-genome data. Two statements follow, and only the first is supportable.
**Supportable:** a whole-genome benchmark was run end to end. **Not
supportable:** the whole-genome-specific value path was exercised.

The counts in this section are from the end-to-end run and are not commensurable
with the technical-validation counts above; the two runs differ in scope and
their totals are not subsets of one another.

## Result qualification — benchmark snapshot, 12 September 2026

The qualification states below are read directly from the engine's own output.
They are not a retrospective assessment of the results, and they were not
assigned after the fact.

| Qualification state | Count | Share |
|---|---|---|
| UNRESOLVED | **33** | 91.7% |
| REFUSED | **3** | 8.3% |
| EVIDENCE_QUALIFIED | 0 | 0% |
| PARTIALLY_SUPPORTED | 0 | 0% |
| ASSUMPTION_DRIVEN | 0 | 0% |
| SYNTHETIC | 0 | 0% |
| **Total results** | **36** | 100% |

| | |
|---|---|
| Bookable results | **0** |
| Headline-eligible results | **0** |

These 36 qualification results are a different unit of account from the 21
findings that reached the economic engine; no ratio between the two is
meaningful and none is computed here.

**The finding is the refusal.** The engine processed a real benchmark genome and
declined to publish a single dollar figure, because every result failed at least
one evidence gate. That is the governance layer behaving exactly as designed —
the same behaviour the test suite asserts and the same rule stated in the
methods: absence of evidence is never reported as evidence, and unresolved is
never coerced to zero.

A number would have been easy to produce. Weakening any one gate would have
yielded a publishable figure, and would have destroyed the only thing this run
demonstrates: that the gates hold against a real genome, not only against
fixtures built to test them.

## Health-economic model status — development status, 14 September 2026

GenomeLens executes an economic evaluation pathway for **36 of 36** evaluated
result rows (100%).

The economic layer distinguishes between findings that support a numerical
decision model and findings whose correct output is conditional, reproductive,
non-applicable, evidence-limited, or explicitly refused. **Missing evidence is
not represented as a zero-dollar benefit.**

| | |
|---|---|
| Economic function coverage | **36 / 36 result rows** |
| Distinct economic pathways behind those rows | **17** |
| Function errors | **0** |
| Rows dropped without accounting | **0** |
| Placeholder `$0` values used to represent missing economics | **0** |
| **Pathways carrying a defensible monetary output** | **4 of 8 with published economics** |

**Rows are not pathways.** Each pathway emits one row per economic quantity, so
a pathway generally appears twice. Counting rows as independent models would
overstate the size of the economics by roughly a factor of two, and both units
are therefore carried separately throughout.

### Disposition breakdown

| Disposition | Rows |
|---|---|
| Conditional — pharmacogenomic, contingent on an exposure that has not occurred | **12** |
| Awaiting comparator-model completion — monogenic | **10** |
| Evidence gap — a named parameter is unsourced, or an effect was sought and not found | **5** |
| Scenario-modelled — two credible comparators answering different questions | 3 |
| Explicitly refused under an evidence gate | 3 |
| Not economically applicable | 2 |
| Reproductive arm — value accrues to a separate beneficiary | 1 |

Two of these counts moved for reasons worth stating, because both moved
in the direction of a *worse-looking* number that is a truer one.

*Conditional fell and evidence gap rose.* Two pharmacogenomic pathways were
reclassified out of "conditional" after their intervention effects were checked
against primary sources. Conditional asserts that real value awaits a trigger.
For these two, randomised trials had looked for the effect and not found it — so
they are waiting on evidence, not on a prescription, and saying otherwise would
promise a payoff those trials refuted.

*A mechanism-routing defect was corrected.* Ten rows had been classified under a
category that did not describe them. The correction moved them into a category
that names a missing comparator, which is the honest state and a less flattering
one.

### Monetized outputs

Monetary outputs are now reported in two tiers, which cannot be added together
and are kept apart for that reason.

| Tier | Pathways | Quantity |
|---|---|---|
| Per-patient conditional net benefit | **2** | **$987 total**, per patient tested |
| Signed verdict from a published ratio alone | **2** | one clears its threshold, one misses it |

The per-patient tier requires six conditions, all disqualifying: both cost and
health components present, a result whose interval excludes zero, an
intervention effect that was measured rather than refuted, a realistic
usual-care comparator, a stated denominator, and no double-counting against
another pathway. Four of eight candidate pathways fail at least one.

The second tier exists because a published cost-effectiveness ratio cannot yield
a per-patient net benefit, but it can yield a signed, scale-free one. **One of
the two verdicts is negative** — a pathway that misses its threshold by a wide
margin — and it is retained. A portfolio that reports only its positive results
is not reporting its evidence.

**Value is conditional, not realised.** Nobody in this configuration is taking
any of these drugs, so realised value is zero, and that zero is a modelled
result rather than an unknown. Converting conditional value into an expected
value would require a lifetime exposure probability, and no such figure was
found published for any of these drugs.

**Two honest qualifications on the coverage figure.**

*Exhaustiveness is a design property, not a measurement.* A row that cannot be
classified is routed to an explicit refusal category rather than to an
unclassified bucket, because a result nobody can classify is itself an error and
naming it one is more useful than leaving a hole in the accounting. "Zero
unclassified" therefore follows from the design and is reported as such. What is
measured, and what is worth reporting, is that the refusal category stays small
and holds only genuine refusals.

*This is a development configuration, not the benchmark genome run.* The 36 rows
here and the 36 results in the benchmark snapshot above are different units of
account over different inputs; the equal totals are coincidental. No ratio
between them is meaningful and none is computed.

**A functioning economic model is not a validated dollar estimate.** GenomeLens
separates model execution, evidence qualification and monetization rather than
publishing unsupported values to increase coverage. Coverage above means 36 of
36 rows are economically **classified**. It does not mean 36 of 36 are
economically **validated**, and it never will: most of these pathways should
not carry a dollar figure on the evidence that exists.

**What the validation pass changed, and in which direction.** Registered
intervention effects were checked against primary sources rather than against
the engine's own tables. Three of five were not what the field requires — a
relative reduction in the modelled event from genotype-guided management against
usual care — and two of those asserted effects that randomised trials had looked
for and failed to find. An independent systematic review of 47 published
economic evaluations reaches the same verdicts on the same pathways, which is
the closest thing to an external check available here.

Two citations were also found not to support the values attached to them: one
named a paper that does not exist as an economic evaluation, and one named a
screening analysis that holds treatment constant across its arms and therefore
carries no management evidence at all. Both were corrected and the affected
figures are now marked unsourced.


## Row-level validation — all 36

Every finding carries its own terminal validation decision. No aggregate
substitutes for a row, because an aggregate can show that *some* findings were
checked and never that a *particular* one was.

| Terminal status | Rows |
|---|---|
| Validated refusal — the reason for refusing is scientifically established | **14** |
| Validated conditional — a real model awaiting a trigger that has not occurred | **12** |
| Validated scenario — credible comparators answering different questions | **7** |
| Validated not economically applicable | 2 |
| Validated reproductive cost-consequence | 1 |

Nothing remains partial, unresolved, unknown or pending. **A validated refusal is
a finished result, not a deferral** — the question is whether the reason survives
contact with the source, not whether the answer is favourable.

**Two counts that a summary would hide.** Fourteen refusals is the largest single
category, and that is the honest shape of this evidence base rather than a
shortfall. And the model-gap count reached zero **without reclassification**: the
ten affected rows resolved on established reasons — six because the treatment is
indicated on the phenotype without a risk calculation at all, so the genomic
finding crosses no threshold the phenotype has not; four because published
management arms exist and are not transferable across jurisdiction, price year
and discount rates.

### A required invariant this build does not meet

An internal standard required zero unsourced material parameters. **Twenty-three
of thirty-six rows still contain at least one model default**, and that is
asserted as a failing count rather than resolved by relabelling. A regression
test pins the number so it cannot quietly be claimed as zero.

The parameter driving most of it is lifetime drug-exposure probability, and the
work done on it illustrates why the count is hard to move.

### Lifetime exposure: reconstructed, and still an interval

No lifetime exposure probability is published for any drug modelled here. Rather
than reuse a round default, five of eight pathways now carry a bracketed
derivation from sourced annual utilization. The reconstruction brackets rather
than chooses, because the literature does not resolve the distinction it turns
on:

> Annual **use prevalence** counts everyone taking a drug this year, including
> people who started years ago. Ever-exposure accumulates **initiations**.
> Treating one as the other compounds continuing users year after year.

| Pathway | Annual | Reconstructed lifetime bracket | Width |
|---|---|---|---|
| Antiplatelet | 1.37% | 1.4% – 42.4% | 31× |
| Anticoagulant | 0.53% | 0.5% – 19.0% | 36× |
| Statin | 2.14% | 2.1% – 57.9% | 27× |
| Antidepressant | 13.2% | 13.2% – 99.7% | 7.5× |
| Immunosuppressant | 0.28% | 0.3% – 10.5% | 38× |

**The width is the finding.** A bracket spanning thirtyfold constrains without
determining, which is why no ex-ante expected value appears anywhere in this
engine: weighting a conditional value by a probability known only to within
thirtyfold produces precision that is entirely fictional.

What it does establish is narrower and still worth having — every registered
default falls *inside* its own bracket. That does not validate the defaults, since
an interval that wide contains almost any plausible value. It rules out the
specific failure of a default sitting outside the range its own evidence permits,
which was not previously known.

**Three pathways returned a displayed zero and are recorded as not found.** The
source shows `0` alongside "no data available", meaning the survey could not
estimate that drug — not that nobody receives it. Reading those as zero exposure
would convert an unmeasured quantity into a measured absence, which is the same
error, in a different register, as reading a missing genotype as a reference
call.

## Where GenomeLens declined to add value

Two controls exist specifically to stop the total rising for reasons that are
not evidence. Both are reported because a system that only ever counts upward
should not be trusted with the count.

**One clinical outcome receives one economic credit.** Several findings in a
single genome can resolve to the same intervention preventing the same class of
event — a polygenic score, several single-gene variants and a second score can
all reduce to lowering one person's cholesterol to prevent one atherosclerotic
event. Booked independently they summed to more than **$1.1 million** of
prevented disease for someone with one cardiovascular system. Credited once per
distinct intervention-and-outcome pair, the figure is **$260,000**, so roughly
**$930,000 was double counting**. The discriminator is the pair, not the gene:
two findings can share a gene and be economically distinct, and two unrelated
genes can collide completely.

The suppressed findings are still reported. A score that confirms what a
single-gene variant already indicated raises confidence and may change intensity
or timing. What it is not is a second prevented event.

**Sequencing is bought once.** Each pharmacogenomic pathway carried its own test
cost, and across the registered genes those summed to **$2,100** charged to one
person for information that arrived in a single sequencing run. Each pathway's
arithmetic was locally correct and the portfolio total was wrong, which is why no
single model caught it. Acquisition is now allocated once at portfolio level,
and only confirmatory testing — the rarest case — is charged per finding.

## Where GenomeLens refused to guess

Three representative stops, paraphrased. Internal identifiers, parameter values
and gate internals are not published.

**Missing comparator.** A pharmacogenomic entry prices an averted adverse event
rather than the difference between the genotype-guided route and realistic usual
care. With no usual-care arm registered, there is no incremental quantity to
compute — only a gross one wearing an incremental label. The result is
**UNRESOLVED**. This is the attribution rule applied to a live result: what
sequencing added, not what correct treatment is worth.

**Missing survival evidence.** No pathway in this run carries registered
survival evidence. Life-years gained is therefore unavailable rather than
assumed, and no mortality benefit is imputed from effect estimates that were
never registered for that purpose. The result is **UNRESOLVED**.

**Verified equation, unverified linkage.** A cardiovascular risk equation is
verified, and its economic linkage is not. Three specific pieces are missing:
the risk-reclassification analysis is absent, the treatment threshold is
unregistered, and the composition from relative to absolute effect is assumed
rather than sourced. A verified predictor does not license an unverified
economic chain built on top of it, so the calculation is **REFUSED** rather than
completed with the missing pieces supplied by default.

## Limitations

**Truth validation covers benchmark-defined regions only.** Nine requested loci
fall outside truth scope and are excluded from every rate reported here. They
are not counted as failures and not counted as successes.

**Technical validation is not clinical validation.** The concordance figures
establish that the pipeline reads a benchmark genome correctly and consistently.
They establish no clinical efficacy, no causal benefit and no economic effect.

**Model coverage is incomplete.** Not every pathway carries survival evidence,
and the qualification counts above reflect that directly.

**Unresolved results are retained as unresolved by design.** They are not
pending items awaiting a default value. A future release may resolve some of
them with sourced evidence; none will be resolved by relaxing a gate.

**Upstream file authenticity is unverified.** This release publishes no upstream
checksums, so the integrity of the downloaded benchmark files is asserted by the
producer and not independently confirmed here.

**Legacy value-of-information figures are excluded.** Figures produced by the
legacy value-of-information path are not evidence-qualified and appear nowhere
in this document.

**Proprietary implementation detail is intentionally omitted.** Parameters,
thresholds, routing rules and module structure are not published, here or
elsewhere in this repository.

## Reproducibility

**The benchmark is public.** The GIAB / NIST v4.2.1 truth set for HG002 on
GRCh38 is downloadable from the GIAB/NIST FTP release. Reproducing the
technical-validation section requires the truth set, its index, and the
high-confidence region definition from the same release.

**The region definition is an opt-in input.** GenomeLens accepts the producer's
confident-region definition as an optional input to the resolution layer.
Without it, behaviour is unchanged: absence stays uninformative, the scoreable
denominator stays at the explicit records, and nothing is inferred to fill the
gap. Supplying it is what licenses the confident-reference state, and the effect
of supplying it is reported above as a denominator change, not as a concordance
change.

**Denominators are always stated.** Every rate in this document names the
population it is computed over. A rate over loci found and a rate over loci
requested differ by exactly the failure mode being measured — loci that were
never evaluable. Reporting only the first makes unevaluable loci disappear from
the arithmetic, which is why the requested total, the scoreable total and the
excluded categories are all carried separately and sum to 490.

**Invocation is not published.** This document describes capability, not
commands. Command-line options, module paths and internal file names are not
part of this repository, in keeping with the rest of its documentation.
