# Limitations

Published in detail because a system that states its boundaries is easier to evaluate than one that does not.

## Status

Research software under active development. **Not a diagnostic device. Not medical advice. Not clinically validated. Not payer-approved.** No regulatory clearance is held or claimed.

## Version scope

Published results describe specific builds: the frozen Deterministic Benchmark v1 build (engine `4eaf729`), the earlier build behind the trio snapshots, and the v1.11 release candidate, whose figures are stated with its build identifier in [`VALIDATION.md`](VALIDATION.md#v111-release-candidate). A result describes the build that produced it and no other.

## Validation scope

The external benchmark checks genotype concordance against a public truth set within a stated scope. The public evidence does not include an independent validation of the health-economic model, its parameters or its conclusions. Technical reproducibility is not clinical effectiveness, payer approval or external health-economic validation.

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

## Hypothetical values

A value whose best-ranked model is one built here (no published model of the decision ranks first for the person) is built from components and assumptions, labeled hypothetical and assumption-based, and shown as a scenario. It is not a supported estimate and is not counted in any evidence-supported total. A published model's outputs are reported as that source's results, under the source's own population, comparator, horizon and assumptions; GenomeLens reporting it does not validate the source.

## Stated age

Public genomes carry no age, so every public demonstration uses a standardized age, stated on every page: 45 for every genome, except the demonstration of a published-model decision, which states the age its source reports (30). Every age-dependent value is therefore a scenario for that age, not a prediction for the donor. Where a source reports one age and the report states another, the value is recomputed at the stated age by a model built here, labeled hypothetical, and the source's value for its own age is shown beside it; the source's value is not carried to the other age as the value.

## Callset limits

The demonstration genomes are whole-genome sequences, so the sequencing reaches the genes the reports discuss. What limits a statement about a gene is the callset the analysis starts from: a variants-only file records only where a person differs from the reference, so silence in it is not evidence of a reference sequence. Reports state, gene by gene, what the callset settles. A gene shown as not checked is a limit of the callset, not a negative finding.

## What a badge does not say

An **evidence-qualified** badge says that every model input is sourced and the value does not rest on an age that is stated rather than observed; it is a statement about provenance. It does not say that the treatment effect behind the value is established: a decision whose inputs are all cited can rest on an observational or contested effect. Where two models in the platform disagree about an effect, the decision's limitation text says so; the genotype-guided antiplatelet decision for CYP2C19 loss-of-function carriers is one (its one-year model uses an observational effect, and a randomized trial did not establish it). A scenario whose own effect estimate is documented as biased toward a larger effect is shown as a scenario and is never counted as a defensible alternative in the alternative-scenario uplift.

## Calls are shown as filed

A finding is shown as the reference databases file it. The reports do not print a ClinVar record's review status, a call's zygosity or, for a mitochondrial variant, its heteroplasmy level, so a "pathogenic" call in an apparently healthy adult (an X-linked or a mitochondrial variant, for example) is a research-genome call to be confirmed, not a diagnosis. Inheritance is stated where the carrier catalog holds it and printed as "uncertain" where it does not. A reproductive row prices the test that would resolve a risk, and a known-variant test is priced for a susceptibility allele (factor V Leiden, an IRF6 cleft-risk allele) only because the catalog treats it as a reproductive finding; a mild or modifier allele (HFE H63D) is classified as the catalog classifies it. None of this is a clinical interpretation.

## The colorectal surveillance example: two comparisons

Two different questions are kept apart for hereditary colorectal cancer. The **primary decision-relevant comparison** is the decision GenomeLens values: acting on a returned finding against not having it. Its value is the best-ranked model of that decision for the person, ranked by the validity of its evidence and its fit to the person's age, not by whether it is published. At the source's own age (30) that is the published analysis, restated to a per-carrier value. At 45 it is a model built here, labeled hypothetical and compared with no surveillance; the one published model that gives a value at 45 (Bennette 2015, from its printed figures restated to 2019 dollars) is shown beside it but not counted, because its penetrance is not corrected for ascertainment. The no-surveillance comparator is not what an unidentified person receives (average-risk screening, which starts at 45), so it overstates the increment a carrier faces, most at that age, by an amount that is not estimated (no average-risk arm is built); the model also credits colorectal cancer deaths averted only, with no morbidity or stage-specific treatment cost, which understates it, again by an amount not estimated, so the direction of the net bias is not known. Bennette's value is 2.2 times the counted one and the reasons are not separated: its penetrance is not corrected for ascertainment, which is why it is not counted, but its registered cumulative risk by 70 (40%) is not higher than the counted model's for an MSH2 man (50%), so penetrance does not by itself explain the gap. The effect size behind the counted value is partly an assumption: the hazard ratio for colonoscopy every 2 years (0.535) is the fitted 3-yearly 0.585 less the source's assumed 0.05 for each year shorter, about 12% more risk reduction than the fit. The interval printed beside the value leaves comparator uncertainty out (each printed interval is labeled parameter uncertainty only). The two primary values (at 30 and at 45) use different comparators and different models and are not comparable; the difference between them is not an age effect.

The **secondary published-strategy comparison** is a published analysis of a different question: which colonoscopy schedule is best for carriers who are already known. It has no "not identified" arm. Its optimal schedule gains health at a cost against the next efficient schedule on the source's own frontier, and gains much more, at lower cost, against the least intensive schedule it lists (current management for a family history, which the source itself rates dominated). Those are two comparisons of one schedule against different comparators; for MLH1 the net monetary benefits are 4,030 and 225,191 USD, and for MSH2, the carrier in the demonstration reports, 1,892 and 91,531 USD (more than fortyfold apart in both genes), and neither is the value of identifying a carrier. Both are recomputed from the source's printed totals by the scripts in [`reproduction/`](../reproduction/), in Python and in R (the table is the publisher's and is not shipped; the folder says how to supply it); one printed ratio in that table does not reproduce from the printed totals, and the kit says so. The source states 3% discounting, but its printed QALYs are larger than 3% discounting can give over its stated horizon, so the convention behind the printed totals is unresolved and every net monetary benefit here is arithmetic on the totals as printed. The source assumes perfect adherence. No interval is reported for a published result the source prints none for.

## Presentation limits of the demonstration reports

These are known and are not corrected in this release; none changes a figure.

- **Empty cells and pages.** A report can carry a row or a page with no figure. The all-findings page lists a population testing-strategy average as "per patient tested" with no dollar amount beside it (that average is shown beside a genome's values, never added to them), and the pages that would hold a pooled disease-model result print "Not available" when no finding of the genome enters that model (for example, several pages of demonstration 2). The layout does not yet suppress those pages.
- **Two tests with similar words.** The first page of the published-model demonstration says its inputs are supported, and later pages say the value is not evidence-qualified. The first is about how the inputs are sourced; the later ones are about the stated age, which is an assumption for every public genome.
- **Phenotype vocabulary.** The pharmacogenomic module names the CYP3A5 phenotypes "non-expressor", "intermediate expressor" and "expressor"; the CPIC guideline's term for the middle group is intermediate metabolizer.
- **The HEOR dossier and the six-page report are made by separate exports.** The public HEOR technical dossier and the featured public six-page report are not covered by the byte-for-byte fidelity check that the demonstration reports pass: a separate export prints a fixed list of fields, and its check is that every decimal, and every whole number of 100 or more apart from years and axis marks, in either document also appears in that genome's public reports, the public documents or the reproduction kit's recorded output (the public suite re-runs it). They omit the model's inputs (the six-page report's input table is replaced by an evidence summary), so the dossier's CHEERS 2022 items that call for them are marked partly or not reported.
- **Rows without an interval.** In the demonstration reports an alternative model or cascade row shown beside a counted value (for example, at age 45, the published rows of Bennette and of Guzauskas and the cascade rows) is not sampled: it shows no interval, and no sentence says why.

## Third-party material

The reproduction kit does not ship the printed totals of the published table whose arithmetic it recomputes (Kastrinos 2021, Table 2): the table belongs to its authors and publisher, the paper is an NIH author manuscript that is free to read but not openly licensed, and the kit tells the reader how to supply it. The documents quote a few printed results of that and other analyses, cited, so that the comparisons can be followed. The public genomes belong to their providers (Harvard Personal Genome Project participants who released their data, Genome in a Bottle, Illumina Platinum Genomes, the 1000 Genomes Project) and are used under each program's published terms for public data. Those terms and the participants' consent are the providers' and are not altered by this repository, and public availability does not by itself remove consent or privacy obligations; the reports attach inferred findings to the genomes as unconfirmed research demonstrations, not as clinical results.
