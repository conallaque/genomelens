# How to read a GenomeLens report

A GenomeLens economic report answers one question for one genome: **for each clinical decision that a finding could change, what would that decision be worth, against what comparator, and how well does the evidence support the answer?** This guide explains the pieces of the report so a number is never read without the words beside it. It describes the public demonstration reports; it does not describe how the production engine reaches them. The public demonstration reports in this repository are the allowlist export of one genome's results, 13 to 21 pages long. The documents in the next table are produced for each genome but are private; they are described here so that the numbers in the public export can be read in context.

## The set of documents

Each genome produces a six-page economic report, the detailed results behind it, and three short documents written for different readers.

| Document | Reader | What it holds | Published here |
|---|---|---|---|
| Economic report (six pages) | Anyone | The headline result and what it means; every finding by class; the calculation behind the headline; sensitivity and uncertainty; assumptions and limitations; evidence provenance | A public edition for demonstration 2 (`demonstrations/v1.11/featured/`; its input table and the detailed results are not in it); the full report is private |
| Detailed results | A reviewer | Every result by category, with its sources | No (private) |
| Payer brief | A health plan | The lead result, who it accrues to, a scale illustration per 1,000 members, what a plan should not read into it | No (private) |
| HEOR technical dossier | A health-economics analyst | A model card for every modeled value; the reporting checklist (CHEERS 2022); the validation record | A public edition for demonstration 1 (`demonstrations/v1.11/dossiers/`, without the input tables); the full edition is private |
| Genomics summary | A sequencing or genomics organization | What the genome yielded, from finding to value; carrier, cascade and population information; what the callset does and does not settle gene by gene | No (private) |

Each of the three short documents opens with a cover sheet: the headline, three key numbers, the main caveat, and a rights notice.

## The headline number

The headline is a **net monetary benefit (NMB)** for one decision: the health the decision gains, valued at a stated threshold per quality-adjusted life-year (QALY), minus the change in cost.

`NMB = (QALYs gained × threshold) − change in cost`

It is a health-economic value, not cash saved, and not a statement of what the genome is worth. The threshold is printed beside it. A positive NMB does not mean the decision saves money; the change in cost is printed separately.

## Where the number comes from: the model basis

Every modeled value says what kind of model stands behind it.

- **Published model.** The counted value is a published cost-effectiveness analysis's printed output for this decision at the person's age, restated by stated arithmetic where the source prints another unit (per person screened, for example). The page states which strategies are contrasted.
- **Secondary published-strategy comparison.** Where a published analysis compares strategies among people who are already identified (not identification against none), it is shown beside the decision as a labeled comparison with each comparator named. It is never the decision's value and never in a total.
- **Hypothetical, assumption-based.** The best-ranked model of the decision for the person is one built here. The value is built from published components and stated assumptions, is labeled hypothetical, and is shown as a scenario; a published value that is not counted (for another age, or with a known bias) is shown beside it.

Where several models describe one decision, one is counted and the others are listed beside it with the reason each was not counted (the "why each alternative was not selected" lines). The choice follows the validity of the evidence and its fit to the person, not whether a model is published.

## Classes of result

A genome's results are sorted into classes, and values in different classes are never added.

| Class | Meaning |
|---|---|
| Primary modeled | The value of the lead decision, at the stated age |
| Once the decision arises | The value if the decision comes up, not weighted by how likely it is to |
| Alternative model | Another model of the same decision, shown for comparison and not counted |
| Cascade testing (relatives) | Value that accrues to relatives, a different beneficiary |
| Conditional | Applies only in a stated context |
| Testing strategy (population) | The value of a testing strategy for a population, not this genome |
| Negative | A modeled decision that loses value |
| True zero | A measured null: the evidence finds no incremental benefit (a trial or a guideline; the interval is printed where the source gives one) |
| Unquantified or unresolved | No supported figure, with the reason stated |

A true zero is a finding. An unquantified result is not a zero: the report says why there is no figure.

## Badges

Beside each figure a badge states how well supported it is: **evidence-qualified** (every model input is sourced and the value does not rest on an age that is stated rather than observed), **partially supported**, **assumption-driven scenario** (it rests on a stated assumption, commonly the age), **limited** (the model would not be evidence-supported even at an observed age), or **unresolved** (in a page footer, a computed figure whose applicability to this donor is not established). A "supported calculation at the assumed age" box lists only the decisions whose every input except the stated age is sourced; an empty list is reported as empty, not as zero. When every decision in the box is a measured null, the box carries the same true-zero badge as those findings, because a measured null does not depend on the age. **Evidence-qualified records provenance, not the strength of the effect evidence:** an input taken from an observational study and one taken from a randomized trial are both "sourced", and the private editions' model cards show which source each input has.

## The stated age

These reports use one standardized age (45) for every genome, because public genomes carry no age; the one exception is the demonstration of a published-model decision, which states the age its source reports (30) so that the published result applies as printed. The age is a stated assumption, labeled on every page, and every age-dependent value is a scenario. It is not the donor's age and the results are not predictions for the donor.

## Reading a model card

A model card documents one modeled value. In the public HEOR technical dossier (`demonstrations/v1.11/dossiers/`) a card gives the decision and its comparator, the population, the perspective (whose costs count), the time horizon, discounting, price year, its key assumptions and its sources, and the value's basis (a published model's printed outputs, or a hypothetical scenario built here) where one is labeled. The private editions' cards also list each input with its value, evidence tier and source; those are not published. The reporting checklist (CHEERS 2022) in the dossier maps each of its items to where this document and the repository address it, and says which are partly addressed or not reported; it records where an item is addressed and is not a claim of compliance.

## Uncertainty

Where a source prints a range, the draws use it (rescaled to a per-carrier value where the source prints another unit); the 5th to 95th percentile is printed beside the value and is parameter uncertainty only. Where no range is registered for an input, it is held at its point, rather than inventing one; an alternative model shown beside a counted value is not sampled and prints no interval. The public demonstration pages print the cost per QALY gained beside each modeled headline value, with the QALYs gained and the change in cost it is calculated from and the threshold it is compared to, or state that the value dominates its comparator (dominance is a status, not a ratio).

## What the reports do not claim

Calls come from research genomes and have not been confirmed by a clinical laboratory. Modeled values depend on the population, comparator, time horizon and inputs the source model used. The reports are a research and decision-support product, not medical advice.
