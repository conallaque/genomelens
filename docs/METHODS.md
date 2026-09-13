# Methods — philosophy

What the model is trying to be right about, and what it refuses to guess.
Parameterization, thresholds and routing rules are not published.

## Fail closed

Every stage refuses rather than approximates. Unsupported input is not coerced
onto a plausible default; an unresolved parameter does not become a zero; an
un-orientable allele does not become a non-carrier.

The rule that follows from this: **absence of evidence is never reported as
evidence.** A locus with no population-frequency lookup is not "rare". A
genotype with no sourced trigger probability has no expected value. A pathway
with no comparator has no attribution.

## Unresolved is not zero

Zero is a claim — that acting here is worth nothing. Unresolved is a different
claim — that the evidence does not yet support a number. A model that emits
zero for both is systematically overconfident in exactly the places it knows
least.

Unresolved pathways are counted, reported, and each names the specific
parameter, identity or attribution problem blocking it.

## Attribution against usual care

The question is not "what is correct treatment worth?" but "what did
sequencing add over the route that would have happened anyway?" A finding that
usual care would have caught has real clinical value and little *genomic*
value. Reporting the first as the second is the most common way genomic
economics overstates itself.

Total management value and genomic-attributable value are carried under
separate names so that no caller can print one while meaning the other.

## Conditional versus expected

A conditional payoff is what a decision is worth **if** its triggering
circumstance occurs. An expected present value is that payoff weighted by the
probability the circumstance occurs, and discounted. They differ by the
reciprocal of the exposure probability — often a factor of fifty or more.

Both are reported. Neither is allowed to wear the other's label, and the two
are never summed.

## Beneficiary separation

Value accruing to a payer, a relative, or a prospective child is not the
sequenced individual's gain. Beneficiaries are tracked per pathway and never
pooled into a single figure attributed to one party. A total spanning several
beneficiaries is labeled as spanning them.

## One decision, one value

A clinical decision reachable through two genomic routes is one decision. It is
counted once, and its value is claimed once.

## Observed, unobserved, and absent

A locus that was examined and matched the reference is not the same as a locus
that was never examined, and neither is the same as a locus the assay could not
resolve. Different input representations express these states differently: a
block-compressed callset records a reference-matching locus as silence inside a
reference block, while a variant-only callset records the same silence for a
locus it says nothing about at all.

Reading absence literally collapses these states into one, and the collapse is
not neutral — it converts an unexamined locus into a negative result. The
resolution layer therefore carries the state explicitly rather than inferring it
from whether a row exists, and what licenses a confident-reference reading is a
producer-supplied region definition, not the shape of the file.

## Representation invariance

The same biology encoded differently must be read the same way. This is
asserted by running both encodings through the pipeline and comparing results,
not by assuming the resolver is correct.

## Reproducibility

Every rendered figure traces to a computed payload; every payload figure traces
to a named parameter with a provenance tier. Renderers are formatters and are
structurally prevented from computing. A render produced from a modified
working tree says so on its face.

## Reproductive value, and whose it is

An earlier version of this page said reproductive outcomes are not priced at
all, on the reasoning that attaching a figure to an affected birth prices a
prospective child. That position has changed and the reasoning behind it has
not.

The reproductive pathway is now valued, under a **separate beneficiary arm**.
What is priced is the decision — partner carrier screening and the reproductive
options it opens — not a child, and not an outcome for a child. The distinction
that carries the weight is structural rather than rhetorical: the value accrues
to a prospective child, is reported as a prospective child's, and is never
summed into what the sequenced person's genome is worth to them. A total
spanning beneficiaries would answer a question nobody asked.

Pricing it required one probability the model had previously refused to guess —
how likely the conditioning circumstance is at all. It is now taken from a
national fertility survey, with its vintage and the direction of its likely
error both recorded, because a probability that multiplies a payoff carries its
bias straight into the result.

**Pricing the pathway did not settle what the genome added to it, and the two
are reported separately.** Carrier screening for this condition has been
recommended for all couples considering pregnancy since 2002, so routine
preconception care realises an unmeasured share of the same value without any
genome involved. What sequencing adds is the residual — and published screening
uptake is wide enough, and inconsistent enough in what it counts, that the
residual cannot currently be pinned to a number. So this pathway carries a
modeled expected value and an open attribution at the same time. The expected
value is labeled as the value of the *pathway*, not of the genome, everywhere
it is reported, and it is excluded from the genome-attributable total by
construction rather than by convention.

The attribution is open, not closed. It would be settled by a measurement that
does not appear to exist yet in published form: the proportion of patients in
routine care who actually complete carrier screening for this condition, with
the denominator defined as patients seen rather than patients offered. If that
figure is published, the residual becomes calculable and this pathway gains an
attributed value. Until then the expected value here is a modeled quantity
about a decision, and the share of it created by sequencing is unmeasured —
which is a statement about the evidence available, not a claim that the
question is unanswerable.

## What is not modeled

Some pathways remain deliberately unvalued and are reported as decisions rather
than omissions. Where a pathway is blocked, the block names the specific
measurement that would end it — and where no measurement could, it says that
instead of reading as unfinished work. One family-cascade input is of that kind:
how many relatives a particular person can reach is a fact about one family, not
a population rate, and substituting a national average would be wrong in a
specific direction for every person it was applied to.
