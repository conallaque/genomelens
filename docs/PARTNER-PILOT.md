# Partner pilot

A narrow first evaluation for a sequencing provider. Research and product
evaluation — not a clinical study, not a payer submission.

## What the pilot asks

Deliberately **not** "does sequencing save money at population scale?" That
question needs evidence a pilot cannot produce. Instead:

1. Can sequencing outputs be processed consistently?
2. Which decision-relevant pathways does a genome yield?
3. Can current, future, conditional and unpriced value be told apart?
4. How much value survives comparison with usual care?
5. Does beneficiary separation hold?
6. How completely is each panel callable — and where is it not?
7. Can results aggregate across a cohort without double counting, and without
   presenting individual estimates as population return?
8. Which evidence gaps most often prevent valuation?

## Shape

```
INPUT    deidentified WGS-derived variant data on a supported build,
         with a study identifier; optional clinical context where permitted

PROCESS  genomic interpretation → decision-pathway mapping →
         usual-care comparison → economic attribution →
         evidence qualification

OUTPUT   per-genome economic intelligence
         cohort-level pathway and economic-state summaries
         callability and evidence coverage
         an evidence-gap map naming what would unlock each unvalued pathway
```

## What success means

Not a large dollar figure. Processing success rate; module callability;
pathway traceability; the share of pathways carrying complete economic
semantics; the share of valued pathways naming a beneficiary and a comparator;
a duplicate-benefit rate of zero; an invented-value rate of zero;
representation-equivalence; reconciliation; zero privacy violations in
shareable output.

## The most useful output

The evidence-gap map. Each unvalued pathway names the specific parameter,
identity or attribution problem blocking it — which makes each one a candidate
for a study, a data partnership, or a linkage a sequencing provider may already
be positioned to supply.

## Synthetic demonstrations

[Partner overview](../assets/partner-overview.pdf) ·
[Pilot concept](../assets/pilot-onepager.pdf) ·
[Example payload](../examples/genomelens_partner_v1.example.json)

All generated from synthetic genomes. No prevalence shown in them is an
epidemiological estimate.

## Limitations

Research and pilot evaluation. Not diagnostic. Not a payer submission. Not a
regulatory evaluation. No claim of population return on investment. No current
partner deployment.
